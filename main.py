from fastapi import FastAPI, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from sqlalchemy import or_
from contextlib import asynccontextmanager  # <--- Thư viện quản lý ngữ cảnh mới
import models
import database
import genome

# --- CẤU HÌNH ---
# Dùng file gốc không nén (vì bạn không cài được pysam/bgzip)
FASTA_PATH = "data/R570.fasta"
genome_reader = None


# --- LIFESPAN (THAY THẾ ON_EVENT) ---
@asynccontextmanager
async def lifespan(app: FastAPI):
    # 1. Code chạy khi Server KHỞI ĐỘNG (Startup)
    global genome_reader
    try:
        print(f"🔄 Đang tải dữ liệu Genome từ: {FASTA_PATH}...")
        genome_reader = genome.GenomeReader(FASTA_PATH)
        print("✅ Genome data loaded successfully!")
    except Exception as e:
        print(f"❌ Lỗi tải Genome: {e}")
        # Không kill app, nhưng sẽ báo lỗi nếu gọi API sequence

    yield  # Điểm phân cách giữa Bật và Tắt

    # 2. Code chạy khi Server TẮT (Shutdown)
    print("🛑 Server đang tắt. Dọn dẹp tài nguyên...")
    # Nếu thư viện có hàm close() thì gọi ở đây. pyfaidx tự đóng nên không cần.


# --- KHỞI TẠO APP VỚI LIFESPAN ---
app = FastAPI(
    title="Sugarcane R570 Genome API",
    lifespan=lifespan  # <--- Đăng ký lifespan vào đây
)


# --- API ENDPOINTS ---

@app.get("/")
def read_root():
    return {"status": "Online", "system": "Sugarcane Genome R570"}


@app.get("/genes/search")
def search_genes(
        q: str = Query(None, description="Tìm theo từ khóa (ID, mô tả)"),
        chrom: str = Query(None, description="Tên nhiễm sắc thể (VD: Sh_205k03)"),
        start: int = None,
        end: int = None,
        limit: int = 10,
        db: Session = Depends(database.get_db)
):
    query = db.query(models.Gene)

    # 1. Nếu có tên nhiễm sắc thể -> Filter ngay
    if chrom:
        query = query.filter(models.Gene.chromosome == chrom)

    # 2. Nếu có tọa độ -> Filter tiếp (Tìm gen nằm đè lên vùng này)
    if start and end:
        # Logic giao thoa (Overlap):
        # (Gen.Start <= Vùng.End) AND (Gen.End >= Vùng.Start)
        query = query.filter(
            models.Gene.start <= end,
            models.Gene.end >= start
        )

    # 3. Nếu có từ khóa -> Filter tiếp
    if q:
        search_fmt = f"%{q}%"
        query = query.filter(or_(
            models.Gene.gene_id.like(search_fmt),
            models.Gene.description.like(search_fmt)
        ))

    results = query.limit(limit).all()
    return {"count": len(results), "data": results}


@app.get("/genes/sequence")
def get_gene_sequence(
        gene_id: str = Query(..., description="Nhập ID của gen vào đây (VD: Sh_205k03_g000010)"),
        db: Session = Depends(database.get_db)
):
    """
    Lấy trình tự DNA.
    Input: Query Param ?gene_id=...
    """
    # 1. Tìm thông tin gen trong Database
    gene = db.query(models.Gene).filter(models.Gene.gene_id == gene_id).first()

    if not gene:
        raise HTTPException(status_code=404, detail=f"Gene ID '{gene_id}' not found")

    # 2. Lấy trình tự từ file FASTA
    if not genome_reader:
        raise HTTPException(status_code=500, detail="Genome system not ready")

    seq = genome_reader.get_sequence(gene.chromosome, gene.start, gene.end)

    if not seq:
        raise HTTPException(status_code=404, detail="Sequence not found in file")

    return {
        "gene_id": gene.gene_id,
        "location": f"{gene.chromosome}:{gene.start}-{gene.end}",
        "length": len(seq),
        "sequence": seq
    }