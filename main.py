from fastapi import FastAPI, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from sqlalchemy import or_
from contextlib import asynccontextmanager  # <--- Thư viện quản lý ngữ cảnh mới
import models
import database
import genome
import crispr
import crispor_engine

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


@app.post("/tools/crispr/design")
def design_crispr_guides(
        sequence: str = None,
        gene_id: str = None,
        db: Session = Depends(database.get_db)
):
    """
    Công cụ thiết kế gRNA.
    Người dùng có thể nhập trực tiếp Sequence HOẶC nhập Gene ID để hệ thống tự lấy sequence.
    """
    target_seq = ""

    # Case 1: Nhập Gene ID
    if gene_id:
        gene = db.query(models.Gene).filter(models.Gene.gene_id == gene_id).first()
        if not gene:
            raise HTTPException(404, "Gene ID not found")
        if not genome_reader:
            raise HTTPException(500, "Genome system not ready")
        target_seq = genome_reader.get_sequence(gene.chromosome, gene.start, gene.end)

    # Case 2: Nhập Sequence trực tiếp
    elif sequence:
        target_seq = sequence
    else:
        raise HTTPException(400, "Phải cung cấp sequence hoặc gene_id")

    # Gọi thuật toán tìm target
    results = crispr.find_crispr_targets(str(target_seq))

    return {
        "gene_id": gene_id,
        "input_length": len(target_seq),
        "candidates_found": len(results),
        "guides": results
    }


@app.post("/tools/crispor")
def run_crispor_tool(
        gene_id: str = None,
        sequence: str = None,
        db: Session = Depends(database.get_db)
):
    """
    Chạy thuật toán mô phỏng CRISPOR.
    Input: Gene ID hoặc Sequence thô.
    """
    target_seq = ""

    print(f"DEBUG: Nhận được gene_id = '{gene_id}'")

    # 1. Lấy sequence
    if gene_id:
        gene = db.query(models.Gene).filter(models.Gene.gene_id == gene_id).first()
        if not gene: raise HTTPException(404, "Gene not found")

        # Cần lấy rộng ra +/- 100bp để thiết kế Primer
        padding = 100
        target_seq = genome_reader.get_sequence(gene.chromosome, gene.start - padding, gene.end + padding)
    elif sequence:
        target_seq = sequence
    else:
        raise HTTPException(400, "Missing input")

    # 2. Chạy thuật toán
    results = crispor_engine.run_crispor_analysis(str(target_seq))

    return {
        "gene_id": gene_id,
        "input_length": len(target_seq),
        "guides_found": len(results),
        "top_guides": results  # Trả về top 10 tốt nhất
    }