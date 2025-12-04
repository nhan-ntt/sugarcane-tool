from typing import List

from fastapi import FastAPI, Depends, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from sqlalchemy.orm import Session
from sqlalchemy import or_
from contextlib import asynccontextmanager
import os

# Import các module nội bộ
import models
import database
import genome
import crispor_engine

# --- CẤU HÌNH TOÀN CỤC ---
# Khởi tạo Manager để quản lý nhiều bộ gen cùng lúc
genome_manager = genome.GenomeManager()


# --- LIFESPAN (QUẢN LÝ VÒNG ĐỜI SERVER) ---
@asynccontextmanager
async def lifespan(app: FastAPI):
    print("🔄 [STARTUP] Đang khởi động hệ thống...")

    # 1. Kết nối DB để lấy danh sách các Genome đã đăng ký
    db = database.SessionLocal()
    try:
        # Kiểm tra xem bảng genomes đã có chưa, nếu chưa thì bỏ qua
        # (Tránh lỗi lần đầu chạy chưa import)
        if engine_has_table("genomes"):
            genomes = db.query(models.Genome).all()
            print(f"📂 Tìm thấy {len(genomes)} bộ gen trong Database.")

            # 2. Load từng file FASTA vào RAM (Index)
            for g in genomes:
                print(f"   -> Loading: {g.id} ({g.fasta_path})")
                genome_manager.load_genome(g.id, g.fasta_path)
        else:
            print("⚠️ Bảng 'genomes' chưa tồn tại. Vui lòng chạy script import_data.py trước.")

    except Exception as e:
        print(f"❌ Lỗi khởi động: {e}")
    finally:
        db.close()

    yield  # --- Server chạy tại đây ---

    print("🛑 [SHUTDOWN] Server đang tắt. Giải phóng tài nguyên...")
    # Pyfaidx tự động đóng file handle nên không cần code thêm


# Hàm phụ trợ kiểm tra bảng
def engine_has_table(table_name):
    from sqlalchemy import inspect
    ins = inspect(database.engine)
    return ins.has_table(table_name)


# --- KHỞI TẠO APP ---
app = FastAPI(
    title="Sugarcane Multi-Genome API",
    version="2.0",
    lifespan=lifespan
)

# Cấu hình CORS (Cho phép Frontend gọi vào)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# --- API ENDPOINTS ---

@app.get("/")
def read_root():
    return {
        "status": "Online",
        "system": "Sugarcane Multi-Genome System",
        "loaded_genomes": list(genome_manager.readers.keys())
    }


@app.get("/genome/search")
def search_genes(
        genome: str = Query(..., description="ID bộ gen (VD: R570, AP85)"),
        q: str = Query(None, description="Từ khóa: ID gen hoặc mô tả"),
        chrom: str = Query(None, description="Tên nhiễm sắc thể"),
        start: int = None,
        end: int = None,
        limit: int = 10,
        db: Session = Depends(database.get_db)
):
    """
    Tìm kiếm gen trong một bộ gen cụ thể.
    """
    # 1. Lọc theo Genome (Bắt buộc)
    query = db.query(models.Gene).filter(models.Gene.genome_id == genome)

    # 2. Lọc theo Chromosome
    if chrom:
        query = query.filter(models.Gene.chromosome == chrom)

    # 3. Lọc theo Vùng (Overlap logic)
    if start and end:
        query = query.filter(
            models.Gene.start <= end,
            models.Gene.end >= start
        )

    # 4. Lọc theo Từ khóa
    if q:
        search_fmt = f"%{q}%"
        query = query.filter(or_(
            models.Gene.gene_id.like(search_fmt),
            models.Gene.description.like(search_fmt)
        ))

    results = query.limit(limit).all()

    return {
        "genome": genome,
        "count": len(results),
        "data": results
    }


@app.get("/genome/sequence")
def get_sequence(
        genome: str = Query(..., description="ID bộ gen (VD: R570)"),
        gene_id: str = Query(..., description="ID của gen"),
        db: Session = Depends(database.get_db)
):
    """
    Lấy trình tự DNA của gen từ file FASTA tương ứng.
    """
    # 1. Tìm thông tin trong DB
    gene = db.query(models.Gene).filter(
        models.Gene.genome_id == genome,
        models.Gene.gene_id == gene_id
    ).first()

    if not gene:
        raise HTTPException(404, detail=f"Không tìm thấy gen '{gene_id}' trong bộ '{genome}'")

    # 2. Lấy sequence từ Manager
    # Lưu ý: Manager tự chọn đúng file FASTA dựa trên `genome` ID
    try:
        seq = genome_manager.get_sequence(genome, gene.chromosome, gene.start, gene.end)
    except ValueError as e:
        raise HTTPException(500, detail=str(e))  # Lỗi chưa load genome

    if not seq:
        raise HTTPException(404, detail="Không đọc được sequence từ file (Check tọa độ/Fasta)")

    return {
        "genome": genome,
        "gene_id": gene.gene_id,
        "location": f"{gene.chromosome}:{gene.start}-{gene.end}",
        "length": len(seq),
        "sequence": seq
    }


# Định nghĩa Body cho Request
class BatchSequenceRequest(BaseModel):
    genome: str
    gene_ids: List[str]


@app.post("/genome/sequence/batch")
def get_sequence_batch_post(
        payload: BatchSequenceRequest,
        db: Session = Depends(database.get_db)
):
    """
    Lấy chi tiết sequence cho danh sách gene (Batch).
    Input: JSON Body { "genome": "R570", "gene_ids": ["ID1", "ID2"] }
    """
    # 1. Truy vấn Database 1 lần duy nhất (Tối ưu SQL)
    genes = db.query(models.Gene).filter(
        models.Gene.genome_id == payload.genome,
        models.Gene.gene_id.in_(payload.gene_ids)
    ).all()

    results = []
    found_ids = set()

    # 2. Duyệt qua các gen tìm thấy để lấy Sequence từ file
    for gene in genes:
        found_ids.add(gene.gene_id)
        try:
            # Lấy sequence từ File FASTA thông qua Manager
            seq = genome_manager.get_sequence(payload.genome, gene.chromosome, gene.start, gene.end)

            results.append({
                "gene_id": gene.gene_id,
                "found": True,
                "location": f"{gene.chromosome}:{gene.start}-{gene.end}",
                "length": len(seq) if seq else 0,
                "sequence": seq  # <--- Dữ liệu quan trọng nhất đây
            })
        except Exception as e:
            results.append({
                "gene_id": gene.gene_id,
                "found": False,
                "error": f"Lỗi đọc file: {str(e)}"
            })

    # 3. Báo cáo các ID không tìm thấy (Missing)
    # Để client biết ID nào bị sai
    requested_set = set(payload.gene_ids)
    missing_ids = requested_set - found_ids

    for mid in missing_ids:
        results.append({
            "gene_id": mid,
            "found": False,
            "error": "Gene ID not found in Database"
        })

    # 4. Trả về kết quả chi tiết
    return {
        "genome": payload.genome,
        "total_requested": len(payload.gene_ids),
        "total_found": len(genes),
        "data": results
    }


@app.post("/tools/crispor")
def run_crispor_tool(
        genome: str = Query(..., description="Chọn bộ gen để lấy sequence (VD: R570)"),
        gene_id: str = None,
        sequence: str = None,
        db: Session = Depends(database.get_db)
):
    """
    Chạy công cụ CRISPOR (Doench '16 + CFD + Primer3).
    Hỗ trợ Multi-Genome với đường dẫn Index động.
    """

    # 1. TÌM ĐƯỜNG DẪN INDEX CỦA GENOME NÀY (BƯỚC QUAN TRỌNG MỚI)
    # Lấy thông tin bộ gen từ DB
    genome_info = db.query(models.Genome).filter(models.Genome.id == genome).first()
    if not genome_info:
        raise HTTPException(404, detail=f"Genome '{genome}' chưa được hỗ trợ hoặc chưa import.")

    # Logic: Biến đường dẫn FASTA thành đường dẫn INDEX
    # Giả sử DB lưu: "data/R570/R570.fasta"
    # Ta cần tạo ra: "/mnt/d/.../data/R570/R570_index"
    abs_fasta_path = os.path.abspath(genome_info.fasta_path)

    # 2. Chuyển đổi sang đường dẫn Index
    # VD: ...\R570.fasta -> ...\R570_index
    abs_index_path = abs_fasta_path.replace(".fasta", "_index")

    # 3. "Phiên dịch" sang chuẩn WSL (Linux)
    # Bước A: Đổi dấu gạch chéo ngược (\) thành xuôi (/)
    wsl_path = abs_index_path.replace("\\", "/")

    # Bước B: Đổi ổ đĩa (D: -> /mnt/d)
    # Lưu ý: Sửa chữ cái ổ đĩa cho đúng máy bạn (thường là c: hoặc d:)
    if wsl_path.lower().startswith("d:"):
        wsl_path = "/mnt/d" + wsl_path[2:]
    elif wsl_path.lower().startswith("c:"):
        wsl_path = "/mnt/c" + wsl_path[2:]

    print(f"DEBUG: Đường dẫn WSL chuẩn -> {wsl_path}")  # In ra để kiểm tra

    # 2. XÁC ĐỊNH SEQUENCE MỤC TIÊU
    target_seq = ""

    # Case A: Dùng Gene ID
    if gene_id:
        gene = db.query(models.Gene).filter(
            models.Gene.genome_id == genome,
            models.Gene.gene_id == gene_id
        ).first()

        if not gene:
            raise HTTPException(404, "Gene not found")

        # Lấy rộng ra +/- 100bp để thiết kế Primer
        padding = 100
        try:
            target_seq = genome_manager.get_sequence(
                genome,
                gene.chromosome,
                gene.start - padding,
                gene.end + padding
            )
        except Exception as e:
            raise HTTPException(500, detail=f"Lỗi đọc Fasta từ Manager: {e}")

    # Case B: Dùng Sequence thô
    elif sequence:
        target_seq = sequence
    else:
        raise HTTPException(400, "Cần cung cấp gene_id hoặc sequence")

    # 3. CHẠY ENGINE VỚI THAM SỐ MỚI
    # Truyền thêm real_index_path vào đây
    try:
        results = crispor_engine.run_crispor_analysis(str(target_seq), wsl_path)
    except Exception as e:
        print(f"Lỗi Engine: {e}")
        # Nếu lỗi (ví dụ chưa build index), trả về list rỗng nhưng không crash
        results = []

    return {
        "genome": genome,
        "index_used": wsl_path,  # Trả về để bạn debug xem đường dẫn đúng không
        "gene_id": gene_id,
        "input_length": len(target_seq),
        "guides_found": len(results),
        "top_guides": results[:20]
    }