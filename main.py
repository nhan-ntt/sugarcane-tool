from typing import List, Optional
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
# KHÔNG dùng genome_reader nữa
genome_manager = genome.GenomeManager()


# Hàm phụ trợ kiểm tra bảng tồn tại
def engine_has_table(table_name):
    from sqlalchemy import inspect
    ins = inspect(database.engine)
    return ins.has_table(table_name)


# --- LIFESPAN (QUẢN LÝ VÒNG ĐỜI SERVER) ---
@asynccontextmanager
async def lifespan(app: FastAPI):
    print("🔄 [STARTUP] Đang khởi động hệ thống...")

    # 1. Kết nối DB để lấy danh sách các Genome đã đăng ký
    db = database.SessionLocal()
    try:
        # Kiểm tra xem bảng genomes đã có chưa
        if engine_has_table("genomes"):
            genomes = db.query(models.Genome).all()
            print(f"📂 Tìm thấy {len(genomes)} bộ gen trong Database.")

            # 2. Load từng file FASTA vào RAM (Index)
            for g in genomes:
                print(f"   -> Loading: {g.id}")
                # Load đủ 3 loại file: Genomic, CDS, Protein
                genome_manager.load_genome(g.id, g.fasta_path, g.cds_path, g.protein_path)
        else:
            print("⚠️ Bảng 'genomes' chưa tồn tại. Vui lòng chạy script import_data.py trước.")

    except Exception as e:
        print(f"❌ Lỗi khởi động: {e}")
    finally:
        db.close()

    yield  # --- Server chạy tại đây ---

    print("🛑 [SHUTDOWN] Server đang tắt. Giải phóng tài nguyên...")


# --- KHỞI TẠO APP ---
app = FastAPI(
    title="Sugarcane Multi-Genome API",
    version="2.0",
    lifespan=lifespan
)

# Cấu hình CORS
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
        "loaded_genomes": list(genome_manager.datasets.keys())
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

    # 2. Lọc theo Nhiễm sắc thể
    if chrom:
        query = query.filter(models.Gene.chromosome == chrom)

    # 3. Lọc theo Vùng
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
        genome: str = Query(..., description="ID bộ gen"),
        gene_id: str = Query(..., description="ID của gen"),
        type: str = Query("genomic", description="genomic, cds, protein, flank"),
        db: Session = Depends(database.get_db)
):
    """
    Lấy trình tự DNA/Protein lẻ tẻ.
    """
    gene = db.query(models.Gene).filter(
        models.Gene.genome_id == genome,
        models.Gene.gene_id == gene_id
    ).first()

    if not gene: raise HTTPException(404, detail="Not found")

    seq = genome_manager.get_data(genome, type, gene.gene_id, gene.chromosome, gene.start, gene.end)
    return {"genome": genome, "gene": gene_id, "type": type, "sequence": seq}


@app.get("/genome/gene_detail")
def get_gene_detail(
        genome: str = Query(..., description="ID bộ gen"),
        gene_id: str = Query(..., description="ID của gen"),
        db: Session = Depends(database.get_db)
):
    """
    Lấy TRỌN BỘ thông tin (Full Detail) cho trang chi tiết.
    """
    gene = db.query(models.Gene).filter(
        models.Gene.genome_id == genome,
        models.Gene.gene_id == gene_id
    ).first()

    if not gene: raise HTTPException(404, detail="Gene not found")

    # Lấy annotation
    try:
        anno = db.query(models.Annotation).filter(models.Annotation.gene_id == gene_id).first()
    except:
        anno = None

    # Lấy sequences
    seq_genomic = genome_manager.get_data(genome, 'genomic', gene.gene_id, gene.chromosome, gene.start, gene.end)
    seq_cds = genome_manager.get_data(genome, 'cds', gene.gene_id)
    seq_protein = genome_manager.get_data(genome, 'protein', gene.gene_id)
    seq_flank = genome_manager.get_data(genome, 'flank', gene.gene_id, gene.chromosome, gene.start, gene.end)

    return {
        "basic_info": {
            "gene_id": gene.gene_id,
            "genome": genome,
            "chromosome": gene.chromosome,
            "location": f"{gene.start}-{gene.end}",
            "strand": gene.strand,
            "description": gene.description
        },
        "sequences": {
            "genomic": seq_genomic,
            "cds": seq_cds,
            "protein": seq_protein,
            "flank_upstream": seq_flank
        },
        "annotations": {
            "swissprot": anno.seed_ortholog if anno else None,
            "go_terms": anno.go_terms if anno else None,
            "kegg_pathways": anno.kegg_pathways if anno else None,
            "pfam": anno.pfam_domains if anno else None,
            "eggnog_desc": anno.description if anno else None
        }
    }


# --- BATCH SEQUENCE API ---
class BatchSequenceRequest(BaseModel):
    genome: str
    gene_ids: List[str]


@app.post("/genome/sequence/batch")
def get_sequence_batch_post(
        payload: BatchSequenceRequest,
        db: Session = Depends(database.get_db)
):
    genes = db.query(models.Gene).filter(
        models.Gene.genome_id == payload.genome,
        models.Gene.gene_id.in_(payload.gene_ids)
    ).all()

    results = []
    found_ids = set()

    for gene in genes:
        found_ids.add(gene.gene_id)
        try:
            seq = genome_manager.get_data(payload.genome, 'genomic', gene.gene_id, gene.chromosome, gene.start,
                                          gene.end)
            results.append({
                "gene_id": gene.gene_id,
                "found": True,
                "location": f"{gene.chromosome}:{gene.start}-{gene.end}",
                "length": len(seq) if seq else 0,
                "sequence": seq
            })
        except Exception as e:
            results.append({"gene_id": gene.gene_id, "found": False, "error": str(e)})

    # Check missing
    for mid in set(payload.gene_ids) - found_ids:
        results.append({"gene_id": mid, "found": False, "error": "Not in DB"})

    return {
        "genome": payload.genome,
        "total_requested": len(payload.gene_ids),
        "total_found": len(genes),
        "data": results
    }


# --- CRISPOR TOOL (QUAN TRỌNG) ---
@app.post("/tools/crispor")
def run_crispor_tool(
        genome: str = Query(..., description="Chọn bộ gen"),
        gene_id: str = None,
        sequence: str = None,
        db: Session = Depends(database.get_db)
):
    """
    CRISPOR Engine: Tìm gRNA + Off-target Bowtie2 + Primer3
    """
    # 1. Tìm đường dẫn Index (WSL Path)
    genome_info = db.query(models.Genome).filter(models.Genome.id == genome).first()
    if not genome_info:
        raise HTTPException(404, detail=f"Genome '{genome}' chưa được hỗ trợ.")

    # Xử lý đường dẫn Windows -> WSL cho Bowtie2
    abs_fasta_path = os.path.abspath(genome_info.fasta_path)
    abs_index_path = abs_fasta_path.replace(".fasta", "_index")  # Quy ước tên index

    wsl_path = abs_index_path.replace("\\", "/")
    if wsl_path.lower().startswith("d:"):
        wsl_path = "/mnt/d" + wsl_path[2:]
    elif wsl_path.lower().startswith("c:"):
        wsl_path = "/mnt/c" + wsl_path[2:]

    print(f"DEBUG: WSL Index Path -> {wsl_path}")

    # 2. Lấy sequence
    target_seq = ""
    if gene_id:
        gene = db.query(models.Gene).filter(models.Gene.genome_id == genome, models.Gene.gene_id == gene_id).first()
        if not gene: raise HTTPException(404, "Gene not found")

        # Lấy rộng ra 100bp để thiết kế Primer
        try:
            # Lưu ý: get_data trả về string, ta cần tính lại tọa độ để lấy padding
            # Hoặc gọi trực tiếp manager nếu cần tùy biến
            # Ở đây ta gọi manager để lấy padding
            target_seq = genome_manager.get_data(genome, 'genomic', chrom=gene.chromosome, start=gene.start - 100,
                                                 end=gene.end + 100)
        except Exception as e:
            raise HTTPException(500, detail=f"Lỗi đọc Fasta: {e}")

    elif sequence:
        target_seq = sequence
    else:
        raise HTTPException(400, "Thiếu input gene_id hoặc sequence")

    # 3. Chạy Engine
    try:
        results = crispor_engine.run_crispor_analysis(str(target_seq), wsl_path)
    except Exception as e:
        print(f"Lỗi Engine: {e}")
        results = []

    return {
        "genome": genome,
        "index_used": wsl_path,
        "input_length": len(target_seq) if target_seq else 0,
        "guides_found": len(results),
        "top_guides": results[:20]
    }