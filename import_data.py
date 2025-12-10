import sys
import os
import argparse
import urllib.parse

# Thêm thư mục hiện tại vào sys.path để import được module trong app
sys.path.append(os.getcwd())

# Import từ app
from database import SessionLocal, engine
from models import Base, Gene, Genome


def run_import(genome_id, gff_path, fasta_path, cds_path=None, protein_path=None):
    """
    Hàm import dữ liệu gen và metadata genome.
    """
    print(f"🚀 Bắt đầu import cho bộ gen: {genome_id}")

    # 1. Khởi tạo Database
    Base.metadata.create_all(bind=engine)
    session = SessionLocal()

    try:
        # 2. Xử lý thông tin Genome (Metadata)
        existing_genome = session.query(Genome).filter(Genome.id == genome_id).first()

        if not existing_genome:
            print(f"➕ Đang tạo mới Genome Metadata: {genome_id}")
            new_genome = Genome(
                id=genome_id,
                name=f"Genome {genome_id}",
                fasta_path=fasta_path,  # Đường dẫn Genomic (.fna)
                cds_path=cds_path,  # Đường dẫn CDS (.cds.fna) - Mới
                protein_path=protein_path  # Đường dẫn Protein (.faa) - Mới
            )
            session.add(new_genome)
        else:
            print(f"ℹ️ Genome {genome_id} đã tồn tại. Đang cập nhật đường dẫn file...")
            # Cập nhật lại đường dẫn nếu chạy lại script
            existing_genome.fasta_path = fasta_path
            if cds_path: existing_genome.cds_path = cds_path
            if protein_path: existing_genome.protein_path = protein_path

        session.commit()

        # 3. Đọc file GFF và nạp dữ liệu Gene
        print(f"📂 Đang đọc file GFF: {gff_path}...")

        # Xóa dữ liệu cũ của genome này (để tránh duplicate nếu import lại)
        deleted_rows = session.query(Gene).filter(Gene.genome_id == genome_id).delete()
        if deleted_rows > 0:
            print(f"🧹 Đã dọn dẹp {deleted_rows} gen cũ của {genome_id} trước khi import mới.")
            session.commit()

        batch = []
        count = 0

        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')

                # Chỉ xử lý dòng gene
                if len(parts) < 9 or parts[2] != 'gene': continue

                attr_str = parts[8]
                attrs = {}
                for item in attr_str.strip().split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        attrs[k] = urllib.parse.unquote(v)

                # Ưu tiên lấy ID, nếu không có lấy Name
                g_id = attrs.get('ID', attrs.get('Name', f'unknown_{count}'))

                gene = Gene(
                    gene_id=g_id,
                    genome_id=genome_id,
                    chromosome=parts[0],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    strand=parts[6],
                    description=attrs.get('Note', attrs.get('description', ''))
                )
                batch.append(gene)
                count += 1

                # Bulk Insert mỗi 5000 dòng
                if len(batch) >= 5000:
                    session.bulk_save_objects(batch)
                    session.commit()
                    batch = []
                    print(f"   -> Đã import {count} gen...")

        # Commit phần còn lại
        if batch:
            session.bulk_save_objects(batch)
            session.commit()

        print(f"✅ HOÀN TẤT! Tổng cộng {count} gen đã được lưu vào Database.")

    except FileNotFoundError as e:
        print(f"❌ Lỗi: Không tìm thấy file - {e}")
    except Exception as e:
        print(f"❌ Lỗi hệ thống: {e}")
        session.rollback()
    finally:
        session.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tool Import dữ liệu Gen mía vào Database")

    # Các tham số bắt buộc
    parser.add_argument("--genome", required=True, help="ID định danh bộ gen (VD: R570, E_coli)")
    parser.add_argument("--gff", required=True, help="Đường dẫn file GFF3")
    parser.add_argument("--fasta", required=True, help="Đường dẫn file Genomic FASTA (.fna)")

    # Các tham số tùy chọn (Mới thêm)
    parser.add_argument("--cds", help="Đường dẫn file CDS FASTA (.cds.fna)", default=None)
    parser.add_argument("--protein", help="Đường dẫn file Protein FASTA (.faa)", default=None)

    args = parser.parse_args()

    run_import(
        genome_id=args.genome,
        gff_path=args.gff,
        fasta_path=args.fasta,
        cds_path=args.cds,
        protein_path=args.protein
    )