import sys
import os
import argparse  # <--- Thêm thư viện này để nhận tham số dòng lệnh
import urllib.parse

sys.path.append(os.getcwd())
from database import SessionLocal, engine
from models import Base, Gene, Genome  # Nhớ import model Genome


def run_import(gff_path, genome_id, fasta_path):
    print(f"🚀 Bắt đầu import cho bộ gen: {genome_id}")

    # 1. Tạo bảng (Nếu chưa có)
    Base.metadata.create_all(bind=engine)
    session = SessionLocal()

    # 2. Đăng ký Genome vào bảng 'genomes' nếu chưa có
    existing_genome = session.query(Genome).filter(Genome.id == genome_id).first()
    if not existing_genome:
        print(f"➕ Đang tạo mới Genome Metadata: {genome_id}")
        new_genome = Genome(id=genome_id, name=f"Genome {genome_id}", fasta_path=fasta_path)
        session.add(new_genome)
        session.commit()

    # 3. Đọc file GFF và nạp Gene
    print(f"📂 Đang đọc file GFF: {gff_path}...")
    batch = []
    count = 0

    try:
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'gene': continue

                attr_str = parts[8]
                attrs = {}
                for item in attr_str.strip().split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        attrs[k] = urllib.parse.unquote(v)

                gene = Gene(
                    gene_id=attrs.get('ID', attrs.get('Name', f'unknown_{count}')),
                    genome_id=genome_id,  # <--- GÁN GENOME ID VÀO ĐÂY
                    chromosome=parts[0],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    strand=parts[6],
                    description=attrs.get('Note', attrs.get('description', ''))
                )
                batch.append(gene)
                count += 1

                if len(batch) >= 5000:
                    session.bulk_save_objects(batch)
                    session.commit()
                    batch = []
                    print(f"   -> Đã import {count} gen...")

        if batch:
            session.bulk_save_objects(batch)
            session.commit()
        print(f"✅ HOÀN TẤT! Tổng cộng {count} gen cho bộ {genome_id}.")

    except FileNotFoundError:
        print(f"❌ Không tìm thấy file: {gff_path}")
    finally:
        session.close()


if __name__ == "__main__":
    # Cấu hình tham số dòng lệnh
    parser = argparse.ArgumentParser(description="Import GFF data for a specific genome")
    parser.add_argument("--genome", required=True, help="ID của bộ gen (VD: R570)")
    parser.add_argument("--gff", required=True, help="Đường dẫn file GFF3")
    parser.add_argument("--fasta", required=True, help="Đường dẫn file FASTA")

    args = parser.parse_args()

    run_import(args.gff, args.genome, args.fasta)