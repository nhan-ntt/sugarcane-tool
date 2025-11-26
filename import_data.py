import sys
import os
import urllib.parse

sys.path.append(os.getcwd())
from database import SessionLocal, engine
from models import Base, Gene


def run_import():
    print("⏳ Creating Database...")
    Base.metadata.create_all(bind=engine)
    session = SessionLocal()

    # ĐƯỜNG DẪN ĐÃ ĐƯỢC CHUẨN HÓA
    gff_path = "data/R570.gff3"

    print(f"📂 Reading {gff_path}...")
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
                    print(f"Imported {count} genes...")

        if batch:
            session.bulk_save_objects(batch)
            session.commit()
        print("✅ DONE!")
    except FileNotFoundError:
        print(f"❌ Error: Could not find {gff_path}. Check your folder structure.")
    finally:
        session.close()


if __name__ == "__main__":
    run_import()