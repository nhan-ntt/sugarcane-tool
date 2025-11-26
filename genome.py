from pyfaidx import Fasta
import os


class GenomeManager:
    def __init__(self):
        # Dictionary lưu các object Fasta: {'R570': FastaObject, 'AP85': FastaObject}
        self.readers = {}

    def load_genome(self, genome_id: str, fasta_path: str):
        """Đăng ký một bộ gen mới vào hệ thống"""
        if os.path.exists(fasta_path):
            # Pyfaidx chỉ load index nhẹ, không tốn RAM
            self.readers[genome_id] = Fasta(fasta_path)
            print(f"✅ Loaded genome: {genome_id}")
        else:
            print(f"❌ File not found: {fasta_path}")

    def get_sequence(self, genome_id: str, chrom: str, start: int, end: int):
        """Lấy sequence từ bộ gen cụ thể"""
        if genome_id not in self.readers:
            raise ValueError(f"Genome {genome_id} chưa được load.")

        try:
            return str(self.readers[genome_id][chrom][start - 1:end])
        except KeyError:
            return None