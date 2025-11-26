from pyfaidx import Fasta
import os


class GenomeReader:
    def __init__(self, fasta_path: str):
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"Missing genome file: {fasta_path}")

        # Pyfaidx tự động xử lý file index (.fai)
        # Nếu chưa có index, nó sẽ tự tạo (mất vài giây lần đầu)
        self.fasta = Fasta(fasta_path)

    def get_sequence(self, chrom: str, start: int, end: int):
        """
        Lấy trình tự DNA.
        Input: Tọa độ 1-based (theo GFF/Sinh học)
        """
        try:
            # Pyfaidx cho phép slicing như chuỗi Python: [start:end]
            # Lưu ý: Pyfaidx dùng 0-based indexing nội bộ nhưng API rất linh hoạt.
            # Cách an toàn nhất để khớp với GFF (1-based, inclusive) là:
            seq = self.fasta[chrom][start - 1:end]
            return str(seq)
        except KeyError:
            return None  # Không tìm thấy Chromosome
        except Exception as e:
            print(f"Error reading sequence: {e}")
            return None