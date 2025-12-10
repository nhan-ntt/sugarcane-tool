from pyfaidx import Fasta
import os
import re


class GenomeManager:
    def __init__(self):
        self.datasets = {}

    def load_genome(self, genome_id: str, fasta_path: str, cds_path: str = None, protein_path: str = None):
        """
        Load genome với bộ lọc header (chỉ lấy phần ID trước dấu cách).
        """
        self.datasets[genome_id] = {}

        # Hàm làm sạch header: ">ID description..." -> "ID"
        clean_key_func = lambda x: x.split()[0].strip()

        if os.path.exists(fasta_path):
            self.datasets[genome_id]['genomic'] = Fasta(fasta_path, key_function=clean_key_func)
            print(f"✅ [{genome_id}] Loaded Genomic")

        if cds_path and os.path.exists(cds_path):
            self.datasets[genome_id]['cds'] = Fasta(cds_path, key_function=clean_key_func)
            print(f"✅ [{genome_id}] Loaded CDS")

        if protein_path and os.path.exists(protein_path):
            self.datasets[genome_id]['protein'] = Fasta(protein_path, key_function=clean_key_func)
            print(f"✅ [{genome_id}] Loaded Protein")

    def _smart_search(self, fasta_obj, gene_id):
        """
        Hàm tìm kiếm siêu thông minh:
        Input: SoffiXsponR570.10Ag000100.v2.1
        Thử lần lượt:
        1. Chính xác: SoffiXsponR570.10Ag000100.v2.1
        2. Bỏ 'v':    SoffiXsponR570.10Ag000100.2.1
        3. Isoform 1: SoffiXsponR570.10Ag000100.1 (Khả năng trúng cao nhất với file của bạn)
        4. Base ID:   SoffiXsponR570.10Ag000100
        """
        # 1. Thử tìm chính xác 100%
        if gene_id in fasta_obj:
            return str(fasta_obj[gene_id])

        # 2. Tìm Base ID (ID gốc)
        # Logic: Cắt bỏ các phần đuôi .v2.1, .1, .t1 ... để lấy cái gốc "SoffiXsponR570.10Ag000100"
        # Regex này tìm chuỗi trước các dấu chấm suffix cuối cùng
        # VD: ID.v2.1 -> ID
        # VD: ID.1    -> ID
        base_match = re.match(r'^(.*?)\.(v?\d+(\.\d+)*)$', gene_id)

        if base_match:
            base_id = base_match.group(1)  # Lấy phần "SoffiXsponR570.10Ag000100"

            # Chiến thuật A: Thử Base ID trần
            if base_id in fasta_obj:
                return str(fasta_obj[base_id])

            # Chiến thuật B: Thử ghép đuôi phổ biến (.1, .2, .t1)
            # Vì file Fasta của bạn có dạng ID.1, ID.2 -> ta thử ID.1 trước (phổ biến nhất)
            candidates = [
                f"{base_id}.1",  # Thử .1 (Isoform 1)
                f"{base_id}.2",  # Thử .2
                f"{base_id}.3",
                f"{base_id}.v1",  # Thử version
                f"{base_id}.v2",
                gene_id.replace(".v", ".")  # Thử bỏ chữ 'v': .v2.1 -> .2.1
            ]

            for cand in candidates:
                if cand in fasta_obj:
                    return str(fasta_obj[cand])

        # 3. Fallback cuối cùng: Nếu ID quá lạ, thử cắt dần từ đuôi lên
        # ID.a.b.c -> ID.a.b -> ID.a
        temp_id = gene_id
        while '.' in temp_id:
            temp_id = temp_id.rsplit('.', 1)[0]
            if temp_id in fasta_obj:
                return str(fasta_obj[temp_id])

        # Nếu vẫn không thấy thì chịu thua
        raise KeyError(f"ID '{gene_id}' not found (even with smart search).")

    def get_data(self, genome_id: str, type: str, gene_id: str = None, chrom: str = None, start: int = 0, end: int = 0):
        if genome_id not in self.datasets:
            return None
        dataset = self.datasets[genome_id]

        try:
            # --- Genomic & Flank ---
            if type == 'genomic':
                return str(dataset['genomic'][chrom][start - 1:end])

            elif type == 'flank':
                flank_len = 2000
                p_end = start - 1
                p_start = max(0, p_end - flank_len)
                return str(dataset['genomic'][chrom][p_start:p_end])

            # --- CDS & Protein (Dùng Smart Search) ---
            elif type == 'cds':
                if 'cds' not in dataset: return None
                try:
                    return self._smart_search(dataset['cds'], gene_id)
                except KeyError:
                    return None

            elif type == 'protein':
                if 'protein' not in dataset: return None
                try:
                    return self._smart_search(dataset['protein'], gene_id)
                except KeyError:
                    return None

        except Exception as e:
            # print(f"Error extracting sequence: {e}") # Bật cái này nếu muốn debug log
            return None

    def get_sequence(self, genome_id, chrom, start, end):
        return self.get_data(genome_id, 'genomic', chrom=chrom, start=start, end=end)