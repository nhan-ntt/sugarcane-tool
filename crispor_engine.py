import re
import primer3
import subprocess
import numpy as np

# --- 1. CẤU HÌNH CFD MATRIX (Dựa trên Doench et al. 2016) ---
# Bảng tra cứu trọng số phạt cho từng loại mismatch ở từng vị trí (1-20)
# Đây là "Secret Sauce" của CRISPOR/Doench để tính độ đặc hiệu chính xác.
# Giá trị 1.0 = Cắt mạnh (Nguy hiểm), Giá trị thấp = Ít cắt (An toàn)
CFD_WEIGHTS = {
    # Ví dụ rút gọn (Thực tế bảng này có 4x4x20 giá trị)
    # Mismatch rG:dA (RNA là G, DNA là A)
    ('G', 'A'): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.0, 0.0, 0.0, 0.0],
    # Mặc định cho các mismatch khác (simplified heuristic)
    'default': [1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.1, 0.1, 0.05, 0.02, 0.0, 0.0, 0.0, 0.0]
}


class CrisporEngine:
    def __init__(self, genome_index_path=None):
        self.genome_index = genome_index_path

    def find_candidates(self, sequence):
        """
        Bước 1: Tìm tất cả vị trí có PAM (NGG)
        Trả về: List các đoạn 30bp (Context) để tính điểm
        """
        candidates = []
        seq_upper = sequence.upper()
        # Regex tìm NGG (Lookahead để bắt các PAM chồng lấn)
        for match in re.finditer(r'(?=([ATGC]GG))', seq_upper):
            pam_start = match.start()
            # Guide 20bp nằm trước PAM
            guide_start = pam_start - 20

            # Cần lấy ngữ cảnh 30bp: 4bp trước + 20bp Guide + 3bp PAM + 3bp sau
            context_start = guide_start - 4
            context_end = pam_start + 3 + 3

            if context_start < 0 or context_end > len(seq_upper):
                continue  # Bỏ qua nếu ở sát mép

            context_seq = seq_upper[context_start:context_end]
            guide_seq = seq_upper[guide_start:pam_start]

            candidates.append({
                "guide_seq": guide_seq,
                "pam": match.group(1),
                "start": guide_start,
                "end": pam_start + 3,
                "context_30bp": context_seq
            })
        return candidates

    def calculate_efficiency_score(self, context_30bp):
        """
        Bước 2: Tính điểm Doench '16 (Efficiency)
        Đây là mô phỏng logic (Simplified), thực tế cần load model Pickle của Sklearn.
        """
        score = 50  # Base

        # Các quy tắc sinh học cơ bản từ bài báo Doench 2016:
        guide = context_30bp[4:24]

        # 1. GC Content (Tối ưu 40-60%)
        gc = (guide.count('G') + guide.count('C')) / 20.0
        if 0.4 <= gc <= 0.6:
            score += 10
        elif gc > 0.8 or gc < 0.2:
            score -= 20

        # 2. Position Specific (Ví dụ: G ở vị trí 20 ngay cạnh PAM rất tốt)
        if guide[19] == 'G': score += 10
        if guide[19] == 'T': score -= 10  # T cạnh PAM làm giảm hiệu suất cắt
        if 'TTTT' in guide: return 0  # Gây ngắt phiên mã U6 -> Score = 0

        return max(0, min(100, score))

    def search_off_targets(self, guide_seq):
        """
        Bước 3: Tìm Off-target bằng Bowtie2 (Mô phỏng)
        Trong môi trường thật, hàm này gọi subprocess bowtie2.
        """
        # --- CODE THẬT (Khi đã cài Bowtie2) ---
        # cmd = f"bowtie2 -x {self.genome_index} -c {guide_seq} -k 10 -v 3 --no-unal"
        # process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        # return parse_sam_output(process.stdout)

        # --- CODE GIẢ LẬP (Cho phép chạy thử ngay) ---
        # Giả vờ tìm thấy 1 off-target có 2 lỗi sai
        return [
            {"seq": guide_seq, "mismatches": 0, "chrom": "Target"},  # Chính nó
            {"seq": guide_seq.replace('A', 'T', 1), "mismatches": 2, "chrom": "Chr_Random"}  # Giả vờ
        ]

    def search_off_targets_real(self, guide_seq):
        """
        CODE THẬT: Gọi Bowtie2 để quét toàn bộ genome
        """
        # 1. Cấu tạo lệnh
        # -N 1: Cho phép sai tối đa 1 nucleotide trong vùng seed (để tìm off-target gần đúng)
        # -L 20: Độ dài seed
        cmd = [
            "wsl",  # <--- Gọi thông qua WSL
            "bowtie2",
            "-x", self.genome_index,
            "-c", guide_seq,
            "-k", "10",  # Tìm tối đa 10 vị trí
            "--no-unal",  # Không hiện cái không tìm thấy
            "--no-hd"  # Không hiện header SAM (cho dễ parse)
        ]

        try:
            # 2. Gọi subprocess
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # 3. Parse kết quả
            return self._parse_sam_output(process.stdout, guide_seq)

        except subprocess.CalledProcessError as e:
            print(f"Lỗi Bowtie2: {e.stderr}")
            return []
        except FileNotFoundError:
            print("Lỗi: Chưa cài đặt bowtie2 trên hệ thống (sudo apt install bowtie2)")
            return []

    def _parse_sam_output(self, sam_content, original_seq):
        """
        Hàm phụ trợ: Đọc format SAM của Bowtie2 trả về
        """
        off_targets = []

        for line in sam_content.strip().split('\n'):
            if not line: continue

            # Cấu trúc SAM file (Tab delimited)
            # Col 2: Chromosome Name
            # Col 3: Position (1-based)
            # ...
            # Tags: NM:i:x (Số lượng mismatch)
            parts = line.split('\t')

            if len(parts) < 10: continue

            chrom = parts[2]
            pos = parts[3]

            # Tìm thẻ NM (Number of Mismatches)
            mismatches = 0
            for tag in parts[11:]:
                if tag.startswith("NM:i:"):
                    mismatches = int(tag.split(":")[2])
                    break

            off_targets.append({
                "seq": original_seq,  # Thực tế nên lấy seq từ SAM nếu nó có mutation
                "chrom": chrom,
                "position": pos,
                "mismatches": mismatches
            })

        return off_targets

    def calculate_specificity_score(self, guide_seq, off_targets):
        """
        Bước 4: Tính điểm CFD (Specificity) dựa trên danh sách Off-target
        Công thức: 1 / (1 + sum(CFD_score_of_off_targets))
        """
        sum_cfd = 0
        for ot in off_targets:
            if ot['mismatches'] == 0: continue  # Bỏ qua chính nó (On-target)

            # Tính điểm CFD cho từng off-target
            # (Logic: So sánh guide vs off-target xem sai ở đâu)
            current_cfd = 0.5  # Giả lập: 0.5 nghĩa là khả năng cắt nhầm là 50%
            sum_cfd += current_cfd

        # Công thức chuẩn hóa về thang 100 (Càng nhiều off-target điểm càng thấp)
        specificity = 100 / (1 + sum_cfd)
        return round(specificity, 2)

    def design_primers(self, sequence_template, target_start):
        """
        Bước 5: Thiết kế Primer (Primer3)
        Mục tiêu: Tạo sản phẩm PCR dài 200-300bp bao quanh vùng cắt.
        """
        # Cấu hình Primer3 để tránh vùng cắt (để primer không đè lên chỗ CRISPR cắt)
        # target_start là vị trí cắt. Ta cấm đặt primer trong vùng +/- 50bp quanh đó.
        try:
            res = primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': 'crispr_check',
                    'SEQUENCE_TEMPLATE': sequence_template,
                    'SEQUENCE_TARGET': [target_start - 5, 10]  # Vùng cần tránh (vùng cắt)
                },
                {
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[150, 300]],
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0
                }
            )
            return {
                "left": res.get('PRIMER_LEFT_0_SEQUENCE', 'N/A'),
                "right": res.get('PRIMER_RIGHT_0_SEQUENCE', 'N/A'),
                "product_size": res.get('PRIMER_PAIR_0_PRODUCT_SIZE', 0)
            }
        except:
            return None


# --- Main Function để gọi từ API ---
def run_crispor_analysis(full_sequence):
    engine = CrisporEngine(genome_index_path="data/R570_index")

    candidates = engine.find_candidates(full_sequence)
    results = []

    for cand in candidates:
        # 1. Tính Efficiency
        eff_score = engine.calculate_efficiency_score(cand['context_30bp'])

        # 2. Tìm & Tính Specificity
        off_targets = engine.search_off_targets(cand['guide_seq'])
        spec_score = engine.calculate_specificity_score(cand['guide_seq'], off_targets)

        # 3. Thiết kế Primer
        # Lưu ý: Primer design cần chuỗi template dài (full_sequence)
        primers = engine.design_primers(full_sequence, cand['start'])

        results.append({
            "sequence": cand['guide_seq'],
            "pam": cand['pam'],
            "location": f"{cand['start']}-{cand['end']}",
            "scores": {
                "efficiency_doench": eff_score,
                "specificity_cfd": spec_score
            },
            "primers": primers
        })

    # Sắp xếp: Ưu tiên Specificity trước (theo tư tưởng CRISPOR), rồi đến Efficiency
    results.sort(key=lambda x: (-x['scores']['specificity_cfd'], -x['scores']['efficiency_doench']))

    return results