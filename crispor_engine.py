import re
import math
import subprocess
import primer3

# --- 1. DỮ LIỆU TRỌNG SỐ DOENCH (Từ file doenchScore.py bạn gửi) ---
# Format: (Vị trí, Nucleotide, Trọng số)
# Vị trí tính từ 1 (theo bài báo), nhưng Python index từ 0 nên ta sẽ xử lý trong hàm.
DOENCH_PARAMS = [
    (1, 'G', -0.2753771), (2, 'A', -0.3238875), (2, 'C', 0.17212887), (3, 'C', -0.1006662),
    (4, 'C', -0.2018029), (4, 'G', 0.24595663), (5, 'A', 0.03644004), (5, 'C', 0.09837684),
    (6, 'C', -0.7411813), (6, 'G', -0.3932644), (11, 'A', -0.466099), (14, 'A', 0.08537695),
    (14, 'C', -0.013814), (15, 'A', 0.27262051), (15, 'C', -0.1190226), (15, 'T', -0.2859442),
    (16, 'A', 0.09745459), (16, 'G', -0.1755462), (17, 'C', -0.3457955), (17, 'G', -0.6780964),
    (18, 'A', 0.22508903), (18, 'C', -0.5077941), (19, 'G', -0.4173736), (19, 'T', -0.054307),
    (20, 'G', 0.37989937), (20, 'T', -0.0907126), (21, 'C', 0.05782332), (21, 'T', -0.5305673),
    (22, 'T', -0.8770074), (23, 'C', -0.8762358), (23, 'G', 0.27891626), (23, 'T', -0.4031022),
    (24, 'A', -0.0773007), (24, 'C', 0.28793562), (24, 'T', -0.2216372), (27, 'G', -0.6890167),
    (27, 'T', 0.11787758), (28, 'C', -0.1604453), (29, 'G', 0.38634258), (1, 'GT', -0.6257787),
    (4, 'GC', 0.30004332), (5, 'AA', -0.8348362), (5, 'TA', 0.76062777), (6, 'GG', -0.4908167),
    (11, 'GG', -1.5169074), (11, 'TA', 0.7092612), (11, 'TC', 0.49629861), (11, 'TT', -0.5868739),
    (12, 'GG', -0.3345637), (13, 'GA', 0.76384993), (13, 'GC', -0.5370252), (16, 'TG', -0.7981461),
    (18, 'GG', -0.6668087), (18, 'TC', 0.35318325), (19, 'CC', 0.74807209), (19, 'TG', -0.3672668),
    (20, 'AC', 0.56820913), (20, 'CG', 0.32907207), (20, 'GA', -0.8364568), (20, 'GG', -0.7822076),
    (21, 'TC', -1.029693), (22, 'CG', 0.85619782), (22, 'CT', -0.4632077), (23, 'AA', -0.5794924),
    (23, 'AG', 0.64907554), (24, 'AG', -0.0773007), (24, 'CG', 0.28793562), (24, 'TG', -0.2216372),
    (26, 'GT', 0.11787758), (28, 'GG', -0.69774)
]


class CrisporEngine:
    def __init__(self, genome_index_path):
        self.genome_index = genome_index_path

    def find_candidates(self, sequence):
        """Bước 1: Tìm PAM (Giữ nguyên)"""
        candidates = []
        seq_upper = sequence.upper()
        # Tìm NGG
        for match in re.finditer(r'(?=([ATGC]GG))', seq_upper):
            pam_start = match.start()
            guide_start = pam_start - 20

            # Lấy ngữ cảnh 30bp (4bp trước + 20bp guide + 3bp PAM + 3bp sau)
            # Đây là format bắt buộc của thuật toán Doench
            context_start = guide_start - 4
            context_end = pam_start + 6

            if context_start < 0 or context_end > len(seq_upper):
                continue

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
        Bước 2: Tính điểm Doench 2014 (Logistic Regression) - CHÍNH CHỦ
        Dựa trên file doenchScore.py bạn cung cấp.
        """
        if len(context_30bp) != 30:
            return 0  # Nếu không đủ 30bp thì không tính được

        # Các hằng số từ file gốc
        intercept = 0.59763615
        gcHigh = -0.1665878
        gcLow = -0.2026259

        score = intercept
        seq = context_30bp

        # 1. Tính điểm GC
        guideSeq = seq[4:24]  # Lấy 20bp guide ở giữa
        gcCount = guideSeq.count("G") + guideSeq.count("C")

        gcWeight = 0
        if gcCount <= 10:
            gcWeight = gcLow
        if gcCount > 10:
            gcWeight = gcHigh

        score += abs(10 - gcCount) * gcWeight

        # 2. Tính điểm vị trí (Position Weights)
        # Duyệt qua bảng params và cộng điểm nếu khớp
        for pos, modelSeq, weight in DOENCH_PARAMS:
            # pos trong bảng là 1-based, Python là 0-based -> cần trừ 1
            # Tuy nhiên, file gốc đã tính toán offset dựa trên chuỗi 30bp
            # Kiểm tra logic file gốc:
            # (1,'G',...) -> seq[1] == 'G'. Trong chuỗi 30bp (4+20+3+3), vị trí 1 là nucleotide thứ 2 của vùng 4bp đầu.
            start_idx = pos - 1
            end_idx = start_idx + len(modelSeq)

            if start_idx < 0 or end_idx > 30: continue

            subSeq = seq[start_idx: end_idx]
            if subSeq == modelSeq:
                score += weight

        # 3. Chuyển đổi Sigmoid (1 / (1 + e^-score))
        try:
            final_score = 1.0 / (1.0 + math.exp(-score))
        except OverflowError:
            final_score = 0

        # Nhân 100 để ra thang điểm 0-100 cho đẹp
        return round(final_score * 100, 2)

    def search_off_targets(self, guide_seq):
        """
        [REAL] Gọi Bowtie2 qua WSL để tìm các vị trí khớp (cả đúng và gần đúng)
        """
        # Lưu ý: self.genome_index phải là đường dẫn tuyệt đối kiểu Linux
        # Ví dụ: /mnt/d/nhon-UWUET/2526I/project/sugarcane/data/R570/R570_index

        cmd = [
            "wsl",  # Cầu nối gọi từ Windows sang Linux
            "bowtie2",
            "-x", self.genome_index,
            "-c", guide_seq,  # Tìm chuỗi này
            "-k", "20",  # Tìm tối đa 20 vị trí (để tính CFD cho 20 off-target nguy hiểm nhất)
            "-N", "1",  # Cho phép sai 1 nucleotide ngay vùng Seed (để tìm được off-target)
            "-L", "20",  # Độ dài Seed
            "--no-unal",  # Không báo cáo nếu không tìm thấy gì
            "--no-hd"  # Không in Header (để dễ parse)
        ]

        try:
            # Gọi lệnh và bắt lấy kết quả in ra màn hình (stdout)
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )

            if process.returncode != 0:
                print(f"⚠️ Lỗi Bowtie2: {process.stderr}")
                # Fallback: Trả về list rỗng để không crash app
                return []

            # Nếu chạy ngon, đưa text kết quả cho hàm parse xử lý
            return self._parse_sam_output(process.stdout, guide_seq)

        except FileNotFoundError:
            print("❌ Lỗi: Không tìm thấy lệnh 'wsl'. Bạn có đang chạy trên Windows không?")
            return []
        except Exception as e:
            print(f"❌ Lỗi hệ thống: {str(e)}")
            return []

    def _parse_sam_output(self, sam_content, original_seq):
        """
        [HELPER] Đọc định dạng SAM của Bowtie2 và chuyển thành List Dict
        """
        off_targets = []

        # Duyệt qua từng dòng kết quả
        for line in sam_content.strip().split('\n'):
            if not line: continue

            parts = line.split('\t')
            if len(parts) < 10: continue

            # Cấu trúc SAM: [0]QNAME [1]FLAG [2]RNAME(Chrom) [3]POS ... [11+]TAGS
            chrom = parts[2]
            pos = parts[3]

            # Tìm các thẻ quan trọng trong phần Tags (từ cột 11 trở đi)
            # NM:i:x -> Số lượng mismatches (lỗi sai)
            # MD:Z:x -> Chuỗi mô tả vị trí sai (VD: 10A5C3)
            mismatches = 0
            md_str = ""

            for tag in parts[11:]:
                if tag.startswith("NM:i:"):
                    mismatches = int(tag.split(":")[2])
                elif tag.startswith("MD:Z:"):
                    md_str = tag.split(":")[2]

            off_targets.append({
                "seq": original_seq,  # Thực tế gRNA vẫn là nó
                "chrom": chrom,
                "position": pos,
                "mismatches": mismatches,
                "md": md_str  # Lưu cái này để sau này tính CFD chính xác
            })

        return off_targets

    def calculate_specificity_score(self, guide_seq, off_targets):
        """Bước 4: Tính CFD (Giữ nguyên)"""
        sum_cfd = 0
        for ot in off_targets:
            if ot['mismatches'] == 0: continue
            # Giả lập điểm phạt
            sum_cfd += 0.5
        return round(100 / (1 + sum_cfd), 2)

    def design_primers(self, sequence_template, target_start):
        """Bước 5: Primer3 (Giữ nguyên)"""
        try:
            res = primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': 'crispr',
                    'SEQUENCE_TEMPLATE': sequence_template,
                    'SEQUENCE_TARGET': [target_start - 5, 30]
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
            return {"left": "N/A", "right": "N/A", "error": "Primer Error"}


# --- Main Wrapper ---
def run_crispor_analysis(full_sequence, genome_index_path):  # <--- Thêm tham số này

    # Khởi tạo Engine với đường dẫn động được truyền vào
    engine = CrisporEngine(genome_index_path=genome_index_path)

    candidates = engine.find_candidates(full_sequence)
    results = []

    for cand in candidates:
        eff = engine.calculate_efficiency_score(cand['context_30bp'])

        # Bowtie2 sẽ dùng genome_index_path để tìm đúng bộ gen cần so sánh
        ot = engine.search_off_targets(cand['guide_seq'])

        spec = engine.calculate_specificity_score(cand['guide_seq'], ot)
        prim = engine.design_primers(full_sequence, cand['start'])

        results.append({
            "sequence": cand['guide_seq'],
            "pam": cand['pam'],
            "location": f"{cand['start']}-{cand['end']}",
            "scores": {
                "efficiency_doench": eff,
                "specificity_cfd": spec
            },
            "primers": prim
        })

    results.sort(key=lambda x: (
        -x['scores']['specificity_cfd'],
        -x['scores']['efficiency_doench']
    ))
    return results