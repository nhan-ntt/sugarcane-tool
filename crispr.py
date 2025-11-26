import re


def calculate_gc_content(sequence):
    """Tính tỷ lệ % G và C trong chuỗi"""
    if not sequence: return 0
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    return ((g_count + c_count) / len(sequence)) * 100


def find_crispr_targets(sequence: str, pam: str = "NGG", spacer_len: int = 20):
    """
    Tìm các vị trí có thể cắt bằng CRISPR-Cas9
    Input: Chuỗi DNA
    Output: Danh sách các gRNA candidates
    """
    sequence = sequence.upper()
    targets = []

    # 1. Xây dựng Regex cho PAM
    # NGG -> [ATGC]GG
    pam_regex = pam.replace("N", "[ATGC]")

    # Tìm tất cả vị trí xuất hiện PAM (lookahead để tìm chồng lấn)
    # (?=(...)) giúp tìm overlap matches
    for match in re.finditer(f"(?=({pam_regex}))", sequence):
        pam_start = match.start()
        pam_end = pam_start + len(pam)
        pam_seq = match.group(1)  # Lấy chuỗi PAM thực tế (vd: AGG, TGG)

        # 2. Lấy Spacer (20bp trước PAM)
        spacer_start = pam_start - spacer_len
        if spacer_start < 0:
            continue  # Bỏ qua nếu ở đầu chuỗi không đủ độ dài

        spacer_seq = sequence[spacer_start:pam_start]

        # 3. Tính toán chỉ số sinh học
        gc_content = calculate_gc_content(spacer_seq)

        # Lọc cơ bản (GC quá thấp hoặc quá cao thường kém hiệu quả)
        if 30 <= gc_content <= 80:
            targets.append({
                "position": f"{spacer_start}-{pam_end}",
                "spacer": spacer_seq,
                "pam": pam_seq,
                "gc_content": round(gc_content, 2),
                "score": "High" if 40 <= gc_content <= 60 else "Medium"
            })

    return targets