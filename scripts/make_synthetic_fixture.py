#!/usr/bin/env python3
"""產生一組**完全合成**的最小輸入，讓 inter_sub_mod 能實際跑完一次。

這組 fixture 要解決的問題
-------------------------
公開 repo 沒有、也不應該有真實的 tumor BAM／reference／somatic VCF，所以任何人
（包括接手的新人與 AI）clone 下來之後都只能編譯與跑單元測試，永遠看不到
「輸入長什麼樣、輸出長什麼樣」。這一步缺口讓人無法確認自己的環境是否正確，
也無法在改動程式後判斷輸出格式有沒有被改壞。

**這組資料沒有任何生物學意義。** 它的唯一用途是證明「這條路徑跑得通」，
以及展示輸入輸出的形狀。不要拿它的數字做任何科學解讀。

為什麼參數要這樣選
------------------
inter_sub_mod 有幾個預設值會把不合格的 read 濾掉，合成資料若沒對齊這些門檻，
會「成功跑完但輸出全空」—— 那比報錯更難除錯，因為看起來像程式壞了：

    --min-read-length  預設 1000 bp   → 每條 read 長度設 1500
    --min-mapq         預設 20        → MAPQ 設 60
    --min-base-quality 預設 20        → base quality 設 40
    --min-common-coverage 預設 3      → 每條 read 至少覆蓋 3 個共同 CpG

參考序列刻意做成高 CpG 密度，讓甲基化矩陣不會稀疏到無法計算距離。

用法
----
    python3 scripts/make_synthetic_fixture.py                 # 產生到預設目錄
    python3 scripts/make_synthetic_fixture.py -o /tmp/myfix   # 指定輸出目錄

產生的檔案（皆不進 git，需要時重新產生即可）：

    synthetic.fa / .fai          參考序列
    synthetic.tumor.bam / .bai   帶 MM/ML 甲基化標籤與 HP 標籤的 tumor reads
    synthetic.somatic.vcf.gz     兩個合成 somatic SNV（附 tabix index）
    RUN.sh                       可直接執行的示範命令
"""

from __future__ import annotations

import argparse
import random
import subprocess
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    sys.exit("需要 pysam：pip install pysam")

CONTIG = "chr1"
CONTIG_LEN = 8000
READ_LEN = 1500          # 必須 > --min-read-length 預設 1000
N_READS = 24             # 一半帶 ALT、一半帶 REF，形成兩個可分辨的群
MAPQ = 60                # 必須 >= --min-mapq 預設 20
BASE_Q = 40              # 必須 >= --min-base-quality 預設 20
READ_START_BASE = 2600   # 每條 read 的起點基準，彼此錯開幾十 bp
READ_START_STAGGER = 6   # 錯開的檔數（實際起點 = BASE + (r % STAGGER) * 20）
SNV_POS = [3000, 3800]   # 0-based；兩個合成 somatic 位點，必須落在所有 read 的覆蓋內
SEED = 20260817          # 固定亂數種子 → 每次產生的 fixture 完全相同


def build_reference(rng: random.Random) -> str:
    """造一段高 CpG 密度的序列。

    每 20 bp 放一個 CG，其餘位置隨機。密度夠高，read 之間才會有足夠的
    common CpG 讓距離矩陣算得出來（--min-common-coverage 預設 3）。
    """
    bases = []
    for i in range(CONTIG_LEN):
        if i % 20 == 0 and i + 1 < CONTIG_LEN:
            bases.append("C")
        elif i % 20 == 1:
            bases.append("G")
        else:
            bases.append(rng.choice("ACGT"))
    return "".join(bases)


def encode_mm_ml(read_seq: str, rng: random.Random, methyl_prob: float):
    """把甲基化狀態編成 SAM 的 MM/ML 標籤。

    MM:Z:C+m?,<skip counts>;  —— 每個數字代表「跳過幾個未標記的 C 才到下一個標記的 C」
    ML:B:C,<0-255 機率>       —— 與 MM 逐一對應的機率值

    🔴 那個 `?` 是必要的，不是可有可無的裝飾。
    `src/core/MethylationParser.cpp` 比對的目標字串就是字面的 `"C+m?"`；
    若 MM 標籤寫成 `C+m.`（表示未列出的鹼基確定未修飾）或裸 `C+m`，
    InterSubMod 會找不到目標區塊而回傳空的 methylation calls ——
    **程式正常結束、不報錯，但 CpG 數量是 0**。這個 fixture 第一版就踩了這個坑。
    使用自己資料的人若發現 "Total CpG sites found: 0"，第一件事就是檢查
    `samtools view <bam> | head -1` 的 MM 標籤是不是 `C+m?` 這個 flavor。

    這裡對每個 CpG 的 C 依 methyl_prob 決定甲基化與否，機率值直接反映該決定，
    所以下游看到的 beta 值會集中在兩端，兩群 read 可分辨。
    """
    deltas, probs = [], []
    skipped = 0
    for i, base in enumerate(read_seq):
        if base != "C":
            continue
        is_cpg = i + 1 < len(read_seq) and read_seq[i + 1] == "G"
        if not is_cpg:
            skipped += 1
            continue
        if rng.random() < methyl_prob:
            deltas.append(skipped)
            probs.append(rng.randint(200, 255))   # 高機率甲基化
            skipped = 0
        else:
            deltas.append(skipped)
            probs.append(rng.randint(0, 40))      # 低機率甲基化
            skipped = 0
    if not deltas:
        return None, None
    return "C+m?," + ",".join(str(d) for d in deltas) + ";", probs


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out-dir",
                    default=str(Path(__file__).resolve().parents[1] / "tests/fixtures/synthetic"),
                    help="輸出目錄")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    # 幾何自檢：每個 SNV 都必須落在**所有** read 的覆蓋範圍內。
    # 不檢查的話，位點落在少數 read 之外只會讓那些 read 靜默地不帶 ALT，
    # 產出一份看起來正常、實際上不是預期設計的 fixture。
    last_start = READ_START_BASE + (READ_START_STAGGER - 1) * 20
    covered_lo, covered_hi = last_start, READ_START_BASE + READ_LEN
    for p in SNV_POS:
        if not (covered_lo <= p < covered_hi):
            sys.exit(f"設定錯誤：SNV 位置 {p} 不在所有 read 的共同覆蓋區間 "
                     f"[{covered_lo}, {covered_hi}) 內。請調整 SNV_POS 或 READ_LEN。")
    if last_start + READ_LEN > CONTIG_LEN:
        sys.exit(f"設定錯誤：最後一條 read 會超出 contig 長度 {CONTIG_LEN}。")

    # ---- 參考序列 ----
    ref_seq = build_reference(rng)
    fa = out / "synthetic.fa"
    with fa.open("w") as fh:
        fh.write(f">{CONTIG}\n")
        for i in range(0, len(ref_seq), 60):
            fh.write(ref_seq[i:i + 60] + "\n")
    subprocess.run(["samtools", "faidx", str(fa)], check=True)

    # ---- tumor BAM ----
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    sam_path = out / "synthetic.tumor.unsorted.bam"
    alt_base = {p: "T" if ref_seq[p] != "T" else "A" for p in SNV_POS}

    with pysam.AlignmentFile(str(sam_path), "wb", header=header) as bam:
        for r in range(N_READS):
            # 讓每條 read 都完整覆蓋兩個 SNV 位點，read 之間起點稍微錯開
            start = READ_START_BASE + (r % READ_START_STAGGER) * 20
            seq = list(ref_seq[start:start + READ_LEN])

            # 前半 read 帶 ALT（模擬帶 somatic 變異的那一群），後半維持 REF
            carries_alt = r < N_READS // 2
            if carries_alt:
                for p in SNV_POS:
                    seq[p - start] = alt_base[p]

            seq_str = "".join(seq)
            # 兩群給不同的甲基化傾向，讓輸出矩陣不是常數
            mm, ml = encode_mm_ml(seq_str, rng, 0.75 if carries_alt else 0.25)

            a = pysam.AlignedSegment()
            a.query_name = f"synth_read_{r:03d}"
            a.query_sequence = seq_str
            a.flag = 0 if r % 2 == 0 else 16
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = MAPQ
            a.cigartuples = [(0, READ_LEN)]          # 全 match，不引入 indel
            a.query_qualities = pysam.qualitystring_to_array(chr(BASE_Q + 33) * READ_LEN)
            tags = [("HP", 1 if carries_alt else 2, "i")]
            if mm:
                tags += [("MM", mm, "Z"), ("ML", ml)]
            a.set_tags(tags)
            bam.write(a)

    sorted_bam = out / "synthetic.tumor.bam"
    subprocess.run(["samtools", "sort", "-o", str(sorted_bam), str(sam_path)], check=True)
    subprocess.run(["samtools", "index", str(sorted_bam)], check=True)
    sam_path.unlink()

    # ---- somatic VCF ----
    vcf_plain = out / "synthetic.somatic.vcf"
    with vcf_plain.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={CONTIG},length={CONTIG_LEN}>\n")
        fh.write('##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,Description="Synthetic, no biological meaning">\n')
        fh.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for p in SNV_POS:
            fh.write(f"{CONTIG}\t{p + 1}\t.\t{ref_seq[p]}\t{alt_base[p]}\t60\tPASS\tSYNTHETIC\n")
    subprocess.run(["bgzip", "-f", str(vcf_plain)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_plain) + ".gz"], check=True)

    # ---- 示範執行腳本 ----
    run_sh = out / "RUN.sh"
    run_sh.write_text(f"""#!/usr/bin/env bash
# 用合成 fixture 跑一次 inter_sub_mod。
# 這組資料沒有生物學意義，只驗證「這條路徑跑得通」與展示輸入輸出的形狀。
set -euo pipefail
HERE="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
BIN="${{ISM_BIN:-$REPO/build/bin/inter_sub_mod}}"

[ -x "$BIN" ] || {{ echo "找不到 $BIN，請先建置：cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"; exit 1; }}

"$BIN" \\
    --tumor-bam "$HERE/synthetic.tumor.bam" \\
    --reference "$HERE/synthetic.fa" \\
    --vcf       "$HERE/synthetic.somatic.vcf.gz" \\
    --output-dir "${{1:-$HERE/output}}" \\
    --window-size 1000 \\
    --threads 2

echo
echo "輸出在 ${{1:-$HERE/output}}/ —— 每個 region 一個目錄，含 metadata.txt / reads.tsv /"
echo "methylation.csv / distance_matrix_*.csv，格式說明見 .claude/rules/output-structure.md"
""")
    run_sh.chmod(0o755)

    print(f"合成 fixture 已產生於 {out}")
    for f in sorted(out.iterdir()):
        print(f"  {f.stat().st_size:>10,} B  {f.name}")
    print(f"\n跑一次：bash {run_sh}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
