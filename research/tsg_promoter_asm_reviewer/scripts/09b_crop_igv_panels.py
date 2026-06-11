#!/usr/bin/env python3
"""
Re-crop IGV PNG using actual PanelLayout dividerFractions from session XML:
  dividerFractions="0.0726, 0.2331, 0.3822, 0.5287, 0.6420, 0.7771, 0.7885"

Panel y-ranges (for 3776 px total height):
  A: 0      - 274    (DataPanel: HCC1395BL methyl phase VCFs)
  B: 274    - 880    (HCC1395 Tmode tagged ClairS pileup BAM TUMOR)
  C: 880    - 1443   (HCC1395BL ONT 5khz normal MOD BAM, grouped by HP)
  D: 1443   - 1996   (HCC1395 Tumor Dorado tagged BAM, grouped by PHASE)  <-- KEY
  E: 1996   - 2424   (HCC1395BL Dorado normal tagged BAM)
  F: 2424   - 2934   (HC SEQC2 BED + sSNV + sINDEL truth)
  G: 2934   - 2978   (junctions)
  H: 2978   - 3776   (Reference + ClairS/DeepSomatic VCFs + RefSeq + fp_cluster)
"""
import os
from PIL import Image, ImageDraw, ImageFont

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
os.makedirs(f"{BASE}/figures/igv_panels", exist_ok=True)

for view in ["full", "zoom"]:
    src = f"{BASE}/figures/brca2_igv_zar1l_session_{view}.png"
    img = Image.open(src)
    W, H = img.size
    print(f"\n{view}: {W}x{H}")

    fractions = [0.0726, 0.2331, 0.3822, 0.5287, 0.6420, 0.7771, 0.7885]
    boundaries = [0] + [int(f * H) for f in fractions] + [H]

    panels = [
        ("00_top_methyl_VCF_phase_calls", boundaries[0], boundaries[1]),
        ("01_HCC1395_Tumor_ClairS_BAM",   boundaries[1], boundaries[2]),
        ("02_HCC1395BL_Normal_HKU_HP",    boundaries[2], boundaries[3]),
        ("03_HCC1395_Tumor_Dorado_PHASE_KEY", boundaries[3], boundaries[4]),
        ("04_HCC1395BL_Normal_Dorado",    boundaries[4], boundaries[5]),
        ("05_SEQC2_truth_tracks",         boundaries[5], boundaries[6]),
        ("06_RefSeq_DeepSomatic",         boundaries[7], boundaries[8]),
    ]

    for name, y0, y1 in panels:
        if y1 - y0 < 30:
            continue
        crop = img.crop((0, y0, W, y1))
        out = f"{BASE}/figures/igv_panels/igv_{view}_{name}.png"
        crop.save(out)
        print(f"  -> {out} ({W}x{y1-y0})")

# Build composite for the KEY panels: Tumor Dorado PHASE-grouped (HP1/HP2/HP1-1)
print("\n[Build comparison composite: HCC1395BL Normal (top) vs HCC1395 Tumor Dorado (bottom)]")
normal_p = Image.open(f"{BASE}/figures/igv_panels/igv_zoom_02_HCC1395BL_Normal_HKU_HP.png")
tumor_p = Image.open(f"{BASE}/figures/igv_panels/igv_zoom_03_HCC1395_Tumor_Dorado_PHASE_KEY.png")
W = 1150
H = normal_p.height + tumor_p.height + 80  # 80 for divider label
composite = Image.new("RGB", (W, H), "white")
composite.paste(normal_p.convert("RGB"), (0, 0))
# Label divider
draw = ImageDraw.Draw(composite)
divider_y = normal_p.height + 10
draw.rectangle([0, divider_y, W, divider_y + 60], fill="#fef3c7")
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 16)
except:
    font = ImageFont.load_default()
draw.text((10, divider_y + 20), "↓ HCC1395 TUMOR Dorado tagged (grouped by PHASE = HP1 / HP2 / HP1-1) ↓",
          fill="#b45309", font=font)
composite.paste(tumor_p.convert("RGB"), (0, normal_p.height + 80))
composite.save(f"{BASE}/figures/igv_panels/composite_normal_vs_tumor_zoom.png")
print(f"  -> composite_normal_vs_tumor_zoom.png ({W}x{H})")

print("\nDone.")
