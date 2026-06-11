#!/usr/bin/env python3
"""
Crop IGV PNG (1150 x 3776) into 4 focused sub-panels with annotations:
  A — Top: 軌道 header + phased methyl VCF + HCC1395 ClairS small BAM (~0-580)
  B — Mid-upper: HCC1395BL ONT 5khz tagged BAM (normal, grouped by HP) (~580-1380)
  C — Mid: HCC1395 Tumor Dorado tagged BAM (KEY — grouped by PHASE = HP1/HP2/HP1-1/HP2-1) (~1380-2680)
  D — Bottom: HCC1395BL Dorado normal + SEQC2 truth + reference + RefSeq (~2680-3776)
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

    # Approximate panel boundaries (visual inspection)
    # Adjusted to capture meaningful slices
    panels = [
        ("A_top_VCF_ClairS_BAM", 0, 600),
        ("B_HCC1395BL_5khz_normal", 600, 1400),
        ("C_HCC1395_Tumor_Dorado_KEY", 1400, 2700),
        ("D_HCC1395BL_Dorado_truth_ref", 2700, H),
    ]

    for name, y0, y1 in panels:
        crop = img.crop((0, y0, W, y1))
        out = f"{BASE}/figures/igv_panels/igv_{view}_{name}.png"
        crop.save(out)
        print(f"  -> {out} ({W}x{y1-y0})")

# Build composite: side-by-side full + zoom for one key panel (C = tumor)
print("\n[Build composite figure]")
img_full = Image.open(f"{BASE}/figures/igv_panels/igv_full_C_HCC1395_Tumor_Dorado_KEY.png")
img_zoom = Image.open(f"{BASE}/figures/igv_panels/igv_zoom_C_HCC1395_Tumor_Dorado_KEY.png")
W = 1150
H = max(img_full.height, img_zoom.height)
composite = Image.new("RGB", (W*2 + 20, H), "white")
composite.paste(img_full.convert("RGB"), (0, 0))
composite.paste(img_zoom.convert("RGB"), (W + 20, 0))
composite.save(f"{BASE}/figures/igv_panels/composite_tumor_full_vs_zoom.png")
print(f"  -> composite_tumor_full_vs_zoom.png ({W*2+20}x{H})")

print("\nAll panel crops written.")
