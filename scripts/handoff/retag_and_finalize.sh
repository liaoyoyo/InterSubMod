#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-only
#
# 重新標記全基因組並整併成最終交付品。**可續跑**、**與 session 解耦**。
#
#   setsid nohup bash retag_and_finalize.sh > <LOG> 2>&1 &
#
# 為什麼要 setsid：這支要跑 ~5 小時，而背景任務會隨啟動它的 session 一起被 SIGTERM
# 掉（2026-08-13 實測：跑到 chr8 時 session 換掉，整個工作被 Terminated）。
# setsid 讓它自成 process group，不再收到 session 的終止訊號。
#
# 為什麼要可續跑：同樣是那次中斷。判斷「這條做完了沒」用 **.bam 與 .receipt 都在**，
# 不能只看 .bam —— 中斷點的那條會留下寫到一半的 BAM 但沒有 receipt，
# 只看 .bam 會把半殘檔當成完成品，那是最危險的失敗方式（安靜且下游看不出來）。

set -uo pipefail

LL=${LL:-/big7_disk/liaoyoyo2001/LongLineage}
STAGE=${STAGE:-/bip7_disk/liaoyoyo2001/lineage_out/HCC1395_v3}
FINAL=${FINAL:-/bip7_disk/liaoyoyo2001/HCC1395_final}
OLD=${OLD:-/bip7_disk/liaoyoyo2001/lineage_out/HCC1395_v1}
PREV=${PREV:-/bip7_disk/liaoyoyo2001/lineage_out/HCC1395_v2}
SAMPLE=${SAMPLE:-HCC1395}
EXPECT=${EXPECT:-40859727}

SRC=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam
SIDE=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/HCC1395.read_tags.tsv.gz
TOPO=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples/HCC1395/HCC1395.topology.jsonl

CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
mkdir -p "$STAGE"/{bam,logs}

echo "=== STAGE 1／3：逐染色體標記（可續跑）$(date -Is) ==="
for c in $CHROMS; do
    b="$STAGE/bam/${SAMPLE}.${c}.lineage_tagged.bam"
    r="$STAGE/bam/${SAMPLE}.${c}.tag_bam.receipt.json"
    if [[ -f "$b" && -f "$r" ]]; then
        echo "  $c SKIP（已完成）"; continue
    fi
    # 只有 .bam 沒 .receipt = 上次中斷在這條，那個 BAM 是半殘的，刪掉重做
    [[ -f "$b" && ! -f "$r" ]] && { echo "  $c 半殘，重做"; rm -f "$b" "$b.bai"; }
    t0=$SECONDS
    if "$LL/build_lineage_migrate/bin/longlineage-tag-bam" \
        --in-bam "$SRC" --sidecar "$SIDE" \
        --assignments "$OLD/assign/${SAMPLE}.all.read_lineage_assignments.tsv.gz" \
        --lineage-paths "$OLD/paths/${SAMPLE}.unit_lineage_paths.tsv.gz" \
        --topology "$TOPO" --region "$c" --threads 16 \
        --out-bam "$b" --receipt "$r" \
        > "$STAGE/logs/tag_bam.${c}.log" 2>&1; then
        echo "  $c OK ($((SECONDS-t0))s)"
    else
        echo "  $c FAIL —— 見 $STAGE/logs/tag_bam.${c}.log"; exit 1
    fi
done

echo "=== STAGE 2／3：合併與驗證 $(date -Is) ==="
parts=(); for c in $CHROMS; do parts+=("$STAGE/bam/${SAMPLE}.${c}.lineage_tagged.bam"); done
M="$STAGE/bam/${SAMPLE}.lineage_tagged.bam"
if [[ ! -f "$M" ]]; then
    samtools merge -@ 16 -f -o "$M" "${parts[@]}" 2> "$STAGE/logs/merge.log" \
        || { echo "merge FAIL"; exit 1; }
fi
[[ -f "$M.bai" ]] || samtools index -@ 16 "$M"

GOT=$(samtools idxstats "$M" | awk '{s+=$3+$4} END{print s}')
echo "  read 數 got=$GOT expect=$EXPECT"
[[ "$GOT" == "$EXPECT" ]] || { echo "FAIL：read 數不符，全部保留不動"; exit 1; }
samtools quickcheck -v "$M" || { echo "FAIL quickcheck"; exit 1; }
for t in lv lg ln ls; do
    n=$(samtools view "$M" chr21 2>/dev/null | head -50000 | grep -c "${t}:" || true)
    [[ "$n" -gt 0 ]] || { echo "FAIL：標籤 $t 沒寫出"; exit 1; }
    echo "  標籤 $t OK（chr21 前 50k 條命中 $n）"
done
echo "  ✓ 四道閘全過"

echo "=== STAGE 3／3：整併成最終交付品 $(date -Is) ==="
mkdir -p "$FINAL"/{bam,lineage}
mv "$M" "$M.bai" "$FINAL/bam/"
mv "$STAGE"/bam/*.tag_bam.receipt.json "$FINAL/bam/"
for d in assign paths bam_pre_lca_receipts; do
    [[ -d "$OLD/$d" ]] && cp -a "$OLD/$d" "$FINAL/lineage/"
done
for f in "${parts[@]}"; do rm -f "$f" "$f.bai"; done

echo "  最終 BAM：$(ls -la "$FINAL/bam/${SAMPLE}.lineage_tagged.bam" | awk '{printf "%.1f GiB", $5/2^30}')"
echo "  重建配方：$(du -sh "$FINAL/lineage" 2>/dev/null | cut -f1)"
echo
echo "DONE $(date -Is)"
echo "下一步（需人工確認後執行）："
echo "  1. 重跑 ISM → $FINAL/ism/"
echo "  2. 重建儀表板 → $FINAL/dashboard/"
echo "  3. build_final_delivery.py 產四個入口檔"
echo "  4. canonical symlink + verify_final_delivery.sh"
echo "  5. 確認無誤後刪 $PREV 與 $STAGE"
