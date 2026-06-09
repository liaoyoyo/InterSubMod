#!/bin/bash
set -x
BT="$(cd "$(dirname "$0")" && pwd)"; export TMPDIR=/big7_disk/liaoyoyo2001/tmp; cd "$BT"
echo "### ASMS (Raineri 2024) repo 搜尋 ###"
gh search repos asms methylation --limit 10 2>&1 | head -20 || echo "gh-search-fail"
for u in EmanueleRaineri/asms eraineri/asms cnag-crg/asms; do echo "try $u:"; gh repo view "$u" 2>&1 | head -3; done
echo "### cvlr 再確認 ###"; gh repo view EmanueleRaineri/cvlr 2>&1 | head -3
echo "### DSS via conda (bioconductor) ###"
conda create -y -p "$BT/r_dss_env" -c bioconda -c conda-forge bioconductor-dss bioconductor-bsseq 2>&1 | tail -15
"$BT/r_dss_env/bin/Rscript" -e 'suppressMessages(library(DSS)); cat("DSS_OK\n")' 2>&1 | tail -3
echo "REST_INSTALL_DONE"
