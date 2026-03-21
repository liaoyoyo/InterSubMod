假設 ID: H004
假設內容: Quality_Score >= 50 AND CramersV > 0.15 組合 rescue
來源: Phase 2 結論：QS 有正向信號但門檻過高，嘗試加 CramersV 條件維持 precision
Pipeline Track: TO
Scale: S1 pilot
修改層級: Tier 1 (Python 模擬，使用既有 rescue_joined_features.tsv)
預期方向: 可能比 QS>=50 alone 差（CramersV≈0 in TO candidates）但 precision 更高
