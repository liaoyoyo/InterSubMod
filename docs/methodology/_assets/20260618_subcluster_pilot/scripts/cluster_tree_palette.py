#!/usr/bin/env python3
"""[演化樹 tree-aware 配色 — 清楚定義 / 強對比版] phylo coarse_label → hex。取代扁平 CLUSTER_COL。

================================ 規則定義(SoT) ================================
標籤語意: coarse_label = 遞迴二元樹路徑。segs=label.split('-')。segs[0]=根群(1/2/3);
  之後每位 ∈{1,2}=該層分裂左/右子(1=左/2=右);深度 depth=len(segs)。例 '1-2-1'=根1→右子→左子(深3)。

配色三維(全 deterministic,同 label 跨 region 顏色固定):
  (1) 色相 HUE — **第一次分裂(segs[1]) 直接跳到「不同色相錨點」**(明顯對比, 非同色深淺):
        root1 左子=magenta / 右子=cyan ; root2 左=rose / 右=violet-magenta。→ 兩大支一眼分得開。
        之後更深的分裂(segs[2:]) = 在錨點鄰域內**微幅 hue 漂移**(左−/右+, 逐層減半)→ 表親可分、仍歸該大支。
  (2) 亮度 L — 隨**深度遞減**(深=暗);冗餘且 luminance 故色盲安全;兄弟 -2 微亮。
  (3) 飽和 S — 固定高飽和(鮮明)。
  例外(other/outlier/edge/'-') = 中性灰 #9ca3af。

護欄(硬規則, cluster 欄與 HP 欄同框並列):
  錨點色相只用 cluster-reserved 區 {magenta, rose, cyan},**脫離 HP 藍(0.58)/紫(0.78)/teal(0.47)、
  脫離 ALT 紅/REF 琥珀/T-N 橘綠**。唯一近鄰 = cyan ↔ HP3-teal(HP3=unphase somatic 罕見、不同欄,可接受+標註)。

⚠ 對比天花板(誠實): 扣掉 HP+保留軸後,真正「明顯不同」的色相只剩 ~3 區(magenta / rose / cyan)。
  兩大支用 magenta vs cyan = 最大對比;3 支以上靠 rose + 亮度。要再多需放寬護欄(讓 cluster 借用 HP 色)。
==============================================================================
"""
import colorsys

ROOT_BASE = {"1": 0.86, "2": 0.95, "3": 0.50}                 # 根群單獨出現時的色相
ANCHOR = {("1", "1"): 0.86, ("1", "2"): 0.50,                 # root1: 左 magenta / 右 cyan(最大對比)
          ("2", "1"): 0.95, ("2", "2"): 0.78,                 # root2(罕見): rose / violet-magenta
          ("3", "1"): 0.50, ("3", "2"): 0.86}
L_BY_DEPTH = {1: 0.62, 2: 0.55, 3: 0.46, 4: 0.38, 5: 0.31}
SAT = 0.68
DRIFT0 = 0.05                                                  # 深層微幅 hue 漂移起始(逐層減半)
GREY = "#9ca3af"


def _hex(h, l, s):
    r, g, b = colorsys.hls_to_rgb(h % 1.0, max(0.0, min(1.0, l)), s)
    return "#%02x%02x%02x" % (round(r * 255), round(g * 255), round(b * 255))


def cluster_tree_color(label, scheme="B"):
    if label in ("other", "outlier", "edge", "-", "", None):
        return GREY
    segs = str(label).split("-")
    root = segs[0]
    if len(segs) == 1:
        hue = ROOT_BASE.get(root, 0.86)
    else:
        hue = ANCHOR.get((root, segs[1]), ROOT_BASE.get(root, 0.86))   # 第一分裂=不同錨點(明顯)
        w = DRIFT0
        for s in segs[2:]:                                              # 深層=錨點鄰域微漂移(表親可分)
            w *= 0.5
            hue += w if s == "2" else -w
    depth = len(segs)
    L = L_BY_DEPTH.get(depth, 0.25)
    if depth > 1 and segs[-1] == "2":
        L += 0.04
    return _hex(hue, min(L, 0.80), SAT)


if __name__ == "__main__":
    labels = ["1", "1-1", "1-2", "1-1-1", "1-1-2", "1-2-1", "1-2-2", "1-1-1-1", "1-1-1-2", "2", "2-1", "other"]
    for L in labels:
        print(f"  {L:<11} {cluster_tree_color(L)}")
    print("\n第一分裂(應明顯不同): 1-1 vs 1-2 =", cluster_tree_color("1-1"), cluster_tree_color("1-2"), "(magenta vs cyan)")
    print("同支子分裂(亮度+微漂): 1-1-1 vs 1-1-2 =", cluster_tree_color("1-1-1"), cluster_tree_color("1-1-2"))
    print("表親(不同支): 1-1-1 vs 1-2-1 =", cluster_tree_color("1-1-1"), cluster_tree_color("1-2-1"))
    print("HP(不動): HP1=#60a5fa HP2=#c084fc HP3=#14b8a6 | ALT=#dc2626")
