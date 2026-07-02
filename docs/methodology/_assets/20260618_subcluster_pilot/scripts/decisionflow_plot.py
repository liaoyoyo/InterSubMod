#!/usr/bin/env python3
"""判別流程 4 統計圖 (讀 decisionflow_summary.json; tumor/merged 對照)。
數字全由 summary JSON 注入 (不手打)。輸出 figs_decisionflow/*.png。
fig1 5-態組成 / fig2 minority-size 分佈(門檻依據) / fig3 門檻敏感度 / fig4 切不出 mean-shift 分解。
註: 統計分佈圖非單位點甲基熱圖, 不適用 dual-panel 規則。"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"]=["Droid Sans Fallback","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False
A=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGS=f"{A}/figs_decisionflow"; os.makedirs(FIGS,exist_ok=True)
S=json.load(open(f"{A}/decisionflow_summary.json"))
SETS=[k for k in ("tumor","merged") if k in S]
ST=["S1","S2","S3","S4","S5"]
LAB={"S1":"①不可驗證","S2":"②1群/無訊號","S3":"③監督可分\n(mean-shift)","S4":"④可分未對齊\n(epiallele?雜訊?)","S5":"⑤確認真結構\n(對齊生物軸)"}
COL={"S1":"#9ca3af","S2":"#64748b","S3":"#d97706","S4":"#7c3aed","S5":"#059669"}
SETNAME={"tumor":"tumor-only (subclone 向)","merged":"merged tumor+normal (cis-control 對照)"}

# fig1: 5-態組成 (水平堆疊 %, tumor vs merged)
fig,ax=plt.subplots(figsize=(11,3.2))
for yi,nm in enumerate(SETS):
    st=S[nm]["states_T4"]; tot=S[nm]["n_total"]; left=0
    for s in ST:
        v=st[s]; w=100*v/tot
        ax.barh(yi,w,left=left,color=COL[s],edgecolor="white",height=0.62)
        if w>3.0: ax.text(left+w/2,yi,f"{s}\n{v}\n{w:.1f}%",ha="center",va="center",fontsize=8,color="white",fontweight="bold")
        left+=w
ax.set_yticks(range(len(SETS))); ax.set_yticklabels([SETNAME[n] for n in SETS],fontsize=9)
ax.set_xlim(0,100); ax.set_xlabel("位點比例 %（valid 群門檻 T=4）"); ax.set_title("HCC1395 全基因組判別流程 — 5 態組成（tumor vs merged 分開）",fontsize=11,fontweight="bold")
hand=[plt.Rectangle((0,0),1,1,color=COL[s]) for s in ST]
ax.legend(hand,[LAB[s].replace("\n"," ") for s in ST],loc="upper center",bbox_to_anchor=(0.5,-0.32),ncol=5,fontsize=7.5,frameon=False)
plt.tight_layout(); plt.savefig(f"{FIGS}/fig1_states.png",dpi=120,bbox_inches="tight"); plt.close()

# fig2: minority-size 分佈 (門檻依據) — aligned vs unaligned, tumor
fig,axes=plt.subplots(1,len(SETS),figsize=(5.6*len(SETS),3.4),squeeze=False)
for ci,nm in enumerate(SETS):
    ax=axes[0][ci]; md=S[nm]["minority_dist"]
    ks=[str(i) for i in range(1,12)]+["12+"]
    ha=md["hist_aligned"]; hu=md["hist_unaligned"]
    va=[ha.get(k,0) for k in ks]; vu=[hu.get(k,0) for k in ks]
    x=np.arange(len(ks))
    ax.bar(x,va,color="#059669",label=f"對齊(→⑤) n={md['aligned'].get('n',0)}")
    ax.bar(x,vu,bottom=va,color="#7c3aed",label=f"未對齊(→④) n={md['unaligned'].get('n',0)}")
    for T,c in ((3,"#dc2626"),(4,"#2563eb"),(5,"#9333ea")):
        ax.axvline(T-1-0.5,ls="--",lw=1.2,color=c,alpha=0.8)
        ax.text(T-1-0.5,ax.get_ylim()[1]*0.96,f"T={T}",color=c,fontsize=7.5,ha="center",va="top")
    ax.set_xticks(x); ax.set_xticklabels(ks,fontsize=7.5)
    ax.set_xlabel("最小 valid 群大小 (minority)"); ax.set_ylabel("位點數")
    ax.set_title(f"{nm}: 切群最小群分佈\n(中位 對齊{md['aligned'].get('median')} / 未對齊{md['unaligned'].get('median')})",fontsize=9.5)
    ax.legend(fontsize=7.5,frameon=False)
plt.suptitle("minority valid 群大小分佈 → valid≥4 門檻依據（T 線=門檻；左側=被當 outlier 的群）",fontsize=11,fontweight="bold")
plt.tight_layout(); plt.savefig(f"{FIGS}/fig2_minority.png",dpi=120,bbox_inches="tight"); plt.close()

# fig3: 門檻敏感度 T=3/4/5
fig,axes=plt.subplots(1,len(SETS),figsize=(5.2*len(SETS),3.4),squeeze=False)
for ci,nm in enumerate(SETS):
    ax=axes[0][ci]; ts=S[nm]["threshold_sensitivity"]
    Ts=sorted(int(t) for t in ts);
    s4=[ts[str(t)]["S4"] for t in Ts]; s5=[ts[str(t)]["S5"] for t in Ts]
    ar=[ts[str(t)]["split_align_rate"] for t in Ts]
    x=np.arange(len(Ts)); w=0.38
    ax.bar(x-w/2,s5,w,color="#059669",label="⑤ 對齊")
    ax.bar(x+w/2,s4,w,color="#7c3aed",label="④ 未對齊")
    ax.set_xticks(x); ax.set_xticklabels([f"T={t}" for t in Ts]); ax.set_ylabel("位點數")
    ax2=ax.twinx(); ax2.plot(x,ar,"o-",color="#dc2626",lw=1.6); ax2.set_ylabel("對齊率 %",color="#dc2626")
    ax2.set_ylim(0,100)
    for xi,a in zip(x,ar): ax2.text(xi,a+3,f"{a:.0f}%",color="#dc2626",fontsize=8,ha="center")
    ax.set_title(f"{nm}: 門檻敏感度",fontsize=9.5); ax.legend(fontsize=7.5,frameon=False,loc="upper right")
plt.suptitle("valid 群門檻 T 敏感度 — ④/⑤ 計數 + 對齊率（T↑→小群更易被當 outlier）",fontsize=11,fontweight="bold")
plt.tight_layout(); plt.savefig(f"{FIGS}/fig3_threshold.png",dpi=120,bbox_inches="tight"); plt.close()

# fig4: 切不出位點的 mean-shift 分解 (S2 / S3_loc / S3_disp)
fig,ax=plt.subplots(figsize=(9,3.0))
sub=["S2","S3_loc","S3_disp"]; subc={"S2":"#64748b","S3_loc":"#059669","S3_disp":"#d97706"}
subl={"S2":"②真無訊號\n(PERMANOVA 不顯著)","S3_loc":"③location-clean\n(真均值差,PERMDISP淨)","S3_disp":"③dispersion\n(散開度混淆)"}
for yi,nm in enumerate(SETS):
    st=S[nm]["states_T4"]; vals={"S2":st["S2"],"S3_loc":st["S3_loc"],"S3_disp":st["S3_disp"]}
    tot=sum(vals.values()); left=0
    for s in sub:
        v=vals[s]; w=100*v/tot if tot else 0
        ax.barh(yi,w,left=left,color=subc[s],edgecolor="white",height=0.6)
        if w>4: ax.text(left+w/2,yi,f"{v}\n{w:.1f}%",ha="center",va="center",fontsize=8,color="white",fontweight="bold")
        left+=w
ax.set_yticks(range(len(SETS))); ax.set_yticklabels([SETNAME[n] for n in SETS],fontsize=9)
ax.set_xlim(0,100); ax.set_xlabel("「切不出」位點中的比例 %")
ax.set_title("「切不出 ≠ 沒訊號」— 切不出位點的 a-priori mean-shift 分解（T=4）",fontsize=11,fontweight="bold")
hand=[plt.Rectangle((0,0),1,1,color=subc[s]) for s in sub]
ax.legend(hand,[subl[s].replace("\n"," ") for s in sub],loc="upper center",bbox_to_anchor=(0.5,-0.34),ncol=3,fontsize=7.5,frameon=False)
plt.tight_layout(); plt.savefig(f"{FIGS}/fig4_meanshift.png",dpi=120,bbox_inches="tight"); plt.close()
print("WROTE 4 figs to",FIGS)
for nm in SETS:
    s=S[nm]; print(f"{nm}: precond {s['pct_precond']}% | T4 {s['states_T4']} | split-align {s['split_align_rate']}% | cantsplit-meanshift {s['nonsplit_meanshift_rate']}%")
