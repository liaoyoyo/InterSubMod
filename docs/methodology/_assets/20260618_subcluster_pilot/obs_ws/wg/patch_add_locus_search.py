#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
patch_add_locus_search.py — add a locus-query box to the WG observation dashboard.

Surgically injects (no regeneration; the 11MB dashboard has no committed generator):
  1. CSS for the search box / help popover  (before </style>)
  2. a "位點查詢" input group as the FIRST control in the .ctrl bar
  3. a locus-query parser (chr:range / wildcard / regex) + one clause in matches()
  4. 'fq' into the input-listener array + clear button + help popover + out-block hide

Query grammar (matched against each case's `chr:pos`):
  chr1:100-1000   range (inclusive); 100- / -1000 one-sided; supports 18M/100k & , _
  chr1:1??        wildcard  ? = one digit, * = any digits (e.g. 18*)
  chr1*           chromosome wildcard (chr1, chr10-19)
  chr1:109267889  exact position ;  chr1  whole chromosome
  /chr2:18[0-9]/  or  re:^chr2:18   regex against "chr:pos"
  space/comma separated terms = OR

Usage: python3 patch_add_locus_search.py            # writes <dash>.patched
       python3 patch_add_locus_search.py --inplace  # backup + overwrite in place
"""
import sys, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, "20260618_HCC1395_observation_dashboard_FULL.standalone.html")

CSS = (
".cg input[type=text]{width:250px;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:11.5px}"
"#fq.bad{border-color:var(--no);color:#fca5a5}"
"#fqx{padding:4px 7px;font-size:11px}"
".qx{font-size:10.5px;color:var(--warn)}"
".qhelp{cursor:pointer;color:var(--ac);border:1px solid var(--ac);border-radius:50%;width:17px;height:17px;display:inline-flex;align-items:center;justify-content:center;font-size:11px;font-weight:700;user-select:none}"
"#qpop{display:none;position:absolute;top:34px;left:60px;z-index:60;background:#0b1222;border:1px solid var(--ac);border-radius:8px;padding:10px 13px;max-width:440px;font-size:11.5px;line-height:1.85;color:var(--fg);box-shadow:0 8px 28px #000b}"
"#qpop code{background:#1e293b;color:#fcd34d;padding:0 4px;border-radius:3px}"
"#qpop b{color:var(--ac)}"
)

CTRL_GROUP = (
' <div class="cg" style="position:relative"><label>位點查詢</label>'
'<input type="text" id="fq" placeholder="chr1:100-1000・chr1:1??・chr2:18M-19M・/chr1[0-9]:/" autocomplete="off" spellcheck="false">'
'<button id="fqx" title="清除查詢">✕</button>'
'<span class="qhelp" id="qhelp" title="查詢語法說明">?</span>'
'<span class="qx" id="fqn"></span><div id="qpop"></div></div>\n'
)

# parser block — raw string. Regex LITERALS use single backslash; strings that
# BECOME regexes use double backslash ('\\d'). Both must appear verbatim in output.
PARSER = r"""
/* ===== locus query (added by patch_add_locus_search.py) : chr:range / wildcard / regex ===== */
function _qnum(t){t=String(t).trim().replace(/[,_]/g,'');if(!t)return NaN;
 var m=t.match(/^(\d+(?:\.\d+)?)([kKmMgG]?)$/);if(!m)return NaN;
 var v=parseFloat(m[1]),u=m[2].toLowerCase();
 if(u==='k')v*=1e3;else if(u==='m')v*=1e6;else if(u==='g')v*=1e9;return Math.round(v);}
function _qesc(s){return s.replace(/[.*+?^${}()|[\]\\]/g,function(m){return '\\'+m;});}
function _qwild(p){var re='^',i,ch;for(i=0;i<p.length;i++){ch=p[i];
 if(ch==='?')re+='\\d';else if(ch==='*')re+='\\d*';else if(ch>='0'&&ch<='9')re+=ch;else re+=_qesc(ch);}
 return new RegExp(re+'$');}
function _qchrnorm(s){return String(s).replace(/^chr/i,'').toUpperCase();}
function _qchr(pat,chr){var c=_qchrnorm(chr);pat=_qchrnorm(pat);
 if(pat===''||pat==='*')return true;
 if(!/[*?]/.test(pat))return c===pat;
 var re='^',i,ch;for(i=0;i<pat.length;i++){ch=pat[i];
  if(ch==='?')re+='.';else if(ch==='*')re+='.*';else re+=_qesc(ch);}
 try{return new RegExp(re+'$').test(c);}catch(e){return false;}}
function _qpos(p){p=p.trim();if(p===''||p==='*')return {kind:'any'};
 var rm=p.match(/^(.*?)(?:\.\.|-|~)(.*)$/);
 if(rm&&!/[*?]/.test(p)&&(rm[1].trim()!==''||rm[2].trim()!=='')){
  var a=_qnum(rm[1]),b=_qnum(rm[2]);
  return {kind:'range',lo:isNaN(a)?-Infinity:a,hi:isNaN(b)?Infinity:b};}
 if(/[*?]/.test(p))return {kind:'wild',re:_qwild(p)};
 var n=_qnum(p);if(!isNaN(n))return {kind:'eq',v:n};
 return {kind:'bad'};}
function _qterm(t){
 var rm=t.match(/^\/(.*)\/([i]?)$/);
 if(t.slice(0,3).toLowerCase()==='re:')rm=[t,t.slice(3),''];
 if(rm){try{return {ok:true,type:'re',re:new RegExp(rm[1],rm[2]||'')};}catch(e){return {ok:false,type:'re'};}}
 var ix=t.indexOf(':');
 if(ix>=0)return {ok:true,type:'cp',chr:t.slice(0,ix),pos:_qpos(t.slice(ix+1))};
 if(/^(chr)?([0-9]{1,2}|[xXyYmM])[*?]*$/.test(t))return {ok:true,type:'chr',chr:t};
 return {ok:true,type:'sub',s:t.toLowerCase()};}
function parseQuery(q){q=(q||'').trim();if(!q)return null;
 var terms=q.split(/[\s,]+/).filter(Boolean).map(_qterm);
 return {ok:terms.every(function(t){return t.ok;}),terms:terms};}
function _qhit(t,c){var loc=(c.chr||'chr1')+':'+c.pos;
 if(t.type==='re')return t.re.test(loc);
 if(t.type==='sub')return loc.toLowerCase().indexOf(t.s)>=0||String(c.id||'').toLowerCase().indexOf(t.s)>=0;
 if(t.type==='chr')return _qchr(t.chr,c.chr||'chr1');
 if(t.type==='cp'){if(!_qchr(t.chr,c.chr||'chr1'))return false;var P=t.pos,pos=c.pos;
  if(P.kind==='any')return true;
  if(P.kind==='range')return pos>=P.lo&&pos<=P.hi;
  if(P.kind==='wild')return P.re.test(String(pos));
  if(P.kind==='eq')return pos===P.v;return false;}
 return false;}
var _QS=null,_QC=null;
function queryMatch(c){var raw=($('fq')?$('fq').value:'')||'';
 if(raw!==_QS){_QS=raw;_QC=parseQuery(raw);}
 if(!_QC)return true;
 var anyok=_QC.terms.some(function(t){return t.ok;});
 if(!anyok)return true;
 var i,t;for(i=0;i<_QC.terms.length;i++){t=_QC.terms[i];if(t.ok&&_qhit(t,c))return true;}
 return false;}
"""

ENH = r"""['fb','fs','ft','fl','fc','fg','cf','fv','sk','sd','tv','ts','fq'].forEach(id=>$(id).addEventListener('input',()=>{page=0;
 if(id==='tv')$('tvl').textContent=parseFloat($('tv').value).toFixed(2);
 if(id==='ts')$('tsl').textContent=parseFloat($('ts').value).toFixed(2);
 if(id==='fq')refreshQueryUI();draw();}));
function refreshQueryUI(){var raw=$('fq').value,q=parseQuery(raw);_QS=raw;_QC=q;
 var bad=!!(q&&q.terms.length&&q.terms.every(function(t){return !t.ok;}));
 $('fq').classList.toggle('bad',bad);
 $('fqn').textContent=raw?(bad?'語法?':''):'';}
$('fqx').addEventListener('click',function(){$('fq').value='';refreshQueryUI();page=0;draw();$('fq').focus();});
$('qhelp').addEventListener('click',function(e){e.stopPropagation();var p=$('qpop');p.style.display=(p.style.display==='block')?'none':'block';});
document.addEventListener('click',function(){$('qpop').style.display='none';});
$('qpop').addEventListener('click',function(e){e.stopPropagation();});
$('qpop').innerHTML='<b>位點查詢語法</b>（比對 <code>chr:pos</code>）<br>'+
 '<code>chr1:100-1000</code> 範圍（含端點）<br>'+
 '<code>chr2:18M-19M</code> 支援 k / M / g、千分位 , _ <br>'+
 '<code>chr1:100-</code> ・ <code>chr1:-1000</code> 單邊開放<br>'+
 '<code>chr1:1??</code> 萬用：<code>?</code>=一位數、<code>*</code>=任意位數（<code>18*</code>）<br>'+
 '<code>chr1:109267889</code> 精確　<code>chr1</code> 整條　<code>chr1*</code> chr1/chr10-19<br>'+
 '<code>/chr2:18[0-9]/</code> ・ <code>re:^chr2:18</code> 正規表示式<br>'+
 '多條件用空白或逗號分隔 ＝ <b>OR</b>（<code>chr2:18M-19M chr8</code>）';
var _draw0=draw;draw=function(){_draw0();if(($('fq').value||'').trim()){var og=$('outgrid');if(og)og.innerHTML='<div class="mut">位點查詢啟用中 — 被篩除區塊已隱藏（清除查詢以檢視近似落選項）</div>';}};
window.draw=draw;
refreshQueryUI();draw();"""


def patch(s):
    if 'id="fq"' in s:
        raise SystemExit("ABORT: dashboard already contains id=\"fq\" (already patched).")

    def repl(hay, old, new, label):
        n = hay.count(old)
        if n != 1:
            raise SystemExit(f"ABORT: anchor [{label}] found {n} times (expected 1).")
        return hay.replace(old, new)

    # 1. CSS before </style>
    s = repl(s, "</style>", CSS + "</style>", "css/</style>")
    # 2. control-bar search group as first .cg
    anchor_ctrl = '<div class="ctrl">\n <div class="cg"><label>桶</label>'
    s = repl(s, anchor_ctrl, '<div class="ctrl">\n' + CTRL_GROUP + ' <div class="cg"><label>桶</label>', "ctrl/first-group")
    # 3. parser block + matches() clause
    anchor_m = "function matches(c){\n"
    s = repl(s, anchor_m, PARSER + "\nfunction matches(c){\n if(!queryMatch(c))return false;\n", "matches/clause")
    # 4. listener array + enhancements (replace the original tail)
    old_tail = (
        "['fb','fs','ft','fl','fc','fg','cf','fv','sk','sd','tv','ts'].forEach(id=>$(id).addEventListener('input',()=>{page=0;\n"
        " if(id==='tv')$('tvl').textContent=parseFloat($('tv').value).toFixed(2);\n"
        " if(id==='ts')$('tsl').textContent=parseFloat($('ts').value).toFixed(2);draw();}));\ndraw();"
    )
    s = repl(s, old_tail, ENH, "listeners/tail")
    return s


def main():
    inplace = "--inplace" in sys.argv
    with open(F, encoding="utf-8") as fh:
        s = fh.read()
    out = patch(s)
    if inplace:
        bak = F + ".prepatch_bak"
        if not os.path.exists(bak):
            shutil.copy2(F, bak)
        with open(F, "w", encoding="utf-8") as fh:
            fh.write(out)
        print(f"[ok] patched IN PLACE: {F} ({len(out)//1024} KB)  backup: {bak}")
    else:
        op = F + ".patched"
        with open(op, "w", encoding="utf-8") as fh:
            fh.write(out)
        print(f"[ok] wrote {op} ({len(out)//1024} KB)  (use --inplace to overwrite)")


if __name__ == "__main__":
    main()
