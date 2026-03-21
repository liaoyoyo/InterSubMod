#!/usr/bin/env python3
"""
將 InterSubMod 研究週報轉成可重複生成的 PPTX。

使用方式:
python scripts/analysis/build_weekly_report_pptx.py \
  --config docs/references/manual/assets/20260311_weekly_report_pptx_config_01.json
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_AUTO_SHAPE_TYPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt


SLIDE_W = 13.333
SLIDE_H = 7.5


def hex_rgb(value: str) -> RGBColor:
    value = value.replace("#", "").strip()
    return RGBColor(int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16))


def load_config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def set_bg(slide, color: str) -> None:
    fill = slide.background.fill
    fill.solid()
    fill.fore_color.rgb = hex_rgb(color)


def add_text(
    slide,
    text: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    font_size: float = 18,
    color: str = "1B1F24",
    bold: bool = False,
    font_name: str = "Calibri",
    align: str = "left",
    valign: MSO_ANCHOR = MSO_ANCHOR.TOP,
    italic: bool = False,
    margin: float = 0.05,
    line_spacing: float | None = None,
) -> None:
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    box.text_frame.clear()
    box.text_frame.vertical_anchor = valign
    box.text_frame.word_wrap = True
    box.text_frame.margin_left = Inches(margin)
    box.text_frame.margin_right = Inches(margin)
    box.text_frame.margin_top = Inches(margin)
    box.text_frame.margin_bottom = Inches(margin)
    p = box.text_frame.paragraphs[0]
    p.alignment = {"left": PP_ALIGN.LEFT, "center": PP_ALIGN.CENTER, "right": PP_ALIGN.RIGHT}[align]
    if line_spacing is not None:
        p.line_spacing = line_spacing
    run = p.add_run()
    run.text = text
    run.font.size = Pt(font_size)
    run.font.name = font_name
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = hex_rgb(color)


def add_bullets(
    slide,
    bullets: Iterable[str],
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    font_size: float = 18,
    color: str = "1B1F24",
    bullet_color: str = "C76B50",
    font_name: str = "Calibri",
    level: int = 0,
) -> None:
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.clear()
    tf.word_wrap = True
    tf.margin_left = Inches(0.02)
    tf.margin_right = Inches(0.02)
    tf.margin_top = Inches(0.02)
    tf.margin_bottom = Inches(0.02)
    for idx, item in enumerate(bullets):
        p = tf.paragraphs[0] if idx == 0 else tf.add_paragraph()
        p.level = level
        p.alignment = PP_ALIGN.LEFT
        p.bullet = True
        p.space_after = Pt(6)
        run = p.add_run()
        run.text = item
        run.font.size = Pt(font_size)
        run.font.name = font_name
        run.font.color.rgb = hex_rgb(color)
        if hasattr(p, "_pPr") and p._pPr is not None:
            pass


def add_rect(slide, x: float, y: float, w: float, h: float, *, fill: str, line: str | None = None, radius: bool = False, transparency: int = 0):
    shape_type = MSO_AUTO_SHAPE_TYPE.ROUNDED_RECTANGLE if radius else MSO_AUTO_SHAPE_TYPE.RECTANGLE
    shape = slide.shapes.add_shape(shape_type, Inches(x), Inches(y), Inches(w), Inches(h))
    shape.fill.solid()
    shape.fill.fore_color.rgb = hex_rgb(fill)
    shape.fill.transparency = transparency
    if line is None:
        shape.line.fill.background()
    else:
        shape.line.color.rgb = hex_rgb(line)
        shape.line.width = Pt(1)
    return shape


def wrap_path_for_slide(path: str, max_len: int = 46, max_lines: int = 4) -> str:
    parts = [token for token in re.split(r"([/ ])", path) if token]
    lines = []
    current = ""
    for part in parts:
        segment = part
        if not current:
            current = segment.lstrip()
            continue
        if len(current) + len(segment) <= max_len:
            current += segment
        else:
            lines.append(current.rstrip())
            current = segment.lstrip()
    if current:
        lines.append(current.rstrip())
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        if len(lines[-1]) >= max_len - 1:
            lines[-1] = lines[-1][: max_len - 1] + "…"
        else:
            lines[-1] += " …"
    return "\n".join(lines)


def add_title(slide, title: str, subtitle: str, theme: dict, *, dark: bool = False) -> None:
    fg = theme["light_text"] if dark else theme["text"]
    muted = theme["light_muted"] if dark else theme["muted"]
    add_text(slide, title, 0.8, 0.55, 11.7, 0.8, font_size=28, color=fg, bold=True, font_name=theme["font_title"], margin=0)
    add_text(slide, subtitle, 0.8, 1.3, 11.3, 0.5, font_size=14, color=muted, font_name=theme["font_body"], margin=0)


def add_footer(slide, theme: dict, page_num: int, total_pages: int, footer_text: str, *, dark: bool = False) -> None:
    color = theme["light_muted"] if dark else "4B5563"
    add_text(slide, footer_text, 0.8, 6.78, 9.0, 0.3, font_size=9.5, color=color, margin=0)
    add_text(slide, f"{page_num}/{total_pages}", 12.0, 6.75, 0.55, 0.3, font_size=11, color=color, align="right", margin=0)


def style_table(table, theme: dict, *, header_fill: str, header_font: str, body_fill: str = "FFFFFF", zebra_fill: str | None = None, font_size: int = 12) -> None:
    rows = len(table.rows)
    cols = len(table.columns)
    for r in range(rows):
        for c in range(cols):
            cell = table.cell(r, c)
            cell.fill.solid()
            if r == 0:
                cell.fill.fore_color.rgb = hex_rgb(header_fill)
            elif zebra_fill and r % 2 == 0:
                cell.fill.fore_color.rgb = hex_rgb(zebra_fill)
            else:
                cell.fill.fore_color.rgb = hex_rgb(body_fill)
            for paragraph in cell.text_frame.paragraphs:
                paragraph.alignment = PP_ALIGN.LEFT if c == 0 else PP_ALIGN.CENTER
                for run in paragraph.runs:
                    run.font.name = theme["font_body"]
                    run.font.size = Pt(font_size)
                    if r == 0:
                        run.font.bold = True
                        run.font.color.rgb = hex_rgb(header_font)
                    else:
                        run.font.color.rgb = hex_rgb(theme["text"])


def add_table(slide, rows: list[list[str]], x: float, y: float, w: float, h: float, theme: dict, *, header_fill: str, body_fill: str = "FFFFFF", zebra_fill: str | None = None, font_size: int = 12, col_widths: list[float] | None = None):
    table = slide.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y), Inches(w), Inches(h)).table
    if col_widths:
        total = sum(col_widths)
        for idx, width in enumerate(col_widths):
            table.columns[idx].width = Inches(w * (width / total))
    else:
        col_width = Inches(w / len(rows[0]))
        for idx in range(len(rows[0])):
            table.columns[idx].width = col_width
    row_h = Inches(h / len(rows))
    for idx in range(len(rows)):
        table.rows[idx].height = row_h
    for r, row in enumerate(rows):
        for c, value in enumerate(row):
            table.cell(r, c).text = str(value)
    style_table(table, theme, header_fill=header_fill, header_font=theme["light_text"], body_fill=body_fill, zebra_fill=zebra_fill, font_size=font_size)
    return table


def create_grouped_bar_chart(output_path: Path, title: str, categories: list[str], series: list[dict], theme: dict, *, ylim: tuple[float, float], ylabel: str) -> None:
    plt.figure(figsize=(10, 4.6), dpi=200)
    x = range(len(categories))
    width = 0.22 if len(series) >= 3 else 0.3
    offsets = [(-width), 0, width, 2 * width]
    for idx, item in enumerate(series):
        values = item["values"]
        pos = [v + offsets[idx] for v in x]
        plt.bar(pos, values, width=width * 0.9, color=f"#{item['color']}", label=item["name"])
        for px, val in zip(pos, values):
            plt.text(px, val + 0.001, f"{val:.4f}", ha="center", va="bottom", fontsize=9)
    plt.xticks(list(x), categories, fontsize=10)
    plt.yticks(fontsize=10)
    plt.ylim(*ylim)
    plt.ylabel(ylabel, fontsize=11)
    plt.title(title, fontsize=13, weight="bold")
    plt.grid(axis="y", linestyle="--", alpha=0.25)
    plt.legend(frameon=False, fontsize=9, ncol=len(series), loc="upper center", bbox_to_anchor=(0.5, 1.15))
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close()


def create_horizontal_bar_chart(output_path: Path, panel_specs: list[dict], theme: dict) -> None:
    fig, axes = plt.subplots(1, len(panel_specs), figsize=(11.2, 4.4), dpi=200, sharex=False)
    if len(panel_specs) == 1:
        axes = [axes]
    for ax, spec in zip(axes, panel_specs):
        labels = [row["name"] for row in spec["rows"]]
        values = [row["f1"] for row in spec["rows"]]
        colors = [f"#{row['color']}" for row in spec["rows"]]
        ypos = list(range(len(labels)))[::-1]
        ax.barh(ypos, values[::-1], color=colors[::-1], height=0.62)
        ax.set_yticks(ypos)
        ax.set_yticklabels(labels[::-1], fontsize=9)
        ax.set_xlim(spec["xlim"][0], spec["xlim"][1])
        ax.set_title(spec["title"], fontsize=11, weight="bold")
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        for y, val in zip(ypos, values[::-1]):
            ax.text(val + 0.001, y, f"{val:.4f}", va="center", fontsize=9)
    fig.suptitle("Paired Raw Pileup/Full vs Final Outputs", fontsize=13, weight="bold")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close()


def create_timeline_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["bg"])
    add_title(slide, "研究推進時間軸", "依日期整理本週主線進展與收斂節點", theme)
    y = 2.2
    x_positions = [0.85, 3.95, 7.05, 10.15]
    add_rect(slide, 1.25, 4.95, 10.35, 0.08, fill=theme["accent"], line=None, radius=False)
    for x, item in zip(x_positions, cfg["deck"]["timeline"]):
        add_rect(slide, x, y, 2.2, 2.55, fill="FFFFFF", line=theme["line"], radius=True)
        add_rect(slide, x + 0.12, y + 0.12, 0.6, 0.35, fill=theme["dark"], line=None, radius=True)
        add_text(slide, item["date"], x + 0.18, y + 0.15, 0.9, 0.2, font_size=11, color=theme["light_text"], bold=True, margin=0)
        add_text(slide, item["title"], x + 0.14, y + 0.56, 1.92, 0.58, font_size=13, color=theme["text"], bold=True, margin=0)
        add_text(slide, item["result"], x + 0.14, y + 1.18, 1.9, 1.12, font_size=11.0, color="4F5B67", margin=0)
    add_rect(slide, 0.8, 6.05, 12.0, 0.45, fill="FFFFFF", line=theme["line"], radius=True)
    add_text(slide, "時間軸重點：3/05 先排除 mixed 主線，3/07 跑通 paired / TO，3/08~09 確立 candidate-specific rescue，3/11 完成 phase 2 與 annotation 回接。", 1.0, 6.2, 11.4, 0.2, font_size=10.3, color="4F5B67", margin=0)


def create_presentation(cfg: dict) -> Presentation:
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)
    prs.core_properties.author = cfg["meta"]["author"]
    prs.core_properties.title = cfg["meta"]["deck_title"]
    prs.core_properties.subject = cfg["meta"]["deck_subtitle"]
    return prs


def add_title_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["dark"])
    add_rect(slide, 0, 0, 13.333, 7.5, fill=theme["dark"], line=None)
    add_rect(slide, 8.9, 0, 4.433, 7.5, fill=theme["accent"], line=None, transparency=12)
    add_rect(slide, 0.7, 0.65, 1.45, 0.35, fill=theme["accent2"], line=None, radius=True)
    add_text(slide, cfg["meta"]["kicker"], 0.83, 0.72, 2.5, 0.2, font_size=11, color=theme["light_text"], bold=True, margin=0)
    add_text(slide, cfg["meta"]["deck_title"], 0.8, 1.35, 8.2, 1.25, font_size=28, color=theme["light_text"], bold=True, font_name=theme["font_title"], margin=0, line_spacing=1.1)
    add_text(slide, cfg["meta"]["deck_subtitle"], 0.84, 2.78, 7.3, 0.55, font_size=15, color=theme["light_muted"], font_name=theme["font_body"], margin=0)
    add_text(slide, cfg["meta"]["summary_line"], 0.84, 3.55, 7.15, 0.95, font_size=18, color=theme["light_text"], bold=True, font_name=theme["font_body"], margin=0)
    add_text(slide, "來源週報：20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md", 0.84, 5.72, 7.0, 0.35, font_size=10.5, color=theme["light_muted"], font_name="Consolas", margin=0)
    add_rect(slide, 9.15, 1.55, 3.35, 3.55, fill="F6EFE6", line=None, radius=True)
    add_text(slide, "本頁導讀", 9.45, 1.84, 1.2, 0.2, font_size=14, color=theme["text"], bold=True, margin=0)
    add_bullets(
        slide,
        [
            "四象限資料：5kHz / DORADO × paired / TO",
            "phase 2：paired raw pileup/full 對照已完成",
            "annotation layer 已正式接回分析輸出層",
        ],
        9.35,
        2.22,
        2.7,
        1.7,
        font_size=11.5,
        color=theme["text"],
    )
    add_text(
        slide,
        f"報告整理時間：{cfg['meta']['generated_date']}\n輸出模式：{cfg['meta']['audience_mode']}\n報告人：{cfg['meta']['author']}",
        9.42,
        4.18,
        2.8,
        0.86,
        font_size=10.8,
        color=theme["muted"],
        margin=0,
    )


def add_key_findings_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["bg"])
    add_title(slide, "本週重點結論", "先看四條證據鏈，快速掌握可直接沿用的研究口徑", theme)
    cards = cfg["deck"]["key_findings"]
    positions = [(0.8, 1.95), (6.95, 1.95), (0.8, 4.32), (6.95, 4.32)]
    for (x, y), card in zip(positions, cards):
        add_rect(slide, x, y, 5.55, 2.08, fill="FFFFFF", line=theme["line"], radius=True)
        add_rect(slide, x + 0.18, y + 0.2, 0.48, 0.48, fill=card["accent"], line=None, radius=True)
        add_text(slide, card["title"], x + 0.8, y + 0.18, 4.4, 0.35, font_size=16, color=theme["text"], bold=True, margin=0)
        add_text(slide, card["body"], x + 0.18, y + 0.76, 5.05, 1.02, font_size=12.0, color="4F5B67", margin=0)


def add_framework_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["bg_alt"])
    add_title(slide, "背景、問題與四象限", "這一週的所有結果都用同一套三層框架與四象限資料集解讀", theme)
    add_rect(slide, 0.8, 1.9, 3.7, 4.82, fill="FFFFFF", line=theme["line"], radius=True)
    add_text(slide, "三層分工", 1.02, 2.12, 1.8, 0.28, font_size=16, color=theme["text"], bold=True, margin=0)
    layer_y = 2.6
    for idx, layer in enumerate(cfg["deck"]["framework"]["layers"]):
        add_rect(slide, 1.0, layer_y + idx * 1.18, 3.25, 0.9, fill=layer["fill"], line=None, radius=True)
        add_text(slide, layer["title"], 1.18, layer_y + 0.08 + idx * 1.18, 2.8, 0.2, font_size=12.8, color=layer["title_color"], bold=True, margin=0)
        add_text(slide, layer["body"], 1.18, layer_y + 0.34 + idx * 1.18, 2.82, 0.38, font_size=10.0, color=layer["body_color"], margin=0)
    add_text(slide, cfg["deck"]["framework"]["problem_statement"], 1.0, 6.2, 3.15, 0.32, font_size=10.2, color="4F5B67", margin=0)

    add_rect(slide, 4.8, 1.9, 7.7, 4.82, fill="FFFFFF", line=theme["line"], radius=True)
    add_text(slide, "資料四象限與可比性", 5.05, 2.12, 2.8, 0.28, font_size=16, color=theme["text"], bold=True, margin=0)
    q_positions = [(5.05, 2.62), (8.55, 2.62), (5.05, 4.48), (8.55, 4.48)]
    for (x, y), quad in zip(q_positions, cfg["deck"]["framework"]["quadrants"]):
        add_rect(slide, x, y, 3.3, 1.48, fill=quad["fill"], line=None, radius=True)
        add_text(slide, quad["title"], x + 0.16, y + 0.12, 2.8, 0.24, font_size=13, color=quad["title_color"], bold=True, margin=0)
        add_text(slide, quad["body"], x + 0.16, y + 0.42, 2.95, 0.72, font_size=10.5, color=quad["body_color"], margin=0)
    add_text(slide, cfg["deck"]["framework"]["scope_note"], 5.05, 6.02, 6.95, 0.5, font_size=10.1, color="4F5B67", margin=0)


def add_goal_status_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["bg"])
    add_title(slide, "主要目標與完成度", "這一週的主線任務已從實驗探索推進到制度化與 annotation 層落地", theme)
    rows = [["主目標", "子目標", "完成度", "本週結論"]]
    for item in cfg["deck"]["goals"]:
        rows.append([item["goal"], item["subgoal"], item["completion"], item["conclusion"]])
    add_table(slide, rows, 0.65, 1.82, 12.0, 3.95, theme, header_fill=theme["dark"], body_fill="FFFFFF", zebra_fill="F4F0E8", font_size=9.6, col_widths=[2.1, 2.45, 0.9, 4.65])
    y = 5.95
    for idx, item in enumerate(cfg["deck"]["goal_summary_cards"]):
        x = 0.8 + idx * 3.95
        add_rect(slide, x, y, 3.55, 0.7, fill=item["fill"], line=None, radius=True)
        add_text(slide, item["title"], x + 0.15, y + 0.08, 3.0, 0.18, font_size=11.5, color=item["title_color"], bold=True, margin=0)
        add_text(slide, item["body"], x + 0.15, y + 0.29, 3.0, 0.2, font_size=9.9, color=item["body_color"], margin=0)


def add_paired_rerun_slide(slide, cfg: dict, theme: dict, chart_path: Path) -> None:
    set_bg(slide, theme["bg_alt"])
    add_title(slide, "paired pure 正式 rerun", "5kHz 仍是最有壓力的主樣本；DORADO 則是穩定交叉驗證基準", theme)
    add_rect(slide, 0.8, 1.8, 7.35, 4.95, fill="FFFFFF", line=theme["line"], radius=True)
    slide.shapes.add_picture(str(chart_path), Inches(1.0), Inches(2.05), width=Inches(6.95), height=Inches(3.75))
    add_text(slide, "圖意：X 軸是 paired dataset，Y 軸是 F1。每組三根柱分別代表 ClairS、LongPhase-S、InterSubMod。", 1.0, 5.88, 6.7, 0.32, font_size=10.2, color="4F5B67", margin=0)
    add_rect(slide, 8.35, 1.8, 4.15, 2.1, fill="FFFFFF", line=theme["line"], radius=True)
    add_text(slide, "5kHz paired", 8.3, 2.05, 2.0, 0.22, font_size=14, color=theme["text"], bold=True, margin=0)
    add_text(slide, "ClairS 0.8443 → LongPhase-S 0.8522 → InterSubMod 0.8532", 8.3, 2.35, 3.8, 0.34, font_size=11.5, color=theme["muted"], margin=0)
    add_text(slide, "意義：在最 noisy 主樣本上，甲基與 phase 整合有實際正增益。", 8.3, 2.86, 3.75, 0.34, font_size=11.5, color=theme["text"], margin=0)
    add_rect(slide, 8.35, 4.08, 4.15, 2.1, fill="FFFFFF", line=theme["line"], radius=True)
    add_text(slide, "DORADO paired", 8.3, 4.27, 2.0, 0.22, font_size=14, color=theme["text"], bold=True, margin=0)
    add_text(slide, "ClairS 0.8565 → LongPhase-S 0.8592 → InterSubMod 0.8590", 8.3, 4.57, 3.85, 0.34, font_size=11.5, color=theme["muted"], margin=0)
    add_text(slide, "意義：跨平台方向可驗證，但同一規則不保證可攜。", 8.3, 5.08, 3.75, 0.34, font_size=11.5, color=theme["text"], margin=0)


def add_to_rescue_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["bg"])
    add_title(slide, "tumor-only candidate-specific rescue", "TO 下已正式證明：甲基訊號不是全負效果，但仍未超過最佳 caller-only gate", theme)
    left = cfg["deck"]["to_rescue"]["hcc1395_5khz_to"]
    right = cfg["deck"]["to_rescue"]["hcc1395_dorado_to"]
    for x, data in [(0.8, left), (6.85, right)]:
        add_rect(slide, x, 1.9, 5.55, 4.65, fill="FFFFFF", line=theme["line"], radius=True)
        add_text(slide, data["label"], x + 0.2, 2.12, 2.8, 0.25, font_size=15, color=theme["text"], bold=True, margin=0)
        rows = [
            ["規則", "TP", "FP", "ΔF1"],
            [data["caller_only"]["rule"], data["caller_only"]["tp"], data["caller_only"]["fp"], data["caller_only"]["delta_f1"]],
            [data["methyl_support"]["rule"], data["methyl_support"]["tp"], data["methyl_support"]["fp"], data["methyl_support"]["delta_f1"]],
            [data["label_support"]["rule"], data["label_support"]["tp"], data["label_support"]["fp"], data["label_support"]["delta_f1"]],
        ]
        add_table(slide, rows, x + 0.18, 2.55, 5.08, 2.0, theme, header_fill=theme["dark"], body_fill="FFFFFF", zebra_fill="F4F0E8", font_size=9.1, col_widths=[3.35, 0.56, 0.56, 0.95])
        add_text(slide, data["interpretation"], x + 0.2, 4.72, 4.92, 0.74, font_size=10.8, color="4F5B67", margin=0)
        add_rect(slide, x + 0.2, 5.72, 4.98, 0.44, fill=data["note_fill"], line=None, radius=True)
        add_text(slide, data["note"], x + 0.34, 5.84, 4.55, 0.18, font_size=10.2, color=data["note_color"], bold=True, margin=0)


def add_phase2_slide(slide, cfg: dict, theme: dict, chart_path: Path) -> None:
    set_bg(slide, theme["bg_alt"])
    add_title(slide, "phase 2：paired raw pileup / full 對照", "paired 問題下，raw caller 有 recall 空間，但 precision 代價過高，不適合直接當最終輸出", theme)
    add_rect(slide, 0.8, 1.8, 7.45, 4.95, fill="FFFFFF", line=theme["line"], radius=True)
    slide.shapes.add_picture(str(chart_path), Inches(1.0), Inches(2.0), width=Inches(7.0), height=Inches(3.85))
    add_text(slide, "圖意：X 軸為 F1，Y 軸為不同 caller 層級來源。左圖是 5kHz paired，右圖是 DORADO paired。", 1.0, 5.98, 6.9, 0.3, font_size=10.0, color="4F5B67", margin=0)
    add_rect(slide, 8.45, 1.8, 4.05, 4.95, fill="FFFFFF", line=theme["line"], radius=True)
    add_text(slide, "phase 2 結論", 8.45, 2.08, 2.0, 0.22, font_size=15, color=theme["text"], bold=True, margin=0)
    add_bullets(
        slide,
        cfg["deck"]["phase2"]["bullets"],
        8.62,
        2.45,
        3.25,
        2.9,
        font_size=11.0,
        color=theme["text"],
    )
    add_rect(slide, 8.62, 5.68, 3.2, 0.42, fill=theme["accent2"], line=None, radius=True)
    add_text(slide, cfg["deck"]["phase2"]["callout"], 8.8, 5.79, 2.85, 0.18, font_size=9.8, color=theme["light_text"], bold=True, margin=0)


def add_feature_roles_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["bg"])
    add_title(slide, "特徵定位與意義", "這一週已把主要特徵從『到處試規則』收斂成可制度化的角色分工", theme)
    positions = [(0.8, 1.9), (6.6, 1.9), (0.8, 3.5), (6.6, 3.5), (0.8, 5.1), (6.6, 5.1)]
    for (x, y), row in zip(positions, cfg["deck"]["feature_roles"]):
        add_rect(slide, x, y, 5.55, 1.42, fill="FFFFFF", line=theme["line"], radius=True)
        add_text(slide, row["feature"], x + 0.18, y + 0.1, 5.0, 0.2, font_size=11.8, color=theme["text"], bold=True, margin=0)
        add_text(slide, row["role"], x + 0.18, y + 0.34, 5.0, 0.16, font_size=10.2, color=theme["accent"], bold=True, margin=0)
        add_text(slide, f"意義：{row['meaning']}", x + 0.18, y + 0.60, 5.0, 0.24, font_size=9.8, color="4F5B67", margin=0)
        add_text(slide, f"限制：{row['limit']}", x + 0.18, y + 0.90, 5.0, 0.24, font_size=9.8, color=theme["muted"], margin=0)


def add_annotation_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["bg_alt"])
    add_title(slide, "annotation layer 決策", "phase 2 的 finer interval 已回接到 annotation 層，但目前仍不應直接改成 hard keep / hard filter", theme)
    positions = [(0.8, 1.95), (0.8, 3.55), (4.9, 1.95), (4.9, 3.55)]
    for (x, y), row in zip(positions, cfg["deck"]["annotation_policy"]["rows"]):
        add_rect(slide, x, y, 3.75, 1.3, fill="FFFFFF", line=theme["line"], radius=True)
        add_text(slide, row["dataset"], x + 0.18, y + 0.12, 3.0, 0.18, font_size=12.2, color=theme["text"], bold=True, margin=0)
        add_text(slide, row["policy"], x + 0.18, y + 0.4, 3.1, 0.18, font_size=10.4, color=theme["accent"], bold=True, margin=0)
        add_text(slide, f"ΔF1 vs baseline：{row['delta']}", x + 0.18, y + 0.68, 2.2, 0.18, font_size=10.0, color="4F5B67", margin=0)
        add_text(slide, row["interpretation"], x + 0.18, y + 0.92, 3.2, 0.2, font_size=9.9, color=theme["muted"], margin=0)
    cards = cfg["deck"]["annotation_policy"]["cards"]
    for idx, card in enumerate(cards):
        x = 9.15
        y = 1.95 + idx * 1.38
        add_rect(slide, x, y, 3.15, 1.05, fill=card["fill"], line=None, radius=True)
        add_text(slide, card["title"], x + 0.18, y + 0.14, 2.6, 0.2, font_size=12.5, color=card["title_color"], bold=True, margin=0)
        add_text(slide, card["body"], x + 0.18, y + 0.42, 2.55, 0.42, font_size=10.5, color=card["body_color"], margin=0)
    add_rect(slide, 0.8, 5.55, 11.7, 0.68, fill="FFFFFF", line=theme["line"], radius=True)
    add_text(slide, cfg["deck"]["annotation_policy"]["bottom_note"], 1.0, 5.77, 11.2, 0.22, font_size=10.3, color="4F5B67", margin=0)


def add_next_steps_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["dark"])
    add_title(slide, "待補齊事項與下一步", "現在最合理的做法是先把 annotation 接到輸出層，再決定是否改核心規則", theme, dark=True)
    add_rect(slide, 0.78, 1.9, 3.85, 4.72, fill="F4EFE7", line=None, radius=True)
    add_text(slide, "目前尚缺", 1.05, 2.15, 1.5, 0.22, font_size=15, color=theme["text"], bold=True, margin=0)
    add_bullets(slide, cfg["deck"]["next_steps"]["gaps"], 0.98, 2.55, 3.28, 3.6, font_size=10.8, color=theme["text"])

    add_rect(slide, 4.74, 1.9, 3.85, 4.72, fill="FFFFFF", line=None, radius=True)
    add_text(slide, "下一步", 5.2, 2.15, 1.5, 0.22, font_size=15, color=theme["text"], bold=True, margin=0)
    add_bullets(slide, cfg["deck"]["next_steps"]["actions"], 4.94, 2.55, 3.28, 3.6, font_size=10.8, color=theme["text"])

    add_rect(slide, 8.7, 1.9, 4.0, 4.72, fill=theme["accent"], line=None, radius=True)
    add_text(slide, "現階段不建議", 9.36, 2.15, 2.0, 0.22, font_size=15, color=theme["light_text"], bold=True, margin=0)
    add_bullets(slide, cfg["deck"]["next_steps"]["not_yet"], 8.98, 2.55, 3.28, 3.6, font_size=10.8, color=theme["light_text"])


def add_paths_slide(slide, cfg: dict, theme: dict) -> None:
    set_bg(slide, theme["bg"])
    add_title(slide, "複查路徑與再生成方式", "本簡報對應的來源、數據與客製化設定檔都可直接回查與重生", theme)
    add_rect(slide, 0.8, 1.82, 4.15, 4.75, fill="FFFFFF", line=theme["line"], radius=True)
    add_rect(slide, 5.12, 1.82, 3.55, 4.75, fill="FFFFFF", line=theme["line"], radius=True)
    add_rect(slide, 8.84, 1.82, 3.7, 4.75, fill="FFFFFF", line=theme["line"], radius=True)
    add_text(slide, "核心文件", 1.0, 2.08, 1.5, 0.22, font_size=15, color=theme["text"], bold=True, margin=0)
    add_text(slide, "主要數據", 5.32, 2.08, 1.5, 0.22, font_size=15, color=theme["text"], bold=True, margin=0)
    add_text(slide, "重生與版本", 9.04, 2.08, 1.8, 0.22, font_size=15, color=theme["text"], bold=True, margin=0)
    for idx, item in enumerate(cfg["deck"]["source_paths"][:3]):
        y = 2.45 + idx * 1.2
        add_text(slide, item["label"], 1.0, y, 1.8, 0.18, font_size=11.5, color=theme["accent"], bold=True, margin=0)
        add_text(slide, wrap_path_for_slide(item["path"], 30, 4), 1.0, y + 0.24, 3.55, 0.82, font_size=8.4, color="4F5B67", font_name="Consolas", margin=0)
    for idx, item in enumerate(cfg["deck"]["source_paths"][3:5]):
        y = 2.45 + idx * 1.55
        add_text(slide, item["label"], 5.32, y, 2.2, 0.18, font_size=10.6, color=theme["accent"], bold=True, margin=0)
        add_text(
            slide,
            wrap_path_for_slide(item["path"], 24, 5),
            5.32,
            y + 0.22,
            2.85,
            0.9,
            font_size=7.7,
            color="4F5B67",
            font_name="Consolas",
            margin=0,
        )
    version_name = Path(cfg["meta"]["output_pptx"]).with_suffix(".version.json").name
    add_text(slide, "版本確認檔", 9.04, 2.45, 1.8, 0.18, font_size=10.6, color=theme["accent"], bold=True, margin=0)
    add_text(slide, version_name, 9.04, 2.7, 2.9, 0.3, font_size=8.8, color="4F5B67", font_name="Consolas", margin=0)
    add_text(slide, "deck 入口", 9.04, 3.2, 1.2, 0.18, font_size=10.6, color=theme["accent"], bold=True, margin=0)
    add_text(slide, Path(cfg["meta"]["output_pptx"]).name, 9.04, 3.45, 3.0, 0.44, font_size=8.8, color="4F5B67", font_name="Consolas", margin=0)
    add_text(slide, "重生指令", 9.04, 4.08, 1.4, 0.18, font_size=10.6, color=theme["accent"], bold=True, margin=0)
    add_rect(slide, 9.04, 4.38, 3.05, 1.24, fill=theme["dark"], line=None, radius=True)
    add_text(
        slide,
        "python build_weekly_report_pptx.py --config\n20260311_weekly_report_pptx_config_01.json",
        9.24,
        4.58,
        2.55,
        0.78,
        font_size=7.7,
        color=theme["light_text"],
        font_name="Consolas",
        margin=0,
    )
    add_text(
        slide,
        "完整來源路徑請回查主週報與 version.json；本頁以可讀、可口頭報告的摘要形式呈現。",
        8.98,
        5.86,
        3.15,
        0.36,
        font_size=8.8,
        color="4F5B67",
        margin=0,
    )


def build_assets(cfg: dict) -> dict:
    assets_dir = Path(cfg["meta"]["assets_dir"])
    ensure_dir(assets_dir)
    theme = cfg["theme"]
    paired_chart = assets_dir / "paired_rerun_f1.png"
    create_grouped_bar_chart(
        paired_chart,
        "Paired Pure Rerun: F1 by Method",
        cfg["charts"]["paired_rerun"]["categories"],
        cfg["charts"]["paired_rerun"]["series"],
        theme,
        ylim=(0.84, 0.862),
        ylabel="F1",
    )
    phase2_chart = assets_dir / "phase2_paired_raw_benchmark.png"
    create_horizontal_bar_chart(phase2_chart, cfg["charts"]["phase2_raw"], theme)
    return {"paired_chart": paired_chart, "phase2_chart": phase2_chart}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    cfg = load_config(args.config)
    theme = cfg["theme"]
    assets = build_assets(cfg)
    prs = create_presentation(cfg)

    slide_builders = [
        lambda s: add_title_slide(s, cfg, theme),
        lambda s: add_key_findings_slide(s, cfg, theme),
        lambda s: add_framework_slide(s, cfg, theme),
        lambda s: add_goal_status_slide(s, cfg, theme),
        lambda s: create_timeline_slide(s, cfg, theme),
        lambda s: add_paired_rerun_slide(s, cfg, theme, assets["paired_chart"]),
        lambda s: add_to_rescue_slide(s, cfg, theme),
        lambda s: add_phase2_slide(s, cfg, theme, assets["phase2_chart"]),
        lambda s: add_feature_roles_slide(s, cfg, theme),
        lambda s: add_annotation_slide(s, cfg, theme),
        lambda s: add_next_steps_slide(s, cfg, theme),
        lambda s: add_paths_slide(s, cfg, theme),
    ]

    total_pages = len(slide_builders)
    for idx, builder in enumerate(slide_builders, start=1):
        slide = prs.slides.add_slide(prs.slide_layouts[6])
        builder(slide)
        add_footer(slide, theme, idx, total_pages, cfg["meta"]["footer"], dark=idx in {1, 11})

    output_path = Path(cfg["meta"]["output_pptx"])
    ensure_dir(output_path.parent)
    prs.save(output_path)
    print(f"saved_pptx={output_path}")
    print(f"assets_dir={cfg['meta']['assets_dir']}")


if __name__ == "__main__":
    main()
