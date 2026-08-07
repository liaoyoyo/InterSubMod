#!/usr/bin/env python3
"""零依賴 PNG 編碼器（只用 stdlib zlib + struct）。

為什麼不用 ism_heatmap_std.svg_dual_panel()：實測單一位點的 inline SVG 是
1,533 KB（每個 cell 一個 <rect>，即使做了 run-length 合併）。同一份資料寫成
1 px = 1 cell 的 PNG 只要 19 KB —— 差 62 倍。genome-scale 下前者不可行。

顏色一律 import scripts/ism_heatmap_std.py 的函式，不自己另立一套色盤
（該檔是本專案的視覺規約 SoT，palette_drift_check.py 以它為準）。
"""
from __future__ import annotations

import struct
import zlib


def write_png(path, pixels, width: int, height: int) -> int:
    """pixels: list of (r,g,b) tuples，長度需為 width*height。回傳寫入 bytes。"""
    raw = bytearray()
    for y in range(height):
        raw.append(0)                       # filter type 0（None）
        base = y * width
        for x in range(width):
            r, g, b = pixels[base + x]
            raw += bytes((r, g, b))

    def chunk(tag: bytes, data: bytes) -> bytes:
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)   # 8-bit RGB
    body = (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", ihdr)
            + chunk(b"IDAT", zlib.compress(bytes(raw), 9))
            + chunk(b"IEND", b""))
    with open(path, "wb") as fh:
        fh.write(body)
    return len(body)


def hex_to_rgb(h: str):
    h = h.lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    return (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))


class Canvas:
    """極簡像素畫布。1 px = 1 個資料格，HTML 端用 image-rendering:pixelated 放大。"""

    def __init__(self, width: int, height: int, bg="#11161f"):
        self.w = width
        self.h = height
        self.px = [hex_to_rgb(bg)] * (width * height)

    def set(self, x: int, y: int, rgb):
        if 0 <= x < self.w and 0 <= y < self.h:
            self.px[y * self.w + x] = rgb

    def rect(self, x: int, y: int, w: int, h: int, rgb):
        for yy in range(y, min(y + h, self.h)):
            if yy < 0:
                continue
            base = yy * self.w
            for xx in range(max(x, 0), min(x + w, self.w)):
                self.px[base + xx] = rgb

    def vline(self, x: int, y0: int, y1: int, rgb):
        for yy in range(max(y0, 0), min(y1, self.h)):
            self.set(x, yy, rgb)

    def save(self, path) -> int:
        return write_png(path, self.px, self.w, self.h)
