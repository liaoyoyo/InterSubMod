#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Font configuration module for matplotlib CJK support.

This module provides consistent font configuration for all plotting scripts
in InterSubMod, ensuring proper display of Chinese characters.

Usage:
    from font_config import configure_matplotlib_fonts
    configure_matplotlib_fonts()
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
import os


def configure_matplotlib_fonts():
    """
    Configure matplotlib to use CJK-compatible fonts.

    This function should be called at the beginning of any plotting script
    to ensure proper rendering of Chinese, Japanese, and Korean characters.

    Font priority:
    1. Noto Sans CJK TC (Traditional Chinese, installed in Docker)
    2. Microsoft JhengHei (Windows)
    3. PingFang TC (macOS)
    4. WenQuanYi Micro Hei (Linux fallback)
    5. DejaVu Sans (final fallback)
    """
    # CJK font priority list
    cjk_fonts = [
        'Noto Sans CJK TC',      # Docker / Linux
        'Noto Sans CJK SC',      # Simplified Chinese fallback
        'Microsoft JhengHei',    # Windows Traditional Chinese
        'Microsoft YaHei',       # Windows Simplified Chinese
        'PingFang TC',           # macOS Traditional Chinese
        'PingFang SC',           # macOS Simplified Chinese
        'WenQuanYi Micro Hei',   # Linux fallback
        'SimHei',                # Windows fallback
        'DejaVu Sans',           # Final fallback
        'sans-serif'
    ]

    # Suppress font warnings
    warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

    # Configure matplotlib
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = cjk_fonts
    plt.rcParams['axes.unicode_minus'] = False  # Fix minus sign display

    # Set default figure parameters
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 150
    plt.rcParams['figure.figsize'] = [10, 8]

    # Use Agg backend for headless environments (Docker)
    if os.environ.get('MPLBACKEND') == 'Agg' or not os.environ.get('DISPLAY'):
        mpl.use('Agg')


def get_available_cjk_font():
    """
    Find the first available CJK font on the system.

    Returns:
        str: Name of the available CJK font, or 'sans-serif' if none found.
    """
    from matplotlib.font_manager import fontManager

    cjk_fonts = [
        'Noto Sans CJK TC',
        'Noto Sans CJK SC',
        'Microsoft JhengHei',
        'PingFang TC',
        'WenQuanYi Micro Hei',
    ]

    available_fonts = set(f.name for f in fontManager.ttflist)

    for font in cjk_fonts:
        if font in available_fonts:
            return font

    return 'sans-serif'


# Auto-configure when imported
configure_matplotlib_fonts()


if __name__ == '__main__':
    # Test font configuration
    print("Testing font configuration...")

    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(8, 6))

    x = np.linspace(0, 10, 100)
    y = np.sin(x)

    ax.plot(x, y, label='正弦波 (Sine Wave)')
    ax.set_title('InterSubMod 中文字體測試')
    ax.set_xlabel('X 軸標籤')
    ax.set_ylabel('Y 軸標籤')
    ax.legend()
    ax.grid(True)

    output_path = '/tmp/font_test.png'
    plt.savefig(output_path, bbox_inches='tight')
    print(f"Test plot saved to: {output_path}")

    available_font = get_available_cjk_font()
    print(f"Available CJK font: {available_font}")
