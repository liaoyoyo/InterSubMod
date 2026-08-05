#!/usr/bin/env python3
"""Deprecated compatibility entrypoint for read-AF tree ordering.

The old name used CCF for an uncorrected family-specific read ALT fraction.
Use read_af_tree_ordering_multisample.py and its explicit CLI instead.
"""

from read_af_tree_ordering_multisample import main


if __name__ == "__main__":
    main()
