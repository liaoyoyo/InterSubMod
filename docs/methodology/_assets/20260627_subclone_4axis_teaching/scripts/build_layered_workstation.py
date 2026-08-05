#!/usr/bin/env python3
"""Compatibility entrypoint for the canonical-v5 layered workstation renderer.

The maintained implementation lives in build_layered_workstation_v5.py.  Keep
this filename so existing repository commands fail closed into the same
hash-bound renderer instead of silently rebuilding a historical interface.
"""

from build_layered_workstation_v5 import main


if __name__ == "__main__":
    main()
