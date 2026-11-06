#!/usr/bin/env python3
"""Compatibility shim. The command is `popseq`; the logic is in popseq/cli.py.

Kept because specification 007 and the README advertised this path first. It
forwards unchanged, including the exit code.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from popseq.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
