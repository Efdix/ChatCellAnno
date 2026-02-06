"""
ChatCellAnno - PySide6 GUI Entry Point

This file is now a wrapper around the modularized GUI package.
See chatcellanno/gui/ for the implementation.
"""

import sys
import os

# Ensure the package is importable
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

from chatcellanno.gui import run_app

if __name__ == "__main__":
    run_app()
