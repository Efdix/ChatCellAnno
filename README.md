# ChatCellAnno: Multimodal AI-Powered Single-Cell Annotation Assistant

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)

ChatCellAnno is a desktop application for single-cell annotation that combines marker genes, optional expression matrices, enrichment evidence, and visual context for LLM-assisted cluster labeling.

[中文文档](README_zh.md)

## Features

- Multimodal context: marker genes, optional UMAP image context, and optional expression matrix.
- Enrichment-assisted prompting: integrates GO/KEGG and related knowledge to reduce annotation ambiguity.
- Flexible interaction modes: built-in browser flow plus configurable API settings.
- Structured parsing: turns LLM responses into cluster annotations.
- Code export: generates Scanpy (Python) and Seurat (R) annotation snippets.
- Modular GUI: PySide6-based layout with separate core logic and UI modules.

## Project Structure

- `chatcellanno/core.py`: workflow orchestration (`generate` and `parse`).
- `chatcellanno/extractor.py`: marker extraction from CSV/TSV formats.
- `chatcellanno/prompt.py`: prompt construction.
- `chatcellanno/parser.py`: response parsing and code generation.
- `chatcellanno/enrichment.py`: enrichment analysis integration.
- `chatcellanno/gui/`: main window, UI blocks, and worker threads.
- `gui.py`: desktop app entry point.

## Installation

### Requirements

- Python 3.10+

### Quick Start

```bash
conda create -n chatcellanno python=3.10
conda activate chatcellanno
pip install pandas pyperclip gseapy biopython PySide6 PySide6-WebEngine
python gui.py
```

## Build (Windows)

```powershell
./build.ps1
```

The executable is generated under `dist/`.

## Typical Workflow

1. Load marker file (`.csv` or `.tsv`).
2. Optionally add enrichment and image/matrix context.
3. Generate a prompt file in Step 2, then either drag it to the built-in browser chat (manual mode) or send it automatically with API mode.
4. Parse the LLM response and export code/report in Step 4 (`Parsing & Report`).

## Development

- Run tests from `tests/` to verify parser and enrichment behavior.
- Keep large runtime data (`database/`, `genome/`, `results/`) out of source control.

## License

MIT. See [LICENSE](LICENSE).
