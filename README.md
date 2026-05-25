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

### Complete Environment Setup and Installation

It is recommended to use `conda` to create an isolated environment to prevent conflicts with other python projects.

```bash
# 1. Create and activate a conda virtual environment
conda create -n chatcellanno python=3.10 -y
conda activate chatcellanno

# 2. Install all required dependencies from requirements.txt
pip install -r requirements.txt
```

If you don't have the source code cloned yet, you can also install the fundamental dependencies manually:
```bash
pip install pandas pyperclip gseapy biopython matplotlib PySide6 PyInstaller
```

### Direct Run

Once the environment is properly configured, run the following command to launch the GUI:
```bash
python gui.py
```

## Build Standalone Executable (Windows)

To package the application into a single, standalone `.exe` file (allowing anyone to download and run it directly without installing a Python environment):

```powershell
./build.ps1
```

After the build is successfully completed, a single executable file `ChatCellAnno.exe` will be generated in the `dist/` folder. You can directly share this `.exe` file with other users.

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
