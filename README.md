# ChatCellAnno: AI-Powered Single-Cell Annotation Assistant

ChatCellAnno is a privacy-first, desktop application designed to bridge the gap between single-cell analysis data (Scanpy/Seurat) and Large Language Models (LLMs). Utilizing a 'Human-in-the-loop' approach and a 'Zero-API' design philosophy, it ensures your data remains local and secure while leveraging the power of modern AI.

[中文文档 (Chinese Documentation)](README_zh.md)

## Core Features

- **GUI Interface**: Modern, drag-and-drop interface powered by PySide6.
- **Model Agnostic**: Compatible with any LLM (Copilot, DeepSeek, ChatGPT, Claude, etc.) via clipboard interaction.
- **Privacy First**: No API keys required. No network requests are made by the application core.
- **Intelligent Parsing**: Automatically detects Scanpy and Seurat marker file formats.
- **Functional Enrichment**: Built-in ORA enrichment analysis (Online via Enrichr or Local via GMT files). Uses database evidence (GO, KEGG, CellMarker) to minimize AI hallucinations.
- **Automated Code Generation**: Parses AI responses to generate executable Python (Scanpy) or R (Seurat) code for cell type annotation.

## Architecture

- **Core Logic**: The \chatcellanno/\ package handles data extraction, prompt engineering, and response parsing.
- **GUI**: \gui.py\ provides the user interface with an integrated QtWebEngine browser.
- **Zero-API**: All data transfer between the app and the LLM is handled exclusively via the system clipboard.

## Installation & Run

### Prerequisites
- Python 3.10+
- Dependencies: Pandas, PySide6, PySide6-WebEngine, Pyperclip

### Quick Start

1. Create a virtual environment (recommended):
   ```bash
   conda create -n chatcellanno python=3.10
   conda activate chatcellanno
   ```

2. Install dependencies:
   ```bash
   pip install pandas pyperclip PySide6 PySide6-WebEngine
   ```

3. Run the application:
   ```bash
   python gui.py
   ```
## Workflow

1. **Load Data**: Drag and drop your marker gene file (.csv, .tsv) into Step 1.
2. **Enrichment (Optional)**: Select Online or Local enrichment in Step 2. You can drag a folder with database files into the app to load local libraries.
3. **Configure & Generate**: Set species and analysis mode in Step 3. Click 'Generate & Copy Prompt'.
4. **AI Browser**: Paste and send the prompt to an LLM using the integrated browser.
5. **Parse & Export**: Copy the AI's markdown table response, paste it into Step 4, and click 'Process AI Output' to get the final R/Python code.

## License

MIT License
