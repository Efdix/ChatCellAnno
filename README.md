# ChatCellAnno: AI-Powered Single-Cell Annotation Assistant

ChatCellAnno is a privacy-first, desktop application designed to bridge the gap between single-cell analysis data (Scanpy/Seurat) and Large Language Models (LLMs). Utilizing a 'Human-in-the-loop' approach and a 'Zero-API' design philosophy, it ensures your data remains local and secure while leveraging the power of modern AI.

[中文文档 (Chinese Documentation)](README_zh.md)

## Core Features

- **GUI Interface**: Modern, drag-and-drop interface powered by PySide6.
- **Model Agnostic**: Compatible with any LLM (Copilot, DeepSeek, ChatGPT, Claude, etc.) via clipboard interaction.
- **Privacy First**: No API keys required. No network requests are made by the application core.
- **Intelligent Parsing**: Automatically detects Scanpy (`names`, `group`) and Seurat (`gene`, `cluster`) marker file formats.
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

1. **Load Data**: Drag and drop your marker gene file (.csv, .tsv) into the application.
2. **Configure**: Enter species, tissue type, and other parameters to refine the context.
3. **Generate Prompt**: Click 'Generate & Copy Prompt'.
4. **Interact with AI**: Paste the prompt into the built-in browser (or any external LLM interface).
5. **Process Response**: Copy the AI's full response and paste it back into ChatCellAnno.
6. **Get Code**: Click 'Process AI Output' to generate the annotation code.

## License

MIT License
