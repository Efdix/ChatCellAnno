# ChatCellAnno: AI-Powered Single-Cell Annotation Assistant 🧬🤖

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)

**ChatCellAnno** is a privacy-first, desktop application designed to bridge the gap between single-cell analysis data (Scanpy/Seurat) and Large Language Models (LLMs). Utilizing a **'Human-in-the-loop'** approach and a **'Zero-API'** design philosophy, it ensures your data remains local and secure while leveraging the power of modern AI.

[中文文档 (Chinese Documentation)](SOFTWARE_DESCRIPTION_CN.md)

---

## ✨ Core Features

- **🎨 Modern GUI**: Intuitive, drag-and-drop interface powered by PySide6.
- **🧠 Model Agostic**: Compatible with **any** LLM (Copilot, DeepSeek, ChatGPT, Claude, etc.) via smart clipboard interaction.
- **🔒 Privacy First**: No API keys required. No network requests are made by the application core.
- **📊 Intelligent Parsing**: Automatically detects Scanpy (`.csv`) and Seurat (`.tsv`) marker file formats.
- **🔬 Functional Enrichment**: Built-in ORA enrichment analysis (Online via Enrichr or Local via GMT files). Uses database evidence (GO, KEGG, CellMarker) to minimize AI hallucinations.
- **🧬 Genome Analysis**: Visualize and export multi-species protein alignments with MEGA-style coloring.
- **💻 Automated Code Generation**: Parses AI responses to generate executable Python (Scanpy) or R (Seurat) code for cell type annotation.

---

## 🏗 Architecture

- **Core Logic**: The `chatcellanno/` package handles data extraction, prompt engineering, and response parsing.
- **GUI**: `gui.py` provides the user interface with an integrated QtWebEngine browser.
- **Zero-API**: All data transfer between the app and the LLM is handled exclusively via the system clipboard.

---

## 🚀 Installation & Run

### Prerequisites
- Python 3.8+
- Dependencies: Pandas, PySide6, PySide6-WebEngine, Pyperclip, Biopython

### Quick Start

1. **Create a virtual environment (recommended):**
   ```bash
   conda create -n chatcellanno python=3.10
   conda activate chatcellanno
   ```

2. **Install dependencies:**
   ```bash
   pip install pandas pyperclip PySide6 PySide6-WebEngine biopython gseapy
   ```

3. **Run the application:**
   ```bash
   python gui.py
   ```

### 📦 Building Executable (.exe)

To distribute ChatCellAnno as a standalone `.exe` file on Windows:

1. Install PyInstaller:
   ```bash
   pip install pyinstaller auto-py-to-exe
   ```

2. Run the provided build script:
   ```powershell
   ./build.ps1
   ```
   *Or manually:*
   ```bash
   pyinstaller ChatCellAnno.spec
   ```
3. The executable will be found in `dist/ChatCellAnno.exe`.

---

## 流程 Workflow

1. **Step 1: Load Data**: Drag and drop your marker gene file (`.csv`, `.tsv`) into the drop zone.
   - **(Optional)** Load a **Gene Expression Matrix** (Cluster-by-Gene mean expression) to provide focus on specific marker candidates.
2. **Step 2: Enrichment (Optional)**: Select Online or Local enrichment to cross-validate cell identities.
3. **Step 3: Configure & Generate**: 
   - Set species and analysis mode. 
   - **(New!) Visual Context**: Paste a screenshot of your UMAP/t-SNE plot to give the AI spatial context.
   - Click **'Generate & Copy Prompt'**.
4. **Step 4: AI Interaction**: Paste the prompt into your favorite LLM (e.g., DeepSeek, ChatGPT). If you added a UMAP image, click "Copy Image" and paste it into the chat too.
5. **Step 5: Parse & Export**: Copy the AI's markdown table response, paste it back, and click **'Process AI Output'**.

---

## 📄 License

MIT License. See [LICENSE](LICENSE) for more information.
