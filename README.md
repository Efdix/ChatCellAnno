# ChatCellAnno: Multimodal AI-Powered Single-Cell Annotation Assistant 🧬🤖

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)

**ChatCellAnno** is a powerful **Multimodal & Hybrid-Interaction** desktop application designed to bridge the gap between complex single-cell omics data (Scanpy/Seurat) and frontier Large Language Models (LLMs). By integrating **Marker Genes**, **Expression Matrices**, and **Visual Context (UMAP/TSNE)**, it provides a high-fidelity grounding for AI-assisted cell type identification while maintaining a flexible "Human-in-the-loop" workflow.

[中文文档 (Chinese Documentation)](SOFTWARE_DESCRIPTION_CN.md)

---

## ✨ Core Features

- **🎨 Multimodal Grounding**: Go beyond simple gene lists. Feed your AI with **Visual Context (UMAP images)** and **Expression Matrices** to resolve ambiguous clusters.
- **🔗 Hybrid Interaction Modes**: 
  - **API Mode**: Directly connect to LLM providers (DeepSeek, OpenAI, SiliconFlow) for high-speed automated pipelines.
  - **Browser Mode**: Use the integrated secure browser for a guided, manual-chat experience with privacy control.
- **🔬 RAG-Enhanced Prompting**: Built-in **Functional Enrichment (GSEApy)** automatically injects biological pathway evidence (GO/KEGG) into prompts to minimize AI hallucinations.
- **🧩 Hot-Swap Plugin System**: Extend the software's capabilities (e.g., Sequence Alignment, CellChat Integration) via a dynamic drag-and-drop plugin architecture.
- **📄 One-Click Analysis Reports**: Automatically generate comprehensive Markdown reports summarizing annotations, enrichment results, and project metadata.
- **💻 Automated Code Generation**: Seamlessly convert AI text responses into executable **Python (Scanpy)** or **R (Seurat)** code snippets.
- **🔒 Flexible Data Persistence**: Portable configuration management designed to run from any folder, preserving your settings across updates.

---

## 🏗 Architecture

- **Core Logic**: The `chatcellanno/` package handles multi-source data extraction, expert-level prompt engineering, and structured response parsing.
- **GUI Engine**: A highly modular PySide6 interface featuring a split-pane design for concurrent data inspection and AI interaction.
- **Plugin Management**: Decoupled framework for runtime function injection without modifying the core codebase.

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

## 🚀 Workflow

1. **Step 1: Data & Enrichment Loading**: 
   - **Differential Gene List (Required)**: Drag and drop your marker gene files (`.csv`, `.tsv`) into the designated area.
   - **(New!) Assisting AI with Enrichment (Optional)**: Check this option to perform local or online (Enrichr) enrichment. Bioinformatic evidence will be automatically written into the prompt to reduce hallucinations.
   - **UMAP/t-SNE Screenshot (Optional)**: Paste a cluster visualization screenshot to provide crucial spatial context.

2. **Step 2: Prompt Configuration**:
   - Set the species, tissue, and output mode (Concise or Detailed).
   - Click **'Generate & Copy Prompt'**. If a screenshot is added, you will be prompted to copy the image simultaneously.

3. **Step 3: AI Interaction**:
   - Paste the generated prompt and image into your favorite LLM (e.g., DeepSeek, ChatGPT, Claude).
   - Use the **Built-in Browser** panel to handle conversations directly within the app.

4. **Step 4: Parse & Export**:
   - Copy the Markdown table response from the AI and paste it into the parsing area.
   - Click **'Process AI Output'** to obtain executable **Python (Scanpy)** or **R (Seurat)** annotation scripts.

## 🛠 Workspace Management

- **Flexible Tabs**: All functional panels (Enrichment, Browser, Genome) on the right can be **closed** or **rearranged** via drag-and-drop.
- **Restore Panels**: If a panel is accidentally closed, reopen it through the **'Show Panels'** menu under the top-right gear icon.

---

## 📄 License

MIT License. See [LICENSE](LICENSE) for more information.
