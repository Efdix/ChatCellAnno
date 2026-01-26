# ChatCell AI - Copilot Instructions

ChatCellAnno is a privacy-first, "human-in-the-loop" Python library and desktop application for single-cell annotation. It bridges Scanpy/Seurat data with any LLM (Copilot, ChatGPT, DeepSeek) using the **System Clipboard** as the API, strictly adhering to a "Zero-API" design philosophy.

## 🏗 Project Architecture

### Core Logic (`chatcellanno/`)
*   **Entry Point**: `annotate_cell_types` in [chatcellanno/core.py](chatcellanno/core.py) orchestrates the "Generate" and "Parse" phases.
*   **Data Extraction**: [chatcellanno/extractor.py](chatcellanno/extractor.py) extracts markers from CSV/TSV.
    *   **Scanpy Mode**: Looks for columns `names`/`gene`/`symbol` (gene) and `group`/`cluster` (cluster).
    *   **Seurat Mode**: Looks for `gene`/`feature` and `cluster`.
*   **Prompt Engineering**: [chatcellanno/prompt.py](chatcellanno/prompt.py) formats markers into deterministic text templates and pushes to clipboard.
*   **Response Parsing**: [chatcellanno/parser.py](chatcellanno/parser.py) maps raw text response lines back to clusters.

### GUI & Distribution
*   **Entry Point**: [gui.py](gui.py) is the Tkinter interface.
    *   *Note*: Handles `sys.path` injection to allow running as a script (dev) vs frozen exe (prod).
*   **Clipboard**: Central data bus. No network requests are made by the app.
*   **Drag & Drop**: Uses `windnd` (Windows only) for better UX, wrapped in try-except for portability.

## 🚀 Critical Developer Workflows

### 1. Build Process
Use the PowerShell script to generate the standalone executable.
```powershell
./build.ps1
```
*   **Tool**: PyInstaller
*   **Configuration**: `--onefile --windowed --hidden-import pandas ...`
*   **Output**: `dist/ChatCellAnno.exe`

### 2. Testing & Verification
Since there is no external API to mock, testing focuses on **Parser Resilience**.
*   **Parser Logic**: verify `parser.py` handles markdown bloat (```json ... ```) and empty lines correctly.
*   **Mocking**: Use synthetic DataFrames similar to manual construction (e.g. `pd.DataFrame({'cluster': ...})`).

## 💻 Codebase Conventions

### 1. Zero-API Design Pattern
*   **Constraint**: NEVER implement direct calls to OpenAI/Anthropic/Google APIs.
*   **Mechanism**: All data transfer must happen via `pyperclip.copy()` and user pasting.
*   **Reasoning**: Privacy, model-agnosticism, and "serverless" operation.

### 2. Prompt Engineering Protocol
Templates in [chatcellanno/prompt.py](chatcellanno/prompt.py) must be **Deterministic**.
*   **Constraint**: The LLM must be instructed to return *exactly one line per cluster*.
*   **Warning**: If the Prompt allows multi-line descriptions per cluster, the `parser.py` logic (line-to-cluster mapping) will fail.
*   **Format**: "Concise" (Type only) or "Detailed" (Type | Reasoning).

### 3. Dependency Management
*   **GUI Imports**: `gui.py` must remain runnable as a standalone script `python gui.py` without installing the package. Do not remove the `sys.path.append` block.
*   **Platform Specifics**: `windnd` is Windows-specific. Always guard it with `try-except ImportError`.

## 📂 Key File Map
*   [chatcellanno/core.py](chatcellanno/core.py): Workflow controller.
*   [chatcellanno/extractor.py](chatcellanno/extractor.py): Pandas logic for Scanpy/Seurat formats.
*   [gui.py](gui.py): Tkinter frontend (Main Thread).
*   [build.ps1](build.ps1): Release build script.

