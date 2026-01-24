# ChatCell AI Instructions

ChatCell is a privacy-first, "human-in-the-loop" Python library for single-cell annotation that uses the system clipboard to bridge Scanpy data with any LLM (Copilot, ChatGPT, etc.), adhering to a "Zero-API" design philosophy.

## 🏗 Project Architecture

### Core Logic
- **Entry Point**: `annotate_cell_types` in [chatcellanno/core.py](chatcellanno/core.py).
- **Data Flow**:
  1.  **Extract**: [chatcellanno/extractor.py](chatcellanno/extractor.py) pulls top marker genes from `adata.uns['rank_genes_groups']` or CSV/TSV.
  2.  **Prompt**: [chatcellanno/prompt.py](chatcellanno/prompt.py) formats markers into a prompt string and pushes it to the **OS Clipboard** (`pyperclip`).
  3.  **Parse**: [chatcellanno/parser.py](chatcellanno/parser.py) ingests raw text response and maps annotations back to clusters.

### GUI & Distribution
- **GUI**: [gui.py](gui.py) provides a Tkinter-based interface for non-programmers.
- **Build System**: [build.ps1](build.ps1) uses PyInstaller to bundle the application.
- **Distribution**: No user installation required; users run the standalone EXE.

### Data Integration
- **Input**: CSV/TSV files containing marker genes.
- **Support**: Tidy format (Sequence/Scanpy output) and Wide format (Matrix).

## 💻 Developer Protocols

### 1. Zero-API Design Pattern
- **Constraint**: Do NOT add direct API calls (e.g., OpenAI SDK). The library must remain agnostic.
- **Mechanism**: Use the clipboard for data transfer. 
  - **Output**: `pyperclip.copy(prompt_string)`.

### 2. Prompt Engineering Protocol
Templates in [chatcellanno/prompt.py](chatcellanno/prompt.py) must be deterministic:
- **Structure**: Instruct models to return *exactly* one line per cluster.
- **Modes**:
  - `concise`: Returns "ClusterX: CellType" only.
  - `detailed`: Returns "ClusterX: CellType | Recommended Markers | Functions | Ranks".
- **Forbidden**: No intro/outro text; No Markdown headers.

### 3. Parser Resilience
- **Strategy**: Strip Markdown blocks, ignore empty lines, and warn if `len(response_lines) != len(clusters)`.
- **Fallback**: Pad with "Unknown" or truncate to maintain alignment.

### 4. Testing & Verification
Verify changes using mock objects in [examples/gen_data.py](examples/gen_data.py) or custom snippets:
```python
import anndata
import pandas as pd
import numpy as np
obs = pd.DataFrame({'leiden': ['0', '1', '0', '1']}, index=['c1', 'c2', 'c3', 'c4'])
adata = anndata.AnnData(np.random.rand(4, 10), obs=obs)
adata.uns['rank_genes_groups'] = {
    'names': pd.DataFrame([['CD3D', 'CD79A'], ['CD3E', 'CD19']], columns=['0', '1']),
    'params': {'groupby': 'leiden'}
}
```

## 📂 File Structure Key
- [chatcellanno/core.py](chatcellanno/core.py): Workflow orchestration.
- [chatcellanno/prompt.py](chatcellanno/prompt.py): Templates and clipboard I/O.
- [chatcellanno/parser.py](chatcellanno/parser.py): Response cleaning and mapping.
- [chatcellanno/extractor.py](chatcellanno/extractor.py): Scanpy/TSV data extraction.

