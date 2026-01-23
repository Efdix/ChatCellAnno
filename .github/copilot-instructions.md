# ChatCell AI Instructions

ChatCell is a privacy-first, "human-in-the-loop" Python library for single-cell annotation that uses the system clipboard to bridge Scanpy data with any LLM (Copilot, ChatGPT, etc.), adhering to a "Zero-API" design philosophy.

## 🏗 Project Architecture

### Core Logic
- **Entry Point**: `chatcell.annotate_cell_types` in [chatcell/core.py](chatcell/core.py).
- **Data Flow**:
  1.  **Extract**: [chatcell/extractor.py](chatcell/extractor.py) pulls top marker genes from `adata.uns['rank_genes_groups']`.
  2.  **Prompt**: [chatcell/prompt.py](chatcell/prompt.py) formats markers into a prompt string and pushes it to the **OS Clipboard**.
  3.  **Parse**: [chatcell/parser.py](chatcell/parser.py) ingests the raw text response (pasted by user/script) and maps annotations back to clusters.

### Scanpy Integration
- **Prerequisite**: Users MUST run `sc.tl.rank_genes_groups` before calling ChatCell.
- **State Management**: The clustering key is strictly retrieved from `adata.uns['rank_genes_groups']['params']['groupby']`.
- **Output**: Writes cell type annotations to `adata.obs['chatcell_annotation']` and optional metadata to `adata.obs['chatcell_extra_info']`.

## 💻 Developer Patterns & Workflows

### 1. Zero-API Design Pattern
- **Constraint**: Do NOT add OpenAI, Anthropic, or other REST API calls. The library must remain agnostic and privacy-centric.
- **Mechanism**: The "API" is the **OS Clipboard** (`pyperclip`).
  - **Output (Generate)**: `pyperclip.copy(prompt_string)`.
  - **Input (Parse)**: User manually passes text to `step='parse', response_text='...'`.

### 2. Prompt Engineering Protocol
All prompts in [chatcell/prompt.py](chatcell/prompt.py) must be deterministic to ensure parsing success:
- **Structure**: Instruct models to return *exactly* one line per cluster.
- **Modes**:
  - `concise`: Returns "CellType" only.
  - `evidence`: Returns "CellType | Supported by: gene1, gene2".
  - `recommendation`: Returns "CellType | Recommended Markers: geneA, geneB".
- **Forbidden**: Do not allow intro text ("Here are your annotations..."), markdown headers, or numbered lists unless stripping logic handles them.

### 3. Parser Resilience
- The parser in [chatcell/parser.py](chatcell/parser.py) must handle "chatty" LLM outputs.
- **Strategy**: Strip Markdown code blocks (` ``` `), ignore empty lines, and warn (don't crash) if `len(response_lines) != len(clusters)`.
- **Fallback**: If array lengths mismatch, pad with "Unknown" or truncate to avoid alignment errors.

### 4. Testing & Verification
- Since there is no full CI/CD suite, verify changes using **mock AnnData objects**.
- **Mock Template**:
  ```python
  import anndata
  import pandas as pd
  import numpy as np
  
  # 1. Mock Data & Obs
  obs = pd.DataFrame({'leiden': ['0', '1', '0', '1']}, index=['c1', 'c2', 'c3', 'c4'])
  adata = anndata.AnnData(np.random.rand(4, 10), obs=obs)
  
  # 2. Mock 'rank_genes_groups' Structure
  adata.uns['rank_genes_groups'] = {
      'names': pd.DataFrame([['CD3D', 'CD79A'], ['CD3E', 'CD19']], columns=['0', '1']),
      'scores': np.zeros((2, 2)),
      'pvals_adj': np.zeros((2, 2)),
      'params': {'groupby': 'leiden'}
  }
  ```

## 📂 File Structure Key
- [chatcell/core.py](chatcell/core.py): Workflow orchestration (`annotate_cell_types`).
- [chatcell/prompt.py](chatcell/prompt.py): Prompt generation text templates & clipboard operations.
- [chatcell/parser.py](chatcell/parser.py): LLM response parsing, validation, and cleaning.
- [chatcell/extractor.py](chatcell/extractor.py): Interface for reading Scanpy `uns` data structures.
