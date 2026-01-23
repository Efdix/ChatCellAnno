# GitCell AI Instructions

GitCell is a privacy-first, "human-in-the-loop" Python library for single-cell annotation that uses the user's local GitHub Copilot Chat interface rather than external APIs.

## 🏗 Project Architecture

### Core Logic
- **Entry Point**: `gitcell.annotate_cell_types` in `gitcell/core.py`.
- **Data Flow**:
  1.  `extractor.py`: Extracts marker genes from `adata.uns['rank_genes_groups']`.
  2.  `prompt.py`: Formats markers into a prompt and puts it on the system **Clipboard**.
  3.  `parser.py`: Ingests Copilot's raw text response (pasted by user/script) and maps it to clusters.

### Scanpy Integration
- The library is designed to work *on top* of existing Scanpy objects (`anndata`).
- **Prerequisite**: Users must run `sc.tl.rank_genes_groups` before calling GitCell.
- **State Management**: Annotation metadata (like the grouping key) is retrieved from `adata.uns['rank_genes_groups']['params']['groupby']`.
- **Output**: Writes result to `adata.obs['gitcell_annotation']`.

## 💻 Developer Patterns & Workflows

### 1. Zero-API Design Pattern
- **Do not** add OpenAI/Anthropic API calls. The core philosophy is "Bring Your Own Copilot" (via the IDE).
- **Mechanism**: The "API" is the **OS Clipboard**.
  - Output: `pyperclip.copy(prompt_string)`
  - Input: User manual paste into `step='parse'` arguments.

### 2. Prompt Engineering Constraints
- The `parser.py` logic is intentionally simple (newline splitting).
- **Rule**: All prompts generated in `prompt.py` MUST instruct Copilot to:
  - Return **strictly** one line per cluster.
  - Exclude introductory text ("Here is the list...").
  - Exclude Markdown formatting that isn't a simple code block.
- **Handling Mismatches**: The parser warns if `len(response_lines) != len(clusters)`. When modifying prompts, always verify output format stability.

### 3. Error Handling
- **Clipboard**: `pyperclip` can fail on headless systems. Always wrap copy operations in `try/except` and fall back to printing the prompt to stdout.
- **Scanpy Data**: Always validate `adata.uns['rank_genes_groups']` exists before proceeding.

### 4. Testing & Verification
- Since there is no automated test suite yet, create mock `anndata` objects for verification:
  ```python
  import anndata
  import pandas as pd
  import numpy as np
  
  # Minimal mock
  obs = pd.DataFrame({'leiden': ['0', '1', '0', '1']}, index=['c1', 'c2', 'c3', 'c4'])
  adata = anndata.AnnData(np.random.rand(4, 10), obs=obs)
  # Mock the rank_genes_groups structure in adata.uns manually if skipping sc.tl call
  ```

## 📂 File Structure Key
- `gitcell/core.py`: Workflow orchestration.
- `gitcell/prompt.py`: Prompt generation & Clipboard ops.
- `gitcell/parser.py`: Text processing & Validation.
- `gitcell/extractor.py`: Scanpy data extraction interface.
