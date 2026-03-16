# ChatCellAnno AI - Copilot Instructions

ChatCellAnno is a multimodal Python library and desktop application for single-cell annotation.

## Project Architecture

### Core Logic (`chatcellanno/`)
- Entry point: `annotate_cell_types` in `chatcellanno/core.py` orchestrates generate/parse workflow.
- Data extraction: `chatcellanno/extractor.py` handles Scanpy/Seurat marker tables from CSV/TSV.
- Prompt engineering: `chatcellanno/prompt.py` composes multimodal prompts from markers, tissue/species, enrichment hints, and optional visual/matrix context.
- Response parsing: `chatcellanno/parser.py` parses LLM responses and exports annotation code for Scanpy/Seurat.
- Enrichment: `chatcellanno/enrichment.py` provides optional biological evidence for prompting.

### GUI and Distribution
- GUI entry: `gui.py` -> `chatcellanno/gui/` modules.
- Main window: `chatcellanno/gui/main_window.py`.
- UI composition: `chatcellanno/gui/ui_blocks.py`.
- Worker threads: `chatcellanno/gui/workers.py`.
- Build pipeline: `build.ps1` and `ChatCellAnno.spec`.

## Developer Workflow

### Build

```powershell
./build.ps1
```

### Test Focus

- Parser robustness against diverse markdown/table outputs.
- Extractor compatibility with varied marker column names.
- Enrichment flow behavior for local and online sources.

## Code Conventions

- Keep parser output format deterministic and testable.
- Keep GUI non-blocking for long tasks via worker/thread mechanisms.
- Preserve compatibility of both browser-driven and API-driven query flows.
- Keep optional dependencies and plugin failures graceful.

## Key Files

- `chatcellanno/core.py`
- `chatcellanno/extractor.py`
- `chatcellanno/prompt.py`
- `chatcellanno/parser.py`
- `chatcellanno/enrichment.py`
- `chatcellanno/gui/main_window.py`
- `build.ps1`
