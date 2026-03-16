# Contributing

## Environment

1. Create and activate environment.
2. Install dependencies.
3. Run local checks before commit.

```bash
conda create -n chatcellanno python=3.10
conda activate chatcellanno
pip install -e .[dev,gui]
```

## Development Rules

- Keep UI logic in `chatcellanno/gui/` and core logic in `chatcellanno/` modules.
- Add or update tests under `tests/` when changing parser/extractor/enrichment behavior.
- Avoid committing runtime data folders (`results/`, `chatcellanno_web_data/`, `database/`, `genome/`).
- Keep i18n strings centralized in `chatcellanno/config.py`.

## Suggested Pre-Commit Check

```bash
python -m compileall chatcellanno gui.py
pytest -q
```

## Pull Request Checklist

- [ ] Scope is focused and minimal.
- [ ] README/README_zh updated if behavior changed.
- [ ] New settings and defaults documented.
- [ ] No generated files or local runtime artifacts included.
