# Release Checklist

## Repository Hygiene

- [ ] `git status` contains only intended changes.
- [ ] `.gitignore` excludes build/runtime/local data.
- [ ] No secrets in tracked files.

## Functional Verification

- [ ] App starts via `python gui.py`.
- [ ] Core workflow runs: load marker -> generate -> parse -> export code.
- [ ] Tests pass (`pytest -q`).

## Build Verification

- [ ] Windows build succeeds with `./build.ps1`.
- [ ] Output exists under `dist/`.

## Documentation

- [ ] English and Chinese README are consistent.
- [ ] Version and license info are accurate.
