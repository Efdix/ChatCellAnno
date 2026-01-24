# ChatCell EXE Build Script
Write-Host "Checking dependencies..." -ForegroundColor Cyan

# Install requirements if missing
pip install pyinstaller windnd pyperclip pandas

Write-Host "Building ChatCell GUI Executable..." -ForegroundColor Green

# Build command:
# --onefile: Bundle into a single EXE
# --windowed: No console window when running
# --name: Name of the output file
# --hidden-import: Ensure some data libraries are hooked correctly (scanpy dependencies)
# --collect-all: Scanpy/Anndata sometimes need all their metadata

pyinstaller --noconfirm --onefile --windowed `
    --name "ChatCellAnno" `
    --hidden-import "pandas" `
    --hidden-import "pyperclip" `
    --hidden-import "windnd" `
    "gui.py"

Write-Host "Build complete! Check the 'dist' folder for ChatCellAnno.exe." -ForegroundColor Yellow
