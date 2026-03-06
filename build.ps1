# ChatCellAnno EXE Build Script
Write-Host "Checking dependencies..." -ForegroundColor Cyan

# Note: Ensure you are in the correct Conda environment (e.g., conda activate chatcellanno)

Write-Host "Building ChatCellAnno GUI Executable..." -ForegroundColor Green

# Ensure icon path is absolute to avoid PyInstaller issues
$ScriptPath = Split-Path -Parent $MyInvocation.MyCommand.Path
$IconPath = Join-Path $ScriptPath "app_icon.ico"

# Build command:
# --clean: Clear cache
# --onedir: Bundle into a directory (faster launch and build than --onefile)
# --windowed: No console window when running

pyinstaller --noconfirm --clean --windowed `
    --name "ChatCellAnno" `
    --icon "$IconPath" `
    --add-data "database;database" `
    --hidden-import "PySide6.QtWebEngineWidgets" `
    --hidden-import "pandas" `
    --hidden-import "pyperclip" `
    --hidden-import "gseapy" `
    --exclude-module "PySide6.QtQml" `
    --exclude-module "tkinter" `
    "gui.py"

Write-Host "Build Complete! Check the 'dist' folder." -ForegroundColor Green

