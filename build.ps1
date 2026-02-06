# ChatCellAnno EXE Build Script
Write-Host "Checking dependencies..." -ForegroundColor Cyan

# Note: Ensure you are in the correct Conda environment (e.g., conda activate chatcellanno)

Write-Host "Building ChatCellAnno GUI Executable..." -ForegroundColor Green

# Ensure icon path is absolute to avoid PyInstaller issues
$ScriptPath = Split-Path -Parent $MyInvocation.MyCommand.Path
$IconPath = Join-Path $ScriptPath "app_icon.ico"

# Build command:
# --clean: Clear cache
# --onefile: Bundle into a single EXE
# --windowed: No console window when running
# --name: Name of the output file
# --icon: Path to icon file
# --hidden-import: Ensure dynamic imports are found

pyinstaller --noconfirm --clean --onefile --windowed `
    --name "ChatCellAnno" `
    --icon "$IconPath" `
    --add-data "database;database" `
    --hidden-import "PySide6.QtWebEngineWidgets" `
    --hidden-import "pandas" `
    --hidden-import "pyperclip" `
    --hidden-import "windnd" `
    --hidden-import "gseapy" `
    "gui.py"

Write-Host "Build Complete! Check the 'dist' folder." -ForegroundColor Green

