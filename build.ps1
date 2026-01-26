# ChatCell EXE Build Script
Write-Host "Checking dependencies..." -ForegroundColor Cyan

# Note: Ensure you are in the correct Conda environment (e.g., conda activate chatcellanno)

Write-Host "Building ChatCell GUI Executable..." -ForegroundColor Green

# Build command:
# --onefile: Bundle into a single EXE
# --windowed: No console window when running
# --name: Name of the output file
# --hidden-import: Ensure dynamic imports are found

pyinstaller --noconfirm --onefile --windowed `
    --name "ChatCellAnno" `
    --hidden-import "pandas" `
    --hidden-import "pyperclip" `
    --hidden-import "windnd" `
    "gui.py"

Write-Host "Build Complete! Check the 'dist' folder." -ForegroundColor Green

