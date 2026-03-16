from pathlib import Path
from setuptools import setup, find_packages

ROOT = Path(__file__).parent
README = (ROOT / "README.md").read_text(encoding="utf-8")

setup(
    name="chatcellanno",
    version="0.1.0",
    description="ChatCellAnno: LLM-assisted Single-cell Annotation Tool",
    long_description=README,
    long_description_content_type="text/markdown",
    author="ChatCellAnno Team",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "pyperclip",
        "gseapy",
        "biopython",
    ],
    extras_require={
        "gui": ["PySide6", "PySide6-WebEngine"],
        "dev": ["pytest", "pytest-mock"],
    },
    entry_points={
        "console_scripts": [
            "chatcellanno-gui=chatcellanno.gui:run_app",
        ],
    },
    python_requires=">=3.10",
)
