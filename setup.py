from setuptools import setup, find_packages

setup(
    name="chatcellanno",
    version="0.1.0",
    description="ChatCellAnno: LLM-assisted Single-cell Annotation Tool",
    author="ChatCellAnno Team",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "pyperclip",
    ],
    python_requires=">=3.8",
)
