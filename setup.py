from setuptools import setup, find_packages

setup(
    name="chatcell",
    version="0.1.0",
    description="ChatCell: LLM-assisted Single-cell Annotation Tool",
    author="ChatCell Team",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "pyperclip",
    ],
    python_requires=">=3.8",
)
