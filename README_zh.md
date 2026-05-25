# ChatCellAnno: AI 驱动的单细胞多模态注释助手

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)

ChatCellAnno 是一个用于单细胞注释的桌面应用。它将 marker 基因、可选表达矩阵、可选视觉上下文和功能富集证据整合到同一流程，辅助 LLM 输出更稳定的聚类注释结果。

[English Documentation](README.md)

## 核心功能

- 多模态上下文：支持 marker、图像和表达矩阵联合输入。
- 富集增强提示词：集成 GO/KEGG 等证据，降低注释歧义。
- 灵活交互：支持内置浏览器流程和可配置 API 参数。
- 结构化解析：将模型回复映射为 cluster -> cell type。
- 代码导出：生成可直接用于 Scanpy/Seurat 的注释代码片段。
- 模块化 GUI：PySide6 分层结构，便于维护和扩展。

## 项目结构

- `chatcellanno/core.py`：主流程编排（生成和解析）。
- `chatcellanno/extractor.py`：CSV/TSV marker 数据提取。
- `chatcellanno/prompt.py`：提示词构建。
- `chatcellanno/parser.py`：回复解析与代码生成。
- `chatcellanno/enrichment.py`：功能富集分析对接。
- `chatcellanno/gui/`：主窗口、UI 区块、后台线程。
- `gui.py`：桌面应用入口。

## 安装与运行

### 环境要求

- Python 3.10+

### 完整环境配置与安装

推荐使用 `conda` 创建隔离环境，确保包版本不会与其它项目冲突。

```bash
# 1. 创建并激活 conda 虚拟环境
conda create -n chatcellanno python=3.10 -y
conda activate chatcellanno

# 2. 从 requirements.txt 安装所有必要的依赖包
pip install -r requirements.txt
```

如果没有下载源码库，也可通过手动执行安装核心依赖：
```bash
pip install pandas pyperclip gseapy biopython matplotlib PySide6 PyInstaller
```

### 直接运行

环境配置好后，运行以下命令启动 GUI：
```bash
python gui.py
```

## Windows 独立环境打包

如果想将应用打包为可以直接运行的独立 `.exe` 文件（方便其他人直接下载双击运行，无需安装 Python 环境）：

```powershell
./build.ps1
```

构建成功后，在 `dist/` 目录下将生成一个独立的单文件应用程序 `ChatCellAnno.exe`。您可以将该 `.exe` 文件直接发送给他人使用。

## 标准使用流程

1. 载入 marker 文件（`.csv` 或 `.tsv`）。
2. 按需启用富集分析并补充图像/表达矩阵上下文。
3. 在步骤 2 生成提示词文件：浏览器模式可拖拽到内置对话框，API 模式可自动发送。
4. 在步骤 4（解析与报告）解析模型回复并导出 Scanpy/Seurat 代码与报告。

## 开发说明

- 在 `tests/` 下运行测试以验证解析和富集链路。
- `database/`、`genome/`、`results/` 等运行数据建议不纳入版本控制。

## 许可证

MIT。详见 [LICENSE](LICENSE)。
