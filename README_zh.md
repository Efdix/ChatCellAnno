# ChatCellAnno: AI 驱动的单细胞注释助手

ChatCellAnno 是一款注重隐私的桌面应用程序，旨在连接单细胞分析数据 (Scanpy/Seurat) 与大语言模型 (LLM)。它采用“人机回环 (Human-in-the-loop)”和“零 API (Zero-API)”的设计理念，确保您的数据安全且无需配置复杂的 API 密钥。

[English Documentation](README.md)

## 核心功能

- **图形界面**: 基于 PySide6 的现代化拖拽式界面，简单易用。
- **模型无关**: 支持任意 LLM (Copilot, DeepSeek, ChatGPT, Claude 等)，通过系统剪贴板进行交互。
- **隐私优先**: 无需 API Key，核心程序不发起任何网络请求，确保数据隐私。
- **智能解析**: 自动识别 Scanpy 和 Seurat 格式的 Marker 文件。
- **功能富集分析**: 集成了本地和在线 (Enrichr) 富集分析功能。通过数据库（GO, KEGG, CellMarker 等）提供的客观证据，辅助 AI 减少幻觉并提高注释权威性。
- **代码生成**: 能够将 AI 的自然语言回复转换为可执行的 Python (Scanpy) 或 R (Seurat) 代码，直接用于细胞类型注释。

## 架构设计

- **核心逻辑**: `chatcellanno/` 目录负责数据提取、提示词构建和结果解析。
- **用户界面**: `gui.py` 提供了基于 PySide6 的图形界面，集成了 QtWebEngine 浏览器。
- **零 API 交互**: 应用程序与 LLM 之间的数据传输完全通过系统剪贴板完成，解耦了特定模型的 API 依赖。

## 安装与运行

### 依赖项
- Python 3.10+
- Pandas, PySide6, PySide6-WebEngine, Pyperclip

### 快速开始
1. 创建虚拟环境 (推荐):
   ```bash
   conda create -n chatcellanno python=3.10
   conda activate chatcellanno
   ```

2. 安装依赖:
   ```bash
   pip install pandas pyperclip PySide6 PySide6-WebEngine
   ```

3. 运行程序:
   ```bash
   python gui.py
   ```

### 构建可执行文件

1. 确保已安装 PyInstaller:
   ```bash
   pip install pyinstaller
   ```
2. 运行构建脚本 (PowerShell):
   ```powershell
   ./build.ps1
   ```
   或者手动构建:
   ```bash
   pyinstaller ChatCellAnno.spec
   ```
3. 可执行文件将在 `dist/ChatCellAnno.exe` 中生成。


## 使用流程

1. **步骤 1: 加载数据**: 将您的 Marker 基因文件 (.csv, .tsv) 拖入软件界面指定区域。
2. **步骤 2: 富集分析 (可选)**: 选择在线或本地富集分析以交叉验证细胞身份。
3. **步骤 3: 提示词生成**:
   - 设置物种、组织及模式（简洁/详细）。
   - **(新功能!) 视觉上下文 (Visual Context)**: 从您的 Scanpy/Seurat 结果中截图 UMAP/t-SNE 聚类图，直接粘贴至软件中（Ctrl+V）。
   - 点击 **“生成并复制提示词”**。
4. **步骤 4: AI 交互**: 将生成的提示词粘贴到您的大模型对话框 (如 ChatGPT, DeepSeek)。如果您添加了 UMAP 图像，请点击 **“复制图像”** 并将其一并发送给 AI。
5. **步骤 5: 解析与导出**: 复制 AI 返回的 Markdown 表格，粘贴回软件，点击 **“解析 AI 输出”** 生成注释结果。

## 许可证

MIT License
