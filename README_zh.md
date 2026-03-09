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

1. **步骤 1: 数据准备与富集辅助 (Data & Enrichment)**:
   - **差异基因列表 (Differential Gene List)**: 将您的 Marker 基因文件 (.csv, .tsv) 拖入上方指定区域（必填）。
   - **(新!) 辅助富集分析 (Enrichment to Assist AI)**: 勾选此项以通过本地或在线 (Enrichr) 工具进行功能富集。富集结果将被客观地写入提示词中，帮助 AI 减少幻觉并提供生物学证据。
   - **UMAP/t-SNE 截图 (Optional)**: 从您的结果中截图可视化图表，点击框体或 Ctrl+V 粘贴。这能为 AI 提供重要的空间聚类上下文。

2. **步骤 2: 提示词配置 (Prompt Configuration)**:
   - 设置物种、组织类型及输出格式（简洁/详细）。
   - 点击 **“生成并复制提示词”**。如果添加了图片，软件会提示您同时复制图像。

3. **步骤 3: AI 交互 (AI Interaction)**:
   - 将提示词和截图一并发送给您的大模型对话框 (如 ChatGPT, DeepSeek, Claude)。
   - 程序内置了浏览器面板，方便您直接在软件内完成对话。

4. **步骤 4: 解析结果 (Parse AI Output)**:
   - 复制 AI 返回的 Markdown 表格，粘贴回软件的解析区域，点击 **“解析 AI 输出”**。
   - 软件将生成可用于 Scanpy (Python) 或 Seurat (R) 的代码片段，直接完成数据集注释。

## 界面管理

- **可定制布局**: 右侧的功能面板（功能富集、内置浏览器、基因组工具）支持**拖拽换位**、**关闭**。
- **恢复面板**: 如果误关了面板，可通过顶部齿轮图标的 **“显示面板 (Show Panels)”** 菜单重新打开。

## 许可证

MIT License
