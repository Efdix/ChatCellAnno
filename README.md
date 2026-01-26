# 🧬 ChatCellAnno: 通用 AI 单细胞注释助手 (GUI版)

**ChatCellAnno** 是一个轻量级、无需安装、开箱即用的 Windows 桌面程序，旨在成为连接您的 **单细胞分析数据 (Seurat/Scanpy)** 与**任意大语言模型 (LLM)** 之间的桥梁。

无论您使用的是 **GitHub Copilot, DeepSeek, ChatGPT (OpenAI), Claude, 豆包, 元宝** 还是本地部署的模型，ChatCell 都能通过**剪贴板**作为通用接口，协助您根据 Marker 基因快速完成细胞类型注释。

## ✨ 核心特性

*   **🖥️ 图形化界面 (GUI)**: 专为不熟悉代码的用户设计，拖拽文件即可生成 AI 提示词。
*   **🌐 模型无关 (Model Agnostic)**: 不绑定任何特定的 AI 模型。只要它是能聊天的 AI，就能用 ChatCellAnno。
*   **🔒 隐私优先 / Zero-API**: 软件本身不发起任何网络请求。您的数据完全掌握在您手中，通过复制粘贴进行交互，无需配置复杂的 API Key。
*   **⚡ 智能文件识别**: 
    *   **Scanpy (Python)**: 支持 `sc.get.rank_genes_groups_df` 导出的表格 (`names`, `group`)。
    *   **Seurat (R)**: 支持 `FindAllMarkers` 导出的表格 (`gene`, `cluster`)。
*   **📝 全流程支持**:
    *   **生成 Prompt**: 一键生成高质量的细胞注释提示词。
    *   **生成代码**: 将 AI 的回答复制回软件，自动生成可执行的 Python (Scanpy) 或 R (Seurat) 代码，直接用于重命名聚类。

## 🚀 快速开始 (无需安装)

1.  **下载**: 直接下载 `ChatCellAnno.exe` (在 Release 页面或 dist 文件夹中)。
2.  **运行**: 双击打开程序。
3.  **使用步骤**:
    *   **Step 1**: 将您的 Marker 表格文件拖入。选择您的数据来源 (**Scanpy** 或 **Seurat**)。
    *   **Step 2**: 设置物种 (如 Human) 和组织 (如 PBMC)。
    *   **Step 3**: 选择输出模式 (Concise 或 Detailed)。
    *   **Step 4**: 点击 **"Generate Prompt"**。提示词已自动复制。
    *   **Step 5**: 前往 AI 聊天界面 (ChatGPT/Claude 等)，粘贴并发送。
    *   **Step 6**: **复制 AI 的回答**，粘贴回软件的 "Step 4: Parse AI Response" 区域。
    *   **Step 7**: 点击 **"Process & Generate Code"**。软件会自动生成对应的 Python 或 R 代码，助您一键完成注释。

## 📄 数据准备指南

ChatCellAnno 严格遵循 Scanpy 和 Seurat 的标准输出格式。

**1. Scanpy 用户 (Python):**
使用 `sc.get.rank_genes_groups_df(adata, None)` 导出的 CSV。
必须包含列: `names` (基因名), `group` 或 `cluster` (聚类号)。

**2. Seurat 用户 (R):**
使用 `FindAllMarkers(seurat_obj)` 导出的 CSV。
必须包含列: `gene` (基因名), `cluster` (聚类号)。

*(注: 不再支持旧版的宽矩阵格式)*

## 🛠️ 开发者指南 (源码运行/自行构建)

如果您是开发者并希望修改源码：

1.  **环境配置**:
    ```bash
    # 使用 Mamba/Conda 创建环境
    mamba create -n chatcellanno python=3.9 -y
    mamba activate chatcellanno
    
    # 安装依赖
    mamba install pandas pyinstaller openpyxl -y
    pip install windnd
    pip install pyperclip
    ```

2.  **运行 GUI**:
    ```bash
    python gui.py
    ```

3.  **构建 EXE**:
    ```bash
    ./build.ps1
    # 或者
    pyinstaller --noconfirm --onefile --windowed --name "ChatCellAnno" --hidden-import "pandas" --hidden-import "pyperclip" --hidden-import "windnd" "gui.py"
    ```

## 📜 许可证

MIT License



