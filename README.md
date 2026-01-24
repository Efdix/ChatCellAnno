# 🧬 ChatCellAnno: 您的通用 AI 单细胞注释助手 (GUI版)

**ChatCellAnno** 是一个轻量级、无需安装、开箱即用的 Windows 桌面程序，旨在成为连接您的 **单细胞分析数据 (Seurat/Scanpy)** 与**任意大语言模型 (LLM)** 之间的桥梁。

无论您使用的是 **GitHub Copilot, DeepSeek, ChatGPT (OpenAI), Claude, 豆包, 元宝** 还是本地部署的模型，ChatCell 都能通过**剪贴板**作为通用接口，协助您根据 Marker 基因快速完成细胞类型注释。

## ✨ 核心特性

*   **🖥️ 图形化界面 (GUI)**: 专为不熟悉代码的用户设计，拖拽文件即可生成 AI 提示词。
*   **🌐 模型无关 (Model Agnostic)**: 不绑定任何特定的 AI 模型。只要它是能聊天的 AI，就能用 ChatCell。
*   **🔒 隐私优先 / Zero-API**: 软件本身不发起任何网络请求。您的数据完全掌握在您手中，通过复制粘贴进行交互，无需配置复杂的 API Key。
*   **⚡ 智能文件识别**: 支持直接读取 Seurat 或 Scanpy 导出的差异基因表格 (`.csv`, `.tsv`, `.txt`)。
*   **📝 两种输出模式**:
    *   **Concise (简洁模式)**: 仅输出细胞类型名称，方便快速浏览。
    *   **Detailed (详细模式)**: 输出推荐 Marker、基因功能解析及原始 Rankings，提供完整的注释证据链。

## 🚀 快速开始 (无需安装)

1.  **下载**: 直接下载 `ChatCellAnno.exe` (通常在 Release 页面或 dist 文件夹中)。
2.  **运行**: 双击打开程序。
3.  **使用步骤**:
    *   **Step 1**: 将您的 Marker 基因表格文件 (CSV/TSV) 拖入窗口，或点击 Browse 选择。
    *   **Step 2**: 设置物种 (如 Human/Mouse) 和组织来源 (如 PBMC, Liver)。选择输出模式 (Concise 或 Detailed)。
    *   **Step 3**: 点击 **"Generate Prompt"** 按钮。
    *   **Step 4**: 此时提示词已自动复制到剪贴板。前往您的 AI 聊天界面 (ChatGPT/Claude/DeepSeek)，按下 `Ctrl+V` 粘贴并发送。
    *   **Step 5**: 阅读 AI 返回的专业的细胞类型注释结果。

## 📄 数据准备指南

ChatCell 需要一份包含 Marker 基因的表格文件。无论您使用 Seurat 还是 Scanpy，只要导出包含 `gene` (基因名) 和 `cluster` (分组) 列的 CSV 文件即可。

**示例文件格式 (Tidy 格式):**
```csv
gene,cluster,avg_log2FC,p_val_adj
CD14,0,2.5,0.0
LYZ,0,2.1,0.0
CD3D,1,3.4,0.0
CD3E,1,3.1,0.0
...
```

**或者宽矩阵格式 (列名=Cluster名):**
```csv
Cluster0,Cluster1
CD14,CD3D
LYZ,CD3E
...
```

## 🛠️ 开发者指南 (源码运行/自行构建)

如果您是开发者并希望修改源码：

1.  **环境配置**:
    ```bash
    # 使用 Mamba/Conda 创建环境
    mamba create -n chatcell python=3.9 -y
    mamba activate chatcell
    
    # 安装依赖
    mamba install pandas pyperclip windnd pyinstaller openpyxl -y
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



