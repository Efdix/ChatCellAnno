# 🧬 ChatCell: 您的通用 AI 单细胞注释助手

**ChatCell** 是一个轻量级的 Python 库，旨在成为连接您的 **Scanpy** 数据分析与**任意大语言模型 (LLM)** 之间的桥梁。

无论您使用的是 **GitHub Copilot, DeepSeek, ChatGPT (OpenAI), Claude, 豆包, 元宝** 还是本地部署的模型，ChatCell 都能通过**剪贴板**作为通用接口，协助您根据 Marker 基因快速完成细胞类型注释。

## ✨ 核心特性

*   **🌐 模型无关 (Model Agnostic)**: 不绑定任何特定的 AI 模型。只要它是能聊天的 AI，就能用 ChatCell。
*   **🔒 隐私优先 / 无需 API**: 库本身不发起任何网络请求。您的数据完全掌握在您手中，通过复制粘贴进行交互，无需配置复杂的 API Key。
*   **⚡ Scanpy 原生**: 无缝对接 `anndata` 对象和 `rank_genes_groups` 结果。
*   **📋 极简交互**: 生成优化过的 Prompt 自动复制及结果自动解析，让手动流程如自动化般顺滑。

## 🛠️ 安装

```bash
git clone https://github.com/Efdix/ChatCell.git
cd ChatCell
pip install -e .
```

## 🚀 快速开始

### 1. 准备工作
确保您已经用 Scanpy 处理了数据并运行了差异表达分析：

```python
import scanpy as sc
import chatcell

# ... 加载并处理您的 adata ...
# 计算 Marker 基因 (关键步骤!)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
```

### 2. 生成 Prompt (步骤 1)
运行此命令提取 Marker 并生成适合 LLM 的 Prompt。它会自动将 Prompt 复制到您的剪贴板。

```python
# 返回 prompt 字符串并复制到剪贴板
chatcell.annotate_cell_types(
    adata, 
    step="generate", 
    species="Human", 
    tissue="PBMC",
    mode="concise" # 模式: concise (简洁), evidence (含证据), recommendation (智能推荐)
)
```

**输出:**
```text
✅ Prompt has been COPIED to your clipboard!
👉 Go to your AI Chat (Copilot, DeepSeek, ChatGPT) and press Ctrl+V (Paste), then Enter.
```

### 3. 与 AI 对话
将 Prompt 粘贴到您喜欢的 AI 聊天界面中（例如 DeepSeek 网页版、VS Code Copilot、ChatGPT）。AI 会回复一个细胞类型列表。

**AI 回复示例 (DeepSeek/Copilot):**
```text
CD4+ Naive T
CD14+ Monocyte
B cell
NK cell
...
```

### 4. 应用注释 (步骤 2)
复制 AI 的回复内容，并传回 ChatCell 以更新您的 `adata`。

```python
ai_response = """
CD4+ Naive T
CD14+ Monocyte
B cell
NK cell
"""

# 解析文本并添加 'chatcell_annotation' 列到 adata.obs
chatcell.annotate_cell_types(
    adata, 
    step="parse", 
    response_text=ai_response
)

# 可视化
sc.pl.umap(adata, color='chatcell_annotation')
```

## 🧠 高级模式

通过 `mode` 参数，您可以获得更丰富的分析结果：

*   **`mode="concise"` (默认)**: 仅获取细胞类型名称。
*   **`mode="evidence"`**: 让 AI 列出支持该判断的 Marker 基因证据。
    *   结果存入: `adata.obs['chatcell_extra_info']`
*   **`mode="recommendation"`**: 让 AI 推荐列表中缺失但有助于确认身份的 Marker 基因。
    *   结果存入: `adata.obs['chatcell_extra_info']`

## 📦 依赖要求
*   Python >= 3.8
*   scanpy
*   anndata
*   pyperclip
*   pandas

## 📄 许可证
MIT



