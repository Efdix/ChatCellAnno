# 🧬 ChatCell: 您的通用 AI 单细胞注释助手

**ChatCell** 是一个轻量级的 Python 库，旨在成为连接您的 **Scanpy** 数据分析与**任意大语言模型 (LLM)** 之间的桥梁。

无论您使用的是 **GitHub Copilot, DeepSeek, ChatGPT (OpenAI), Claude, 豆包, 元宝** 还是本地部署的模型，ChatCell 都能通过**剪贴板**作为通用接口，协助您根据 Marker 基因快速完成细胞类型注释。

## ✨ 核心特性

*   **🌐 模型无关 (Model Agnostic)**: 不绑定任何特定的 AI 模型。只要它是能聊天的 AI，就能用 ChatCell。
*   **🔒 隐私优先 / 无需 API**: 库本身不发起任何网络请求。您的数据完全掌握在您手中，通过复制粘贴进行交互，无需配置复杂的 API Key。
*   **⚡ 双模式支持**: 
    *   **Scanpy 原生**: 无缝对接 `anndata` 对象和 `rank_genes_groups` 结果。
    *   **通用文件模式**: 支持直接读取 CSV/TSV 格式的 Marker 表进行分析。
*   **📋 极简交互**: 生成优化过的 Prompt 自动复制及结果自动解析，让手动流程如自动化般顺滑。

## 🛠️ 安装

pip 安装较慢时，推荐先用 conda/mamba 安装核心依赖：

```bash
# 1. 创建环境并安装依赖 (速度更快)
conda create -n chatcell -c conda-forge python=3.9 scanpy pandas anndata pyperclip -y
# 或者: mamba create -n chatcell -c conda-forge python=3.9 scanpy pandas anndata pyperclip -y

conda activate chatcell

# 2. 安装 ChatCell
git clone https://github.com/Efdix/ChatCell.git
cd ChatCell
pip install -e .
```

## 📖 使用教程

ChatCell 提供两种使用方式：直接基于 Scanpy 对象 (`adata`) 或基于 Marker 基因表格文件。

### 方式一：Scanpy 对象工作流 (`anndata`)

**1. 准备数据**
确保您已经运行了差异表达分析：
```python
import scanpy as sc
import chatcell

# ... 加载并处理您的 adata ...
# 计算 Marker 基因 (关键步骤!)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
```

**2. 生成 Prompt (Step 1)**
此命令提取 Marker 并将优化后的 Prompt 复制到剪贴板。
```python
chatcell.annotate_cell_types(
    adata=adata,        # 传入 AnnData 对象
    step="generate", 
    species="Human", 
    tissue="PBMC",
    mode="concise"      # 模式: concise, evidence, recommendation
)
```
*输出: `✅ Prompt has been COPIED to your clipboard!`*

**3. 与 AI 对话**
将 Prompt 粘贴给 AI（Copilot, DeepSeek, ChatGPT...）。AI 会回复细胞类型列表。

**4. 应用注释 (Step 2)**
复制 AI 的回复内容，传回 ChatCell 更新 `adata`。
```python
ai_response = """
CD4+ Naive T
CD14+ Monocyte
...
"""

chatcell.annotate_cell_types(
    adata=adata,
    step="parse", 
    response_text=ai_response
)

# 结果已自动写入: data.obs['chatcell_annotation']
sc.pl.umap(adata, color='chatcell_annotation')
```

---

### 方式二：通用文件工作流 (CSV/TSV)

如果您没有 `adata` 对象，只有一个包含 Marker 基因的表格（列名=簇名，列值=基因列表），也可以使用 ChatCell。

**1. 准备文件 (markers.tsv)**
```csv
Cluster0    Cluster1
CD14        CD3D
LYZ         CD3E
...         ...
```

**2. 生成 Prompt (Step 1)**
```python
import chatcell

chatcell.annotate_cell_types(
    marker_file="markers.tsv",  # 传入文件路径
    step="generate",
    species="Mouse",
    tissue="Brain"
)
```

**3. 解析结果 (Step 2)**
```python
ai_response = "..." # 从 AI 处复制

annotations, extra_info = chatcell.annotate_cell_types(
    marker_file="markers.tsv",
    step="parse",
    response_text=ai_response
)

print(annotations)
# Output: {'Cluster0': 'Microglia', 'Cluster1': 'T Cell'}
```

## 🧠 高级模式 (Prompt Engineering)

通过 `mode` 参数，您可以获得更丰富的分析结果，Prompt 已经过针对性优化：

*   **`mode="concise"` (默认)**: 
    *   仅获取标准细胞类型名称 (Cell Ontology)。
    *   适合后续直接用于自动标注。
*   **`mode="evidence"`**: 
    *   让 AI 列出支持该判断的 Marker 基因证据。
    *   格式: `Cell Type | Supported by: gene1, gene2`
*   **`mode="recommendation"`**: 
    *   让 AI 推荐列表中缺失但有助于确认身份的 Marker 基因。
    *   格式: `Cell Type | Recommended Markers: geneA, geneB`

## 📂 示例代码

在 `examples/` 目录下有完整的运行示例。无需 `cd` 进入文件夹，直接在项目根目录运行即可：

*   生成测试数据: `python examples/gen_data.py`
*   Scanpy 流程: `python examples/run_with_adata.py`
*   TSV 流程: `python examples/run_with_tsv.py`

## 📦 依赖要求

*   pandas
*   pyperclip
*   anndata
*   scanpy (可选，仅用于方式一)



