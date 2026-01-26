# ChatCellAnno 架构与源码学习指南

本文档旨在帮助开发者了解 ChatCellAnno 的内部架构、设计理念以及各模块的具体实现。
如果你想学习 "如何编写一个不需要 API Key 就能利用大模型能力的 Python 工具"，这份文档就是为你准备的。

## 1. 项目概述 (Overview)

**ChatCellAnno** 是一个用于单细胞测序数据注释的工具。它不需要调用 OpenAI 或 Anthropic 的付费 API，而是利用系统的 **剪贴板 (Clipboard)** 作为数据传输的桥梁，连接本地 Python 程序和网页版大模型 (如 ChatGPT, Copilot, DeepSeek)。

### 核心理念 (Core Philosophy)
1.  **Zero-API (无 API)**: 不依赖任何第三方服务的 API Key，完全免费使用。
2.  **Clipboard as Transport (剪贴板即传输)**: 数据流向是 `本地程序 -> 剪贴板 -> 用户粘贴给 AI -> AI 回复 -> 用户复制 -> 剪贴板 -> 本地程序`。
3.  **Human-in-the-loop (人机交互)**: 强调人工审核，而不是全自动黑盒。
4.  **Actionable Code (代码生成)**: 不仅仅给出一个名字，而是直接生成可以运行的 Python/R 代码，真正打通分析流程的最后一公里。

---

## 2. 目录结构说明 (Directory Structure)

```text
ChatCellAnno/
├── chatcellanno/           # 核心逻辑包
│   ├── __init__.py         # 包初始化
│   ├── core.py             # [核心] 总控制器，增加 source 参数和代码生成返回
│   ├── extractor.py        # [数据] 区分处理 Scanpy/Seurat 格式，移除 Matrix 格式支持
│   ├── prompt.py           # [生成] 负责构建发给 AI 的提示词 (Prompt Engineering)
│   ├── parser.py           # [解析] 负责解析回复并生成 Python/R 代码
├── examples/               # 示例数据
│   └── example.csv
├── gui.py                  # [界面] 新增 Step 4: Parse & Generate Code 区域
├── build.ps1               # 打包脚本 (PyInstaller)
├── setup.py                # 安装脚本
└── README.md
```

---

## 3. 模块详解与学习顺序 (Module Breakdown)

建议按照以下顺序阅读源码，以便循序渐进地理解：

### 第一步：理解数据入口 `extractor.py`
这个模块展示了如何处理生物信息学中常见的 CSV 数据。
- **关键点**: **数据源区分**。代码现在有了 `source='scanpy'` 和 `source='seurat'` 分支。
- **变更**: 不再支持猜测列名，而是严格遵循 Scanpy (`names`, `group`) 和 Seurat (`gene`, `cluster`) 的标准字段名。这是为了保证后续代码生成的准确性。

### 第二步：理解提示词工程 `prompt.py`
这是本项目最 "AI" 的部分。
- **关键点**: 如何构造一个能让 AI 稳定输出格式化数据的 Prompt。
- **技巧**: 注意 `generate_annotation_prompt` 函数中，我们明确要求 AI "Return exactly X lines" (返回确定的行数)，并给出 Example。这是一种 **Few-Shot Prompting (少样本提示)** 技巧。

### 第三步：理解控制中心 `core.py`
它是连接各个器官的大脑。
- **关键点**: `annotate_cell_types` 函数现在不仅要传递 `step`，还要传递 `source`。
- **设计模式**: 这是一个典型的 "Facade Pattern" (外观模式)。它现在在 `step='parse'` 时返回两个值：`Dataframe` (结果表) 和 `str` (生成代码)。

### 第四步：理解结果解析 `parser.py`
当 AI 把结果通过剪贴板还给你时，如何处理？
- **关键点**: `parse_llm_response` (解析) 和 `generate_annotation_code` (代码生成)。
- **代码生成**: 这是一个新功能。根据用户选择的 `source`，生成 Python 字典映射代码 (Scanpy) 或 R 语言向量命名代码 (Seurat)。这体现了 "Output as Code" 的思想。
- (**对齐逻辑**: 同前，依赖行号一致性)

### 第五步：理解用户界面 `gui.py`
最后看 GUI，了解如何将上述逻辑封装成普通用户可用的软件。
- **界面变更**: 
    - 增加了 "Source Type" 单选框。
    - 界面高度增加，新增了 "Step 4" 区域，允许用户粘贴 AI 回复并获取代码。
- **交互逻辑**: 用户在 Step 4 点击按钮 -> 调用 core -> core 调用 parser -> 返回代码 -> 显示在 GUI 文本框。

---

## 4. 关键数据流 (Data Flow)

### 流程一：生成 Prompt (Step 1: Generate)
1.  **用户** 在 GUI 选择 CSV 文件，**并选择数据来源 (Scanpy/Seurat)**。
2.  `gui.py` 调用 `core.annotate_cell_types(step='generate', source='...')`。
3.  `core` 调用 `extractor.extract_markers_from_file` -> 严格检查列名并提取。
4.  `core` 调用 `prompt.generate_annotation_prompt` -> 格式化为 Prompt。
5.  `prompt` 模块调用 `pyperclip.copy()` 将文本写入剪贴板。

### 流程二：解析与代码生成 (Step 2: Parse) (GUI 已支持)
1.  **用户** 将 ChatGPT 的回复复制，粘贴到 GUI 的 **Step 4** 文本框。
2.  用户点击 "Process & Generate Code"。
3.  `gui.py` 调用 `core.annotate_cell_types(step='parse', response_text='...', source='...')`。
4.  `core` **再次**调用 `extractor` 读取原文件 -> 获取 Cluster 列表（作为校验基准）。
5.  `core` 调用 `parser.parse_llm_response` -> 得到 `{ "0": "T Cell" }`。
6.  `core` 调用 `parser.generate_annotation_code` -> 生成 `adata.obs['cell_type'] = ...` 或 `seurat_obj <- RenameIdents...`。
7.  `gui.py` 将生成的代码显示在界面下方供用户复制。

---

## 5. 开发者建议 (Tips for Developers)

1.  **调试技巧**: 
    - 不要每次都运行 GUI。可以在项目根目录下创建一个 `debug.py`，直接调用 `chatcellanno.core` 中的函数进行调试。
    - 模拟 AI 回复：在开发 Parser 时，不需要真的去问 AI，可以手动写一个字符串变量来模拟 AI 的输出，测试解析器的健壮性。

2.  **扩展方向**:
    - 目前只支持 CSV/TSV，可以尝试添加 `.h5ad` (Anndata) 支持。
    - Prompt 目前是英文的，可以尝试修改 `prompt.py` 支持中文 Prompt。

3.  **Python 特性**:
    - 本项目大量使用了 Type Hints (类型提示)，如 `def func(a: str) -> bool:`。这对于大型项目维护非常有帮助。

希望这份文档能帮助你快速上手并理解 ChatCellAnno 的代码逻辑！
