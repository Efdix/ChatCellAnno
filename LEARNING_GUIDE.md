# 🎓 ChatCellAnno 源码开发学习指南

这份指南将带你从零开始，解构 **ChatCellAnno** 是如何开发出来的。无论你是生物信息学研究者想学 Python 开发，还是 Python 开发者想了解如何利用 LLM 赋能传统软件，这份教程都适合你。

---

## 🧭 学习路线图 (Roadmap)

我们将开发过程拆解为 6 个阶段，这也是 ChatCellAnno 实际的演进路径：

1.  **理念验证 (MVP)**: 实现 "Python <-> 剪贴板 <-> AI" 的通路。
2.  **数据层 (Data)**: 处理生信领域的 Pandas 数据框 (Scanpy/Seurat)。
3.  **提示词工程 (Prompt)**: 学习如何"控制" AI 的输出格式。
4.  **业务逻辑 (Core)**: 串联提取、生成、解析全流程。
5.  **界面开发 (GUI)**: 使用 Tkinter 打造图形界面。
6.  **工程化 (Build)**: 使用 PyInstaller 打包成无需安装的 EXE。

---

## Phase 1: 核心机制 - 剪贴板交互

**目标**: 不申请 OpenAI API Key，也能让 Python 程序和 AI 聊天。

**核心库**: `pyperclip`

**原理**:
操作系统有一个共享的内存区域叫"剪贴板"。
- 当你按下 `Ctrl+C`，系统把文本存入剪贴板。
- 当你按下 `Ctrl+V`，系统从剪贴板读出文本。
Python 可以通过 `pyperclip` 模拟这两个动作。

**代码原型**:
```python
import pyperclip
import time

# 1. 发送数据给 AI
prompt = "请解释什么是单细胞测序？"
pyperclip.copy(prompt)
print("已复制到剪贴板，请去网页版 AI 粘贴并提问。")

# 2. 等待用户把 AI 的回答复制回来
input("等你把 AI 的回答复制后，按回车键...")
response = pyperclip.paste()
print("Python 收到了 AI 的回复:", response)
```
> **思考**: 这就是 ChatCellAnno 的 "Zero-API" 核心。虽然土，但极其有效且免费。

---

## Phase 2: 数据处理层 (Extractor)

**目标**: 从复杂的生信分析结果 (CSV) 中提取出 AI 需要的信息 (Marker 基因)。

**核心库**: `pandas`

**挑战**: Scanpy 和 Seurat 导出的表格格式不同。
- **Scanpy**: 列名通常是 `names` (基因) 和 `group` (聚类)。
- **Seurat**: 列名通常是 `gene` (基因) 和 `cluster` (聚类)。

**实现逻辑 (参见 `chatcellanno/extractor.py`)**:
我们需要编写一个能够"自适应"的函数。
```python
def extract_markers(df, source):
    if source == 'scanpy':
        # 找 names 和 group 列
        ...
    elif source == 'seurat':
        # 找 gene 和 cluster 列
        ...
    # 无论来源如何，最后都统一转成字典格式:
    # { '0': 'CD3D, CD8A', '1': 'CD19, MS4A1' }
```

---

## Phase 3: 提示词工程 (Prompt Engineering)

**目标**: 让 AI 不仅能回答，还要按照我们**指定的格式**回答，方便程序解析。

**核心文件**: `chatcellanno/prompt.py`

**关键技巧**:
1.  **Role Play (角色扮演)**: "你是一个专家..."
2.  **Few-Shot (少样本提示)**: 给它一个 Example。
    > *不要只说 "返回格式化数据"，要给它看 "Cluster0: T Cell | CD3D"*
3.  **Constraint (强约束)**: "Return exactly one line per cluster" (每一类只返回一行)。

**Prompt 模板示例**:
```text
Identify cell types...

IMPORTANT: Return exactly 5 lines.
Format: ClusterID: CellType

Example:
Cluster0: T Cell
Cluster1: B Cell

Task Data:
Cluster0: CD3D, ...
```

---

## Phase 4: 解析器开发 (Parser)

**目标**: 把 AI 返回的自然语言文本，变回 Python 的字典或代码。

**核心文件**: `chatcellanno/parser.py`

**难点**: AI 即使被约束了，有时也会罗嗦。比如它可能会输出：
> "Sure! Here is the annotation:
> ```
> Cluster0: T Cell
> ```"

**解决方案**:
1.  **清洗**: 去掉 "Sure!", "Here is...", Markdown 标记 (```)。
2.  **映射**: 因为我们在 Prompt 里要求了"顺序一致"，所以我们可以简单粗暴地：
    - 第 1 行有效文本 -> 对应 Cluster 0
    - 第 2 行有效文本 -> 对应 Cluster 1

---

## Phase 5: GUI 开发 (Tkinter)

**目标**: 给 Python 脚本穿上一层"衣服"，让不懂代码的人也能点点点。

**核心文件**: `gui.py`

**为什么选 Tkinter?**
- Python 内置，无需安装额外库 (PyQt 需要 pip install，打包体积大)。
- 足够简单，适合这种工具类软件。

**架构模式**:
- **UI 线程**: 负责显示按钮、响应点击。
- **Core 逻辑**: 按钮点击后调用 `annotate_cell_types`。

**关键控件**:
- `tk.Button`: 按钮。
- `ttk.Entry` / `tk.StringVar`: 输入框和变量绑定。
- `windnd` (第三方): 实现文件拖拽进窗口的功能 (Windows 特有黑科技)。

---

## Phase 6: 打包发布 (PyInstaller)

**目标**: 把 `.py` 文件变成 `.exe` 文件，这样别人的电脑上没有 Python 也能运行。

**构建脚本**: `build.ps1`

**PyInstaller 重要参数**:
- `--onefile`: 把所有依赖 (Pandas, Numpy 等) 压进一个 exe 文件。
- `--windowed`: 运行时不显示黑色的命令行窗口。
- `--hidden-import`: 显式告诉打包器要把 `pandas` 等库打进去 (有时候自动分析会漏掉)。

---

## 🎯 系统架构总结

当你阅读源码时，请脑补这个流程图：

1. **GUI (gui.py)**: 用户拖入 CSV -> 获取文件路径。
2. **Core (core.py)**: 读取文件 -> 调用 `extractor.py` 提取 Marker 字典。
3. **Core**: 将 Marker 字典交给 `prompt.py` -> 生成 Prompt 字符串。
4. **GUI**: 提示用户 Prompt 已复制。
   --- (用户去和 AI 聊天，把结果复制回来) ---
5. **GUI**: 用户点击 "Generate Code" -> 把剪贴板内容传给 `core.py`。
6. **Core**: 调用 `parser.py` 清洗文本 -> 生成 Python/R 代码。
7. **GUI**: 显示最终代码。

## 💡 动手实践建议

1. **修改 Prompt**: 尝试在 `chatcellanno/prompt.py` 里修改 `detailed` 模式的模板，看看 AI 的回答会有什么变化。
2. **支持中文**: 现在的 Prompt 是英文的，试着把它改成中文指令，让 AI 用中文解释细胞类型。
3. **Debug GUI**: 在 `gui.py` 里加一行 `print("Button Clicked")`，运行代码，看看控制台的输出，理解事件驱动编程。

祝你学习愉快！🚀
