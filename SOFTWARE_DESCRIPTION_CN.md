# 软件说明书 (Software Description) for ChatCellAnno

## 1. 软件基本信息
- **软件名称**: ChatCellAnno 单细胞智能注释辅助系统
- **软件版本**: V1.0.0
- **开发语言**: Python 3.10+
- **运行平台**: Windows 10/11 (主要), macOS, Linux
- **依赖框架**: PySide6 (GUI), Pandas (数据处理), Pyperclip (系统剪贴板交互)

## 2. 软件概述
ChatCellAnno 是一款专为单细胞转录组测序（scRNA-seq）数据分析设计的辅助软件。它采用“人在回路”（Human-in-the-loop）的设计理念，旨在填补传统的单细胞分析流程（如 Scanpy, Seurat）与现代大语言模型（LLM，如 ChatGPT, DeepSeek, Claude）之间的鸿沟。

本软件首创“零API（Zero-API）”架构，无需用户配置任何 API Key，也无需联网上传敏感数据，仅通过系统剪贴板即可实现与任意 LLM 的无缝交互，极大降低了科研人员使用 AI 进行细胞类型注释的门槛，同时确保了数据的绝对隐私安全。

## 3. 运行环境要求
### 3.1 硬件环境
- **CPU**: Intel Core i5 或同级以上处理器
- **内存**: 8GB RAM 以上（建议 16GB 以处理大型单细胞数据）
- **硬盘**: 至少 500MB 可用空间
- **显示器**: 分辨率 1920x1080 及以上

### 3.2 软件环境
- **操作系统**: Windows 10 / 11 64位
- **Python环境**: Python 3.8 及以上版本
- **浏览器**: 系统需安装 Edge 或 Chrome 浏览器（用于内部 WebEngine 渲染）

## 4. 软件主要功能
### 4.1 标记基因提取与解析
- 支持导入 Scanpy (.csv) 和 Seurat (.tsv/.txt) 格式的 Marker 基因差异分析结果表。
- 自动识别 Marker 表中的关键列（如 Gene Name, Cluster ID, LogFoldChange 等）。
- 支持用户自定义筛选 Top N 个高变基因用于后续分析。

### 4.2 功能富集分析 (Functional Enrichment)
- **在线模式**: 自动调用 Enrichr API 进行 GO (Gene Ontology) 和 KEGG 通路富集分析。
- **离线模式**: 支持导入本地 GMT 格式数据库（如 CellMarker 2.0, PanglaoDB），在无网环境下进行超几何分布检验。
- **辅助决策**: 将富集结果作为“提示线索（Hints）”注入到 AI Prompt 中，有效减少 AI 幻觉（Hallucination）。

### 4.3 智能提示词工程 (Prompt Engineering)
- 内置针对生物信息学优化的 Prompt 模板。
- **视觉上下文 (Visual Context)**: [新] 支持用户粘贴 UMAP 或 t-SNE 降维聚类图的截图。软件会自动提取图像数据，并引导用户将其作为多模态输入（Multi-modal Input）提供给大模型，辅助 AI 结合空间分布模式进行细胞类型推断。
- 支持“简洁模式”与“详细推理模式”切换。
- 自动构建包含角色设定、任务描述、数据约束的结构化提示词，一键复制到剪贴板。

### 4.4 跨物种基因组比对 (Genome Alignment)
- 集成多物种同源基因搜索模块，支持 Human, Mouse 及其他非模式生物（如 Flying Squirrel）。
- **可视化比对**: 提供 MEGA 风格的序列比对视图，支持彩色氨基酸显示和位点标尺。
- **专业算法**: 内置 Biopython 的星状比对（Star Alignment）逻辑，支持 Gap 处理和保守性分析。
- **导出功能**: 支持将比对后的多物种序列导出为标准 FASTA 格式。

### 4.5 代码生成与导出
- 解析 AI 返回的 Markdown 表格，自动映射回原始 Cluster ID。
- 根据用户选择的下游平台（Scanpy 或 Seurat），自动生成可直接运行的 Python 或 R 注释代码。

## 5. 使用流程
1.  **数据导入**: 将差异基因文件拖入“Step 1”区域，软件自动解析。
2.  **富集分析 (可选)**: 在“Step 2”选择数据库进行预分析，获取功能线索。
3.  **生成提示词**: 在“Step 3”配置物种和组织信息，点击生成 Prompt 并自动复制。
4.  **AI 交互**: 在右侧内置浏览器打开任意 AI 助手（如 DeepSeek），粘贴 Prompt 并发送。
5.  **结果解析**: 复制 AI 返回的表格内容，粘贴回“Step 4”输入框。
6.  **代码导出**: 点击处理按钮，在“Step 5”获取最终注释代码。

## 6. 技术特点与创新点
1.  **隐私优先**: 彻底摒弃 API 调用，数据交互完全本地化/剪贴板化，杜绝测序数据泄露风险。
2.  **模型无关性 (Model Agnostic)**: 不绑定特定大模型，用户可随意切换 GPT-4, Claude 3.5, DeepSeek-V3 等最新模型。
3.  **多模态验证**: 结合传统的统计学富集分析与 AI 语义推理，显著提高了细胞类型注释的准确度。
4.  **交互式基因组学**: 在单细胞注释流程中无缝集成了基因组序列比对功能，支持从转录组到基因组的跨维度分析。
