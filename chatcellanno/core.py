"""
核心逻辑模块 (Core Logic Module)

本模块作为 ChatCellAnno 的主要入口点，协调各个子模块（提取器、提示词生成器、解析器）的工作。
它实现了 "Generate" (生成提示词) 和 "Parse" (解析 AI 回复) 这两个核心步骤。
"""

from .extractor import extract_markers_from_file, extract_markers_from_df
from .prompt import generate_annotation_prompt
from .parser import parse_llm_response, generate_annotation_code
# 尝试导入富集分析模块 (Try importing enrichment module)
try:
    from .enrichment import perform_enrichment
    HAS_ENRICHMENT = True
except ImportError:
    HAS_ENRICHMENT = False

import pandas as pd
import os

def annotate_cell_types(
    marker_file: str,
    step: str = "generate",
    response_text: str = None,
    species: str = "Human",
    tissue: str = "PBMC",
    top_n: int = 10,
    mode: str = "concise",
    source: str = "scanpy",
    exclude_types: str = "",
    use_enrichment: bool = False,
    enrichment_db: str = "GO_Biological_Process_2021",
    enrichment_hints: dict = None,
    visual_context: str = None,
    check_expression_file: str = None
):
    """
    ChatCellAnno 标注工作流管理器 (Workflow Controller).
    
    参数 (Parameters):
    - marker_file: Marker 基因文件路径 (.csv/.tsv)
    - step: 当前执行步骤 ("generate" 生成提示词, "parse" 解析回复)
    - response_text: AI 返回的文本内容 (仅在 "parse" 步骤需要)
    - species: 物种名称 (如 Human, Mouse)
    - tissue: 组织类型 (如 PBMC, Brain)
    - top_n: 每个聚类选取的 Top Marker 数量
    - mode: 提示词详细程度 ("concise" 仅类型, "detailed" 包含理由)
    - source: 数据来源格式 ("scanpy" 或 "seurat")
    - exclude_types: 排除的细胞类型关键字
    - use_enrichment: 是否启用功能富集分析
    - enrichment_db: 富集分析数据库路径或名称
    - enrichment_hints: 预生成的富集分析提示信息 (可选)
    - visual_context: 视觉上下文描述 (如 "UMAP plot shows cluster 0 is isolated")
    - check_expression_file: 关键基因表达矩阵文件路径 (可选)
    """
    
    # 检查输入文件是否存在 (仅在 generate 阶段是必须的)
    if step == "generate" and not os.path.exists(marker_file):
         raise FileNotFoundError(f"Marker file not found: {marker_file}")

    # 第一阶段：提取 Marker 并生成 Prompt (Phase 1: Extraction & Generation)
    if step == "generate":
        # 1. 提取用于构造 Prompt 的 Top Markers (Extract top markers for prompt)
        markers_for_prompt = extract_markers_from_file(marker_file, top_n=top_n)

        # 1.1 处理可选的表达矩阵 (Process Optional Expression Matrix)
        # 提供基因表达强度信息可以帮助 AI 区分高度相近的亚群
        expression_matrix_text = None
        if check_expression_file and os.path.exists(check_expression_file):
            try:
                # 读取 CSV/TSV 表达矩阵 (Read matrix data)
                if check_expression_file.endswith('.csv'):
                    df_expr = pd.read_csv(check_expression_file, index_col=0)
                else:
                    df_expr = pd.read_table(check_expression_file, index_col=0)
                
                # 尝试转换为 Markdown 表格格式 (Format to Markdown table)
                try:
                    expression_matrix_text = df_expr.to_markdown()
                except ImportError:
                    # 如果未安装 tabulate，回退到简单的列表格式 (Fallback to TSV if tabulate missing)
                    expression_matrix_text = "Cluster\t" + "\t".join(df_expr.columns) + "\n"
                    for idx, row in df_expr.iterrows():
                        expression_matrix_text += f"{idx}\t" + "\t".join([str(round(v, 4)) for v in row.values]) + "\n"
            except Exception as e:
                print(f"Error reading expression matrix: {e}")
        
        # 1.2 执行功能富集分析 (Perform Pathway Enrichment Analysis)
        full_enrichment_data = None
        if use_enrichment and HAS_ENRICHMENT and not enrichment_hints:
            try:
                # 提取更多的 Marker (如 Top 100) 以提高富集分析的统计效力 (Extract more markers for better enrichment power)
                markers_for_enrich = extract_markers_from_file(marker_file, top_n=100)
                
                # 准备富集分析输入格式 (Prepare enrichment input)
                enrich_input = {k: v.split(", ") for k,v in markers_for_enrich.items()}
                
                # 调用富集模块 (Execute GSEApy/Enrichr logic)
                full_enrichment_data = perform_enrichment(
                    enrich_input, 
                    species=species, 
                    database_path=enrichment_db,
                    top_term_n=3
                )
                
                # 提取富集提示词 (Extract hints for prompt injection)
                enrichment_hints = {k: v['hints'] for k, v in full_enrichment_data.items()}
            except Exception as e:
                print(f"Enrichment Analysis skipped due to error: {e}")
        
        # 1.3 构造多模态提示词 (Construct multimodal prompt for LLM)
        # 整合 Markers, 物种, 组织, 富集信息, 视觉描述和表达矩阵
        prompt = generate_annotation_prompt(
            markers_for_prompt, 
            species, 
            tissue, 
            mode=mode, 
            exclude_types=exclude_types,
            enrichment_hints=enrichment_hints,
            visual_context=visual_context,
            expression_matrix=expression_matrix_text
        )
        return prompt, full_enrichment_data

    # 第二阶段：解析 LLM 的回复 (Phase 2: Parsing LLM Response)
    elif step == "parse":
        # 验证回复文本是否存在 (Check if response exists)
        if not response_text:
            raise ValueError("For step='parse', you must provide 'response_text'.")
        
        # 尝试从原始 Marker 文件获取聚类 ID 列表 (Try to get cluster IDs from original file)
        cluster_names = None
        if marker_file and os.path.exists(marker_file):
            markers = extract_markers_from_file(marker_file, top_n=top_n)
            cluster_names = list(markers.keys())

        # 2.1 调用解析规则处理 AI 文本 (Call parser to extract cell types)
        # annotations: {cluster_id: cell_type}
        # extra_info: {cluster_id: reasoning/other}
        annotations, extra_info = parse_llm_response(response_text, cluster_names)
        
        # 如果初始没有 cluster_names (即用户直接输入回复且未读文件)，则从解析出的字典更新
        # (Update cluster list if it was empty)
        if cluster_names is None:
            cluster_names = list(annotations.keys())
        
        # 2.2 生成配套分析代码 (Generate annotation code snippet)
        # 根据 source (Scanpy/Seurat) 生成对应的 Python 或 R 代码
        code_snippet = generate_annotation_code(annotations, source)
        
        # 2.3 构建格式化的结果表格 (Create formatted result DataFrame)
        results = []
        for cluster in cluster_names:
            cell_type = annotations.get(cluster, "Unknown")
            reasoning = extra_info.get(cluster, "")
            results.append({
                "Cluster": cluster,
                "Cell Type": cell_type,
                "Reasoning": reasoning
            })
        
        df_results = pd.DataFrame(results)
        return df_results, code_snippet
