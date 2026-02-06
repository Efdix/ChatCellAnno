"""
核心逻辑模块 (Core Logic Module)

本模块作为 ChatCellAnno 的主要入口点，协调各个子模块（提取器、提示词生成器、解析器）的工作。
它实现了 "Generate" (生成提示词) 和 "Parse" (解析 AI 回复) 这两个核心步骤。
"""

from .extractor import extract_markers_from_file, extract_markers_from_df
from .prompt import generate_annotation_prompt
from .parser import parse_llm_response, generate_annotation_code
# Try importing enrichment module
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
    enrichment_hints: dict = None
):
    """
    ChatCellAnno 标注工作流管理器。
    """
    
    # 检查输入文件是否存在 (仅在 generate 阶段是必须的)
    if step == "generate" and not os.path.exists(marker_file):
         raise FileNotFoundError(f"Marker file not found: {marker_file}")

    # 第一步：提取 Marker 并生成 Prompt (Extract Markers)
    if step == "generate":
        # 1. 提取用于 Prompt 的 Markers
        markers_for_prompt = extract_markers_from_file(marker_file, top_n=top_n)
        
        full_enrichment_data = None
        if use_enrichment and HAS_ENRICHMENT and not enrichment_hints:
            try:
                # 2. 提取用于 Enrich 的 Markers (top 100)
                markers_for_enrich = extract_markers_from_file(marker_file, top_n=100)
                
                enrich_input = {k: v.split(", ") for k,v in markers_for_enrich.items()}
                
                full_enrichment_data = perform_enrichment(
                    enrich_input, 
                    species=species, 
                    database_path=enrichment_db,
                    top_term_n=3
                )
                
                enrichment_hints = {k: v['hints'] for k, v in full_enrichment_data.items()}
            except Exception as e:
                print(f"Enrichment Analysis skipped due to error: {e}")
        
        # 根据提取的数据生成给 LLM 的提示词
        prompt = generate_annotation_prompt(
            markers_for_prompt, 
            species, 
            tissue, 
            mode=mode, 
            exclude_types=exclude_types,
            enrichment_hints=enrichment_hints
        )
        return prompt, full_enrichment_data

    # 第二步：解析 LLM 的回复 (Parse LLM Response)
    elif step == "parse":
        if not response_text:
            raise ValueError("For step='parse', you must provide 'response_text'.")
        
        # 尝试获取 cluster names (如果文件存在)
        cluster_names = None
        if marker_file and os.path.exists(marker_file):
            markers = extract_markers_from_file(marker_file, top_n=top_n)
            cluster_names = list(markers.keys())

        # 调用解析器处理 LLM 返回的文本
        # 如果 cluster_names 为 None，解析器将尝试从表格中提取 ID
        annotations, extra_info = parse_llm_response(response_text, cluster_names)
        
        # 如果没有 cluster_names (即没读文件且从 response 提取)，则从 annotations 更新
        if cluster_names is None:
            cluster_names = list(annotations.keys())
        
        # 生成代码片段
        code_snippet = generate_annotation_code(annotations, source)
        
        # 创建结果 DataFrame (Create a result DataFrame)
        results = []
        for cluster in cluster_names:
            row = {
                'Cluster': cluster,
                'Annotation': annotations.get(cluster, 'Unknown'),  # 获取标注结果，如果解析失败默认为 Unknown
            }
            # 如果存在额外信息 (detailed 模式)，也添加到结果中
            if extra_info:
                row['Extra_Info'] = extra_info.get(cluster, '')
            results.append(row)
            
        result_df = pd.DataFrame(results)
        print("✅ Annotations Parsed Successfully. (解析成功)")
        
        # 返回 DataFrame 和 代码
        return result_df, code_snippet
    
    else:
        raise ValueError("Invalid step. Use 'generate' or 'parse'.")


