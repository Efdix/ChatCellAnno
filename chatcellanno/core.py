"""
核心逻辑模块 (Core Logic Module)

本模块作为 ChatCell 的主要入口点，协调各个子模块（提取器、提示词生成器、解析器）的工作。
它实现了 "Generate" (生成提示词) 和 "Parse" (解析 AI 回复) 这两个核心步骤。
"""

from .extractor import extract_markers_from_file
from .prompt import generate_annotation_prompt
from .parser import parse_llm_response, generate_annotation_code
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
    source: str = "scanpy"
):
    """
    ChatCell 标注工作流管理器。
    
    参数:
        marker_file: 包含标记基因的 CSV/TSV 文件路径。
        step: 当前步骤，'generate' (生成提示词) 或 'parse' (解析回复)。
        response_text: AI 返回的文本内容 (仅用于 parse 步骤)。
        species: 物种名称 (例如 Human, Mouse)，用于构建提示词。
        tissue: 组织名称 (例如 PBMC, Lung)，用于构建提示词。
        top_n: 每个聚类提取的标记基因数量。
        mode: 提示词模式，'concise' (简洁) 或 'detailed' (详细)。
        source: 数据来源，'scanpy' 或 'seurat'。
    """
    
    # 检查输入文件是否存在
    if not os.path.exists(marker_file):
         raise FileNotFoundError(f"Marker file not found: {marker_file}")

    # 第一步：提取 Marker 并生成 Prompt (Extract Markers)
    if step == "generate":
        # 从文件中提取 Marker 数据，返回字典格式 {cluster_key: "gene1, gene2..."}
        markers = extract_markers_from_file(marker_file, top_n=top_n, source=source)
        # 根据提取的数据生成给 LLM 的提示词
        return generate_annotation_prompt(markers, species, tissue, mode=mode)

    # 第二步：解析 LLM 的回复 (Parse LLM Response)
    elif step == "parse":
        if not response_text:
            raise ValueError("For step='parse', you must provide 'response_text'.")
        
        # 我们需要聚类名称来验证解析结果的一致性
        # 重新提取 Marker 仅仅是为了获取 keys (确保顺序与生成时一致)
        # 注意：这里假设文件在两步之间没有发生变化
        markers = extract_markers_from_file(marker_file, top_n=top_n, source=source)
        cluster_names = list(markers.keys())

        # 调用解析器处理 LLM 返回的文本
        annotations, extra_info = parse_llm_response(response_text, cluster_names)
        
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


