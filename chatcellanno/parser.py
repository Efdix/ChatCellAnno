"""
解析模块 (Response Parser Module)

负责处理从 LLM (Copilot, ChatGPT 等) 复制回来的文本。
将非结构化的 AI 回复文本转换为结构化的数据 (字典)，以便后续映射回原始 Cluster。
"""

import warnings

def parse_llm_response(response_text: str, cluster_names: list) -> tuple:
    """
    解析 LLM 的原始文本回复。
    
    逻辑说明：
    1. 清理：去除 Markdown 代码块标记 (```) 和空行。
    2. 验证：检查清洗后的行数是否与 Cluster 数量一致。
    3. 映射：将每一行文本与 Cluster ID 对应起来。
    
    参数:
        response_text: 从 AI 聊天窗口复制的文本。
        cluster_names: Cluster ID 列表 (用于验证和映射)。
        
    返回:
        tuple: (annotations_dict, extra_info_dict)
               annotations_dict: {cluster_name: annotation}
               extra_info_dict:  {cluster_name: extra_info} (如果没有额外信息则为空)
    """
    # 1. 清理 Markdown 代码块 (Cleanup Markdown Code blocks)
    lines = response_text.strip().split('\n')
    clean_lines = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # 跳过代码块标记
        if line.startswith('```'):
            continue
        clean_lines.append(line)
        
    # 2. 长度验证 (Validation)
    # 我们假设 AI 能够严格遵守 "一行对应一个 Cluster" 的指令
    if len(clean_lines) != len(cluster_names):
        warnings.warn(
            f"⚠️ Length Mismatch! (长度不匹配) Expected {len(cluster_names)} annotations, got {len(clean_lines)}.\n"
            f"Mapping might be incorrect. Please check the AI output."
        )
        # 优雅处理不匹配情况 (截断或填充)
        if len(clean_lines) > len(cluster_names):
            # 如果 AI 返回多了，取前 N 行
            clean_lines = clean_lines[:len(cluster_names)]
        else:
            # 如果 AI 返回少了，用 "Unknown" 填充
            clean_lines += ["Unknown"] * (len(cluster_names) - len(clean_lines))
            
    # 3. 映射 (Map)
    annotations = {}
    extra_info_map = {}
    
    for i, cluster in enumerate(cluster_names):
        text = clean_lines[i]
        
        # 检查是否存在分隔符 '|' (Check for pipe separator)
        # 可以在 prompt.py 中看到，detailed 模式会使用 '|' 分隔
        if '|' in text:
            parts = text.split('|', 1)
            main_anno = parts[0].strip() # 细胞类型名称
            extra = parts[1].strip()     # 额外详细信息
            
            annotations[cluster] = main_anno
            extra_info_map[cluster] = extra
        else:
            annotations[cluster] = text
            # 没有额外信息
            
    return annotations, extra_info_map

def generate_annotation_code(annotations: dict, source: str) -> str:
    """
    根据解析的注释结果，生成对应的 Python (Scanpy) 或 R (Seurat) 代码。
    
    参数:
        annotations: 注释字典 {cluster_id: cell_type}
        source: 'scanpy' or 'seurat'
        
    返回:
        str: 对应的代码片段
    """
    source = source.lower()
    code_snippet = ""
    
    if source == 'scanpy':
        code_snippet += "# Python / Scanpy Annotation Code\n"
        code_snippet += "new_cluster_names = {\n"
        for cluster, cell_type in annotations.items():
            # 转义单引号
            safe_type = cell_type.replace("'", "\\'")
            # 尝试判断 cluster 是否为数字，如果是数字字符串，保留引号较为安全
            code_snippet += f"    '{cluster}': '{safe_type}',\n"
        code_snippet += "}\n\n"
        code_snippet += "# Assuming 'adata' is your AnnData object and 'leiden' is your cluster key\n"
        code_snippet += "# 假设 'adata' 是您的对象，聚类列为 'leiden' (请根据实际情况修改: louvain, cluster等)\n"
        code_snippet += "adata.obs['cell_type'] = adata.obs['leiden'].map(new_cluster_names)\n"
        
    elif source == 'seurat':
        code_snippet += "# R / Seurat Annotation Code\n"
        code_snippet += "new.cluster.ids <- c(\n"
        
        r_assignments = []
        for cluster, cell_type in annotations.items():
            # 转义双引号
            safe_type = cell_type.replace('"', '\\"')
            # R 语言命名向量: c("0" = "T Cell", "1" = "B Cell")
            r_assignments.append(f'  "{cluster}" = "{safe_type}"')
            
        code_snippet += ",\n".join(r_assignments)
        code_snippet += "\n)\n\n"
        code_snippet += "# Rename Idents\n"
        code_snippet += "seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)\n"
        code_snippet += "seurat_obj$cell_type <- Idents(seurat_obj)\n"
        
    else:
        code_snippet = "# Unknown source format. Please choose 'scanpy' or 'seurat'."
        
    return code_snippet
