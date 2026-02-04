"""
解析模块 (Response Parser Module)

负责处理从 LLM (Copilot, ChatGPT 等) 复制回来的文本。
将非结构化的 AI 回复文本转换为结构化的数据 (字典)，以便后续映射回原始 Cluster。
"""

import warnings

def parse_llm_response(response_text: str, cluster_names: list = None) -> tuple:
    """
    解析 LLM 的原始文本回复。
    
    逻辑说明：
    1. 清理：去除 Markdown 代码块标记 (```) 和空行。
    2. 验证：如果提供了 cluster_names，检查行数是否一致。
    3. 映射：将每一行文本与 Cluster ID 对应起来。
    
    参数:
        response_text: 从 AI 聊天窗口复制的文本。
        cluster_names: Cluster ID 列表 (可选)。如果不提供，尝试从表格第一列解析 ID。
        
    返回:
        tuple: (annotations_dict, extra_info_dict)
               annotations_dict: {cluster_name: annotation}
               extra_info_dict:  {cluster_name: extra_info}
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
        # 跳过 Markdown 表格的标题行和分割线
        if "|" in line:
            # 简单判断是否是表头或分割线: 包含多个横杠或包含 "Cluster" 关键字且在第一行
            if "---" in line:
                continue
            if "Cluster" in line and len(clean_lines) == 0:
                # 提取数据的逻辑需要进一步细化，这里先保留该行进行后续过滤
                # 如果第一行包含 "Cluster" 且包含更多列名，通常是表头
                parts = [p.strip() for p in line.split("|") if p.strip()]
                if "Cell Type" in parts or "Annotation" in parts:
                    continue
        
        clean_lines.append(line)
        
    # 2. 映射 (Map)
    annotations = {}
    extra_info_map = {}

    # 模式 A: 已知 Cluster List (Strict Mode)
    if cluster_names:
        # 长度验证
        if len(clean_lines) != len(cluster_names):
            warnings.warn(
                f"⚠️ Length Mismatch! (长度不匹配) Expected {len(cluster_names)} annotations, got {len(clean_lines)}.\n"
                f"Mapping might be incorrect. Please check the AI output."
            )
            if len(clean_lines) > len(cluster_names):
                clean_lines = clean_lines[:len(cluster_names)]
            else:
                clean_lines += ["Unknown"] * (len(cluster_names) - len(clean_lines))
        
        for i, cluster in enumerate(cluster_names):
            text = clean_lines[i]
            
            # 尝试解析表格行
            if '|' in text:
                parts = [p.strip() for p in text.split('|') if p.strip()]
                # 假设格式: | ClusterID | CellType | Reason... |
                # 如果有 cluster_names，我们主要用顺序映射，但也检查一下内容
                if len(parts) >= 2:
                    main_anno = parts[1]
                    extra = " | ".join(parts[2:]) if len(parts) > 2 else ""
                    annotations[cluster] = main_anno
                    extra_info_map[cluster] = extra
                else:
                    annotations[cluster] = text.strip('| ')
            else:
                annotations[cluster] = text

    # 模式 B: 未知 Cluster List (Infer Mode)
    else:
        for line in clean_lines:
            if '|' in line:
                parts = [p.strip() for p in line.split('|') if p.strip()]
                if len(parts) >= 2:
                    # 假设第一列是 ID
                    cluster_id = parts[0]
                    main_anno = parts[1]
                    extra = " | ".join(parts[2:]) if len(parts) > 2 else ""
                    annotations[cluster_id] = main_anno
                    extra_info_map[cluster_id] = extra
            else:
                # 无法解析 ID，只能忽略
                pass

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
