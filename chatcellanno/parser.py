"""
响应解析模块 (Response Parser Module)

负责处理从 LLM (Copilot, ChatGPT, DeepSeek 等) 对话框中复制回来的预测结果。
本模块的核心职能是将非结构化的、可能带有 Markdown 格式的 AI 回复文本，
精准地转换为 Python 字典结构，从而实现聚类 ID (Cluster ID) 与细胞类型 (Cell Type) 的自动映射。
"""

import warnings

def parse_llm_response(response_text: str, cluster_names: list = None) -> tuple:
    """
    解析 LLM 的原始文本回复并提取关键注释信息 (Parse LLM response and extract annotations).
    
    算法逻辑 (Algorithm Logic):
    1. 清理 (Cleanup): 自动识别并去除 Markdown 表格边框、代码块标记 (```) 以及空行。
    2. 校验 (Validation): 针对已知的聚类列表进行长度校验，确保每一行都能对应到一个物理聚类。
    3. 映射 (Mapping): 
       - 优先从表格列中提取 Cluster ID 和 Cell Type。
       - 如果表格格式不规范，则采用“行对齐”模式进行强行映射。
    
    参数 (Parameters):
    - response_text: AI 生成的回复全文。
    - cluster_names: 预期的聚类 ID 列表 (如 ['0', '1', '2'])。
    
    返回 (Returns):
    - tuple: (annotations_dict, extra_info_dict)
             annotations_dict: {聚类ID: 细胞类型名称}
             extra_info_dict:  {聚类ID: AI 提供的生物学理由/备注}
    """
    # 1. 预处理：清洗 Markdown 代码块与表格噪声 (Pre-processing: Clean Markdown & Table noise)
    lines = response_text.strip().split('\n')
    clean_lines = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # 跳过 Markdown 围栏代码块标记 (Skip code block markers)
        if line.startswith('```'):
            continue
        # 智能过滤 Markdown 表格的表头行和装饰线 (Filter table headers and separators)
        if "|" in line:
            # 如果包含连续横线，判定为表格分割线 (e.g., |---|---|)
            if "---" in line:
                continue
            # 如果是首行且包含关键字，判定为表头 (e.g., | Cluster | Cell Type |)
            if "Cluster" in line and len(clean_lines) == 0:
                parts = [p.strip() for p in line.split("|") if p.strip()]
                if "Cell Type" in parts or "Annotation" in parts:
                    continue
        
        clean_lines.append(line)
        
    # 2. 核心映射逻辑 (Core Mapping Logic)
    annotations = {}
    extra_info_map = {}

    # 场景 A: 存在预期的聚类列表 (Strict/Ordered Mode)
    # 此模式下，即使 AI 没有在表格中写明 ID，我们也会按顺序将其填入对应的聚类。
    if cluster_names:
        # 长度容错处理 (Length Mismatch Handling)
        if len(clean_lines) != len(cluster_names):
            warnings.warn(
                f"⚠️ Length Mismatch! (长度不匹配) Expected {len(cluster_names)} annotations, got {len(clean_lines)}.\n"
                f"Mapping might be incorrect. Please check the AI output."
            )
            # 补齐或截断，防止程序崩溃 (Pad or truncate to prevent crash)
            if len(clean_lines) > len(cluster_names):
                clean_lines = clean_lines[:len(cluster_names)]
            else:
                clean_lines += ["Unknown"] * (len(cluster_names) - len(clean_lines))
        
        for i, cluster in enumerate(cluster_names):
            text = clean_lines[i]
            
            # 尝试从 Markdown 表格行中提取字段内容 (Extract fields from Markdown row)
            if '|' in text:
                parts = [p.strip() for p in text.split('|') if p.strip()]
                # 标准格式假设: | ClusterID | CellType | Reasoning |
                if len(parts) >= 2:
                    # 优先取第二列作为细胞类型 (Prioritize 2nd col as cell type)
                    main_anno = parts[1]
                    # 将第三列及以后的内容合并作为“详细理由” (Merge remaining cols as extra info)
                    extra = " | ".join(parts[2:]) if len(parts) > 2 else ""
                    annotations[cluster] = main_anno
                    extra_info_map[cluster] = extra
                else:
                    # 如果只有一列有效内容，直接作为类型 (Single col fallback)
                    annotations[cluster] = text.strip('| ')
            else:
                # 非表格行，直接整行映射 (Plain text fallback)
                annotations[cluster] = text

    # 场景 B: 未知聚类列表 (Infer/Heuristic Mode)
    # 当用户跳过 Marker 文件直接粘贴 AI 结果时，尝试从文本内容中反演 Cluster ID。
    else:
        for line in clean_lines:
            if '|' in line:
                parts = [p.strip() for p in line.split('|') if p.strip()]
                if len(parts) >= 2:
                    # 假设首列即为 Cluster ID (Assume 1st col is the key)
                    cluster_id = parts[0]
                    main_anno = parts[1]
                    extra = " | ".join(parts[2:]) if len(parts) > 2 else ""
                    annotations[cluster_id] = main_anno
                    extra_info_map[cluster_id] = extra
            else:
                # 无法识别的行暂不处理
                pass
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
