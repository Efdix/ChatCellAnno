"""
数据提取模块 (Data Extraction Module)

负责从 CSV/TSV 文件或 Pandas DataFrame 中提取标记基因 (Marker Genes)。
支持两种主要的数据格式：
1. Tidy/Long Format (长宽表): 常见的分析软件输出格式 (如 Scanpy, Seurat)，包含 'cluster' 和 'gene' 列。
2. Wide/Matrix Format (宽表): 每一列代表一个 Cluster，列中的值为 Gene Names。
"""

import pandas as pd
import os

def extract_markers_from_df(df: pd.DataFrame, source: str = 'scanpy', top_n: int = 10):
    """
    从 Pandas DataFrame 中提取标记基因。
    
    参数:
        df: 输入的 DataFrame 数据。
        source: 数据来源 ('scanpy' 或 'seurat')。
        top_n: 每个 Cluster 提取前 N 个基因。
        
    返回:
        dict: {cluster_name: "gene1, gene2, ..."}
    """
    processed_markers = {}
    source = source.lower()
    
    cluster_col = None
    gene_col = None

    # 1. Scanpy Logic
    # 对应 sc.get.rank_genes_groups_df 输出
    # 通常期望列: 'names' (基因), 'group' (聚类)
    if source == 'scanpy':
        possible_gene_cols = ['names', 'gene', 'symbol']
        possible_cluster_cols = ['group', 'cluster', 'leiden', 'louvain']
        
        gene_col = next((c for c in df.columns if c.lower() in possible_gene_cols), None)
        cluster_col = next((c for c in df.columns if c.lower() in possible_cluster_cols), None)
        
        if not gene_col or not cluster_col:
            raise ValueError(f"Scanpy data requires columns for genes ({possible_gene_cols}) and clusters ({possible_cluster_cols}). Found: {list(df.columns)}")

    # 2. Seurat Logic
    # 对应 FindAllMarkers 输出
    # 通常期望列: 'gene', 'cluster'
    elif source == 'seurat':
        possible_gene_cols = ['gene', 'feature']
        possible_cluster_cols = ['cluster']
        
        gene_col = next((c for c in df.columns if c.lower() in possible_gene_cols), None)
        cluster_col = next((c for c in df.columns if c.lower() in possible_cluster_cols), None)
        
        if not gene_col or not cluster_col:
            raise ValueError(f"Seurat data requires columns for genes ({possible_gene_cols}) and clusters ({possible_cluster_cols}). Found: {list(df.columns)}")
            
    else:
        raise ValueError("Source must be 'scanpy' or 'seurat'.")

    # 执行提取
    # 按 Cluster 分组 (Group by cluster)
    grouped = df.groupby(cluster_col)
    for name, group in grouped:
        # 获取该组的前 N 个基因
        genes = group[gene_col].astype(str).values[:top_n]
        # 将基因列表转换为逗号分隔的字符串
        processed_markers[str(name)] = ", ".join(genes)
            
    return processed_markers

def extract_markers_from_file(file_path: str, source: str = 'scanpy', top_n: int = 10, sep: str = None):
    """
    从 CSV 或 TSV 文件中提取标记基因。
    作为 `extract_markers_from_df` 的文件输入包装器。
    
    参数:
        file_path: 文件路径。
        source: 数据来源 ('scanpy' 或 'seurat')。
        top_n: 提取数量。
        sep: 分隔符 (如果为 None 则自动根据文件扩展名推断)。
    """
    if sep is None:
        if file_path.endswith('.csv'):
            sep = ','
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            sep = '\t'
        else:
            sep = ',' # 默认回退值
            
    try:
        # 读取文件到 Pandas DataFrame
        df = pd.read_csv(file_path, sep=sep)
        return extract_markers_from_df(df, source=source, top_n=top_n)
    except Exception as e:
        raise ValueError(f"Failed to read marker file (读取文件失败): {e}")


