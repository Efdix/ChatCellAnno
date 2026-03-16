"""
数据提取模块 (Data Extraction Module)

负责从 CSV/TSV 文件或 Pandas DataFrame 中提取标记基因 (Marker Genes)。
本模块具备自动识别主流单细胞分析软件（如 Scanpy 和 Seurat）输出格式的能力。
"""

import pandas as pd
import os

def extract_markers_from_df(df: pd.DataFrame, top_n: int = 10):
    """
    自动检测并从 DataFrame 中提取 Marker 基因 (Auto-detect markers from DF).
    兼容 Scanpy 和 Seurat 的标准输出格式。
    
    参数 (Parameters):
    - df: 包含 Marker 信息的 Pandas DataFrame
    - top_n: 每个聚类保留的前 N 个基因 (用于提示词构造)
    """
    processed_markers = {}
    
    # 启发式搜索 Scanpy/Seurat 常见的列名 (Scan columns based on common biological data conventions)
    possible_gene_cols = ['names', 'gene', 'symbol', 'feature'] # 基因名列
    possible_cluster_cols = ['group', 'cluster', 'leiden', 'louvain'] # 聚类/分组列
    
    # 识别列名 (Identify relevant columns)
    gene_col = next((c for c in df.columns if c.lower() in possible_gene_cols), None)
    cluster_col = next((c for c in df.columns if c.lower() in possible_cluster_cols), None)
    
    # 如果找不到列，抛出详细错误 (Raise error if column structure is unrecognized)
    if not gene_col or not cluster_col:
        raise ValueError(f"Could not identify gene/cluster columns. \nFound: {list(df.columns)}\nNeed one of {possible_gene_cols} and {possible_cluster_cols}")

    # 执行分组提取 (Execute extraction by grouping)
    # 按 Cluster (聚类) 分组并提取对应的基因列表
    grouped = df.groupby(cluster_col)
    for name, group in grouped:
        # 获取指定数量的基因并转换为逗号分隔的字符串 (Extract N genes and join as string)
        genes = group[gene_col].astype(str).values[:top_n]
        processed_markers[str(name)] = ", ".join(genes)
            
    return processed_markers

def extract_markers_from_file(file_path: str, top_n: int = 10, sep: str = None):
    """
    从本地文件中自动识别并提取数据 (Read markers from a local file).
    
    参数 (Parameters):
    - file_path: 文件路径 (.csv, .tsv, .txt)
    - top_n: 提取数量
    - sep: 分隔符 (可选，默认为按扩展名判断)
    """
    # 自动推断分隔符 (Infer separator by file extension)
    if sep is None:
        if file_path.endswith('.csv'):
            sep = ','
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            sep = '\t'
        else:
            sep = ',' # Default fallback
            
    try:
        # 加载数据 (Load data)
        df = pd.read_csv(file_path, sep=sep)
        return extract_markers_from_df(df, top_n=top_n)
    except Exception as e:
        # 捕获并重传错误以供 GUI 显示 (Relay errors for UI handling)
        raise ValueError(f"Failed to read marker file (读取文件失败): {e}")


