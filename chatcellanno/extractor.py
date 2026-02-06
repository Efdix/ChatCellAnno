"""
数据提取模块 (Data Extraction Module)

负责从 CSV/TSV 文件或 Pandas DataFrame 中提取标记基因 (Marker Genes)。
支持两种主要的数据格式：
1. Tidy/Long Format (长宽表): 常见的分析软件输出格式 (如 Scanpy, Seurat)，包含 'cluster' 和 'gene' 列。
2. Wide/Matrix Format (宽表): 每一列代表一个 Cluster，列中的值为 Gene Names。
"""

import pandas as pd
import os

def extract_markers_from_df(df: pd.DataFrame, top_n: int = 10):
    """
    Auto-detect and extract marker genes from DataFrame.
    Supports both Scanpy and Seurat standard outputs.
    """
    processed_markers = {}
    
    # Combined search for common columns in Scanpy/Seurat
    possible_gene_cols = ['names', 'gene', 'symbol', 'feature']
    possible_cluster_cols = ['group', 'cluster', 'leiden', 'louvain']
    
    gene_col = next((c for c in df.columns if c.lower() in possible_gene_cols), None)
    cluster_col = next((c for c in df.columns if c.lower() in possible_cluster_cols), None)
    
    if not gene_col or not cluster_col:
        raise ValueError(f"Could not identify gene/cluster columns. \nFound: {list(df.columns)}\nNeed one of {possible_gene_cols} and {possible_cluster_cols}")

    # 执行提取
    # 按 Cluster 分组 (Group by cluster)
    grouped = df.groupby(cluster_col)
    for name, group in grouped:
        # Get specified top_n genes
        genes = group[gene_col].astype(str).values[:top_n]
        processed_markers[str(name)] = ", ".join(genes)
            
    return processed_markers

def extract_markers_from_file(file_path: str, top_n: int = 10, sep: str = None):
    """
    Auto-detect and extract from file.
    """
    if sep is None:
        if file_path.endswith('.csv'):
            sep = ','
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            sep = '\t'
        else:
            sep = ',' 
            
    try:
        df = pd.read_csv(file_path, sep=sep)
        return extract_markers_from_df(df, top_n=top_n)
    except Exception as e:
        raise e
    except Exception as e:
        raise ValueError(f"Failed to read marker file (读取文件失败): {e}")


