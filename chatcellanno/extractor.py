import pandas as pd
import os

def extract_markers_from_df(df: pd.DataFrame, top_n: int = 10):
    """
    Extracts markers from a pandas DataFrame.
    Supports two formats:
    1. Tidy/Long format: Columns like 'gene', 'cluster', 'avg_log2FC' (standard Seurat/Scanpy output table).
    2. Wide/Matrix format: Columns are cluster names, values are genes.
    """
    processed_markers = {}
    
    # Check for Tidy/Long format
    # Common column names for clustering results
    cols = [c.lower() for c in df.columns]
    
    # Identify key columns
    cluster_col = next((c for c in df.columns if c.lower() in ['cluster', 'group', 'leiden', 'louvain']), None)
    gene_col = next((c for c in df.columns if c.lower() in ['gene', 'names', 'symbol', 'feature']), None)
    
    # If we found both 'cluster' and 'gene' columns, treat as Long format
    if cluster_col and gene_col:
        # Standardize sorting if possible (optional: relies on input being pre-sorted or having a stats column)
        # We assume the input file is already sorted by significance/fold-change, respecting input order.
        
        # Group by cluster and extract top N genes
        grouped = df.groupby(cluster_col)
        for name, group in grouped:
            # Get top N genes from this group
            # Assuming the file is already sorted, we just take the first N
            genes = group[gene_col].astype(str).values[:top_n]
            processed_markers[str(name)] = ", ".join(genes)
            
    else:
        # Fallback to Wide/Matrix format (Columns = Clusters)
        for col in df.columns:
            # Take top n values, convert to string
            genes = df[col].dropna().astype(str).tolist()
            group_genes = genes[:top_n]
            processed_markers[str(col)] = ", ".join(group_genes)
            
    return processed_markers

def extract_markers_from_file(file_path: str, top_n: int = 10, sep: str = None):
    """
    Extracts markers from a CSV or TSV file.
    """
    if sep is None:
        if file_path.endswith('.csv'):
            sep = ','
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            sep = '\t'
        else:
            sep = ',' # Default fallback
            
    try:
        # Read the file
        df = pd.read_csv(file_path, sep=sep)
        return extract_markers_from_df(df, top_n)
    except Exception as e:
        raise ValueError(f"Failed to read marker file: {e}")


