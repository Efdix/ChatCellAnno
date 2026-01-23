import pandas as pd
import scanpy as sc
import anndata

def extract_markers_from_scanpy(adata: anndata.AnnData, top_n: int = 10, uns_key: str = 'rank_genes_groups'):
    """
    Extracts top N markers from Scanpy's rank_genes_groups results.
    
    Args:
        adata: Annotated Data Matrix (AnnData)
        top_n: Number of top markers to extract per cluster
        uns_key: Key in adata.uns where rank_genes_groups results are stored
        
    Returns:
        dict: A dictionary mapping cluster names to comma-separated marker strings.
    """
    if uns_key not in adata.uns:
        raise ValueError(f"'{uns_key}' not found in adata.uns. Please run sc.tl.rank_genes_groups first.")
    
    results = adata.uns[uns_key]
    
    processed_markers = {}
    
    try:
        # Handling structured arrays common in scanpy or DataFrames
        names = pd.DataFrame(results['names'])
        
        # Get group names from columns
        groups = names.columns.tolist()
        
        # Note: We currently only use 'names' to get the gene list.
        # Future versions might use 'scores' or 'pvals_adj' to filter, but for now we trust the ranking.
            
        for group in groups:
            # Take top N
            group_genes = names[group].values[:top_n]
            
            # Convert to comma separated string
            marker_str = ", ".join([str(g) for g in group_genes])
            processed_markers[group] = marker_str
            
    except Exception as e:
        raise ValueError(f"Failed to parse rank_genes_groups structure: {str(e)}")
        
    return processed_markers

def extract_markers_from_df(df: pd.DataFrame, top_n: int = 10):
    """
    Extracts markers from a pandas DataFrame where columns are groups and rows are genes.
    """
    processed_markers = {}
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
        else:
            sep = '\t'
            
    try:
        df = pd.read_csv(file_path, sep=sep)
        return extract_markers_from_df(df, top_n)
    except Exception as e:
        raise ValueError(f"Failed to read marker file: {e}")

