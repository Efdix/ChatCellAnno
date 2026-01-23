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
        
        scores = pd.DataFrame(results['scores']) # or logfoldchanges
        pvals = pd.DataFrame(results['pvals_adj'])
        
        # Check if 'logfoldchanges' exists, prefer it over scores for formatting, though extraction logic is similar
        if 'logfoldchanges' in results:
            lfc = pd.DataFrame(results['logfoldchanges'])
        else:
            lfc = scores
            
        for group in groups:
            # Create a temporary DF for this group to sort/filter if needed
            # Scanpy usually already sorts them, but we take top N just in case
            group_genes = names[group].values[:top_n]
            
            # Filter? Usually scanpy results are already ranked. 
            # We assume the user has run rank_genes_groups with appropriate method (e.g. wilcoxon)
            
            # Convert to comma separated string
            marker_str = ", ".join([str(g) for g in group_genes])
            processed_markers[group] = marker_str
            
    except Exception as e:
        raise ValueError(f"Failed to parse rank_genes_groups structure: {str(e)}")
        
    return processed_markers
