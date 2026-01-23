import anndata
import pandas as pd
import numpy as np
import os

def generate_example_data():
    print("Generating realistic PBMC example data...")
    
    # Define real markers for 5 clusters (approx 20-30 genes each)
    # These are common PBMC markers
    marker_dict = {
        '0': ["IL7R", "CD3D", "CD3E", "LTB", "CD4", "LDHB", "TPT1", "TMSB10", "TRAC", "GEM", "CD27", "CD2", "ITGB1", "PTPRC", "CD28", "CD69", "TRBC1", "TRBC2", "LCK", "FYN", "JUNB", "S100A4", "NOSIP"],
        '1': ["CD14", "LYZ", "S100A8", "S100A9", "CST3", "FCN1", "LST1", "AIF1", "S100A12", "TYROBP", "FCER1G", "CD68", "CD163", "CSF1R", "CCR2", "IL1B", "CX3CR1", "SERPINA1", "CSTA", "SPI1", "MNDA", "CTSS", "FCN1"],
        '2': ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "CD19", "BANK1", "PAX5", "BLK", "FCRL1", "FCRL2", "IGHM", "IGHD", "IGKC", "IGLC1", "TNFRSF13C", "CD22", "CD37", "CD40", "RALGPS2", "MEF2C", "HVCN1"],
        '3': ["GNLY", "NKG7", "CST7", "GZMB", "KLRB1", "CD247", "FCGR3A", "GZMA", "PRF1", "FGFBP2", "KLRD1", "SPON2", "GZMH", "KLRF1", "HOPX", "CLIC3", "CTSW", "IL2RB", "NCR1", "XCL1", "XCL2"],
        '4': ["FCER1A", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CLEC10A", "CD1C", "CD1E", "PLD4", "ITGAX", "CCR6", "IRF4", "TMEM176B", "CD1B", "PKIB", "NDRG2", "CPVL", "LTA", "VIM"]
    }
    
    clusters = list(marker_dict.keys())
    # Find max length to pad
    max_len = max(len(genes) for genes in marker_dict.values())
    
    # Create the 'names' structured array (or DF)
    # Pad with 'NaN' or just repeat last gene (Scanpy usually has fixed shape)
    # For simplicity, we just take the first 20 genes for the strict structure if exists, 
    # but let's try to be flexible and use exact lists via a DataFrame that might have NaNs (ChatCell handles string conversion)
    
    # Let's just crop to min_len to keep it rectangular and simple for this mock
    min_len = min(len(genes) for genes in marker_dict.values())
    print(f"Using top {min_len} genes per cluster for the mock structure.")
    
    data_dict = {k: v[:min_len] for k, v in marker_dict.items()}
    names_df = pd.DataFrame(data_dict)
    
    # 1. Create a mock AnnData object
    # We don't need real counts for ChatCell's extraction part if we mock .uns directly
    # But let's make it consistent.
    n_obs = 500
    n_vars = 100 # minimal vars
    obs = pd.DataFrame({
        'leiden': np.random.choice(clusters, n_obs)
    }, index=[f'cell_{i}' for i in range(n_obs)])
    obs['leiden'] = obs['leiden'].astype('category')
    
    # Dummy X
    adata = anndata.AnnData(X=np.zeros((n_obs, n_vars)), obs=obs)
    
    # Mock rank_genes_groups
    adata.uns['rank_genes_groups'] = {}
    adata.uns['rank_genes_groups']['params'] = {'groupby': 'leiden', 'method': 'wilcoxon', 'reference': 'rest', 'use_raw': False}
    
    adata.uns['rank_genes_groups']['names'] = names_df.to_records(index=False)
    
    # Mock other keys required by scanpy structure (optional but good for stability)
    # Just zeros or randoms
    dummy_matrix = np.zeros((min_len, len(clusters)))
    adata.uns['rank_genes_groups']['scores'] = pd.DataFrame(dummy_matrix, columns=clusters).to_records(index=False)
    adata.uns['rank_genes_groups']['pvals'] = pd.DataFrame(dummy_matrix, columns=clusters).to_records(index=False)
    adata.uns['rank_genes_groups']['pvals_adj'] = pd.DataFrame(dummy_matrix, columns=clusters).to_records(index=False)
    adata.uns['rank_genes_groups']['logfoldchanges'] = pd.DataFrame(dummy_matrix, columns=clusters).to_records(index=False)
    
    # Save AnnData
    adata_path = "example_realistic_adata.h5ad"
    adata.write(adata_path)
    print(f"✅ Saved {adata_path}")
    
    # 2. Save TSV
    tsv_path = "example_realistic_markers.tsv"
    names_df.to_csv(tsv_path, sep='\t', index=False)
    print(f"✅ Saved {tsv_path}")
    
    print("\nSimulated Biology:")
    print("Cluster 0: T cells markers")
    print("Cluster 1: Monocytes markers")
    print("Cluster 2: B cells markers")
    print("Cluster 3: NK cells markers")
    print("Cluster 4: Dendritic cells markers")

if __name__ == "__main__":
    generate_example_data()
