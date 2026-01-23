from .extractor import extract_markers_from_scanpy
from .prompt import generate_annotation_prompt
from .parser import parse_llm_response
import anndata

def annotate_cell_types(
    adata: anndata.AnnData,
    step: str = "generate",
    response_text: str = None,
    species: str = "Human",
    tissue: str = "PBMC",
    top_n: int = 10,
    mode: str = "concise"
):
    """
    Workflow manager for ChatCell annotation.
    
    Args:
        adata: AnnData object
        step: 'generate' or 'parse'
        response_text: AI's response (for parse step)
        species: Species validation
        tissue: Tissue context
        top_n: Number of markers
        mode: 'concise', 'evidence', or 'recommendation'
    
    Usage:
    
    # Step 1: Generate Prompt
    chatcell.annotate_cell_types(adata, step="generate", mode="evidence")
    
    # ... User pastes to AI Chat, gets response ...
    
    # Step 2: Parse and Apply
    ai_output = "..."
    chatcell.annotate_cell_types(adata, step="parse", response_text=ai_output)
    """
    
    # Extract Markers
    try:
        markers = extract_markers_from_scanpy(adata, top_n=top_n)
    except Exception as e:
        print("❌ Error extracting markers from adata. Did you run `sc.tl.rank_genes_groups`?")
        raise e
        
    cluster_names = list(markers.keys())
    
    if step == "generate":
        return generate_annotation_prompt(markers, species, tissue, mode=mode)
        
    elif step == "parse":
        if not response_text:
            raise ValueError("For step='parse', you must provide 'response_text'.")
            
        annotations, extra_info = parse_llm_response(response_text, cluster_names)
        
        # Apply to AnnData
        # We assume the user wants to add this to .obs?
        # Typically annotations are mapped by cluster ID.
        
        # 1. Create a mapping dictionary
        # annotations is {cluster_id: cell_type}
        
        # 2. Map to a new column
        # Assuming the cluster column used for rank_genes_groups is available. 
        # Scanpy stores the grouping key in adata.uns['rank_genes_groups']['params']['groupby']
        
        try:
            groupby_key = adata.uns['rank_genes_groups']['params']['groupby']
            new_col_name = f"chatcell_annotation"
            
            # Map the annotations
            # Convert annotations dict to a mapper that handles the types correctly
            adata.obs[new_col_name] = adata.obs[groupby_key].map(annotations).astype('category')
            print(f"✅ Annotations added to adata.obs['{new_col_name}']")
            
            # 3. Handle extra info (evidence/recommendation) if present
            if extra_info:
                extra_col_name = f"chatcell_extra_info"
                adata.obs[extra_col_name] = adata.obs[groupby_key].map(extra_info).astype('object')
                print(f"✅ Extra info added to adata.obs['{extra_col_name}']")
            
            return adata.obs[new_col_name]
            
        except KeyError:
            print("⚠️ Could not automatically detect 'groupby' column to apply annotations directly to .obs.")
            print("Returning the annotation dictionary instead.")
            if extra_info:
                return annotations, extra_info
            return annotations
            
    else:
        raise ValueError("Invalid step. Use 'generate' or 'parse'.")
