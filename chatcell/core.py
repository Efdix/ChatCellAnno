from .extractor import extract_markers_from_scanpy, extract_markers_from_file
from .prompt import generate_annotation_prompt
from .parser import parse_llm_response
import anndata
import pandas as pd
import os

def annotate_cell_types(
    adata: anndata.AnnData = None,
    marker_file: str = None,
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
        adata: AnnData object (optional if marker_file is provided)
        marker_file: Path to a CSV/TSV file containing markers (optional if adata is provided)
        step: 'generate' or 'parse'
        response_text: AI's response (for parse step)
        species: Species validation
        tissue: Tissue context
        top_n: Number of markers
        mode: 'concise', 'evidence', or 'recommendation'
    
    Usage:
    
    # Method 1: Scanpy Object
    chatcell.annotate_cell_types(adata=adata, step="generate")
    
    # Method 2: Marker File
    chatcell.annotate_cell_types(marker_file="markers.tsv", step="generate")
    """
    
    # Extract Markers
    if step == "generate":
        markers = {}
        if adata is not None:
            try:
                markers = extract_markers_from_scanpy(adata, top_n=top_n)
            except Exception as e:
                print("❌ Error extracting markers from adata. Did you run `sc.tl.rank_genes_groups`?")
                raise e
        elif marker_file is not None:
             if not os.path.exists(marker_file):
                 raise FileNotFoundError(f"Marker file not found: {marker_file}")
             markers = extract_markers_from_file(marker_file, top_n=top_n)
        else:
            raise ValueError("You must provide either 'adata' or 'marker_file' for step='generate'.")
            
        return generate_annotation_prompt(markers, species, tissue, mode=mode)
        
    elif step == "parse":
        if not response_text:
            raise ValueError("For step='parse', you must provide 'response_text'.")
        
        # We need cluster names to validate the parse
        cluster_names = []
        
        # If adata is present, we try to get names from it to be safe, but strictly we need the keys from the prompt generation.
        # But here we might not know them if we don't re-extract. 
        # Ideally, the user provides the same context. 
        # Let's try to infer cluster names from the source again if available.
        
        if adata is not None:
            # We assume the user has the same adata state
             markers = extract_markers_from_scanpy(adata, top_n=top_n)
             cluster_names = list(markers.keys())
        elif marker_file is not None:
             markers = extract_markers_from_file(marker_file, top_n=top_n)
             cluster_names = list(markers.keys())
        else:
             # If neither is provided, we can't fully validate length match against original clusters easily
             # unless we change the API to require cluster_names.
             # For now, let's warn and try to rely on the parser to just return what it finds if we can't validate.
             # Actually, parse_llm_response REQUIRES cluster_names to map row 1 to cluster A.
             raise ValueError("For step='parse', you must provide 'adata' or 'marker_file' again to align the AI response to the correct clusters.")

        annotations, extra_info = parse_llm_response(response_text, cluster_names)
        
        # Apply to AnnData if available
        if adata is not None:
            try:
                # We assume the user wants to add this to .obs?
                # Typically annotations are mapped by cluster ID.
                # Scanpy stores the grouping key in adata.uns['rank_genes_groups']['params']['groupby']
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
            # If no adata (e.g. file mode), just return the dicts
            print("✅ Parsed annotations returned as dictionary.")
            if extra_info:
                return annotations, extra_info
            return annotations
            
    else:
        raise ValueError("Invalid step. Use 'generate' or 'parse'.")

