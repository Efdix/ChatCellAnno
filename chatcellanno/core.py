from .extractor import extract_markers_from_file
from .prompt import generate_annotation_prompt
from .parser import parse_llm_response
import pandas as pd
import os

def annotate_cell_types(
    marker_file: str,
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
        marker_file: Path to a CSV/TSV file containing markers.
        step: 'generate' or 'parse'
        response_text: AI's response (for parse step)
        species: Species validation
        tissue: Tissue context
        top_n: Number of markers per cluster
        mode: 'concise', 'evidence', or 'recommendation'
    """
    
    if not os.path.exists(marker_file):
         raise FileNotFoundError(f"Marker file not found: {marker_file}")

    # Extract Markers
    if step == "generate":
        markers = extract_markers_from_file(marker_file, top_n=top_n)
        return generate_annotation_prompt(markers, species, tissue, mode=mode)

        
    elif step == "parse":
        if not response_text:
            raise ValueError("For step='parse', you must provide 'response_text'.")
        
        # We need cluster names to validate the parse
        # Re-extract to ensure alignment
        markers = extract_markers_from_file(marker_file, top_n=top_n)
        cluster_names = list(markers.keys())

        annotations, extra_info = parse_llm_response(response_text, cluster_names)
        
        # Create a result DataFrame
        results = []
        for cluster in cluster_names:
            row = {
                'Cluster': cluster,
                'Annotation': annotations.get(cluster, 'Unknown'),
            }
            if extra_info:
                row['Extra_Info'] = extra_info.get(cluster, '')
            results.append(row)
            
        result_df = pd.DataFrame(results)
        print("✅ Annotations Parsed Successfully.")
        return result_df
    
    else:
        raise ValueError("Invalid step. Use 'generate' or 'parse'.")


