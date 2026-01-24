import pyperclip
import platform

def copy_to_clipboard(text: str):
    """
    Attempts to copy text to the system clipboard.
    """
    try:
        pyperclip.copy(text)
        return True
    except Exception as e:
        print(f"Warning: Could not copy to clipboard automatically: {e}")
        return False

def generate_annotation_prompt(
    markers: dict,
    species: str = "Human",
    tissue: str = "PBMC",
    mode: str = "concise",
    auto_copy: bool = True
) -> str:
    """
    Generates the prompt for the LLM (Copilot, ChatGPT, DeepSeek, etc.).
    
    Args:
        markers: Dictionary of cluster -> markers
        species: Species name
        tissue: Tissue name
        mode: Output mode. Options: 'concise' (default), 'detailed'
        auto_copy: Whether to copy to clipboard automatically
    """
    
    marker_lines = []
    for cluster, gene_str in markers.items():
        marker_lines.append(f"{cluster}: {gene_str}")
    
    marker_block = "\n".join(marker_lines)
    num_clusters = len(markers)
    
    base_instruction = f"Identify cell types of {species} {tissue} cells using the following markers separately for each row.\n" \
                       f"You MUST use standardized cell type names from the Cell Ontology (CL)."

    mode_instruction = ""
    example_block = ""

    if mode == "concise":
        mode_instruction = "Provide the cell type name for each cluster. Format: ClusterX: Cell Type."
        example_block = """Example Output:
Cluster0: CD4+ T cell
Cluster1: B cell
Cluster2: CD14+ Monocyte"""
    elif mode == "detailed":
        mode_instruction = "For each cluster, provide the Cell Type followed by a detailed explanation. Include recommended markers, their functions, and their ranks in the original list."
        example_block = """Example Output:
Cluster0: CD4+ T Cell | Recommended Markers: CD3D, CD4 | Marker Functions: CD3D (T-cell receptor complex), CD4 (Helper T-cell marker) | Ranks: CD3D (Rank 1), CD4 (Rank 2)
Cluster1: B Cell | Recommended Markers: CD19, MS4A1 | Marker Functions: CD19 (B-cell activation), MS4A1 (B-cell receptor signaling) | Ranks: CD19 (Rank 3), MS4A1 (Rank 5)"""
    else:
        raise ValueError("Invalid mode. Choose from 'concise', 'detailed'.")

    prompt = f"""{base_instruction}

{mode_instruction}

IMPORTANT: Return exactly {num_clusters} lines, one for each row. 
Do not use Markdown header or code blocks. Just plain text lines.
Do not add "Here is the list" or "Sure".

{example_block}

---
Task Data:
{marker_block}
"""

    print("=" * 80)
    print("🤖 ChatCell: AI Prompt Generated")
    print("=" * 80)
    
    if auto_copy:
        success = copy_to_clipboard(prompt)
        if success:
            print("✅ Prompt has been COPIED to your clipboard!")
            print("👉 Go to your AI Chat (Copilot, DeepSeek, ChatGPT) and press Ctrl+V (Paste), then Enter.")
        else:
            print("📋 Please copy the prompt below manually:")
            print(prompt)
    else:
        print(prompt)
        
    print("=" * 80)
    
    return prompt
