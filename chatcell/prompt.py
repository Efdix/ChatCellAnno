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
        mode: Output mode. Options: 'concise' (default), 'evidence', 'recommendation'
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
    if mode == "concise":
        mode_instruction = "Only provide the cell type name. Do not include any numbers or extra annotations before the name. " \
                           "Do not include any explanatory text, introductory phrases, or descriptions."
    elif mode == "evidence":
        mode_instruction = "For each row, provide the Cell Type followed by a pipe character '|' and then a list of the specific genes from my provided list that most strongly support this decision.\n" \
                           "Format: Cell Type Name | Supported by: gene1, gene2"
    elif mode == "recommendation":
        mode_instruction = "For each row, provide the Cell Type followed by a pipe character '|' and then recommend 2-3 canonical marker genes that are MISSING from my list but would confirm this cell type.\n" \
                           "Format: Cell Type Name | Recommended Markers: gene1, gene2"
    else:
        raise ValueError("Invalid mode. Choose from 'concise', 'evidence', 'recommendation'.")

    prompt = f"""{base_instruction}
{mode_instruction}
IMPORTANT: Return exactly {num_clusters} lines, one for each row. Do not use Markdown header or code blocks if not necessary, just plain text lines.

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
