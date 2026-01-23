import chatcell
import os
import pandas as pd

def main():
    marker_file = "example_realistic_markers.tsv"
    if not os.path.exists(marker_file):
        print(f"Please run gen_data.py first to generate {marker_file}!")
        return

    print("--- Running from TSV File ---")
    
    print("\n--- Step 1: Generate Prompt ---")
    # This will copy prompt to clipboard
    prompt = chatcell.annotate_cell_types(
        marker_file=marker_file,
        step="generate",
        species="Human",
        tissue="PBMC",
        mode="evidence"
    )
    
    print("\n" + "="*60)
    print("👉 ACTION REQUIRED: ")
    print("1. Go to your AI Chat (ChatGPT/Claude/DeepSeek).")
    print("2. Paste (Ctrl+V) the prompt.")
    print("3. Copy the AI's response.")
    print("="*60 + "\n")

    print("(For this demo, you can either paste the real AI response below,")
    print(" or just press ENTER to use a pre-calculated simulation.)")
    
    print("\nWaiting for your input (Paste AI response, end with empty line or Ctrl+D/Z):")
    
    lines = []
    try:
        while True:
            line = input()
            if not line:
                break
            lines.append(line)
    except EOFError:
        pass
        
    user_input = "\n".join(lines).strip()
    
    if user_input:
        ai_response = user_input
        print("\nUsing YOUR pasted response.")
    else:
        # Realistic simulated response with evidence
        ai_response = """CD4+ T Cell | Supported by: CD3D, CD4, IL7R
CD14+ Monocyte | Supported by: CD14, LYZ
B Cell | Supported by: MS4A1, CD79A
NK Cell | Supported by: GNLY, NKG7
Dendritic Cell | Supported by: FCER1A, HLA-DPB1"""
        print("\n[Using Simulated Response based on Ground Truth]")
        print(ai_response)
    
    print("\n--- Step 2: Parse Results ---")
    # Note: When using file input, 'parse' returns the dictionaries directly
    annotations, evidences = chatcell.annotate_cell_types(
        marker_file=marker_file,
        step="parse",
        response_text=ai_response
    )
    
    print("\n--- Parsed Annotations ---")
    print("Annotations:", annotations)
    print("Evidence:", evidences)
    
    # You can now save this to a CSV if you want
    results_df = pd.DataFrame({
        'Cluster': list(annotations.keys()),
        'CellType': list(annotations.values()),
        'Evidence': [evidences.get(k, "") for k in annotations.keys()]
    })
    print("\n--- Results DataFrame ---")
    print(results_df)

if __name__ == "__main__":
    main()
