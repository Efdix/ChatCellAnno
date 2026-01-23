import scanpy as sc
import chatcell
import os

def main():
    adata_path = "example_realistic_adata.h5ad"
    if not os.path.exists(adata_path):
        print(f"Please run gen_data.py first to generate {adata_path}!")
        return

    print("--- Loading Realistic PBMC AnnData ---")
    adata = sc.read_h5ad(adata_path)
    print(f"Loaded data with {adata.n_obs} cells and {len(adata.obs['leiden'].unique())} clusters.")
    
    print("\n--- Step 1: Generate Prompt ---")
    # This will copy prompt to clipboard
    prompt = chatcell.annotate_cell_types(
        adata=adata,
        step="generate",
        species="Human",
        tissue="PBMC",
        mode="concise"
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
    
    # Simple multi-line input simulation
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
        # Realistic simulated response based on the markers in gen_data.py
        # 0: T cells, 1: Mono, 2: B cells, 3: NK, 4: DC
        ai_response = """CD4+ T Cell
CD14+ Monocyte
B Cell
NK Cell
Dendritic Cell"""
        print("\n[Using Simulated Response based on Ground Truth]")
        print(ai_response)
    
    print("\n--- Step 2: Parse and Apply ---")
    chatcell.annotate_cell_types(
        adata=adata,
        step="parse",
        response_text=ai_response
    )
    
    print("\n--- Result in adata.obs (First 10 cells) ---")
    print(adata.obs[['leiden', 'chatcell_annotation']].head(10))
    
    print("\n✅ Demo Complete. In real analysis, you would now run:")
    print("   sc.pl.umap(adata, color='chatcell_annotation')")

if __name__ == "__main__":
    main()
