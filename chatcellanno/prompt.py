"""
提示词生成模块 (Prompt Engineering Module)

负责构建高质量的 Prompt，指导 LLM 进行细胞注释。
还包括了与系统剪贴板的交互功能，方便用户直接粘贴到 AI 对话框。
"""

import pyperclip
import platform

def copy_to_clipboard(text: str):
    """
    尝试将文本复制到系统剪贴板。
    使用 `pyperclip` 库，兼容 Windows, Mac, Linux。
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
    生成用于 LLM (Copilot, ChatGPT, DeepSeek 等) 的提示词。
    
    参数:
        markers: 字典 {cluster -> markers_string}
        species: 物种名称 (Human, Mouse...)
        tissue: 组织名称 (PBMC, Liver...)
        mode: 输出模式
             - 'concise': 简洁模式，只输出细胞类型名称。
             - 'detailed': 详细模式，包含推荐 Marker 和功能解释。
        auto_copy: 是否自动复制到剪贴板。
    """
    
    # 构建 Marker 数据块
    marker_lines = []
    for cluster, gene_str in markers.items():
        marker_lines.append(f"{cluster}: {gene_str}")
    
    marker_block = "\n".join(marker_lines)
    num_clusters = len(markers)
    
    # 基础指令 (Base Instruction)
    base_instruction = f"Identify cell types of {species} {tissue} cells using the following markers separately for each row.\n" \
                       f"You MUST use standardized cell type names from the Cell Ontology (CL)."

    mode_instruction = ""
    example_block = ""

    # 根据模式设置具体的格式指令
    if mode == "concise":
        mode_instruction = "Provide the cell type name for each cluster. Format: ClusterX: Cell Type."
        example_block = """Example Output:
Cluster0: CD4+ T cell
Cluster1: B cell
Cluster2: CD14+ Monocyte"""
    elif mode == "detailed":
        mode_instruction = "For each cluster, provide the Cell Type followed by a detailed explanation. Include recommended markers, their functions, and their ranks in the original list."
        # 使用对应的 Example，确保 AI 理解 '|' 分隔符
        example_block = """Example Output:
Cluster0: CD4+ T Cell | Recommended Markers: CD3D, CD4 | Marker Functions: CD3D (T-cell receptor complex), CD4 (Helper T-cell marker) | Ranks: CD3D (Rank 1), CD4 (Rank 2)
Cluster1: B Cell | Recommended Markers: CD19, MS4A1 | Marker Functions: CD19 (B-cell activation), MS4A1 (B-cell receptor signaling) | Ranks: CD19 (Rank 3), MS4A1 (Rank 5)"""
    else:
        raise ValueError("Invalid mode. Choose from 'concise', 'detailed'.")

    # 组装最终 Prompt
    # 强调返回行数必须与 Cluster 数量一致，这对于后续解析非常关键
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
            print("✅ Prompt has been COPIED to your clipboard! (已复制到剪贴板)")
            print("👉 Go to your AI Chat (Copilot, DeepSeek, ChatGPT) and press Ctrl+V (Paste), then Enter.")
        else:
            print("📋 Please copy the prompt below manually:")
            print(prompt)
    else:
        print(prompt)
        
    print("=" * 80)
    
    return prompt
