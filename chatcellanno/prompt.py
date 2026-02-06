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
    exclude_types: str = "",
    enrichment_hints: dict = None,
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
        exclude_types: 需要排除的细胞类型（逗号分隔字符串）。
        enrichment_hints: 字典 {cluster -> [hint1, hint2...]}，可选的富集分析结果。
        auto_copy: 是否自动复制到剪贴板。
    """
    
    # 构建 Marker 数据块
    # 如果有 enrichment_hints，我们需要改变格式来包含这些信息
    marker_lines = []
    
    use_hints = enrichment_hints is not None and len(enrichment_hints) > 0

    for cluster, gene_str in markers.items():
        if use_hints and str(cluster) in enrichment_hints:
            hints = enrichment_hints[str(cluster)]
            # 将提示信息拼接到同一行，或者作为子项
            # 为了表格解析器方便，我们尽量保持行结构，建议将 Hints 放在方括号里
            # 例如: Cluster0: Genes... [Hints: T cell activation(P=...), ...]
            hints_str = "; ".join(hints)
            line = f"{cluster}: Markers: {gene_str} | Functional Hints from Database: [{hints_str}]"
            marker_lines.append(line)
        else:
            # 原有的格式
            marker_lines.append(f"{cluster}: {gene_str}")
    
    marker_block = "\n".join(marker_lines)
    num_clusters = len(markers)
    
    # 基础指令 (Base Instruction)
    base_instruction = f"Identify cell types of {species} {tissue} cells using the following markers separately for each row.\n" \
                       f"You MUST use standardized cell type names from the Cell Ontology (CL)."

    if use_hints:
        base_instruction += "\nUse the provided 'Functional Hints' (ORA Enrichment Results) as strong evidence to cross-validate gene markers."

    # 排除项指令 (Exclude Instruction)
    exclude_instruction = ""
    if exclude_types and exclude_types.strip():
        exclude_instruction = f"IMPORTANT: The following cell types are known NOT to be present in this sample. Do NOT identify any cluster as: {exclude_types}."

    mode_instruction = ""
    example_block = ""

    # 根据模式设置具体的格式指令
    if mode == "concise":
        mode_instruction = "For each cluster, provide the result in a Markdown Table format. The table must have columns: 'Cluster' and 'Cell Type'."
        example_block = """Example Output:
| Cluster | Cell Type |
| :--- | :--- |
| Cluster0 | CD4+ T Cell |
| Cluster1 | B Cell |
| Cluster2 | CD14+ Monocyte |"""
    elif mode == "detailed":
        mode_instruction = "For each cluster, provide the result in a Markdown Table format. The table must have columns: 'Cluster', 'Cell Type', 'Recommended Markers', and 'Reasoning/Functions'."
        # 使用对应的 Example，确保 AI 输出表格
        example_block = """Example Output:
| Cluster | Cell Type | Recommended Markers | Reasoning/Functions |
| :--- | :--- | :--- | :--- |
| Cluster0 | CD4+ T Cell | CD3D, CD4 | CD3D is a T-cell receptor complex component... |
| Cluster1 | B Cell | CD19, MS4A1 | CD19 is a classic B-cell marker... |"""
    else:
        raise ValueError("Invalid mode. Choose from 'concise', 'detailed'.")

    # 组装最终 Prompt
    # 强调返回行数必须与 Cluster 数量一致，这对于后续解析非常关键
    prompt = f"""{base_instruction}

{exclude_instruction}

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
    print("🤖 ChatCellAnno: AI Prompt Generated")
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
