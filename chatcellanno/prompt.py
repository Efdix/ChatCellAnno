"""
提示词生成模块 (Prompt Engineering Module)

负责构建高质量的 Prompt，指导 LLM 进行细胞注释。
包含 Role-Play、Few-Shot Learning 和 Chain-of-Thought (CoT) 等策略。
"""

import pyperclip
import platform
from typing import Dict, List, Optional

class PromptBuilder:
    """
    用于构建复杂 LLM 提示词的建造者类。
    支持链式调用，分模块组装 Prompt。
    """
    
    def __init__(self):
        self._parts = []
    
    def add_role(self, role_desc: str) -> 'PromptBuilder':
        """添加角色设定"""
        self._parts.append(f"## Role\n{role_desc}\n")
        return self
    
    def add_task(self, task_desc: str) -> 'PromptBuilder':
        """添加任务描述"""
        self._parts.append(f"## Task\n{task_desc}\n")
        return self
    
    def add_constraints(self, constraints: List[str]) -> 'PromptBuilder':
        """添加约束条件"""
        cons_str = "\n".join([f"- {c}" for c in constraints])
        self._parts.append(f"## Constraints\n{cons_str}\n")
        return self
    
    def add_data(self, data_title: str, content: str) -> 'PromptBuilder':
        """添加核心数据"""
        self._parts.append(f"## Data ({data_title})\n```text\n{content}\n```\n")
        return self
    
    def add_output_format(self, format_desc: str) -> 'PromptBuilder':
        """添加输出格式要求"""
        self._parts.append(f"## Output Format\n{format_desc}\n")
        return self
        
    def build(self) -> str:
        """生成最终的提示词字符串"""
        return "\n".join(self._parts)


def copy_to_clipboard(text: str) -> bool:
    """
    尝试将文本复制到系统剪贴板。
    """
    try:
        pyperclip.copy(text)
        return True
    except Exception as e:
        print(f"Warning: Could not copy to clipboard automatically: {e}")
        return False

def generate_annotation_prompt(
    markers: Dict[str, str],
    species: str = "Human",
    tissue: str = "PBMC",
    mode: str = "concise",
    exclude_types: str = "",
    enrichment_hints: Optional[Dict[str, List[str]]] = None,
    visual_context: Optional[str] = None,
    expression_matrix: Optional[str] = None,
    auto_copy: bool = False
) -> str:
    """
    构建专业的单细胞注释 Prompt。
    
    参数:
        visual_context: 视觉上下文描述
        expression_matrix: 含有某些关键基因在各个 cluster 中表达量的矩阵 (Markdown 格式)
    """
    
    # 1. 准备数据块
    data_lines = []
    
    # Check if enrichment_hints is actually a dict and has content
    use_hints = isinstance(enrichment_hints, dict) and len(enrichment_hints) > 0

    for cluster_id, gene_list_str in markers.items():
        base_line = f"Cluster {cluster_id}: {gene_list_str}"
        
        # 注入富集分析线索 (RAG - Retrieval Augmented Generation)
        # Handle both string keys and int keys for cluster_id
        c_key_str = str(cluster_id)
        c_key_int = int(cluster_id) if str(cluster_id).isdigit() else cluster_id

        hints = []
        if use_hints:
            if c_key_str in enrichment_hints:
                hints = enrichment_hints[c_key_str]
            elif c_key_int in enrichment_hints:
                hints = enrichment_hints[c_key_int]
        
        if hints:
            # Join hints, limit to top 5 to avoid token overflow
            hints_text = "; ".join(hints[:5]) 
            base_line += f"\n   [Reference Hints]: {hints_text}"
            
        data_lines.append(base_line)
    
    data_content = "\n".join(data_lines)
    
    # 2. 构建 Prompt
    builder = PromptBuilder()
    
    # Role
    builder.add_role(
        f"You are an expert Cell Biologist and Bioinformatician specializing in {species} {tissue} single-cell RNA-seq analysis."
    )
    
    # Task
    task_desc = f"Identify the cell types for the following clusters based on their top marker genes.\nThe tissue of origin is **{species} {tissue}**."
    
    if visual_context:
        task_desc += f"\n\n**Visual Context Provided**: I have attached an image of the {visual_context} plot. Please examine the spatial relationships and clustering patterns in the image to cross-validate your annotations (e.g., adjacent clusters likely share lineage)."
        
    builder.add_task(task_desc)
    
    # Constraints
    constraints = [
        "Use standard, accepted cell ontology nomenclature.",
        "Be specific (e.g., 'CD8+ T cell' instead of just 'T cell' if markers support it).",
        "If markers are ambiguous or low quality, label as 'Unknown' or 'Low Quality'.",
    ]
    
    if exclude_types:
        constraints.append(f"Do NOT use these cell types: {exclude_types}.")
        
    if mode == "concise":
        constraints.append("Keep the reasoning brief.")
    else:
        constraints.append("Provide a short biological reasoning for each cell type based on the markers.")

    builder.add_constraints(constraints)
    
    # Output Format
    # 强制要求用户先输出摘要表格，再进行详细推理 (Enforce Summary Table followed by Detailed Analysis)
    if mode == "detailed":
        format_desc = (
            "Please structure your response as follows:\n\n"
            "### 1. Summary Table\n"
            "Provide a Markdown table with columns: | Cluster | Cell Type | Reasoning Summary |\n\n"
            "### 2. Comprehensive Analysis (per Cluster)\n"
            "For EACH cluster, provide a detailed paragraph (4-6 sentences) explaining the inference. "
            "Your explanation MUST collaboratively synthesize all available evidence:\n"
            "- How the **Top Marker Genes** define the cell identity.\n"
            "- How the **Reference Hints (Enrichment)** confirm the biological process.\n"
            "- How the **Expression Matrix** (if provided) quantitatively supports the cluster's uniqueness.\n"
            "- How the **Visual Patterns (UMAP/t-SNE)** (if provided) reflect the cluster's spatial or lineage relationships.\n"
            "Finally, explain why this identification makes sense in the context of **{species} {tissue}**."
        ).format(species=species, tissue=tissue)
    else:
        format_desc = (
            "Output ONLY a Markdown table with the following columns:\n"
            "| Cluster | Cell Type | Reasoning |\n"
            "Do not include any other text before or after the table."
        )
    builder.add_output_format(format_desc)
    
    # Data
    builder.add_data("Top Marker Genes per Cluster", data_content)
    
    if expression_matrix:
        matrix_description = (
            "The following matrix shows the mean expression of specific genes of interest across each cluster. "
            "Please use this expression data as additional evidence to refine your cell type identification. "
            "Rows represent clusters and columns represent genes of interest."
        )
        builder.add_data("Gene Expression Matrix (Additional Evidence)", f"{matrix_description}\n\n{expression_matrix}")
    
    # 3. 完成组装并处理输出 (Phase 3: Final execution & clipboard)
    prompt = builder.build()
    
    # 打印一些状态信息到控制台 (Print status to console)
    print("=" * 80)
    print("🤖 ChatCellAnno: AI Prompt Generated via PromptBuilder")
    print("=" * 80)
    
    # 如果开启自动复制，将结果存入系统剪贴板 (Auto-copy to system clipboard if enabled)
    if auto_copy:
        success = copy_to_clipboard(prompt)
        if success:
            print("✅ Prompt has been COPIED to your clipboard! (已复制到剪贴板)")
            print("👉 Go to your AI Chat (Copilot, DeepSeek, ChatGPT) and press Ctrl+V (Paste), then Enter.")
        else:
            # 兼容模式：如果剪贴板调用失败，则在控制台打印全文 (Fallback: print full prompt)
            print("📋 Please copy the prompt below manually:")
            print(prompt)
    else:
        print(prompt)
        
    print("=" * 80)
    
    return prompt
