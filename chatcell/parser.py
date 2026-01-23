import warnings

def parse_llm_response(response_text: str, cluster_names: list) -> tuple:
    """
    Parses the raw text response from the LLM (Copilot, ChatGPT, DeepSeek, etc.) into dictionaries.
    
    Args:
        response_text: The string copied from the AI Chat.
        cluster_names: A list of cluster names (strings).
        
    Returns:
        tuple: (annotations_dict, extra_info_dict)
               annotations_dict: {cluster_name: annotation}
               extra_info_dict: {cluster_name: extra_info} or empty if simple mode
    """
    # 1. Cleanup Markdown Code blocks
    lines = response_text.strip().split('\n')
    clean_lines = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('```'):
            continue
        clean_lines.append(line)
        
    # 2. Validation
    if len(clean_lines) != len(cluster_names):
        warnings.warn(
            f"⚠️ Length Mismatch! Expected {len(cluster_names)} annotations, got {len(clean_lines)}.\n"
            f"Mapping might be incorrect. Please check the AI output."
        )
        # Handle mismatch gracefully (pad or truncate)
        if len(clean_lines) > len(cluster_names):
            clean_lines = clean_lines[:len(cluster_names)]
        else:
            clean_lines += ["Unknown"] * (len(cluster_names) - len(clean_lines))
            
    # 3. Map
    annotations = {}
    extra_info_map = {}
    
    for i, cluster in enumerate(cluster_names):
        text = clean_lines[i]
        
        # Check for pipe separator
        if '|' in text:
            parts = text.split('|', 1)
            main_anno = parts[0].strip()
            extra = parts[1].strip()
            
            annotations[cluster] = main_anno
            extra_info_map[cluster] = extra
        else:
            annotations[cluster] = text
            # No extra info
            
    return annotations, extra_info_map
