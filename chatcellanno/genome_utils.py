import os
import re

def find_gene_sequences(gene_name, genome_dir="genome"):
    """
    Search for gene name in FAA files.
    Optimized for files where sequence names (IDs) or headers contain the gene symbol.
    """
    all_records = []
    if not os.path.exists(genome_dir):
        return all_records

    # Filter and sort files (Human first, then Mouse, then others)
    faa_files = [f for f in os.listdir(genome_dir) if f.endswith((".faa", ".fasta", ".fa"))]
    def sort_key(name):
        name_l = name.lower()
        if "human" in name_l: return 0
        if "mouse" in name_l: return 1
        return 2
    faa_files.sort(key=sort_key)
    
    # Pre-compile patterns: 
    # 1. Exact ID match (e.g. >PRRX1 ...)
    # 2. Key-value matches (gene=PRRX1)
    # 3. Word match anywhere in header
    symbol_patterns = [
        re.compile(rf"^>{re.escape(gene_name)}(\s|$)", re.I),
        re.compile(rf"\[gene={re.escape(gene_name)}\]", re.I),
        re.compile(rf"Gene[ _]symbol={re.escape(gene_name)}\b", re.I),
        re.compile(rf"gene={re.escape(gene_name)}\b", re.I),
        re.compile(rf"\b{re.escape(gene_name)}\b", re.I)
    ]

    for faa_file in faa_files:
        filepath = os.path.join(genome_dir, faa_file)
        header, seq = _search_in_file(filepath, symbol_patterns)
        if header:
            all_records.append(_format_record(faa_file, header, seq))

    return all_records

def _search_in_file(filepath, patterns):
    found_header = None
    sequence_lines = []
    is_extracting = False
    
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if is_extracting: break
                    if any(p.search(line) for p in patterns):
                        found_header = line
                        is_extracting = True
                elif is_extracting:
                    sequence_lines.append(line)
    except: pass
    return found_header, "".join(sequence_lines)

def _format_record(filename, header, sequence):
    # Regex to strip common bioinformatics suffixes and extensions case-insensitively
    # This will remove .Renamed, .raw, .faa, .fasta etc.
    import re
    cleaned_name = re.sub(r'(\.renamed|\.raw|\.faa|\.fasta|\.fa|\.txt)', '', filename, flags=re.I)
    
    # Capitalize for display (Human, Mouse...)
    species = cleaned_name.replace("_", " ").title()
    return {
        "species": species,
        "header": header,
        "sequence": sequence,
        "filename": filename
    }

def format_mega_style(records, wrap=60):
    """
    Format sequences into MEGA-style alignment view.
    If wrap is None or 0, it produces one long line per species.
    """
    if not records: return ""
    
    ref_seq = records[0]['sequence']
    output = []
    
    # Determine the longest name for alignment
    name_width = max(len(r['species']) for r in records) + 2
    
    # If wrap is None, treat it as one giant block
    effective_wrap = wrap if wrap and wrap > 0 else len(ref_seq)
    
    # Break into chunks
    for start in range(0, len(ref_seq), effective_wrap):
        end = min(start + effective_wrap, len(ref_seq))
        
        # Ruler/Coord
        if wrap:
            output.append(f"{' ' * name_width} {start + 1}")
        
        for i, rec in enumerate(records):
            seq = rec['sequence']
            current_chunk = seq[start:end]
            
            # Display logic
            display_chunk = ""
            if i == 0:
                display_chunk = current_chunk
            else:
                # Compare with reference
                ref_chunk = ref_seq[start:end]
                for j in range(len(current_chunk)):
                    if j < len(ref_chunk) and ref_chunk[j] == current_chunk[j]:
                        display_chunk += "."
                    else:
                        display_chunk += current_chunk[j]
                    
            output.append(f"{rec['species'].ljust(name_width)} {display_chunk}")
        
        if wrap:
            output.append("") # Gap between blocks
        
    return "\n".join(output)

def format_alignment_html(records):
    """
    Generate a beautiful HTML-based alignment view.
    Includes color-coded amino acids and highlighted conservation.
    """
    if not records: return ""
    
    # Simple AA Color Map (Clustal-like)
    color_map = {
        'A': '#80a0f0', 'R': '#f01505', 'N': '#00ff00', 'D': '#c048c0',
        'C': '#f08080', 'Q': '#00ff00', 'E': '#c048c0', 'G': '#f09048',
        'H': '#1505ff', 'I': '#80a0f0', 'L': '#80a0f0', 'K': '#f01505',
        'M': '#80a0f0', 'F': '#80a0f0', 'P': '#ffff00', 'S': '#f09048',
        'T': '#f09048', 'W': '#80a0f0', 'Y': '#1505ff', 'V': '#80a0f0',
        '-': '#e0e0e0'
    }

    ref_seq = records[0]['sequence']
    name_width = max(len(r['species']) for r in records)
    
    html = [
        "<style>",
        "table { border-collapse: collapse; font-family: 'Consolas', monospace; font-size: 13px; }",
        "th { text-align: left; padding: 4px 10px; background: #f2f2f2; color: #666; position: sticky; left: 0; z-index: 1; }",
        ".aa { width: 18px; height: 22px; text-align: center; border: 0.5px solid #eee; }",
        ".match { color: #ccc; }",
        ".coord { font-size: 10px; color: #999; border-bottom: 1px solid #ddd; }",
        "</style>",
        "<table>"
    ]

    # Add Ruler Row (Every 10)
    html.append("<tr><th>Position</th>")
    for j in range(len(ref_seq)):
        if (j + 1) % 10 == 0 or j == 0:
            html.append(f"<td class='coord' colspan='1'>{j+1}</td>")
        else:
            html.append("<td class='coord'></td>")
    html.append("</tr>")

    for i, rec in enumerate(records):
        html.append(f"<tr><th>{rec['species']}</th>")
        seq = rec['sequence']
        
        for j in range(len(ref_seq)):
            char = seq[j] if j < len(seq) else "-"
            
            # Match Logic
            if i > 0 and char != "-" and j < len(ref_seq) and char == ref_seq[j]:
                html.append(f"<td class='aa match'>.</td>")
            else:
                bg = color_map.get(char.upper(), "#ffffff")
                color = "white" if bg != "#ffff00" and bg != "#ffffff" else "black"
                html.append(f"<td class='aa' style='background-color:{bg}; color:{color};'>{char}</td>")
        
        html.append("</tr>")
        
    html.append("</table>")
    return "".join(html)

def align_records(records):
    """
    Perform Multiple Sequence Alignment (MSA) using Biopython.
    If Biopython is not installed, falls back to a simple heuristic.
    """
    if len(records) < 2: return records
    
    # Try Biopython (Star Alignment)
    try:
        return _run_biopython_alignment(records)
    except:
        # Fallback to simple heuristic if biopython is missing or fails
        return _run_heuristic_alignment(records)

def _run_biopython_alignment(records):
    """
    Use Biopython's PairwiseAligner to perform a Star Alignment.
    The first sequence is the center (reference).
    """
    try:
        from Bio import Align
    except ImportError:
        raise Exception("Biopython not installed")

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -1

    ref_rec = records[0]
    ref_seq = ref_rec['sequence']
    
    new_records = [ref_rec.copy()]
    
    for i in range(1, len(records)):
        target_rec = records[i].copy()
        target_seq = target_rec['sequence']
        
        # Get the best alignment
        alignments = aligner.align(ref_seq, target_seq)
        best = alignments[0]
        
        # Extract the aligned sequences (with gaps)
        # best[0] is ref with gaps, best[1] is target with gaps
        aligned_ref = str(best[0])
        aligned_target = str(best[1])
        
        # Re-using the logic:
        # Fix the target sequence length to match the reference sequence 
        # so they can be compared vertically in the grid.
        res_target = ""
        # We walk the alignment. We only care about positions corresponding to the reference
        for r_char, t_char in zip(aligned_ref, aligned_target):
            if r_char != "-":
                res_target += t_char
            else:
                # Target has an insertion; we skip to keep coordinate stability vs Reference
                continue
        
        # Padding if target ended prematurely
        if len(res_target) < len(ref_seq):
            res_target += "-" * (len(ref_seq) - len(res_target))
            
        target_rec['sequence'] = res_target[:len(ref_seq)]
        new_records.append(target_rec)
        
    return new_records

def _run_heuristic_alignment(records):
    import difflib
    ref_rec = records[0]
    ref_seq = ref_rec['sequence']
    new_records = [ref_rec.copy()]
    
    for i in range(1, len(records)):
        target_rec = records[i].copy()
        target_seq = target_rec['sequence']
        s = difflib.SequenceMatcher(None, ref_seq, target_seq)
        new_target = ""
        for tag, i1, i2, j1, j2 in s.get_opcodes():
            if tag == 'equal':
                new_target += target_seq[j1:j2]
            elif tag == 'replace':
                diff_len = (i2-i1) - (j2-j1)
                new_target += target_seq[j1:j2]
                if diff_len > 0: new_target += "-" * diff_len
            elif tag == 'delete':
                new_target += "-" * (i2-i1)
            elif tag == 'insert':
                pass 
        if len(new_target) < len(ref_seq):
            new_target += "-" * (len(ref_seq) - len(new_target))
        target_rec['sequence'] = new_target[:len(ref_seq)]
        new_records.append(target_rec)
    return new_records

def save_fasta_combined(records, output_path):
    import difflib
    ref_rec = records[0]
    ref_seq = ref_rec['sequence']
    new_records = [ref_rec.copy()]
    
    for i in range(1, len(records)):
        target_rec = records[i].copy()
        target_seq = target_rec['sequence']
        s = difflib.SequenceMatcher(None, ref_seq, target_seq)
        new_target = ""
        for tag, i1, i2, j1, j2 in s.get_opcodes():
            if tag == 'equal':
                new_target += target_seq[j1:j2]
            elif tag == 'replace':
                diff_len = (i2-i1) - (j2-j1)
                new_target += target_seq[j1:j2]
                if diff_len > 0: new_target += "-" * diff_len
            elif tag == 'delete':
                new_target += "-" * (i2-i1)
            elif tag == 'insert':
                pass 
        if len(new_target) < len(ref_seq):
            new_target += "-" * (len(ref_seq) - len(new_target))
        target_rec['sequence'] = new_target[:len(ref_seq)]
        new_records.append(target_rec)
    return new_records

def get_fasta_string(records):
    """Generate a combined FASTA string for all records."""
    lines = []
    for r in records:
        lines.append(f">{r['species']} | {r['header'].lstrip('>')}")
        # Break sequence into 80-char lines for standard FASTA
        seq = r['sequence']
        for i in range(0, len(seq), 80):
            lines.append(seq[i:i+80])
        lines.append("") # Gap between records
    return "\n".join(lines)

def save_fasta_combined(records, output_path):
    """Save all records to a single FASTA file."""
    fasta_content = get_fasta_string(records)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(fasta_content)
