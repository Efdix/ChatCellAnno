"""
功能富集分析模块 (Functional Enrichment Module)

利用 gseapy 库对标记基因列表进行 Over-Representation Analysis (ORA) 富集分析。
此模块为 ChatCellAnno 提供"双引擎验证"能力，通过数据库结果修正 AI 幻觉。

依赖:
    - gseapy
    - pandas
"""

import pandas as pd
import gseapy as gp
import os
import matplotlib.pyplot as plt

def perform_enrichment(gene_dict: dict, species: str = "Human", database_path: str = "", top_term_n: int = 3, out_dir: str = "results/enrichment", is_local: bool = True):
    """
    对每个 Cluster 的基因列表执行 ORA 富集分析（支持本地数据库和 Enrichr 在线数据库）。
    """
    enrichment_summary = {}
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if is_local:
        if not os.path.exists(database_path):
            print(f"❌ Error: Local database file not found at {database_path}")
            return {str(c): {"hints": [f"Error: Database not found"], "table_path": None, "plot_path": None} for c in gene_dict.keys()}
        # Enforce absolute path for gseapy to avoid resolution issues
        database_path = os.path.abspath(database_path)

    mode_str = "Local" if is_local else "Online (Enrichr)"
    print(f"📊 Starting {mode_str} Enrichment Analysis using {database_path}...")

    for cluster, genes in gene_dict.items():
        if isinstance(genes, str):
            genes = [g.strip() for g in genes.split(',')]
        
        # Standardize genes to Uppercase for matching against GMT/Enrichr libraries
        # Most databases (CellMarker, GO, KEGG) use all-caps gene symbols.
        genes = [g.upper() for g in genes if g]
        
        if len(genes) < 3:
            enrichment_summary[str(cluster)] = {"hints": ["Not enough genes"], "table_path": None, "plot_path": None}
            continue

        try:
            enr = gp.enrich(
                gene_list=genes,
                gene_sets=database_path,
                background=None, 
                outdir=None,
                no_plot=True
            )
            
            res_df = enr.results
            if res_df is None or res_df.empty:
                enrichment_summary[str(cluster)] = {"hints": ["No significant terms"], "table_path": None, "plot_path": None}
                continue
                
            # If the only result is an error message or similar (sometimes gseapy does this)
            if "Term" not in res_df.columns:
                enrichment_summary[str(cluster)] = {"hints": ["Invalid results format"], "table_path": None, "plot_path": None}
                continue

            # 1. 保存原始结果表格
            table_filename = f"cluster_{cluster}_enrichment.csv"
            table_path = os.path.join(out_dir, table_filename)
            res_df.to_csv(table_path, index=False)

            # 2. 生成可视化图片
            plot_filename = f"cluster_{cluster}_plot.png"
            plot_path = os.path.join(out_dir, plot_filename)
            
            # 使用 gseapy 的 dotplot 绘图，如果没有结果则跳过
            try:
                # 只针对前 10 个显著的条目绘图
                from gseapy import dotplot
                ax = dotplot(res_df, 
                             column="Adjusted P-value", 
                             x='Combined Score', # 或者用 Combined Score
                             size=10, 
                             top_term=10, 
                             figsize=(8, 6), 
                             title=f"Cluster {cluster} Enrichment",
                             ofname=plot_path)
                plt.close('all') # 释放内存
            except Exception as plot_err:
                print(f"Plotting failed for cluster {cluster}: {plot_err}")
                plot_path = None

            # 3. 提取 AI 提示词 Hint
            res_df = res_df.sort_values(by="Adjusted P-value")
            hints = []
            for i in range(min(top_term_n, len(res_df))):
                row = res_df.iloc[i]
                term = row["Term"]
                pval = row["Adjusted P-value"]
                hints.append(f"{term} (P={pval:.1e})")
            
            enrichment_summary[str(cluster)] = {
                "hints": hints,
                "table_path": os.path.abspath(table_path),
                "plot_path": os.path.abspath(plot_path) if plot_path else None,
                "full_df": res_df.head(10) # 也可以存个简版用于 UI 展示
            }
            
        except Exception as e:
            print(f"⚠️ Enrichment failed for Cluster {cluster}: {e}")
            enrichment_summary[str(cluster)] = {"hints": [f"Error: {str(e)[:20]}"], "table_path": None, "plot_path": None}

    return enrichment_summary
