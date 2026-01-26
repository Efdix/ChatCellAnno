"""
图形用户界面 (GUI) 模块

基于 Tkinter 构建的桌面应用程序入口。
它为 chatcellanno 库提供了一个可视化操作界面，支持文件选择（或拖拽）、参数配置、一键生成 Prompt 并复制等功能。
主要使用了 `core.annotate_cell_types` 来驱动业务逻辑。

依赖项:
- tkinter (Python 内置 GUI 库)
- windnd (Windows 拖拽支持，可选)
- chatcellanno (本项目核心库)
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import pandas as pd
import pyperclip
import sys

# 尝试导入核心逻辑
try:
    from chatcellanno.core import annotate_cell_types
except ImportError:
    # 如果直接作为脚本运行 (python gui.py)，可能需要添加路径才能找到 chatcellanno 包
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from chatcellanno.core import annotate_cell_types

# 尝试导入 windnd 以支持拖拽文件 (Windows Only)
try:
    import windnd
    DND_SUPPORT = True
except ImportError:
    DND_SUPPORT = False

class ChatCellApp:
    def __init__(self, root):
        self.root = root
        self.root.title("ChatCell AI - Single Cell Annotation (Zero-API)")
        self.root.geometry("800x950") # 增加高度以容纳新功能
        
        # 定义 UI 绑定的变量 (UI Variables)
        self.file_path = tk.StringVar()
        self.species = tk.StringVar(value="Human")
        self.tissue = tk.StringVar(value="PBMC")
        self.top_n = tk.StringVar(value="10")
        self.mode = tk.StringVar(value="concise")
        self.source = tk.StringVar(value="scanpy") # 新增 Source 变量
        
        # 初始化界面布局
        self.setup_ui()
        
        # 注册拖拽事件
        if DND_SUPPORT:
            windnd.hook_dropfiles(self.root, func=self.on_drop)
            
    def setup_ui(self):
        """配置并放置所有 UI 组件"""
        style = ttk.Style()
        style.configure("TLabel", font=("Arial", 10))
        style.configure("TButton", font=("Arial", 10))
        
        padding = {'padx': 10, 'pady': 5}
        
        # --- 滚动区域设置 (Main Scrollable Area) ---
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill="both", expand=True)
        
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Save canvas for scrolling
        self.canvas = canvas
        
        # Bind MouseWheel for scrolling anywhere
        # Windows: <MouseWheel>, Linux: Button-4/5
        self.root.bind_all("<MouseWheel>", self._on_mousewheel)
        
        # ----------------------------------------
        
        # 1. 区域：文件选择 (File Section)
        file_frame = ttk.LabelFrame(scrollable_frame, text="Step 1: Data & Source", padding=(10, 10))
        file_frame.pack(fill="x", **padding)
        
        # 文件路径
        ttk.Label(file_frame, text="File Path:").grid(row=0, column=0, sticky="w")
        self.file_entry = ttk.Entry(file_frame, textvariable=self.file_path, width=50)
        self.file_entry.grid(row=0, column=1, **padding)
        btn_browse = ttk.Button(file_frame, text="Browse", command=self.browse_file)
        btn_browse.grid(row=0, column=2, **padding)

        # 数据来源 (Source)
        ttk.Label(file_frame, text="Source Type:").grid(row=1, column=0, sticky="w")
        source_frame = ttk.Frame(file_frame)
        source_frame.grid(row=1, column=1, sticky="w", **padding)
        ttk.Radiobutton(source_frame, text="Scanpy (Python)", variable=self.source, value="scanpy").pack(side="left", padx=5)
        ttk.Radiobutton(source_frame, text="Seurat (R)", variable=self.source, value="seurat").pack(side="left", padx=5)
        
        if DND_SUPPORT:
            ttk.Label(file_frame, text="(Drag & Drop supported)", foreground="gray").grid(row=2, column=1, sticky="w")
        else:
            ttk.Label(file_frame, text="(Install 'windnd' for drag-and-drop)", foreground="gray").grid(row=2, column=1, sticky="w")

        # 2. 区域：参数配置 (Parameters Section)
        param_frame = ttk.LabelFrame(scrollable_frame, text="Step 2: Configuration", padding=(10, 10))
        param_frame.pack(fill="x", **padding)
        
        # Species, Tissue, Top N
        ttk.Label(param_frame, text="Species:").grid(row=0, column=0, sticky="w")
        ttk.Entry(param_frame, textvariable=self.species).grid(row=0, column=1, sticky="w", **padding)
        
        ttk.Label(param_frame, text="Tissue:").grid(row=1, column=0, sticky="w")
        ttk.Entry(param_frame, textvariable=self.tissue).grid(row=1, column=1, sticky="w", **padding)
        
        ttk.Label(param_frame, text="Top Markers:").grid(row=2, column=0, sticky="w")
        ttk.Entry(param_frame, textvariable=self.top_n).grid(row=2, column=1, sticky="w", **padding)
        
        # Mode
        ttk.Label(param_frame, text="Output Mode:").grid(row=3, column=0, sticky="w")
        mode_radio_frame = ttk.Frame(param_frame)
        mode_radio_frame.grid(row=3, column=1, sticky="w", **padding)
        ttk.Radiobutton(mode_radio_frame, text="Concise (Type only)", variable=self.mode, value="concise").pack(side="left", padx=5)
        ttk.Radiobutton(mode_radio_frame, text="Detailed", variable=self.mode, value="detailed").pack(side="left", padx=5)
        
        # 3. 区域：生成 Prompt (Generate Section)
        gen_frame = ttk.LabelFrame(scrollable_frame, text="Step 3: Generate Prompt (Copy to LLM)", padding=(10, 10))
        gen_frame.pack(fill="both", expand=True, **padding)
        
        ttk.Button(gen_frame, text="Generate Prompt", command=self.generate_prompt).pack(anchor="w", pady=5)
        
        self.output_text = tk.Text(gen_frame, height=8, font=("Consolas", 9))
        self.output_text.pack(fill="both", expand=True)
        ttk.Button(gen_frame, text="Copy Prompt", command=self.copy_to_clipboard).pack(pady=5)
        
        # 4. 区域：解析结果 (Parse Section)
        parse_frame = ttk.LabelFrame(scrollable_frame, text="Step 4: Parse AI Response & Generate Code", padding=(10, 10))
        parse_frame.pack(fill="both", expand=True, **padding)
        
        ttk.Label(parse_frame, text="Paste AI Response here:").pack(anchor="w")
        self.response_text = tk.Text(parse_frame, height=6, font=("Consolas", 9))
        self.response_text.pack(fill="x", pady=5)
        
        ttk.Button(parse_frame, text="Process & Generate Code", command=self.process_response).pack(pady=5)
        
        ttk.Label(parse_frame, text="Generated Annotation Code:").pack(anchor="w")
        self.code_output = tk.Text(parse_frame, height=10, font=("Consolas", 9), bg="#f0f0f0")
        self.code_output.pack(fill="both", expand=True, pady=5)
        ttk.Button(parse_frame, text="Copy Code", command=self.copy_code).pack(pady=5)

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Select Data File",
            filetypes=[("CSV/TSV files", "*.csv;*.tsv;*.txt"), ("All files", "*.*")]
        )
        if filename:
            self.file_path.set(filename)

    def on_drop(self, files):
        """处理文件拖拽事件"""
        if files:
            path = files[0].decode('gbk') if isinstance(files[0], bytes) else files[0]
            self.file_path.set(path)

    def generate_prompt(self):
        """Step 3: 生成 Prompt"""
        path = self.file_entry.get().strip()
        if not path:
            messagebox.showerror("Error", "Please select a file first.")
            return
        
        if not os.path.exists(path):
            messagebox.showerror("Error", "File does not exist.")
            return

        try:
            self.output_text.delete(1.0, tk.END)
            self.output_text.insert(tk.END, "Processing...\n")
            self.root.update_idletasks()
            
            top_n = int(self.top_n.get())
            source_type = self.source.get()
            
            # 使用 core.annotate_cell_types 生成 Prompt
            prompt = annotate_cell_types(
                marker_file=path,
                species=self.species.get(),
                tissue=self.tissue.get(),
                top_n=top_n,
                mode=self.mode.get(),
                step="generate",
                source=source_type # 传递 source
            )
            
            self.output_text.delete(1.0, tk.END)
            self.output_text.insert(tk.END, prompt)
            
            pyperclip.copy(prompt)
            messagebox.showinfo("Success", "Prompt generated and copied to clipboard!")
            
        except Exception as e:
            self.output_text.delete(1.0, tk.END)
            messagebox.showerror("Error", f"Failed to generate prompt:\n{str(e)}")

    def copy_to_clipboard(self):
        content = self.output_text.get(1.0, tk.END).strip()
        if content:
            pyperclip.copy(content)
            messagebox.showinfo("Copied", "Prompt copied.")

    def process_response(self):
        """Step 4: 解析回复并生成代码"""
        path = self.file_entry.get().strip()
        if not path:
             messagebox.showerror("Error", "Marker file path missing (needed for cluster alignment).")
             return

        ai_response = self.response_text.get(1.0, tk.END).strip()
        if not ai_response:
            messagebox.showerror("Error", "Please paste AI response first.")
            return
            
        try:
            top_n = int(self.top_n.get())
            source_type = self.source.get()
            
            # 调用 core 解析
            result_df, code_snippet = annotate_cell_types(
                marker_file=path,
                step="parse",
                response_text=ai_response,
                top_n=top_n,
                source=source_type
            )
            
            # 显示代码
            self.code_output.delete(1.0, tk.END)
            self.code_output.insert(tk.END, code_snippet)
            
            messagebox.showinfo("Success", "Code generated successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to parse:\n{str(e)}")
            
    def copy_code(self):
        content = self.code_output.get(1.0, tk.END).strip()
        if content:
            pyperclip.copy(content)
            messagebox.showinfo("Copied", "Code copied.")

    def _on_mousewheel(self, event):
        """
        处理鼠标滚轮事件 (Handle mouse wheel scrolling)。
        
        该方法绑定到全局 ("<MouseWheel>")，因此鼠标在应用任何位置滚动时，
        都会触发主 Canvas 的滚动，解决了鼠标必须悬停在特定区域才能滚动的问题。
        
        Windows 系统下 event.delta 也是 120 的倍数。
        -1 * (delta/120) 用于将滚动方向转换为 Tkinter yview_scroll 所需的单位。
        """
        # Windows uses delta, dividing by 120 is standard
        # A negative value means scrolling down
        self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")

if __name__ == "__main__":
    root = tk.Tk()
    app = ChatCellApp(root)
    root.mainloop()
