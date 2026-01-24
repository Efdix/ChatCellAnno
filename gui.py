import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import pandas as pd
import pyperclip
import sys

# Import core logic from chatcell
try:
    from chatcellanno.core import annotate_cell_types
except ImportError:
    # If running as script, ensure chatcell is in path
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from chatcellanno.core import annotate_cell_types

# Try to import windnd for drag and drop support
try:
    import windnd
    DND_SUPPORT = True
except ImportError:
    DND_SUPPORT = False

class ChatCellApp:
    def __init__(self, root):
        self.root = root
        self.root.title("ChatCell AI - Single Cell Annotation (Zero-API)")
        self.root.geometry("700x850")
        
        self.file_path = tk.StringVar()
        self.species = tk.StringVar(value="Human")
        self.tissue = tk.StringVar(value="PBMC")
        self.top_n = tk.StringVar(value="10")
        self.mode = tk.StringVar(value="concise")
        
        self.setup_ui()
        
        if DND_SUPPORT:
            windnd.hook_dropfiles(self.root, func=self.on_drop)
            
    def setup_ui(self):
        style = ttk.Style()
        style.configure("TLabel", font=("Arial", 10))
        style.configure("TButton", font=("Arial", 10))
        
        padding = {'padx': 10, 'pady': 5}
        
        # File Section
        file_frame = ttk.LabelFrame(self.root, text="Step 1: Select Data File (.csv, .tsv)", padding=(10, 10))
        file_frame.pack(fill="x", **padding)
        
        ttk.Label(file_frame, text="File Path:").grid(row=0, column=0, sticky="w")
        self.file_entry = ttk.Entry(file_frame, textvariable=self.file_path, width=60)
        self.file_entry.grid(row=0, column=1, **padding)
        
        btn_browse = ttk.Button(file_frame, text="Browse", command=self.browse_file)
        btn_browse.grid(row=0, column=2, **padding)
        
        if DND_SUPPORT:
            ttk.Label(file_frame, text="(Drag and drop files onto this window)", foreground="gray").grid(row=1, column=1, sticky="w")
        else:
            ttk.Label(file_frame, text="(Install 'windnd' for drag-and-drop support)", foreground="gray").grid(row=1, column=1, sticky="w")
            
        # Parameters Section
        param_frame = ttk.LabelFrame(self.root, text="Step 2: Configuration", padding=(10, 10))
        param_frame.pack(fill="x", **padding)
        
        # Species
        ttk.Label(param_frame, text="Species:").grid(row=0, column=0, sticky="w")
        ttk.Entry(param_frame, textvariable=self.species).grid(row=0, column=1, sticky="w", **padding)
        
        # Tissue
        ttk.Label(param_frame, text="Tissue:").grid(row=1, column=0, sticky="w")
        ttk.Entry(param_frame, textvariable=self.tissue).grid(row=1, column=1, sticky="w", **padding)
        
        # Top N
        ttk.Label(param_frame, text="Top Markers:").grid(row=2, column=0, sticky="w")
        ttk.Entry(param_frame, textvariable=self.top_n).grid(row=2, column=1, sticky="w", **padding)
        
        # Mode
        ttk.Label(param_frame, text="Output Mode:").grid(row=3, column=0, sticky="w")
        mode_radio_frame = ttk.Frame(param_frame)
        mode_radio_frame.grid(row=3, column=1, sticky="w", **padding)
        
        ttk.Radiobutton(mode_radio_frame, text="Concise (Cell Type only)", variable=self.mode, value="concise").pack(side="left", padx=5)
        ttk.Radiobutton(mode_radio_frame, text="Detailed (Type | Markers | Functions | Ranks)", variable=self.mode, value="detailed").pack(side="left", padx=5)
        
        # Action Section
        btn_generate = ttk.Button(self.root, text="Generate Prompt", command=self.generate_prompt)
        btn_generate.pack(pady=10)
        
        # Output Section
        output_frame = ttk.LabelFrame(self.root, text="Step 3: Result (Copy to LLM)", padding=(10, 10))
        output_frame.pack(fill="both", expand=True, **padding)
        
        self.output_text = tk.Text(output_frame, height=15, font=("Consolas", 10))
        self.output_text.pack(fill="both", expand=True, side="left")
        
        scrollbar = ttk.Scrollbar(output_frame, command=self.output_text.yview)
        scrollbar.pack(side="right", fill="y")
        self.output_text.config(yscrollcommand=scrollbar.set)
        
        btn_copy = ttk.Button(self.root, text="Copy to Clipboard", command=self.copy_to_clipboard)
        btn_copy.pack(pady=10)
        
        # Footer
        footer = ttk.Label(self.root, text="Human-in-the-loop AI annotation for Scanpy", foreground="gray")
        footer.pack(side="bottom", pady=5)

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Select Data File",
            filetypes=[("CSV/TSV files", "*.csv;*.tsv;*.txt"), ("All files", "*.*")]
        )
        if filename:
            self.file_path.set(filename)

    def on_drop(self, files):
        if files:
            # windnd returns binary strings on some versions or list of bytes
            path = files[0].decode('gbk') if isinstance(files[0], bytes) else files[0]
            self.file_path.set(path)

    def generate_prompt(self):
        path = self.file_entry.get().strip()
        if not path:
            messagebox.showerror("Error", "Please select a file first.")
            return
        
        if not os.path.exists(path):
            messagebox.showerror("Error", "File does not exist.")
            return

        try:
            self.output_text.delete(1.0, tk.END)
            self.output_text.insert(tk.END, "Loading data and generating prompt... Please wait.\n")
            self.root.update_idletasks()
            
            top_n = int(self.top_n.get())
            
            # Use core.annotate_cell_types
            prompt = annotate_cell_types(
                marker_file=path,
                species=self.species.get(),
                tissue=self.tissue.get(),
                top_n=top_n,
                mode=self.mode.get(),
                step="generate"
            )
            
            self.output_text.delete(1.0, tk.END)
            self.output_text.insert(tk.END, prompt)
            
            # Also copy to clipboard automatically as per library design
            pyperclip.copy(prompt)
            messagebox.showinfo("Success", "Prompt generated and copied to clipboard!")
            
        except Exception as e:
            self.output_text.delete(1.0, tk.END)
            messagebox.showerror("Error", f"Failed to generate prompt:\n{str(e)}")

    def copy_to_clipboard(self):
        content = self.output_text.get(1.0, tk.END).strip()
        if content:
            pyperclip.copy(content)
            messagebox.showinfo("Copied", "Result copied to clipboard.")

if __name__ == "__main__":
    root = tk.Tk()
    app = ChatCellApp(root)
    root.mainloop()
