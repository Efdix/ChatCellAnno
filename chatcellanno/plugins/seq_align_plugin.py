from PySide6.QtWidgets import QAction, QDialog, QVBoxLayout, QTextEdit, QPushButton, QLabel

class SequenceAlignmentPlugin:
    """
    A plugin for basic sequence alignment features.
    """
    def __init__(self, main_window):
        self.main_window = main_window

    def init_ui(self):
        # Add to the gear menu or create a new menu if possible
        # Since main_window has self.gear_menu, let's attach to it for now
        if hasattr(self.main_window, 'gear_menu'):
            align_action = self.main_window.gear_menu.addAction("🧬 Sequence Alignment")
            align_action.triggered.connect(self.show_dialog)
            
    def show_dialog(self):
        dialog = QDialog(self.main_window)
        dialog.setWindowTitle("Sequence Alignment")
        dialog.resize(500, 400)
        
        layout = QVBoxLayout(dialog)
        
        layout.addWidget(QLabel("Paste Sequence 1:"))
        seq1_edit = QTextEdit()
        layout.addWidget(seq1_edit)
        
        layout.addWidget(QLabel("Paste Sequence 2:"))
        seq2_edit = QTextEdit()
        layout.addWidget(seq2_edit)
        
        btn = QPushButton("Align")
        layout.addWidget(btn)
        
        result_label = QLabel("Result will appear here...")
        layout.addWidget(result_label)
        
        def do_align():
            s1 = seq1_edit.toPlainText().strip()
            s2 = seq2_edit.toPlainText().strip()
            if not s1 or not s2:
                result_label.setText("Please provide both sequences.")
                return
            
            # Very basic identical match demonstration
            if s1 == s2:
                result_label.setText("100% Match! (Exact Match)")
            else:
                # simple edit distance/overlap or just a placeholder message for the plugin demo
                result_label.setText(f"Sequences differ. Lengths: {len(s1)} vs {len(s2)}")
                
        btn.clicked.connect(do_align)
        
        dialog.exec()

def register_plugin(main_window):
    return SequenceAlignmentPlugin(main_window)
