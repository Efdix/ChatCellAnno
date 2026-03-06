from PySide6.QtWidgets import (
    QGroupBox, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QLabel,
    QComboBox, QRadioButton, QButtonGroup, QToolButton, QTextEdit
)
from PySide6.QtCore import Qt

def build_data_selection_ui(main_window):
    """Builds Step 1: Data Selection UI"""
    data_group = QGroupBox(main_window.config.T("step1"))
    data_layout = QVBoxLayout()
    
    file_path_layout = QHBoxLayout()
    main_window.path_edit = QLineEdit()
    main_window.path_edit.setPlaceholderText(main_window.config.T("placeholder_path"))
    browse_btn = QPushButton(main_window.config.T("browse"))
    browse_btn.clicked.connect(main_window.browse_file)
    file_path_layout.addWidget(main_window.path_edit)
    file_path_layout.addWidget(browse_btn)
    data_layout.addLayout(file_path_layout)

    # Optional Matrix Row
    data_layout.addWidget(QLabel(main_window.config.T("optional_matrix")))
    matrix_path_layout = QHBoxLayout()
    main_window.matrix_path_edit = QLineEdit()
    main_window.matrix_path_edit.setPlaceholderText(main_window.config.T("matrix_hint"))
    main_window.matrix_path_edit.setToolTip("Optional CSV file: mean expression of key genes per cluster")
    browse_matrix_btn = QPushButton(main_window.config.T("browse"))
    browse_matrix_btn.clicked.connect(main_window.browse_matrix_file)
    matrix_path_layout.addWidget(main_window.matrix_path_edit)
    matrix_path_layout.addWidget(browse_matrix_btn)
    data_layout.addLayout(matrix_path_layout)

    # --- Visual Context (Moved from Step 3 to Step 1) ---
    main_window.visual_box = QGroupBox(main_window.config.T("visual_context"))
    visual_layout = QVBoxLayout()
    
    btn_img_layout = QHBoxLayout()
    main_window.btn_paste_img = QPushButton(main_window.config.T("paste_img"))
    main_window.btn_paste_img.clicked.connect(main_window.paste_image_from_clipboard)
    main_window.btn_clear_img = QPushButton(main_window.config.T("clear_img"))
    main_window.btn_clear_img.clicked.connect(main_window.clear_image)
    btn_img_layout.addWidget(main_window.btn_paste_img)
    btn_img_layout.addWidget(main_window.btn_clear_img)
    
    main_window.img_type_combo = QComboBox()
    main_window.img_type_combo.addItems(["UMAP Plot", "t-SNE Plot", "Spatial Plot", "Heatmap"])
    
    main_window.img_preview = QLabel(main_window.config.T("no_img_loaded"))
    main_window.img_preview.setAlignment(Qt.AlignCenter)
    main_window.img_preview.setStyleSheet("border: 1px dashed #aaa; padding: 10px; color: #888; background: #fdfdfd;")
    main_window.img_preview.setMinimumHeight(100)
    main_window.img_preview.setMaximumHeight(130)
    main_window.img_data = None 
    
    visual_layout.addLayout(btn_img_layout)
    visual_layout.addWidget(main_window.img_type_combo)
    visual_layout.addWidget(main_window.img_preview)

    main_window.visual_box.setLayout(visual_layout)
    data_layout.addWidget(main_window.visual_box)
    # ----------------------------------------------------
    
    data_group.setLayout(data_layout)
    return data_group
