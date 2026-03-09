from PySide6.QtWidgets import (
    QGroupBox, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QLabel,
    QComboBox, QRadioButton, QButtonGroup, QToolButton, QCheckBox, QWidget
)
from PySide6.QtCore import Qt

def build_data_selection_ui(main_window):
    """
    构建步骤 1：数据选择与多模态输入界面
    """
    data_group = QGroupBox(main_window.config.T("step1"))
    data_layout = QVBoxLayout()
    
    # 1.1 核心 Marker 文件选择 (Core Marker File Selection)
    data_layout.addWidget(QLabel(main_window.config.T("marker_list_label")))
    file_path_layout = QHBoxLayout()
    main_window.path_edit = QLineEdit()
    main_window.path_edit.setPlaceholderText(main_window.config.T("placeholder_path"))
    browse_btn = QPushButton(main_window.config.T("browse"))
    browse_btn.clicked.connect(main_window.browse_file)
    file_path_layout.addWidget(main_window.path_edit)
    file_path_layout.addWidget(browse_btn)
    data_layout.addLayout(file_path_layout)
    
    # 1.2 可选：功能富集分析 (Enrichment Analysis inside Step 1)
    main_window.use_enrichment_cb = QCheckBox(main_window.config.T("enrich_assist_label"))
    main_window.use_enrichment_cb.setChecked(False)
    data_layout.addWidget(main_window.use_enrichment_cb)
    
    main_window.enrich_widget = QWidget()
    enrich_layout = QVBoxLayout(main_window.enrich_widget)
    enrich_layout.setContentsMargins(10, 0, 0, 0) # Indent slightly
    
    enrich_layout.addWidget(QLabel(main_window.config.T("enrich_source")))
    source_layout = QHBoxLayout()
    main_window.rb_online = QRadioButton(main_window.config.T("online"))
    main_window.rb_local = QRadioButton(main_window.config.T("local"))
    main_window.rb_online.setChecked(True) 
    main_window.source_group_enrich = QButtonGroup(main_window)
    main_window.source_group_enrich.addButton(main_window.rb_online)
    main_window.source_group_enrich.addButton(main_window.rb_local)
    source_layout.addWidget(main_window.rb_online)
    source_layout.addWidget(main_window.rb_local)
    enrich_layout.addLayout(source_layout)

    main_window.enrich_db_combo = QComboBox()
    db_select_layout = QHBoxLayout()
    db_select_layout.addWidget(main_window.enrich_db_combo)
    
    main_window.load_db_btn = QToolButton()
    main_window.load_db_btn.setText("📁")
    main_window.load_db_btn.setToolTip("Load custom database folder")
    main_window.load_db_btn.clicked.connect(main_window.browse_db_folder)
    main_window.load_db_btn.hide()
    db_select_layout.addWidget(main_window.load_db_btn)
    
    enrich_layout.addWidget(QLabel("Select GeneSet Database:"))
    enrich_layout.addLayout(db_select_layout)
    
    main_window.custom_db_label = QLabel("")
    main_window.custom_db_label.setStyleSheet("color: #666; font-size: 10px;")
    enrich_layout.addWidget(main_window.custom_db_label)

    main_window.rb_local.toggled.connect(main_window.update_enrich_db_list)
    main_window.rb_online.toggled.connect(main_window.update_enrich_db_list)
    main_window.update_enrich_db_list()
    
    if main_window.config.last_db_path:
        main_window.load_databases_from_folder(main_window.config.last_db_path)

    main_window.run_enrich_btn = QPushButton(main_window.config.T("gen_btn_enrich"))
    main_window.run_enrich_btn.setStyleSheet("background-color: #673AB7; color: white; font-weight: bold; height: 35px;")
    main_window.run_enrich_btn.clicked.connect(main_window.start_enrichment_thread)
    enrich_layout.addWidget(main_window.run_enrich_btn)
    
    main_window.enrich_widget.setVisible(False)
    main_window.use_enrichment_cb.toggled.connect(main_window.enrich_widget.setVisible)
    data_layout.addWidget(main_window.enrich_widget)

    # 1.3 可选：表达矩阵输入 (Optional Matrix Input)
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

    # 1.4 视觉上下文模块 (Visual Context -> Text/Image analysis)
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
    
    data_group.setLayout(data_layout)
    return data_group

