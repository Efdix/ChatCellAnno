import sys
import os
import pyperclip
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QLabel, QLineEdit, QPushButton, QRadioButton, QButtonGroup, 
    QGroupBox, QTextEdit, QFileDialog, QMessageBox, QScrollArea,
    QSplitter, QComboBox, QDialog, QFormLayout, QToolButton,
    QListWidget, QListWidgetItem, QAbstractItemView, QMenu, QWidgetAction,
    QTabWidget, QProgressBar
)
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtCore import Qt, QUrl
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebEngineCore import QWebEngineProfile, QWebEnginePage, QWebEngineSettings

# Relative imports within the package structure
# Assuming this file is at chatcellanno/gui/main_window.py
# So we go up to chatcellanno then down to config
from ..config import ConfigManager
from .workers import EnrichmentWorker

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.config = ConfigManager()
        
        # Local Database Path Management
        if hasattr(sys, "_MEIPASS"):
            base_dir = sys._MEIPASS
            # If internal, we might have bundled database in root or chatcellanno/database
        else:
            base_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
            
        # Fallback to script dir if database not found next to executable
        db_check = os.path.join(base_dir, "database")
        if not os.path.exists(db_check):
             # Try going up one level from this file location (chatcellanno/gui/.. -> chatcellanno)
             # Actually, if running from source, base_dir usually is the root 'ChatCellAnno'
             # Let's trust the logic that worked in gui.py, but adapted if __file__ is inside chatcellanno/gui
             base_dir_pkg = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
             if os.path.exists(os.path.join(base_dir_pkg, "database")):
                 base_dir = base_dir_pkg

        self.db_dir = os.path.join(base_dir, "database")
        self.db_path_map = {} 

        self.online_db_list = [
            "GO_Biological_Process_2025",
            "GO_Cellular_Component_2025",
            "GO_Molecular_Function_2025",
            "KEGG_2026",
            "CellMarker_2024",
            "PanglaoDB_Augmented_2021"
        ]
        
        # We don't need load_config here as ConfigManager does it on init
        
        # Configure WebEngine Profile for Persistence using ConfigManager path
        self.web_profile = QWebEngineProfile("ChatCellAnnoProfile", self)
        self.web_profile.setPersistentStoragePath(self.config.config_paths["web_storage"])
        self.web_profile.setCachePath(self.config.config_paths["web_storage"])
        
        # State variables
        self.current_marker_file = ""
        self.enrichment_data = None
        self.all_possible_services = self.config.chat_services.copy()
        
        self.setWindowTitle(self.config.T("title"))
        if os.path.exists("app_icon.ico"):
            self.setWindowIcon(QIcon("app_icon.ico"))

        # Global stylesheet
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f5f5f7; 
            }
            QWidget {
                font-family: "Segoe UI", "Microsoft YaHei", sans-serif;
                font-size: 10pt;
            }
            QGroupBox {
                background-color: white;
                border: 1px solid #dcdcdc;
                border-radius: 8px;
                margin-top: 10px;
                padding-top: 15px; 
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 0 5px;
                color: #333;
                font-weight: bold;
                font-size: 15px;
            }
            QPushButton {
                background-color: #007aff; 
                color: white; 
                border-radius: 5px; 
                padding: 6px 12px;
                border: none;
            }
            QPushButton:hover {
                background-color: #0060c0;
            }
            QPushButton:pressed {
                background-color: #004080;
            }
            QLineEdit, QTextEdit, QComboBox {
                border: 1px solid #ccc;
                border-radius: 4px;
                padding: 4px;
                background-color: white;
                selection-background-color: #007aff;
            }
            QLineEdit:focus, QTextEdit:focus, QComboBox:focus {
                border: 1px solid #007aff;
            }
            QListWidget {
                border: none;
                background-color: transparent;
            }
            QListWidget::item {
                background-color: white;
                border: 1px solid #ddd;
                border-radius: 4px;
                margin: 2px;
                padding: 5px;
                color: #333;
            }
            QListWidget::item:hover {
                 background-color: #f0f8ff;
                 border-color: #007aff;
            }
            QScrollArea {
                border: none;
                background-color: transparent;
            }
        """)

        self.setup_ui()
        self.setAcceptDrops(True) 
        
        # Set stretch factors to equal values (1:1)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 1)

    def setup_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)

        self.splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(self.splitter)

        # --- Left Panel: Controls ---
        left_container = QWidget()
        left_layout = QVBoxLayout(left_container)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        left_layout.addWidget(scroll)
        
        control_widget = QWidget()
        control_layout = QVBoxLayout(control_widget)
        scroll.setWidget(control_widget)

        # 1. Data Selection
        data_group = QGroupBox(self.config.T("step1"))
        data_layout = QVBoxLayout()
        
        file_path_layout = QHBoxLayout()
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText(self.config.T("placeholder_path"))
        browse_btn = QPushButton(self.config.T("browse"))
        browse_btn.clicked.connect(self.browse_file)
        file_path_layout.addWidget(self.path_edit)
        file_path_layout.addWidget(browse_btn)
        data_layout.addLayout(file_path_layout)
        data_group.setLayout(data_layout)
        control_layout.addWidget(data_group)

        # 2. Enrichment Analysis
        enrich_group = QGroupBox(self.config.T("step2"))
        enrich_layout = QVBoxLayout()
        
        enrich_layout.addWidget(QLabel(self.config.T("enrich_source")))
        source_layout = QHBoxLayout()
        self.rb_online = QRadioButton(self.config.T("online"))
        self.rb_local = QRadioButton(self.config.T("local"))
        self.rb_online.setChecked(True) 
        self.source_group_enrich = QButtonGroup(self)
        self.source_group_enrich.addButton(self.rb_online)
        self.source_group_enrich.addButton(self.rb_local)
        source_layout.addWidget(self.rb_online)
        source_layout.addWidget(self.rb_local)
        enrich_layout.addLayout(source_layout)

        self.enrich_db_combo = QComboBox()
        db_select_layout = QHBoxLayout()
        db_select_layout.addWidget(self.enrich_db_combo)
        
        self.load_db_btn = QToolButton()
        self.load_db_btn.setText("📁")
        self.load_db_btn.setToolTip("Load custom database folder")
        self.load_db_btn.clicked.connect(self.browse_db_folder)
        self.load_db_btn.hide()
        db_select_layout.addWidget(self.load_db_btn)
        
        enrich_layout.addWidget(QLabel("Select GeneSet Database:"))
        enrich_layout.addLayout(db_select_layout)
        
        self.custom_db_label = QLabel("")
        self.custom_db_label.setStyleSheet("color: #666; font-size: 10px;")
        enrich_layout.addWidget(self.custom_db_label)

        self.rb_local.toggled.connect(self.update_enrich_db_list)
        self.rb_online.toggled.connect(self.update_enrich_db_list)
        self.update_enrich_db_list()

        self.run_enrich_btn = QPushButton(self.config.T("gen_btn_enrich"))
        self.run_enrich_btn.setStyleSheet("background-color: #673AB7; color: white; font-weight: bold; height: 35px;")
        self.run_enrich_btn.clicked.connect(self.start_enrichment_thread)
        enrich_layout.addWidget(self.run_enrich_btn)

        enrich_group.setLayout(enrich_layout)
        control_layout.addWidget(enrich_group)

        # 3. Prompt Configuration
        gen_group = QGroupBox(self.config.T("step3"))
        gen_layout = QVBoxLayout()
        
        self.species_edit = QLineEdit(self.config.session_params.get("species", "Human"))
        self.tissue_edit = QLineEdit(self.config.session_params.get("tissue", "PBMC"))
        self.top_n_edit = QLineEdit(self.config.session_params.get("top_n", "10"))
        
        self.species_edit.textChanged.connect(self.save_config_state)
        self.tissue_edit.textChanged.connect(self.save_config_state)
        self.top_n_edit.textChanged.connect(self.save_config_state)

        gen_layout.addWidget(QLabel(self.config.T("species")))
        gen_layout.addWidget(self.species_edit)
        gen_layout.addWidget(QLabel(self.config.T("tissue")))
        gen_layout.addWidget(self.tissue_edit)
        gen_layout.addWidget(QLabel(self.config.T("top_n")))
        gen_layout.addWidget(self.top_n_edit)
        
        self.exclude_edit = QLineEdit(self.config.session_params.get("exclude", ""))
        self.exclude_edit.setPlaceholderText("e.g. Neuron, Astrocyte")
        self.exclude_edit.textChanged.connect(self.save_config_state)
        
        gen_layout.addWidget(QLabel(self.config.T("exclude")))
        gen_layout.addWidget(self.exclude_edit)

        gen_layout.addWidget(QLabel(self.config.T("prompt_mode")))
        self.mode_group = QButtonGroup(self)
        mode_layout = QHBoxLayout()
        self.rb_concise = QRadioButton(self.config.T("concise"))
        self.rb_detailed = QRadioButton(self.config.T("detailed"))
        
        if self.config.session_params.get("prompt_mode") == "detailed":
            self.rb_detailed.setChecked(True)
        else:
            self.rb_concise.setChecked(True)
            
        self.mode_group.addButton(self.rb_concise)
        self.mode_group.addButton(self.rb_detailed)
        mode_layout.addWidget(self.rb_concise)
        mode_layout.addWidget(self.rb_detailed)
        gen_layout.addLayout(mode_layout)

        self.gen_btn = QPushButton(self.config.T("gen_btn"))
        self.gen_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; height: 30px;")
        self.gen_btn.clicked.connect(self.generate_prompt)
        self.prompt_display = QTextEdit()
        self.prompt_display.setMaximumHeight(80)
        self.prompt_display.setFontFamily("Consolas")
        
        gen_layout.addWidget(self.gen_btn)
        gen_layout.addWidget(self.prompt_display)
        gen_group.setLayout(gen_layout)
        control_layout.addWidget(gen_group)

        # 4. Parse Response
        parse_group = QGroupBox(self.config.T("step4"))
        parse_layout = QVBoxLayout()
        
        self.response_input = QTextEdit()
        self.response_input.setPlaceholderText(self.config.T("placeholder_resp"))
        self.response_input.setMaximumHeight(80)
        
        parse_btn = QPushButton(self.config.T("process_btn"))
        parse_btn.setStyleSheet("background-color: #2196F3; color: white; font-weight: bold; height: 30px;")
        parse_btn.clicked.connect(self.process_response)
        
        clean_btn = QPushButton(self.config.T("clean"))
        clean_btn.setFixedWidth(60)
        clean_btn.clicked.connect(lambda: self.response_input.clear())
        
        self.source_group = QButtonGroup(self)
        platform_layout = QHBoxLayout()
        self.rb_scanpy = QRadioButton(self.config.T("scanpy"))
        self.rb_scanpy.setChecked(True)
        self.rb_seurat = QRadioButton(self.config.T("seurat"))
        self.source_group.addButton(self.rb_scanpy)
        self.source_group.addButton(self.rb_seurat)
        platform_layout.addWidget(QLabel("Target Platform:"))
        platform_layout.addWidget(self.rb_scanpy)
        platform_layout.addWidget(self.rb_scanpy) # Wait, adding rb_scanpy twice? Ah, no the original code added rb_seurat second.
        # Fix:
        platform_layout.replaceWidget(platform_layout.itemAt(2).widget(), self.rb_seurat)
        
        btn_row = QHBoxLayout()
        btn_row.addWidget(parse_btn)
        btn_row.addWidget(clean_btn)
        
        parse_layout.addWidget(self.response_input)
        parse_layout.addLayout(platform_layout)
        parse_layout.addLayout(btn_row)
        parse_group.setLayout(parse_layout)
        control_layout.addWidget(parse_group)

        # 5. Export Code
        export_group = QGroupBox(self.config.T("step5"))
        export_layout = QVBoxLayout()

        self.code_display = QTextEdit()
        self.code_display.setReadOnly(True)
        self.code_display.setFontFamily("Consolas")
        export_layout.addWidget(QLabel(self.config.T("code_label")))
        export_layout.addWidget(self.code_display)
        
        copy_code_btn = QPushButton(self.config.T("copy_code"))
        copy_code_btn.clicked.connect(self.copy_code)
        export_layout.addWidget(copy_code_btn)
        
        export_group.setLayout(export_layout)
        control_layout.addWidget(export_group)

        # --- Right Panel ---
        right_panel_container = QWidget()
        right_panel_layout = QVBoxLayout(right_panel_container)
        
        toolbar_top = QHBoxLayout()
        toolbar_top.addStretch()

        self.pin_btn = QToolButton()
        self.pin_btn.setText("📌")
        self.pin_btn.setCheckable(True)
        self.pin_btn.setToolTip("Always on Top")
        self.pin_btn.clicked.connect(self.toggle_topmost_btn)
        toolbar_top.addWidget(self.pin_btn)

        self.gear_btn = QToolButton()
        self.gear_btn.setText("⚙️")
        self.gear_btn.setToolTip(self.config.T("settings"))
        self.gear_menu = QMenu(self)
        
        self.gear_menu.setStyleSheet("""
            QMenu { 
                menu-scrollable: 1; 
                background-color: white; 
                border: 1px solid #ccc;
                border-radius: 6px;
                padding: 4px 0;
            }
            QMenu::item { 
                padding: 6px 25px 6px 25px; 
                text-align: left; 
                border-radius: 4px;
                margin: 0 4px;
            }
            QMenu::item:selected { 
                background-color: #e0e6f8; 
                color: black; 
            }
            QMenu::separator {
                height: 1px;
                background: #e0e0e0;
                margin: 4px 10px;
            }
        """)

        startup_widget = QWidget()
        startup_layout = QHBoxLayout(startup_widget)
        startup_layout.setContentsMargins(15, 2, 15, 2)
        startup_layout.setSpacing(10)
        
        label = QLabel(self.config.T("default_page"))
        label.setStyleSheet("border: none; background: transparent;")
        startup_layout.addWidget(label)
        
        self.default_page_combo = QComboBox()
        self.default_page_combo.setContextMenuPolicy(Qt.CustomContextMenu)
        self.default_page_combo.customContextMenuRequested.connect(self.on_startup_combo_context_menu)
        self.default_page_combo.currentTextChanged.connect(self.on_default_service_changed)
        startup_layout.addWidget(self.default_page_combo)
        
        startup_action = QWidgetAction(self)
        startup_action.setDefaultWidget(startup_widget)
        self.gear_menu.addAction(startup_action)

        self.gear_menu.addSeparator()

        lang_menu = QMenu(self.config.T("language"), self)
        zh_action = lang_menu.addAction(self.config.T("lang_zh"))
        zh_action.triggered.connect(lambda: self.set_language("zh"))
        en_action = lang_menu.addAction(self.config.T("lang_en"))
        en_action.triggered.connect(lambda: self.set_language("en"))
        
        self.gear_menu.addMenu(lang_menu)

        self.gear_menu.addSeparator()

        guide_action = self.gear_menu.addAction(self.config.T("guide"))
        guide_action.triggered.connect(self.show_help)
        
        self.gear_btn.setMenu(self.gear_menu)
        self.gear_btn.setPopupMode(QToolButton.InstantPopup)
        toolbar_top.addWidget(self.gear_btn)
        
        right_panel_layout.addLayout(toolbar_top)

        self.right_tabs = QTabWidget()
        
        # Tab 1: Enrichment Results
        enrichment_container = QWidget()
        enrichment_layout = QVBoxLayout(enrichment_container)

        self.enrich_status_container = QWidget()
        status_layout = QVBoxLayout(self.enrich_status_container)
        self.status_label = QLabel(self.config.T("enrich_running"))
        self.status_label.setStyleSheet("color: #673AB7; font-weight: bold;")
        self.enrich_progress = QProgressBar()
        self.enrich_progress.setRange(0, 100)
        status_layout.addWidget(self.status_label)
        status_layout.addWidget(self.enrich_progress)
        self.enrich_status_container.hide()
        enrichment_layout.addWidget(self.enrich_status_container)
        
        enrich_scroll = QScrollArea()
        enrich_scroll.setWidgetResizable(True)
        enrich_scroll_content = QWidget()
        self.enrich_detail_layout = QVBoxLayout(enrich_scroll_content)
        
        self.enrich_placeholder = QLabel("No enrichment results yet. Generate a prompt with enrichment enabled to see data here.")
        self.enrich_placeholder.setWordWrap(True)
        self.enrich_placeholder.setAlignment(Qt.AlignCenter)
        self.enrich_detail_layout.addWidget(self.enrich_placeholder)
        
        enrich_scroll.setWidget(enrich_scroll_content)
        enrichment_layout.addWidget(enrich_scroll)
        
        self.right_tabs.addTab(enrichment_container, self.config.T("tab_enrichment"))

        # Tab 2: AI Browser
        browser_container = QWidget()
        browser_layout = QVBoxLayout(browser_container)

        initial_url = self.config.chat_services.get(self.config.default_service, "https://chat.deepseek.com")
        
        browser_header = QHBoxLayout()
        self.url_bar = QLineEdit(initial_url)
        self.url_bar.returnPressed.connect(self.navigate_to_url)
        go_btn = QPushButton(self.config.T("go"))
        go_btn.clicked.connect(self.navigate_to_url)
        
        refresh_btn = QPushButton(self.config.T("refresh"))
        refresh_btn.clicked.connect(lambda: self.browser.reload())

        browser_header.addWidget(QLabel(self.config.T("browser_label")))
        browser_header.addWidget(self.url_bar)
        browser_header.addWidget(go_btn)
        browser_header.addWidget(refresh_btn)
        
        bookmarks_bar_layout = QHBoxLayout()

        self.bookmarks_list = QListWidget()
        self.bookmarks_list.setFlow(QListWidget.LeftToRight)
        self.bookmarks_list.setFixedHeight(35)
        self.bookmarks_list.setDragDropMode(QAbstractItemView.InternalMove)
        self.bookmarks_list.setWrapping(False)
        self.bookmarks_list.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.bookmarks_list.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.bookmarks_list.setContextMenuPolicy(Qt.CustomContextMenu)
        self.bookmarks_list.customContextMenuRequested.connect(self.on_bookmarks_context_menu)
        self.bookmarks_list.itemClicked.connect(self.on_bookmark_clicked)
        self.bookmarks_list.model().rowsMoved.connect(self.on_bookmarks_reordered)
        
        def clear_selection(event):
            item = self.bookmarks_list.itemAt(event.position().toPoint())
            if not item:
                self.bookmarks_list.clearSelection()
            QListWidget.mousePressEvent(self.bookmarks_list, event)
        self.bookmarks_list.mousePressEvent = clear_selection

        self.bookmarks_list.setStyleSheet("""
            QListWidget { background: transparent; border: none; }
            QListWidget::item { 
                background: #f0f0f0; 
                border: 1px solid #ccc; 
                border-radius: 4px; 
                margin: 2px; 
                padding: 4px 8px; 
            }
            QListWidget::item:hover { background: #e0e0e0; }
            QListWidget::item:selected { background: #d0d0d0; border-color: #999; }
        """)
        
        self.refresh_service_buttons()

        add_web_btn = QPushButton(self.config.T("add_web"))
        add_web_btn.setFixedWidth(60)
        add_web_btn.clicked.connect(self.show_add_service_dialog)
        
        bookmarks_bar_layout.addWidget(self.bookmarks_list)
        bookmarks_bar_layout.addSpacing(5) 
        bookmarks_bar_layout.addWidget(add_web_btn)
        bookmarks_bar_layout.addStretch() 
        
        self.browser = QWebEngineView()
        page = QWebEnginePage(self.web_profile, self.browser)
        self.browser.setPage(page)
        
        self.browser.settings().setAttribute(QWebEngineSettings.WebAttribute.JavascriptCanAccessClipboard, True)
        self.browser.settings().setAttribute(QWebEngineSettings.WebAttribute.JavascriptCanOpenWindows, True)
        self.browser.setUrl(QUrl(initial_url))
        
        page.featurePermissionRequested.connect(self.on_permission_requested)
        
        browser_layout.addLayout(browser_header)
        browser_layout.addLayout(bookmarks_bar_layout)
        browser_layout.addWidget(self.browser)
        
        self.right_tabs.addTab(browser_container, self.config.T("tab_browser"))

        right_panel_layout.addWidget(self.right_tabs)

        self.splitter.addWidget(left_container)
        self.splitter.addWidget(right_panel_container)
        self.splitter.setCollapsible(0, False)
        self.splitter.setCollapsible(1, False)

        self.sync_and_save() 

    def on_permission_requested(self, url, feature):
        if feature in (QWebEnginePage.PermissionType.ClipboardReadWrite, 
                       QWebEnginePage.PermissionType.ClipboardWrite):
            self.sender().setFeaturePermission(url, feature, QWebEnginePage.PermissionPolicy.PermissionGrantedByUser)
        else:
            self.sender().setFeaturePermission(url, feature, QWebEnginePage.PermissionPolicy.PermissionDeniedByUser)

    def update_enrich_db_list(self):
        self.enrich_db_combo.clear()
        if self.rb_online.isChecked():
            self.enrich_db_combo.addItems(self.online_db_list)
            self.load_db_btn.hide()
        else:
            if not self.db_path_map:
                self.enrich_db_combo.addItem(self.config.T("local_db_hint"))
            else:
                self.enrich_db_combo.addItems(list(self.db_path_map.keys()))
            self.load_db_btn.show()

    def save_config_state(self):
        if hasattr(self, "species_edit"):
            self.config.session_params["species"] = self.species_edit.text()
            self.config.session_params["tissue"] = self.tissue_edit.text()
            self.config.session_params["top_n"] = self.top_n_edit.text()
            self.config.session_params["exclude"] = self.exclude_edit.text()
            self.config.session_params["prompt_mode"] = "concise" if self.rb_concise.isChecked() else "detailed"
        self.config.save_config()

    def set_language(self, lang):
        if self.config.language == lang:
            return
        self.config.language = lang
        self.config.save_config()
        QMessageBox.information(self, self.config.T("restart_required"), self.config.T("language_changed").format(lang))

    def on_default_service_changed(self, text):
        if text == self.config.T("web_plus"):
            self.show_add_service_dialog()
            self.default_page_combo.blockSignals(True)
            self.default_page_combo.setCurrentText(self.config.default_service)
            self.default_page_combo.blockSignals(False)
            return

        self.config.default_service = text
        self.config.save_config()
    
    def toggle_topmost_btn(self, checked):
        if checked:
            self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)
            self.pin_btn.setStyleSheet("background-color: #87CEEB; border-radius: 4px;") 
        else:
            self.setWindowFlags(self.windowFlags() & ~Qt.WindowStaysOnTopHint)
            self.pin_btn.setStyleSheet("background-color: transparent;") 
        self.show()

    def refresh_service_buttons(self):
        self.bookmarks_list.clear()
        total_width = 0
        fm = self.bookmarks_list.fontMetrics()
        for name, url in self.config.chat_services.items():
            item = QListWidgetItem(name)
            item.setData(Qt.UserRole, url)
            item.setTextAlignment(Qt.AlignCenter)
            self.bookmarks_list.addItem(item)
            item_w = fm.horizontalAdvance(name) + 25 
            total_width += item_w
        self.bookmarks_list.setMaximumWidth(total_width + 10)

    def on_bookmark_clicked(self, item):
        url = item.data(Qt.UserRole)
        self.browser.setUrl(QUrl(url))
        self.url_bar.setText(url)

    def browse_db_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Database Folder")
        if folder:
            self.load_databases_from_folder(folder)

    def load_databases_from_folder(self, folder):
        gmt_files = [f for f in os.listdir(folder) if f.endswith(('.gmt', '.txt'))]
        if not gmt_files:
            QMessageBox.warning(self, "No Databases Found", "No .gmt or .txt files found in the selected folder.")
            return
            
        self.db_dir = folder
        self.db_path_map = {f: f for f in gmt_files}
        self.custom_db_label.setText(f"Folder: ...{folder[-30:]}")
        self.rb_local.setChecked(True)
        self.update_enrich_db_list()
        QMessageBox.information(self, "Success", f"Loaded {len(gmt_files)} databases from folder.")

    def on_bookmarks_reordered(self, parent, start, end, destination, row):
        new_services = {}
        for i in range(self.bookmarks_list.count()):
            item = self.bookmarks_list.item(i)
            name = item.text()
            url = item.data(Qt.UserRole)
            new_services[name] = url
        self.config.chat_services = new_services
        self.config.save_config()

    def on_bookmarks_context_menu(self, pos):
        item = self.bookmarks_list.itemAt(pos)
        if not item: return
        name = item.text()
        url = item.data(Qt.UserRole)

        menu = QMenu()
        edit_action = menu.addAction(self.config.T("edit_fmt").format(name))
        delete_action = menu.addAction(self.config.T("delete_fmt").format(name))
        
        action = menu.exec(self.bookmarks_list.mapToGlobal(pos))
        
        if action == edit_action:
            self.show_add_service_dialog(edit_mode=True, old_name=name, old_url=url)
        elif action == delete_action:
            if len(self.config.chat_services) <= 1:
                QMessageBox.warning(self, self.config.T("warning"), self.config.T("cannot_delete_last"))
                return
            confirm = QMessageBox.question(self, self.config.T("confirm_delete"), self.config.T("delete_confirmation").format(name), QMessageBox.Yes | QMessageBox.No)
            if confirm == QMessageBox.Yes:
                del self.config.chat_services[name]
                if name in self.all_possible_services:
                    del self.all_possible_services[name]
                if self.config.default_service == name:
                    self.config.default_service = list(self.config.chat_services.keys())[0]
                self.sync_and_save()

    def sync_and_save(self):
        self.default_page_combo.blockSignals(True)
        self.default_page_combo.clear()
        services = list(self.config.chat_services.keys())
        self.default_page_combo.addItems(services)
        self.default_page_combo.addItem(self.config.T("web_plus"))

        if self.config.default_service in services:
             self.default_page_combo.setCurrentText(self.config.default_service)
        else:
             if services:
                self.config.default_service = services[0]
                self.default_page_combo.setCurrentText(self.config.default_service)
             
        self.default_page_combo.blockSignals(False)
        self.refresh_service_buttons()
        self.config.save_config()

    def on_startup_combo_context_menu(self, pos):
        name = self.default_page_combo.currentText()
        if not name or name not in self.config.chat_services: return
        url = self.config.chat_services[name]

        menu = QMenu()
        edit_action = menu.addAction(self.config.T("edit_fmt").format(name))
        delete_action = menu.addAction(self.config.T("delete_fmt").format(name))
        
        action = menu.exec(self.default_page_combo.mapToGlobal(pos))
        
        if action == edit_action:
            self.show_add_service_dialog(edit_mode=True, old_name=name, old_url=url)
        elif action == delete_action:
            if len(self.config.chat_services) <= 1:
                QMessageBox.warning(self, self.config.T("warning"), self.config.T("cannot_delete_last"))
                return
            confirm = QMessageBox.question(self, self.config.T("confirm_delete"), self.config.T("delete_confirmation").format(name), QMessageBox.Yes | QMessageBox.No)
            if confirm == QMessageBox.Yes:
                del self.config.chat_services[name]
                if self.config.default_service == name:
                    self.config.default_service = list(self.config.chat_services.keys())[0]
                self.sync_and_save()

    def show_add_service_dialog(self, edit_mode=False, old_name="", old_url=""):
        dialog = QDialog(self)
        dialog.setWindowTitle(self.config.T("edit_web_link") if edit_mode else self.config.T("add_web_link"))
        layout = QFormLayout(dialog)
        
        name_edit = QLineEdit(old_name)
        url_edit = QLineEdit(old_url if edit_mode else "https://")
        
        layout.addRow(self.config.T("web_name"), name_edit)
        layout.addRow(self.config.T("url"), url_edit)
        
        btn_box = QHBoxLayout()
        add_btn = QPushButton(self.config.T("save") if edit_mode else self.config.T("add"))
        cancel_btn = QPushButton(self.config.T("cancel"))
        btn_box.addWidget(add_btn)
        btn_box.addWidget(cancel_btn)
        layout.addRow(btn_box)
        
        def do_add():
            name = name_edit.text().strip()
            url = url_edit.text().strip()
            if name and url.startswith("http"):
                if edit_mode and old_name != name:
                    if old_name in self.config.chat_services: del self.config.chat_services[old_name]
                    if old_name in self.all_possible_services: del self.all_possible_services[old_name]
                
                self.config.chat_services[name] = url
                self.all_possible_services[name] = url
                
                self.config.default_service = name
                self.sync_and_save()
                dialog.accept()
            else:
                QMessageBox.warning(dialog, self.config.T("input_error"), self.config.T("input_error_msg"))
        
        add_btn.clicked.connect(do_add)
        cancel_btn.clicked.connect(dialog.reject)
        dialog.exec()

    def start_enrichment_thread(self):
        path = self.path_edit.text().strip()
        if not os.path.exists(path):
            QMessageBox.critical(self, self.config.T("error"), self.config.T("file_not_found"))
            return
            
        is_local = self.rb_local.isChecked()
        db_name = self.enrich_db_combo.currentText()
        
        if is_local:
            db_file = self.db_path_map.get(db_name, "")
            enrich_db_path = os.path.join(self.db_dir, db_file)
            if not db_file or not os.path.exists(enrich_db_path):
                QMessageBox.critical(self, self.config.T("error"), f"Local database file not found: {db_file}\nSearched in: {self.db_dir}")
                return
        else:
            enrich_db_path = db_name
            
        self.run_enrich_btn.setEnabled(False)
        self.enrich_status_container.show()
        self.status_label.setText(self.config.T("enrich_running"))
        self.enrich_progress.setValue(0)
        
        self.right_tabs.setCurrentIndex(0)
        self.clear_enrichment_display()
        
        self.worker = EnrichmentWorker(path, self.species_edit.text(), enrich_db_path, is_local=is_local)
        self.worker.progress.connect(self.enrich_progress.setValue)
        self.worker.finished.connect(self.on_enrichment_finished)
        self.worker.error.connect(self.on_enrichment_error)
        self.worker.start()

    def on_enrichment_finished(self, data):
        self.enrichment_data = data
        self.run_enrich_btn.setEnabled(True)
        self.enrich_status_container.hide() 
        self.update_enrichment_display(data)
        QMessageBox.information(self, self.config.T("success"), "Enrichment complete! AI hints are ready for Step 3.")

    def on_enrichment_error(self, err_msg):
        self.run_enrich_btn.setEnabled(True)
        self.enrich_status_container.hide() 
        QMessageBox.critical(self, self.config.T("error"), f"Enrichment failed: {err_msg}")

    def clear_enrichment_display(self):
        for i in reversed(range(self.enrich_detail_layout.count())): 
            item = self.enrich_detail_layout.itemAt(i)
            if item.widget():
                item.widget().setParent(None)
        
        placeholder = QLabel(self.config.T("enrich_running")) 
        placeholder.setAlignment(Qt.AlignCenter)
        self.enrich_detail_layout.addWidget(placeholder)

    def update_enrichment_display(self, enrichment_data):
        for i in reversed(range(self.enrich_detail_layout.count())): 
            item = self.enrich_detail_layout.itemAt(i)
            if item.widget():
                item.widget().setParent(None)
        
        for cluster, data in enrichment_data.items():
            c_group = QGroupBox(f"Cluster {cluster}")
            c_layout = QVBoxLayout()
            
            c_layout.addWidget(QLabel("<b>[AI Prompt Hints]</b>"))
            hints_text = "; ".join(data['hints'])
            c_layout.addWidget(QLabel(hints_text))
            
            if data.get('table_path'):
                c_layout.addWidget(QLabel("<b>[Database Results]</b>"))
                btn_open = QPushButton(f"Open Full Result Table (CSV)")
                def make_open_func(p):
                    return lambda: os.startfile(p)
                btn_open.clicked.connect(make_open_func(data['table_path']))
                c_layout.addWidget(btn_open)
            
            if data.get('plot_path') and os.path.exists(data['plot_path']):
                c_layout.addWidget(QLabel("<b>[Enrichment Visualization]</b>"))
                img_label = QLabel()
                pixmap = QPixmap(data['plot_path'])
                img_label.setPixmap(pixmap.scaledToWidth(500, Qt.SmoothTransformation))
                c_layout.addWidget(img_label)

            c_group.setLayout(c_layout)
            self.enrich_detail_layout.addWidget(c_group)
        
        self.enrich_detail_layout.addStretch()

    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Marker File", "", "CSV/TSV Files (*.csv *.tsv *.txt)")
        if file_path:
            self.path_edit.setText(file_path)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        files = [u.toLocalFile() for u in event.mimeData().urls()]
        if not files:
            return
            
        path = files[0]
        if os.path.isdir(path):
            self.load_databases_from_folder(path)
        else:
            self.path_edit.setText(path)

    def navigate_to_url(self):
        url = self.url_bar.text()
        if not url.startswith("http"):
            url = "https://" + url
        self.browser.setUrl(QUrl(url))

    def generate_prompt(self):
        path = self.path_edit.text().strip()
        if not os.path.exists(path):
            QMessageBox.critical(self, self.config.T("error"), self.config.T("file_not_found"))
            return
        
        try:
            mode = "concise" if self.rb_concise.isChecked() else "detailed"
            
            enrichment_hints = None
            if self.enrichment_data:
                enrichment_hints = {k: v['hints'] for k, v in self.enrichment_data.items()}

            # Import core relative logic
            from chatcellanno.core import annotate_cell_types
            
            prompt, _ = annotate_cell_types(
                marker_file=path,
                step="generate",
                species=self.species_edit.text(),
                tissue=self.tissue_edit.text(),
                top_n=int(self.top_n_edit.text()),
                mode=mode,
                exclude_types=self.exclude_edit.text(),
                use_enrichment=False, 
                enrichment_hints=enrichment_hints
            )
            
            self.prompt_display.setText(prompt)
            pyperclip.copy(prompt)

            self.right_tabs.setCurrentIndex(1)

            QMessageBox.information(self, self.config.T("success"), self.config.T("prompt_copied"))
        except Exception as e:
            QMessageBox.critical(self, self.config.T("error"), f"Failed: {str(e)}")

    def process_response(self):
        path = self.path_edit.text().strip()
        response = self.response_input.toPlainText().strip()
        
        if not response:
            QMessageBox.warning(self, self.config.T("warning"), self.config.T("paste_warning"))
            return

        if not path or not os.path.exists(path):
            reply = QMessageBox.question(
                self, 
                self.config.T("no_file_title"), 
                self.config.T("no_file_msg"),
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
            path = None

        try:
            source = "scanpy" if self.rb_scanpy.isChecked() else "seurat"
            
            from chatcellanno.core import annotate_cell_types

            result_df, code_snippet = annotate_cell_types(
                marker_file=path,
                step="parse",
                response_text=response,
                source=source,
                top_n=int(self.top_n_edit.text())
            )
            
            self.last_annotations = result_df
            self.code_display.setText(code_snippet)
            
            QMessageBox.information(self, self.config.T("success"), "AI response processed! Code generated in Step 5.")
        except Exception as e:
            QMessageBox.critical(self, self.config.T("error"), f"Parsing failed: {str(e)}")

    def copy_code(self):
        pyperclip.copy(self.code_display.toPlainText())
        QMessageBox.information(self, "Copied", "Code copied to clipboard.")

    def show_help(self):
        msg = """
        ChatCellAnno User Guide:
        
        Step 1 (Data Source): 
        - Load your marker file (CSV/TSV/TXT) from Scanpy/Seurat.
        
        Step 2 (Enrichment): 
        - [Optional] Run functional enrichment analysis (Online via Enrichr or Local GMT files).
        - Drag a folder into the software to load offline databases.
        - Database results help AI minimize hallucinations.
        
        Step 3 (Prompt Configuration): 
        - Configure species, tissue, and mode (Concise/Detailed).
        - Click 'Generate & Copy Prompt' - it's automatically copied to your clipboard.
        
        Step 4 (AI Interaction):
        - Use the AI Browser on the right to visit ChatGPT, DeepSeek, or Claude.
        - Paste (Ctrl+V) and Send.
        - Copy the AI's MARKDOWN TABLE response.
        - Select your target platform (Scanpy/Seurat) and click 'Process AI Output'.
        
        Step 5 (Export):
        - Copy and run the generated Python/R code to apply annotations to your data.
        """
        QMessageBox.information(self, self.config.T("guide"), msg.strip())
