"""
ChatCellAnno - PySide6 GUI Entry Point

This refactored GUI uses PySide6 for a modern interface and includes 
an embedded web browser (QtWebEngine) to access LLMs directly within the app.

Architecture:
- MainWindow: Main container using QMainWindow.
- ControlPanel: Left side with all parameters and buttons.
- BrowserPanel: Right side with an embedded QWebEngineView.
"""

import sys
import os
import json
import pyperclip
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QLabel, QLineEdit, QPushButton, QRadioButton, QButtonGroup, 
    QGroupBox, QTextEdit, QFileDialog, QMessageBox, QTableWidget, 
    QTableWidgetItem, QHeaderView, QCheckBox, QTabWidget, QScrollArea,
    QSplitter, QComboBox, QDialog, QFormLayout, QToolButton, QStyle,
    QListWidget, QListWidgetItem, QAbstractItemView, QMenu, QWidgetAction
)
from PySide6.QtGui import QIcon, QAction
from PySide6.QtCore import Qt, QUrl
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebEngineCore import QWebEngineProfile, QWebEnginePage, QWebEngineSettings

# Core logic imported lazily/later or protected
try:
    # Try importing initially to check presence, but we might rely on 
    # lazy import inside method if we want speed, though Python caches it.
    # For now, keep it global as moving it requires refactoring all usages.
    from chatcellanno.core import annotate_cell_types
except ImportError:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from chatcellanno.core import annotate_cell_types

# Late import pandas only when needed if possible, 
# but it's used in extract_markers so we'll leave it 
# or move it to a function if it's the bottleneck.
# gui.py seems to use pandas in 'process_ai_response' indirectly via core?
# Let's check where 'pandas' (pd) is used in gui.py.
# Based on previous read, 'import pandas as pd' is at line 18.
import pandas as pd

TRANSLATIONS = {
    "en": {
        "title": "ChatCellAnno",
        "step1": "Step 1: Data & Source",
        "placeholder_path": "Drag & drop file here or browse...",
        "browse": "Browse",
        "scanpy": "Scanpy (Python)",
        "seurat": "Seurat (R)",
        "step2": "Step 2: Configuration",
        "species": "Species:",
        "tissue": "Tissue:",
        "top_n": "Top Markers:",
        "exclude": "Exclude Types (optional):",
        "prompt_mode": "Prompt Mode:",
        "concise": "Concise",
        "detailed": "Detailed",
        "step3": "Step 3: Generate Prompt",
        "gen_btn": "Generate & Copy Prompt",
        "step4": "Step 4: Result & Code",
        "placeholder_resp": "Paste AI response here...",
        "process_btn": "Process AI Output",
        "code_label": "Generated Code:",
        "copy_code": "Copy Code",
        "default_page": "Default Page:",
        "browser_label": "Browser:",
        "go": "Go",
        "refresh": "Refresh",
        "guide": "User Guide",
        "settings": "Settings",
        "add_web": "+ Web",
        "language": "Language",
        "confirm_delete": "Confirm Delete",
        "delete_msg": "Are you sure you want to delete '{}'?",
        "last_service_msg": "Cannot delete the last service.",
        "success": "Success",
        "prompt_copied": "Prompt copied! Paste it into the browser on the right.",
        "error": "Error",
        "file_not_found": "File not found!",
        "warning": "Warning",
        "paste_warning": "Please paste the AI response first.",
        "no_file_title": "No Marker File",
        "no_file_msg": "No marker file selected.\n\nDo you want to generate code using ONLY the AI response? (Ensure the AI output contains Cluster IDs)",
        "edit_web": "Edit Web Link",
        "add_web_title": "Add Web Link",
        "web_name": "Web Name:",
        "web_url": "URL:",
        "save": "Save",
        "add": "Add",
        "cancel": "Cancel",
        "input_error": "Input Error",
        "input_error_msg": "Please provide a valid name and URL (starting with http/https).",
        "restart_msg": "Please restart the application to apply language changes.",
        "edit_fmt": "Edit '{}'",
        "delete_fmt": "Delete '{}'",
        "cannot_delete_last": "Cannot delete the last service.",
        "delete_confirmation": "Are you sure you want to delete '{}'?",
        "edit_web_link": "Edit Web Link",
        "add_web_link": "Add Web Link",
        "url": "URL:",
        "restart_required": "Restart Required",
        "language_changed": "Language changed to {}. Please restart the application.",
        "web_plus": "+ Web",
        "toggle_language": "Toggle Language",
        "open_in_new": "Open in New Window",
        "initial_page": "Default Page",
        "chat_services": "Chat Services",
        "bookmarks": "Bookmarks",
        "lang_zh": "Chinese",
        "lang_en": "English"
    },
    "zh": {
        "title": "ChatCellAnno - 单细胞注释助手",
        "step1": "步骤 1: 数据源",
        "placeholder_path": "拖放文件或点击浏览...",
        "browse": "浏览",
        "scanpy": "Scanpy (Python)",
        "seurat": "Seurat (R)",
        "step2": "步骤 2: 参数配置",
        "species": "物种:",
        "tissue": "组织:",
        "top_n": "Top Markers:",
        "exclude": "排除类型 (可选):",
        "prompt_mode": "提示词模式:",
        "concise": "简洁",
        "detailed": "详细",
        "step3": "步骤 3: 生成提示词",
        "gen_btn": "生成并复制提示词",
        "step4": "步骤 4: 结果转换",
        "placeholder_resp": "在此粘贴 AI 回复...",
        "process_btn": "解析 AI 回复",
        "code_label": "生成的代码:",
        "copy_code": "复制代码",
        "default_page": "默认主页:",
        "browser_label": "浏览器:",
        "go": "前往",
        "refresh": "刷新",
        "guide": "用户指南",
        "settings": "设置",
        "add_web": "+ Web",
        "language": "语言 / Language",
        "confirm_delete": "确认删除",
        "delete_msg": "确定要删除 '{}' 吗?",
        "last_service_msg": "无法删除最后一个服务。",
        "success": "成功",
        "prompt_copied": "提示词已复制！请粘贴到右侧浏览器中。",
        "error": "错误",
        "file_not_found": "未找到文件！",
        "warning": "警告",
        "paste_warning": "请先粘贴 AI 回复内容。",
        "no_file_title": "未选择 Marker 文件",
        "no_file_msg": "未检测到 Marker 文件。\n\n是否仅根据 AI 回复生成代码？(请确保 AI 回复中包含 Cluster ID)",
        "edit_web": "编辑网页链接",
        "add_web_title": "添加网页链接",
        "web_name": "名称:",
        "web_url": "链接 (URL):",
        "save": "保存",
        "add": "添加",
        "cancel": "取消",
        "input_error": "输入错误",
        "input_error_msg": "请输入有效的名称和网址 (需以 http/https 开头)。",
        "restart_msg": "语言设置已更改，请重启软件生效。",
        "edit_fmt": "编辑 '{}'",
        "delete_fmt": "删除 '{}'",
        "cannot_delete_last": "无法删除最后一个服务。",
        "delete_confirmation": "您确定要删除 '{}' 吗？",
        "edit_web_link": "编辑链接",
        "add_web_link": "添加链接",
        "url": "URL：",
        "restart_required": "需要重启",
        "language_changed": "语言已切换为 {}。请重启应用以生效。",
        "web_plus": "+ 网站",
        "toggle_language": "切换语言",
        "open_in_new": "在新窗口打开",
        "initial_page": "默认主页",
        "chat_services": "聊天服务",
        "bookmarks": "书签",
        "lang_zh": "中文",
        "lang_en": "English"
    }
}

class ChatCellAnnoApp(QMainWindow):
    def __init__(self):
        super().__init__()
        
        # Load Config First to get Language
        self.config_paths = {
            "config": "config.json",
        }
        self.all_possible_services = {} 
        self.chat_services = {
            "DeepSeek": "https://chat.deepseek.com",
            "ChatGPT": "https://chatgpt.com", 
            "Claude": "https://claude.ai",
            "Copilot": "https://copilot.microsoft.com",
            "Doubao": "https://www.doubao.com/chat/"
        }
        self.default_service = "DeepSeek"
        self.language = "en"
        self.session_params = {"species": "Human", "tissue": "PBMC", "top_n": "10", "exclude": ""}
        
        self.load_config() # Updates self.language
        
        self.setWindowTitle(self.tr("title"))
        
        # Apply Global Styling for GroupBox Titles
        self.setStyleSheet("""
            QGroupBox::title {
                font-size: 15px;
                font-weight: bold;
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 0 3px;
            }
        """)

        # 1. Default Size (Moderate)
        self.resize(1280, 800)
        
        # 2. Configurable Chat Services logic
        # Use user's home directory or current directory for config to ensure write permissions
        # and persistence across EXE updates
        user_config_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.config_paths = {
            "config": os.path.join(user_config_dir, "config.json"),
            "web_storage": os.path.join(user_config_dir, "chatcellanno_web_data")
        }
        
        # Configure WebEngine Profile for Persistence
        self.web_profile = QWebEngineProfile("ChatCellAnnoProfile", self)
        self.web_profile.setPersistentStoragePath(self.config_paths["web_storage"])
        self.web_profile.setCachePath(self.config_paths["web_storage"])
        
        # Default fallback services (Clean: Only DeepSeek)
        self.all_possible_services = {
            "DeepSeek": "https://chat.deepseek.com"
        }
        self.chat_services = {"DeepSeek": "https://chat.deepseek.com"}
        self.default_service = "DeepSeek"
        
        # Default session params
        self.session_params = {
            "species": "Human",
            "tissue": "PBMC",
            "top_n": "10",
            "exclude": ""
        }
        
        self.load_config()
        
        # State variables
        self.current_marker_file = ""
        
        self.setWindowTitle(self.T("title"))
        if os.path.exists("app_icon.ico"):
            self.setWindowIcon(QIcon("app_icon.ico"))

        # Global stylesheet for a more modern look
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
        self.setAcceptDrops(True) # Enable native drag and drop
        
        # 3. Default 50/50 Split
        # Set stretch factors to equal values (1:1)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 1)
        # Note: splitter sizes are set in setup_ui but will be managed by layout stretch afterwards

    def setup_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)


        # Use a splitter for resizable panels
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
        data_group = QGroupBox(self.T("step1"))
        data_layout = QVBoxLayout()
        
        file_path_layout = QHBoxLayout()
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText(self.T("placeholder_path"))
        browse_btn = QPushButton(self.T("browse"))
        browse_btn.clicked.connect(self.browse_file)
        file_path_layout.addWidget(self.path_edit)
        file_path_layout.addWidget(browse_btn)
        data_layout.addLayout(file_path_layout)

        self.source_group = QButtonGroup(self)
        source_rb_layout = QHBoxLayout()
        self.rb_scanpy = QRadioButton(self.T("scanpy"))
        self.rb_scanpy.setChecked(True)
        self.rb_seurat = QRadioButton(self.T("seurat"))
        self.source_group.addButton(self.rb_scanpy)
        self.source_group.addButton(self.rb_seurat)
        source_rb_layout.addWidget(self.rb_scanpy)
        source_rb_layout.addWidget(self.rb_seurat)
        data_layout.addLayout(source_rb_layout)
        data_group.setLayout(data_layout)
        control_layout.addWidget(data_group)

        # 2. Configuration
        conf_group = QGroupBox(self.T("step2"))
        conf_layout = QVBoxLayout()
        
        self.species_edit = QLineEdit(self.session_params.get("species", "Human"))
        self.tissue_edit = QLineEdit(self.session_params.get("tissue", "PBMC"))
        self.top_n_edit = QLineEdit(self.session_params.get("top_n", "10"))
        self.exclude_edit = QLineEdit(self.session_params.get("exclude", ""))
        self.exclude_edit.setPlaceholderText("e.g. Neuron, Astrocyte")
        
        # Connect changes to save
        self.species_edit.textChanged.connect(self.save_config)
        self.tissue_edit.textChanged.connect(self.save_config)
        self.top_n_edit.textChanged.connect(self.save_config)
        self.exclude_edit.textChanged.connect(self.save_config)

        conf_layout.addWidget(QLabel(self.T("species")))
        conf_layout.addWidget(self.species_edit)
        conf_layout.addWidget(QLabel(self.T("tissue")))
        conf_layout.addWidget(self.tissue_edit)
        conf_layout.addWidget(QLabel(self.T("top_n")))
        conf_layout.addWidget(self.top_n_edit)
        conf_layout.addWidget(QLabel(self.T("exclude")))
        conf_layout.addWidget(self.exclude_edit)

        conf_layout.addWidget(QLabel(self.T("prompt_mode")))
        self.mode_group = QButtonGroup(self)
        mode_layout = QHBoxLayout()
        self.rb_concise = QRadioButton(self.T("concise"))
        self.rb_detailed = QRadioButton(self.T("detailed"))
        self.rb_detailed.setChecked(True)
        self.mode_group.addButton(self.rb_concise)
        self.mode_group.addButton(self.rb_detailed)
        mode_layout.addWidget(self.rb_concise)
        mode_layout.addWidget(self.rb_detailed)
        conf_layout.addLayout(mode_layout)
        conf_group.setLayout(conf_layout)
        control_layout.addWidget(conf_group)

        # 3. Prompt Generation
        gen_group = QGroupBox(self.T("step3"))
        gen_layout = QVBoxLayout()
        gen_btn = QPushButton(self.T("gen_btn"))
        gen_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; height: 30px;")
        gen_btn.clicked.connect(self.generate_prompt)
        self.prompt_display = QTextEdit()
        self.prompt_display.setMaximumHeight(100)
        self.prompt_display.setFontFamily("Consolas")
        
        gen_layout.addWidget(gen_btn)
        gen_layout.addWidget(self.prompt_display)
        gen_group.setLayout(gen_layout)
        control_layout.addWidget(gen_group)

        # 4. Parse Response
        parse_group = QGroupBox(self.T("step4"))
        parse_layout = QVBoxLayout()
        
        self.response_input = QTextEdit()
        self.response_input.setPlaceholderText(self.T("placeholder_resp"))
        self.response_input.setMaximumHeight(100)
        
        parse_btn = QPushButton(self.T("process_btn"))
        parse_btn.setStyleSheet("background-color: #2196F3; color: white; font-weight: bold; height: 30px;")
        parse_btn.clicked.connect(self.process_response)
        
        # Code Output area
        self.code_display = QTextEdit()
        self.code_display.setReadOnly(True)
        self.code_display.setFontFamily("Consolas")
        code_container = QWidget()
        code_vbox = QVBoxLayout(code_container)
        code_vbox.addWidget(QLabel(self.T("code_label")))
        code_vbox.addWidget(self.code_display)
        copy_code_btn = QPushButton(self.T("copy_code"))
        copy_code_btn.clicked.connect(self.copy_code)
        code_vbox.addWidget(copy_code_btn)

        parse_layout.addWidget(self.response_input)
        parse_layout.addWidget(parse_btn)
        parse_layout.addWidget(code_container)
        parse_group.setLayout(parse_layout)
        control_layout.addWidget(parse_group)

        # --- Right Panel: Browser ---
        browser_container = QWidget()
        browser_layout = QVBoxLayout(browser_container)
        
        # 0. Browser Toolbar (Top Bar)
        toolbar_top = QHBoxLayout()
        toolbar_top.addStretch()

        self.pin_btn = QToolButton()
        self.pin_btn.setText("📌")
        self.pin_btn.setCheckable(True)
        self.pin_btn.setToolTip("Always on Top")
        self.pin_btn.clicked.connect(self.toggle_topmost_btn)
        toolbar_top.addWidget(self.pin_btn)

        # Gear Menu for Settings
        self.gear_btn = QToolButton()
        self.gear_btn.setText("⚙️")
        self.gear_btn.setToolTip(self.T("settings"))
        self.gear_menu = QMenu(self)
        
        # Adjust menu style for alignment and hover color. 
        # Standard QMenu::item padding is roughly 4px vertical, 20-30px horizontal.
        # We also style the separator.
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

        # -- Startup Setting (using QWidgetAction) --
        # We reuse the padding logic to try and match QMenu::item visually.
        # QMenu padding (left) is usually where the icon/check goes. 
        # For a QWidgetAction, we need to manually offset or accept it's a widget.
        startup_widget = QWidget()
        startup_layout = QHBoxLayout(startup_widget)
        # Margins: Left/Right 20px to match text indent roughly, Top/Bottom 4px 
        startup_layout.setContentsMargins(15, 2, 15, 2)
        startup_layout.setSpacing(10)
        
        label = QLabel(self.T("default_page"))
        label.setStyleSheet("border: none; background: transparent;")
        startup_layout.addWidget(label)
        
        self.default_page_combo = QComboBox()
        # Items populated in sync_and_save()
        self.default_page_combo.setContextMenuPolicy(Qt.CustomContextMenu)
        self.default_page_combo.customContextMenuRequested.connect(self.on_startup_combo_context_menu)
        # Connect currentTextChanged instead of activated to get the text directly
        self.default_page_combo.currentTextChanged.connect(self.on_default_service_changed)
        startup_layout.addWidget(self.default_page_combo)
        
        startup_action = QWidgetAction(self)
        startup_action.setDefaultWidget(startup_widget)
        self.gear_menu.addAction(startup_action)

        self.gear_menu.addSeparator()

        # -- Language Sub-menu --
        lang_menu = QMenu(self.T("language"), self)
        
        zh_action = lang_menu.addAction(self.T("lang_zh"))
        zh_action.triggered.connect(lambda: self.set_language("zh"))
        
        en_action = lang_menu.addAction(self.T("lang_en"))
        en_action.triggered.connect(lambda: self.set_language("en"))
        
        self.gear_menu.addMenu(lang_menu)

        self.gear_menu.addSeparator()

        # -- Guide Action --
        guide_action = self.gear_menu.addAction(self.T("guide"))
        guide_action.triggered.connect(self.show_help)
        
        self.gear_btn.setMenu(self.gear_menu)
        self.gear_btn.setPopupMode(QToolButton.InstantPopup)
        toolbar_top.addWidget(self.gear_btn)
        
        browser_layout.addLayout(toolbar_top)

        initial_url = self.all_possible_services.get(self.default_service, "https://chat.deepseek.com")
        
        browser_header = QHBoxLayout()
        self.url_bar = QLineEdit(initial_url)
        self.url_bar.returnPressed.connect(self.navigate_to_url)
        go_btn = QPushButton(self.T("go"))
        go_btn.clicked.connect(self.navigate_to_url)
        
        # Add Refresh Button next to Go
        refresh_btn = QPushButton(self.T("refresh"))
        refresh_btn.clicked.connect(lambda: self.browser.reload())

        browser_header.addWidget(QLabel(self.T("browser_label")))
        browser_header.addWidget(self.url_bar)
        browser_header.addWidget(go_btn)
        browser_header.addWidget(refresh_btn)
        
        # Service Buttons Toolbar (Bookmarks bar with drag-reorder)
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
        
        # Allow clearing selection when clicking on empty space
        original_mouse_press = self.bookmarks_list.mousePressEvent
        def clear_selection(event):
            item = self.bookmarks_list.itemAt(event.position().toPoint())
            if not item:
                self.bookmarks_list.clearSelection()
            original_mouse_press(event)
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

        add_web_btn = QPushButton(self.T("add_web"))
        add_web_btn.setFixedWidth(60)
        # add_web_btn.clicked.connect(self.show_add_service_dialog) # Use T for logic or keep same?
        add_web_btn.clicked.connect(self.show_add_service_dialog)
        
        bookmarks_bar_layout.addWidget(self.bookmarks_list)
        bookmarks_bar_layout.addSpacing(5) # Small gap
        bookmarks_bar_layout.addWidget(add_web_btn)
        bookmarks_bar_layout.addStretch() # Push both to the left
        
        # Initialize WebEngine with Profile
        self.browser = QWebEngineView()
        page = QWebEnginePage(self.web_profile, self.browser)
        self.browser.setPage(page)
        
        # Allow JavaScript to access clipboard (fixes Copy button in web apps)
        self.browser.settings().setAttribute(QWebEngineSettings.WebAttribute.JavascriptCanAccessClipboard, True)
        self.browser.settings().setAttribute(QWebEngineSettings.WebAttribute.JavascriptCanOpenWindows, True)

        self.browser.setUrl(QUrl(initial_url))
        
        # Handle permission requests (Clipboard)
        page.featurePermissionRequested.connect(self.on_permission_requested)
        
        browser_layout.addLayout(browser_header)
        browser_layout.addLayout(bookmarks_bar_layout)
        browser_layout.addWidget(self.browser)

        # Add both to splitter
        self.splitter.addWidget(left_container)
        self.splitter.addWidget(browser_container)
        # Ensure 50/50 split initially
        self.splitter.setCollapsible(0, False)
        self.splitter.setCollapsible(1, False)

        # Initial combo population
        self.sync_and_save() 


    # --- Slots & Logic ---
    
    def on_permission_requested(self, url, feature):
        # Grant clipboard reading/writing permissions
        if feature in (QWebEnginePage.PermissionType.ClipboardReadWrite, 
                       QWebEnginePage.PermissionType.ClipboardWrite):
            self.sender().setFeaturePermission(url, feature, QWebEnginePage.PermissionPolicy.PermissionGrantedByUser)
        else:
            self.sender().setFeaturePermission(url, feature, QWebEnginePage.PermissionPolicy.PermissionDeniedByUser)

    def T(self, key):
        return TRANSLATIONS.get(self.language, TRANSLATIONS["en"]).get(key, key)

    def load_config(self):
        self.default_service = "DeepSeek"
        if os.path.exists(self.config_paths["config"]):
            try:
                with open(self.config_paths["config"], "r", encoding="utf-8") as f:
                    config = json.load(f)
                    self.default_service = config.get("default_service", "DeepSeek")
                    self.session_params = config.get("session_params", self.session_params)
                    custom_services = config.get("chat_services")
                    self.language = config.get("language", "en") # Load language
                    if custom_services:
                        self.chat_services = custom_services
            except:
                pass
        self.all_possible_services = self.chat_services.copy()

    def save_config(self):
        # Update session params from UI if initialized
        if hasattr(self, "species_edit"):
            self.session_params["species"] = self.species_edit.text()
            self.session_params["tissue"] = self.tissue_edit.text()
            self.session_params["top_n"] = self.top_n_edit.text()
            self.session_params["exclude"] = self.exclude_edit.text()

        config = {
            "default_service": self.default_service,
            "session_params": self.session_params,
            "chat_services": self.chat_services,
            "language": self.language
        }
        try:
            with open(self.config_paths["config"], "w", encoding="utf-8") as f:
                json.dump(config, f, indent=4)
        except:
            pass

    def set_language(self, lang):
        if self.language == lang:
            return
        self.language = lang
        self.save_config()
        QMessageBox.information(self, self.T("restart_required"), self.T("language_changed").format(lang))

    def on_default_service_changed(self, text):
        if text == self.T("web_plus"):
            self.show_add_service_dialog()
            # Revert to valid selection in UI until new service is added
            self.default_page_combo.blockSignals(True)
            self.default_page_combo.setCurrentText(self.default_service)
            self.default_page_combo.blockSignals(False)
            return

        self.default_service = text
        self.save_config()
    
    def toggle_topmost_btn(self, checked):
        if checked:
            self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)
            self.pin_btn.setStyleSheet("background-color: #87CEEB;") # Visual cue
        else:
            self.setWindowFlags(self.windowFlags() & ~Qt.WindowStaysOnTopHint)
            self.pin_btn.setStyleSheet("")
        self.show()

    def refresh_service_buttons(self):
        self.bookmarks_list.clear()
        total_width = 0
        fm = self.bookmarks_list.fontMetrics()
        for name, url in self.chat_services.items():
            item = QListWidgetItem(name)
            item.setData(Qt.UserRole, url)
            item.setTextAlignment(Qt.AlignCenter)
            self.bookmarks_list.addItem(item)
            # Calculate width based on text + padding + margin
            item_w = fm.horizontalAdvance(name) + 25 
            total_width += item_w
        
        # Set a reasonable maximum width so it doesn't take up the whole bar 
        # unless there are many items. This pushes the + Web button to the left.
        self.bookmarks_list.setMaximumWidth(total_width + 10)

    def on_bookmark_clicked(self, item):
        url = item.data(Qt.UserRole)
        self.browser.setUrl(QUrl(url))
        self.url_bar.setText(url)

    def on_bookmarks_reordered(self, parent, start, end, destination, row):
        # Update self.chat_services base on new order
        new_services = {}
        for i in range(self.bookmarks_list.count()):
            item = self.bookmarks_list.item(i)
            name = item.text()
            url = item.data(Qt.UserRole)
            new_services[name] = url
        self.chat_services = new_services
        self.save_config()

    def on_bookmarks_context_menu(self, pos):
        item = self.bookmarks_list.itemAt(pos)
        if not item: return
        name = item.text()
        url = item.data(Qt.UserRole)

        menu = QMenu()
        edit_action = menu.addAction(self.T("edit_fmt").format(name))
        delete_action = menu.addAction(self.T("delete_fmt").format(name))
        
        action = menu.exec(self.bookmarks_list.mapToGlobal(pos))
        
        if action == edit_action:
            self.show_add_service_dialog(edit_mode=True, old_name=name, old_url=url)
        elif action == delete_action:
            if len(self.chat_services) <= 1:
                QMessageBox.warning(self, self.T("warning"), self.T("cannot_delete_last"))
                return
            confirm = QMessageBox.question(self, self.T("confirm_delete"), self.T("delete_confirmation").format(name), QMessageBox.Yes | QMessageBox.No)
            if confirm == QMessageBox.Yes:
                del self.chat_services[name]
                if name in self.all_possible_services:
                    del self.all_possible_services[name]
                if self.default_service == name:
                    self.default_service = list(self.chat_services.keys())[0]
                self.sync_and_save()

    def sync_and_save(self):
        """Update UI elements and save to config.json"""
        self.default_page_combo.blockSignals(True)
        self.default_page_combo.clear()
        services = list(self.chat_services.keys())
        self.default_page_combo.addItems(services)
        self.default_page_combo.addItem(self.T("web_plus"))

        if self.default_service in services:
             self.default_page_combo.setCurrentText(self.default_service)
        else:
             if services:
                self.default_service = services[0]
                self.default_page_combo.setCurrentText(self.default_service)
             
        self.default_page_combo.blockSignals(False)
        
        self.refresh_service_buttons()
        self.save_config()

    def on_startup_combo_context_menu(self, pos):
        name = self.default_page_combo.currentText()
        if not name or name not in self.chat_services: return
        url = self.chat_services[name]

        menu = QMenu()
        edit_action = menu.addAction(self.T("edit_fmt").format(name))
        delete_action = menu.addAction(self.T("delete_fmt").format(name))
        
        action = menu.exec(self.default_page_combo.mapToGlobal(pos))
        
        if action == edit_action:
            self.show_add_service_dialog(edit_mode=True, old_name=name, old_url=url)
        elif action == delete_action:
            if len(self.chat_services) <= 1:
                QMessageBox.warning(self, self.T("warning"), self.T("cannot_delete_last"))
                return
            confirm = QMessageBox.question(self, self.T("confirm_delete"), self.T("delete_confirmation").format(name), QMessageBox.Yes | QMessageBox.No)
            if confirm == QMessageBox.Yes:
                del self.chat_services[name]
                if self.default_service == name:
                    self.default_service = list(self.chat_services.keys())[0]
                self.sync_and_save()

    def show_add_service_dialog(self, edit_mode=False, old_name="", old_url=""):
        dialog = QDialog(self)
        dialog.setWindowTitle(self.T("edit_web_link") if edit_mode else self.T("add_web_link"))
        layout = QFormLayout(dialog)
        
        name_edit = QLineEdit(old_name)
        url_edit = QLineEdit(old_url if edit_mode else "https://")
        
        layout.addRow(self.T("web_name"), name_edit)
        layout.addRow(self.T("url"), url_edit)
        
        btn_box = QHBoxLayout()
        add_btn = QPushButton(self.T("save") if edit_mode else self.T("add"))
        cancel_btn = QPushButton(self.T("cancel"))
        btn_box.addWidget(add_btn)
        btn_box.addWidget(cancel_btn)
        layout.addRow(btn_box)
        
        def do_add():
            name = name_edit.text().strip()
            url = url_edit.text().strip()
            if name and url.startswith("http"):
                # If editing and name changed, remove old
                if edit_mode and old_name != name:
                    if old_name in self.chat_services: del self.chat_services[old_name]
                    if old_name in self.all_possible_services: del self.all_possible_services[old_name]
                
                self.chat_services[name] = url
                self.all_possible_services[name] = url
                
                self.default_service = name
                self.sync_and_save()
                dialog.accept()
            else:
                QMessageBox.warning(dialog, self.T("input_error"), self.T("input_error_msg"))
        
        add_btn.clicked.connect(do_add)
        cancel_btn.clicked.connect(dialog.reject)
        dialog.exec()

    def toggle_topmost(self, state):
        # Deprecated logic kept for compatibility if needed, but UI removed
        pass

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
        if files:
            self.path_edit.setText(files[0])

    def navigate_to_url(self):
        url = self.url_bar.text()
        if not url.startswith("http"):
            url = "https://" + url
        self.browser.setUrl(QUrl(url))

    def generate_prompt(self):
        path = self.path_edit.text().strip()
        if not os.path.exists(path):
            QMessageBox.critical(self, "Error", "File not found!")
            return
        
        try:
            source = "scanpy" if self.rb_scanpy.isChecked() else "seurat"
            mode = "concise" if self.rb_concise.isChecked() else "detailed"
            
            prompt = annotate_cell_types(
                marker_file=path,
                step="generate",
                species=self.species_edit.text(),
                tissue=self.tissue_edit.text(),
                top_n=int(self.top_n_edit.text()),
                mode=mode,
                source=source,
                exclude_types=self.exclude_edit.text()
            )
            self.prompt_display.setText(prompt)
            pyperclip.copy(prompt)
            QMessageBox.information(self, "Success", "Prompt copied! Paste it into the browser on the right.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed: {str(e)}")

    def process_response(self):
        path = self.path_edit.text().strip()
        response = self.response_input.toPlainText().strip()
        
        if not response:
            QMessageBox.warning(self, "Warning", "Please paste the AI response first.")
            return

        # If no file path, ask user for confirmation to use text-only parsing
        if not path:
            reply = QMessageBox.question(
                self, 
                "No Marker File", 
                "No marker file selected.\n\nDo you want to generate code using ONLY the AI response? (Ensure the AI output contains Cluster IDs)",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
            
        try:
            source = "scanpy" if self.rb_scanpy.isChecked() else "seurat"
            # Now only using the code part
            _, code = annotate_cell_types(
                marker_file=path,
                step="parse",
                response_text=response,
                top_n=int(self.top_n_edit.text()),
                source=source
            )
            
            self.code_display.setText(code)
            QMessageBox.information(self, "Success", "Processing complete! Check the code output tab.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Parsing failed: {str(e)}")

    def copy_code(self):
        pyperclip.copy(self.code_display.toPlainText())
        QMessageBox.information(self, "Copied", "Code copied to clipboard.")

    def show_help(self):
        msg = """
        ChatCellAnno Guide:
        1. Select your marker file (Scanpy/Seurat).
        2. Set Parameters (Species, Tissue, Exclude list).
        3. Click 'Generate Prompt'. It's auto-copied!
        4. Use the internal browser on the right. Login to ChatGPT/DeepSeek.
        5. Paste (Ctrl+V) and Send.
        6. Copy AI's table response and paste it into Step 4 box (bottom left).
        7. Click 'Process AI Output' to see results and get code.
        """
        QMessageBox.information(self, "User Guide", msg.strip())

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ChatCellAnnoApp()
    window.show()
    sys.exit(app.exec())
