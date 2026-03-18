import sys
import os
import pyperclip
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QLabel, QLineEdit, QPushButton, QRadioButton, QButtonGroup, 
    QGroupBox, QTextEdit, QFileDialog, QMessageBox, QScrollArea,
    QSplitter, QComboBox, QDialog, QFormLayout, QToolButton,
    QListWidget, QListWidgetItem, QAbstractItemView, QMenu, QWidgetAction,
    QTabWidget, QProgressBar, QApplication
)
from PySide6.QtGui import QIcon, QPixmap, QImage, QImageReader
from PySide6.QtCore import Qt, QUrl, QMimeData
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebEngineCore import QWebEngineProfile, QWebEnginePage, QWebEngineSettings

# Relative imports within the package structure
# Assuming this file is at chatcellanno/gui/main_window.py
# So we go up to chatcellanno then down to config
from ..config import ConfigManager
from .workers import EnrichmentWorker
from ..genome_utils import find_gene_sequences, format_mega_style

class DragOutLabel(QLabel):
    """
    文件拖拽输出控件 (File Drag-Out Component).
    
    该类继承自 QLabel，允许用户将生成的提示词文件或报告直接从应用程序界面
    拖拽到外部环境（如浏览器窗口、文件夹或桌面）。
    """
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.file_paths = [] # 存储可供拖拽的文件路径列表

    def setFiles(self, paths):
        """设置当前关联的文件路径 (Set associated file paths)"""
        self.file_paths = paths

    def mousePressEvent(self, event):
        """鼠标按下事件：触发系统的 QDrag 拖拽操作 (Trigger system drag-and-drop)"""
        if event.button() == Qt.LeftButton and self.file_paths:
            from PySide6.QtGui import QDrag
            from PySide6.QtCore import QMimeData, QUrl
            drag = QDrag(self)
            mime = QMimeData()
            # 将本地路径转换为 URL 格式 (Convert local paths to QUrls)
            urls = [QUrl.fromLocalFile(p) for p in self.file_paths if os.path.exists(p)]
            if urls:
                mime.setUrls(urls)
                drag.setMimeData(mime)
                # 执行“复制”类型的拖拽 (Execute copy-action drag)
                drag.exec_(Qt.CopyAction)

class MainWindow(QMainWindow):
    """
    主窗体类 (Main Application Window).
    
    这是 ChatCellAnno 的核心 UI 控制器，负责：
    1. 管理多步骤的单细胞分析流 (Marker 提取 -> 富集 -> AI 对话 -> 代码生成)。
    2. 集成 QtWebEngine 浏览器以便与在线 LLM 交互。
    3. 管理 API 调用、插件系统以及本地数据库路径。
    """
    def __init__(self):
        super().__init__()
        
        # 初始化配置管理器 (Initialize global config)
        self.config = ConfigManager()
        
        # 路径资产管理 (Resource Path Management)
        # 自动识别是运行在 PyInstaller 打包后的环境还是原始 Python 环境
        if hasattr(sys, "_MEIPASS"):
            base_dir = sys._MEIPASS # 打包后的临时解压路径
        else:
            base_dir = os.path.dirname(os.path.abspath(sys.argv[0])) # 脚本所在目录
            
        # 数据库目录校验逻辑 (Database directory validation)
        db_check = os.path.join(base_dir, "database")
        if not os.path.exists(db_check):
             # 针对源码开发模式的路径修正 (Path correction for source dev mode)
             base_dir_pkg = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
             if os.path.exists(os.path.join(base_dir_pkg, "database")):
                 base_dir = base_dir_pkg

        self.db_dir = os.path.join(base_dir, "database") # 通路富集数据库目录
        self.genome_dir = os.path.join(base_dir, "genome") # 基因对比序列目录
        self.db_path_map = {} 

        # 默认在线数据库列表 (Online DB identifiers for Enrichr)
        self.online_db_list = [
            "GO_Biological_Process_2025",
            "GO_Cellular_Component_2025",
            "GO_Molecular_Function_2025",
            "KEGG_2026",
            "CellMarker_2024",
            "PanglaoDB_Augmented_2021"
        ]
        
        # 配置 WebEngine 存储持久化 (Configure WebEngine persistence)
        # 确保浏览器中的 AI 账号登录状态、Cookie 不会因为应用关闭而丢失
        self.web_profile = QWebEngineProfile("ChatCellAnnoProfile", self)
        self.web_profile.setPersistentStoragePath(self.config.config_paths["web_storage"])
        self.web_profile.setCachePath(self.config.config_paths["web_storage"])
        
        # 应用状态变量 (Application State Variables)
        self.current_marker_file = "" # 当前加载的标记基因文件
        self.enrichment_data = None    # 存储最近一次富集分析的结果
        self.all_possible_services = self.config.chat_services.copy() # AI 服务列表
        
        # 设置窗体基础属性 (Set basic window properties)
        self.setWindowTitle(self.config.T("title"))
        if os.path.exists("app_icon.ico"):
            self.setWindowIcon(QIcon("app_icon.ico"))

        # 全局样式表定义 (Global Stylesheet Architecture)
        # 模仿现代桌面 UI 设计风格 (Apple/Windows 11 Style)
        self.setStyleSheet("""
            QMainWindow { background-color: #f5f5f7; }
            QWidget { font-family: "Segoe UI", "Microsoft YaHei", sans-serif; font-size: 10pt; }
            QGroupBox {
                background-color: white; border: 1px solid #dcdcdc;
                border-radius: 8px; margin-top: 10px; padding-top: 15px; 
            }
            QGroupBox::title {
                subcontrol-origin: margin; subcontrol-position: top left;
                padding: 0 5px; color: #333; font-weight: bold; font-size: 15px;
            }
            QPushButton {
                background-color: #007aff; color: white; 
                border-radius: 5px; padding: 6px 12px; border: none;
            }
            QPushButton:hover { background-color: #0060c0; }
            QPushButton:pressed { background-color: #004080; }
            QLineEdit, QTextEdit, QComboBox {
                border: 1px solid #ccc; border-radius: 4px;
                padding: 4px; background-color: white; selection-background-color: #007aff;
            }
            QLineEdit:focus, QTextEdit:focus, QComboBox:focus { border: 1px solid #007aff; }
            QListWidget { border: none; background-color: transparent; }
            QListWidget::item {
                background-color: white; border: 1px solid #ddd;
                border-radius: 4px; margin: 2px; padding: 5px; color: #333;
            }
            QListWidget::item:hover { background-color: #f0f8ff; border-color: #007aff; }
            QScrollArea { border: none; background-color: transparent; }
        """)

        # 启动 UI 构建 (Start UI Setup)
        self.setup_ui()
        # 允许窗体接受拖放事件 (Enable drag-and-drop for file inputs)
        self.setAcceptDrops(True) 
        
        # 设置左右面板拉伸比例为 1:1 (Set balanced 1:1 splitter stretch)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 1)

    def setup_ui(self):
        """
        组装程序主界面 (Assemble the Main UI Layout).
        采用 QSplitter 分为左侧控制台和右侧浏览器/结果区。
        """
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)

        # 使用分割器实现左右可调节布局 (Use splitter for adjustable layouts)
        self.splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(self.splitter)

        # --- 左侧面板：交互控制区 (Left Panel: Controls) ---
        left_container = QWidget()
        left_layout = QVBoxLayout(left_container)
        
        # 滚动区域支持，防止小屏幕下内容溢出 (Scroll area for responsiveness)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        left_layout.addWidget(scroll)
        
        control_widget = QWidget()
        control_layout = QVBoxLayout(control_widget)
        scroll.setWidget(control_widget)

        # 1. 数据选择与多模态输入模块 (Data Selection & Multimodal Inputs)
        # 引用重构后的 UI Block (Importing refactored modular UI blocks)
        from .ui_blocks import build_data_selection_ui
        data_group = build_data_selection_ui(self)
        control_layout.addWidget(data_group)



        # 3. Prompt Configuration (Step 2)
        gen_group = QGroupBox(self.config.T("step2"))
        gen_layout = QVBoxLayout()
        
        # 3.1 Basic Settings
        self.species_edit = QLineEdit(self.config.session_params.get("species", "Human"))
        self.tissue_edit = QLineEdit(self.config.session_params.get("tissue", "PBMC"))
        self.top_n_edit = QLineEdit(self.config.session_params.get("top_n", "10"))
        
        # 信号连接：当用户修改物种、组织或 Top N 时，自动保存状态到当前会话 (Update config state on text change)
        self.species_edit.textChanged.connect(self.save_config_state)
        self.tissue_edit.textChanged.connect(self.save_config_state)
        self.top_n_edit.textChanged.connect(self.save_config_state)

        # 布局：物种与组织输入 (Species & Tissue inputs)
        gen_layout.addWidget(QLabel(self.config.T("species")))
        gen_layout.addWidget(self.species_edit)
        gen_layout.addWidget(QLabel(self.config.T("tissue")))
        gen_layout.addWidget(self.tissue_edit)
        gen_layout.addWidget(QLabel(self.config.T("top_n")))
        gen_layout.addWidget(self.top_n_edit)
        
        # --- 交互策略配置 (Query Strategy Setup) ---
        # 排除特定的已知细胞类型，以引导 AI 产出更准确的结果 (Exclude specific known cell types)
        self.exclude_edit = QLineEdit(self.config.session_params.get("exclude", ""))
        self.exclude_edit.setPlaceholderText("e.g. Neuron, Astrocyte")
        self.exclude_edit.textChanged.connect(self.save_config_state)
        
        gen_layout.addWidget(QLabel(self.config.T("exclude")))
        gen_layout.addWidget(self.exclude_edit)

        # 提示词精简度选择 (Prompt detail mode: Concise vs Detailed)
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

        # 3.2 交互模式切换：浏览器 vs 直接 API (Hybrid Mode: Browser vs API Direct)
        gen_layout.addWidget(QLabel(self.config.T("query_mode")))
        self.query_mode_group = QButtonGroup(self)
        query_mode_layout = QHBoxLayout()
        self.rb_query_browser = QRadioButton(self.config.T("browser_query")) # 浏览器手动对话 (Manual Chat)
        self.rb_query_api = QRadioButton(self.config.T("api_query"))         # API 自动请求 (API Automation)
        self.rb_query_browser.setChecked(True)
        self.query_mode_group.addButton(self.rb_query_browser)
        self.query_mode_group.addButton(self.rb_query_api)
        query_mode_layout.addWidget(self.rb_query_browser)
        query_mode_layout.addWidget(self.rb_query_api)
        gen_layout.addLayout(query_mode_layout)

        # 3.3 动态交互控制区 (Dynamic Action Area)
        # 包含两种模式下的按钮和拖放区域，会根据 Query Mode 自动切换显示 (Visibility toggled dynamically)
        
        # 浏览器模式下的操作区：生成提示词文件并提供“拖拽输出” (Browser Mode Actions)
        self.browser_action_widget = QWidget()
        browser_layout = QVBoxLayout(self.browser_action_widget)
        browser_layout.setContentsMargins(0, 0, 0, 0)
        
        self.gen_btn = QPushButton(self.config.T("gen_prompt_file"))
        self.gen_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; height: 35px;")
        self.gen_btn.clicked.connect(self.generate_prompt)
        browser_layout.addWidget(self.gen_btn)
        
        # 拖拽区：生成的 Prompt 文件可以通过此 Label 拖拽到网页版 AI 对话框 (Drag files directly to Chat)
        self.drag_area_group = QGroupBox(self.config.T("drag_info"))
        self.drag_area_layout = QVBoxLayout()
        self.drag_area_info = DragOutLabel(self.config.T("drag_info"))
        self.drag_area_info.setAlignment(Qt.AlignCenter)
        self.drag_area_info.setStyleSheet("border: 2px dashed #4CAF50; border-radius: 5px; padding: 15px; background: #f0fff0; color: #2e7d32; font-weight: bold;")
        self.drag_area_info.setMinimumHeight(100)
        self.drag_area_layout.addWidget(self.drag_area_info)
        self.drag_area_group.setLayout(self.drag_area_layout)
        self.drag_area_group.hide() # 初始隐藏，文件生成后显示 (Hidden initially)
        browser_layout.addWidget(self.drag_area_group)
        
        # API 模式下的操作区：直连推理服务器 (API Mode Actions)
        self.api_action_widget = QWidget()
        api_layout = QHBoxLayout(self.api_action_widget)
        api_layout.setContentsMargins(0, 0, 0, 0)
        
        self.btn_call_api = QPushButton("🚀 " + self.config.T("query_btn"))
        self.btn_call_api.setStyleSheet("background-color: #673AB7; color: white; font-weight: bold; height: 35px;")
        self.btn_call_api.clicked.connect(self.call_api_directly)
        
        self.btn_api_settings = QPushButton("⚙️")
        self.btn_api_settings.clicked.connect(self.show_api_settings_dialog)
        self.btn_api_settings.setFixedWidth(40)
        self.btn_api_settings.setToolTip(self.config.T("api_settings"))
        
        api_layout.addWidget(self.btn_call_api)
        api_layout.addWidget(self.btn_api_settings)
        
        gen_layout.addWidget(self.browser_action_widget)
        gen_layout.addWidget(self.api_action_widget)

        # 动态切换控制逻辑 (Toggle UI visibility based on Query Mode radio buttons)
        def on_query_mode_changed():
            is_api = self.rb_query_api.isChecked()
            self.api_action_widget.setVisible(is_api)
            self.browser_action_widget.setVisible(not is_api)

        self.rb_query_browser.toggled.connect(on_query_mode_changed)
        self.rb_query_api.toggled.connect(on_query_mode_changed)
        on_query_mode_changed()

        # 提示词预览显示区 (Live Prompt Preview Window)
        self.prompt_display = QTextEdit()
        self.prompt_display.setMaximumHeight(80)
        self.prompt_display.setFontFamily("Consolas")
        
        gen_layout.addWidget(self.prompt_display)
        gen_group.setLayout(gen_layout)
        control_layout.addWidget(gen_group)

        # 4. 析 AI 回复模块 (Step 3: AI Response Processing)
        parse_group = QGroupBox(self.config.T("step3"))
        parse_layout = QVBoxLayout()
        
        # 用户在此粘贴 AI 生成的 Markdown 表格 (Paste AI response here)
        self.response_input = QTextEdit()
        self.response_input.setPlaceholderText(self.config.T("placeholder_resp"))
        self.response_input.setMaximumHeight(80)
        
        # 执行脚本解析与映射按钮 (Execution button for parser logic)
        parse_btn = QPushButton(self.config.T("process_btn"))
        parse_btn.setStyleSheet("background-color: #2196F3; color: white; font-weight: bold; height: 30px;")
        parse_btn.clicked.connect(self.process_response)
        
        clean_btn = QPushButton(self.config.T("clean"))
        clean_btn.setFixedWidth(60)
        clean_btn.clicked.connect(lambda: self.response_input.clear())
        
        # 目标代码框架选择 (Target analysis platform selector)
        self.source_group = QButtonGroup(self)
        platform_layout = QHBoxLayout()
        self.rb_scanpy = QRadioButton(self.config.T("scanpy")) # Python 生态
        self.rb_scanpy.setChecked(True)
        self.rb_seurat = QRadioButton(self.config.T("seurat")) # R 语言生态
        self.source_group.addButton(self.rb_scanpy)
        self.source_group.addButton(self.rb_seurat)
        platform_layout.addWidget(QLabel("Target Platform:"))
        platform_layout.addWidget(self.rb_scanpy)
        platform_layout.addWidget(self.rb_seurat)
        
        btn_row = QHBoxLayout()
        btn_row.addWidget(parse_btn)
        btn_row.addWidget(clean_btn)
        
        parse_layout.addWidget(self.response_input)
        parse_layout.addLayout(platform_layout)
        parse_layout.addLayout(btn_row)
        parse_group.setLayout(parse_layout)
        control_layout.addWidget(parse_group)

        # 5. 代码导出与报告生成模块 (Step 4: Code Export & Final Report)
        export_group = QGroupBox(self.config.T("step4"))
        export_layout = QVBoxLayout()

        # 显示生成的注释映射代码 (Display mapping code snippets)
        self.code_display = QTextEdit()
        self.code_display.setReadOnly(True)
        self.code_display.setFontFamily("Consolas")
        export_layout.addWidget(QLabel(self.config.T("code_label")))
        export_layout.addWidget(self.code_display)
        
        copy_code_btn = QPushButton(self.config.T("copy_code"))
        copy_code_btn.clicked.connect(self.copy_code)
        
        # 一键导出 Markdown 格式的详细分析报告 (Export Analysis Report button)
        self.gen_report_btn = QPushButton("📝 " + self.config.T("gen_report"))
        self.gen_btn_style = "background-color: #FF9800; color: white; font-weight: bold; height: 35px;"
        self.gen_report_btn.setStyleSheet(self.gen_btn_style)
        self.gen_report_btn.clicked.connect(self.export_analysis_report)
        
        btns_layout = QHBoxLayout()
        btns_layout.addWidget(copy_code_btn)
        btns_layout.addWidget(self.gen_report_btn)
        export_layout.addLayout(btns_layout)
        
        export_group.setLayout(export_layout)
        control_layout.addWidget(export_group)

        # --- 右侧面板：内嵌浏览器与数据展示区 (Right Panel: Web & Visualization Area) ---
        right_panel_container = QWidget()
        right_panel_layout = QVBoxLayout(right_panel_container)
        
        # 顶部工具条：包含置顶、设置等按钮 (Top control toolbar)
        toolbar_top = QHBoxLayout()
        toolbar_top.addStretch()

        # 窗口置顶切换 (Stay-on-top toggle)
        self.pin_btn = QToolButton()
        self.pin_btn.setText("📌")
        self.pin_btn.setCheckable(True)
        self.pin_btn.setToolTip("Always on Top")
        self.pin_btn.clicked.connect(self.toggle_topmost_btn)
        toolbar_top.addWidget(self.pin_btn)

        # 全局设置菜单按钮 (Global settings gear-menu)
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

        # --- 配置多语言切换菜单 (Multilingual Support Menu) ---
        # 支持在运行时动态切换界面语言，无需手动修改配置文件
        lang_menu = QMenu(self.config.T("language"), self)
        zh_action = lang_menu.addAction(self.config.T("lang_zh")) # 中文
        zh_action.triggered.connect(lambda: self.set_language("zh"))
        en_action = lang_menu.addAction(self.config.T("lang_en")) # 英文 (English)
        en_action.triggered.connect(lambda: self.set_language("en"))
        
        self.gear_menu.addMenu(lang_menu)

        self.gear_menu.addSeparator()

        # 添加面板控制菜单
        panels_menu = QMenu(self.config.T("panels_control"), self)
        
        enrichment_action = panels_menu.addAction(self.config.T("tab_enrichment"))
        enrichment_action.triggered.connect(lambda: self.show_panel_tab(self.tab_enrichment_widget, "tab_enrichment"))
        
        browser_action = panels_menu.addAction(self.config.T("tab_browser"))
        browser_action.triggered.connect(lambda: self.show_panel_tab(self.tab_browser_widget, "tab_browser"))
        
        genome_action = panels_menu.addAction(self.config.T("tab_genome"))
        genome_action.triggered.connect(lambda: self.show_panel_tab(self.tab_genome_widget, "tab_genome"))

        self.gear_menu.addMenu(panels_menu)

        self.gear_menu.addSeparator()

        # 添加用户指南/帮助文档入口 (User Guide Action)
        guide_action = self.gear_menu.addAction(self.config.T("guide"))
        guide_action.triggered.connect(self.show_help)
        
        # --- 注入可扩展插件系统 (Plugin System Architecture Injection) ---
        # 该部分实现了软件的“热插拔”特性，允许在不重启的情况下加载外部功能模块
        self.gear_menu.addSeparator()
        
        def show_plugin_manager():
            """弹出插件管理对话框，允许用户安装或管理 .py 插件 (Show Plugin Manager Dialog)"""
            from chatcellanno.plugins.manager_gui import PluginManagerDialog
            # 检查 plugin_manager 实例是否存在 (Verify plugin manager instance)
            dlg = PluginManagerDialog(self, hasattr(self, 'plugin_manager') and self.plugin_manager)
            dlg.exec()
            
        # 安全获取语言包中的插件管理字样 (Safely fetch i18n label)
        try:
            plugin_text = self.config.T("plugin_manager")
        except:
            plugin_text = "插件管理 (Plugins)"
            
        plugin_action = self.gear_menu.addAction("🧩 " + plugin_text)
        plugin_action.triggered.connect(show_plugin_manager)
        
        # 尝试初始化插件管理器并自动加载已存在的插件 (Initialize and auto-load plugins)
        try:
            from chatcellanno.plugins import PluginManager
            # 建立主窗体与插件管理器的双向绑定关系 (Bi-directional binding)
            self.plugin_manager = PluginManager(self)
            self.plugin_manager.load_plugins() # 扫描 plugins 目录
            self.plugin_manager.init_plugins() # 执行插件的初始化钩子 (Initial hooks)
        except Exception as e:
            # 插件加载失败不应导致主程序崩溃，仅在控制台记录错误 (Log error without crashing)
            print(f"Failed to load plugins: {e}")
        # -------------------------------
        
        self.gear_btn.setMenu(self.gear_menu)
        self.gear_btn.setPopupMode(QToolButton.InstantPopup)
        toolbar_top.addWidget(self.gear_btn)
        
        right_panel_layout.addLayout(toolbar_top)

        # --- 右侧功能标签页系统 (Right-side Tabbed Dashboard System) ---
        # 包含：1. 富集结果看板 2. AI 智能浏览器 3. 基因序列比对中心
        self.right_tabs = QTabWidget()
        self.right_tabs.setTabsClosable(True)
        self.right_tabs.setMovable(True)
        self.right_tabs.tabCloseRequested.connect(self.close_right_tab)
        
        # 标签页 1：功能富集可视化区 (Tab 1: Functional Enrichment Visualization)
        enrichment_container = QWidget()
        enrichment_layout = QVBoxLayout(enrichment_container)

        # 任务进度状态条 (Status Container for active enrichment tasks)
        self.enrich_status_container = QWidget()
        status_layout = QVBoxLayout(self.enrich_status_container)
        self.status_label = QLabel(self.config.T("enrich_running"))
        self.status_label.setStyleSheet("color: #673AB7; font-weight: bold;")
        self.enrich_progress = QProgressBar() # 进度条反馈 (Real-time progress feedback)
        self.enrich_progress.setRange(0, 100)
        status_layout.addWidget(self.status_label)
        status_layout.addWidget(self.enrich_progress)
        self.enrich_status_container.hide() # 默认隐藏，在 worker 运行时开启 (Visible during task)
        enrichment_layout.addWidget(self.enrich_status_container)
        
        # 富集详细结果滚动区 (Scrollable area for detailed term results)
        enrich_scroll = QScrollArea()
        enrich_scroll.setWidgetResizable(True)
        enrich_scroll_content = QWidget()
        self.enrich_detail_layout = QVBoxLayout(enrich_scroll_content)
        
        # 初始占位提示文字 (Initial placeholder text)
        self.enrich_placeholder = QLabel("No enrichment results yet. Generate a prompt with enrichment enabled to see data here.")
        self.enrich_placeholder.setWordWrap(True)
        self.enrich_placeholder.setAlignment(Qt.AlignCenter)
        self.enrich_detail_layout.addWidget(self.enrich_placeholder)
        
        enrich_scroll.setWidget(enrich_scroll_content)
        enrichment_layout.addWidget(enrich_scroll)
        
        self.tab_enrichment_widget = enrichment_container
        self.right_tabs.addTab(enrichment_container, self.config.T("tab_enrichment"))

        # 标签页 2：AI 智能对话浏览器 (Tab 2: Integrated AI Chat Browser)
        # 该模块核心在于 QWebEngineView，通过它用户可以直接与 ChatGPT/DeepSeek 对话
        browser_container = QWidget()
        browser_layout = QVBoxLayout(browser_container)

        # 确定初始加载页面 (Determination of initial URL from config)
        initial_url = self.config.chat_services.get(self.config.default_service, "https://chat.deepseek.com")
        
        # 浏览器地址栏区域 (Browser Navigation Toolbar)
        browser_header = QHBoxLayout()
        self.url_bar = QLineEdit(initial_url)
        self.url_bar.returnPressed.connect(self.navigate_to_url) # 回车触发表单提交/跳转
        go_btn = QPushButton(self.config.T("go"))
        go_btn.clicked.connect(self.navigate_to_url)
        
        # 刷新页面按钮 (Reload button)
        refresh_btn = QPushButton(self.config.T("refresh"))
        refresh_btn.clicked.connect(lambda: self.browser.reload())

        browser_header.addWidget(QLabel(self.config.T("browser_label")))
        browser_header.addWidget(self.url_bar)
        browser_header.addWidget(go_btn)
        browser_header.addWidget(refresh_btn)
        
        # 书签/快捷服务栏 (Bookmarks / Quick Service Bar)
        # 允许用户通过单击快速切换不同的 AI 服务提供商 (Quick switch between provider models)
        bookmarks_bar_layout = QHBoxLayout()

        self.bookmarks_list = QListWidget()
        self.bookmarks_list.setFlow(QListWidget.LeftToRight) # 水平排列布局 (Horizontal list layout)
        self.bookmarks_list.setFixedHeight(35)
        self.bookmarks_list.setDragDropMode(QAbstractItemView.InternalMove) # 允许用户手动拖拽调整书签顺序
        self.bookmarks_list.setWrapping(False)
        self.bookmarks_list.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.bookmarks_list.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.bookmarks_list.setContextMenuPolicy(Qt.CustomContextMenu) # 支持右键菜单 (Custom right-click menu)
        
        # 连接书签相关的信号槽 (Signal-slot connections for bookmarks)
        self.bookmarks_list.customContextMenuRequested.connect(self.on_bookmarks_context_menu)
        self.bookmarks_list.itemClicked.connect(self.on_bookmark_clicked)
        self.bookmarks_list.model().rowsMoved.connect(self.on_bookmarks_reordered)
        
        # 点击空白处清除选中的辅助逻辑 (Auxiliary logic to clear selection)
        def clear_selection(event):
            item = self.bookmarks_list.itemAt(event.position().toPoint())
            if not item:
                self.bookmarks_list.clearSelection()
            QListWidget.mousePressEvent(self.bookmarks_list, event)
        self.bookmarks_list.mousePressEvent = clear_selection

        # 书签项的美化样式设置 (Bookmark item styling)
        self.bookmarks_list.setStyleSheet("""
            QListWidget { background: transparent; border: none; }
            QListWidget::item { 
                background: #f0f0f0; border: 1px solid #ccc; 
                border-radius: 4px; margin: 2px; padding: 4px 8px; color: #333; 
            }
            QListWidget::item:hover { background: #e0e0e0; }
            QListWidget::item:selected { background: #d0d0d0; border-color: #999; }
        """)
        
        # 从配置文件加载并刷新书签列表 (Load and refresh bookmarks from config manager)
        self.refresh_service_buttons()

        # 添加新网页/书签的快捷入口 (Quick add button for new services)
        add_web_btn = QPushButton(self.config.T("add_web"))
        add_web_btn.setFixedWidth(60)
        add_web_btn.clicked.connect(self.show_add_service_dialog)
        
        bookmarks_bar_layout.addWidget(self.bookmarks_list)
        bookmarks_bar_layout.addSpacing(5) 
        bookmarks_bar_layout.addWidget(add_web_btn)
        bookmarks_bar_layout.addStretch() 
        
        # --- 核心浏览器引擎组件 (Core Web Engine Component) ---
        # 集成安全沙箱机制，允许用户在应用内完成复杂的 AI 登录和对话
        self.browser = QWebEngineView()
        # 将浏览器页面绑定到具有持久性存储的 WebEngineProfile (Bind page to persistent profile)
        page = QWebEnginePage(self.web_profile, self.browser)
        self.browser.setPage(page)
        
        # 核心设置：开启 Javascript 访问剪贴板权限，以配合系统级的复制工作流 (Crucial: Enable Clipboard access)
        self.browser.settings().setAttribute(QWebEngineSettings.WebAttribute.JavascriptCanAccessClipboard, True)
        self.browser.settings().setAttribute(QWebEngineSettings.WebAttribute.JavascriptCanOpenWindows, True)
        self.browser.setUrl(QUrl(initial_url))
        
        # 权限请求处理回调 (Handle feature permission requests like clipboard reading)
        page.featurePermissionRequested.connect(self.on_permission_requested)
        
        browser_layout.addLayout(browser_header)
        browser_layout.addLayout(bookmarks_bar_layout)
        browser_layout.addWidget(self.browser)
        
        self.tab_browser_widget = browser_container
        self.right_tabs.addTab(browser_container, self.config.T("tab_browser"))

        # --- 标签页 3：基因对比中心 (Tab 3: Genome Comparison / Multiple Sequence Alignment) ---
        # 本模块支持在各物种间搜索同源基因并以 MEGA 风格配色展示序列
        genome_container = QWidget()
        genome_layout = QVBoxLayout(genome_container)
        
        # 基因搜索栏架构 (Gene Search Bar Architecture)
        genome_search_layout = QHBoxLayout()
        self.genome_search_edit = QLineEdit()
        self.genome_search_edit.setPlaceholderText(self.config.T("placeholder_gene"))
        self.genome_search_edit.returnPressed.connect(self.search_genome) # 回车直接搜索 (Enter triggers search)
        genome_search_btn = QPushButton(self.config.T("search"))
        genome_search_btn.clicked.connect(self.search_genome)
        genome_search_layout.addWidget(self.genome_search_edit)
        genome_search_layout.addWidget(genome_search_btn)
        genome_layout.addLayout(genome_search_layout)
        
        # 自定义滚动展示区，用于动态生成比对卡片 (Dynamic area for alignment cards)
        self.genome_results_scroll = QScrollArea()
        self.genome_results_scroll.setWidgetResizable(True)
        self.genome_results_content = QWidget()
        self.genome_results_layout = QVBoxLayout(self.genome_results_content)
        self.genome_results_scroll.setWidget(self.genome_results_content)
        genome_layout.addWidget(self.genome_results_scroll)
        
        self.tab_genome_widget = genome_container
        self.right_tabs.addTab(genome_container, self.config.T("tab_genome"))

        # 将所有的 UI 结构体压入分割器 (Push components into the master splitter)
        right_panel_layout.addWidget(self.right_tabs)
        self.splitter.addWidget(left_container)
        self.splitter.addWidget(right_panel_container)
        
        # 禁用分割器折叠功能以保持 UI 稳定性 (Disable splitter collapse for layout stability)
        self.splitter.setCollapsible(0, False)
        self.splitter.setCollapsible(1, False)

        # 启动时的状态同步 (Final state sync and data saving on init completion)
        self.sync_and_save() 

    def on_permission_requested(self, url, feature):
        """
        处理 WebEngine 发起的系统级权限请求 (Handle browser permission requests).
        对于单细胞分析，最重要的权限是“剪贴板读写”以支持结果直传 AI。
        """
        if feature in (QWebEnginePage.PermissionType.ClipboardReadWrite, 
                       QWebEnginePage.PermissionType.ClipboardWrite):
            # 自动授予用户请求的剪贴板权限 (Automatically grant clipboard permissions)
            self.sender().setFeaturePermission(url, feature, QWebEnginePage.PermissionPolicy.PermissionGrantedByUser)
        else:
            # 对于摄像头、麦克风等不相关权限，默认不自动开启以保证隐私安全 (Privacy priority)
            self.sender().setFeaturePermission(url, feature, QWebEnginePage.PermissionPolicy.PermissionDeniedByUser)

    def update_enrich_db_list(self):
        """
        根据当前选择的模式（在线/本地）更新富集数据库下拉列表 (Update db list based on mode).
        """
        self.enrich_db_combo.clear()
        if self.rb_online.isChecked():
            # 在线模式：加载常用的公共基因集 (Load common public gene sets)
            self.enrich_db_combo.addItems(self.online_db_list)
            self.load_db_btn.hide()
        else:
            # 本地模式：显示用户已加载的本地 GMT 文件 (Show loaded local GMT files)
            if not self.db_path_map:
                self.enrich_db_combo.addItem(self.config.T("local_db_hint"))
            else:
                self.enrich_db_combo.addItems(list(self.db_path_map.keys()))
            self.load_db_btn.show()

    def save_config_state(self):
        """
        同步 UI 状态到配置对象并持久化到本地 JSON 文件 (Sync UI state to config object).
        涵盖了物种、组织、Top-N 基因数量和排除词等关键参数。
        """
        if hasattr(self, "species_edit"):
            # 将 UI 上的当前输入值存入 Python 持久化对象 (Capture UI inputs)
            self.config.session_params["species"] = self.species_edit.text()
            self.config.session_params["tissue"] = self.tissue_edit.text()
            self.config.session_params["top_n"] = self.top_n_edit.text()
            self.config.session_params["exclude"] = self.exclude_edit.text()
            # 记录用户选择的 Prompt 模式 (Record current prompting strategy mode)
            self.config.session_params["prompt_mode"] = "concise" if self.rb_concise.isChecked() else "detailed"
        self.config.save_config() # 写入 config.json


    def close_right_tab(self, index):
        self.right_tabs.removeTab(index)

    def show_panel_tab(self, widget, label_key):
        for i in range(self.right_tabs.count()):
            if self.right_tabs.widget(i) == widget:
                self.right_tabs.setCurrentIndex(i)
                return
        self.right_tabs.addTab(widget, self.config.T(label_key))
        self.right_tabs.setCurrentWidget(widget)

    def set_language(self, lang):
        """
        切换应用程序界面语言 (Switch interface language).
        :param lang: 'zh' 或 'en' (Language code)
        """
        if self.config.language == lang:
            return
        self.config.language = lang
        self.config.save_config()
        # 提示用户重启生效 (Notify user that restart is required for full i18n application)
        QMessageBox.information(self, self.config.T("restart_required"), self.config.T("language_changed").format(lang))

    def on_default_service_changed(self, text):
        """
        当用户修改默认启动页服务时触发的状态保存 (Handle change of startup service).
        """
        if text == self.config.T("web_plus"):
            # 如果选择“添加更多”，则弹出添加对话框 (Trigger 'Add New Service' dialog)
            self.show_add_service_dialog()
            self.default_page_combo.blockSignals(True)
            self.default_page_combo.setCurrentText(self.config.default_service)
            self.default_page_combo.blockSignals(False)
            return

        # 更新配置对象的默认值并保存 (Update and save default preference)
        self.config.default_service = text
        self.config.save_config()
    
    def toggle_topmost_btn(self, checked):
        """
        控制窗口是否始终置顶 (Toggle always-on-top window property).
        此功能便于用户一边查看数据文档，一边操作 AI 界面。
        """
        if checked:
            self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)
            self.pin_btn.setStyleSheet("background-color: #87CEEB; border-radius: 4px;") 
        else:
            self.setWindowFlags(self.windowFlags() & ~Qt.WindowStaysOnTopHint)
            self.pin_btn.setStyleSheet("background-color: transparent;") 
        # 修改 Flag 后必须重新调用 show() 才能在 Windows 上立即生效
        self.show()

    def refresh_service_buttons(self):
        """
        重新渲染浏览器书签栏的图标/文字按钮 (Re-render bookmark bar items).
        会根据配置中的 chat_services 字典进行循环构建。
        """
        self.bookmarks_list.clear() # 先清空旧的按钮 (Clear existing)
        total_width = 0
        fm = self.bookmarks_list.fontMetrics()
        for name, url in self.config.chat_services.items():
            item = QListWidgetItem(name)
            item.setData(Qt.UserRole, url) # 将 URL 隐藏在 UserRole 中 (Embed URL data)
            item.setTextAlignment(Qt.AlignCenter)
            self.bookmarks_list.addItem(item)
            # 动态计算按钮所需宽度 (Calculate width dynamically based on text length)
            item_w = fm.horizontalAdvance(name) + 25 
            total_width += item_w
        self.bookmarks_list.setMaximumWidth(total_width + 10)

    def on_bookmark_clicked(self, item):
        """
        响应点击书签的操作，让浏览器加载目标 URL (Load target URL on bookmark click).
        """
        url = item.data(Qt.UserRole)
        self.browser.setUrl(QUrl(url))
        self.url_bar.setText(url)

    def browse_db_folder(self):
        """
        打开文件夹选择器，让用户导入本地自定义 GMT 数据库 (Open directory picker for GMT db).
        """
        folder = QFileDialog.getExistingDirectory(self, "Select Database Folder")
        if folder:
            self.load_databases_from_folder(folder)

    def load_databases_from_folder(self, folder):
        """
        扫描指定文件夹中的 GMT/TXT 文件并加载为本地富集数据库 (Scan and load local databases).
        """
        if not os.path.exists(folder):
            return
            
        # 只保留符合基因集格式的文件扩展名 (Filter for gene set format extensions)
        gmt_files = [f for f in os.listdir(folder) if f.endswith(('.gmt', '.txt'))]
        if not gmt_files:
            return
            
        self.db_dir = folder
        self.config.last_db_path = folder  # 记录路径，下次启动自动定位 (Persistence for UX)
        self.config.save_config()
        
        # 将文件名映射到其路径，供下拉菜单使用 (Map filename to path for combo box)
        self.db_path_map = {f: f for f in gmt_files}
        self.custom_db_label.setText(f"Folder: ...{folder[-30:]}")
        self.rb_local.setChecked(True) # 自动切换到本地模式 (Auto-switch to local mode)
        self.update_enrich_db_list()

    def on_bookmarks_reordered(self, parent, start, end, destination, row):
        """
        当用户通过拖拽改变书签顺序时同步更新配置文件 (Sync config after drag-and-drop reordering).
        """
        new_services = {}
        for i in range(self.bookmarks_list.count()):
            item = self.bookmarks_list.item(i)
            name = item.text()
            url = item.data(Qt.UserRole)
            new_services[name] = url
        self.config.chat_services = new_services
        self.config.save_config() # 持久化新的排序顺序 (Persist new sequence)

    def on_bookmarks_context_menu(self, pos):
        """
        浏览器书签栏的右键菜单：支持编辑、删除操作 (Context menu for bookmarks).
        """
        item = self.bookmarks_list.itemAt(pos)
        if not item: return
        name = item.text()
        url = item.data(Qt.UserRole)

        menu = QMenu()
        edit_action = menu.addAction(self.config.T("edit_fmt").format(name))
        delete_action = menu.addAction(self.config.T("delete_fmt").format(name))
        
        action = menu.exec(self.bookmarks_list.mapToGlobal(pos))
        
        if action == edit_action:
            # 开启编辑模式弹出框 (Open dialog in edit mode)
            self.show_add_service_dialog(edit_mode=True, old_name=name, old_url=url)
        elif action == delete_action:
            # 防止用户删光所有服务导致的崩溃风险 (Prevent deleting all services)
            if len(self.config.chat_services) <= 1:
                QMessageBox.warning(self, self.config.T("warning"), self.config.T("cannot_delete_last"))
                return
            confirm = QMessageBox.question(self, self.config.T("confirm_delete"), self.config.T("delete_confirmation").format(name), QMessageBox.Yes | QMessageBox.No)
            if confirm == QMessageBox.Yes:
                del self.config.chat_services[name]
                if hasattr(self, 'all_possible_services') and name in self.all_possible_services:
                    del self.all_possible_services[name]
                if self.config.default_service == name:
                    # 如果删除的是默认服务，重置为列表第一个 (Fallback if default deleted)
                    self.config.default_service = list(self.config.chat_services.keys())[0]
                self.sync_and_save()

    def sync_and_save(self):
        """
        全量同步 UI 状态到 Config 并刷新所有关联的视图 (Full sync and refresh).
        """
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
        self.refresh_service_buttons() # 刷新书签栏
        self.config.save_config() # 保存到磁盘

    def on_startup_combo_context_menu(self, pos):
        """
        设置面板中启动页下拉框的右键快捷操作 (Context menu for settings startup combo).
        """
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
        """
        弹出对话框，允许用户添加新的 AI 对话网页或修改现有书签 (Dialog for adding/editing AI service bookmarks).
        :param edit_mode: 是否处于编辑模式 (Whether in edit mode)
        :param old_name: 编辑前的服务名称 (Original name before edit)
        :param old_url: 编辑前的 URL (Original URL before edit)
        """
        dialog = QDialog(self)
        dialog.setWindowTitle(self.config.T("edit_web_link") if edit_mode else self.config.T("add_web_link"))
        layout = QFormLayout(dialog)
        
        # 输入框：服务名称和 URL (Inputs for service name and URL)
        name_edit = QLineEdit(old_name)
        url_edit = QLineEdit(old_url if edit_mode else "https://")
        
        layout.addRow(self.config.T("web_name"), name_edit)
        layout.addRow(self.config.T("url"), url_edit)
        
        # 按钮栏：保存与取消 (Buttons for Save and Cancel)
        btn_box = QHBoxLayout()
        add_btn = QPushButton(self.config.T("save") if edit_mode else self.config.T("add"))
        cancel_btn = QPushButton(self.config.T("cancel"))
        btn_box.addWidget(add_btn)
        btn_box.addWidget(cancel_btn)
        layout.addRow(btn_box)
        
        def do_add():
            """
            执行添加/修改逻辑，并进行基础校验 (Perform validation and addition logic).
            """
            name = name_edit.text().strip()
            url = url_edit.text().strip()
            # 校验名称不为空且 URL 格式基本正确 (Basic validation for non-empty name and HTTP protocol)
            if name and url.startswith("http"):
                if edit_mode and old_name != name:
                    # 如果名称改变，删除旧条目 (Delete old entries if name changed during edit)
                    if old_name in self.config.chat_services: del self.config.chat_services[old_name]
                    if hasattr(self, 'all_possible_services') and old_name in self.all_possible_services: 
                        del self.all_possible_services[old_name]
                
                # 更新字典并刷新配置 (Update dictionary and refresh persistence)
                self.config.chat_services[name] = url
                if hasattr(self, 'all_possible_services'):
                    self.all_possible_services[name] = url
                
                # 设置当前新增/编辑的为默认页 (Set the new/edited one as default for focus)
                self.config.default_service = name
                self.sync_and_save()
                dialog.accept()
            else:
                QMessageBox.warning(dialog, self.config.T("input_error"), self.config.T("input_error_msg"))
        
        add_btn.clicked.connect(do_add)
        cancel_btn.clicked.connect(dialog.reject)
        dialog.exec()

    def start_enrichment_thread(self):
        """
        启动后台线程执行 GSEA/Over-representation 分析 (Launch background thread for enrichment analysis).
        该过程耗时较长，通过 QThread 避免阻塞主界面。
        """
        path = self.path_edit.text().strip()
        if not os.path.exists(path):
            QMessageBox.critical(self, self.config.T("error"), self.config.T("file_not_found"))
            return
            
        is_local = self.rb_local.isChecked()
        db_name = self.enrich_db_combo.currentText()
        
        # 确定 GMT 数据库的物理路径 (Resolve physical path of the GMT database)
        if is_local:
            db_file = self.db_path_map.get(db_name, "")
            enrich_db_path = os.path.join(self.db_dir, db_file)
            if not db_file or not os.path.exists(enrich_db_path):
                QMessageBox.critical(self, self.config.T("error"), f"Local database file not found: {db_file}\nSearched in: {self.db_dir}")
                return
        else:
            # 在线模式直接透传数据库名称 (Online mode passes db name to engine)
            enrich_db_path = db_name
            
        # UI 状态切换：禁用按钮，显示进度条 (UI state toggle: disable button, show progress)
        self.run_enrich_btn.setEnabled(False)
        self.enrich_status_container.show()
        self.status_label.setText(self.config.T("enrich_running"))
        self.enrich_progress.setValue(0)
        
        # 切换到富集结果标签页以便用户观察 (Auto-switch to enrichment tab)
        self.show_panel_tab(self.tab_enrichment_widget, 'tab_enrichment')
        self.clear_enrichment_display()
        
        # 初始化 worker 线程并连接信号槽 (Initialize worker and connect signals/slots)
        from chatcellanno.gui.workers import EnrichmentWorker # 动态导入避免循环依赖
        self.worker = EnrichmentWorker(path, self.species_edit.text(), enrich_db_path, is_local=is_local)
        self.worker.progress.connect(self.enrich_progress.setValue)
        self.worker.finished.connect(self.on_enrichment_finished)
        self.worker.error.connect(self.on_enrichment_error)
        self.worker.start()

    def on_enrichment_finished(self, data):
        """
        富集分析成功后的回调 (Callback on successful enrichment completion).
        :param data: 包含各 Cluster 富集项及其可视化路径的字典 (Enrichment result dictionary).
        """
        self.enrichment_data = data
        self.run_enrich_btn.setEnabled(True)
        self.enrich_status_container.hide() 
        self.update_enrichment_display(data) # 渲染可视化卡片
        QMessageBox.information(self, self.config.T("success"), "Enrichment complete! AI hints are ready for Step 3.")

    def on_enrichment_error(self, err_msg):
        """
        富集分析失败后的错误处理 (Error handling for enrichment failure).
        """
        self.run_enrich_btn.setEnabled(True)
        self.enrich_status_container.hide() 
        QMessageBox.critical(self, self.config.T("error"), f"Enrichment failed: {err_msg}")

    def clear_enrichment_display(self):
        """
        清空右侧面板中累积的旧富集卡片 (Clear old enrichment result cards).
        """
        for i in reversed(range(self.enrich_detail_layout.count())): 
            item = self.enrich_detail_layout.itemAt(i)
            if item.widget():
                item.widget().setParent(None)
        
        # 添加加载中占位符 (Add 'Running...' placeholder)
        placeholder = QLabel(self.config.T("enrich_running")) 
        placeholder.setAlignment(Qt.AlignCenter)
        self.enrich_detail_layout.addWidget(placeholder)

    def update_enrichment_display(self, enrichment_data):
        """
        将分析得到的 GSEA/Enrichr 数据渲染为精美的 UI 组合框 (Render enrichment results as GUI group boxes).
        包含了：1. AI 提示语 2. 详细数据表格链接 3. 彩色可视化图表
        """
        # 1. 首先彻底清理布局中的现有组件 (Thoroughly clean current layout)
        for i in reversed(range(self.enrich_detail_layout.count())): 
            item = self.enrich_detail_layout.itemAt(i)
            if item.widget():
                item.widget().setParent(None)
        
        # 2. 遍历每个 Cluster 构建结果卡片 (Iterate per cluster to build cards)
        for cluster, data in enrichment_data.items():
            c_group = QGroupBox(f"Cluster {cluster}")
            c_layout = QVBoxLayout()
            
            # 渲染翻译后的 AI 参考术语 (Render translated AI reference terms)
            c_layout.addWidget(QLabel("<b>[AI Prompt Hints]</b>"))
            hints_text = "; ".join(data['hints'])
            c_layout.addWidget(QLabel(hints_text))
            
            # 提供直接打开详细 CSV 文件的按钮 (Provide button to open raw CSV table)
            if data.get('table_path'):
                c_layout.addWidget(QLabel("<b>[Database Results]</b>"))
                btn_open = QPushButton(f"Open Full Result Table (CSV)")
                # 使用闭包捕获当前路径 (Capture current path via closure)
                def make_open_func(p):
                    return lambda: os.startfile(p)
                btn_open.clicked.connect(make_open_func(data['table_path']))
                c_layout.addWidget(btn_open)
            
            # 如果存在可视化图表，则在 UI 中嵌入缩略图 (Embed visualization plot if available)
            if data.get('plot_path') and os.path.exists(data['plot_path']):
                c_layout.addWidget(QLabel("<b>[Enrichment Visualization]</b>"))
                img_label = QLabel()
                pixmap = QPixmap(data['plot_path'])
                # 缩放图片以适应侧边栏宽度，保持平滑抗锯齿 (Scale image with smooth transformation)
                img_label.setPixmap(pixmap.scaledToWidth(500, Qt.SmoothTransformation))
                c_layout.addWidget(img_label)

            c_group.setLayout(c_layout)
            self.enrich_detail_layout.addWidget(c_group)
        
        # 添加底部弹簧以便卡片向上对齐 (Add stretch to align cards to the top)
        self.enrich_detail_layout.addStretch()

    def browse_file(self):
        """
        调用文件选择对话框获取 Marker 基因文件 (File picker for marker genes).
        """
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Marker File", "", "CSV/TSV Files (*.csv *.tsv *.txt)")
        if file_path:
            self.path_edit.setText(file_path)

    def browse_matrix_file(self):
        """
        调用文件选择对话框获取全量表达矩阵文件 (File picker for expression matrix).
        用于验证 marker 基因在不同 cluster 中的特异性表达情况。
        """
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Expression Matrix File", "", "CSV/TSV Files (*.csv *.tsv *.txt)")
        if file_path:
            self.matrix_path_edit.setText(file_path)

    def dragEnterEvent(self, event):
        """
        实现拖拽进入回调，检查是否有合法的 URL（文件路径）被拖入 (Handle drag enter).
        """
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        """
        处理文件释放动作，实现“拖入即加载”的智能识别逻辑 (Handle drop file with smart auto-detection).
        可以自动区分文件夹（GMT数据库）、Marker文件和表达量矩阵。
        """
        files = [u.toLocalFile() for u in event.mimeData().urls()]
        if not files:
            return
            
        path = files[0] # 获取拖入的第一个文件路径 (Get first file path)
        
        if os.path.isdir(path):
            # 如果拖入的是目录，尝试将其作为本地数据库加载 (Load folder as local database)
            self.load_databases_from_folder(path)
        else:
            # 智能检测：根据内容初步推断是 Matrix 还是 Marker (Content-based smart detection)
            is_matrix = False
            try:
                # 读取前几行进行特征分析 (Analyze first few lines)
                with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                    header = f.readline().lower()
                    first_line = f.readline()
                
                # 逻辑 1：普通表达矩阵通常包含大量基因列，或者包含 'Cluster_' 标识 (Matrix markers)
                if 'cluster_' in first_line.lower() or 'cluster_' in header:
                    is_matrix = True
                # 逻辑 2：如果列数极多（>15）且不含常见的差异分析字段 (Check column count vs marker fields)
                elif ',' in header and len(header.split(',')) > 15:
                    if not any(x in header for x in ['p_val', 'logfold', 'adj']):
                        is_matrix = True
            except:
                pass # 解析失败时采用保守策略

            if is_matrix:
                # 自动填入表达矩阵输入框 (Auto-fill matrix input)
                self.matrix_path_edit.setText(path)
            else:
                # 默认填入 Marker 文件输入框 (Default to marker input)
                self.path_edit.setText(path)

    def navigate_to_url(self):
        """
        让内置浏览器跳转到地址栏指定的 URL (Trigger browser navigation).
        支持自动补全 https:// 协议头 (Auto-prepend https).
        """
        url = self.url_bar.text()
        if url and not url.startswith("http"):
            url = "https://" + url
        self.browser.setUrl(QUrl(url))

    def browse_visual_image(self):
        """从本地文件浏览器选择图片文件"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Image", "", "Images (*.png *.xpm *.jpg *.jpeg *.bmp *.tif *.tiff *.gif *.webp)"
        )
        if file_path:
            self.img_path_edit.setText(file_path)
            pixmap = QPixmap(file_path)
            if not pixmap.isNull():
                self.set_preview_image(pixmap)

    def img_preview_drag_enter(self, event):
        """处理图片预览框的拖入事件"""
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                if url.toLocalFile().lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff', '.gif', '.webp')):
                    event.acceptProposedAction()
                    return

    def img_preview_drop(self, event):
        """处理图片预览框的放下事件"""
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            self.img_path_edit.setText(file_path)
            pixmap = QPixmap(file_path)
            if not pixmap.isNull():
                self.set_preview_image(pixmap)
                break

    def set_preview_image(self, pixmap):
        """统一设置预览图和状态"""
        self.img_data = pixmap
        # 使用固定的最大预览高度 130 确保比例协调
        scaled_pixmap = self.img_data.scaledToHeight(130, Qt.SmoothTransformation)
        self.img_preview.setPixmap(scaled_pixmap)
        self.img_preview.setText("") # 清除占位文字

    def generate_prompt(self):
        """
        核心逻辑：根据用户的所有输入参数，工厂化生产 AI 提示词 (Orchestrate prompt generation factory).
        整合了：1. Marker 基因 2. 富集分析结果 3. 视觉上下文 4. 全量矩阵校验
        """
        path = self.path_edit.text().strip()
        if not os.path.exists(path):
            QMessageBox.critical(self, self.config.T("error"), self.config.T("file_not_found"))
            return
        
        try:
            mode = "concise" if self.rb_concise.isChecked() else "detailed"
            
            # 提取富集分析产生的 AI 关键词提示 (Extract high-level functional hints)
            enrichment_hints = None
            if hasattr(self, 'use_enrichment_cb') and self.use_enrichment_cb.isChecked() and self.enrichment_data:
                enrichment_hints = {k: v['hints'] for k, v in self.enrichment_data.items()}

            # 提取视觉背景类型（UMAP/TSNE 等）(Visual modality context)
            visual_context_str = None
            if self.img_data:
                visual_context_str = "UMAP/t-SNE Plot"
                if hasattr(self, "img_path_edit") and self.img_path_edit.text():
                    visual_context_str += f" ({os.path.basename(self.img_path_edit.text())})"
                
            # 调用 backend core 逻辑生成多模态 Prompt (Call core logic for multimodal fusion)
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
                enrichment_hints=enrichment_hints,
                visual_context=visual_context_str,
                check_expression_file=self.matrix_path_edit.text().strip()
            )
            
            # 在预览窗口展示结果 (Show finalized text in UI)
            self.prompt_display.setText(prompt)
            
            # 自动切换到 AI 浏览器标签页 (Auto-switch to interaction tab)
            self.show_panel_tab(self.tab_browser_widget, 'tab_browser')
            
            # --- 多模态导出策略：文件暂存与拖拽支持 (Export strategy with file temp storage) ---
            # 建立软件执行目录下的 Temp_Prompts 文件夹，方便用户在浏览器中使用文件上传模式
            export_dir = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "Temp_Prompts")
            # 每次生成前清理陈旧数据 (Clean old files before new export)
            if os.path.exists(export_dir):
                import shutil
                try: shutil.rmtree(export_dir) 
                except: pass
            os.makedirs(export_dir, exist_ok=True)
            
            file_paths = []
            
            # 1. 保存提示词文本文件 (Save prompt text file)
            text_path = os.path.join(export_dir, "1_Annotation_Prompt.txt")
            with open(text_path, 'w', encoding='utf-8') as f:
                f.write(prompt)
            file_paths.append(text_path)
            
            # 2. 如果存在视觉上下文，保存图片文件 (Save visual context image if existing)
            if self.img_data:
                img_path = os.path.join(export_dir, "2_Visual_Context.png")
                self.img_data.toImage().save(img_path, "PNG")
                file_paths.append(img_path)
            
            # 3. 配置拖拽组件：将文件路径绑定到“绿色拖拽区” (Configure drag-out area with file paths)
            self.drag_area_info.setFiles(file_paths)
            self.drag_area_info.setText(f"📁 已就绪 ({len(file_paths)}个文件)\n在此处按住鼠标左键\n并拖拽进浏览器对话框")
            self.drag_area_group.show() # 显示拖拽手柄 (Show drag handle)
            
            msg = f"已生成拉拽文件！\n\n文件保存于: {export_dir}\n\n使用方式：按住下方绿色区域并拖动至AI聊天界面。"
            
            # 弹出成功提示 (Show success notification)
            QMessageBox.information(self, self.config.T("success"), msg)
        except Exception as e:
            QMessageBox.critical(self, self.config.T("error"), f"Failed: {str(e)}")

    def show_api_settings_dialog(self):
        """
        弹出 API 配置对话框，允许用户输入 Key 和自定义 Endpoint (API configuration dialog).
        支持 SiliconFlow、OpenAI、DeepSeek 等兼容 OpenAI 格式的后端。
        """
        dialog = QDialog(self)
        dialog.setWindowTitle(self.config.T("api_settings"))
        layout = QFormLayout(dialog)
        
        # 构建输入控件 (Build input widgets)
        api_key_edit = QLineEdit(self.config.api_settings.get("api_key", ""))
        api_key_edit.setEchoMode(QLineEdit.Password) # 隐藏 API Key (Password echo mode)
        api_base_edit = QLineEdit(self.config.api_settings.get("api_base", "https://api.siliconflow.cn/v1"))
        api_model_edit = QLineEdit(self.config.api_settings.get("api_model", "Qwen/Qwen2.5-72B-Instruct"))
        api_model_edit.setToolTip("e.g. Qwen/Qwen2.5-72B-Instruct (Text) or Pro/Qwen/Qwen2-VL-72B-Instruct (Vision/Image)")
        
        layout.addRow(self.config.T("api_key"), api_key_edit)
        layout.addRow(self.config.T("api_base"), api_base_edit)
        layout.addRow(self.config.T("api_model"), api_model_edit)
        
        # 保存按钮逻辑 (Save button logic)
        btn_box = QHBoxLayout()
        save_btn = QPushButton(self.config.T("save"))
        save_btn.clicked.connect(lambda: dialog.accept())
        btn_box.addWidget(save_btn)
        
        layout.addRow(btn_box)
        
        if dialog.exec() == QDialog.Accepted:
            # 更新持久化配置 (Update and persist configuration)
            self.config.api_settings["api_key"] = api_key_edit.text().strip()
            self.config.api_settings["api_base"] = api_base_edit.text().strip()
            self.config.api_settings["api_model"] = api_model_edit.text().strip()
            self.config.save_config()

    def call_api_directly(self):
        """
        核心 API 调用逻辑：跳过浏览器手动操作，直接通过 Python 通信获取 AI 结果 (Direct API calling bypass).
        支持多模态请求：将生成的 Prompt 和 Base64 编码的图片一并发送。
        """
        # 1. 检查 API Key 是否已设置 (Pre-flight check for API key)
        if not self.config.api_settings.get("api_key"):
            QMessageBox.warning(self, self.config.T("warning"), "Please set your API Key first.")
            self.show_api_settings_dialog()
            if not self.config.api_settings.get("api_key"):
                return
                
        # 2. 首先在后台静默生成 Prompt (Quietly generate prompt backend)
        path = self.path_edit.text().strip()
        if not os.path.exists(path):
            QMessageBox.critical(self, self.config.T("error"), self.config.T("file_not_found"))
            return
            
        try:
            mode = "concise" if self.rb_concise.isChecked() else "detailed"
            enrichment_hints = {k: v['hints'] for k, v in self.enrichment_data.items()} if self.enrichment_data else None
            visual_context_str = "UMAP/t-SNE Plot" if self.img_data else None
                
            from chatcellanno.core import annotate_cell_types
            prompt, _ = annotate_cell_types(
                marker_file=path, step="generate", species=self.species_edit.text(),
                tissue=self.tissue_edit.text(), top_n=int(self.top_n_edit.text()),
                mode=mode, exclude_types=self.exclude_edit.text(), use_enrichment=False, 
                enrichment_hints=enrichment_hints, visual_context=visual_context_str,
                check_expression_file=self.matrix_path_edit.text().strip()
            )
            
            # 3. 如果存在视觉上下文，转换为 Base64 编码 (Convert images to base64 for multimodal API)
            base64_img = None
            if self.img_data:
                from PySide6.QtCore import QByteArray, QBuffer, QIODevice
                import base64
                ba = QByteArray()
                buffer = QBuffer(ba)
                buffer.open(QIODevice.WriteOnly)
                self.img_data.toImage().save(buffer, "PNG")
                base64_img = base64.b64encode(ba.data()).decode('utf-8')
            
            # 4. 初始化 API Worker 后台线程 (Initialize API worker thread)
            from chatcellanno.gui.workers import ApiWorker 
            self.api_worker = ApiWorker(
                prompt=prompt,
                api_key=self.config.api_settings["api_key"],
                base_url=self.config.api_settings["api_base"],
                model=self.config.api_settings["api_model"],
                base64_image=base64_img
            )
            self.api_worker.finished.connect(self.on_api_success)
            self.api_worker.error.connect(self.on_api_error)
            
            # UI 交互：禁用按钮并展示状态提示 (Disable button and update status bar)
            self.btn_call_api.setEnabled(False)
            self.btn_call_api.setText("Loading...")
            
            self.statusBar().showMessage(self.config.T("api_calling"))
            self.api_worker.start() # 启动线程 (Launch)
            
        except Exception as e:
            QMessageBox.critical(self, self.config.T("error"), f"Failed: {str(e)}")
            
    def on_api_success(self, response_text):
        """
        API 调用成功后的 UI 回调 (Callback on API response success).
        """
        self.btn_call_api.setEnabled(True)
        self.btn_call_api.setText("🚀 " + self.config.T("use_api"))
        self.statusBar().showMessage(self.config.T("api_success"), 5000)
        
        # 将结果填充到解析区输入框 (Populate result into response area)
        self.response_input.setText(response_text)
        
        # 自动切换到分析结果页 (Auto-switch to step 4)
        self.show_panel_tab(self.tab_browser_widget, 'tab_browser')
        
        QMessageBox.information(self, self.config.T("success"), "API call successful! Now click 'Process AI Output'.")
        
    def on_api_error(self, error_msg):
        """
        API 调用失败时的错误提示 (Handle API error response).
        """
        self.btn_call_api.setEnabled(True)
        self.btn_call_api.setText("🚀 " + self.config.T("use_api"))
        self.statusBar().showMessage("API Error", 5000)
        QMessageBox.critical(self, self.config.T("error"), error_msg)

    def process_response(self):
        """
        执行解析流程：将 AI 返回的文本注释重新映射到单细胞对象 (Execute LLM parsing pipeline).
        该过程会识别格式、映射 Cluster 并生成对应的 Python/R 处理代码。
        """
        path = self.path_edit.text().strip()
        response = self.response_input.toPlainText().strip()
        
        # 基础校验 (Pre-parsing check)
        if not response:
            QMessageBox.warning(self, self.config.T("warning"), self.config.T("paste_warning"))
            return

        if not path or not os.path.exists(path):
            # 允许无本地文件时的预览解析 (Allow parsing without source file for preview)
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
            
            # 调用 backend parser 逻辑 (Call core logic for parsing)
            from chatcellanno.core import annotate_cell_types

            # result_df: 数据帧格式的结果映射 (Mapping DF)
            # code_snippet: 生成的自动重命名代码块 (Generated rename code)
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
        """
        将 Step 5 生成的代码片段快速复制到剪贴板 (Copy generated code to clipboard).
        用户可以直接在 Jupyter Lab 或 RStudio 中粘贴执行以应用注释。
        """
        pyperclip.copy(self.code_display.toPlainText())
        QMessageBox.information(self, "Copied", "Code copied to clipboard.")

    def export_analysis_report(self):
        """
        全量报告生成功能：整合标注结果、富集分析和视觉背景 (Generate comprehensive analysis report).
        采用 Markdown 格式导出，支持图文混排，便于整理科学论文或实验记录。
        """
        try:
            # 1. 采集关键指标数据 (Gather key data streams)
            results_text = self.response_input.toPlainText().strip()
            if not results_text:
                QMessageBox.warning(self, self.config.T("warning"), "No annotation results found. Please run Step 3/4 first.")
                return

            # 构建富集结果汇总信息 (Construct enrichment summary info)
            enrichment_info = ""
            if hasattr(self, "enrichment_data") and self.enrichment_data:
                # 选取每个 Cluster 的前 5 个显著项进行罗列 (Top 5 significant terms per cluster)
                enrichment_info += "### 🧬 Pathway Enrichment (GSEApy)\n"
                for cluster, df in self.enrichment_data.items():
                    top_terms = df.sort_values("Adjusted P-value").head(5)["Term"].tolist()
                    enrichment_info += f"- **Cluster {cluster}**: " + ", ".join(top_terms) + "\n"
            else:
                # 如果未进行富集分析，展示占位符 (Placeholder if enrichment skipped)
                enrichment_info = "*(Pathway enrichment data not available)*"

            # 2. 应用报告模板填充内容 (Format with Markdown template)
            from datetime import datetime
            template = self.config.get("report_template", "")
            if not template:
                # 默认后备模板：包含标题、日期、注释结果和生物学背景 (Default fallback template)
                template = "# Single-Cell Analysis Report\nDate: {date}\n\n## 🏷 Annotation Results\n{results}\n\n## 🧬 Biological Context\n{enrichment}"

            report_content = template.format(
                date=datetime.now().strftime("%Y-%m-%d %H:%M"),
                results=results_text,
                enrichment=enrichment_info
            )

            # 3. 弹出保存对话框 (Standard file save dialog)
            from PySide6.QtWidgets import QFileDialog
            file_path, _ = QFileDialog.getSaveFileName(
                self, 
                self.config.T("gen_report"), 
                f"Cell_Analysis_Report_{datetime.now().strftime('%m%d_%H%M')}.md",
                "Markdown (*.md);;Text (*.txt)"
            )

            if file_path:
                # 写入本地文件 (Write to local disk)
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(report_content)
                QMessageBox.information(self, "Success", f"Report saved to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, self.config.T("error"), f"Report failed: {str(e)}")

    def search_genome(self):
        """
        执行基因组跨物种搜索与多序列比对预览 (Perform cross-species genome search & MSA).
        该功能允许用户根据基因名称在本地序列库中检索并以 MEGA 配色方案预览序列差异。
        """
        gene_name = self.genome_search_edit.text().strip()
        if not gene_name:
            return
            
        # 调用底层序列检索引擎 (Call sequence search engine)
        from chatcellanno.genome_utils import find_gene_sequences
        results = find_gene_sequences(gene_name, getattr(self, "genome_dir", ""))
        
        # 排序策略：人类序列排首位作为参考，随后是小鼠及其他 (Sort: Human -> Mouse -> Others)
        def sort_key(r):
            name_l = r['species'].lower()
            if "human" in name_l: return (0, r['species'])
            if "mouse" in name_l: return (1, r['species'])
            return (2, r['species'])
        results.sort(key=sort_key)
        
        # UI 清空逻辑：彻底删除旧的结果卡片和布局弹簧 (Thorough UI clean)
        while self.genome_results_layout.count():
            item = self.genome_results_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
                
        if not results:
            # 未找到结果时的视觉反馈 (No results found feedback)
            label = QLabel(self.config.T("no_results"))
            label.setAlignment(Qt.AlignCenter)
            self.genome_results_layout.addWidget(label)
            return
            
        # 1. 进行动态序列比对与格式化 (Perform relative alignment & HTML formatting)
        from chatcellanno.genome_utils import align_records, format_alignment_html, save_fasta_combined, get_fasta_string
        
        # 对选中的序列进行对齐操作 (MSA processing)
        aligned_results = align_records(results)
        
        # 2. 构建渲染结果容器 (Build rendering container)
        alignment_group = QGroupBox(f"Sequence Comparison for '{gene_name}'")
        alignment_layout = QVBoxLayout()
        
        # 生成基于 HTML 的高级着色视图 (High-fidelity HTML color view)
        html_content = format_alignment_html(aligned_results)
        view_advanced = QTextEdit()
        view_advanced.setReadOnly(True)
        view_advanced.setHtml(html_content)
        view_advanced.setLineWrapMode(QTextEdit.NoWrap) # 允许横向滚动查看长序列
        
        alignment_layout.addWidget(view_advanced)
        
        # 3. 底部操作按钮：复制与导出 (Bottom actions: Copy & Export)
        btn_layout = QHBoxLayout()
        
        copy_all_btn = QPushButton("Copy Combined FASTA")
        fasta_text = get_fasta_string(results)
        copy_all_btn.clicked.connect(lambda: pyperclip.copy(fasta_text))
        
        download_btn = QPushButton("Export Combined FASTA")
        def do_download():
            save_path, _ = QFileDialog.getSaveFileName(self, "Export Combined FASTA", f"{gene_name}_homologs.fa", "FASTA Files (*.fa *.fasta)")
            if save_path:
                try:
                    save_fasta_combined(results, save_path)
                    QMessageBox.information(self, "Success", f"Saved to {save_path}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to save: {e}")
        
        download_btn.clicked.connect(do_download)
        
        btn_layout.addWidget(copy_all_btn)
        btn_layout.addWidget(download_btn)
        btn_layout.addStretch() # 向左对齐按钮 (Left align buttons)
        
        alignment_layout.addLayout(btn_layout)
        alignment_group.setLayout(alignment_layout)
        self.genome_results_layout.addWidget(alignment_group)
        
        # 添加底部弹簧以便卡片向上对齐 (Add bottom stretch)
        self.genome_results_layout.addStretch()

    def show_help(self):
        """
        提供完整的软件使用指南，涵盖从数据加载到结果导出的全流程 (Comprehensive user guide dialog).
        """
        msg = """
        ChatCellAnno (单细胞智能注释助手) 使用指南:
        
        步骤 1 (数据准备): 
        - 加载您的 Marker 基因文件 (CSV/TSV/TXT)，支持 Scanpy 或 Seurat 格式。
        
        步骤 2 (富集分析): 
        - [可选] 运行功能富集。您可以关联本地数据库文件夹（直接拖入目录）。
        - 这一步产生的生物学关键词能有效防止 AI 标注时的“幻觉”现象。
        
        步骤 3 (生成提示词): 
        - 设置物种、组织环境等参数，一键生成多模态 Prompt。提示词会自动存入剪贴板。
        
        步骤 4 (AI 交互):
        - 利用内置浏览器直接在 ChatGPT/DeepSeek 中提问。
        - 粘贴刚刚生成的 Prompt，并复制 AI 返回的 MARKDOWN 表格。
        - 之后点击“处理 AI 输出”即可完成数据回传。
        
        步骤 5 (同步数据):
        - 复制生成的 R/Python 代码片段，在您的分析脚本中运行即可应用翻译结果。
        """
        QMessageBox.information(self, self.config.T("guide"), msg.strip())
