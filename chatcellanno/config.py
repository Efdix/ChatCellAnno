import json
import os
import sys

TRANSLATIONS = {
    "en": {
        "title": "ChatCellAnno",
        "step1": "Step 1: Data Source",
        "placeholder_path": "Drag & drop file here or browse...",
        "browse": "Browse",
        "scanpy": "Scanpy (Python)",
        "seurat": "Seurat (R)",
        "step2": "Step 2: Enrichment Analysis",
        "species": "Species:",
        "tissue": "Tissue:",
        "top_n": "Top Markers:",
        "exclude": "Exclude Types (optional):",
        "prompt_mode": "Prompt Mode:",
        "concise": "Concise",
        "detailed": "Detailed",
        "step3": "Step 3: Prompt Configuration",
        "gen_btn": "Generate & Copy Prompt",
        "gen_btn_enrich": "Run Enrichment Analysis",
        "step4": "Step 4: AI Response",
        "placeholder_resp": "Paste AI response here...",
        "process_btn": "Process AI Output",
        "step5": "Step 5: Export Code",
        "code_label": "Generated Code:",
        "copy_code": "Copy Code",
        "default_page": "Default Page:",
        "browser_label": "Browser:",
        "go": "Go",
        "refresh": "Refresh",
        "enrich_running": "Enrichment Analysis in Progress...",
        "clean": "Clean",
        "enrich_complete": "Analysis Complete!",
        "local": "Local Database",
        "online": "Online (Enrichr)",
        "enrich_source": "Enrichment Source:",
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
        "lang_en": "English",
        "tab_browser": "AI Browser",
        "tab_enrichment": "Enrichment Results",
        "local_db_hint": "Please click the folder button to select database directory"
    },
    "zh": {
        "title": "ChatCellAnno - 单细胞注释助手",
        "step1": "步骤 1: 数据源",
        "placeholder_path": "拖放文件或点击浏览...",
        "browse": "浏览",
        "scanpy": "Scanpy (Python)",
        "seurat": "Seurat (R)",
        "step2": "步骤 2: 功能富集分析",
        "species": "物种:",
        "tissue": "组织:",
        "top_n": "Top Markers:",
        "exclude": "排除类型 (可选):",
        "prompt_mode": "提示词模式:",
        "concise": "简洁",
        "detailed": "详细",
        "step3": "步骤 3: 提示词配置",
        "gen_btn": "生成并复制提示词",
        "gen_btn_enrich": "运行富集分析",
        "step4": "步骤 4: AI 回答结果",
        "placeholder_resp": "在此粘贴 AI 回复...",
        "process_btn": "处理 AI 输出",
        "step5": "步骤 5: 导出代码",
        "code_label": "生成的代码:",
        "copy_code": "复制代码",
        "default_page": "默认主页:",
        "browser_label": "浏览器:",
        "go": "前往",
        "refresh": "刷新",
        "enrich_running": "正在进行富集分析...",
        "clean": "清空",
        "enrich_complete": "分析完成！",
        "local": "本地数据库",
        "online": "在线 (Enrichr)",
        "enrich_source": "富集数据源:",
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
        "lang_en": "English",
        "tab_browser": "AI 浏览器",
        "tab_enrichment": "功能富集结果",
        "local_db_hint": "请点击右侧文件夹按钮选择本地数据库目录"
    }
}

class ConfigManager:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ConfigManager, cls).__new__(cls)
            cls._instance.initialize()
        return cls._instance

    def initialize(self):
        # Determine paths
        # Use user's home directory or current directory for config to ensure write permissions
        # and persistence across EXE updates
        user_config_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.config_paths = {
            "config": os.path.join(user_config_dir, "config.json"),
            "web_storage": os.path.join(user_config_dir, "chatcellanno_web_data")
        }
        
        self.language = "en"
        self.chat_services = {
            "DeepSeek": "https://chat.deepseek.com",
            "ChatGPT": "https://chatgpt.com", 
            "Claude": "https://claude.ai",
        }
        self.default_service = "DeepSeek"
        self.session_params = {
            "species": "Human",
            "tissue": "PBMC",
            "top_n": "10",
            "exclude": "",
            "prompt_mode": "concise"
        }
        
        self.load_config()

    def load_config(self):
        if os.path.exists(self.config_paths["config"]):
            try:
                with open(self.config_paths["config"], "r", encoding="utf-8") as f:
                    config = json.load(f)
                    self.default_service = config.get("default_service", "DeepSeek")
                    loaded_params = config.get("session_params", {})
                    self.session_params.update(loaded_params)
                    custom_services = config.get("chat_services")
                    self.language = config.get("language", "en") 
                    if custom_services:
                        self.chat_services = custom_services
            except:
                pass

    def save_config(self):
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
            
    def T(self, key):
        return TRANSLATIONS.get(self.language, TRANSLATIONS["en"]).get(key, key)
