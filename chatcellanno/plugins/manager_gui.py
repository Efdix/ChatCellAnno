import os
import shutil
import sys
from PySide6.QtWidgets import QDialog, QVBoxLayout, QLabel, QListWidget, QMessageBox
from PySide6.QtCore import Qt

class PluginManagerDialog(QDialog):
    def __init__(self, parent, plugin_manager):
        super().__init__(parent)
        self.plugin_manager = plugin_manager
        self.setWindowTitle("🧩 插件管理 (Plugin Manager)")
        self.resize(400, 300)
        self.setAcceptDrops(True)
        
        layout = QVBoxLayout(self)
        
        # Instruction label
        info_label = QLabel("⬇️ 将 *.py 插件文件拖拽到此窗口即可完成安装\n\n已安装的插件：")
        info_label.setAlignment(Qt.AlignCenter)
        info_label.setStyleSheet("padding: 10px; color: #555; border: 2px dashed #aaa;")
        layout.addWidget(info_label)
        
        # List of plugins
        self.list_widget = QListWidget()
        layout.addWidget(self.list_widget)
        
        self.refresh_list()

    def refresh_list(self):
        self.list_widget.clear()
        
        # Get external plugin dir
        exe_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.ext_plugin_dir = os.path.join(exe_dir, "plugins")
        
        if not os.path.exists(self.ext_plugin_dir):
            try:
                os.makedirs(self.ext_plugin_dir, exist_ok=True)
            except Exception as e:
                print(f"Error creating plugin directory: {e}")
                
        # List all plugin files in the directory
        if os.path.exists(self.ext_plugin_dir):
            for filename in os.listdir(self.ext_plugin_dir):
                if filename.endswith(".py") and not filename.startswith("__"):
                    self.list_widget.addItem(f"📄 {filename}")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        urls = event.mimeData().urls()
        installed_count = 0
        
        for url in urls:
            local_path = url.toLocalFile()
            if local_path.endswith(".py"):
                filename = os.path.basename(local_path)
                dest_path = os.path.join(self.ext_plugin_dir, filename)
                try:
                    shutil.copy2(local_path, dest_path)
                    installed_count += 1
                except Exception as e:
                    QMessageBox.critical(self, "错误", f"拷贝插件 {filename} 失败:\n{str(e)}")
            else:
                QMessageBox.warning(self, "警告", "只支持挂载以 .py 结尾的 Python 插件脚本文件。")
                
        if installed_count > 0:
            self.refresh_list()
            QMessageBox.information(self, "安装成功", f"成功安装 {installed_count} 个插件！\n提示：请重启软件以加载新插件。")
