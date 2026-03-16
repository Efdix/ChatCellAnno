import os
import sys
import importlib.util

class PluginManager:
    def __init__(self, main_window):
        self.main_window = main_window
        self.plugins = []
        
    def load_plugins(self):
        # 1. Internal/Built-in plugins directory
        plugin_dirs = [os.path.dirname(__file__)]
        
        # 2. External plugins directory (next to executable/script)
        exe_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        ext_plugin_dir = os.path.join(exe_dir, "plugins")
        if os.path.exists(ext_plugin_dir) and ext_plugin_dir not in plugin_dirs:
            plugin_dirs.append(ext_plugin_dir)
        else:
            try:
                os.makedirs(ext_plugin_dir, exist_ok=True)
                plugin_dirs.append(ext_plugin_dir)
            except:
                pass
                
        for plugin_dir in plugin_dirs:
            if not os.path.exists(plugin_dir):
                continue
            for filename in os.listdir(plugin_dir):
                if filename.endswith(".py") and not filename.startswith("__"):
                    module_name = f"chatcellanno.plugins_dynamic.{filename[:-3]}"
                    spec = importlib.util.spec_from_file_location(module_name, os.path.join(plugin_dir, filename))
                    if spec and spec.loader:
                        try:
                            module = importlib.util.module_from_spec(spec)
                            spec.loader.exec_module(module)
                            if hasattr(module, 'register_plugin'):
                                plugin = module.register_plugin(self.main_window)
                                self.plugins.append(plugin)
                        except Exception as e:
                            print(f"Error loading plugin {filename}: {e}")
                        
    def init_plugins(self):
        for plugin in self.plugins:
            if hasattr(plugin, 'init_ui'):
                plugin.init_ui()
