
import os
config_path = r'd:\System\Documents\GitHub\ChatCellAnno\chatcellanno\config.py'
with open(config_path, 'r', encoding='utf-8') as f:
    text = f.read()

# English section
if '"query": "Query"' in text:
    text = text.replace('"query": "Query"', 
                       '"query": "Query",\n        "marker_list_label": "Marker Gene List (Required):",\n        "enrich_assist_label": "Perform enrichment to assist AI annotation (Optional)",\n        "panels_control": "Show Tab Panels"')

# Chinese section
if '"copy_img_success": "图像已复制！现在您可以将它粘贴到 AI 的对话框中了。"' in text:
    text = text.replace('"copy_img_success": "图像已复制！现在您可以将它粘贴到 AI 的对话框中了。"', 
                       '"copy_img_success": "图像已复制！现在您可以将它粘贴到 AI 的对话框中了。",\n        "marker_list_label": "差异基因列表 (必选):",\n        "enrich_assist_label": "进行功能富集以辅助 AI 注释 (可选)",\n        "panels_control": "显示面板 (Panels)"')

with open(config_path, 'w', encoding='utf-8') as f:
    f.write(text)
print("Finished i18n update")
