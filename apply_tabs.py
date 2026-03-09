"
import sys

with open(r'd:\System\Documents\GitHub\ChatCellAnno\chatcellanno\gui\main_window.py', 'r', encoding='utf-8') as f:
    text = f.read()

text = text.replace('        self.right_tabs.addTab(enrichment_container, self.config.T(\
tab_enrichment\))',
                    '        self.tab_enrichment_widget = enrichment_container\n        self.right_tabs.addTab(enrichment_container, self.config.T(\tab_enrichment\))')

text = text.replace('        self.right_tabs.addTab(browser_container, self.config.T(\tab_browser\))',
                    '        self.tab_browser_widget = browser_container\n        self.right_tabs.addTab(browser_container, self.config.T(\tab_browser\))')

text = text.replace('        self.right_tabs.addTab(genome_container, self.config.T(\tab_genome\))',
                    '        self.tab_genome_widget = genome_container\n        self.right_tabs.addTab(genome_container, self.config.T(\tab_genome\))')

with open(r'd:\System\Documents\GitHub\ChatCellAnno\chatcellanno\gui\main_window.py', 'w', encoding='utf-8') as f:
    f.write(text)
"
