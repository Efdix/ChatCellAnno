from PySide6.QtCore import QThread, Signal
from chatcellanno.extractor import extract_markers_from_file
from chatcellanno.enrichment import perform_enrichment

class EnrichmentWorker(QThread):
    finished = Signal(dict)
    error = Signal(str)
    progress = Signal(int)

    def __init__(self, marker_file, species, database_path, is_local=True):
        super().__init__()
        self.marker_file = marker_file
        self.species = species
        self.database_path = database_path
        self.is_local = is_local

    def run(self):
        try:
            self.progress.emit(10)
            # 1. Extract markers for enrichment (top 100)
            markers_for_enrich = extract_markers_from_file(self.marker_file, top_n=100)
            
            self.progress.emit(30)
            # 2. Prepare enrichment input
            enrich_input = {k: v.split(", ") for k,v in markers_for_enrich.items()}
            
            # 3. Perform enrichment
            results = perform_enrichment(
                enrich_input, 
                species=self.species, 
                database_path=self.database_path,
                top_term_n=3,
                is_local=self.is_local
            )
            self.progress.emit(100)
            self.finished.emit(results)
        except Exception as e:
            self.error.emit(str(e))
