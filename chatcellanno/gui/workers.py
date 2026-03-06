import json
import urllib.request
from typing import Optional
from PySide6.QtCore import QThread, Signal
from chatcellanno.extractor import extract_markers_from_file
from chatcellanno.enrichment import perform_enrichment

class ApiWorker(QThread):
    finished = Signal(str)
    error = Signal(str)

    def __init__(self, prompt: str, api_key: str, base_url: str, model: str, base64_image: Optional[str] = None):
        super().__init__()
        self.prompt = prompt
        self.api_key = api_key
        # Ensure url ends with /chat/completions for OpenAI compat
        if not base_url.endswith("/chat/completions"):
            if base_url.endswith("/"):
                self.base_url = base_url + "chat/completions"
            else:
                self.base_url = base_url + "/chat/completions"
        else:
            self.base_url = base_url
            
        self.model = model
        self.base64_image = base64_image

    def run(self):
        try:
            # Build OpenAI-Compatible Payload
            messages = []
            
            if self.base64_image:
                # Multimodal request
                user_content = [
                    {"type": "text", "text": self.prompt},
                    {"type": "image_url", "image_url": {"url": f"data:image/png;base64,{self.base64_image}"}}
                ]
                messages.append({"role": "user", "content": user_content})
            else:
                # Text only request
                messages.append({"role": "user", "content": self.prompt})

            data = {
                "model": self.model,
                "messages": messages,
                "temperature": 0.1,  # Low temp for table reasoning
                "max_tokens": 4096,
                "stream": False
            }
            
            req = urllib.request.Request(self.base_url, data=json.dumps(data).encode('utf-8'))
            req.add_header('Content-Type', 'application/json')
            req.add_header('Authorization', f'Bearer {self.api_key}')
            
            try:
                with urllib.request.urlopen(req, timeout=120) as response:
                    result = json.loads(response.read().decode('utf-8'))
                    
                    # Extract the message
                    reply = result['choices'][0]['message']['content']
                    self.finished.emit(reply)
            except urllib.error.HTTPError as e:
                err_body = e.read().decode('utf-8')
                self.error.emit(f"HTTP Error {e.code}: {err_body}")
        
        except Exception as e:
            self.error.emit(f"API Request Failed: {str(e)}")

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
