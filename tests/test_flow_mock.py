import sys
import os
import unittest
from unittest.mock import MagicMock, patch

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import AFTER path update
import chatcellanno.core

class TestEnrichmentFlow(unittest.TestCase):
    
    def test_enrichment_integration(self):
        """
        Test that IF enrichment is enabled and works, does the prompt get updated correctly?
        We mock the actual enrichment calculation to bypass broken gseapy.
        """
        
        # 1. Mock the perform_enrichment function
        mock_enrichment_results = {
            "0": ["GO:0001 (T cell activation) P=1e-5", "GO:0002 (Immune response) P=1e-4"],
            "1": ["GO:0003 (B cell activation) P=1e-6"]
        }
        
        mock_perform = MagicMock(return_value=mock_enrichment_results)
        
        # 2. Patch core module
        # We need to ensure HAS_ENRICHMENT is True and perform_enrichment is our mock
        with patch('chatcellanno.core.HAS_ENRICHMENT', True), \
             patch('chatcellanno.core.perform_enrichment', mock_perform):
            
            print("\n[Test] Running annotate_cell_types with mocked enrichment...")
            
            # 3. Define path to example data
            example_csv = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'examples', 'example.csv'))
            
            # 4. Run the Core Function
            prompt, enrichment_hints = chatcellanno.core.annotate_cell_types(
                step="generate",
                marker_file=example_csv,
                species="Human",
                tissue="PBMC",
                top_n=5,
                use_enrichment=True, # Enable enrichment!
                enrichment_db="GO_Biological_Process_2021"
            )
            
            # 5. Verify the Mock was called
            # We don't verify arguments strictly (too fragile), just that it was called
            mock_perform.assert_called_once()
            print("[Test] perform_enrichment was called successfully.")
            
            # 6. Verify the Prompt Output contains the hints
            print("[Test] Verifying prompt content...")
            self.assertIn("Functional Hints from Database", prompt)
            self.assertIn("T cell activation", prompt)
            self.assertIn("GO:0003", prompt)
            
            # 7. Verify enrichment_hints return
            self.assertEqual(len(enrichment_hints), 2)
            self.assertIn("0", enrichment_hints)
            
            print("[Test] Prompt contains correct enrichment hints!")
            print("-" * 50)
            print("Snippet of Generated Prompt:")
            print(prompt.split('\n')[0:15]) # Print first few lines
            print("...")

if __name__ == '__main__':
    unittest.main()
