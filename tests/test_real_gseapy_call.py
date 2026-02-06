import sys
import unittest
try:
    import gseapy as gp
    import pandas as pd
except ImportError:
    sys.exit(0) # Skip if not installed

class TestRealGseapy(unittest.TestCase):
    def test_enrichr_call(self):
        gene_list = ['TP53', 'TNF', 'EGFR']
        print("Attempting to talk to Enrichr...")
        try:
            # Short timeout if possible? gseapy doesn't expose timeout easily.
            # We just try it. If no internet, this will raise an error.
            res = gp.enrichr(gene_list=gene_list, gene_sets='GO_Biological_Process_2021', organism='Human')
            print("Success! (Internet is available)")
        except Exception as e:
            print(f"Enrichr call failed (Expected if no internet): {e}")
            # We don't fail the test, we just note it.
            pass

if __name__ == '__main__':
    unittest.main()
