import sys
try:
    import gseapy
    print(f"SUCCESS: gseapy version {gseapy.__version__} is installed and working.")
except ImportError as e:
    print(f"FAILURE: Could not import gseapy. Error: {e}")
except Exception as e:
    print(f"FAILURE: An unexpected error occurred while importing gseapy: {e}")
