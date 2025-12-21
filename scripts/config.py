"""
Configuration for AI-Enabled Therapeutic Target Discovery Pipeline

Modify BASE_DIR to point to your project directory.
"""

from pathlib import Path
import os

# Auto-detect base directory (parent of scripts folder)
BASE_DIR = Path(__file__).parent.parent.resolve()

# Or set manually via environment variable
if os.environ.get('TARGET_DISCOVERY_DIR'):
    BASE_DIR = Path(os.environ['TARGET_DISCOVERY_DIR'])

# Directory structure
DATA_DIR = BASE_DIR / 'data'
RAW_DIR = DATA_DIR / 'raw'
PROCESSED_DIR = DATA_DIR / 'processed'
EXTERNAL_DIR = DATA_DIR / 'external'
RESULTS_DIR = BASE_DIR / 'results'
TABLES_DIR = RESULTS_DIR / 'tables'
FIGURES_DIR = BASE_DIR / 'figures'

# Create directories if they don't exist
for d in [RAW_DIR, PROCESSED_DIR, EXTERNAL_DIR, TABLES_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Print configuration on import (optional)
if __name__ == '__main__':
    print(f"BASE_DIR: {BASE_DIR}")
    print(f"DATA_DIR: {DATA_DIR}")
    print(f"RESULTS_DIR: {RESULTS_DIR}")
