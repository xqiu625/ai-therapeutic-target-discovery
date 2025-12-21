#!/usr/bin/env python
"""
Download GSE167363 data from GEO.
This dataset contains scRNA-seq data from sepsis patients (survivors + non-survivors) and healthy controls.

Usage:
    python scripts/download_data.py

The script downloads:
    - Raw count matrices (10X format)
    - Sample metadata
"""

import os
import subprocess
import tarfile
import gzip
import shutil
from pathlib import Path

# Configuration
GEO_ACCESSION = "GSE167363"
RAW_DATA_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")

# GEO FTP URLs for GSE167363
# Note: These are the supplementary files containing the 10X data
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167363/suppl/"


def download_from_geo():
    """Download data files from GEO."""

    RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Downloading GSE167363 data...")
    print("=" * 60)

    # Method 1: Use GEOparse if available
    try:
        import GEOparse
        print("Using GEOparse to download metadata...")
        gse = GEOparse.get_GEO(geo=GEO_ACCESSION, destdir=str(RAW_DATA_DIR))
        print(f"Downloaded metadata for {GEO_ACCESSION}")

        # Print sample information
        print("\nSamples in dataset:")
        for gsm_name, gsm in gse.gsms.items():
            title = gsm.metadata.get('title', ['Unknown'])[0]
            print(f"  {gsm_name}: {title}")

    except ImportError:
        print("GEOparse not available, using wget...")
    except Exception as e:
        print(f"GEOparse failed: {e}, trying wget...")

    # Method 2: Direct download of supplementary files
    print("\nDownloading supplementary files...")

    # The main tar file containing all 10X data
    supp_file = f"{GEO_ACCESSION}_RAW.tar"
    supp_url = f"{BASE_URL}{supp_file}"

    tar_path = RAW_DATA_DIR / supp_file

    if not tar_path.exists():
        print(f"Downloading {supp_file}...")
        try:
            subprocess.run([
                "wget", "-c", "-P", str(RAW_DATA_DIR), supp_url
            ], check=True)
        except subprocess.CalledProcessError:
            # Try curl as fallback
            subprocess.run([
                "curl", "-L", "-o", str(tar_path), supp_url
            ], check=True)
    else:
        print(f"{supp_file} already exists, skipping download.")

    # Extract tar file
    if tar_path.exists():
        print(f"Extracting {supp_file}...")
        with tarfile.open(tar_path) as tar:
            tar.extractall(path=RAW_DATA_DIR)
        print("Extraction complete!")

    print("\n" + "=" * 60)
    print("Download complete!")
    print(f"Data saved to: {RAW_DATA_DIR}")


def download_scenic_databases():
    """Download pySCENIC reference databases for human."""

    EXTERNAL_DIR = Path("data/external")
    EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)

    print("\nDownloading pySCENIC databases (human)...")
    print("=" * 60)

    # cisTarget database URLs (hg38)
    databases = {
        "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather":
            "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather":
            "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl":
            "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    }

    # TF list
    tf_url = "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt"

    # Download TF list
    tf_path = EXTERNAL_DIR / "hs_hgnc_tfs.txt"
    if not tf_path.exists():
        print("Downloading TF list...")
        subprocess.run(["wget", "-O", str(tf_path), tf_url], check=True)

    # Download databases (these are large files, ~1-2GB each)
    for filename, url in databases.items():
        filepath = EXTERNAL_DIR / filename
        if not filepath.exists():
            print(f"Downloading {filename}...")
            print("(This may take a while, files are 1-2GB each)")
            subprocess.run(["wget", "-c", "-O", str(filepath), url], check=True)
        else:
            print(f"{filename} already exists, skipping.")

    print("\n" + "=" * 60)
    print("pySCENIC databases downloaded!")


def verify_download():
    """Verify that all required files are present."""

    print("\nVerifying downloaded files...")
    print("=" * 60)

    # Check raw data
    raw_files = list(RAW_DATA_DIR.glob("*.gz"))
    print(f"Found {len(raw_files)} compressed files in {RAW_DATA_DIR}")

    for f in raw_files[:10]:  # Show first 10
        print(f"  - {f.name}")
    if len(raw_files) > 10:
        print(f"  ... and {len(raw_files) - 10} more")

    # Check external databases
    external_dir = Path("data/external")
    if external_dir.exists():
        external_files = list(external_dir.glob("*"))
        print(f"\nFound {len(external_files)} files in {external_dir}")
        for f in external_files:
            size_mb = f.stat().st_size / (1024 * 1024)
            print(f"  - {f.name} ({size_mb:.1f} MB)")


def main():
    """Main function to download all required data."""

    print("=" * 60)
    print("SEPSIS TARGET DISCOVERY - DATA DOWNLOAD")
    print("=" * 60)
    print(f"GEO Accession: {GEO_ACCESSION}")
    print(f"Output directory: {RAW_DATA_DIR}")
    print("=" * 60)

    # Download GEO data
    download_from_geo()

    # Download SCENIC databases (optional, can be slow)
    response = input("\nDownload pySCENIC databases? (y/n, ~5GB total): ")
    if response.lower() == 'y':
        download_scenic_databases()
    else:
        print("Skipping pySCENIC databases. You can download later with:")
        print("  python scripts/download_data.py --scenic-only")

    # Verify
    verify_download()

    print("\n" + "=" * 60)
    print("Data download complete!")
    print("Next step: Run notebook 01_data_preparation.ipynb")
    print("=" * 60)


if __name__ == "__main__":
    import sys

    if "--scenic-only" in sys.argv:
        download_scenic_databases()
    else:
        main()
