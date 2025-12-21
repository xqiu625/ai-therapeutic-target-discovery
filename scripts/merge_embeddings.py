#!/usr/bin/env python
"""
Merge Foundation Model Embeddings

This script combines embeddings from scGPT and Geneformer into a single h5ad file.
Run after both embedding scripts have completed.
"""

import scanpy as sc
import numpy as np
from pathlib import Path
import argparse


def main(args):
    print("="*60)
    print("MERGING FOUNDATION MODEL EMBEDDINGS")
    print("="*60)

    # Paths
    processed_dir = Path(args.processed_dir)
    embeddings_dir = Path(args.embeddings_dir)

    # Load base processed data
    base_file = processed_dir / 'adata_processed.h5ad'
    print(f"\nLoading base data: {base_file}")
    adata = sc.read_h5ad(base_file)
    print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Load scGPT embeddings
    scgpt_file = embeddings_dir / 'scgpt_embeddings.npy'
    if scgpt_file.exists():
        scgpt_emb = np.load(scgpt_file)
        if scgpt_emb.shape[0] == adata.n_obs:
            adata.obsm['X_scgpt'] = scgpt_emb
            is_placeholder = np.all(scgpt_emb == 0)
            print(f"Loaded scGPT embeddings: {scgpt_emb.shape} {'(placeholder)' if is_placeholder else ''}")
        else:
            print(f"WARNING: scGPT embeddings shape mismatch ({scgpt_emb.shape[0]} vs {adata.n_obs})")
    else:
        print(f"scGPT embeddings not found: {scgpt_file}")

    # Load Geneformer embeddings
    gf_file = embeddings_dir / 'geneformer_embeddings.npy'
    if gf_file.exists():
        gf_emb = np.load(gf_file)
        if gf_emb.shape[0] == adata.n_obs:
            adata.obsm['X_geneformer'] = gf_emb
            is_placeholder = np.all(gf_emb == 0)
            print(f"Loaded Geneformer embeddings: {gf_emb.shape} {'(placeholder)' if is_placeholder else ''}")
        else:
            print(f"WARNING: Geneformer embeddings shape mismatch ({gf_emb.shape[0]} vs {adata.n_obs})")
    else:
        print(f"Geneformer embeddings not found: {gf_file}")

    # Determine best embedding for downstream analysis
    embeddings_available = []
    for emb_name in ['X_scgpt', 'X_geneformer']:
        if emb_name in adata.obsm and not np.all(adata.obsm[emb_name] == 0):
            embeddings_available.append(emb_name)

    if embeddings_available:
        # Use first available real embedding
        adata.obsm['X_embedding'] = adata.obsm[embeddings_available[0]]
        print(f"\nPrimary embedding set to: {embeddings_available[0]}")
    else:
        # Fall back to PCA
        if 'X_pca' in adata.obsm:
            adata.obsm['X_embedding'] = adata.obsm['X_pca'][:, :50]
            print("\nNo FM embeddings available, using PCA as fallback")
        else:
            print("\nWARNING: No embeddings available!")

    # Summary
    print("\n" + "="*60)
    print("EMBEDDING SUMMARY")
    print("="*60)
    for key in adata.obsm.keys():
        shape = adata.obsm[key].shape
        is_zero = np.all(adata.obsm[key] == 0) if 'X_' in key else False
        status = "(placeholder)" if is_zero else ""
        print(f"  {key}: {shape} {status}")

    # Save merged file
    output_file = processed_dir / 'adata_embeddings.h5ad'
    adata.write(output_file)
    print(f"\nSaved merged embeddings to: {output_file}")

    # Also save to embeddings directory
    adata.write(embeddings_dir / 'adata_embeddings.h5ad')

    print("\nMerge complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge foundation model embeddings')
    parser.add_argument('--processed_dir', type=str, required=True,
                        help='Path to processed data directory')
    parser.add_argument('--embeddings_dir', type=str, required=True,
                        help='Path to embeddings directory')

    args = parser.parse_args()
    main(args)
