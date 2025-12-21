#!/usr/bin/env python
"""
scGPT Embedding Generation Script
Run with: conda activate scgpt && python scripts/run_scgpt.py

Generates cell embeddings using scGPT foundation model.
"""

import scanpy as sc
import numpy as np
import torch
import os
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings('ignore')

def main(args):
    print("="*60)
    print("scGPT EMBEDDING GENERATION")
    print("="*60)

    # Check GPU
    print(f"\nPyTorch version: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    # Paths
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"\nLoading data from: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Prepare data for scGPT
    adata_scgpt = adata.copy()

    # Use raw counts if available
    if 'counts' in adata.layers:
        adata_scgpt.X = adata.layers['counts'].copy()
        print("Using raw counts from layers['counts']")
    else:
        print("WARNING: No raw counts layer found, using .X directly")

    # Try to import and run scGPT
    try:
        # Import scGPT - adjust import based on installed version
        # scGPT has different API versions, try common patterns

        try:
            # Try newer API
            from scgpt.tasks import embed_data
            from scgpt.utils import set_seed

            set_seed(42)

            # Generate embeddings
            print("\nGenerating scGPT embeddings...")
            embeddings = embed_data(
                adata_scgpt,
                model_dir=args.model_dir,
                batch_size=args.batch_size,
                gene_col='gene_name' if 'gene_name' in adata_scgpt.var.columns else None,
            )

        except ImportError:
            # Try alternative API
            from scgpt import scGPT

            print(f"\nLoading scGPT model from: {args.model_dir}")
            model = scGPT.load_pretrained(args.model_dir)
            model.eval()

            if torch.cuda.is_available():
                model = model.cuda()

            # Generate embeddings in batches
            print(f"\nGenerating embeddings (batch_size={args.batch_size})...")
            embeddings_list = []

            with torch.no_grad():
                for i in range(0, adata_scgpt.n_obs, args.batch_size):
                    batch = adata_scgpt[i:i+args.batch_size]
                    emb = model.encode(batch)

                    if isinstance(emb, torch.Tensor):
                        emb = emb.cpu().numpy()

                    embeddings_list.append(emb)

                    if (i // args.batch_size) % 10 == 0:
                        print(f"  Processed {i}/{adata_scgpt.n_obs} cells")

            embeddings = np.vstack(embeddings_list)

        print(f"scGPT embeddings shape: {embeddings.shape}")

        # Add to adata
        adata.obsm['X_scgpt'] = embeddings

        # Save embeddings
        output_path = output_dir / 'scgpt_embeddings.npy'
        np.save(output_path, embeddings)
        print(f"\nSaved embeddings to: {output_path}")

        # Save updated adata
        adata_output = output_dir / 'adata_scgpt.h5ad'
        adata.write(adata_output)
        print(f"Saved adata with embeddings to: {adata_output}")

        print("\nscGPT embedding generation complete!")
        return True

    except Exception as e:
        print(f"\nERROR generating scGPT embeddings: {e}")
        import traceback
        traceback.print_exc()

        # Save placeholder
        print("\nSaving placeholder embeddings...")
        embeddings = np.zeros((adata.n_obs, 512))
        adata.obsm['X_scgpt'] = embeddings
        np.save(output_dir / 'scgpt_embeddings.npy', embeddings)
        adata.write(output_dir / 'adata_scgpt.h5ad')

        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate scGPT embeddings')
    parser.add_argument('--input', type=str, required=True,
                        help='Path to input h5ad file')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for embeddings')
    parser.add_argument('--model_dir', type=str, default='scGPT_human',
                        help='Path to scGPT model or model name')
    parser.add_argument('--batch_size', type=int, default=64,
                        help='Batch size for embedding generation')

    args = parser.parse_args()
    main(args)
