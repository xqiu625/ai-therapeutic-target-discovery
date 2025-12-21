#!/usr/bin/env python
"""
Geneformer Embedding Generation Script
Run with Singularity: singularity exec geneformer.sif python scripts/run_geneformer.py

Generates cell embeddings using Geneformer foundation model.
"""

import scanpy as sc
import numpy as np
import os
from pathlib import Path
import argparse
import warnings
import tempfile
import shutil
warnings.filterwarnings('ignore')


def convert_gene_names_to_ensembl(adata):
    """
    Convert gene names to Ensembl IDs if needed.
    Geneformer requires Ensembl IDs.
    """
    # Check if already Ensembl IDs
    if adata.var_names[0].startswith('ENSG'):
        print("Gene names already in Ensembl format")
        return adata

    # Try to use gene_ids if available
    if 'gene_ids' in adata.var.columns:
        print("Using gene_ids column for Ensembl IDs")
        # Create mapping
        adata.var['original_names'] = adata.var_names.copy()
        adata.var_names = adata.var['gene_ids'].values
        return adata

    print("WARNING: Gene names not in Ensembl format and no gene_ids column found")
    print("Geneformer may not work correctly without Ensembl IDs")
    return adata


def main(args):
    print("="*60)
    print("GENEFORMER EMBEDDING GENERATION")
    print("="*60)

    # Paths
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"\nLoading data from: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Prepare data for Geneformer
    adata_gf = adata.copy()

    # Use raw counts
    if 'counts' in adata.layers:
        adata_gf.X = adata.layers['counts'].copy()
        print("Using raw counts from layers['counts']")
    else:
        print("WARNING: No raw counts layer found")

    # Convert gene names to Ensembl IDs
    adata_gf = convert_gene_names_to_ensembl(adata_gf)

    # Try to import and run Geneformer
    try:
        from geneformer import TranscriptomeTokenizer, EmbExtractor
        print("Geneformer imported successfully")

        # Create temporary directory for tokenized data
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Save adata in format Geneformer expects
            data_dir = tmpdir / 'data'
            data_dir.mkdir()
            adata_gf.write(data_dir / 'data.h5ad')

            # Tokenize
            print("\nTokenizing data...")
            tokenizer = TranscriptomeTokenizer(
                custom_attr_name_dict={'cell_type': 'cell_type'} if 'cell_type' in adata_gf.obs.columns else None,
                model_input_size=2048
            )

            tokenized_dir = tmpdir / 'tokenized'
            tokenized_dir.mkdir()

            tokenizer.tokenize_data(
                data_directory=str(data_dir),
                output_directory=str(tokenized_dir),
                output_prefix='sepsis'
            )
            print("Tokenization complete")

            # Extract embeddings
            print("\nExtracting embeddings...")
            extractor = EmbExtractor(
                model_type='Pretrained',
                num_classes=0,  # For embeddings only
                emb_layer=-1,   # Last layer
                emb_mode='cell',
                forward_batch_size=args.batch_size
            )

            emb_dir = tmpdir / 'embeddings'
            emb_dir.mkdir()

            embeddings = extractor.extract_embs(
                model_directory=args.model_dir,
                input_data_file=str(tokenized_dir / 'sepsis.dataset'),
                output_directory=str(emb_dir),
                output_prefix='geneformer'
            )

            # Convert to numpy if needed
            if hasattr(embeddings, 'values'):
                embeddings = embeddings.values
            embeddings = np.array(embeddings)

        print(f"Geneformer embeddings shape: {embeddings.shape}")

        # Add to adata
        adata.obsm['X_geneformer'] = embeddings

        # Save embeddings
        output_path = output_dir / 'geneformer_embeddings.npy'
        np.save(output_path, embeddings)
        print(f"\nSaved embeddings to: {output_path}")

        # Save updated adata
        adata_output = output_dir / 'adata_geneformer.h5ad'
        adata.write(adata_output)
        print(f"Saved adata with embeddings to: {adata_output}")

        print("\nGeneformer embedding generation complete!")
        return True

    except ImportError as e:
        print(f"\nGeneformer import error: {e}")
        print("Make sure you're running this inside the Geneformer Singularity container")
        return False

    except Exception as e:
        print(f"\nERROR generating Geneformer embeddings: {e}")
        import traceback
        traceback.print_exc()

        # Save placeholder
        print("\nSaving placeholder embeddings...")
        embeddings = np.zeros((adata.n_obs, 512))
        adata.obsm['X_geneformer'] = embeddings
        np.save(output_dir / 'geneformer_embeddings.npy', embeddings)
        adata.write(output_dir / 'adata_geneformer.h5ad')

        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate Geneformer embeddings')
    parser.add_argument('--input', type=str, required=True,
                        help='Path to input h5ad file')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for embeddings')
    parser.add_argument('--model_dir', type=str, default='ctheodoris/Geneformer',
                        help='Path to Geneformer model or HuggingFace model name')
    parser.add_argument('--batch_size', type=int, default=32,
                        help='Batch size for embedding extraction')

    args = parser.parse_args()
    main(args)
