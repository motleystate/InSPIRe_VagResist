#!/usr/bin/env python3

"""
Calculate the Phenotypic Resistance Diversity Index (PRDI).

PRDI quantifies the diversity of predicted antimicrobial resistance
phenotypes represented by acquired antimicrobial resistance genes (ARGs).

(1) This script converts ResFinder gene names to their corresponding accession
    numbers using the mapping dictionary provided in `data/`.

(2) The calculation is based on Faith's phylogenetic diversity formula applied to a
    phenotype-informed ARG tree (PhenoARGTree) also provided in `data/`.

(3) The output is a TSV file containing each SampleID and its corresponding
    calculated PRDI value.
"""

import argparse
import json
import skbio
from pathlib import Path

import pandas as pd
from skbio import TreeNode
from skbio.diversity.alpha import faith_pd


# ---------------------------------------------------------------------
# Resource files
# ---------------------------------------------------------------------

script = Path(__file__).resolve().parent

tree = script / "data" / "PhenoARGTree.newick"
dict = script / "data" / "RF_gene_info_dict.json"


# ---------------------------------------------------------------------
# Load resources
# ---------------------------------------------------------------------

def load_resfinder_dict(dict_file):

    with open(dict_file, "r", encoding="utf-8") as handle:
        rf_dict = json.load(handle)

    return {
        entry["Gene"]: entry["Accession"]
        for entry in rf_dict
    }


def load_tree(tree_file):

    pheno_tree = TreeNode.read(tree_file)
    pheno_tree = pheno_tree.root_at(pheno_tree)

    return pheno_tree


# ---------------------------------------------------------------------
# Prepare ARG abundance table
# ---------------------------------------------------------------------

def convert_genes_to_accessions(arg_table, gene_to_accession):
    """
    Convert ResFinder gene names to accession numbers used in PhenoARGTree.
    """

    table = arg_table.copy()

    table.columns = table.columns.map(gene_to_accession)

    # ARGs sharing the same accession are collapsed by summing their counts,
    # as they share the same predicted phenotypes.
    table = table.T.groupby(level=0).sum().T

    # Match accession formatting used in PhenoARGTree
    table.columns = (
        table.columns
        .str.replace(r"^NC_", "NC ", regex=True)
        .str.replace(r"^NG_", "NG ", regex=True)
    )

    return table


# ---------------------------------------------------------------------
# Check accession IDs
# ---------------------------------------------------------------------

def check_tree_ids(arg_table, pheno_tree):

    tree_ids = {tip.name for tip in pheno_tree.tips()}
    arg_ids = set(arg_table.columns)

    missing_ids = sorted(arg_ids - tree_ids)

    if missing_ids:

        print(
            f"Error: {len(missing_ids)} accession(s) were not found "
            "in PhenoARGTree:"
        )

        for accession in missing_ids:
            print(f"  - {accession}")

        raise ValueError(
            "PRDI cannot be calculated because some accession IDs are absent from "
            "PhenoARGTree. Please check that the accession IDs exist in the tree "
            "and are correctly formatted."
        )

# ---------------------------------------------------------------------
# PRDI calculation
# ---------------------------------------------------------------------

def calculate_prdi(arg_table, pheno_tree):

    taxa_ids = arg_table.columns.tolist()
    results = {}

    for sample_id in arg_table.index:

        abundances = (arg_table.loc[sample_id, taxa_ids].astype(float).to_numpy())

        # Samples without detected ARGs are assigned a PRDI of 0.
        if abundances.sum() == 0:
            results[sample_id] = 0.0

        else:
            results[sample_id] = faith_pd(
                abundances,
                taxa_ids,
                pheno_tree,
                validate=True
            )

    return pd.Series(results, name="PRDI")


# ---------------------------------------------------------------------
# Command-line
# ---------------------------------------------------------------------

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Calculate PRDI from an ARG abundance table."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input ARG abundance table (TSV)."
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output PRDI table (TSV)."
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    args = parse_arguments()

    # Load ARG abundance table
    print("Loading ARG abundance table...")
    arg_table = pd.read_csv(
        args.input,
        sep="\t",
        index_col=0
    )

    # Load PRDI resources
    gene_to_accession = load_resfinder_dict(dict)

    print("Loading PhenoARGTree (this may take a few minutes)...")
    pheno_tree = load_tree(tree)
    print("\nPhenoARGTree loaded successfully.")

    # Convert ResFinder gene names to accession numbers
    arg_table_accessions = convert_genes_to_accessions(
        arg_table,
        gene_to_accession
    )

    # Check accession IDs
    check_tree_ids(
        arg_table_accessions,
        pheno_tree
    )
    print("\nAll accession IDs were found in PhenoARGTree.")

    # Calculate PRDI
    print("\nCalculating PRDI (this may take a few minutes)...")
    prdi = calculate_prdi(
        arg_table_accessions,
        pheno_tree
    )

    # Save results
    prdi.to_frame().to_csv(
        args.output,
        sep="\t",
        index_label="SampleID"
    )

    print(
        f"\nPRDI calculated for {len(prdi)} samples : {args.output}"
    )
    print(
        "\nCite us: Boutouchent N, Baud A, Tazi A, et al. "
        "Bacterial community structure shapes the vaginal resistome during pregnancy. "
        "npj Biofilms and Microbiomes (2026). "
        "https://doi.org/10.1038/s41522-026-01144-y"
    )

if __name__ == "__main__":
    main()