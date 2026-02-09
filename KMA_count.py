###################################### HELLO! ###########################################################

# This script post-process KMA mapping results to quantify genes frequency across multiple samples.

# ------------
# Input
# ------------
# The base directory must contain one subdirectory per sample. Each sample directory must include:
#     resfinder_kma/
#         kma_*.frag.gz   Read-level alignments/hits
#         kma_*.res       Per-template summary statistics (coverage/identity, etc.)

# -----------
# Steps 
# -----------
# 1) Gene filtering (template-level quality control)
#    Genes (templates) are retained if their estimated coverage and identity in the corresponding
#    *.res file meet user-defined thresholds (default: coverage >= 60%, identity >= 80%).
#    These thresholds are commonly used in the resistome literature, but should be adapted to the
#    organism, sequencing protocol, and research question.

# 2) Read-pair aggregation and assignment
#    Hits are aggregated by read pair identifier (R1/R2). For each read pair:
#    - if all retained hits map to a single gene, the pair contributes 1 count to that gene;
#    - if the pair maps to multiple genes, the gene with the highest KMA bit score is selected
#      (ties are tracked in the summary report).

# 3) Per-sample gene counting
#    Gene counts are computed per sample and aggregated into a sample-by-gene matrix.

# ---------
# Outputs
# ---------

# 1) Gene abundance table (CSV)
#    A sample-by-gene table where values correspond to read-pair counts assigned to each gene.

# 2) Mapping summary report (CSV)
#    Per-sample summary including numbers of read pairs processed, retained, and the frequency
#    of multi-hit assignments and score ties.

############################################################################################################

# ---------------------------------
# Parsing, filtering and read count 
# ---------------------------------

import os
import argparse
import glob
import gzip
import argparse
import pandas as pd
from collections import defaultdict

# Gene filtering from the *.res file based on coverage and identity thresholds.
def filter_genes(res_file, coverage_threshold, ID_threshold):
    """
    *.frag.gz contains the following information:
      - read sequence (column 0)
      - number of equivalent templates (column 1)
      - mapping score (bit score, column 2)
      - start position (column 3)
      - end position (column 4)
      - selected template (gene, column 5)
      - read ID (column 6)
      - additional information (column 7)
    """
    gene_quality = {}
    with open(res_file, "r") as res:
        for line in res:
            if line.startswith("#"):
                continue
            parts = line.split()
            gene = parts[0]
            coverage = float(parts[5])
            identity = float(parts[4])
            q_value = float(parts[9])
            p_value = float(parts[10])
            if coverage >= coverage_threshold and identity >= ID_threshold:
                gene_quality[gene] = (q_value, p_value)
    return gene_quality

# read filtering from the *.frag.gz file by keeping genes in gene_quality dictionary.
def filter_reads(frag_file, gene_quality):
    counted_reads = set()
    filtered_reads = []
    read_to_hits = defaultdict(list)
    total_reads_before = set()
    multiple_hit_reads = 0
    same_score_count = 0

    with gzip.open(frag_file, "rt") as frag:
        for line in frag:
            parts = line.split()
            read_id = parts[6]
            gene = parts[5]
            if gene not in gene_quality:
                continue
            bit_score = float(parts[2])
            paired_read_id = read_id.rsplit(" ", 1)[0]
            total_reads_before.add(paired_read_id)
            read_to_hits[paired_read_id].append((gene, bit_score))

    for read_id, hits in read_to_hits.items():
        gene_to_score = {}
        for gene, score in hits:
            if gene in gene_to_score:
                if score > gene_to_score[gene]:
                    gene_to_score[gene] = score
            else:
                gene_to_score[gene] = score

        if len(gene_to_score) == 1:
            best_gene = next(iter(gene_to_score))
        else:
            best_gene = None
            best_score = float("-inf")
            for gene, score in gene_to_score.items():
                if score > best_score:
                    best_score = score
                    best_gene = gene

            if sum(1 for s in gene_to_score.values() if s == best_score) > 1:
                same_score_count += 1
            multiple_hit_reads += 1

        if read_id not in counted_reads:
            filtered_reads.append((read_id, best_gene))
            counted_reads.add(read_id)

    return filtered_reads, multiple_hit_reads, same_score_count, len(total_reads_before), len(counted_reads)


def count_genes(filtered_reads):
    gene_counts = defaultdict(int)
    for _, gene in filtered_reads:
        gene_counts[gene] += 1
    return gene_counts

# ---------------------------
# CLI / main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(description="Post-processing KMA outputs to build a gene frequency table.")
    parser.add_argument("-p", "--base_path", required=True, help="Path to the directory with sample subfolders.")
    parser.add_argument("-cov", "--coverage", type=float, default=60, help="Coverage threshold (%) (default: 60).")
    parser.add_argument("-ID", "--identity", type=float, default=80, help="Identity threshold (%) (default: 80).")
    parser.add_argument("--out_table", required=True)
    parser.add_argument("--out_report", required=True)
    args = parser.parse_args()

    base_path = args.base_path
    cov = args.coverage
    ID = args.identity

    print("Starting ResFinder_count.py!")
    print(f"Base path: {base_path}")
    print(f"Thresholds: coverage ≥ {cov}%, identity ≥ {ID}%")

    all_gene_counts = defaultdict(lambda: defaultdict(int))
    all_genes = set()
    report_data = []

    for sample_path in glob.glob(os.path.join(base_path, "*")):
        if not os.path.isdir(sample_path):
            continue

        sample_name = os.path.basename(sample_path)
        print(f"\nProcessing sample: {sample_name}")

        frag_files = glob.glob(os.path.join(sample_path, "resfinder_kma", "kma_*.frag.gz"))

        total_multi_hits = 0
        total_same_score = 0
        total_reads_before = 0
        total_reads_after = 0

        for frag_file in frag_files:
            res_file = os.path.join(sample_path,"resfinder_kma",
                                    f"kma_{os.path.basename(frag_file).replace('kma_', '').replace('.frag.gz', '')}.res")

            gene_quality = filter_genes(res_file, coverage_threshold=cov, ID_threshold=ID)
            filtered_reads, multi_hit_reads, same_score_reads, reads_before, reads_after = filter_reads(
                frag_file, gene_quality)
            gene_counts = count_genes(filtered_reads)

            total_multi_hits += multi_hit_reads
            total_same_score += same_score_reads
            total_reads_before += reads_before
            total_reads_after += reads_after

            for gene, count in gene_counts.items():
                all_gene_counts[sample_name][gene] = count
                all_genes.add(gene)

        report_data.append([sample_name, total_reads_before, total_reads_after, total_multi_hits, total_same_score])
        print("Done.")

    rows = []
    for sample_name in all_gene_counts:
        rows.append([sample_name] + [all_gene_counts[sample_name].get(gene, 0) for gene in all_genes])

    df = pd.DataFrame(rows, columns=["SampleID"] + list(all_genes))
    df.to_csv(args.out_table, index=False)
    print(f"\nGene frequency table saved: {args.out_table}")

    report_df = pd.DataFrame(report_data, 
                             columns=["SampleID", "Total_Reads_Before", "Total_Reads_After", "Reads_Multi_Hits", "Multi_Hit_Same_Score"],)
    report_df.to_csv(args.out_report, index=False)
    print(f"Mapping report saved: {args.out_report}")

if __name__ == "__main__":
    main()