import os 
from snakemake.io import glob_wildcards

configfile: "config.yaml"

SOURCE_DIR   = config["paths"]["source_dir"]
OUTPUT_MAPPING_DIR  = config["paths"]["output_dir"]
DB_RF = config["paths"]["db_resfinder"]
COUNT_OUTDIR = config["paths"]["count_outdir"]

coverage = int(config["kma_count"]["coverage"])
identity = int(config["kma_count"]["identity"])

Samples, = glob_wildcards(os.path.join(config["paths"]["source_dir"], "{sample}_R1_trimmed_nohost.fastq.gz"))

print ("Welcome to MoonCrater pipeline !")
rule all:
    input:
        "done.txt"

rule run_resfinder: 
    input: 
        trimmed_nohost_R1 = os.path.join(SOURCE_DIR, "{sample}_R1_trimmed_nohost.fastq.gz"),
        trimmed_nohost_R2 = os.path.join(SOURCE_DIR, "{sample}_R2_trimmed_nohost.fastq.gz")
    output:
        out_RF = directory(os.path.join(OUTPUT_MAPPING_DIR, "{sample}"))
    threads: 8

    shell: 
        """
        echo "[MoonCrater] Running ResFinder/KMA for sample: {wildcards.sample}..."
        mkdir -p {OUTPUT_MAPPING_DIR}/{wildcards.sample}
        run_resfinder.py \
            -ifq {input.trimmed_nohost_R1} {input.trimmed_nohost_R2} \
            -acq \
            -s "Other" \
            -db_res {DB_RF} \
            -o {output.out_RF}
        echo "[MoonCrater] Finished ResFinder/KMA for sample: {wildcards.sample}"
        """

rule read_count:
    input: 
        expand(os.path.join(OUTPUT_MAPPING_DIR, "{sample}"), sample=Samples),
        script = "KMA_count.py" 
    output:
        table = os.path.join(COUNT_OUTDIR,f"gene_abundance_table_Cov{coverage}_ID{identity}.csv"),
        report = os.path.join(COUNT_OUTDIR, f"read_mapping_report.csv")
    threads: 4
    shell: 
        """
        echo "[MoonCrater] Aggregating KMA results across samples"
        echo "[MoonCrater] Thresholds applied: coverage ≥ {coverage}%, identity ≥ {identity}%" 
        mkdir -p {COUNT_OUTDIR}
        python {input.script} -p "{OUTPUT_MAPPING_DIR}" -cov {coverage} -ID {identity} --out_table "{output.table}" --out_report "{output.report}"

        """
rule run_completed:
    input: 
        table = os.path.join(COUNT_OUTDIR,f"gene_abundance_table_Cov{coverage}_ID{identity}.csv"),
        report = os.path.join(COUNT_OUTDIR, f"read_mapping_report.csv")
    output:
        "done.txt"
    shell:
        """
        touch {output}
        echo "[MoonCrater] We hope this pipeline was useful for your analyses!"
        echo "[MoonCrater] Questions or suggestions? please contact : Nassim Boutouchent (nassim.boutouchent@chu-rouen.fr)"
        """