from datetime import datetime
import os
from snakemake.io import glob_wildcards

configfile: "config.yaml"

SOURCE_DIR = config["input_data"]["source_dir"]
OUTPUT_DIR = config["paths"].get("output_dir", "mooncrater_output")
OUTPUT_MAPPING_DIR = os.path.join(OUTPUT_DIR, "ResFinder_mapping_out")
COUNT_OUTDIR = os.path.join(OUTPUT_DIR, "count_out")

DB_RF = config["paths"]["db_resfinder"]


START_FILE = ".mooncrater_start_time"
def return_runtime(start, end):
    string = []
    delta_seconds = (end - start).seconds  # omitting microseconds
    d = delta_seconds // 86400
    if d > 0:
        delta_seconds -= (d * 86400)
        string += [f"{d}d"]
    h = delta_seconds // 3600
    if h > 0:
        delta_seconds -= (h * 3600)
        string += [f"{h}h"]
    m = delta_seconds // 60
    if m > 0:
        delta_seconds -= (m * 60)
        string += [f"{m}m"]
    string += [f"{delta_seconds}s"]
    return " ".join(string)


coverage = int(config["kma_count"]["coverage"])
identity = int(config["kma_count"]["identity"])


# How samples and paired reads will be detected with wildcards.
__input_suffix = config["input_data"].get("suffix", ".fastq.gz")
__input_prefix = config["input_data"].get("pair_prefix", "_")
if config["input_data"].get("paired", True):
    (SAMPLES,) = glob_wildcards("{dir}/{{sample}}{prefix}1{suffix}".format(dir=SOURCE_DIR, prefix=__input_prefix, suffix=__input_suffix))
    __run_resfinder_input = expand("{dir}/{sample}{prefix}{{read}}{suffix}".format(dir=SOURCE_DIR, sample="{{sample}}", prefix=__input_prefix, suffix=__input_suffix), read=["1", "2"])
else:
    (SAMPLES,) = glob_wildcards("{dir}/{{sample}}{suffix}".format(dir=__input_dir, suffix=__input_suffix))
    __run_resfinder_input = "{dir}/{sample}{suffix}".format(dir=SOURCE_DIR, sample="{sample}", suffix=__input_suffix)


__all_output = os.path.join(OUTPUT_DIR, "report.txt")
rule all:
    input:
        __all_output


onstart:
    print("Welcome to MoonCrater pipeline!", flush=True)
    with open(START_FILE, 'w') as f:
        f.write(datetime.now().isoformat())


rule run_resfinder:
    input:
        __run_resfinder_input
    output:
        out_RF = directory(os.path.join(OUTPUT_MAPPING_DIR, "{sample}"))
    params:
        input_format = config["input_data"]["format"]
    threads: 8
    run:
        if params.input_format == "fasta":
            input_arg = "-ifa"
        elif params.input_format == "fastq":
            input_arg = "-ifq"
        else:
            raise ValueError("format of input data files should be 'fasta' or 'fastq'")
        print(f"[MoonCrater] Running ResFinder/KMA for sample: {wildcards.sample}...")
        shell(r'run_resfinder.py {input_arg} {input} -acq -s "Other" -db_res {DB_RF} -o {output.out_RF}')
        print(f"[MoonCrater] Finished ResFinder/KMA for sample: {wildcards.sample}")

__read_count_outtable = os.path.join(COUNT_OUTDIR, f"gene_abundance_table_Cov{coverage}_ID{identity}.csv")
__read_count_outreport = os.path.join(COUNT_OUTDIR, "read_mapping_report.csv")
rule read_count:
    input:
        expand(os.path.join(OUTPUT_MAPPING_DIR, "{sample}"), sample=SAMPLES)
    output:
        table = __read_count_outtable,
        report = __read_count_outreport
    params:
        cov = coverage,
        ident = identity
    threads: 4
    shell:
        """
        echo "[MoonCrater] Aggregating KMA results across samples"
        echo "[MoonCrater] Thresholds applied: coverage ≥ {params.cov}%, identity ≥ {params.ident}%" 
        python KMA_count.py -p "{OUTPUT_MAPPING_DIR}" -cov {params.cov} -ID {params.ident} --out_table "{output.table}" --out_report "{output.report}"
        """

rule run_completed:
    input:
        table = __read_count_outtable,
        report = __read_count_outreport
    output:
        __all_output
    run:
        from datetime import datetime
        import time
        import os

        outfile = open(str(output), "w")
        ending_date = datetime.now()
        offset = time.timezone if (time.localtime().tm_isdst == 0) else time.altzone
        offset = offset / -3600
        outfile.write(ending_date.strftime('[MoonCrater] ended at %Y-%m-%d (YMD) %H:%M:%S ')+f"UTC+{offset}")
        if os.path.exists(START_FILE):
            with open(START_FILE) as f:
                starting_date = datetime.fromisoformat(f.read().strip())
            outfile.write(f" after running for {return_runtime(starting_date, ending_date)}\n")
            os.remove(START_FILE)
        else:
            outfile.write("\n")
        outfile.write(f"KMA parameters:\n  - coverage: {coverage}\n  - identity: {identity}\n")
        samples = SAMPLES
        outfile.write(f"Data:\n  - {len(samples)} samples: {', '.join(samples)}\n")
        dt_m = datetime.fromtimestamp(os.path.getmtime(DB_RF)).strftime('%Y-%m-%d %H:%M')
        outfile.write(f"  - ResFinder database: {DB_RF} [{dt_m}]")
        outfile.close()
        print("[MoonCrater] We hope this pipeline was useful for your analyses!")
        print("[MoonCrater] Questions or suggestions? Please contact: Nassim Boutouchent (nassim.boutouchent@chu-rouen.fr)")
