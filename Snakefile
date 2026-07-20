# ==============================================================================
# IsoDecipher Snakemake pipeline
#
# Usage:
#   cd /Users/renegibson/Desktop/githubrepo/IsoMatrix/IsoDecipher
#   snakemake --configfile config/cheng_bcell.yaml --cores 4 -n   # dry run
#   snakemake --configfile config/cheng_bcell.yaml --cores 4       # run
#   snakemake --configfile config/lymphoma.yaml --cores 3          # lymphoma batch
#
# Re-run specific rule:
#   snakemake --configfile config/lymphoma.yaml --cores 3 \
#             --forcerun assign_reads
#
# Per-sample path overrides (see config/lymphoma.yaml):
#   Use `per_sample:` block when samples live in different directories or have
#   non-standard BAM names. Without per_sample, all paths derive from data_dir.
# ==============================================================================

import os

# ── Config ────────────────────────────────────────────────────────────────────
SAMPLES     = config["samples"]
DATA_DIR    = config["data_dir"]
GLOBAL_PAN  = config["global_panel"]
INDICES_PKL = config.get("indices_pkl")
SCRIPTS     = os.path.join(workflow.basedir, "IsoDecipher", "scripts")

BATCH_PAN   = os.path.join(DATA_DIR, config["batch_panel"])
OUT_H5AD    = os.path.join(DATA_DIR, config["output_h5ad"])
BATCH_LOGS  = os.path.join(DATA_DIR, "results", "logs")


# ── Per-sample config helpers ─────────────────────────────────────────────────
def _per(sample, key):
    """Return per-sample config override or None."""
    return config.get("per_sample", {}).get(sample, {}).get(key)

def sample_src_dir(sample):
    """Source dir for BAM + barcodes — may differ from DATA_DIR."""
    return _per(sample, "data_dir") or os.path.join(DATA_DIR, sample)

def bam_path(sample):
    bam = _per(sample, "bam") or config["bam"]
    return os.path.join(sample_src_dir(sample), bam)

def barcodes_path(sample):
    mdir = _per(sample, "matrix_dir") or config.get("matrix_dir", "filtered_feature_bc_matrix")
    return os.path.join(sample_src_dir(sample), mdir, "barcodes.tsv.gz")

def matrix_dir(sample):
    mdir = _per(sample, "matrix_dir") or config.get("matrix_dir", "filtered_feature_bc_matrix")
    return os.path.join(sample_src_dir(sample), mdir)

def insert_lookup(sample):
    """Pre-computed lookup takes priority; else expect validate_insert_size output."""
    pre = _per(sample, "insert_size_lookup")
    return pre or os.path.join(DATA_DIR, sample, "results", "insert_size",
                               f"{sample}_insert_size_panel_lookup.json")

def counts_base(sample):
    return os.path.join(DATA_DIR, sample, "results", "counts", sample)

def h5_path(sample):
    return _per(sample, "h5") or os.path.join(DATA_DIR, sample, "filtered_feature_bc_matrix.h5")


# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        OUT_H5AD


# ── Rule 1: prep_sample_panel ─────────────────────────────────────────────────
rule prep_sample_panel:
    input:
        global_panel = GLOBAL_PAN,
        matrix_dirs  = [matrix_dir(s) for s in SAMPLES],
    output:
        BATCH_PAN,
    log:
        os.path.join(BATCH_LOGS, "prep_panel.log"),
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        python {SCRIPTS}/prep_sample_panel.py \
            --global-panel {input.global_panel} \
            --matrix-dirs  {input.matrix_dirs} \
            --out          {output} \
            > {log} 2>&1
        """


# ── Rule 2: validate_insert_size (per sample, parallelizable) ────────────────
rule validate_insert_size:
    input:
        bam          = lambda wc: bam_path(wc.sample),
        global_panel = GLOBAL_PAN,
    output:
        png    = os.path.join(DATA_DIR, "{sample}", "results", "insert_size",
                              "{sample}_insert_size_panel.png"),
        lookup = os.path.join(DATA_DIR, "{sample}", "results", "insert_size",
                              "{sample}_insert_size_panel_lookup.json"),
    log:
        os.path.join(DATA_DIR, "{sample}", "results", "logs",
                     "{sample}_insert_size.log"),
    params:
        n_reads = config["insert_size"]["n_reads"],
        n_sites = config["insert_size"]["n_sites"],
    shell:
        """
        mkdir -p $(dirname {output.png})
        mkdir -p $(dirname {log})
        python {SCRIPTS}/validate_insert_size.py \
            --bam     {input.bam} \
            --panel   {input.global_panel} \
            --out     {output.png} \
            --n-reads {params.n_reads} \
            --n-sites {params.n_sites} \
            > {log} 2>&1
        """


# ── Rule 3: assign_reads (per sample, parallelizable) ────────────────────────
rule assign_reads:
    input:
        bam      = lambda wc: bam_path(wc.sample),
        barcodes = lambda wc: barcodes_path(wc.sample),
        panel    = BATCH_PAN,
        lookup   = lambda wc: insert_lookup(wc.sample),
    output:
        matrix  = os.path.join(DATA_DIR, "{sample}", "results", "counts",
                               "{sample}_matrix.npz"),
        obs     = os.path.join(DATA_DIR, "{sample}", "results", "counts",
                               "{sample}_obs.txt"),
        var     = os.path.join(DATA_DIR, "{sample}", "results", "counts",
                               "{sample}_var.txt"),
        summary = os.path.join(DATA_DIR, "{sample}", "results", "counts",
                               "{sample}_summary.json"),
    log:
        os.path.join(DATA_DIR, "{sample}", "results", "logs",
                     "{sample}_assign_reads.log"),
    params:
        out_base       = lambda wc: counts_base(wc.sample),
        indices_pkl_arg = f"--indices-pkl {INDICES_PKL}" if INDICES_PKL else "",
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        mkdir -p $(dirname {log})
        python {SCRIPTS}/assign_reads.py \
            --bam                {input.bam} \
            --panel              {input.panel} \
            --barcodes           {input.barcodes} \
            --insert-size-lookup {input.lookup} \
            --out                {params.out_base} \
            {params.indices_pkl_arg} \
            > {log} 2>&1
        """


# ── Rule 4: integrate_samples ────────────────────────────────────────────────
rule integrate_samples:
    input:
        matrices = [
            os.path.join(DATA_DIR, s, "results", "counts", f"{s}_matrix.npz")
            for s in SAMPLES
        ],
        h5_files    = [h5_path(s) for s in SAMPLES],
        batch_panel = BATCH_PAN,
    output:
        OUT_H5AD,
    log:
        os.path.join(BATCH_LOGS, "integrate.log"),
    params:
        samples   = " ".join(SAMPLES),
        data_dir  = DATA_DIR,
        suffix    = lambda wc: "{sample}/results/counts/{sample}",
        h5_paths  = " ".join(f"{s}={h5_path(s)}" for s in SAMPLES),
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        python {SCRIPTS}/integrate_samples.py \
            --samples   {params.samples} \
            --data_dir  {params.data_dir} \
            --iso_dir   {params.data_dir} \
            --suffix    '{params.suffix}' \
            --h5-paths  {params.h5_paths} \
            --panel     {input.batch_panel} \
            --out       {output} \
            > {log} 2>&1
        """
