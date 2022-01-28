# Porthmeus
# 23.10.19
rule fastqc_raw:
    input:
        R1 = "raw/{sample}_R1_001.fastq.gz",
        R2 = "raw/{sample}_R2_001.fastq.gz"
    output:
        R1 = "fastqc/{sample}_R1_001_fastqc.zip",
        R2 = "fastqc/{sample}_R2_001_fastqc.zip"
    conda:
        "../envs/qiime2.yaml"
    params:
        outdir = "fastqc"
    threads: 4
    log:
         "logs/{sample}_fastqc_raw.log"
    shell:
        "fastqc {input} -o {params.outdir} -t {threads} &> {log}"
