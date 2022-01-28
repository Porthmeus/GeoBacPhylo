# Porthmeus
# 22.05.20

# create a aggregated report from the fastqc files


rule multiqc_raw:
    input: expand("fastqc/{sample}_fastqc.zip", sample = samples)
    output: 
        report("reports/multiqc_raw.html",
                caption = "../captions/multiqc_raw.rst",
                category = "Quality control")
    log: "logs/multiqc_raw.log"
    params: 
        outdir = "reports/"
    conda: "../envs/qiime2.yaml"
    shell: "multiqc {input} -n {output} &> {log}" 
