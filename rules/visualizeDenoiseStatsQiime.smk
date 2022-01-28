# Porthmeus
# 26.05.20

rule visualizeDenoiseStatsQiime:
    input: "reports/PE_denoise_{dn}_Stats.qza"
    output: 
        report("reports/PE_denoise_{dn}_Stats.qzv",
                caption = "../captions/denoisingStats.rst",
                category = "Quality control")
    log: "logs/visualizeDenoiseStatsQiime_{dn}.log"
    conda: "../envs/qiime2.yaml"
    shell:
        "qiime metadata tabulate \
                --m-input-file {input} \
                --o-visualization {output} &>{log}"
