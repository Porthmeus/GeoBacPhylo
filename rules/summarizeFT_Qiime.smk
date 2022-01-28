# Porthmeus
# 26.05.20

# get the visualization object for feature table

rule summarizeFT_Qiime:
    input:
        FT = "data/PE_denoise_{dn}_FT.qza",
        meta = "MetaData.csv"
    output: 
        report("reports/PE_denoise_{dn}_sumFT.qzv",
                category = "Quality control",
                caption = "../captions/PE_denoise_sumFT.rst")
    log: "logs/summarizeDenoisedFTQiime_{dn}.log"
    conda: "../envs/qiime2.yaml"
    shell:
        "qiime feature-table summarize \
                --i-table {input.FT} \
                --m-sample-metadata-file {input.meta} \
                --o-visualization {output} &> {log}"
