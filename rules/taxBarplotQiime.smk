# Porthmeus
# 07.06.20

# create an interactive barplot

rule taxBarplotQiime:
    input:
        FT = "data/PE_denoise_{dn}_FT.qza",
        tax = "data/PE_denoise_{dn}_Tax.qza",
        meta = "MetaData.csv"
    output:
        bar = report("results/taxa_bar_{dn}.qzv",
                caption = "../captions/taxBarplot.rst",
                category = "Taxonomy")
    log: "logs/taxBarplotQiime_{dn}.log"
    conda: "../envs/qiime2.yaml"
    shell:
        "qiime taxa barplot \
                --i-table {input.FT} \
                --i-taxonomy {input.tax} \
                --m-metadata-file {input.meta} \
                --o-visualization {output.bar} &> {log}"
