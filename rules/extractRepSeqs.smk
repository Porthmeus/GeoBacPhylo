# Porthmeus
# 26.05.20

# extract the representing sequences of the features

rule extractRepSeqs:
    input: "data/PE_denoise_{dn}_RS.qza"
    output: "reports/PE_denoise_{dn}_seqTab.qzv"
    log: "logs/extractRepSeqs_{dn}.log"
    conda: "../envs/qiime2.yaml"
    shell:
        "qiime feature-table tabulate-seqs \
                --i-data {input} \
                --o-visualization {output} &> {log}"
