# Porthmeus
# 22.05.20

# visualize the results of the trimming

rule SummarizeTrimmingQiime:
    input: "data/{infile}_trim.qza"
    output: "reports/{infile}_trim_summary.qzv"
    log: "logs/SummarizeTrimmingQiime_{infile}.log"
    conda: "../envs/qiime2.yaml"
    shell:
        "qiime demux summarize \
                --i-data {input} \
                --o-visualization {output} &> {log} "
