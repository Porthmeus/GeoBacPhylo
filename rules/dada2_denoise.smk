# Porthmeus
# 22.05.20

rule dada2_denoise:
    input: "data/PE_demux_trim.qza"
    output: 
        seqs = "data/PE_demux_trim_dada2.qza",
        data = "data/PE_dada2_table.qza"
    log: "logs/dada2_denoise.log"
    conda:"../envs/qiime2.yaml"
    shell:
        "qiime dada2 denoise-paired \
                -i-demultiplexed-seqs 
