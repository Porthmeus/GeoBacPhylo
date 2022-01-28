# Porthmeus
# 22.05.20

# trim adapter sequences from the end of the raw reads

rule trimAdapterQiime:
    input: "data/PE_demux.qza"
    output: temp("data/PE_demux_cut.qza")
    log: "logs/trimAdapterQiime.log"
    threads: 4
    conda: "../envs/qiime2.yaml"
    shell:
        "qiime cutadapt trim-paired \
                --verbose \
                --p-cores {threads} \
                --p-adapter-f AGRGTTYGATYMTGGCTCAG \
                --p-adapter-r TGCTGCCTCCCGTAGGAGT \
                --p-match-read-wildcards \
                --p-match-adapter-wildcards \
                --p-minimum-length 20 \
                --i-demultiplexed-sequences {input} \
                --o-trimmed-sequences {output} &> {log}"
