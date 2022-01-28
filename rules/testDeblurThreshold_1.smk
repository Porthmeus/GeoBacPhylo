# Porthmeus
# 23.06.20

# test different trimming length for deblur and its influence on the diversity
rule testDeblurThreshold:
    input:
        trim = "data/PE_demux_joined_trim.qza",
    output:
        RS = "temp/{thr}_deblur_RS.qza",
        FT = "temp/{thr}_deblur_FT.qza",
        stats= "temp/{thr}_deblur_stats.qza",
    threads: 4
    conda: "../envs/qiime2.yaml"
    log: "logs/{thr}_testDeblurThreshold.log"
    shell:
        """
        qiime deblur denoise-16S \
                --i-demultiplexed-seqs {input.trim}\
                --p-trim-length {wildcards.thr}\
                --p-sample-stats \
                --p-jobs-to-start {threads} \
                --o-representative-sequences {output.RS} \
                --o-table {output.FT} \
                --o-stats {output.stats} &> {log}
        """
