# Porthmeus
# 23.06.20

# test different trimming length for deblur and its influence on the diversity
rule testDeblurSilvaThreshold:
    input:
        trim = "data/PE_demux_joined_trim.qza",
        ref = "misc/silva132_refSeqs.qza"
    output:
        RS = "FilterThrl/{thr}_deblurSilva_RS.qza",
        FT = "FilterThrl/{thr}_deblurSilva_FT.qza",
        stats= "FilterThrl/{thr}_deblurSilva_stats.qza",
    threads: 4
    conda: "../envs/qiime2.yaml"
    log: "logs/{thr}_testDeblurSilvaThreshold.log"
    shell:
        """
        qiime deblur denoise-other \
                --i-demultiplexed-seqs {input.trim}\
                --i-reference-seqs {input.ref} \
                --p-trim-length {wildcards.thr}\
                --p-sample-stats \
                --p-jobs-to-start {threads} \
                --o-representative-sequences {output.RS} \
                --o-table {output.FT} \
                --o-stats {output.stats} &> {log}
        """
