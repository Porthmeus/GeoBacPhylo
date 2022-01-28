# Porthmeus
# 23.06.20


# test different quality thresholds for dada2 and its influence on the diversity, first step: filter and denoise
rule testDada2Threshold_1:
    input:
        trim = "data/PE_demux_cut.qza",
    output:
        RS = "FilterThrl/{q}_dada2_RS.qza",
        FT = "FilterThrl/{q}_dada2_FT.qza",
        stats= "FilterThrl/{q}_dada2_stats.qza",
    params:
        trunc_f = 290,#200
        trunc_r = 235,#200
        trim_front= 9
    threads: 4
    conda: "../envs/qiime2.yaml"
    log: "logs/{q}_Dada2Threshold_1.log"
    shell:
        """
        qiime dada2 denoise-paired \
                --verbose \
                --i-demultiplexed-seqs {input.trim} \
                --p-trunc-len-f {params.trunc_f} \
                --p-trunc-len-r {params.trunc_r} \
                --p-trunc-q {wildcards.q} \
                --p-trim-left-f {params.trim_front} \
                --p-trim-left-r {params.trim_front} \
                --p-n-threads {threads} \
                --o-table {output.FT} \
                --o-representative-sequences {output.RS} \
                --o-denoising-stats {output.stats} &> {log}
        """
