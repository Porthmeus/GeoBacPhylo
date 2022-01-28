# Porthmeus
# 25.05.20

# denoising using uncombined reads with dada2

#trunc_thr = [300, 270, 240, 210, 180]
#qualities = [20, 15, 10, 5, 1]

rule dada2NoFiltPairedQiime:
    input: "data/PE_demux_cut.qza"
    output: 
        feat_table = "data/PE_denoise_dada2NoFilt_FT.qza",
        repr_seq = "data/PE_denoise_dada2NoFilt_RS.qza",
        stats = "reports/PE_denoise_dada2NoFilt_Stats.qza",
    log: "logs/dada2NoFiltPairedQiime.log"
    threads: 2
    params:
        trunc_f = 290,#200
        trunc_r = 235,#200
        trim_front= 9,
        trunc_q = 10,
        maxEE = "Inf"
    conda:"../envs/qiime2.yaml"
    shell:
        """
        qiime dada2 denoise-paired \
                --verbose \
                --i-demultiplexed-seqs {input} \
                --p-trunc-len-f {params.trunc_f} \
                --p-trunc-len-r {params.trunc_r} \
                --p-trunc-q {params.trunc_q} \
                --p-trim-left-f {params.trim_front} \
                --p-trim-left-r {params.trim_front} \
                --p-max-ee-f {params.maxEE} \
                --p-max-ee-r {params.maxEE} \
                --p-n-threads {threads} \
                --o-table {output.feat_table} \
                --o-representative-sequences {output.repr_seq} \
                --o-denoising-stats {output.stats} &> {log}
        """
