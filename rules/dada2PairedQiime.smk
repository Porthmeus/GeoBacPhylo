# Porthmeus
# 25.05.20

# denoising using uncombined reads with dada2

trunc_thr = [300, 270, 240, 210, 180]
qualities = [20, 15, 10, 5, 1]


rule dada2PairedQiime:
    input: "data/PE_demux_cut.qza"
    output: 
        feat_table = expand ("data/PE_cut_dada2_f{fwd}_r{rev}_q{qual}_FT.qza",
                fwd = trunc_thr,
                rev = trunc_thr,
                qual = qualities),

        repr_seq = expand("data/PE_cut_dada2_f{fwd}_r{rev}_q{qual}_RS.qza",
                fwd = trunc_thr,
                rev = trunc_thr,
                qual = qualities),

        stats = expand("reports/PE_cut_dada2_f{fwd}_r{rev}_q{qual}_Stats.qza",
                fwd = trunc_thr,
                rev = trunc_thr,
                qual = qualities)
    log: "logs/dada2PairedQiime.log"
    threads: 2
    params:
        trunc_f = "{fwd}", #285,
        trunc_r = "{rev}", #269,
        trim_front= 9,
        trunc_q = "{qual}" #15
    conda:"../envs/qiime2.yaml"
    shell:
        "qiime dada2 denoise-paired \
                --verbose \
                --i-demultiplexed-seqs {input} \
                --p-trunc-len-f {params.trunc_f} \
                --p-trunc-len-r {params.trunc_r} \
                --p-trunc-q {params.trunc_q} \
                --p-trim-left-f {params.trim_front} \
                --p-trim-left-r {params.trim_front} \
                --p-n-threads {threads} \
                --o-table {output.feat_table} \
                --o-representative-sequences {output.repr_seq} \
                --o-denoising-stats {output.stats} &> {log}"
