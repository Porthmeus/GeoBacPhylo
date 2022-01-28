# Porthmeus
# 26.05.20

# denoise the samples using deblur

rule deblurQiime:
    input: "data/PE_demux_joined_trim.qza"
    output:
        FT = "data/PE_denoise_deblur_FT.qza",
        RS = "data/PE_denoise_deblur_RS.qza",
        stats = "reports/PE_denoise_deblur_Stats.qza",
        visStats = report("reports/PE_denoise_deblur_Stats.qzv",
                category = "Quality control",
                caption = "../captions/denoisingStats.rst")
    log: "logs/deblurQiime.log"
    conda: "../envs/qiime2.yaml"
    threads: 4
    shell:
        """
        qiime deblur denoise-16S \
                --i-demultiplexed-seqs {input}\
                --p-trim-length 295\
                --p-sample-stats \
                --p-jobs-to-start {threads} \
                --o-representative-sequences {output.RS} \
                --o-table {output.FT} \
                --o-stats {output.stats} &> {log}
        qiime deblur visualize-stats \
                --i-deblur-stats {output.stats} \
                --o-visualization {output.visStats} &> {log}
        """
