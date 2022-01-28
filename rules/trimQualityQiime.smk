# Porthmeus
# 25.05.20

# Trim the sequences for quality, before assembly to ESVs

rule trimQualityQiime:
    input: "data/PE_demux_cut_joined.qza"
    output:
        data = "data/PE_demux_joined_trim.qza",
        stats= "reports/PE_demux_joined_trim_stats.qza"
    log: "logs/trimQualityQiime.log"
    shell:
        "qiime quality-filter q-score-joined \
                --i-demux {input} \
                --p-min-quality 20 \
                --verbose \
                --o-filtered-sequences {output.data} \
                --o-filter-stats {output.stats} &> {log}"
                #--p-min-length-fraction 0.8 \
