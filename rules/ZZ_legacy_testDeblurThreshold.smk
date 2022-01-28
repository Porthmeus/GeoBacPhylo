# Porthmeus
# 23.06.20

# test different trimming length for deblur and its influence on the diversity
rule testDeblurThreshold:
    input:
        trim = "data/PE_demux_joined_trim.qza",
        meta = "MetaData.csv"
    output:
        RS = temp("temp/temp_deblur_RS.qza"),
        FT = temp("temp/temp_deblur_FT.qza"),
        stats= temp("temp/temp_deblur_stats.qza"),
        DF = "reports/deblurTrimOnAlphaDiv.csv"
    threads: 4
    conda: "../envs/qiime2.yaml"
    log: "logs/testDeblurThreshold.log"
    script: "../scripts/testDeblurThreshold.R"

