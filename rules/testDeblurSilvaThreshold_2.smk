# Porthmeus
# 24.06.20

# test different quality thresholds for dada2 and its influence on the diversity, second step: read the data and calculate diversity measures

rule testDeblurSilvaThreshold_2:
    input:
        meta = "MetaData.csv",
        FT = expand("FilterThrl/{thr}_deblurSilva_FT.qza", thr = trunc_thr)
    output:
        DF = "reports/deblurSilva_ThrlOnAlphaDiv.csv"
    log: "logs/testDeblurSilva2Threshold_2.log"
    conda: "../envs/R.yaml"
    script: "../scripts/testDada2Threshold_2.R"
