# Porthmeus
# 24.06.20

# test different quality thresholds for dada2 and its influence on the diversity, second step: read the data and calculate diversity measures

rule testDeblurThreshold_2:
    input:
        meta = "MetaData.csv",
        FT = expand("FilterThrl/{thr}_deblur_FT.qza", thr = trunc_thr)
    output:
        DF = "reports/deblur_ThrlOnAlphaDiv.csv"
    log: "logs/testDeblur2Threshold_2.log"
    conda: "../envs/R.yaml"
    script: "../scripts/testDada2Threshold_2.R"
