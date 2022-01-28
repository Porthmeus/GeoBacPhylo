# Porthmeus
# 24.06.20

# test different quality thresholds for dada2 and its influence on the diversity, second step: read the data and calculate diversity measures

rule testDada2Threshold_2:
    input:
        meta = "MetaData.csv",
        FT = expand("FilterThrl/{q}_dada2_FT.qza", q = qualities)
    output:
        DF = "reports/dada2_ThrlOnAlphaDiv.csv"
    log: "logs/testDada2Threshold_2.log"
    conda: "../envs/R.yaml"
    script: "../scripts/testDada2Threshold_2.R"
