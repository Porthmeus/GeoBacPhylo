# Porthmeus
# 13.07.20

# test how many reads are filtered with different thrsholds

rule thresholdSurvivorsDeblur:
    input:
        stats = expand("FilterThrl/{thr}_deblur_stats.qza",
                thr = trunc_thr)
    output:
        csv = "reports/thresholdSurvivors_deblur.csv",
        plot = "plots/eval/thresholdSurvivors_deblur.svg"
    log: "logs/thresholdSurvivorsDeblur.log"
    conda: "../envs/R.yaml"
    script: "../scripts/thresholdSurvivors.R"
