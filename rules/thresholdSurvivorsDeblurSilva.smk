# Porthmeus
# 13.07.20

# test how many reads are filtered with different thrsholds

rule thresholdSurvivorsDeblurSilva:
    input:
        stats = expand("FilterThrl/{thr}_deblurSilva_stats.qza",
                thr = trunc_thr)
    output:
        csv = "reports/thresholdSurvivors_deblurSilva.csv",
        plot = "plots/eval/thresholdSurvivors_deblurSilva.svg"
    log: "logs/thresholdSurvivorsDeblurSilva.log"
    conda: "../envs/R.yaml"
    script: "../scripts/thresholdSurvivors.R"
