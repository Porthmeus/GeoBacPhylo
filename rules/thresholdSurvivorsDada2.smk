# Porthmeus
# 13.07.20

# test how many reads are filtered with different thrsholds

rule thresholdSurvivorsDada2:
    input:
        stats = expand("FilterThrl/{thr}_dada2_stats.qza",
                thr = qualities)
    output:
        csv = "reports/thresholdSurvivors_dada2.csv",
        plot = "plots/eval/thresholdSurvivors_dada2.svg"
    log: "logs/thresholdSurvivorsDada2.log"
    conda: "../envs/R.yaml"
    script: "../scripts/thresholdSurvivors.R"
