# Porthmeus
# 13.07.20

# test how many reads are filtered with different thrsholds

rule thresholdSurvivorsDada2NoFilt:
    input:
        stats = expand("FilterThrl/{thr}_dada2NoFilt_stats.qza",
                thr = qualities)
    output:
        csv = "reports/thresholdSurvivors_dada2NoFilt.csv",
        plot = "plots/eval/thresholdSurvivors_dada2NoFilt.svg"
    log: "logs/thresholdSurvivorsDada2NoFilt.log"
    conda: "../envs/R.yaml"
    script: "../scripts/thresholdSurvivors.R"
