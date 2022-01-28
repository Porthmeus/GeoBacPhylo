# Porthmeus
# 25.05.20

# create some plots for measuring the effect from the dada2 filtering

rule dada2ThresholdPlots:
    input: "data/PE_cut_dada2_f{fwd}_r{rev}_q{qual}_FT.qza"
    output: 
        svg = "reports/dada2Thresholds.svg",
        dataframe = "reports/dada2Tresholds.csv"
    log: "logs/dada2ThresholdPlots.log"
    conda: "../envs/qiime2.yaml"
    script:
        "scripts/dada2TresholdPlots.R"
