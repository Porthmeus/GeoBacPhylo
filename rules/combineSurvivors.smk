# Porthmeus
# 20.07.20

# combine the surviving plots

rule combineSurvivors:
    input: 
        survivors = expand("reports/thresholdSurvivors_{dn}.csv",
                dn = dn)
    output:
        plot = report("plots/eval/thresholdSurvivors.svg",
                category = "Filter thresholds",
                caption = "../captions/thresholdSurvivors.rst")
    conda: "../envs/R.yaml"
    log: "logs/combineSurvivors.log"
    script: "../scripts/combineSurvivors.R"
