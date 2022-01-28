# Porthmeus
# 14.07.20

# plot the dependence of alpha diversity
rule thresholdAlphaDiv:
    input:
        tbls = expand("reports/{dn}_ThrlOnAlphaDiv.csv",
                dn = dn)
    output:
        plot = report("plots/eval/ThrlOnAlphaDiv.svg",
                caption = "../captions/thresholdAlphaDiv.rst",
                category = "Filter thresholds")
    conda: "../envs/R.yaml"
    log: "logs/thresholdAlphaDiv.logs"
    script: "../scripts/thresholdAlphaDiv.R"
