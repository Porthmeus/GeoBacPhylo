# Porthmeus
# 16.07.20

# test the differences in beta diversity within and outside a Population ID in dependence of the threshold used

fls = expand("FilterThrl/{ql}_{dada}_FT.qza", ql = qualities, dada = ["dada2","dada2NoFilt"])

fls.extend(expand("FilterThrl/{thr}_{deblur}_FT.qza", thr = trunc_thr, deblur = ["deblurSilva","deblur"]))

rule thesholBetaDiv:
    input:
        files = fls,
        meta = "MetaData.csv"
    output:
        plot = report("plots/eval/ThrlOnBetaDiv.jpg",
                caption = "../captions/thresholdBetaDiv.rst",
                category = "Filter thresholds"),
        plot_effsize = report("plots/eval/ThrlOnBetaDivEffSize.svg",
                caption = "../captions/thresholdBetaDiv_EffSize.rst",
                category = "Filter thresholds")
    log: "logs/thresholdBetaDiv.logs"
    conda: "../envs/R.yaml"
    script: "../scripts/thresholdBetaDiv.R"
