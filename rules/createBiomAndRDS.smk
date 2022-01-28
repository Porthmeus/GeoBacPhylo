# Porthmeus
# 23.06.20

# convolute the data from qiime into a phyloseq object and save it as biom and RDS file
# note that conda is currently not supported as the function depends on a package called
rule createBiomAndRDS:
    input:
        FT = "data/PE_denoise_{dn}_FT.qza",
        tree = "data/PE_denoise_{dn}_tree.qza",
        meta = "MetaData.csv",
        tax = "data/PE_denoise_{dn}_Tax.qza",
        RS = "data/PE_denoise_{dn}_RS.qza"
    output:
        RDS = report("data/PE_denoise_{dn}_physeq.RDS",
                category = "Data",
                caption = "../captions/physeqRDS.rst"),
        biom = report("data/PE_denoise_{dn}_biom.biom",
                category = "Data",
                caption = "../captions/biom.rst")
    log: "logs/createBiomAndRDS_{dn}.log"
    conda: "../envs/R.yaml"
    script: "../scripts/createBiomAndRDS.R"

