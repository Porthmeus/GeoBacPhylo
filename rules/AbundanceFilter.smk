# Porthmeus
# 30.07.20

# filter the otu table for abundance and importance of samples and ASVs

rule AbundanceFilter:
    input:
        physeq = "data/PE_denoise_dada2_physeq.RDS"
    output:
        FT = "data/dada2_AbFilt_FT.csv",
        meta = "data/dada2_AbFilt_meta.csv",
        tax = "data/dada2_AbFilt_Tax.csv",
        fasta = "data/dada2_AbFilt_fasta.fasta"
    conda: "../envs/R.yaml"
    log: "logs/AbundanceFilter.log"
    script: "../scripts/AbundanceFilter.R"

