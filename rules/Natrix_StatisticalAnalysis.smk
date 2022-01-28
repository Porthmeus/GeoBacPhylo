# do the statistical analysis on the abundance data

rule Natrix_StatisticalAnalysis:
    input:
        AbFilt_FT = "data/dada2_AbFilt_FT.csv",
        AbFilt_Tax = "data/dada2_AbFilt_Tax.csv",
        AbFilt_meta = "data/dada2_AbFilt_meta.csv",
        AbFilt_fasta_tree = "data/dada2_AbFilt_fasta_tree.qza",
        MetaLakes = "MetaDataLakes.csv"
    output:
        html = report("Rmarkdown/Natrix_StatisticalAnalysis.html",
                caption = "../captions/Natrix_StatisticalAnalysis.rst",
                category = "Results")
    log: "logs/Natrix_StatisticalAnalysis.log"
    conda : "../envs/R.yaml"
    script:
        "../Rmarkdown/Natrix_StatisticalAnalysis.Rmd"

