# Porthmeus
# 30.07.20

# check the abundance filter applied to the data

rule Natrix_EvaluationAbundanceFilter:
    input:
        physeq = "data/PE_denoise_dada2_physeq.RDS"
    output:
        html = report("Rmarkdown/Natrix_EvaluationAbundanceFilter.html",
                caption = "../captions/Natrix_EvaluationAbundanceFilter.rst",
                category = "Filter thresholds")
    log: "logs/Natrix_EvaluationAbundanceFilter.log"
    conda: "../envs/R.yaml"
    script:
        "../Rmarkdown/Natrix_EvaluationAbundanceFilter.Rmd"
    
    
    # if one wants to generate several output files, one have to hardcode them
    # within the Rmarkdown script and render the html using the following
    # command. Outputs can afterwards be defined in the output section of the
    # snakemake rule to have snakemake check whether all outputs where
    # generated.
    #shell:
    #    '''R -e "rmarkdown::render('Rmarkdown/Natrix_EvaluationAbundanceFilter.Rmd')" 2&> {log}'''
