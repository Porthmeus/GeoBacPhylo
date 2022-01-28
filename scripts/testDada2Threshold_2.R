# Porthmeus
# 23.06.20

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# test different quality and its influence on the diversity
require(phyloseq)
if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R) # developmental package can be installed from github: https://github.com/jbisanz/qiime2R


ftables <- snakemake@input[["FT"]]
meta <- snakemake@input[["meta"]]
# run deblur with different thresholds
dfs <- list()
for(ft in ftables){
    print(ft)
    qual <- strsplit(basename(ft), split="_")[[1]][1]
    phy <- qza_to_phyloseq(features = ft, meta = meta)
    # print(qual)
    dfs[[qual]] <- cbind(estimate_richness(phy,
            measures=c("Observed","Shannon","Simpson")),
        thrl = as.integer(qual))
}

df <- do.call(rbind,dfs)

write.csv(file = snakemake@output[["DF"]], df)
