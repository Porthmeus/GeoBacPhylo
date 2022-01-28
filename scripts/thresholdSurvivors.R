# Porthmeus
# 13.07.20

# compare the different thresholds and what the percent of filtered reads is
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R)
require(ggplot2)
require(cowplot)

#files <- file.path("temp",list.files("temp/", pattern = "deblur_stats.qza"))
for(fl in snakemake@input[["stats"]]){
    thr <- strsplit(basename(fl), split = "_")[[1]][1]
    df <- cbind(as.data.frame(qza_to_phyloseq(fl)@.Data), 
                threshold = thr)
    if(exists("outputDF")){
        outputDF <- rbind(outputDF, df)
    } else {
        outputDF <- df
    }
}

survivors <- c("merged","reads.hit.reference")
sel <- sapply(survivors, function(x) x %in% names(outputDF))
surv <- survivors[sel]
outputDF <- cbind(outputDF, survivors = outputDF[,surv], fractionSurvived = outputDF[,surv]/outputDF[,1])
write.csv(outputDF, file = snakemake@output[["csv"]], row.names =FALSE)

# plot the surviving reads
p_total <- ggplot(outputDF, aes(y=survivors, x = threshold)) + 
    geom_boxplot() +
    theme_bw() +
    ggtitle(strsplit(snakemake@input[["stats"]][1], split= "_")[[1]][2])
p_fraction <- ggplot(outputDF, aes(y = fractionSurvived, x = threshold)) + 
    geom_boxplot()+
    theme_bw() +
    ggtitle(strsplit(snakemake@input[["stats"]][1], split= "_")[[1]][2])
p_both <- plot_grid(p_total, p_fraction, nrow = 2)

ggsave(file = snakemake@output[["plot"]], p_both, width = 8, height=6, units = "in")
