# Porthmeus
# 14.07.20

# get the alpha diverstiy measures of the different algorithms and thresholds and plot them to compare
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")



require(data.table)
require(ggplot2)
require(cowplot)

minMaxNorm <- function(x){
    (x-min(x))/(max(x)-min(x))
}

fls <- snakemake@input[["tbls"]]
#fls <- file.path("reports",list.files("reports", pattern = "OnAlphaDiv"))
tbls<-list()
for(f in fls){
    method <- strsplit(basename(f), split = "_")[[1]][1]
    #method <- basename(f)
    df <- fread(f)
    df[,rel_thrl := minMaxNorm(thrl)]
    tbls[[method]] <- cbind(df, Method= method)
}
df <- do.call(rbind, tbls)
suppressWarnings(df2<-melt(df, id = c(1,-2:0+ncol(df))))

# Plot
p <- ggplot(df2[variable %in% c("Observed","Shannon","Simpson"),],
            aes(x=factor(rel_thrl), y = value, fill = Method)) + 
    geom_boxplot() + 
    facet_grid(variable~., scale = "free_y") +
    xlab("Quality/Trimming threshold (min-max-normalized)")+
    theme_bw()


ggsave(p, file = snakemake@output[["plot"]], width = 8, height = 6, units = "in")
