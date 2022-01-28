# Porthmeus
# 20.07.20

# combine the grafics for the survivor for direct comparisons

# sink to log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# load packages
require(data.table)
require(ggplot2)
require(cowplot)

# small function to normalize the thresholds
minMaxNorm <- function(x){
    (x-min(x))/(max(x)-min(x))
}

#fls <- file.path("../reports/",list.files("../reports/", pattern = "thresholdSurvivors"))
fls <- snakemake@input[["survivors"]]
# load the tables, relativize threshold, add method column, remove unecessary columns
tbls <- list()
for(f in fls){
    mthd <- gsub(".csv","",strsplit(f,split="_")[[1]][2])
    df <- fread(f)
    df[,Method := mthd]
    df[,rel_thrl := minMaxNorm(threshold)]
    sel <- ncol(df) - (0:4)
    tbls[[mthd]] <- df[,..sel]
}
df <- do.call(rbind, tbls)

# create plots and save
p_abs <- ggplot(df, aes(x=as.factor(rel_thrl), y=survivors, fill = Method)) + 
    geom_boxplot() +
    theme_bw()+
    xlab("Quality/Trimming threshold (min-max-normalized)")+
    ylab("Number of reads surviving filtering")
p_rel <- ggplot(df, aes(x=as.factor(rel_thrl), y=fractionSurvived, fill = Method)) + 
    geom_boxplot() +
    theme_bw()+
    xlab("Quality/Trimming threshold (min-max-normalized)")+
    ylab("Fraction of reads surviving filtering")
p_all <- plot_grid(p_abs,p_rel, ncol=1)


ggsave(p_all, file = snakemake@output[["plot"]], width = 8, height = 6, units = "in")