# Porthmeus
# 14.07.20
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R)
require(phyloseq)
require(vegan)
require(data.table)
require(ggplot2)
require(effsize)

# write a small function to form the data frame needed for plotting the beta diversity
withinOutside <- function(otuTable, meta){
    dstnc <- vegdist(t(otuTable), method = "morisita")
    dstnc <- data.table(reshape2::melt(as.matrix(dstnc)))
    df <- merge(dstnc, meta[,.(ID,PopID)], by.x = "Var1", by.y="ID") 
    df <- merge(df, meta[,.(ID,PopID)], by.x = "Var2", by.y="ID")
    sel <- df[["PopID.x"]] == df[["PopID.y"]]
    df[["compareTo"]] <- "outside"
    df[sel,"compareTo"] <- "within"
    return(df)
}

# small function to normalize the thresholds
minMaxNorm <- function(x){
    (x-min(x))/(max(x)-min(x))
}

meta <- fread(snakemake@input[["meta"]])
fls <- snakemake@input[["files"]] # file.path("temp",list.files("temp",pattern ="_FT.qza"))

tbls <- list()

for(fl in fls){
    print(fl)
    qualAlg <- strsplit(basename(fl), split="_")[[1]][1:2]
    df <- read_qza(fl)
    df <- cbind(withinOutside(df$data, meta), thrl = as.integer(qualAlg[1]), Method = qualAlg[2])
    tbls[[paste(qualAlg, collapse = "_")]] <- df
}
df <- do.call(rbind,tbls)
df[,rel_thrl := minMaxNorm(thrl), by=Method]

p <- ggplot(df, aes(fill=compareTo, y=value, x =factor(rel_thrl))) +
    geom_boxplot() +
    facet_grid(Method~.) +
    theme_bw() +
    ylab("Morisita score") +
    xlab("Quality/Trimming threshold (min-max-normalized)")
    
ggsave(p, file = snakemake@output[["plot"]], width = 8, height = 8, units="in")

# calculate the effect size to get a better impression

d_hedge <- lapply(tbls, function(x) c(unlist(cohen.d(value~compareTo, x, hedges.correction=T,na.rm=T)[c("estimate","conf.int")])))
df_hedge <- do.call(rbind, d_hedge)
df_hedge <- cbind(df_hedge, 
    thrl = as.integer(sapply(rownames(df_hedge), function(x) strsplit(x, split = "_")[[1]][1])), 
    Method = sapply(rownames(df_hedge), function(x) strsplit(x, split = "_")[[1]][2])
    )
df_hedge <- data.table(df_hedge)
df_hedge[, rel_thrl := minMaxNorm(as.integer(thrl)), by = Method]

p_hedge <- ggplot(df_hedge, aes(x=rel_thrl, y=as.numeric(estimate), ymin=as.numeric(conf.int.lower), ymax=as.numeric(conf.int.upper), color=Method, fill = Method)) +
    geom_line() +
    geom_point() +
    geom_ribbon(alpha=0.5)

ggsave(p_hedge, file = snakemake@output[["plot_effsize"]], width = 8, height = 4, units = "in")
