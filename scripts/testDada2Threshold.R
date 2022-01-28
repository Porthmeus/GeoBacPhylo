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


# construct the command for the qiime dada2
cmd <- "qiime dada2 denoise-paired"
args <- list(inp = "--i-demultiplexed-seqs",
    inp.val = snakemake@input[["cut"]],
    par1 = "--p-trunc-len-f",
    par1.val = snakemake@params[["trunc_f"]],
    par2 ="--p-trunc-len-r",
    par2.val = snakemake@params[["trunc_r"]],
    par3 = "--p-n-threads",
    par3.val = snakemake@threads[[1]],
    par4 = "--p-trim-left",
    par4.val = snakemake@params[["trim_front"]],
    par5 = "--p-trim-right",
    par5.val = snakemake@params[["trim_front"]],
    out1 = "--o-representative-sequences",
    out1.val = snakemake@output[["RS"]],
    out2 = "--o-table",
    out2.val = snakemake@output[["FT"]],
    out3 = "--o-denoising-stats",
    out3.val = snakemake@output[["stats"]],
    par6 = "--p-trunc-q",
    par6.val = 1)

# run deblur with different thresholds
thrl <- round(seq(10,30, length = 10))
dfs <- list()
for(i in thrl){
    args[["par6.val"]] <- i
    cm <- paste(cmd, paste(args, collapse=" "))
    system(command = cm)
    phy <- qza_to_phyloseq(features = snakemake@output[["FT"]], meta = snakemake@input[["meta"]])
    dfs[[as.character(i)]] <- estimate_richness(phy)
}

df <- do.call(rbind,dfs)
df[["thrl"]] <- sort(rep(thrl,nrow(df)))

write.csv(file = snakemake@output[["DF"]], df)
