# Porthmeus
# 23.06.20

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# test different trimming length and its influence on the diversity
require(phyloseq)
if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R) # developmental package can be installed from github: https://github.com/jbisanz/qiime2R
print(snakemake@input)
print(snakemake@output)
system("qiime tools view --help")
print(Sys.getenv("PATH"))
# construct the command for the qiime deblur
cmd <- "qiime deblur denoise-16S"
args <- list(inp = "--i-demultiplexed-seqs",
    inp.val = snakemake@input[["trim"]],
    par1 = "--p-sample-stats",
    par2 ="--p-jobs-to-start",
    par.val = snakemake@threads[[1]],
    out1 = "--o-representative-sequences",
    out1.val = snakemake@output[["RS"]],
    out2 = "--o-table",
    out2.val = snakemake@output[["FT"]],
    out3 = "--o-stats",
    out3.val = snakemake@output[["stats"]],
    par3 = "--p-trim-length",
    par3.val = -1)

# run deblur with different thresholds
thrl <- round(seq(100,300, length = 10))
dfs <- list()
for(i in thrl){
    args[["par3.val"]] <- i
    cm <- paste(cmd, paste(args, collapse=" "))
    system(command = cm)
    phy <- qza_to_phyloseq(features = snakemake@output[["FT"]], meta = snakemake@input[["meta"]])
    dfs[[as.character(i)]] <- estimate_richness(phy)
}

df <- do.call(rbind,dfs)
df[["thrl"]] <- sort(rep(thrl,nrow(df)))

write.csv(file = snakemake@output[["DF"]], df)
