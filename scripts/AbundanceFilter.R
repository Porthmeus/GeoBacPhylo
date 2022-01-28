# Porthmeus
# 30.07.20

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# use the abundance filters evaluated from the R markdown "EvaluationAbundanceFilter" to filter the OTU table

require(phyloseq)
require(data.table)
print(snakemake@input)

phy <- readRDS(snakemake@input[["physeq"]])

otu <- phy@otu_table

# calculate maximum contribution to a samples information and prevalence
depth <- apply(otu,2,sum)
maxFrac <- apply(otu, 1, function(x) max(x/depth, na.rm=TRUE))
prevalence <- apply(otu, 1, function(x) sum(x>0))

# filter the data for bacteria contribution
sel <- ((prevalence >3 | maxFrac >0.01))
otu2 <- otu[sel,]

# filter the data for sequencing depth
depth <- apply(otu2,2,sum)
sel <- scale(log(depth+1)) < -1*max(scale(log(depth+1)))
otu3 <- otu2[,!sel]

# remove all ASV, which have no sequence information in the leftover samples.
sel <- apply(otu3,1,sum) != 0
otu3 <- otu3[sel,]

sm <- colnames(otu3)
asv <- rownames(otu3)

tax <- phy@tax_table[asv,]
sm_dat <- phy@sam_data[sm,]

# write the files singular to disk
write.csv(file = snakemake@output[["FT"]], otu3)
write.csv(file = snakemake@output[["meta"]], sm_dat)
write.csv(file = snakemake@output[["tax"]], tax)

# write a fasta file from the leftover samples
cat(file = snakemake@output[["fasta"]], "", append = FALSE)

for(nm in asv){
    entry <- paste(">", nm,"\n", tax[nm,"rep.seq"], "\n", sep = "")
    cat(file = snakemake@output[["fasta"]], entry, append =TRUE)
}


