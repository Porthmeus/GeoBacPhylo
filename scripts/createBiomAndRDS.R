# Porthmeus
# 23.06.20

# import the different outputs from qiime2 and generate a phyloseq data set from it, as well as write a biom file

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(phyloseq)
if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R) # developmental package can be installed from github: https://github.com/jbisanz/qiime2R
require(biomformat)

# read the data
physeq <- qza_to_phyloseq(features = snakemake@input[["FT"]],tree = snakemake@input[["tree"]], metadata = snakemake@input[["meta"]])
# add the taxonomy and sequence information
# read tax and reformat
tax <- read_qza(snakemake@input[["tax"]])
colnames(tax$data)[3] <- "Confidence"
tax <- (parse_taxonomy(tax$data))

# read representative sequences for ASV
seq <- read_qza(snakemake@input[["RS"]])
seq <- as.data.frame(seq$data)
colnames(seq)[1] <- "rep.seq"

# merge and reformat
st <- merge( tax, seq, by = 0)
rownames(st) <- st[["Row.names"]] 
st <- st[,-1]

# add to phyloseq object
physeq <- merge_phyloseq(physeq, tax_table(as.matrix(st)))

# save the file as rds data for easy import in R
saveRDS(object = physeq, file = snakemake@output[["RDS"]])

# write the dataset to biom format (without tree information though)
biom <- make_biom(data=physeq@otu_table, sample_metadata = physeq@sam_data, observation_metadata = physeq@tax_table)
write_biom(x = biom, biom_file = snakemake@output[["biom"]])
