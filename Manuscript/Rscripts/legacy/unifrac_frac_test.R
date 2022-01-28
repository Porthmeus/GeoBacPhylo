
otu_frac <- t(t(otu)/apply(otu,2,sum))
phy_frac <- phyloseq(otu_table(otu_frac,
                              taxa_are_rows = TRUE),
                    tax_table(as.matrix(tax)),
                    sample_data(meta),
                    phy_tree(tree))
wunif <- phyloseq::distance(phy_frac, method = "unifrac")
lyotWunif <- umap::umap(as.matrix(wunif), distance = TRUE)$layout
lyotWunif <- data.table(lyotWunif, keep.rownames=TRUE)
test <- merge(lyotWunif, data.table(meta,keep.rownames=TRUE), by = "rn")
p <- ggplot(test, aes(x=V1,y=V2,color = Species)) + geom_point()
p
