# Porthmeus
#11.02.21

require(ggplot2)
require(phyloseq)
require(qiime2R)
require(data.table)
require(jtools)
require(cowplot)

otu <- read.csv("../data/dada2_AbFilt_FT.csv", row.names = 1)
tax <- read.csv("../data/dada2_AbFilt_Tax.csv", row.names =1)
meta <- read.csv("../data/dada2_AbFilt_meta.csv", row.names=1)
tree <- read_qza("../data/dada2_AbFilt_fasta_tree.qza")$data

map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]

# create a phyloseq object
phy <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax)),
                sample_data(meta),
                phy_tree(tree))

# The idea of Jacint was to associate the alpha diversity to the 
alphaDiv <- estimate_richness(phy)
alphaNames <- colnames(alphaDiv)
alphaDiv <- merge(alphaDiv, meta, by=0, all.x=TRUE)
# transform the data to fit normality a little better
alphaDiv[["Simpson"]] <- alphaDiv[["Simpson"]]
alphaDiv[["Chao1"]] <- (alphaDiv[["Chao1"]])


alphaPlots <- list()
for(nm in c("Simpson","Shannon","Chao1")){
    for(group in c("PopID", "Species", "Reproduction")){
        id <- paste(nm, group, sep ="_")
        p <- ggplot(alphaDiv, aes_string(x = group, y= nm, fill = group))+
                    geom_boxplot() +
                    ylab(nm) +
                    xlab(group) +
                    theme_bw() +
                    theme(legend.position = "None",
                          axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) 
        alphaPlots[[id]] <- p
    }
}
pp <- cowplot::plot_grid(plotlist = alphaPlots, rel_widths = rep(c(3,1,1),3))

ggsave(pp, file ="figures/Boxplot_alphaDiv.pdf", width = 7, height = 6)


simpMod <- lm(data=alphaDiv, Simpson ~ PopID+Species)
shanMod <- lm(data=alphaDiv, Shannon ~ PopID+Species+Reproduction)
chaoMod <- lm(data=alphaDiv, Chao1 ~ PopID+Species)

# export the anova and summary tables of the models
simpModTable <- summary(simpMod)
write.csv(file = "tables/Simpson_coefficientStats.csv", simpModTable[["coefficients"]])
simpModAnova <- anova(simpMod)
write.csv(file = "tables/Simpson_Anova.csv", as.data.frame(simpModAnova))

shanModTable <- summary(shanMod)
write.csv(file = "tables/Shannon_coefficientStats.csv", shanModTable[["coefficients"]])
shanModAnova <- anova(shanMod)
write.csv(file = "tables/Shannon_Anova.csv", as.data.frame(shanModAnova))

chaoModTable <- summary(chaoMod)
write.csv(file = "tables/Chao_coefficientStats.csv", chaoModTable[["coefficients"]])
chaoModAnova <- anova(chaoMod)
write.csv(file = "tables/Chao_Anova.csv", as.data.frame(chaoModAnova))


models <- list("Simpson" = simpMod ,"Shannon" = shanMod,"Chao1" = chaoMod)
effecSizePlots <- list()
for(n in names(models)){
    for(vrb in c("PopID","Species")){
    mod  <- models[[n]]
    prd  <- make_predictions(mod, pred = vrb, data = alphaDiv)
    effecSizePlots[[paste(n,vrb,sep="_")]] <-
        ggplot(prd, aes_string(x = vrb, y=n ,ymin = "ymin", ymax = "ymax", color = vrb)) +
        geom_point() +
        geom_errorbar() +
        theme_bw() +
        theme(legend.position = "None",
              axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) 
    }
}

pp_effect <- plot_grid(plotlist=effecSizePlots, ncol = 2, rel_widths = rep(c(3,1),3))
ggsave(pp_effect, file = "figures/Effect_alphaDiv.pdf", width = 6, height = 5)
