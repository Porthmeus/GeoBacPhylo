# Porthmeus 
# 14.12.21

# analysis of the Population M90 and M89 only to test for species specific changes in one sampling site

require(data.table)
require(DESeq2)
if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R)
require(ggplot2)
require(vegan)
require(phyloseq)
require(lme4)
require(lmerTest)
require(car)
require(lsmeans)
require(multcomp)
require(umap)

# load data

AbFilt_FT = "../data/dada2_AbFilt_FT.csv"
AbFilt_Tax = "../data/dada2_AbFilt_Tax.csv"
AbFilt_meta = "../data/dada2_AbFilt_meta.csv"
AbFilt_fasta_tree = "../data/dada2_AbFilt_fasta_tree.qza"
MetaLakes = "../MetaDataLakes.csv"


otu <- read.csv(AbFilt_FT, row.names = 1)
# remove the smallest 5% (<3368 reads) from the samples - ~14
sumRC <- apply(otu,2,sum)
otu <- otu[,sumRC > quantile(sumRC,0.05)]
tax <- read.csv(AbFilt_Tax, row.names=1)
meta <- read.csv(AbFilt_meta, row.names=1)
metaLakes <- read.csv(MetaLakes,row.names=1)

# simplify the reproduction vector
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]

# add the correct naming for the nutritional state of the water bodies
map_nut <- c("Hypereutr.", "Eutr.","Mesoeutr.")
meta[["Nutrient_load"]] <- map_nut[metaLakes[meta$PopID,"NutrientLoad"]]
meta[["Water_body2"]] <- metaLakes[meta$PopID, "Waterbody2"]
meta[["Water_body"]] <- metaLakes[meta$PopID, "waterbody"]
meta[["Water_body"]] <- gsub("River","running",gsub("Lake","standing",meta[["Water_body"]]))
meta[["LocationIDshort"]] <- sapply(meta$LocationID, function(x) strsplit(x, split ="/")[[1]][2])

# reorder the sampling sites, so that they occur in a similar order as no the map
siteOrder <- c("M34","M72","M71","M70","M31","M67","M89","M90","M86","M85","M84","M79","M78","M109","M28","M26","M44","M108","M107","M83","R12")
meta[["SiteID"]] <- factor(meta[["PopID"]], levels = siteOrder)

# remove rare species (lower 25% of the data)
countsPerESV <- apply(otu,1,sum)
thrs <- quantile(countsPerESV,0.25)
summary(countsPerESV)/sum(otu)
plot(ecdf(log10(countsPerESV+1)))
abline(v = log10(thrs))
hist(log10(countsPerESV))
abline(v=log10(thrs))

ESVs <- names(countsPerESV)[countsPerESV > thrs]
otu <- otu[ESVs,]
tax <- tax[ESVs,]

# remove chloroplast ESVs from (probably) algea
sel <- which(tax[["Order"]] != "Chloroplast")
tax <- tax[sel,]
otu <- otu[rownames(tax),]

# create phyloseq object
phy <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax)),
                sample_data(meta))

alpha <- estimate_richness(phy,measures = "Shannon")
alpha <- merge(alpha, meta, by=0)
alpha[["Nutrient_load"]] <- relevel(as.factor(alpha[["Nutrient_load"]]), "Mesoeutr.")
alpha <- data.table(alpha)

#### ALPHA DIVERSITY ####
# select only M90 and M89
alpha <- alpha[PopID %in% c("M89","M90"),]

alpha_Species_plots <- list()
alpha_Reproduction_plots <- list()
alpha_stats <- list()
for(pop in unique(alpha[,PopID])){
    alphaPop <- alpha[PopID == pop,]
# plot alpha diversity on species
    p_alphaSpecies <- ggplot(alphaPop, aes(x = Species, y = Shannon, fill = Species))+
        geom_boxplot() +
        #ggtitle("Alpha diversity") +
        theme_bw() +
        theme(legend.position = "None")


    p_alphaReproduction <- ggplot(alphaPop, aes(x = Reproduction, y = Shannon, fill = Reproduction))+
        geom_boxplot() +
        theme_bw()+
        theme(legend.position = "None")
    alpha_Species_plots[[pop]] <- p_alphaSpecies
    alpha_Reproduction_plots[[pop]] <- p_alphaReproduction
    
    alpha_stats[[paste(pop,"Species","Reproduction", sep = "_")]] <- lm(data = alphaPop, Shannon ~ Species + Reproduction)
    alpha_stats[[paste(pop,"Species", sep = "_")]] <- lm(data = alphaPop, Shannon ~ Species)
    alpha_stats[[paste(pop,"Reproduction", sep = "_")]] <- lm(data = alphaPop, Shannon ~ Reproduction)
}


#### BETA DIVERSITY ####

# add a single count to the OTU matrix, stabilize variance, and shift the values to positive values
set.seed(5245)
otu2 <- otu+1
otu_vst <- varianceStabilizingTransformation(as.matrix(otu2), fitType = "local")
if(min(otu_vst) <0){
    otu_vst <- otu_vst + min(otu_vst)*-1
}

phy_vst <- phyloseq(otu_table(otu_vst,
                              taxa_are_rows = TRUE),
                              tax_table(as.matrix(tax)),
                              sample_data(meta))


# calculate bray curtis
brayCurt <- distance(phy_vst, method = "bray")
smpls <- rownames(as.matrix(brayCurt))
meta <- meta[smpls,]
meta <- data.table(meta, keep.rownames=TRUE)
metaM9089 <- meta[SiteID %in% c("M90","M89"),]

stats <- list()
umap_plots <- list()
pca_plots <- list()
dist_Species_plots <- list()
dist_Reproduction_plots <- list()
for(pop in unique(metaM9089[,SiteID])){
    metaPop <- metaM9089[SiteID == pop,]
    bcPop <- as.dist(as.matrix(brayCurt)[metaPop[,rn],metaPop[,rn]])
    # do the stats
    if(length(unique(metaPop[,Reproduction])) > 1){
        stats[[paste0(pop, "_Species_Reproduction")]]  <- adonis2(data = metaPop, formula = bcPop ~ Species + Reproduction)
        stats[[paste0(pop, "_Reproduction_Species")]] <- adonis2(data = metaPop, formula = bcPop ~  Reproduction + Species)
        stats[[paste0(pop, "_Reproduction")]] <- adonis2(data = metaPop, formula = bcPop ~ Reproduction)
    }
    stats[[paste0(pop, "_Species")]] <- adonis2(data = metaPop, formula = bcPop ~ Species)

    # save the data for the plots
    # start with the UMAP plots
    umap_dat <- umap(as.matrix(bcPop), input = "dist")$layout
    colnames(umap_dat) <- c("UMAP1","UMAP2")
    umap_dat <- data.table(umap_dat, keep.rownames = TRUE)
    umap_dat <- merge(umap_dat, metaPop, by = "rn")
    umap_p <- ggplot(umap_dat, aes(x=UMAP1, y = UMAP2, color = Species, shape = Reproduction)) +
        geom_point(size = 2) +
        #ggtitle("UMAP") +
        theme_bw()
    umap_plots[[pop]] <- umap_p

    # go on with the PCA
    pca_dat <- prcomp(bcPop, scale = TRUE, center = TRUE)
    pca_var <- (pca_dat$sdev)^2/sum(pca_dat$sdev^2)
    pca_dat <- data.table(pca_dat[["x"]][,1:2], keep.rownames = TRUE)
    pca_dat <- merge(pca_dat, metaPop, by = "rn")
    pca_p <- ggplot(pca_dat, aes( x= PC1, y = PC2, color = Species, shape = Reproduction)) +
        geom_point(size = 2) + 
        #ggtitle("PCA") +
        ylab(paste0("PC2 (",round(pca_var[2]*100,2),"%)")) +
        xlab(paste0("PC1 (",round(pca_var[1]*100,2),"%)")) +
        theme_bw()
    pca_p

    pca_plots[[pop]] <- pca_p
    
    # add the distance plots for Species
    brayCurt_long <- as.vector(bcPop)
    n <- as.data.frame(t(combn(rownames(as.matrix(bcPop)),2)))
    brayCurt_long <- cbind(n,brayCurt = brayCurt_long)
    colnames(brayCurt_long) <- c("Sample1","Sample2","Bray_Curtis")
    brayCurt_long <- data.table(brayCurt_long)
    setkey(metaPop, "rn")

    # add the information if within or between population IDs
    brayCurt_long[,Species := c("within","between")[as.integer(metaPop[Sample1, "Species"] != metaPop[Sample2, "Species"])+1]]
    brayCurt_long[,Species := factor(Species, levels = c("within","between"))]


    Species_long <- paste(metaPop[brayCurt_long[,Sample1],Species],metaPop[brayCurt_long[,Sample2], Species], sep = "-")
    n <- combn(sort(unique(metaPop[["Species"]])),2)
    crct <- apply(n,2,paste,collapse="-")
    wrng <- apply(n[c(2,1),],2,paste,collapse="-")
    for(i in 1:ncol(n)){
        Species_long <- gsub(wrng[i],crct[i],Species_long)
    }
    brayCurt_long[,Species_long:=factor(Species_long, levels = c("CIR-CIR","VUL-VUL","OLI-OLI","CIR-VUL","CIR-OLI","OLI-VUL"))]



    p_boxDist <- ggplot(data = brayCurt_long, aes(x=Species_long, y=Bray_Curtis)) +
        geom_boxplot() +
        ylab("Bray-Curtis distance") +
        theme_bw() +
        xlab("") +
        #ggtitle("Beta Diversity")+
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
    dist_Species_plots[[pop]] <- p_boxDist

    # add the distance plots for Reproduction
    brayCurt_long <- as.vector(bcPop)
    n <- as.data.frame(t(combn(rownames(as.matrix(bcPop)),2)))
    brayCurt_long <- cbind(n,brayCurt = brayCurt_long)
    colnames(brayCurt_long) <- c("Sample1","Sample2","Bray_Curtis")
    brayCurt_long <- data.table(brayCurt_long)
    setkey(metaPop, "rn")

    # add the information if within or between population IDs
    brayCurt_long[,Reproduction := c("within","between")[as.integer(metaPop[Sample1, "Reproduction"] != metaPop[Sample2, "Reproduction"])+1]]
    brayCurt_long[,Reproduction := factor(Reproduction, levels = c("within","between"))]
    Reproduction_long <- paste(metaPop[brayCurt_long[,Sample1],Reproduction],metaPop[brayCurt_long[,Sample2], Reproduction], sep = "-")
    n <- combn(sort(unique(metaPop[["Reproduction"]])),2)
    crct <- apply(n,2,paste,collapse="-")
    wrng <- apply(n[c(2,1),],2,paste,collapse="-")
    for(i in 1:ncol(n)){
        Reproduction_long <- gsub(paste0("^",wrng[i],"$"),crct[i],Reproduction_long)
    }
    brayCurt_long[,Reproduction_long:=factor(Reproduction_long)]


    p_boxDist <- ggplot(data = brayCurt_long, aes(x=Reproduction_long, y=Bray_Curtis)) +
        geom_boxplot() +
        ylab("Bray-Curtis distance") +
        theme_bw() +
        xlab("") +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
    dist_Reproduction_plots[[pop]] <- p_boxDist

}


# create the panels
for(pop in unique(metaM9089[,PopID])){
    pp <- cowplot::plot_grid(umap_plots[[pop]],
                             pca_plots[[pop]],
                             dist_Species_plots[[pop]],
                             alpha_Species_plots[[pop]],
                             dist_Reproduction_plots[[pop]],
                             alpha_Reproduction_plots[[pop]],
                             ncol = 2,
                             nrow = 3,
                             rel_heights = c(1.5,1,1),
                             labels = c("A","B","C","E","D","F")
    )
    ggsave(pp, file = paste0("figures/FigSX_",pop,".pdf"), width= 8, height = 8)
}

print("PERMANOVA results:")
print(stats)
print("Linear model results:")
print(lapply(alpha_stats, anova))
