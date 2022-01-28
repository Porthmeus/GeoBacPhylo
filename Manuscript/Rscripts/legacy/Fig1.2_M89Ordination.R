# Porthmeus
# 16.03.21

require(data.table)
require(DESeq2)
if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R)
require(ggplot2)
require(vegan)
require(phyloseq)
require(umap)
require(pals)
require(plotly)


AbFilt_FT = "../data/dada2_AbFilt_FT.csv"
AbFilt_Tax = "../data/dada2_AbFilt_Tax.csv"
AbFilt_meta = "../data/dada2_AbFilt_meta.csv"
AbFilt_fasta_tree = "../data/dada2_AbFilt_fasta_tree.qza"
MetaLakes = "../MetaDataLakes.csv"


otu <- read.csv(AbFilt_FT, row.names = 1)
tax <- read.csv(AbFilt_Tax, row.names=1)
meta <- read.csv(AbFilt_meta, row.names=1)

# simplify the reproduction vector
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]

# add a single count to the OTU matrix, stabilize variance, and shift the values to positive values
otu2 <- otu+1
otu_vst <- varianceStabilizingTransformation(as.matrix(otu2), fitType = "local")
if(min(otu_vst) <0){
    otu_vst <- otu_vst + min(otu_vst)*-1
}

phy_vst <- phyloseq(otu_table(otu_vst,
                              taxa_are_rows = TRUE),
                              tax_table(as.matrix(tax)),
                              sample_data(meta))

# calculate bray-curtis and jaccard distance
brayCurt <- phyloseq::distance(phy_vst, method = "bray")
jac <- phyloseq::distance(phy_vst, method = "jaccard")




# calculate the UMAP plot only for specific pop
for(specPop in c("M89","M90","M28","M44","M70","M72","M79","M107","M108","M109","M26","M44","M67","M83")){
    print(specPop)
    meta2 <- data.table(meta,keep.rownames=TRUE)
    specPop_sel <- meta2[PopID == specPop, rn]
    jac_specPop <- as.matrix(jac)[specPop_sel,specPop_sel]
    #umap_jac_specPop <- umap(jac_specPop,dist=TRUE)[["layout"]]
   # colnames(umap_jac_specPop) <- c("UMAP1_jac","UMAP2_jac")
    brayCurt_specPop <- as.matrix(brayCurt)[specPop_sel,specPop_sel]
   # umap_brayCurt_specPop <- umap(brayCurt_specPop,dist=TRUE)[["layout"]]
   # colnames(umap_brayCurt_specPop) <- c("UMAP1_bray","UMAP2_bray")

   # p_dat_specPop <- merge(meta,umap_jac_specPop,by=0)
   # p_dat_specPop <- merge(p_dat_specPop,umap_brayCurt_specPop,by.x="Row.names",by.y=0)

   # p_jac_specPop <- ggplot(p_dat_specPop, aes(x=UMAP1_jac, y=UMAP2_jac)) +
   #     geom_point(aes(color = Species, shape = Reproduction,size=3)) +
   #     theme_bw() +
   #     ggtitle("Jaccard") +
   #     ylab ("UMAP2") +
   #     xlab ("UMAP1") 
   # fln = file.path("figures",paste0("UMAP_jac_",specPop,".pdf"))
   # ggsave(p_jac_specPop,
   #        file = fln,
   #        width= 7,
   #        height = 6)


   # p_bray_specPop <- ggplot(p_dat_specPop, aes(x=UMAP1_bray, y=UMAP2_bray)) +
   #     geom_point(aes(color = Species, shape = Reproduction, size=3)) +
   #     theme_bw() +
   #     ggtitle("Bray-Curtis") +
   #     ylab ("UMAP2") +
   #     xlab ("UMAP1")

   # fln = file.path("figures",paste0("UMAP_bray_",specPop,".pdf"))
   # ggsave(p_bray_specPop,
   #        file = fln,
   #        width= 7,
   #        height = 6)



    pca_jac <- prcomp(jac_specPop)
    pca_jac[["varExpl"]] <- round(pca_jac[["sdev"]]^2/sum(pca_jac[["sdev"]]^2)*100, digits = 2)
    pca_bray <- prcomp(brayCurt_specPop)
    pca_bray[["varExpl"]] <- round(pca_bray[["sdev"]]^2/sum(pca_bray[["sdev"]]^2)*100, digits = 2)
    pcs_jac <- pca_jac[["x"]]
    colnames(pcs_jac) <- paste(colnames(pcs_jac),"Jaccard",sep="_")



    pcs_bray <- pca_bray[["x"]]
    colnames(pcs_bray) <- paste(colnames(pcs_bray),"Bray",sep="_")
    pca_dat <- merge(pcs_bray,pcs_jac, by=0)
    pca_dat <- merge(pca_dat, meta, by.x="Row.names",by.y=0)



    p_jac <- ggplot(pca_dat, aes(x=PC1_Jaccard, y=PC2_Jaccard, color = Species,shape=LocationID)) +
        geom_point(size=3) +
        ggtitle("Jaccard") +
        ylab(paste0("PC2 (",pca_jac[["varExpl"]][2],"%)")) +
        xlab(paste0("PC1 (",pca_jac[["varExpl"]][1],"%)")) +
        theme_bw()

    fln <- file.path("figures",paste0("PCA_jac_",specPop,".pdf"))
    ggsave(p_jac,
           file = fln,
           width= 4.8,
           height = 4)
    
    p3d <- plot_ly(pca_dat, x = ~PC1_Jaccard, y = ~PC2_Jaccard, z = ~PC3_Jaccard,
                   color = ~Species,
                   colors = c("dark red","dark blue","dark green"),
                   symbol = ~as.factor(LocationID),
                   symbols = c("circle","cross","diamond"))
    p3d <- p3d %>% add_markers()
    p3d <- p3d %>% layout(scene = list( xaxis = list(title = paste0("PC1 (",pca_jac[["varExpl"]][1],"%)")),
                                       yaxis = list(title = paste0("PC2 (",pca_jac[["varExpl"]][2],"%)")),
                                       zaxis = list(title = paste0("PC3 (",pca_jac[["varExpl"]][3],"%)"))),
                          title = "Jaccard")
    fln <- file.path("figures",paste0("PCA3D_jac_",specPop,".html"))
    htmlwidgets::saveWidget(p3d, file = fln, selfcontained = TRUE)

    p_bray <- ggplot(pca_dat, aes(x=PC1_Bray, y=PC2_Bray, color = Species,shape=LocationID)) +
        geom_point(size=3) +
        ggtitle("Bray-Curtis") +
        ylab(paste0("PC2 (",pca_bray[["varExpl"]][2],"%)")) +
        xlab(paste0("PC1 (",pca_bray[["varExpl"]][1],"%)")) +
        theme_bw()

    fln <- file.path("figures",paste0("PCA_bray_",specPop,".pdf"))
    ggsave(p_bray,
           file = fln,
           width= 4.8,
           height = 4)

    p3d <- plot_ly(pca_dat, x = ~PC1_Bray, y = ~PC2_Bray, z = ~PC3_Bray,
                   color = ~Species,
                   colors = c("dark red","dark blue","dark green"),
                   symbol = ~as.factor(LocationID),
                   symbols = c("circle","cross","diamond"))
    p3d <- p3d %>% add_markers()
    p3d <- p3d %>% layout(scene = list( xaxis = list(title = paste0("PC1 (",pca_bray[["varExpl"]][1],"%)")),
                                       yaxis = list(title = paste0("PC2 (",pca_bray[["varExpl"]][2],"%)")),
                                       zaxis = list(title = paste0("PC3 (",pca_bray[["varExpl"]][3],"%)"))),
                          title = "Bray-Curtis")
    fln <- file.path("figures",paste0("PCA3D_bray_",specPop,".html"))
    htmlwidgets::saveWidget(p3d, file = fln, selfcontained = TRUE)
}



pca_jac <- prcomp(jac)
pca_jac[["varExpl"]] <- round(pca_jac[["sdev"]]^2/sum(pca_jac[["sdev"]]^2)*100, digits = 2)
pcs_jac <- pca_jac[["x"]]
colnames(pcs_jac) <- paste(colnames(pcs_jac),"Jaccard",sep="_")

pca_bray <- prcomp(brayCurt)
pca_bray[["varExpl"]] <- round(pca_bray[["sdev"]]^2/sum(pca_bray[["sdev"]]^2)*100, digits = 2)
pcs_bray <- pca_bray[["x"]]
colnames(pcs_bray) <- paste(colnames(pcs_bray),"Bray",sep="_")
pca_dat <- merge(pcs_bray,pcs_jac, by=0)
pca_dat <- merge(pca_dat, meta, by.x="Row.names",by.y=0)

p3d <- plot_ly(pca_dat, x = ~PC1_Jaccard, y = ~PC2_Jaccard, z = ~PC3_Jaccard,
               color = ~PopID,
               #colors = c("dark red","dark blue","dark green"),
               symbol = ~Species,
               symbols = c("cross","circle","diamond"))
p3d <- p3d %>% add_markers()
p3d <- p3d %>% layout(scene = list( xaxis = list(title = paste0("PC1 (",pca_jac[["varExpl"]][1],"%)")),
                                   yaxis = list(title = paste0("PC2 (",pca_jac[["varExpl"]][2],"%)")),
                                   zaxis = list(title = paste0("PC3 (",pca_jac[["varExpl"]][3],"%)"))),
                      title = "Jaccard")
fln <- file.path("figures","PCA3D_jac_all.html")
htmlwidgets::saveWidget(p3d, file = fln, selfcontained = TRUE)

p3d <- plot_ly(pca_dat, x = ~PC1_Bray, y = ~PC2_Bray, z = ~PC3_Bray,
               color = ~Species,
               colors = c("dark red","dark blue","dark green"),
               symbol = ~as.factor(LocationID),
               symbols = c("circle","cross","diamond"))
p3d <- p3d %>% add_markers()
p3d <- p3d %>% layout(scene = list( xaxis = list(title = paste0("PC1 (",pca_bray[["varExpl"]][1],"%)")),
                                   yaxis = list(title = paste0("PC2 (",pca_bray[["varExpl"]][2],"%)")),
                                   zaxis = list(title = paste0("PC3 (",pca_bray[["varExpl"]][3],"%)"))),
                      title = "Bray-Curtis")
fln <- file.path("figures","PCA3D_bray_all.html")
htmlwidgets::saveWidget(p3d, file = fln, selfcontained = TRUE)

