# Porthmeus
# 29.07.21


require(data.table)
require(DESeq2)
if(!("qiime2R" %in% installed.packages()[,"Package"])){
    devtools::install_github("jbisanz/qiime2R")
}
require(qiime2R)
require(ggplot2)
require(vegan)
require(phyloseq)
require(ggVennDiagram)
require(cowplot)

AbFilt_FT = "../data/dada2_AbFilt_FT.csv"
AbFilt_Tax = "../data/dada2_AbFilt_Tax.csv"
AbFilt_meta = "../data/dada2_AbFilt_meta.csv"
AbFilt_fasta_tree = "../data/dada2_AbFilt_fasta_tree.qza"
MetaLakes = "../MetaDataLakes.csv"


otu <- read.csv(AbFilt_FT, row.names = 1)
tax <- read.csv(AbFilt_Tax, row.names=1)
meta <- read.csv(AbFilt_meta, row.names=1)
metaLakes <- read.csv(MetaLakes,row.names=1)

# simplify the reproduction vector
map_vect <- c("SEX","NR","SEX","ASEX","SEX","SEX","SEX","SEX")
names(map_vect) <- unique(meta$ReproductiveMode)
meta[["Reproduction"]] <- map_vect[meta[["ReproductiveMode"]]]

# add the correct naming for the nutritional state of the water bodies
map_nut <- c("Hypereutr.", "Eutr.","Mesoeutr.")
map_nut <- factor(map_nut, levels= map_nut)
meta[["Nutrient_load"]] <- map_nut[metaLakes[meta$PopID,"NutrientLoad"]]
meta[["Water_body2"]] <- metaLakes[meta$PopID, "Waterbody2"]
meta[["Water_body"]] <- metaLakes[meta$PopID, "waterbody"]
meta[["Water_body"]] <- gsub("River","flowing",gsub("Lake","standing",meta[["Water_body"]]))
meta[["LocationIDshort"]] <- sapply(meta$LocationID, function(x) strsplit(x, split ="/")[[1]][2])

# reorder the sampling sites, so that they occur in a similar order as no the map
siteOrder <- c("M34","M72","M71","M70","M31","M67","M89","M90","M86","M85","M84","M79","M78","M109","M28","M26","M44","M108","M107","M83","R12")
meta[["SiteID"]] <- factor(meta[["PopID"]], levels = siteOrder)
meta[["Nutrient_load"]] <- factor(meta[["Nutrient_load"]], levels = rev(map_nut))


# remove rare species (lower 25% of the data)
countsPerESV <- apply(otu,1,sum)
thrs <- quantile(countsPerESV,0.25)
summary(countsPerESV)/sum(otu)
plot(ecdf(log10(countsPerESV)))
abline(v = log10(thrs))
hist(log10(countsPerESV))
abline(v=log10(thrs))

ESVs <- names(countsPerESV)[countsPerESV > thrs]

# remove species wich are not annotated on class level
otu <- otu[ESVs,]
tax <- tax[ESVs,]

# create short tax numbers
len <- 2
shESV <- sapply(rownames(tax), function(x) paste(strsplit(x, split="")[[1]][1:len],collapse=""))
while(sum(duplicated(shESV))){
    len <- len+1
    shESV <- sapply(rownames(tax), function(x) paste(strsplit(x, split="")[[1]][1:len],collapse=""))
}
tax[["shortID"]] <- factor(shESV, levels = shESV[order(tax[["Class"]])])

# add a tax phylum_class anno
tax[["PhyClass"]] <- paste(tax[["Phylum"]], tax[["Class"]], tax[["Order"]], sep = "/")
tax[["ClassOrder"]] <- paste(tax[["Class"]], tax[["Order"]], sep = "/")
subVec <- list("NA" = "uncult.",
               "uncultured bacterium" = "uncult.",
               "uncultured organism" = "uncult.",
               "bacteriales" = "bac.",
               "bacteria" = "bac.",
               "metagenome" = "uncult.",
               "Subgroup " = "subgr.",
               "WS6 (Dojkabacteria)" = "Dojkabac.",
               "SL56 marine group" = "SL56",
               "ales" = ".",
               " clade(Marine group B)" = "",
               " Incertae Sedis" = "",
               " clade" = "",
               " VC2.1 Bac22" = "",
               "Candidatus " = "",
               "JGI 0000069-P22" = "uncult."
)
for(n in names(subVec)){
tax[["PhyClass"]] <- gsub(n,subVec[[n]],tax[["PhyClass"]],fixed=TRUE)
}

# add another category phylum/class were the less abundant bacteria are grouped to "other"
otu_other <- matrix(0, ncol = ncol(otu), nrow = length(unique(tax[["PhyClass"]])), 
                    dimnames = list(unique(tax[["PhyClass"]]), colnames(otu)))
for(n in rownames(otu_other)){
    otu_other[n,] <- colSums(otu[rownames(tax[tax[["PhyClass"]] %in% n,]),])
}
contr <- cumsum(sort(rowSums(otu_other)/sum(otu_other)))

rmPhyClass <- names(contr[!(contr > 0.05)])
tax[["PhyClass2"]] <- tax[["PhyClass"]]
tax[tax[["PhyClass"]] %in% rmPhyClass, "PhyClass2"] <- "_other"
# test the results on ESV level
otu2 <- otu + 1
deseq <- DESeqDataSetFromMatrix(otu2, colData = meta[colnames(otu),], design = ~ Reproduction+Species+Nutrient_load + Water_body)
deseq <- DESeq(deseq)

getMainEffectResults <- function(obj, lfcThreshold = 0, alpha = 0.1){
    # get the summaries for the main effects
    
    # get the variable names of the main effects
    design <- as.character(design(obj))[2]
    main <- gsub(" ","",strsplit(gsub("*","+",design, fixed=TRUE), split = "+", fixed = T)[[1]])
    
    # remove interaction terms (not covered here)
    inter <- grep(":",main, fixed=TRUE)
    if(length(inter) >0){
        main <- main[-inter]
    }
    
    
    resultsSummary <- list()
    for(variable in main){
        # get the combination of levels in the variable
        cntrst <- combn(unique(as.character(colData(obj)[,variable])),2)
        colnames(cntrst) <- apply(cntrst, 2, paste, collapse = "_vs_")
        
        # get the results
        resultTabs <- list()
        for(con in colnames(cntrst)){
            result <- results(obj, contrast = c(variable, cntrst[,con]), lfcThreshold = lfcThreshold, alpha = alpha)
            resultTabs[[con]] <- result
        }

        # create a venn mat
        vennMat <- matrix(0, ncol = ncol(cntrst), nrow = nrow(obj), 
                          dimnames = list(c(rownames(result)), contrasts = colnames(cntrst)))
        vennList <- list()
        for(con in names(resultTabs)){
            res <- data.table(as.data.frame(resultTabs[[con]]), keep.rownames=TRUE)
            genes <- res[abs(log2FoldChange) > lfcThreshold & !is.na(padj) & padj < alpha, rn]
            vennMat[genes,con] <- 1
            vennList[[con]] <- genes
        }

        # save the results into one list
        resultsSummary[[variable]] <- list(vennMat = vennMat, vennList = vennList, results = resultTabs)
    }
    return(resultsSummary)
} 

    
result <- getMainEffectResults(deseq, alpha = 0.05, lfcThreshold = 1)


# get the indicator species for each condition tested
for(valName in names(result)){
    value <- result[[valName]]
    cntrsts <- colnames(value[["vennMat"]])
    features <- unique(as.vector(unlist(strsplit(cntrsts, split="_vs_"))))
    indicators <- list()
    # add those ESVs which are gradually abundant from one condition to the next
    ESVs <- rownames(value[["vennMat"]])[apply(value[["vennMat"]], 1, sum) ==length(cntrsts)]
    for(feat in features){
        indicators[[feat]] <- ESVs
    }
    if(length(cntrsts) > 1){

        # add those ESVs which are indicators for a specific condition
        ESVs2 <- rownames(value[["vennMat"]])[apply(value[["vennMat"]], 1, sum) ==length(cntrsts)-1]
        for(esv in ESVs2){
            cnt <- cntrsts[as.logical(value[["vennMat"]][esv,])]
            cnt <- as.vector((sapply(cnt, function(x) strsplit(x, split = "_vs_")[[1]])))
            cnt <- cnt[duplicated(cnt)]
            indicators[[cnt]] <- c(indicators[[cnt]], esv)
        }
    }
    barTabs <- list()
    for(con in names(indicators)){
        resultTabs <- grep(con, names(value[["results"]]),value = TRUE)
        barTab <- list()
        for(tabname in resultTabs){
            tab <- value[["results"]][[tabname]][indicators[[con]],]

            # correct the direction if necessary
            if(!startsWith(tabname,con)){
                tab[["log2FoldChange"]] <- tab[["log2FoldChange"]]*-1
            }
            tab[["vs."]] <- gsub(paste(paste0(con,"_vs_"),paste0("_vs_",con), sep = "|"),"",tabname)
            tab[["ESV"]] <- rownames(tab)
            barTab[[tabname]] <- tab
        }
        barTab <- as.data.frame(do.call(rbind,barTab))
        barTab <- merge(barTab, tax, by.x = "ESV",by.y=0)
        barTab[["Indicator"]] <- con
        barTabs[[con]] <- barTab
    }
    barTab <- do.call(rbind, barTabs)
    if(valName == "Nutrient_load"){
        barTab[["Indicator"]] <- relevel(factor(barTab[["Indicator"]]),"Mesoeutr.")
        barTab[["vs."]] <- relevel(factor(barTab[["vs."]]),"Mesoeutr.")
    }
    
    # reorder 
    barTab[["shortID"]] <- factor(barTab[["shortID"]], level = unique(barTab[["shortID"]][order(paste0(barTab[["vs."]],barTab[["PhyClass"]],barTab[["log2FoldChange"]]))]))
    
    if(length(cntrsts) == 1){
        barTab <- barTab[!duplicated(barTab[["shortID"]]),]
        p_barplot <- ggplot(barTab, aes(x=log2FoldChange, y = shortID, fill = PhyClass))+
            geom_bar(stat = "identity",position = position_dodge(width = 0.9),width = 0.8, color = "black") +
            theme_bw() +
            ylab("ESV") +
            xlab("Log2 fold change") +
            labs(fill = "Phylum/Class/Order")+
            ggtitle(gsub("vs","vs.",gsub("_"," ",paste0(toupper(substr(cntrsts[[1]],1,1)),substr(cntrsts[[1]],2,nchar(cntrsts[[1]]))))))
    } else {
        p_barplot <- ggplot(barTab, aes(x=log2FoldChange, y = shortID, fill = PhyClass, color = vs.))+
            geom_bar(stat = "identity",position = position_dodge(width = 0.9),width = 0.8) +
            scale_color_grey()+
            theme_bw() +
            ylab("ESV") +
            xlab("Log2 fold change") +
            labs(fill = "Phylum/Class/Order")+
            facet_wrap(~Indicator, scale ="free_y")
            #theme(legend.position = "bottom", legend.direction = "vertical") +
            #guides(fill = guide_legend(title.position = "top"),
            #       color = guide_legend(title.position = "top"))

    }
        p_barplot



    # create a bar/boxplot for the Class composition of each sample/condition
    otu_frac <- t(t(otu+0.1)/colSums(otu))
    otu_tab <- melt(data.table(otu_frac,keep.rownames=TRUE), id.vars="rn")
    otu_tab <- merge(otu_tab, data.table(tax, keep.rownames=TRUE), by = "rn")
    otu_tab <- otu_tab[,.(rel_abundance = mean(value)), by = .(PhyClass2, variable)]
    otu_tab <- merge(otu_tab, data.table(meta, keep.rownames=TRUE), by.x = "variable", by.y = "rn")

    p_boxAbund <- ggplot(otu_tab, aes(x=PhyClass2, y=rel_abundance)) +
        geom_boxplot(aes_string(fill = valName)) +
        scale_y_log10() +
        theme_bw()+
        labs(x ="Phylum/Class/Order", y = "Relative abundance per sample")+
        guides(fill = guide_legend(nrow = 1, title.position = "top"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
              legend.position = "top")
    p_boxAbund


    otu_long <- melt(data.table(otu_frac, keep.rownames=TRUE), id.vars = "rn")
    otu_long <- merge(otu_long, data.table(tax, keep.rownames=TRUE), by = "rn")
    otu_long <- otu_long[,.(abundance = sum(value)), by = .(PhyClass2, Sample = variable)]
    otu_long <- merge(otu_long, data.table(meta, keep.rownames=TRUE), by.x = "Sample", by.y = "rn")
    otu_long <- otu_long[,.(abundance = sum(abundance)), by = .(PhyClass2, Condition = get(valName))]
    otu_long[, rel_abundance := abundance/sum(abundance), by=Condition]
    
    
    if(length(cntrsts) > 1){
        p_barAbund<- ggplot(otu_long, aes(x=Condition, y = rel_abundance, fill = PhyClass2)) +
           geom_bar(stat="identity",color = "black") +
           theme_bw() +
           labs(x= valName, y = "Relative abundance per condition", fill = "Phylum/Class/Order")+
           theme(legend.position = "right",
                 axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1)) +
           guides(fill = guide_legend(ncol = 1, title.position = "top"))

#        p_venn <- ggVennDiagram(value[["vennList"]],  labels = "none",) +
#            scale_fill_gradient(low="white",high="gray80")+
#            scale_color_manual(values=rep("black",length(cntrsts)))

        p_top <- plot_grid(p_barAbund, p_boxAbund, nrow = 1, rel_widths = c(1,1,0.8), labels="AUTO")
        #p_bottom <- plot_grid(p_venn, p_barplot, nrow = 1, rel_width = c(1,3), labels = c("C","D"))
        p_all <- plot_grid(p_top, p_barplot, ncol = 1,labels = c(NA,"C"), rel_heights=c(1,1))

    }else{
        p_barAbund<- ggplot(otu_long, aes(x=Condition, y = rel_abundance, fill = PhyClass2)) +
           geom_bar(stat="identity",color = "black") +
           theme_bw() +
           labs(x= valName, y = "Relative abundance per condition", fill = "Phylum/Class/Order")+
           theme(legend.position = "right",
                 axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1)) +
           guides(fill = guide_legend(ncol = 1, title.position = "top"))

        p_top <- plot_grid(p_barAbund, p_boxAbund, nrow = 1, rel_widths = c(1,1), labels="AUTO")
        p_all <- plot_grid(p_top, p_barplot, ncol = 1,  labels=c(NA,"C"), rel_heights=c(1,1))
    }
    ggsave(p_all, file = paste0("figures/FigSX_",valName,"indicators.pdf"), width = 12, height =10)

}




