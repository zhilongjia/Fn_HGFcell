

################################################################################
# import RNA-Seq read count data
# Note: The fifth column provides the expected read count in each transcript, 
# which can be utilized by tools like EBSeq, DESeq and edgeR for differential 
# expression analysis. ref: https://github.com/bli25broad/RSEM_tutorial

library(dplyr)
gene_count_rpkm_fn <- dir("../data/rawdata_fromWenyanKang/Analysis_Report/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression",
    pattern="gene", full.names=TRUE)

gene_count_list <- list()
for (fn_i in gene_count_rpkm_fn ) {
    gene_count_rpkm <- readr::read_tsv(fn_i) %>% 
        dplyr::select(Symbol, expected_count) %>% 
        dplyr::group_by(Symbol) %>%
        summarise_all(funs(sum))
    
    sample_id <- sub(".gene.fpkm.xls","",basename(fn_i) )
    colnames(gene_count_rpkm)[2] <- sample_id
    
    gene_count_list[[sample_id]] <- gene_count_rpkm
}

sample_count_df <- Reduce(dplyr::inner_join, gene_count_list )

################################################################################
# Normlization
library(edgeR)
library(limma)

# All gene count noramlization
S_raw <- as.matrix(select(sample_count_df, -Symbol) )
rownames(S_raw) <- sample_count_df$Symbol
dge <- DGEList(counts=S_raw)
dge <- calcNormFactors(dge)

pheno <- tibble(sample_ID=colnames(S_raw), 
                    individual=substr(colnames(S_raw), 1,1) ,
                    timepoint=paste0(gsub("[[:alpha:]]","", colnames(S_raw) ), "h") , 
                    state=gsub("[[:alpha:]]*[[:digit:]]","", colnames(S_raw)))

pheno[pheno$sample_ID %in% c("B0", "C0", "D0", "E0", "F0"), "state" ] <- "C"

pheno$individual <- as.factor(pheno$individual)
pheno$state <- factor(pheno$state, levels = c("C", "F"))
pheno$label <- factor(paste(pheno$state, pheno$timepoint, sep="_"))

save(dge, pheno, file="../results/Fn_HGF_Gene_exprs.Rdata")

################################################################################
# DEA for each timepoint
DEA_list <- list()
# limma_res_list <- list()
# DEG_list <- list()
# DEGexp_list <- list()
for (tp in c("2h", "6h", "12h", "24h", "48h") ) {
    
    print(tp)
    #subset data and pheno
    subsample_pheno <- dplyr::filter(pheno, timepoint==tp)
    Expdesign <- model.matrix(~subsample_pheno$individual + subsample_pheno$state)
    
    subsample_dge <- dge[,subsample_pheno$sample_ID]
    
    # DEA
    v <- voom(subsample_dge, Expdesign)
    Expfit1 <- lmFit(v, Expdesign)
    Expfit2 <- eBayes(Expfit1)
    GSE_limma <- topTable(Expfit2, coef=tail(colnames(Expdesign), 1), number=Inf) %>% 
        tibble::rownames_to_column(var=tp )
    
    # limma_res_list[[ tp ]] <- GSE_limma
    
    DEG_limma_filter <- dplyr::filter(GSE_limma, adj.P.Val<=0.05, abs(logFC)>=2 )
    # DEG_list[[tp]] <- DEG_limma_filter[[tp]]
    
    # DEGexp_list[[tp]] <- v$E[DEG_limma_filter[[tp]],]
    # print (paste(nrow(DEGexp_list[[tp]]), "DEGs"))
    
    DEA_list[["limma_res"]][[tp]] <- GSE_limma
    DEA_list[["DEG"]][[tp]] <- DEG_limma_filter[[tp]]
    DEA_list[["DEGexp"]][[tp]] <- v$E[DEG_limma_filter[[tp]],]
}

################################################################################
# gene symbol annotation
sapply(DEA_list$DEG, length)

intersect_time_indep_DEG <- Reduce(intersect, DEA_list$DEG)
union_time_indep_DEG <- Reduce(union, DEA_list$DEG)
source("./gene2KEGG.R")
intersect_time_indep_gene_KEGG_GO <- gene2KEGG( intersect_time_indep_DEG )
union_time_indep_DEgene_KEGG_GO <- gene2KEGG( union_time_indep_DEG )

write.table(intersect_time_indep_gene_KEGG_GO, file="../results/intersect_time_indep_gene_KEGG_GO.txt", sep="\t", quote=FALSE, row.names = FALSE)
write.table(union_time_indep_DEgene_KEGG_GO, file="../results/union_time_indep_DEgene_KEGG_GO.txt", sep="\t", quote=FALSE, row.names = FALSE)


################################################################################
# PCA of DEGs
Expdesign <- model.matrix(~0+ pheno$individual+ pheno$label )
v <- voom(dge,Expdesign)

####################################################
# intersect
DEexprs_time_indep_intersect <- v$E[intersect_time_indep_DEG,]

###################
# heatmap
gplots::heatmap.2(DEexprs_time_indep_intersect, trace = "none", col="greenred", scale="row",
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=TRUE,
                  ylab=paste(nrow(DEexprs_time_indep_intersect), "DE Genes"), margins = c(4, 6))

# order based on time
sample_id_ordered <- dplyr::arrange(pheno, state, as.numeric(sub("h", "",timepoint))) %>% dplyr::select(sample_ID) %>% unlist(use.names = FALSE)
gplots::heatmap.2(DEexprs_time_indep_intersect[,sample_id_ordered], trace = "none", col="greenred", scale="row",
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=FALSE,
                  ylab=paste(nrow(DEexprs_time_indep_intersect), "DE Genes"), margins = c(4, 6))

####################################################
#united DEGs
DEexprs_time_indep_united <- v$E[union_time_indep_DEG,]

gplots::heatmap.2(DEexprs_time_indep_united, trace = "none", col="greenred", scale="row",
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=TRUE, labRow=FALSE,
                  ylab=paste(nrow(DEexprs_time_indep_united), "DE Genes"), margins = c(4, 2))

# order based on time
sample_id_ordered <- dplyr::arrange(pheno, state, as.numeric(sub("h", "",timepoint))) %>% dplyr::select(sample_ID) %>% unlist(use.names = FALSE)
gplots::heatmap.2(DEexprs_time_indep_united[,sample_id_ordered], trace = "none", col="greenred", scale="row",
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=FALSE, labRow=FALSE,
                  ylab=paste(nrow(DEexprs_time_indep_united), "DE Genes"), margins = c(4, 2))



#########################################################
# 2D PCA
DEexprs_time_indep <- DEexprs_time_indep_united

library(ggfortify)


#color group and time
autoplot(prcomp(t(DEexprs_time_indep), scale=F), data=pheno, colour = "label", 
         size = 15, label = TRUE, label.colour="black", ts.colour="black" )

# color individ
autoplot(prcomp(t(DEexprs_time_indep)), data=pheno, colour = "individual", 
         size = 15, label = TRUE, label.colour="black", ts.colour="black" )

# time
autoplot(prcomp(t(DEexprs_time_indep)), data=pheno, colour = "timepoint", 
         size = 15, label = TRUE, label.colour="black", ts.colour="black" )

############################
# 3D-PCA
colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
            "#D55E00", "#CC79A7", "green", "#CC6666", "#9999CC")

names(colors) = unique(pheno$label)

# scaled or not
pc <- prcomp(t(DEexprs_time_indep), scale=FALSE)
pc <- prcomp(t(DEexprs_time_indep), scale=TRUE)

s3d <- scatterplot3d::scatterplot3d(pc$x[,c(1,3, 2)], pch=19, color=colors[pheno$label], 
                                    cex.symbols=5, 
                                    ylab = "",
                                    # ylab=paste0("PC1:\t", formattable::percent(as.data.frame(summary(pc)$importance)[2,"PC1"] )), 
                                    zlab=paste0("PC2:  ", formattable::percent(as.data.frame(summary(pc)$importance)[2,"PC2"] )), 
                                    xlab=paste0("PC1:  ", formattable::percent(as.data.frame(summary(pc)$importance)[2,"PC1"] ))
)
# legend("topleft", legend = names(colors), col = colors, pch = 19)
text(s3d$xyz.convert(pc$x[, c(1,3,2)]), labels = pheno$sample_ID, cex= 0.8, col ="grey4" )
mtext(paste0("PC3:  ", formattable::percent(as.data.frame(summary(pc)$importance)[2,"PC3"] )), side=4,las=2, padj=17, line=-7)

################################################################################
library(Vennerable)

w0 <- Venn(Sets=DEA_list$DEG, Weight = Weight)
# Weights(w0) <- log2(Weights(w0)+1)

w <- compute.Venn(w0, type = "ChowRuskey")

# parameters
gp <- VennThemes(w)
gp[["Face"]][["11111"]]$fill <-  "green"

# setText size
gp[["SetText"]][["Set1"]]$fontsize <- 20
gp[["SetText"]][["Set2"]]$fontsize <- 20
gp[["SetText"]][["Set3"]]$fontsize <- 20
gp[["SetText"]][["Set4"]]$fontsize <- 20
gp[["SetText"]][["Set5"]]$fontsize <- 20


# setname location
set_label <- VennGetSetLabels(w)
set_label[2,"y"] <- -16
set_label[2,"x"] <- -2
set_label[1,"y"] <- 8
set_label[3,"x"] <- -26
set_label[4,"y"] <- -28
set_label[4,"x"] <- -40


w <- VennSetSetLabels(w, set_label)

grid::grid.newpage()
plot(w, gp=gp)

genes_in_each_set <- w0@IntersectionSets

save.image("../results/time_indep_DEA_stage1.RData")


################################################################################
# Pathway & GO analysis of ChowRuskey genesets
library(clusterProfiler)
source("~/Dropbox/zzz.R")

genesEntrezID_in_each_set <- sapply(genes_in_each_set, symbol2entrezID )

genesEntrezID_in_each_set_KEGG <- compareCluster(genesEntrezID_in_each_set, fun='enrichKEGG')
dotplot(genesEntrezID_in_each_set_KEGG)

genesEntrezID_in_each_set_GO <- compareCluster(genesEntrezID_in_each_set, fun='enrichGO', OrgDb='org.Hs.eg.db')
dotplot(genesEntrezID_in_each_set_GO)

# filter entrez gene less than 10
genesEntrezID_in_each_set_filtered <- genesEntrezID_in_each_set[names(which(sapply(genesEntrezID_in_each_set, length)>=10))]
genesEntrezID_in_each_set_filtered_KEGG <- compareCluster(genesEntrezID_in_each_set_filtered, fun='enrichKEGG')
dotplot(genesEntrezID_in_each_set_filtered_KEGG)

genesEntrezID_in_each_set_filtered_GO <- compareCluster(genesEntrezID_in_each_set_filtered, fun='enrichGO', OrgDb='org.Hs.eg.db')
dotplot(genesEntrezID_in_each_set_filtered_GO)


################################################################################
# Pathway & GO analysis of DEG at each timepoint

genesEntrezID_eachtimepoint <- sapply(DEA_list[["DEG"]], symbol2entrezID )

genesEntrezID_eachtimepoint_KEGG <- compareCluster(genesEntrezID_eachtimepoint, fun='enrichKEGG')
dotplot(genesEntrezID_eachtimepoint_KEGG, showCategory=30)

genesEntrezID_eachtimepoint_GO <- compareCluster(genesEntrezID_eachtimepoint, fun='enrichGO', OrgDb='org.Hs.eg.db')
dotplot(genesEntrezID_eachtimepoint_GO, showCategory=30)

save.image("../results/time_indep_DEA_stage2.RData")

################################################################################






