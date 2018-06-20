
load("../results/Fn_HGF_Gene_exprs.Rdata")

################################################################################
# DEA over time
#subset data and pheno
subsample_pheno <- dplyr::bind_rows(dplyr::filter(pheno, state=="F" ),
                                    dplyr::filter(pheno, timepoint=="0h" ) )

subsample_pheno$label <- factor(as.character(subsample_pheno$label), levels=c("C_0h", "F_2h", "F_6h", "F_12h", "F_24h",  "F_48h") )

Expdesign <- model.matrix(~0+ subsample_pheno$label+ subsample_pheno$individual  )
colnames(Expdesign) <- gsub("subsample_pheno\\$individual|subsample_pheno\\$label","", colnames(Expdesign) )

subsample_dge <- dge[,subsample_pheno$sample_ID]

# DEA
v <- voom(subsample_dge, Expdesign)
Expfit1 <- lmFit(v, Expdesign)
cont.matrix <- makeContrasts(F_2h-C_0h, F_6h-F_2h, F_12h-F_6h, F_24h-F_12h, F_48h-F_24h, levels=Expdesign)
fit2 <- contrasts.fit(Expfit1, cont.matrix)
Expfit2 <- eBayes(fit2)

#########################
DEA_list <- list()
for (vs_i in colnames(cont.matrix) ) {
    
    # vs_i <- colnames(cont.matrix)[1]
    print(vs_i)
    
    GSE_limma <- topTable(Expfit2, coef=vs_i, number=Inf) %>% 
        tibble::rownames_to_column(var=vs_i )
    
    DEG_limma_filter <- dplyr::filter(GSE_limma, adj.P.Val<=0.05, abs(logFC)>=2 )

    
    DEA_list[["limma_res"]][[vs_i]] <- GSE_limma
    DEA_list[["DEG"]][[vs_i]] <- DEG_limma_filter[[vs_i]]
    DEA_list[["DEGexp"]][[vs_i]] <- v$E[DEG_limma_filter[[vs_i]],]
}
sapply(DEA_list[["DEG"]], length)

save.image("../results/time_course_DEA_stage1.RData")


################################################################################
# gene symbol annotation
sapply(DEA_list$DEG, length)

for (vs_i in names(DEA_list$limma_res) ){

    write.table(DEA_list$limma_res[[vs_i]][,c(vs_i, "logFC")], 
                file = paste0( "../results/pathview_time_course/DEA_list_limma_res_", gsub(" ","",vs_i) ,"_pathview.txt"), 
                sep="\t", quote=FALSE, row.names = FALSE)
    
}

DEG_logFC_list <- list()
for (tp_i in names(DEA_list$limma_res) ){
    DEG_logFC <- DEA_list$limma_res[[tp_i]][,c(tp_i, "logFC")]
    colnames(DEG_logFC) <- c("gene", tp_i)
    DEG_logFC_list[[tp_i]] <- DEG_logFC
}

DEG_logFC_df <- Reduce(dplyr::inner_join, DEG_logFC_list)
write.table(DEG_logFC_df, 
            file = paste0( "../results/pathview_time_course/DEG_logFC_df.txt"), 
            sep="\t", quote=FALSE, row.names = FALSE)


# intersect_time_indep_DEG <- Reduce(intersect, DEA_list$DEG)
union_time_course_DEG <- Reduce(union, DEA_list$DEG)
source("~/Dropbox/zzz.R")
# intersect_time_indep_gene_KEGG_GO <- gene2KEGG( intersect_time_indep_DEG )
union_time_course_DEgene_KEGG_GO <- gene2KEGG( union_time_course_DEG )

# write.table(intersect_time_indep_gene_KEGG_GO, file="../results/intersect_time_indep_gene_KEGG_GO.txt", sep="\t", quote=FALSE, row.names = FALSE)
write.table(union_time_course_DEgene_KEGG_GO, file="../results/union_time_course_DEgene_KEGG_GO.txt", sep="\t", quote=FALSE, row.names = FALSE)


################################################################################
# PCA of DEGs

####################################################
#united DEGs
DEexprs_time_course_united <- v$E[union_time_course_DEG,]

gplots::heatmap.2(DEexprs_time_course_united, trace = "none", col="greenred", scale="row",
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=TRUE, labRow=FALSE,
                  ylab=paste(nrow(DEexprs_time_course_united), "DE Genes"), margins = c(4, 2))

# order based on time
sample_id_ordered <- dplyr::arrange(subsample_pheno, state, as.numeric(sub("h", "",timepoint))) %>% dplyr::select(sample_ID) %>% unlist(use.names = FALSE)
gplots::heatmap.2(DEexprs_time_course_united[,sample_id_ordered], trace = "none", col="greenred", scale="row",
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=FALSE, labRow=FALSE,
                  ylab=paste(nrow(DEexprs_time_course_united), "DE Genes"), margins = c(4, 2))



#########################################################
# 2D PCA

library(ggfortify)


#color group and time
autoplot(prcomp(t(DEexprs_time_course_united), scale=F), data=subsample_pheno, colour = "label", 
         size = 15, label = TRUE, label.colour="black", ts.colour="black" )

# color individ
autoplot(prcomp(t(DEexprs_time_course_united)), data=subsample_pheno, colour = "individual", 
         size = 15, label = TRUE, label.colour="black", ts.colour="black" )

# time
autoplot(prcomp(t(DEexprs_time_course_united)), data=subsample_pheno, colour = "timepoint", 
         size = 15, label = TRUE, label.colour="black", ts.colour="black" )


################################################################################
library(Vennerable)

w0 <- Venn(Sets=DEA_list$DEG, Weight = Weight)
# Weights(w0) <- log2(Weights(w0)+1)

w <- compute.Venn(w0, type = "ChowRuskey")

# parameters
gp <- VennThemes(w)
# gp[["Face"]][["11111"]]$fill <-  "green"
# 
# # setText size
# gp[["SetText"]][["Set1"]]$fontsize <- 20
# gp[["SetText"]][["Set2"]]$fontsize <- 20
# gp[["SetText"]][["Set3"]]$fontsize <- 20
# gp[["SetText"]][["Set4"]]$fontsize <- 20
# gp[["SetText"]][["Set5"]]$fontsize <- 20
# 
# 
# # setname location
# set_label <- VennGetSetLabels(w)
# set_label[2,"y"] <- -16
# set_label[2,"x"] <- -2
# set_label[1,"y"] <- 8
# set_label[3,"x"] <- -26
# set_label[4,"y"] <- -28
# set_label[4,"x"] <- -40
# 
# 
# w <- VennSetSetLabels(w, set_label)

grid::grid.newpage()
plot(w, gp=gp)

genes_in_each_set <- w0@IntersectionSets

save.image("../results/time_course_DEA_stage2.RData")


################################################################################
# Pathway & GO analysis of ChowRuskey genesets
library(clusterProfiler)
source("~/Dropbox/zzz.R")

genesEntrezID_in_each_set <- sapply(genes_in_each_set, symbol2entrezID )


# filter entrez gene less than 10
genesEntrezID_in_each_set_filtered <- genesEntrezID_in_each_set[names(which(sapply(genesEntrezID_in_each_set, length)>=10))]
genesEntrezID_in_each_set_filtered_KEGG <- compareCluster(genesEntrezID_in_each_set_filtered, fun='enrichKEGG')
dotplot(genesEntrezID_in_each_set_filtered_KEGG)

genesEntrezID_in_each_set_filtered_GO <- compareCluster(genesEntrezID_in_each_set_filtered, fun='enrichGO', OrgDb='org.Hs.eg.db')
dotplot(genesEntrezID_in_each_set_filtered_GO)


################################################################################
# Pathway & GO analysis of DEG at each timepoint comparsion

genesEntrezID_eachtimepoint <- sapply(DEA_list[["DEG"]], symbol2entrezID )

genesEntrezID_eachtimepoint_KEGG <- compareCluster(genesEntrezID_eachtimepoint, fun='enrichKEGG')
dotplot(genesEntrezID_eachtimepoint_KEGG, showCategory=30)

genesEntrezID_eachtimepoint_GO <- compareCluster(genesEntrezID_eachtimepoint, fun='enrichGO', OrgDb='org.Hs.eg.db')
dotplot(genesEntrezID_eachtimepoint_GO, showCategory=30)

save.image("../results/time_course_DEA_stage3.RData")

################################################################################




