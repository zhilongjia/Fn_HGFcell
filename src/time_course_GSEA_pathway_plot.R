
enrichedPathway_list <- list()
for (fn in dir("../results/GSEA", pattern="gsea_report_for_na.*xls", recursive=TRUE, full.names=TRUE) ) {
    # fn <- "../results/GSEA/Fn_HGFcell_12h.GseaPreranked.1533520587759/gsea_report_for_na_neg_1533520587759.xls"
    print (fn)
    tp <- strsplit(strsplit(fn, "/")[[1]][4], "_|\\.")[[1]][3]
    pso_neg <- strsplit(strsplit(fn, "/")[[1]][5], "_")[[1]][5]
    
    tp_dir <- readr::read_tsv(fn)[,c("NAME", "NES", "FDR q-val")]
    # tp_dir <- dplyr::filter(tp_dir, `FDR q-val`<=0.1)
    tp_dir$tp <- tp
    enrichedPathway_list[[paste0(tp, "_", pso_neg)]] <- tp_dir
    
}

enrichedPathway_2h <- rbind(enrichedPathway_list[["2h_pos"]], enrichedPathway_list[["2h_neg"]]) %>% 
    dplyr::rename(NES_2h=NES) %>% dplyr::select(-c(tp, `FDR q-val`))
enrichedPathway_6h <- rbind(enrichedPathway_list[["6h_pos"]], enrichedPathway_list[["6h_neg"]]) %>% 
    dplyr::rename(NES_6h=NES) %>% dplyr::select(-c(tp, `FDR q-val`))
enrichedPathway_12h <- rbind(enrichedPathway_list[["12h_pos"]], enrichedPathway_list[["12h_neg"]]) %>% 
    dplyr::rename(NES_12h=NES) %>% dplyr::select(-c(tp, `FDR q-val`))
enrichedPathway_24h <- rbind(enrichedPathway_list[["24h_pos"]], enrichedPathway_list[["24h_neg"]]) %>% 
    dplyr::rename(NES_24h=NES) %>% dplyr::select(-c(tp, `FDR q-val`))
enrichedPathway_48h <- rbind(enrichedPathway_list[["48h_pos"]], enrichedPathway_list[["48h_neg"]]) %>% 
    dplyr::rename(NES_48h=NES) %>% dplyr::select(-c(tp, `FDR q-val`))

# enrichedPathway_2h <- rbind(enrichedPathway_list[["2h_pos"]], enrichedPathway_list[["2h_neg"]]) %>% 
#     dplyr::rename(NES_2h=NES, FDR_2h=`FDR q-val`) %>% dplyr::select(-tp)
# enrichedPathway_6h <- rbind(enrichedPathway_list[["6h_pos"]], enrichedPathway_list[["6h_neg"]]) %>% 
#     dplyr::rename(NES_6h=NES, FDR_6h=`FDR q-val`) %>% dplyr::select(-tp)
# enrichedPathway_12h <- rbind(enrichedPathway_list[["12h_pos"]], enrichedPathway_list[["12h_neg"]]) %>% 
#     dplyr::rename(NES_12h=NES, FDR_12h=`FDR q-val`) %>% dplyr::select(-tp)
# enrichedPathway_24h <- rbind(enrichedPathway_list[["24h_pos"]], enrichedPathway_list[["24h_neg"]]) %>% 
#     dplyr::rename(NES_24h=NES, FDR_24h=`FDR q-val`) %>% dplyr::select(-tp)
# enrichedPathway_48h <- rbind(enrichedPathway_list[["48h_pos"]], enrichedPathway_list[["48h_neg"]]) %>% 
#     dplyr::rename(NES_48h=NES, FDR_48h=`FDR q-val`) %>% dplyr::select(-tp)

enrichedPathway_timecourse <- Reduce(dplyr::full_join, list(enrichedPathway_2h,
                                                            enrichedPathway_6h,
                                                            enrichedPathway_12h,
                                                            enrichedPathway_24h,
                                                            enrichedPathway_48h))


enrichedPathway_timecourse$NAME <- gsub("kegg_","",lettercase::str_lower_case(enrichedPathway_timecourse$NAME ))
enrichedPathway_timecourse_mat <- as.matrix(enrichedPathway_timecourse[,-1])
rownames(enrichedPathway_timecourse_mat) <- enrichedPathway_timecourse$NAME
colnames(enrichedPathway_timecourse_mat) <- gsub("NES_","",colnames(enrichedPathway_timecourse_mat))

# enrichedPathway_timecourse_mat[enrichedPathway_timecourse_mat<1.2 & 
#                                    enrichedPathway_timecourse_mat>-1.2] =0
# 
# enrichedPathway_timecourse_mat <- enrichedPathway_timecourse_mat[names(which(rowSums(enrichedPathway_timecourse_mat)!=0)),]

pdf("../results/heatmap_pathway.pdf",width=10,height=25)
gplots::heatmap.2(enrichedPathway_timecourse_mat, dendrogram="row",
                  trace = "none", col="greenred", key=T,
                  keysize=1.1, key.title=NA, key.ylab="NES", key.xlab=NA,
                  srtCol=0, cexCol=1.5, cexRow=1.0, Colv=FALSE,
                  lwid = c(1,4), lhei = c(1,12),
                  margins = c(2, 26))
dev.off()


enrichedPathway <- Reduce(rbind, enrichedPathway_list)

enrichedPathway$tp <- factor(enrichedPathway$tp, 
                             levels=c("2h", "6h","12h", "24h", "48h"))

enrichedPathway1 <- dplyr::filter(enrichedPathway, `FDR q-val`<=0.05)

library(ggplot2)

enrichedPathway1 <- enrichedPathway[sample(790, 100),]
ggplot(enrichedPathway, aes(tp, NAME))+ 
    geom_point(aes(color=NES, size = -log10(`FDR q-val`))) +
    # guides(size=FALSE) + 
    # coord_polar("y") +
    scale_color_gradient2(midpoint=0, low="blue", mid="white",
                          high="red", space ="Lab" )



