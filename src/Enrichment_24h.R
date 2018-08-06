

# Enrichment analysis only for 24h

library(clusterProfiler)
source("~/Dropbox/zzz.R")
symbol2entrezID_df <- function(gene_symbols) {
    symbol_ENTREZID <- clusterProfiler::bitr(gene_symbols, fromType="SYMBOL", 
                                             toType="ENTREZID", OrgDb="org.Hs.eg.db")
    return(symbol_ENTREZID)
}


load("../results/time_indep_DEA_stage1.RData")

tp <- "24h"

DEGs <- DEA_list$DEG[[tp]]
limma_res <- DEA_list$limma_res[[tp]]
DEGexp <- DEA_list$DEGexp[[tp]]

write.table(as.data.frame(DEGexp), "../results/24h/DEGexp.tsv", quote=FALSE)
write.table(as.data.frame(limma_res), "../results/24h/limma_res.tsv", quote=FALSE, row.names = FALSE)


##heatmap
col1 <- colorRampPalette(c("blue", "white", "firebrick2"))
col1 <- "greenred"
gplots::heatmap.2(DEGexp, trace = "none", col=col1, scale="row",
                  labRow=F,
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=TRUE,
                  ylab=paste(nrow(DEGexp), "DE Genes"), margins = c(4, 4))



up_dn_Genes <- list()
up_dn_Genes$upGenes <- dplyr::filter(limma_res, adj.P.Val<=0.05, logFC>=2 ) %>% dplyr::select(tp) %>% unlist(use.names = FALSE)
up_dn_Genes$dnGenes <- dplyr::filter(limma_res, adj.P.Val<=0.05, logFC<=-2 ) %>% dplyr::select(tp) %>% unlist(use.names = FALSE)

################################################################################
# Pathway analysis
DEGs_EntrezID <- symbol2entrezID(DEGs)

kk <- enrichKEGG(gene = DEGs_EntrezID, organism = 'hsa', pvalueCutoff = 0.05)
head(kk)
dotplot(kk, showCategory=20)

### GSE
DEGs_Symbol_EntrezID <- symbol2entrezID_df(limma_res$`24h`) %>% 
    dplyr::left_join(limma_res, by=c("SYMBOL"="24h")) %>% 
    dplyr::arrange(-t)
DEGs_EntrezID_gse <- DEGs_Symbol_EntrezID$t
names(DEGs_EntrezID_gse) <- DEGs_Symbol_EntrezID$ENTREZID

kk2 <- gseKEGG(geneList     = DEGs_EntrezID_gse,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
dotplot(kk2, showCategory=20)

###
library(GSEABase)
c2 <- read.gmt("./cogena/inst/extdata/c2.all.v6.1.symbols.gmt")

egmt <- enricher(DEGs, TERM2GENE=c2)
head(egmt)
dotplot(egmt, showCategory=20)

####
# genesEntrezID_in_each_set <- sapply(up_dn_Genes, symbol2entrezID )
# 
# genesEntrezID_in_each_set_KEGG <- compareCluster(genesEntrezID_in_each_set, fun='enrichKEGG')
# dotplot(genesEntrezID_in_each_set_KEGG)


################################################################################
ego_BP <- enrichGO(gene          = DEGs_EntrezID,
                universe      = symbol2entrezID(gene_count_rpkm$Symbol),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

dotplot(ego_BP)
cnetplot(ego_BP)
goplot(ego_BP)

####
ego_MF <- enrichGO(gene          = DEGs_EntrezID,
                   universe      = symbol2entrezID(gene_count_rpkm$Symbol),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

dotplot(ego_MF)
cnetplot(ego_MF)
# goplot(ego_MF)

####
ego_CC <- enrichGO(gene          = DEGs_EntrezID,
                   universe      = symbol2entrezID(gene_count_rpkm$Symbol),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

dotplot(ego_CC, showCategory=20)
cnetplot(ego_CC, showCategory=20)
goplot(ego_CC, showCategory=20)

##########################################
library("cogena")

nClust <- 2:6
ncore <- 5
clMethods <- c("hierarchical","kmeans","pam")

genecl_result <- coExp(DEGexp, nClust=nClust, 
                       clMethods=clMethods, 
                       metric="correlation", 
                       method="complete", 
                       ncore=ncore, 
                       verbose=TRUE)

annoGMT <- "c2.cp.kegg.v6.1.symbols.gmt"
# annoGMT <- "c2.all.v6.1.symbols.gmt"
annoGMT <- "c2.cp.reactome.v6.1.symbols.gmt"; nClust =4; clMethods=c("kmeans")

annofile <- system.file("extdata", annoGMT, package="cogena")

sampleLabel <- factor(substr( colnames(DEGexp), 4,4), levels = c("C", "F"))
names(sampleLabel) <- colnames(DEGexp)

cogena_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

summary(cogena_result)

heatmapCluster(cogena_result, "k", "4",  add2=FALSE, 
               heatmapcol=colorRampPalette(c("blue", "white", "firebrick2")),
               cexCol=1.2, maintitle="cogena: co-expression analysis" )
heatmapPEI(cogena_result, "k", "4", maintitle="cogena: KEGG pathway analysis", 
           add2=FALSE)


save.image("../results/24h/Enrichment_24h_cogena_KEGG.RData")


################################################################################
save.image("../results/24h/Enrichment_24h.RData")



