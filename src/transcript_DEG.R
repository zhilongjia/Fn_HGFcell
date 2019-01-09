
library(dplyr)
gene_count_rpkm_fn <- dir("../data/rawdata_fromWenyanKang/Analysis_Report/Analysis_Report/BGI_result/Quantify/GeneExpression/GeneExpression",
                          pattern="transcript", full.names=TRUE)

gene_count_rpkm_fn <- gene_count_rpkm_fn[grep("12", gene_count_rpkm_fn)]

gene_count_list <- list()
for (fn_i in gene_count_rpkm_fn ) {
    # fn_i <- gene_count_rpkm_fn[1]
    gene_count_rpkm <- readr::read_tsv(fn_i) %>% 
        dplyr::select(transcript_id, expected_count) %>% 
        dplyr::group_by(transcript_id) %>%
        summarise_all(funs(sum))
    
    sample_id <- sub(".transcript.fpkm.xls","",basename(fn_i) )
    colnames(gene_count_rpkm)[2] <- sample_id
    
    gene_count_list[[sample_id]] <- gene_count_rpkm
}

sample_count_df <- Reduce(dplyr::inner_join, gene_count_list )


################################################################################
# Normlization
library(edgeR)
library(limma)

# All gene count noramlization
S_raw <- as.matrix(select(sample_count_df, -transcript_id) )
rownames(S_raw) <- sample_count_df$transcript_id
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

subsample_pheno <- dplyr::filter(pheno, timepoint=="12h")
Expdesign <- model.matrix(~subsample_pheno$individual + subsample_pheno$state)

subsample_dge <- dge[,subsample_pheno$sample_ID]

# DEA
v <- voom(subsample_dge, Expdesign)
Expfit1 <- lmFit(v, Expdesign)
Expfit2 <- eBayes(Expfit1)
GSE_limma <- topTable(Expfit2, coef=tail(colnames(Expdesign), 1), number=Inf)

write.table(rownames(GSE_limma), file="../results/GSE_limma_DEGs.txt", 
            row.names = F, col.names = F, quote=FALSE)

GSE_limma$refseq_raw <- rownames(GSE_limma)
GSE_limma$refseq_id <- sub("\\..*", "", rownames(GSE_limma) )
refseq_symbol_lib <- readr::read_tsv("../data/genesymbol_refseq_lib.txt")

refseq_limma <- dplyr::left_join(GSE_limma, refseq_symbol_lib, by=c("refseq_id"="RefSeq(supplied by NCBI)" ) )

length(which(!is.na(refseq_limma$`Approved Symbol`)))

save.image("~/tmp/transcript_DEG.RData")


library(biomaRt)
listFilters(mart=mart)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
refseq <- sub("\\..*", "", rownames(GSE_limma) )
tmp <- getBM(filters="refseq_mrna", attributes=c("refseq_mrna", "hgnc_symbol"), values=refseq, mart=mart)
tmp2 <- getBM(filters="refseq_mrna_predicted", attributes=c("refseq_mrna_predicted", "hgnc_symbol"), values=refseq, mart=mart)




