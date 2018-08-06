
library(ggfortify)

GE_raw <- read.csv("../data/rawdata_fromWenyanKang/AllSamples.GeneExpression.FPKM.tsv.csv")

GE_all <- GE_raw[,-1]
colnames(GE_all) <- gsub("_FPKM","", colnames(GE_all) )


pheno_allsample <- data.frame(sample_IDraw=character(54), sample_ID=character(54), stringsAsFactors=F)
pheno_allsample$sample_IDraw <- colnames(GE_all)
pheno_allsample$sample_ID <- gsub("_FPKM","", pheno_allsample$sample_IDraw)
pheno_allsample <- dplyr::left_join(pheno_allsample, pheno)

pheno_allsample$label <- factor(as.character(pheno_allsample$label), 
                                levels=c("C_0h", "C_2h", "C_6h",  "C_12h", "C_24h",  "C_48h",
                                         "F_2h", "F_6h", "F_12h", "F_24h", "F_48h"))

#color group and time
autoplot(prcomp(t(GE_all), scale=F), data=pheno_allsample, colour = "label", 
         size = 18, label = TRUE, label.size=4, label.colour="black", 
         jitter=TRUE, ts.colour="black" ) 


GE_PCA_raw <- read.csv("../data/rawdata_fromWenyanKang/T.PCA_resul.csv")
GE_PCA <- GE_PCA[,1:3] %>% dplyr::left_join(pheno_allsample, by=c("X"="sample_IDraw"))



GE_PCA1 <- tidyr::gather(GE_PCA, "Comp", "value", 2:3)
plot(GE_PCA[,2:3])

qplot(data=GE_PCA, x=Comp.1, y=Comp.2, color=label, size=8 )+ scale_colour_brewer(palette = "Set1")

ggplot(GE_PCA, aes(Comp.1, Comp.2, label=label)) + 
    geom_point(aes(colour = label), label=TRUE, size=9 ) +
    geom_text(aes(label=label) ) + 
    labs(x = "PC1(91.68%)", y="PC2(4.04%)")


ggplot(GE_PCA, aes(Comp.1, Comp.2, label=label)) + 
    geom_point(aes(colour = label), size=9 ) + 
    labs(x = "PC1(91.68%)", y="PC2(4.04%)")





