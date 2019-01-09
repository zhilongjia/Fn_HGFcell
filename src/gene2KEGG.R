
################################################################################
# output gene_symbol kegg and GO
gene2KEGG <- function(gene_symbols, db="org.Hs.eg.db") {
    library(GO.db)
    library(KEGG.db)
    # library(reactome.db)
    # library(org.Rn.eg.db)
    library(org.Hs.eg.db)
    library(dplyr)
    
    # to do 20180613
    # species recog
    
    gene_kegg_GO_fulltable <- select(org.Hs.eg.db, keys=gene_symbols, columns = c("SYMBOL","GENENAME", "GO", "PATH"),
                                     keytype="SYMBOL", multiVals="list") %>% 
        dplyr::left_join(toTable(KEGGPATHID2NAME), c("PATH"="path_id")) %>% 
        dplyr::left_join(toTable(GOTERM)[c("go_id", "Term", "Ontology", "Definition")], c("GO"="go_id")) 
    
    gene_kegg_GO_fulltable <- gene_kegg_GO_fulltable[c("SYMBOL", "GENENAME", "PATH", "path_name", "GO", "Term", "Ontology", "Definition")]
    
    # collase duplicated items
    vec2string <- function(x) {
        x <- unique(na.omit(x))
        paste(x, collapse = "; ")}
    
    gene_kegg_GO_table <- dplyr::group_by(gene_kegg_GO_fulltable, SYMBOL) %>% 
        dplyr::summarise_all(funs(vec2string))
    
    # write.table(gene_kegg_GO_table, file="../result/AMS_DEA_Liver_DEG_KEGG_GO.txt", sep="\t", quote=FALSE, row.names = FALSE)
    
}

