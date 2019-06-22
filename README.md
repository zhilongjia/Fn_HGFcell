## Reproducible research

The R code and data to reproduce all the figures and tables in the manuscript, *Time-course Transcriptome Analysis for Drug Repositioning in Fusobacterium nucleatum-infected Human Gingival Fibroblasts*.

## Directories and Files:

    data : 
        rawdata_fromWenyanKang/.../GeneExpression/ : the raw read counts of genes for all the samples.
        B0.gene.fpkm.xls : FPKM of genes for all the samples.
        genesymbol_refseq_lib.txt : the mapping relationship of gene symbol an refseqID.

    documents ：slides for the results. 

Set the src directory as the working dir to run the code. For example, `setwd("./src")`. change as necessary.

    src : the R code
        ./time_indep_DEA.R : the main pipeline, including
            Differential Expression Analysis. 
            get the annotation of the core genes (Table S2)
            Heatmap of Differentially expressed genes (Fig S1)
            PCA (Fig 2B) 
            veen diagram, (Fig 2A) 
            Pathway analysis of DEG at each timepoint (Fig 3) 
            GO Ontology analysis of DEG at each timepoint (Fig S2) 
            Coexpression haetmap analysis (Fig 4A) 
            cogena-based pathway analysis (Fig 4B) 
            cogena-based drug repositioning analysis (Fig 5ABC) 
        ./gene2KEGG.R : function to get the annotation of the input genes
        ./cogena cogena package. revised the ColSideColors design ( ncol=4, horiz=FALSE) for this paper.

        # results not shown in this paper, but perhaps, useful to others.
        ./PCA_allsamples.R : PCA of all samples using RPKM of gene.
        ./Enrichment_24h.R : An analysis for the data in the timepoint 24h
        ./time_course_DEA.R : an analysis using a timepoinit and the timeponit before it.
        ./time_course_GSEA_pathway_plot.R : GSEA pathway analysis.
        ./transcript_DEG.R : Differential transcript analysis.

## Update:

    Released as a Public repository at June 22, 2019.
    Created as a Private repository at Jan 9, 2019.

## Issue

Issuing question at [here](https://github.com/zhilongjia/Fn_HGFcell/issues) or [email](zhilongjia#gmail.com) me, thank you.

            