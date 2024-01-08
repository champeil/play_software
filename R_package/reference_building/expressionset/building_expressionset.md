# this tutorial is for building the expressionset according to [Biobase](https://www.bioconductor.org/packages/release/bioc/html/Biobase.html)
- author: laojp
- time: 2024.01.08
- position: SYSUCC bioinformatic platform

## introduction
- Biobase contains standardized data structures to represent genomic data. The ExpressionSet class is designed to combine several dierent sources of information into a single convenient structure
- expressionset object is mainly contain several parts:
  - assayData: expresssion data matrix
  - phenoData: dataframe to records the meta data
  - featureData: descripe the chips and technologies 
  - experimentData: describe the experiment data

## create the expressionset object
``` R
# the script is for constructing the expressionset objects
# necessary data
- assayData: expression data matrix with the row as gene names and col as sample names
- phenoData: dataframe to record each samples' pheno with the rownames are the same with assayData and col names are pheno name

# script
expressionset2 <- ExpressionSet(
  assayData = matrix,
  phenoData = new(
    "AnnotatedDataFrame",
    data=tibble(
      SampleID=colnames(matrix)
    ) %>%
      dplyr::left_join(clinical_data2,by=c("SampleID"="Tumor_Sample_Barcode")) %>%
      tibble::column_to_rownames(var="SampleID"),
        varMetadata=data.frame(labelDescription=c(colnames(clinical_data2)[-1]),row.names=c(colnames(clinical_data2)[-1]))
      ),
  annotation=set_2
)
```

## reference
- Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, Bravo HC, Davis S, Gatto L, Girke T, Gottardo R, Hahne F, Hansen KD, Irizarry RA, Lawrence M, Love MI, MacDonald J, Obenchain V, Ole's AK, Pag'es H, Reyes A, Shannon P, Smyth GK, Tenenbaum D, Waldron L, Morgan M (2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” Nature Methods, 12(2), 115–121. http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html.
