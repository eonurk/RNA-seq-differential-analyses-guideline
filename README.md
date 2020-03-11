# RNA Seq Differential Analyses Guideline in R
I won't analyze any particular dataset in this tutorial, rather this is just an attempt to draw a roadmap for those who are interested in reproducible research. Any help to this cause is cordially appreciated!

    Note: I will use "%>%" (pipes) from dplyr package throughout this tutorial, so remember to import it first:
```r
library(dplyr) #install.packages("dplyr")
```

## Data Preprocessing
Before starting to analyze any dataset, you should make sure that the quality of it is okay and to do so a preprocessing step is generally required.

First of all, ensure that your raw count matrix is `k x (1+n)` where `k` is the number of genes and `n` is the number of samples that you have and plus `1` is the gene names (or ids). Then, you can re-order the genes according to their standard deviations, so that the most variant genes are at the top:

```r
# Don't include 1st column since it is the gene names (or ids)
count.matrix <- count.matrix [rev(order(apply(count.matrix[,-1], 1, sd))),]
```
and now you can delete the duplicated genes using:

```r
count.matrix <- count.matrix [!duplicated(count.matrix[,1]),]
```
The reason why we first re-order and then delete the duplicated genes is that we don't want to lose significant information by deleting the wrong, uninformative duplicated genes later. 
Now, we can set the row names of count matrix to gene names and delete the gene name column. Since there are no duplicated genes in the matrix anymore setting the gene name column to row names will not throw any errors.

```r
# Set row names to gene names (or IDs)
rownames(count.matrix) <- count.matrix[,1]

# Filter gene names
count.matrix <- count.matrix[,-1]
```
To prevent problems which may happen further down the pipeline, we should make sure that all counts are integers:
```r
# Enforce all counts to be integers
count.matrix <- round(count.matrix, 0)
```
Now that we have count matrix ready, we can filter the lowly expressed genes:
```r
# Filter low-expressed genes
# Keep the genes that have Count-Per-Million more than k = 0.5 in n = 1 libraries
# It is pretty similar to filterByExpr(y, min.count = 0.5) from edgeR [1] but its choice of n is different.
count.matrix <- count.matrix [rowSums(cpm(count.matrix) >= 0.5) >= 1,]
```
Note that you can be more strict about the filtering (higher k and n) depending on your dataset. 

At this point, we have a raw count matrix, which is filtered and sorted according to gene rows' stardard deviations. However, for some differential analyses methods (e.g. limma-trend) we may need normalized count matrix. To do that, we can use `DESeq2` [2] package but first we need a meta data matrix, which helps us to describe our particular problem. This meta data matrix is composed of the features that you want to compare (Status, Age, Gender etc.) per sample. So it is an `n x f` matrix where `n` is the number of samples that we defined as above and `f` is the number of features you want to compare. Having created meta data matrix, we can create our groups, and assign each sample that we have to these groups.

```r
# Create your desired groups
# Remember these parameters (Status, Age, Gender) depends on your own analyses!
group <- paste(meta.data$Status,meta.data$Gender ,sep="_")

# Assign each sample to its group assuming the column names of 
# count matrix are the sample names.
colData <- cbind(colnames(count.matrix), group)
colnames(colData)  = c("sample", "groups")
```
Now we can normalize our count matrix using `DESeq2` package, with `log2(cpm + c)`.
```r
library(DESeq2)

# Create DEseq Object and estimate factors to be used in normalization.
dds <- DESeqDataSetFromMatrix(countData = count.matrix, colData = colData, design = ~ groups)
dds <- estimateSizeFactors(dds)

# Normalize the counts using log2(cpm + c)
# c term is added in order to avoid log(0)
count.matrix.normalized <- log2(counts(dds, normalized=TRUE) + 4)
```
At this point we have both raw and normalized count matrix which could be used down the pipeline. 

One possible problem that you may encounter is the gene name convention complication. You may have the Ensembl Gene IDs or Gene Symbols in row names of your count matrix, however, you might need other or vice versa during the rest of the analyses. (or Entrez IDs)

## Gene Name Conversion
To change the gene naming convention, you can use `biomaRt` package [3]. The following example demonstrates the conversion from `Gene Ensembl ID` to `Gene Symbol`, but you can change it to any other by playing with the parameters:

```r
library(biomaRt)

# Load Biomart DB 
# Notice that this database is for human, to see different DBs:
# mart = useMart("ensembl"); listDatasets(mart)
mart <- useMart("ensembl", 
                host = "useast.ensembl.org",
                dataset = "hsapiens_gene_ensembl")

# Get the mapping matrix for ensembl_gene_id to hgnc_symbol!
mapping <- getBM(mart = mart, 
                 useCache = T,
                 uniqueRows = F,
                 filters = "ensembl_gene_id",
                 values = rownames(count.matrix),
                 attributes = c("ensembl_gene_id","hgnc_symbol"))

convert_Ensembl_to_GeneID <- function (mapping, count.matrix) {
    # Map Ensembl Gene IDs to Gene IDs
    count.matrix <- merge(mapping, count.matrix, 
                          by.x = "ensembl_gene_id", by.y = "row.names", all.y = T)
    
    # Use original Ensembl GeneIDs for non-converted genes
    ix <- which(count.matrix[,2] %in% "" | is.na(count.matrix[,2]))
    count.matrix[ix,2] <- count.matrix[ix,1]
    
    # Order genes according to their standard deviation in decreasing order
    count.matrix <- count.matrix [rev(order(apply(count.matrix[,c(-1,-2)], 1, sd))),]
    
    # Remove duplicated Gene IDs
    count.matrix <- count.matrix [!duplicated(count.matrix[,2]),]
    
    # Re-order genes according to their standard deviation in decreasing order
    count.matrix <- count.matrix [rev(order(apply(count.matrix[,c(-1,-2)], 1, sd))),]
    
    # Make row names Gene IDs
    rownames(count.matrix) <- count.matrix[,2] 
    
    # Filter Ensembl & Gene IDs and Return
    return(count.matrix[,c(-1,-2)])
}

# Convert row names to Gene IDs
count.matrix <- convert_Ensembl_to_GeneID(mapping, count.matrix)
count.matrix.normalized <- convert_Ensembl_to_GeneID(mapping, count.matrix.normalized)
```

    If you completed all the steps above, now you should have your raw and normalized count matrices ready. Yeey!

## Quality Control
The first thing you want to check now is the library sizes that you have, which is the sum of all counts of the genes per sample. This is because later this information will help us to determine the method that we want to use in our differential analyses. You can check your library sizes as follows:

```r
barplot(colSums(count.matrix)/1e6, 
        las = 3, 
        col = 'red', 
        main="Total read counts (millions)")  
```
To see whether the normalization helped eliminating the variation among library sizes, you can check out the distributions of the normalized library sizes:
```r
boxplot(count.matrix.normalized, 
        las  = 3, 
        col  = 'red',
        ylab = 'Normalized expression levels',
        main = 'Distribution of transformed data') 
```
## PCA Plots
PCA plots are suitable for checking whether any feature that you have explains your data in an obvious manner.
Using the meta data matrix that we created earlier, and `ggplot2` package, we can create pretty plots for different features that want to investigate. For instance, we can seperate the samples according to their genders & disease status and represent genders with different colors and disease status with various shapes:

```r
library(ggplot2)

# We don't need to scale or center anymore
# Samples should be in rows
pca <- prcomp(t(count.matrix.normalized))

d <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)
xl <- sprintf("PC 1: %.1f %%", d[1])
yl <- sprintf("PC 2: %.1f %%", d[2])

df <- data.frame(PC1=as.numeric(pca$x[,1]),
                 PC2=as.numeric(pca$x[,2]),
                 Gender = meta.data$Gender,
                 Status = meta.data$Status)

PCA.plot.age_status <- 
    ggplot(df, aes(PC1, PC2, color = Gender, shape = Status)) + 
    geom_point() +  
    labs(x=xl,y=yl) +  
    theme_minimal() + 
    coord_fixed(ratio = 1)
```

## Differential Analyses
We will run the same pipeline with 4 different methods: edgeR, limma-voom, limma-trend and DESeq2. So that you can choose any of them to conduct your own analyses.

First, let's define global thresholds for the analyses, which could be tweaked to your taste later.
```r
FDR.cutoff <- 0.1
LFC.cutoff <- 0
```

I will start with `edgeR` since its pipeline can be used for both `limma` methods as well.

### edgeR 
One thing to notice here is that edgeR needs **raw count matrix** to run the differential analyses, so actually we did not have to create the normalized count matrix to use it (sorry!). The following script creates an `DGEList` object, calculates the normalization factors (using TMM) which are used internally in `edgeR` pipeline, creates a design matrix from our `group` that we created earlier in Data Processing chapter, estimates dispersion among genes and fits a generalized linear model to our data with the design matrix. Yeah it does a lot of things...

```r
# Create DGEList object
y <- DGEList(counts = count.matrix)

# Calculate normalization factors for library sizes with TMM
y <- calcNormFactors(y, method = "TMM")

# Add intercept term for multiple comparisons
design <- model.matrix(~ 0 + group) 
rownames(design) <- colnames(count.matrix)
colnames(design) <- levels(group)

# Estimate dispersion for genes with Bayesian Shrinkage
y <- estimateDisp(y,design)

# Fit the model
fit.glm <- glmQLFit(y,design)
```

Now, all we have to do is to create our contrasts and start checking for important the genes! Keep in mind that this step entirely depends on the feature names that you selected. 

For instance let's say you have the following meta data matrix,

|          	| Status  	| Gender 	|
|----------	|---------	|--------	|
| Sample_1 	| Patient 	| Female 	|
| Sample_2 	| Control 	| Male   	|
| Sample_3 	| Control 	| Female 	|
| Sample_4 	| Patient 	| Male   	|
| Sample_5 	| Patient 	| Male   	|

Then, we can create some contrasts as follows (the order you defined the `group` is also important!):

```r
# Create contrasts for comparison!
contrasts <- makeContrasts(
    
    # Gender comparison for patients    
    MxF_Patient = Patient_Male - Patient_Female,
    
    # Status comparison for males
    PxC_Male = Patient_Male - Control_Male,
    
    levels = fit.glm$design
)
```
Having defined the `contrasts` we can test and create a differentially expressed gene list for each of them:
```r
# Create DE gene list for edgeR
DE.genes.edger <- list()

for (i in seq_len(ncol(contrasts))){
    contrast.name <- colnames(contrasts)[i]
    # Test for the contrast
    qlf <- glmQLFTest(fit.glm,contrast = contrasts[,i])
    # Get the list of all genes
    top.table <- topTags(qlf, n = Inf, p.value = FDR.cutoff)$table
    # Place them into the list
    DE.genes.edger[[contrast.name]] <- top.table[abs(top.table$logFC) >= LFC.cutoff,]
}
```
You can now check out `DE.genes.edger` to see differentially expressed genes for different comparisons!

### limma-trend
For, limma-trend we need normalized counts with our design matrix and while computing the statistics with `eBayes` function we must set `trend = T`:

```r
library(limma)

fit.trend <- lmFit(count.matrix.normalized, design)
fit.trend2 <- eBayes(contrasts.fit(fit.trend, contrasts), trend = T)

# Create DE gene list for limma-trend
DE.genes.trend <- list()

for (i in seq_len(ncol(contrasts))){
    contrast.name <- colnames(contrasts)[i]
    top.table <- topTable(fit.trend2, 
                          coef = colnames(contrasts)[i],
                          p.value = FDR.cutoff,
                          lfc = LFC.cutoff,
                          number = Inf)
    DE.genes.trend[[contrast.name]] <- top.table
}
```
### limma-voom
`voom` uses quantile normalization to eliminate library depth differences between replicates, therefore fits a slightly smoother curve to the mean-variance trend, compared to `limma-trend` pipeline. However, it does not mean that this will increase the number of true positive genes. So, if you have library depth differences, you can consider using voom pipeline.

To see the mean-variance plot you can set `plot=TRUE` in `voom` function.

```r
v <- voom(count.matrix, design, plot=F)
fit.voom <- lmFit(v, design)
fit.voom2 <- eBayes(contrasts.fit(fit.voom, contrasts))
# summary(decideTests(fit.voom2, method="separate", lfc = 0, p.value = 0.1))

# Create DE gene list for limma-voom
DE.genes.voom <- list()

for (i in seq_len(ncol(contrasts))){
    contrast.name <- colnames(contrasts)[i]
    top.table <- topTable(fit.voom2, 
                          coef = colnames(contrasts)[i],
                          p.value = FDR.cutoff, 
                          lfc = LFC.cutoff, 
                          number = Inf)
    DE.genes.voom[[contrast.name]] <- top.table
}
```

### DESeq2
`DESeq2` has a slightly different pipeline than above three. First, you assign the each sample to its corresponding group in a column. Then, you create an object with the raw count matrix and the column you created:

```r
# Create your desired groups
group <- paste(meta.data$Vaccine,meta.data$Response, meta.data$Day ,sep=".")

# Assign each sample to its group
colData <- cbind(colnames(count.matrix), group) %>% as.data.frame()
colnames(colData)  = c("sample", "groups")

# Create DEseq Object
dds <- DESeqDataSetFromMatrix(countData = count.matrix, colData = colData, design = ~ groups)
```
Then, you are ready to go:

```r
dds = DESeq(dds)
```
    Note: It may take longer compared to other methods, so saving the object may help for further analyses.
    
Then, we can just test the results as we did earlier,
```r
# Create DE gene list for DESeq2
DE.genes.deseq2 <- list()

for (i in seq_len(ncol(contrasts))){
    contrast.name <- colnames(contrasts)[i]
    DEseq.contrast <- rownames(contrasts)[contrasts[,i] != 0]
    res <- results(dds, c("groups", DEseq.contrast[2], DEseq.contrast[1]), tidy = T)
    rownames(res) <- res$row
    res.ordered <- res[order(res$pvalue),]
    res.significant <- subset(res.ordered, padj <= FDR.cutoff)
    res.significant <- subset(res.significant, abs(log2FoldChange) >= LFC.cutoff)
    DE.genes.deseq2[[contrast.name]] <- res.significant
}
```

## Conclusion
You can investigate `DE.genes.[method]` to see the differentially expressed genes for different methods. See whether they overlap for your analyses. 

Please note that, this is just a bad recap of all of these methods, they all have really nice user guides explaining all of these steps in more detailed and structured manner. So, if you want to learn more about them, please refer to their user guides.

## Acknowledgement

This tutorial is vastly inspired by [IDEP.90](http://bioinformatics.sdstate.edu/idep/) project. 
Check it out, it is very cool!

## Contribution

You can send pull requests to make your contributions.

I occasionally mess up, so all comments are appreciated!


## Author

- EO Karakaslar

## License

- GNU General Public License v3.0

## References

- [1] Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140. doi: 10.1093/bioinformatics/btp616.

- [2] Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

- [3] Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, 1184–1191.