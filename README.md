# referenceNof1

The *referenceNof1* R package is designed to create robust reference standards using Jaccard Indices to calculate agreement and identify regions in which agreement is maximized among techniques for differential gene expression. 

## Methods

The *referenceNof1* R package uses the following methods to create a concordance matrix to evaluate whether the parameter filters are good or not: 

- DeSeq

- DESeq2

- edgeR

- NOISEq

- DEGseq

The concordance matrix is built using the Jaccard Index (JI) which is given by the following formula: 
 $$JI = \frac{|A\cap B|}{|A \cup B|}$$

## Data
The data provided in this R package is a subset of the MCF7 data (access: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51403), from the following manuscript: 

Liu Y, Zhou J, White KP. RNA-seq differential expression studies: more sequence or more replication? Bioinformatics 2014 Feb 1;30(3):301-4. PMID: 24319002


## Usage

The way referenceNof1 works is by identifying regions of Fold change and expression level cutoffs that decrease statistical noise and ensure greater agreement among techniques for differential gene expression. 

Measuring fold change and log fold change as the ratio difference
between a (S) stimulated and an unstimulated (U) sample: S/U, the referenceNof1 algorithm, takes in as input an expression-level cutoff and a fold-change cutoff, and returns a concordance matrix using JI as the similarity/agreement distance, as well as the median JI per method. Using the code below, we can apply an expression level cutoff of 30, and search for agreement in a fold change region of 1-1.3.

```
require(locfit)
mcf7 = referenceNof1::mcf7
res <- referenceNof1::applyFilter(mcf7, expressionCutoff = 30, fold_change_window = c(0,.3), 1)
```

The concordance matrix, provides the concordance or similarity scores across techniques, which can be accessed from the `res` object:

```
print(res$JaccardMatrix)
```


```{r echo=F}
print(res$JaccardMatrix)

```

and can be visualized into a concordance heatmap:

```
referenceNof1::concordance_heatmap(jaccard_matrix = res$JaccardMatrix)
```

## Identifying Optimal Reference Standard

To identify the region of optimal performance, the referenceNof1 R package has a function tha conducts a grid-search on user-selected parameters. For example, if a bioinformatician wants to determine the best cutoff threshold for lowly-expressed genes, a grid-search like below: 

```
referenceNof1::optimize_parameters(gs.mat = mcf7, expressionCutoffs = c(30,50), FCRegions = list(c(1, 1.5)),FDR.cutoff = 0.1)
```

where referenceNof1 will construct reference standards using a cutoff of 0 (i.e., no cutoff) and 30, and use those two thresholds to determine which one provides better concordance across techniques. 

# Install

To install the referenceNof1 R package, we invite you to install the latest/experimental release on GitHub, via the commands below, as well as the stable release on CRAN upon availibity.

## Github
install_github("samirrachidzaim/referenceNof1")

## Cran (coming shortly)

install.packages('referenceNof1')

