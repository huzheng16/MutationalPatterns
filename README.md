# MutationalPatterns

The MutationalPatterns R package provides a comprehensive set of functions for the extraction and plotting of mutational patterns in Single Nucleotide Variant (SNV) data.

## Installation

1. First install and load devtools package

```{r}
install.packages("devtools")
library(devtools)
```
2. Install and load MutationalPatterns package

```{r}
install_github("CuppenResearch/MutationalPatterns")
library(MutationalPatterns)
```

## Reference genome

1. List all available reference genomes (BSGenome)

```{r}
available.genomes()
```
2. Download and load your reference genome of choice

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
```