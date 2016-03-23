# MutationalPatterns

The MutationalPatterns R package provides a comprehensive set of functions for the extraction and plotting of mutational patterns in Single Nucleotide Variant (SNV) data.

## Getting started

### Installation

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

### Reference genome

1. List all available reference genomes (BSgenome)

  ```{r}
  available.genomes()
  ```
2. Download and load your reference genome of interest

  ```{r}
  source("http://bioconductor.org/biocLite.R")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
  library("BSgenome.Hsapiens.UCSC.hg19")
  ```
  
### Load SNV data

Load a single vcf file
  ```{r}
  vcf = read_vcf("your_file.vcf")
  ```

Load a list of vcf files from directory
  ```{r}
  vcf_file_list = list.files("your_dir", full.names = T)
  sample_names = c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6")
  vcf_list = read_vcf_list(vcf_file_list, sample_names)
  ```

Include relevant metadata in your analysis, e.g. donor id, cell type, age, tissue type, mutant or wild type
  ```{r}
  sample_type = c("mutant", "mutant", "mutant", "wt", "wt", "wt")
  ```

##  Analyses

### Mutation types

Retrieve base substitutions from vcf object as "REF>ALT"
  ```{r}
  get_muts(vcf)
  ```
  
Retrieve base substitution from vcf and converted to the 6 types of base substitution types that are distinguished by convention: C>A, C>G, C>T, T>A, T>C, T>G. For example, if the reference allele is G and the alternative allele is T (G>T), this functions returns the G:C>T:A mutation as a C>A mutation.
  ```{r}
  get_types(vcf)
  ```
  
Retrieve the context (1 base upstream and 1 base downstream) of the positions in the vcf object from the reference genome.
  ```{r}
  get_mut_context(vcf, ref_genome)
  ```

Retrieve the types and context of the base substitution types for all positions in the vcf object. For the base substitutions that are converted to the conventional base substitution types, the reverse complement of the context is returned.
  ```{r}
  get_type_context(vcf, ref_genome)
  ```

Count mutation type occurences for one vcf object
  ```{r}
  type_occurences = count_type_occurences(list(vcf), ref_genome)
  ```

Count mutation type occurences for all samples in a list of vcf objects
  ```{r}
  type_occurences = count_type_occurences(vcf_list)
  ```

Plot mutation spectrum over all samples. Plottes is the mean relative contribution of each of the 6 base substitution types. Error bars indicate standard deviation over all samples. The n indicates the total number of mutations in the set.
  ```{r}
  plot_spectrum(type_occurences)
  ```

Plot mutation spectrum with distinction between C>T at CpG sites
  ```{r}
  plot_spectrum(type_occurences, CT = T)
  ```

Specify 7 colors for spectrum
  ```{r}
  colors = c("black", "red", "yellow", "orange", "green", "brown", "pink")
  plot_spectrum(type_occurences, CT = T)
  ```

Plot spectrum for each sample type separately
  ```{r}
  plot_spectrum(type_occurences_ercc1, by = sample_type, CT = T)
  ```
  ![spectrum](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/spectrum_per_type.pdf)

### Rainfall plot

Define the chromosomes you want to include in your rainfall plot
  ```{r}
  chromosomes = seqnames(BSgenome.Hsapiens.UCSC.hg19)[1:19]
  ```

Make rainfall plot
  ```{r}
  rainfall_plot(vcf_list[[1]], ref_genome = BSgenome.Hsapiens.UCSC.hg19, chromosomes = chromosomes)
  ```
  
  