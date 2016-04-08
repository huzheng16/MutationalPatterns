# MutationalPatterns

The MutationalPatterns R package provides a comprehensive set of flexible functions for easy finding and plotting of mutational patterns in Single Nucleotide Variant (SNV) data.

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
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
  source("http://bioconductor.org/biocLite.R")
  biocLite(ref_genome)
  library(ref_genome, character.only = T)
  ```
  
### Load SNV data

Find package example data
  ```{r}
  vcf_files = list.files(system.file("extdata", package="MutationalPatterns"), full.names = T)
  ```

Load a single vcf file
  ```{r}
  vcf = read_vcf(vcf_files[1], "sample1")
  ```

Load a list of vcf files
  ```{r}
  sample_names = c("colon1", "colon2", "colon3", "intestine1", "intestine2", "intestine3", "liver1", "liver2", "liver3")
  vcfs = read_vcf(vcf_files, sample_names)
  ```

Include relevant metadata in your analysis, e.g. donor id, cell type, age, tissue type, mutant or wild type
  ```{r}
  tissue = c("colon", "colon", "colon", "intestine", "intestine", "intestine", "liver", "liver", "liver")
  ```

##  Analyses

### Mutation types

Retrieve base substitutions from vcf object as "REF>ALT"
  ```{r}
  get_muts(vcfs[[1]])
  ```
  
Retrieve base substitutions from vcf and convert to the 6 types of base substitution types that are distinguished by convention: C>A, C>G, C>T, T>A, T>C, T>G. For example, if the reference allele is G and the alternative allele is T (G>T), this functions returns the G:C>T:A mutation as a C>A mutation.
  ```{r}
  get_types(vcfs[[1]])
  ```
  
Retrieve the context (1 base upstream and 1 base downstream) of the positions in the vcf object from the reference genome.
  ```{r}
  get_mut_context(vcfs[[1]], ref_genome)
  ```

Retrieve the types and context of the base substitution types for all positions in the vcf object. For the base substitutions that are converted to the conventional base substitution types, the reverse complement of the context is returned.
  ```{r}
  get_type_context(vcfs[[1]], ref_genome)
  ```

Count mutation type occurences for one vcf object
  ```{r}
  type_occurences = count_type_occurences(vcfs[[1]], ref_genome)
  ```

Count mutation type occurences for all samples in a list of vcf objects
  ```{r}
  type_occurences = count_type_occurences(vcfs, ref_genome)
  ```

### Mutation spectrum

Plot mutation spectrum over all samples. Plottes is the mean relative contribution of each of the 6 base substitution types. Error bars indicate standard deviation over all samples. The n indicates the total number of mutations in the set.
  ```{r}
  plot_spectrum(type_occurences)
  ```

Plot mutation spectrum with distinction between C>T at CpG sites
  ```{r}
  plot_spectrum(type_occurences, CT = T)
  ```

Plot spectrum without legend
  ```{r}
  plot_spectrum(type_occurences, CT = T)
  ```

  ![spectra1](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/spectra1.png)

Specify 7 colors for spectrum plotting
  ```{r}
  my_colors = c("pink", "orange", "blue", "lightblue", "green", "red", "purple")
  plot_spectrum(type_occurences, CT = T, legend = T, colors = my_colors)
  ```

Plot spectrum for each tissue separately
  ```{r}
  plot_spectrum(type_occurences, by = tissue, CT = T)
  ```
  
  ![spectra2](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/spectra2.png)

### 96 Mutation Profile

### Extract Signatures

### Fit 96 mutation profiles to known signatures  


### Rainfall plot

A rainfall plot visualizes mutation types and intermutation distance. Rainfall plots can be used to visualize the distribution of mutations along the genome or a subset of chromosomes. The y-axis corresponds to the distance of a mutation with the previous mutation and is log10 transformed. Drop-downs from the plots indicate clusters or "hotspots" of mutations.

Make rainfall plot of all autosomal chromosomes
  ```{r}
  # define autosomal chromosomes
  chromosomes = seqnames(get(ref_genome))[1:22]
  # make rainfall plot
  plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes, cex = 1.5)

  ```

  ![rainfall1](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/rainfall1.png)
  
Make rainfall plot of chromosome 1

  ```{r}
  chromosomes = seqnames(get(ref_genome))[1]
  plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes[1], cex = 2)
  ```
  ![rainfall2](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/rainfall2.png)

  