---
name: Omics
topic: Genomics
maintainer: Julie Aubert, Florian Privé and Nathalie Vialaneix 
email: 
version: 2022-08-25
source: 
---

In this task view, we focused on the most important CRAN packages, which have been published more than one year ago and are regularly updated. The task view is structured into main topics:

-   [Annotation](#annotation)
-   [Genomics](#genomics)
-   [Human genetic epidemiology](#human)
-   [Methylation](#methylation)
-   [Transcriptomics](#transcriptomics)
-   [Proteomics](#proteomics)
-   [Metabolomics](#metabolomics)
-   [Bacterial genomics, microbiome and metagenomics](#metagenomics)
-   [Integration of differents omics](#integration)
-   [Specific tasks](#tasks)
-   [Specific application fields](#applications)

Complementary information might also be found in `r view("Phylogenetics")`, `r view("StatisticalGenetics")`, `r view("MachineLearning")` `r view("Agriculture")`. Note that packages covering statistical genetics data (DNA data, haplotype estimation, population structure, genetic epidemiology), phylogenetics and phylogenomics are not covered by the Statistical Genomics task view.

If you think we have missed some important packages in this list, please e-mail the maintainers or submit an issue or pull request in the GitHub repository linked above.

[**Annotation**]{#annotation}

-   `r pkg("AnnotationBustR")` extracts subsequences into FASTA files from GenBank annotations where gene names may vary among accessions.

-   `r pkg("biomartr", priority = "core")` is an interface to the 
    [BioMart](http://www.ensembl.org/info/data/biomart/index.html) database and 
    a standardized way to automate data retrieval from several databases.
    `r pkg("BED")` is an interface for the [Neo4j](https://neo4j.com/) database 
    providing mapping between different identifiers of biological entities.
    `r pkg("toprdata")` provides gene and exon information from Ensembl genome
    builds.
-   `r pkg("valr")` can be used to read and manipulate genome intervals and
    signals.
-   `r pkg("vhica")` can be used to detect horizontal transfers of transposable 
    elements from their divergence compared with regular genes.
-   `r pkg("WebGestaltR", priority = "core")` uses 
    [WebGestalt](http://www.webgestalt.org/) to perform Gene Set Enrichment and 
    Network Topology analysis. `r pkg("VAM")` proposes a gene set testing method
    that is better designed than standard ones to handle single-cell RNA-seq 
    data. `r pkg("tmod")`provides functions for gene set enrichment analysis in 
    transcriptomic and metabolic profiling data. `r pkg("sigora")` corrects 
    over-representation biases in pathway enrichment analysis. 
    `r pkg("SlaPMEG")` performs pathway enrichment analysis for longitudinal
    omics data.
-   `r pkg("topologyGSA")` performs gene expression data tests using given
    pathway information on genes.
-   `r pkg("sonicLength")` estimates the abundance of cell clones from DNA
    fragments. Similarly, `r pkg("scoper")` perform spectral clustering to 
    identify clones from B cell data.

   
[**Genomics**]{#genomics}

-   `r pkg("apex")` contains a collection of tools for the analysis aligned DNA sequences from multiple genes.
-   `r pkg("aroma.cn")` implements several methods for normalizing and analyzing DNA copy-number data.
-   `r pkg("babelgene")` converts between human and non-human gene orthologs/homologs and integrates orthology assertion predictions sourced from multiple databases as compiled by the HGNC Comparison of Orthology Predictions.
-   `r pkg("crispRdesignR")` designs guide sequences for CRISPR/Cas9 genome editing.
-   `r pkg("cumSeg")` estimates the number and location of change points in mean-shift (piecewise constant) models, such as genomic sequences.
-   `r pkg("sigminer")` computes alteration signatures from genomic alteration 
    records.
-   `r pkg("Signac")` is a framework for the analysis and exploration of 
    single-cell chromatin data.

[**Methylation**]{#methylation}

-   `r pkg("TCA")` can deconvolve bulk condition-specific DNA methylation data
    into condition and individual specific methylation levels and detect 
    associations with phenotypes.

[**Transcriptomics**]{#transcriptomics}

*Microarray or other continuous expression data*

-   `r pkg("TailRank")`provides a tail-rank non parametric test for microarray
    datasets. `r pkg("SMVar")` implements structural model for variances for
    differential analysis of gene expression data.
-   `r pkg("ssize.fdr")` provides functions that calculate appropriate sample 
    sizes for gene expression tests.
-   `r pkg("TcGSA")` and `r pkg("TGS")` implement methods for longitudinal 
    gene-expression data analysis. `r pkg("survival666")` implements a method
    to eliminate influence of co-expressed genes in survival analyses.
-   `r pkg("slfm")` performs gene expression analysis with a Bayesian sparse 
    latent factor model.
-   `r pkg("SMDIC")` performs the identification of somatic mutation-driven 
    immune cells from expression data.

*RNA-seq*

-   `r pkg("Tmisc")`is a collection of utility functions to manipulate gene
    expression data.
-   `r pkg("seqgendiff")` provides a framework to simulate RNA-seq data under
    various assumptions and `r pkg("SeqNet")` offers simulations of RNA-seq data
    based on regulatory networks.
-   `r pkg("SIBERG")` implements a method to identify bimodally expressed genes.

*single-cell RNA-seq*

-   `r pkg("Seurat", priority = "core")` contains a collection of functions for 
    single-cell data analysis, including differential analysis that is based on
    the structure of `r pkg("SeuratObject")`. `r pkg("sccore")` and 
    `r pkg("singleCellHaystack")` also provide methods for single-cell 
    differential analysis, the second being based on KL divergence.
-   `r pkg("SoupX")` provides a method to remove mRNA contamination from 
    single-cell RNA-seq data.
-   `r pkg("conos")` and `r pkg("scINSIGHT")` can be used to identify recurrent 
    cell clusters in collections of single-cell RNA-seq datasets obtained in
    various conditions. `r pkg("stochprofML")` models heterogeneity from 
    populations of cells using a mixture of log-normal distributions.
-   `r pkg("scSorter")` assigns cell to known cell types according to marker
    genes and `r pkg("SignacX")` uses neural network to identify cell types.
-   `r pkg("treefit")` and `r pkg("SCORPIUS")` infer cell trajectories from 
    single-cell gene expression data.
-   `r pkg("scBio")` proposes a deconvolution algorithm for bulk RNA-seq data 
    when one single-cell reference sample is available. `r pkg("scMappR")`
    assigns cell-type specificity score to DEGs obtained from bulk RNA-seq data
    by using a deconvolution method.

*???*
- `r pkg("cubfits")` estimates mutation and selection coefficients on synonymous codon bias usage based on models of ribosome overhead cost (ROC), and estimates and predicts protein production rates.
- `r pkg("bda")` implements algorithms for binned data analysis, gene expression data analysis and measurement error models for ordinal data analysis. 





    

[**Proteomics**]{#proteomics}


-   `r pkg("aLFQ")` implements the most commonly used absolute label-free protein abundance estimation methods for LC-MS/MS modes quantifying on either MS1-, MS2-levels or spectral counts together with validation algorithms to enable automated data analysis and error estimation.
-   `r pkg("bio3d")`contains utilities to process, organize and explore protein structure, sequence and dynamics data.

-   `r pkg("wrProteo", priority = "core")` contains a collection of functions 
    for the analysis of mass spectrometry proteomic data.
-   `r pkg("ypssc")` is designed to analyzed outputs of 
    [MaxQuant](https://www.maxquant.org)], which is a quantitative proteomics 
    tool for large mass-spectrometric datasets.


[**Bacterial genomics, microbiome and metagenomics**]{#metagenomics}

-   `r pkg("BarcodingR")` performs species identification using DNA barcodes.\
-   `r pkg("BacArena")` can be used for simulation of organisms living in     communities. Each organism is represented individually and genome scale metabolic models determine the uptake and release of compounds.

[**Integration**]{#integration}

-   `r pkg("solvebio")` is a binding for the 
    [SolveBio](https://www.solvebio.com/) API that is a biomedical knowledge hub
    oriented toward omics data integration.
-   `r pkg("wrMisc")` contains a collection of tools to manipulate omics data
    and to perform various simple statistical analyses (including normalization
    and some statistical tests).
-   `r pkg("CovCombR")` can be used to combine heterogeneous data sets through a
    covariance based method. `r pkg("semmcmc")` also provides a omics 
    integration method based on structural equation modelling (SEM).



[**Specific tasks**]{#tasks}

-   *Multiple testing*:

-   *High dimensional data - regularization*: `r pkg("whitening")` implements 
    whitening methods and CCA for high-dimensional omics data.\
    `r pkg("supclust")` implements methods for supervised clustering of
    potentially many predictor variables (including 'PELORA' and 'WILMA').

-   *Networks*:
    `r pkg("scLink")` and `r pkg("scTenifoldNet")` can be used to infer 
    co-expression networks from single-cell data, the latter including an 
    approach to compare these networks between different conditions. 
    Similarly, `r pkg("scTenifoldKnk")` uses co-expression networks to identify
    differentially regulated genes by a virtual knockout approach and 
    `r pkg("SIMMS")` enables integration of molecular profiles with functional 
    networks (such as PPI networks) to detect biomarkers from survival data.

-   *Visualization*: 
    `r pkg("seqinr")` provides visualizations for biological sequence (DNA and 
    protein) data.
    `r pkg("cRegulome")` builds a 'SQLite' database and enables the 
    visualization of Regulome-Gene Expression Correlations in Cancer and 
    `r pkg("tinyarray")` is dedicated to the visualization of GEO and TCGA 
    expression data.\ 
    `r pkg("valr")` can be used to visualize genome-scale data.\ 
    `r pkg("VALERIE")` enables visualization of alternative splicing event from
    single-cell data.\ 
    `r pkg("statVisual")` implements visulizations for translational medicine 
    and biomarker discovery. Similarly, `r pkg("volcano3D")` provides 3D volcano
    and polar plots that is well suited to visualize biomarker differential 
    analysis results for 3-class problems.\ 
    `r pkg("wilson")` is a web-based tool dedicated to the visualization of 
    multi-omics data in an interactive way. \ 

-   *Missing values*: `r pkg("SurrogateRegression")` performs estimation and 
    inference on a partially missing target outcome while borrowing information 
    from a correlated surrogate outcome.

[**Specific application fields**]{#applications}

-   *Cancer*: `r pkg("tidyestimate")` infers tumor purity from expression data.

-   *Plant breeding*:

### Other ressources

-   [The Bioconductor project](https://www.bioconductor.org/)
-   [The PhyloGenetic CRAN task view](to%20be%20completed)
-   [The Epidemiology CRAN task view](to%20be%20completed)
-   [The StatisticalGenetics CRAN task view](to%20be%20completed)
