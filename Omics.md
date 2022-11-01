---
name: Omics
topic: Genomics
maintainer: Julie Aubert, Florian Privé and Nathalie Vialaneix 
email: 
version: 2022-08-25
source: 
---

In this task view, we focused on the most important CRAN packages, which have been published more than one year ago and are regularly updated. The task view is structured into main topics:

-   [Annotation and databases](#annotation)
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

[**Annotation and Databases**]{#annotation}

-   `r pkg("aliases2entrez")` converts human gene symbols to curated gene entrezID from NCBI database.

-   `r pkg("AnnotationBustR")` extracts subsequences into FASTA files from GenBank annotations where gene names may vary among accessions.

-   `r pkg("BED")` for Biological Entity Dictionary is an interface for the [Neo4j](https://neo4j.com/) database providing mapping between different identifiers of biological entities.

-   `r pkg("biomart")` provides  an interface to the 'BioMart' database and a standardized way to automate genome, proteome, 'RNA', coding sequence ('CDS'), 'GFF', and metagenome retrieval from 'NCBI RefSeq', 'NCBI Genbank', 'ENSEMBL', and 'UniProt' databases.
-   `r pkg("cbioportalR")` browse and query clinical and genomic Data from [cBioPortal](http://www.cbioportal.org/).
-   `r pkg("CePa")` aims to find significant pathways through network topology information.
-   `r pkg("toprdata")` provides gene and exon information from Ensembl genome builds.
 `r pkg("read.gb")` can read records with `.gb` extension form the NCBI Nucleotide database. `r pkg("babelgene")` converts between human and 
    non-human gene orthologs/homologs from multiple databases. 
    `r pkg("riceidconverter")` can be used to convert rice biological 
    identifiers. 
-   `r pkg("RNAsmc")` provides functions to mine, compare and plot RNA secondary
    structure.
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
    
[**Genomics**](#genomics)
*CRISPR*
-   `r pkg("crispRdesignR")` designs guide sequences for CRISPR/Cas9 genome editing.
-   `r pkg("CB2")` provided functions for CRISPR pooled screen analysis using Beta-Binomial Test.

*To be classified*
-   `r pkg("agvgd")` is an extension of the original 'Grantham' distance to multiple sequence alignments to predict "missense" based on the properties of amino acid side chains and protein multiple sequence alignments.

-   `r pkg("alakazam")` provides methods for for high-throughput adaptive immune receptor repertoire sequencing (AIRR-Seq; Rep-Seq) analysis.
-   `r pkg("AnaCoDa")` contains a collection of models to analyze genome scale codon data using a Bayesian framework.
-   `r pkg("apex")` contains a collection of tools for the analysis aligned DNA sequences from multiple genes.
-   `r pkg("aroma.cn")` implements several methods for normalizing and analyzing DNA copy-number data.
-   `r pkg("babelgene")` converts between human and non-human gene orthologs/homologs and integrates orthology assertion predictions sourced from multiple databases as compiled by the HGNC Comparison of Orthology Predictions.
-     `r pkg("BASiNET")` implements a method to classify RNA Sequences using Complex Network Theory.

-     `r pkg("bioseq")` is a toolbox for manipulating biological (DNA, RNA and amino acid) sequences.
- `r pkg("Cascade")` implements a modeling tool allowing gene selection, reverse engineering, and prediction in cascade networks. Some such experimental data are available in the `r pkg("CascadeData")`

-   `r pkg("cumSeg")` estimates the number and location of change points in mean-shift (piecewise constant) models, such as genomic sequences.

-   `r pkg("sigminer")` computes alteration signatures from genomic alteration 
    records and `r pkg("Rediscover")` identify mutually exclusive mutations 
    using a Poisson-Binomial model.
-   `r pkg("Signac")` is a framework for the analysis and exploration of 
    single-cell chromatin data.

[**Methylation**]{#methylation}

-   `r pkg("BiasCorrector")` is a GUI to correct measurement bias in DNA methylation analyses.

-   `r pkg("TCA")` can deconvolve bulk condition-specific DNA methylation data
    into condition and individual specific methylation levels and detect 
    associations with phenotypes.

[**Transcriptomics**]{#transcriptomics}
*Microarray or other continuous expression data*

-   `r pkg("TailRank")`provides a tail-rank non parametric test for microarray
    datasets. `r pkg("SMVar")` implements structural model for variances for
    differential analysis of gene expression data. `r pkg("RNentropy")`
    implements a method based on information theory to detect significant 
    variation in gene expression.  `r pkg("depthTools` (and its 
    [R commander](https://socialsciences.mcmaster.ca/jfox/Misc/Rcmdr/) 
    plugin `r pkg("RcmdrPlugin.depthTools")`) is a collection of statistical
    tools based on data depth for gene expression analysis.
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
-   `r pkg("BisqueRNA")` provides tools to accurately estimate cell type abundances from heterogeneous bulk expression based either on reference-based or marker-based methods.
-   `r pkg("Tmisc")` is a collection of utility functions to manipulate gene
    expression data.
-   `r pkg("seqgendiff")` provides a framework to simulate RNA-seq data under
    various assumptions and `r pkg("SeqNet")` offers simulations of RNA-seq data
    based on regulatory networks.
-   `r pkg("QuasiSeq")` includes function to perform differential analysis of
    RNA-seq data using quasi-Poisson or quasi-negative binomial models.
-   `r pkg("SIBERG")` implements a method to identify bimodally expressed genes.
-   `r pkg("CeRNASeek")` provides several functions to identify and analyse miRNA sponge.

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
- `r pkg("ADAPTS")` constructs cell-type signature matrices using flow sorted or single cell samples and deconvolve bulk gene expression data.  
-   `r pkg("scSorter")` assigns cell to known cell types according to marker
    genes and `r pkg("SignacX")` uses neural network to identify cell types. `r pkg("clustermole")` can be used to identify human and mouse single-Cell transcriptomic data cell type.
-   `r pkg("treefit")` and `r pkg("SCORPIUS")` infer cell trajectories from 
    single-cell gene expression data.
-   `r pkg("scBio")` proposes a deconvolution algorithm for bulk RNA-seq data 
    when one single-cell reference sample is available. `r pkg("scMappR")`
    assigns cell-type specificity score to DEGs obtained from bulk RNA-seq data
    by using a deconvolution method.
-   `r pkg("rPanglaoDB")` can download and merge single-cell RNA-seq data from
    the PanglaoDB https://panglaodb.se/

*???*
- `r pkg("cubfits")` estimates mutation and selection coefficients on synonymous codon bias usage based on models of ribosome overhead cost (ROC), and estimates and predicts protein production rates.
- `r pkg("bda")` implements algorithms for binned data analysis, gene expression data analysis and measurement error models for ordinal data analysis. 
- `r pkg("BClustLonG")` implements a dirichlet process mixture model for clustering longitudinal gene expression data.





    

[**Proteomics**]{#proteomics}

- `r pkg("ampir")` is a toolkit to predict antimicrobial peptides from protein sequences on a genome-wide scale.
-   `r pkg("canprot")`contains data files of published differentially expressed proteins in cancer and cell culture proteomics experiments and tools for calculate chemical metrics.
-   `r pkg("ChemoSpec", priority = "core")` contains a collection of functions for top-down exploratory chemometrics for spectroscopy.
-   `r pkg("compas")` manipulates and analyzes 3-D structural geometry of Protein Data Bank (PDB) files.
-   `r pkg("cp4p")` provides calibration plot for proteomics to check assumptions of FDR (false discovery rate) control procedures and to compute adjusted p-values.
-   `r pkg("bio3d")` contains utilities to process, organize and explore protein
    structure, sequence and dynamics data.
-   `r pkg("aLFQ")` implements the most common absolute label-free protein 
    abundance estimation methods for LC-MS/MS.
-   `r pkg("wrProteo", priority = "core")` contains a collection of functions 
    for the analysis of mass spectrometry proteomic data.
-   `r pkg("ypssc")` is designed to analyzed outputs of 
    [MaxQuant](https://www.maxquant.org)], which is a quantitative proteomics 
    tool for large mass-spectrometric datasets.
-   `r pkg("RPPASAPACE")` provides tools for the analysis of reverse-phase 
    protein arrays.

[**Metabolomics**]{#metabolomics}

-   `r pkg("RAMClustR")` includes a clustering algorithm for mass spectrometry
    metabolomic data.


[**Bacterial genomics, microbiome and metagenomics**]{#metagenomics}

-   `r pkg("BarcodingR")` performs species identification using DNA barcodes.
-   `r pkg("BacArena")` can be used for simulation of organisms living in 
    communities. 

[**Integration**]{#integration}

- `r pkg("ActivePathways")` uses p-value merging to combine gene- or protein-level signals, followed by ranked hypergeometric tests to determine enriched pathways and processes. 

- `r pkg("CovCombR")` can be used to combine heterogeneous data sets through a covariance based method.    

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

-   *Data normalization*: `r pkg("RUVIIIC")` performs normalization using
    negative control variables and replications.\
    `r pkg("rBiasCorrection")` performs bias correction in methylation data.

-   *Multiple testing*: `r pkg("bayefdr")` implements a bayesian estimation and optimisation of expected False Discovery Rate.\
`r pkg("rSEA")` performs simultaneous enrichment
    analysis that controls the FWER.\ 
    `r pkg("RobustRankAggreg")` provides a method for aggregating rank lists,
    especially lists of genes.

-   *High dimensional data - regularization*: `r pkg("whitening")` implements 
    whitening methods and CCA for high-dimensional omics data.\
    `r pkg("supclust")` implements methods for supervised clustering of
    potentially many predictor variables (including 'PELORA' and 'WILMA').

-   *Networks*:
    `r pkg("RGBM")` implements bootstrap based algorithm for network inference
    from microarray and RNA-seq data.
    `r pkg("scLink")` and `r pkg("scTenifoldNet")` can be used to infer 
    co-expression networks from single-cell data, the latter including an 
    approach to compare these networks between different conditions. \ 
    `r pkg("scTenifoldKnk")` uses co-expression networks to identify 
    differentially regulated genes by a virtual knockout approach and 
    `r pkg("SIMMS")` enables integration of molecular profiles with functional 
    networks (such as PPI networks) to detect biomarkers from survival data.
    Similarly, `r pkg("regnet")` provides network-base regularized models to 
    perform variable selection in high-dimensional biological data.

 

-   *Visualization*: 

 `r pkg("BioInsight")` filters and plots the abundance of different RNA biotypes present in a count matrix. 
    `r pkg("chromoMap")` provides interactive genomic visualization of the chromosomes or chromosome regions of any living organism.
    `r pkg("seqinr")` provides visualizations for biological sequence (DNA and 
    protein) data. \ 
    `r pkg("RVA")` is dedicated to visualization of RNA-seq data and especially 
    of results of differential analysis. Similarly, `r pkg("cRegulome")` builds 
    a 'SQLite' database and enables the visualization of Regulome-Gene 
    Expression Correlations in Cancer and `r pkg("tinyarray")` is dedicated to 
    the visualization of GEO and TCGA expression data.\ 
    `r pkg("valr")` can be used to visualize genome-scale data, 
    `r pkg("RCircos")` includes collection of 2D circos plots for genomic 
    visualization, and `r pkg("RIdeogram")` provides functions to display 
    genome-wide data on ideograms.\ 
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

-   *Nutrition*:
`r pkg("BRINDA")` computes a Biomarkers Reflecting Inflammation and Nutritional Determinants of Anemia (BRINDA) adjustment method from the BRINDA multi-agency and multi-country partnership.

### Other ressources

-   [The Bioconductor project](https://www.bioconductor.org/)
