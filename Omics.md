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

-   `r pkg("toprdata")` provides gene and exon information from Ensembl genome
    builds.
-   `r pkg("valr")` can be used to read and manipulate genome intervals and
    signals.
-   `r pkg("vhica")` can be used to detect horizontal transfers of transposable 
    elements from their divergence compared with regular genes.
-   `r pkg("WebGestaltR")` uses [WebGestalt](http://www.webgestalt.org/) to 
    perform Gene Set Enrichment and Network Topology analysis. `r pkg("VAM")`
    proposes a gene set testing method that is better designed than standard 
    ones to handle single-cell RNA-seq data. `r pkg("tmod")`provides functions
    for gene set enrichment analysis in transcriptomic and metabolic profiling
    data.
-   `r pkg("topologyGSA")` performs gene expression data tests using given
    pathway information on genes.

[**Methylation**]{#methylation}

-   `r pkg("TCA")` can deconvolve bulk condition-specific DNA methylation data
    into condition and individual specific methylation levels and detect 
    associations with phenotypes.

[**Transcriptomics**]{#transcriptomics}

-   `r pkg("Tmisc")`is a collection of utility functions to manipulate gene
    expression data.
-   `r pkg("TailRank")`provides a tail-rank non parametric test for microarray
    datasets.
-   `r pkg("TwoPhaseInd")` estimates gene-treatment interactions in randomized
    clinical trials.
-   `r pkg("TcGSA")` and `r pkg("TGS")` implement methods for longitudinal 
    gene-expression data analysis.
-   `r pkg("treefit")` infers cell trajectories (as trees) from single-cell gene
    expression data.

[**Proteomics**]{#proteomics}

-   `r pkg("wrProteo", priority = "core")` contains a collection of functions 
    for the analysis of mass spectrometry proteomic data.
-   `r pkg("ypssc")` is designed to analyzed outputs of 
    [MaxQuant](https://www.maxquant.org)], which is a quantitative proteomics 
    tool for large mass-spectrometric datasets.

[**Integration**]{#integration}

-   `r pkg("wrMisc")` contains a collection of tools to manipulate omics data
    and to perform various simple statistical analyses (including normalization
    and some statistical tests).

[**Specific tasks**]{#tasks}

-   *Multiple testing*:

-   *High dimensional data - regularization*: `r pkg("whitening")` implements 
    whitening methods and CCA for high-dimensional omics data.

-   *Networks*:

-   *Visualization*: `r pkg("tinyarray")` is dedicated to the visualization of
    GEO and TCGA expression data.\ 
    `r pkg("valr")` can be used to visualize genome-scale data. 
    `r pkg("VALERIE")` enables visualization of alternative splicing event from
    single-cell data.\ 
    `r pkg("volcano3D")` provides 3D volcano and polar plots that is well suited
    to visualize biomarker differential analysis results for 3-class problems.\ 
    `r pkg("wilson")` is a web-based tool dedicated to the visualization of 
    multi-omics data in an interactive way. \ 

[**Specific application fields**]{#applications}

-   *Cancer*: `r pkg("tidyestimate")` infers tumor purity from expression data.

-   *Plant breeding*:

### Other ressources

-   [The Bioconductor project](https://www.bioconductor.org/)
-   [The PhyloGenetic CRAN task view](to%20be%20completed)
-   [The Epidemiology CRAN task view](to%20be%20completed)
-   [The StatisticalGenetics CRAN task view](to%20be%20completed)
