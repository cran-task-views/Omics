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
-   [Proteomics](#proteomics)
-   [Transcriptomics](#transcriptomics)
-   [Metabolomics](#metabolomics)
-   [Bacterial genomics, microbiome and metagenomics](#metagenomics)
-   [Integration of differents omics](#integration)
-   [Specific tasks](#tasks)
-   [Specific application fields](#applications)

Complementary information might also be found in `r view("Phylogenetics")`, `r view("StatisticalGenetics")`, `r view("MachineLearning")` `r view("Agriculture")`. Note that packages covering statistical genetics data (DNA data, haplotype estimation, population structure, genetic epidemiology), phylogenetics and phylogenomics are not covered by the Statistical Genomics task view.

If you think we have missed some important packages in this list, please e-mail the maintainers or submit an issue or pull request in the GitHub repository linked above.

[**Annotation**]{#annotation}

-   `r pkg("WebGestaltR")` uses [WebGestalt](http://www.webgestalt.org/) to 
    perform Gene Set Enrichment and Network Topology analysis.

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

-   *Visualization*: `r pkg("wilson")` is a web-based tool dedicated to the
    visualization of multi-omics data in an interactive way.

[**Specific application fields**]{#applications}

-   *Cancer*:

-   *Plant breeding*:

### Other ressources

-   [The Bioconductor project](https://www.bioconductor.org/)
-   [The PhyloGenetic CRAN task view](to%20be%20completed)
-   [The Epidemiology CRAN task view](to%20be%20completed)
-   [The StatisticalGenetics CRAN task view](to%20be%20completed)
