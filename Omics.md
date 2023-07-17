---
name: Omics
topic: Genomics, Proteomics, Metabolomics, Transcriptomics, and Other Omics
maintainer: Julie Aubert, Toby Dylan Hocking, Nathalie Vialaneix
email: julie.aubert@inrae.fr
version: 2023-04-03
source: https://github.com/cran-task-views/Omics/
---

In this task view, we focus on
the most-important CRAN packages dedicated to omics data analysis including those related to annotation and databases. The term omics concerns various biological disciplines ending with the suffix -omics, such as genomics, proteomics, metabolomics, transcriptomics. 
We include CRAN packages which have been published more than one year ago and are regularly updated. Therefore, very recent packages may not appear yet but may be included soon if they are maintained. A few additional packages from [the Bioconductor project](https://www.bioconductor.org/) are
included, not looking for an exhaustive description but pointing to the must-see
packages or to packages focused on omics not well covered by CRAN packages. 
Complementary information might also be found in `r view("Phylogenetics")`, 
`r view("MachineLearning")`, `r view("Agriculture")` or `r view("MissingData")` task views. Note that packages covering
statistical genetics data (DNA data, haplotype estimation, population structure,
genetic epidemiology), phylogenetics and phylogenomics are not covered by this task view.
Still further packages might be found on
[Bioconductor](https://www.bioconductor.org/),
[R-Forge](https://r-forge.r-project.org/), and [GitHub](https://github.com/).

If you think we have missed some important packages in this list, please e-mail
the maintainers or submit an issue or pull request in the GitHub repository 
linked above. 

The task view is structured into main topics:

-   [Annotation and databases](#annotation-and-databases)
-   [Genomics](#genomics)
-   [Transcriptomics](#transcriptomics)
-   [Proteomics](#proteomics)
-   [Metabolomics](#metabolomics)
-   [Other omics](#other-omics)
-   [Multiple omics](#multiple-omics)
-   [Specific tasks](#specific-tasks)
-   [Specific application fields](#specific-application-fields)



### Annotation and databases

-   *Databases*: `r pkg("geneExpressionFromGEO")` can be used to easily download
    a gene expression dataset from [GEO](https://www.ncbi.nlm.nih.gov/geo/).
    `r pkg("cbioportalR")` makes available clinical and genomic data from
    [cBioPortal](http://www.cbioportal.org/). `r pkg("pinfsc50")` provides 
    genomic data for the plant pathogen *Phytophthora infestans*.
-   *Functional annotation*: `r pkg("enrichR")` provides an interface to the
    [Enrichr](https://maayanlab.cloud/Enrichr/) databases. 
    `r pkg("biomartr", priority = "core")` provides  an interface to the 
    [BioMart](https://www.ensembl.org/info/data/biomart/) database. 
    `r pkg("msigdbr")` provides an interface for the 
    [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) for Gene Set 
    Enrichment Analysis. `r pkg("gprofiler2")` provides an interface for 
    [g:Profiler](https://biit.cs.ut.ee/gprofiler/) (including functional 
    enrichment analysis, identifier conversion, or orthology across species).\
    `r pkg("toprdata")` provides gene and exon information from Ensembl genome 
    builds. `r pkg("read.gb")` can read records with `.gb` extension from the 
    NCBI Nucleotide database. `r pkg("Map2NCBI")` provides information on
    markers described by their positions by querying the 
    [NCBI](www.ncbi.nlm.nih.gov/) database.  `r pkg("PSSMCOOL")` contains the 
    computation of various features from Position Specific Scoring Matrix
    (PSSM).
-   *Gene symbol conversions*: `r pkg("aliases2entrez")` converts human gene 
    symbols to curated gene entrezID from [NCBI](https://www.ncbi.nlm.nih.gov/) 
    database. `r pkg("AnnotationBustR")` extracts subsequences into FASTA files  
    from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) annotations where gene 
    names may vary among accessions. `r pkg("BED")` is an interface for the 
    [Neo4j](https://neo4j.com/) database providing a mapping between different 
    identifiers of biological entities. More specifically, 
    `r pkg("riceidconverter")` can be used to convert rice biological 
    identifiers.\
    `r pkg("babelgene")` converts between human and non-human gene 
    orthologs/homologs from multiple databases. 
-   *Gene set/pathway enrichment*: `r pkg("WebGestaltR", priority = "core")`
    uses [WebGestalt](http://www.webgestalt.org/) to perform Gene Set 
    Enrichment and Network Topology analysis. `r pkg("GOxploreR")` provides a
    method to explore the gene ontology. `r pkg("liger")` and `r pkg("GSA")` 
    also propose methods for gene set enrichment analysis. `r pkg("VAM")` 
    proposes a gene set testing method that is better designed than standard 
    ones to handle single-cell RNA-seq data. `r pkg("tmod")` provides functions 
    for gene set enrichment analysis in transcriptomic and metabolic profiling 
    data. `r pkg("sigora")` corrects over-representation biases in pathway 
    enrichment analysis. `r pkg("permPATH")` performs pathway enrichment 
    analysis with permutation tests. `r pkg("CePa")` aims to find significant 
    pathways through network topology information and `r pkg("pathfindR")` 
    (accompanied by the data package `r pkg("pathfindR.data")`) and
    `r pkg("netgsa")` also use network information for enrichment analysis. 
    `r pkg("SlaPMEG")` performs pathway enrichment analysis for longitudinal 
    omics data. `r pkg("topologyGSA")` performs gene expression data tests using 
    given pathway information on genes. `r pkg("diffEnrich")` compares
    functional enrichment between two experimentally-derived groups of genes or
    proteins. `r pkg("DysPIA")` identifies dysregulated pathways based on a 
    pre-ranked gene pair list, using the `r pkg("DysPIAData")`.
    
### Genomics

-   *Sequence manipulation and analysis*: 
    `r pkg("GenomicTools.fileHandler", priority = "core")` is a collection of 
    I/O tools for handling the most commonly used genomic datafile types.
    `r pkg("valr")` can be used to read and manipulate genome intervals and 
    signals. `r pkg("KRIS")` provides various routine functions for 
    bioinformatics analyses of sequences. `r pkg("bioseq")` is a toolbox
    for manipulating biological (DNA, RNA and amino acid) sequences. 
    `r pkg("ftrCOOL")` extracts features from nucleotide and peptide sequences
    and converts them to discrete numbers in order to be used as predictors in 
    machine learning models.\
    `r pkg("agvgd")` predicts "missense" substitutions based on the properties of amino acid
    side chains and protein multiple sequence alignments. `r pkg("RNAsmc")` 
    provides functions to mine, compare and plot RNA secondary structure. 
    `r pkg("vhica")` can be used to detect horizontal transfers of transposable
    elements from their divergence compared with regular genes.\
    `r pkg("AnaCoDa")` contains a collection of models to analyze genome scale
    codon data using a Bayesian framework. `r pkg("apex")` contains a collection
    of tools for the analysis of aligned DNA sequences from multiple genes. 
    `r pkg("BASiNET")` implements a method to classify RNA Sequences using 
    network theory. `r pkg("MINTplates")` is a framework dedicated to
    exploration and annotation of tRNA and other genomic sequences. 
    `r pkg("cumSeg")` estimates the number and location of change points in 
    mean-shift (piecewise constant) models, such as genomic sequences.
-   *Cell analyses*: `r pkg("sonicLength")` estimates the abundance of cell 
    clones from DNA fragments. Similarly, `r pkg("scoper")` performs spectral 
    clustering to identify clones from B cell data. `r pkg("Signac")` is a
    framework for the analysis and exploration of single-cell chromatin data.
    `r pkg("popPCR")` allows the classification of droplets in fluorescence
    activated cell sorting data to estimate DNA target concentration.
-   *Copy number data*: `r pkg("ioncopy")` implements calculation of 
    copy-numbers in amplicon sequencing data. `r pkg("aroma.cn")` implements 
    several methods for normalizing and analyzing  DNA copy-number data and 
    `r pkg("PSCBS")` focuses on the analysis of parent-specific DNA copy-number 
    data.
-   *Mutation and CRISPR*: `r pkg("sigminer")` computes alteration signatures
    from genomic alteration records, `r pkg("Rediscover")` identifies mutually 
    exclusive mutations using a Poisson-Binomial model, and 
    `r pkg("mutSignatures")` identifies mutation signatures from somatic 
    mutational catalogs. More generally, `r pkg("ICAMS")` analyses mutational 
    signatures. Also, based on mutation data, `r pkg("ICBioMark")` identifies 
    immunotherapy biomarkers. `r pkg("pathwayTMB")` estimates tumor mutational
    burden (TMB) with a new pathway-based gene panel. `r pkg("CB2")` provides
    functions for CRISPR pooled screen analysis using Beta-Binomial Test.

### Transcriptomics

#### Microarray or other continuous expression data

-   *Quality control and normalization*: `r pkg("POD")` computes the 
    probability of detection (POD) curve and the limit of detection (LOD). 
    `r pkg("NACHO")` and `r pkg("nanostringr")` provide tools for quality 
    control or normalization of NanoString nCounter data. `r pkg("MiRNAQCD")` 
    is a toolbox for QC control of miRNA expression data.
-   *Integrated tools (with GUI)*: `r pkg("maGUI")` provides a graphical user 
    interface to analyze microarray data (including annotation, tests, or 
    network inference). `r pkg("depthTools")` (and its 
    [R commander](https://socialsciences.mcmaster.ca/jfox/Misc/Rcmdr/) 
    plugin `r pkg("RcmdrPlugin.depthTools")`) is a collection of statistical
    tools based on data depth for gene expression analysis.
    [OOMPA](http://oompa.r-forge.r-project.org/) provides a collection of 
    CRAN packages for microarray and proteomics analysis.
-   *Visualization*: `r pkg("tinyarray")` is dedicated to the visualization of 
    GEO and TCGA expression data. `r pkg("mpm")` provides exploratory graphical 
    analysis of gene expression data with various factorial approaches 
    (including PCA).
-   *Differential analysis*: `r bioc("limma")` is among the most-important 
    packages for microarray differential analysis. In addition, 
    `r pkg("TailRank")` provides a tail-rank non
    parametric test for microarray datasets. `r pkg("SMVar")` implements 
    structural model for variances for differential analysis of gene expression 
    data. `r pkg("RNentropy")` implements a method based on information theory 
    to detect significant variation in gene expression. `r pkg("optBiomarker")` 
    estimates the optimal number of biomarkers in 2-group assays. 
    `r pkg("leapp")` proposes a latent effect adjustment method to counter
    their effects on the rankings of hypotheses in tests. `r pkg("MCMC.qpcr")` 
    implements Bayesian approaches for normalization and differential analysis 
    of qPCR data. From a more global point of view, `r pkg("directPA")` can
    identify combinatorial effects of multiple treatments/conditions on pathways
    from different omics, including microarray and RNA-seq data.\
    `r pkg("ssize.fdr")` provides functions to calculate appropriate sample 
    sizes for gene expression tests and `r pkg("NewmanOmics")` implements Newman
    studentized range statistics for genome scale transcriptomics.
-   *Time series*: `r pkg("TcGSA")` and `r pkg("TGS")` implement methods for 
    longitudinal gene-expression data analysis. `r pkg("survival666")` 
    implements a method to eliminate the influence of co-expressed genes in 
    survival analyses. `r pkg("GeneCycle")` also analyses gene expression time 
    series to detect periodically expressed genes. `r pkg("BClustLonG")` 
    implements a Dirichlet process mixture model for clustering longitudinal 
    gene expression data.
-   *Clustering*: `r pkg("slfm")` performs gene expression analysis with a 
    Bayesian sparse latent factor model. `r pkg("lmQCM")` implements a 
    graph-based method for gene co-expression module discovery.
-   *Prediction*: `r pkg("LPS")` implements the linear prediction score
    approach for gene expression signatures. `r pkg("SMDIC")` performs the 
    identification of somatic mutation-driven immune cells from expression data.
    `r pkg("cubfits")` uses sequence and expression data to predict protein 
    production rates.

#### RNA-seq data

-   *Generic*: `r bioc("DESeq2")` and `r bioc("edgeR")` are among the
    most-important packages for RNA-seq preprocessing and differential analysis.
-   *Mapping*: `r pkg("MAAPER")` assigns 3' RNA-seq reads to polyA sites. 
-   *Deconvolution*: `r pkg("BisqueRNA")` and `r pkg("InteRD")` provide methods 
    to estimate cell type abundances from bulk expression data. `r pkg("scBio")`
    also proposes a deconvolution algorithm for bulk RNA-seq data but requires 
    one single-cell reference sample. More specifically, `r pkg("imsig")` 
    estimates the abundance of immune subpopulation cells in solid tumours.
    Finally, `r pkg("omicwas")` tests association with phenotypes that are 
    cell-type specific in bulk omics experiments.  
-   *Generic tool*: `r pkg("Tmisc")` is a collection of utility functions to 
    manipulate gene expression data.
-   *Simulation*: `r pkg("seqgendiff")` provides a framework to simulate 
    RNA-seq data under various assumptions and `r pkg("SeqNet")` and 
    `r pkg("graphsim")` both offer simulations of RNA-seq data based on 
    regulatory networks.
-   *Differential analysis*: 
    `r pkg("PQLseq")` implements a mixed model to account for 
    population structure in RNA-seq (and other count data) differential 
    analyses and `r pkg("glmmSeq")` also uses a mixed model for repeated 
    measurements in transcriptomic datasets. `r pkg("NBPSeq")` and 
    `r pkg("NBBttest")` are two packages that use (respectively) Negative 
    Binomial and Negative Binomial Beta tests for differential analyses of
    RNA-seq data. `r pkg("DiPALM")` enables differential analysis of
    time-course gene expression in different conditions. `r pkg("RVA")` is 
    dedicated to the visualization of RNA-seq data and especially of the results
    of differential analysis.\
    Similarly, `r pkg("markerpen")` uses penalized PCA for biomarker discovery.
    After differential analysis, `r pkg("DGEobj")` and  `r pkg("DGEobj.utils")` 
    provide a flexible container and a function toolkit to manage and annotate 
    Differential Gene Expression (DGE) analysis results.
-   *Clustering*: `r pkg("HTSCluster")` contains a model based on Poisson 
    mixture to cluster RNA-seq datasets. 
-   *Misc*: `r pkg("SIBERG")` implements a method to identify bimodally 
    expressed genes. `r pkg("CeRNASeek")` provides several functions to 
    identify and analyse miRNA sponge.

#### Single-cell RNA-seq

-   *Generic tools*: `r pkg("Seurat", priority = "core")`, `r pkg("iCellR")` 
    and `r pkg("pagoda2", priority = "core")` contain a collection of functions 
    for single-cell data analysis, including differential analysis (the first 
    one is based on the structure of `r pkg("SeuratObject")`). 
    `r pkg("sccore")` and `r pkg("singleCellHaystack")` also provide methods 
    for single-cell differential analysis, the second being based on KL 
    divergence.
-   *Quality control*: `r pkg("SoupX")` provides a method to remove mRNA 
    contamination from single-cell RNA-seq data.
-   *Cell clustering and annotation*: `r pkg("conos")` and `r pkg("scINSIGHT")` 
    can be used to identify recurrent cell clusters in collections of 
    single-cell RNA-seq datasets obtained in various conditions. 
    `r pkg("stochprofML")` models heterogeneity from populations of cells using 
    a mixture of log-normal distributions. `r pkg("ADAPTS")` constructs 
    cell-type signature matrices using flow sorted or single cell samples and 
    deconvolve bulk gene expression data. `r pkg("FiRE")` can discover rare 
    cells from voluminous single cell expression data.\
   `r pkg("scSorter")` assigns cells to known cell types according to marker
    genes and `r pkg("SignacX")` uses neural networks to identify cell types. 
    `r pkg("clustermole")` can be used to identify human and mouse single-cell 
    transcriptomic data cell type.\
    `r pkg("immunarch")` can be used to analyze specifically T-cells and 
    B-cells.
-   *Biomarker discovery*: `r pkg("DIscBIO")` is a user-friendly 
    multi-algorithmic pipeline for biomarker discovery in single-cell 
    transcriptomics. 
-   *Cell trajectories*: `r pkg("dynwrap")`, `r pkg("treefit")`, and
    `r pkg("SCORPIUS")` infer cell trajectories from single-cell gene 
    expression data. `r pkg("phateR")` can be used to visualize single-cell data
    with trajectories.
-   *Misc*: `r pkg("rPanglaoDB")` can download and merge single-cell RNA-seq 
    data from the [PanglaoDB](https://panglaodb.se/).

### Proteomics

-   *Generic tool*: `r pkg("ChemoSpec", priority = "core")` contains a 
    collection of functions for top-down exploratory chemometrics for various 
    proteomics spectroscopy data.
-   *Relation between proteins and DNA*: `r pkg("geno2proteo")` allows to find 
    the DNA and protein sequences of any given genomic loci using the 
    [ENSEMBL](https://www.ensembl.org/) annotations. 
    `r pkg("PredCRG")` contains a computational model to predict proteins 
    encoded by circadian genes.
-   *Protein structure and sequence*: `r pkg("compas")` manipulates and 
    analyzes 3-D structural geometry of Protein Data Bank (PDB) files. 
    `r pkg("bio3d")` contains tools to process, organize and explore protein 
    structure, sequence and dynamics data. `r pkg("ampir")` is a toolkit to 
    predict antimicrobial peptides from protein sequences on a genome-wide 
    scale. `r pkg("ptm")` contains functions for the analysis of 
    post-translational modifications in proteins. `r pkg("pbm")` contains 
    various models to analyze protein-ligand interactions.
-   *Mass spectrometry proteomics*: `r pkg("wrProteo", priority = "core")` 
    contains a collection of functions for the analysis of mass spectrometry 
    proteomic data. `r pkg("aLFQ", priority = "core")` implements the most common
    absolute label-free protein abundance estimation methods for LC-MS/MS and 
    `r pkg("iq")` also implements protein quantification for mass 
    spectrometry. `r pkg("protViz")` contains functions to visualize and analyze
    small mass spectrometry proteomic datasets.\
    `r pkg("PTXQC")` and `r pkg("ypssc")` are both designed to analyze outputs 
    of [MaxQuant](https://www.maxquant.org) (a quantitative proteomics tool for
    large mass-spectrometric datasets). `r pkg("protti")` can also analyze 
    outputs of MaxQuant in addition to outputs of 
    [Spectronaut](https://biognosys.com/software/spectronaut/) or of
    [Proteome Discoverer](https://www.thermofisher.com/fr/fr/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html).\ 
    `r pkg("prozor")` determines the minimal protein set explaining peptide 
    spectrum matches. `r pkg("net4pg")` implements a method to handle ambiguity
    for the identification in shotgun proteomics.\
    `r pkg("Inflect")` contains functions to compute melt temperature and shifts
    from LC-MS/MS abundance data in TPP experiments.
-   *Imputation*: Several packages implement functions for the imputation of
    proteomics datasets: `r pkg("imputeLCMD")` (left-censored imputation), 
    `r pkg("imp4p")` (including estimation of the missing value mechanism), and
    `r pkg("mi4p")` (multiple imputation). See also the `r view("MissingData")`
    view.
-   *Differential analysis*: `r pkg("canprot")` contains data files of published 
    differentially expressed proteins in cancer and cell culture experiments.
-   *Visualization*: `r pkg("cp4p")` provides calibration plots for proteomics 
    to check assumptions of FDR control procedures and to compute adjusted 
    p-values. `r pkg("LSPFP")` provides plots of peptides from shotgun 
    proteomics analysis.
-   *Prediction*: `r pkg("krm")` is a package dedicated to kernel based 
    regression methods that includes a kernel for protein sequence.
-   *Single-cell protein data (CITE-seq)*: `r pkg("dsb")` contains functions to 
    normalize and denoise droplet single-cell proteomics data.

### Metabolomics

-   *Data*: `r pkg("metaboData")` contains examples of various metabolomics 
    datasets.
-   *Generic tools*: `r pkg("omu")` provides metabolomics analysis tools 
    (including visualization) and `r pkg("lilikoi")` is a set of tools for 
    metabolomics, including normalization and various prediction methods. 
    `r pkg("MetabolomicsBasics")` also includes basic functions (peak picking, 
    deconvolution, ...) to explore metabolomics data.
-   *Identification*: `r pkg("erah")` can be used to identify metabolites by 
    spectral library matching in GC-MS untargeted metabolomics.
-   *Tests and biomarker discovery*: `r pkg("iCARH")` implements a method for 
    the analysis of time course metabolomic datasets, including biomarker 
    discovery. `r pkg("MetabolicSurv")` contains biomarker discovery and 
    survival analysis methods for metabolomics signatures. `r pkg("MetSizeR")` 
    computes sample size for target statistical power in metabolomics 
    experiments.
-   *Clustering*: `r pkg("RAMClustR")` includes a clustering algorithm for mass 
    spectrometry metabolomic data.
-   *Flux data*: `r pkg("fbar")` is a toolkit for flux balance analysis and 
    related metabolic modeling techniques.

### Other omics

-   *Methylation*: `r pkg("BiasCorrector")` is a GUI to correct measurement 
    bias in DNA methylation analyses. `r pkg("rBiasCorrection")` performs bias 
    correction in methylation data. `r pkg("TCA")` can deconvolve bulk 
    condition-specific DNA methylation data into condition and individual 
    specific methylation levels and detect associations with phenotypes.
-   *Hi-C*: `r bioc("InteractionSet")` is a Bioconductor package providing a 
    data structure adapted to Hi-C data. `r bioc("HiTC")` is a package to 
    manipulate and analyze Hi-C data and `r bioc("diffHic")` is focused on 
    differential analysis of Hi-C experiments.
-   *Lipidomics*: `r pkg("interep")` can perform interaction analysis in high 
    dimensional lipidomics datasets with repeated measurements. 
    `r pkg("flippant")` allows the analysis of the activity of the lipid 
    scrambling activity based on a fluorescence assay (dithionite scramblase 
    assay).
-   *Bacterial genomics, microbiome and metagenomics*: `r bioc("phyloseq")` is
    among the most-important packages for the analysis of high-throughput
    microbiome data. `r pkg("BarcodingR")` performs species identification
    using DNA barcodes. `r pkg("enveomics.R")` 
    contains a collection of functions for microbial ecology and other 
    applications of genomics and metagenomics and is the companion package for 
    the [Enveomics](http://enve-omics.ce.gatech.edu/enveomics/) collection.
-   *AIRR-seq / Rep-seq*: `r pkg("alakazam")` provides methods for
    high-throughput adaptive immune receptor repertoire sequencing
    (AIRR-Seq; Rep-Seq) analysis.

### Multiple omics

-   *Generic*: `r bioc("MOFA2")` and `r bioc("mixOmics")` are among the most two
    important packages for integration multi-omics data sets.
-   *Single omics analysis*: `r pkg("wrMisc")` contains a collection of tools 
    to manipulate omics data and to perform various statistical analyses 
    (including normalization and some statistical tests). `r pkg("integIRTy")` 
    provides a method to identify genes that are consistently altered in cancer 
    across different omics.
-   *Pathway analyses*: `r pkg("ActivePathways")` combines p-values to obtain
    enriched pathways and processes. `r pkg("ICDS")` identifies pathways 
    dysfunctional in cancer based on the integration of multiple omics.
-   *Exploratory integrative analyses*: `r pkg("o2plsda", priority = "core")` 
    provides functions to perform O2PLS-DA for multi-omics data integration. 
    Similarly, `r pkg("CovCombR")` can be used to combine heterogeneous data 
    sets through a covariance based method, and
    `r pkg("PMA", priority = "core")` proposes a sparse canonical correlation
    analysis to integrate multiple omics datasets. `r pkg("CovCombR")` and 
    `r pkg("packMBPLSDA")` contain other methods to combine heterogeneous data 
    sets through a covariance based method (the second being based on PLS-DA),
    and `r pkg("MOSS")` provides a sparse SVD-based multi-omics integration 
    method. `r pkg("semmcmc")` also provides an omics integration method based on
    structural equation modelling (SEM), and `r pkg("IMIX")` uses Gaussian
    mixtures for multi-omics data integration. `r pkg("solvebio")` is a binding 
    for the [SolveBio](https://www.solvebio.com/) API, a biomedical 
    knowledge hub oriented toward omics data integration. `r pkg("IntLIM")` 
    integrates two omics using linear modeling. More generally, `r pkg("dnet")` 
    performs integration from different angles including integration with 
    molecular networks, enrichments using ontologies, etc.\
    More specifically, `r pkg("iBATCGH")` integrates transcriptomic and CGH data
    with a Bayesian approach, and `r pkg("desiR")` provides functions for 
    ranking, selecting, and integrating genes, proteins and metabolite data.
-   *Prediction*: `r pkg("prioritylasso")` performs prediction from multiple
    omics using successive Lasso models using different priorities for the 
    omics.
-   *Clustering*: `r pkg("PINSPlus")` implements a robust clustering method
    for omics data integration and `r pkg("LUCIDus")` includes a ML based
    method for integrated clustering.
-   *Meta-analyses*: `r pkg("metaRNASeq")` and `r pkg("metaMA")`  implement
    p-value combination techniques for meta-analysis of RNA-seq data and 
    microarray data respectively. `r pkg("MetaIntegrator")` provides a pipeline
    for the meta-analysis of gene expression data.

### Specific tasks

#### Quality control and normalization

-   `r pkg("lineup")` and `r pkg("lineup2")` include tools for detecting and 
    correcting sample mix-ups in omics data. `r pkg("RUVIIIC")` performs 
    normalization using negative control variables and replications (originally 
    designed for proteomics).
    
#### Peak calling and analysis in data like ChIP-seq and ATAC-seq

-   `r pkg("PeakSegOptimal")` implements a change point detection method based 
    on the Poisson distribution for count data that includes a constraint suited
    for peak calling and `r pkg("PeakSegDisk")` provides a large scale 
    implementation of the method using on-disk storage. `r pkg("PeakError")` 
    computes true and false positive in peak calling with respect to annotated 
    region labels.
    
#### Multiple testing

-   `r pkg("bayefdr")` implements Bayesian estimation and optimization of 
    expected False Discovery Rate and `r pkg("hommel")` includes methods for 
    close testing with Simes' inequality. `r pkg("rSEA")` performs simultaneous 
    enrichment analysis that controls the FWER. More specifically, 
    `r pkg("RobustRankAggreg")` provides a method for aggregating rank lists, 
    especially lists of genes, `r pkg("DiscreteFDR")` implements multiple 
    testing procedures adapted for discrete tests.

#### High-dimensional data regularization

-   `r pkg("whitening")` implements whitening methods and CCA for 
    high-dimensional omics data. `r pkg("mpmi")` uses a kernel smoothing 
    approach for comparison of pairs of variables in large genomic datasets. 
    `r pkg("hsstan")` uses a hierarchical Bayesian approach for biomarker 
    discovery in high-dimensional datasets. `r pkg("supclust")` implements 
    methods for supervised clustering of potentially many predictor variables 
    (including 'PELORA' and 'WILMA').
    
#### Networks

-   *Gene network inference*: `r pkg("WGCNA", priority = "core")` implements 
    gene network inference with correlation based methods and 
    `r pkg("GeneNet", priority = "core")` implements gene network inference with
    Gaussian Graphical Models. `r pkg("RGBM")` implements bootstrap based 
    algorithms for network inference from microarray and RNA-seq data. 
    `r pkg("RNAseqNet")` implements network inference with a log-linear Poisson 
    model and can handle missing individuals. `r pkg("Patterns")` contains tools
    to infer 
    biological networks with approaches designed for single or multiple joint 
    omics. `r pkg("parmigene")` performs network inference with mutual 
    information methods. `r pkg("networkABC")` performs network inference with 
    Approximate Bayesian Computation. `r pkg("JSparO")` implements joint sparse 
    optimization for gene network inference for cell fate conversion. 
    `r pkg("Cascade")` implements a modeling tool allowing gene selection, 
    reverse engineering, and prediction in cascade networks with experimental 
    data provided in `r pkg("CascadeData")`. `r pkg("miic")` can be used to 
    infer causal and non causal networks using information theory.
-   *Single-cell data*: `r pkg("scLink")` and `r pkg("scTenifoldNet")` can be 
    used to infer co-expression networks from single-cell data, the latter 
    including an approach to compare these networks between different 
    conditions.
-   *Differential analysis*: `r pkg("scTenifoldKnk")` uses co-expression 
    networks to identify differentially regulated genes by a virtual knockout 
    approach and `r pkg("SIMMS")` enables integration of molecular profiles with
    functional networks (such as PPI networks) to detect biomarkers from 
    survival data. Similarly, `r pkg("regnet")` provides network-based 
    regularized models to perform variable selection in high-dimensional 
    biological data. More generally, `r pkg("DiffCorr")` implements a method for
    identifying pattern changes between 2 experimental conditions in correlation
    networks and `r pkg("dnapath")` integrates pathway information into the 
    differential network analysis of two gene expression datasets. 
    `r pkg("GANPA")` is a network-based gene weighting algorithm for pathway 
    enrichment analysis and `r pkg("DRaWR")` is a network-based method for 
    ranking genes or properties related to a given gene set.
-   *PPI networks*: `r pkg("prioGene")` can be used to define disease specific 
    PPI networks and to deduce candidate gene prioritization.

#### Clustering

-   `r pkg("EMMIXgene")` implements a mixture model-based approach for the 
    clustering of microarray expression data and `r pkg("ORIClust")` performs
    clustering of short time-course or dose-response microarray gene 
    expressions. `r pkg("OptCirClust")` performs clustering for circular data 
    (like circular DNA or RNA molecules) and `r pkg("adjclust")` performs
    constrained clustering adapted to genomic constraints. The latter 
    includes a wrapper for Hi-C datasets.

#### Visualization

-   *Sequence*: `r pkg("JBrowseR")` provides an R interface to the
    [JBrowse 2](https://jbrowse.org/jb2/) genome browser. `r pkg("chromoMap")` 
    provides interactive genomic visualization of the chromosomes or chromosome 
    regions of any living organism. Similarly, `r pkg("seqinr")` provides 
    visualizations for biological sequence (DNA and protein) data.
-   *Genome scale data*: `r pkg("genoPlotR")` produces various gene or genome 
    map figures ready for publications. Similarly, `r pkg("valr")` can be used
    to visualize genome-scale data, `r pkg("RCircos")` includes a collection of 2D 
    circos plots for genomic visualization, and `r pkg("RIdeogram")` provides 
    functions to display genome-wide data on ideograms. More specifically, 
    `r pkg("gggenes")` draws gene arrow maps with a
    [ggplot2](https://ggplot2.tidyverse.org/) approach. `r pkg("VALERIE")` 
    enables the visualization of alternative splicing events from single-cell 
    data. `r pkg("PACVr")` provides a function to visualize the coverage depth 
    of a complete plastid genome.
-   *Results of biomarker discovery methods*: `r pkg("func2vis")` provides 
    visualization and clean-up of enriched gene ontologies (GO) terms, protein 
    complexes and pathways using 
    [ConsensusPathDB](http://www.consensuspathdb.org/YCPDB). 
    `r pkg("statVisual")` implements visualizations for translational medicine 
    and biomarker discovery. Similarly, `r pkg("volcano3D")` provides 3D volcano
    and polar plots that is suited to visualize biomarker differential 
    analysis results for 3-class problems.
-   *Multiple omics*:  `r pkg("OmicNavigator")` is dedicated to the 
    visualization of omics data. `r pkg("wilson")` is a web-based tool dedicated 
    to the visualization of multi-omics data in an interactive way. 
-   *Misc*: `r pkg("BioInsight")` filters and plots the abundance of different 
    RNA biotypes present in a count matrix. 

#### Missing values

-   `r pkg("SurrogateRegression")` performs estimation and inference on a 
    partially missing target outcome while borrowing information from a 
    correlated surrogate outcome.

### Specific application fields

-   *Cancer*: `r pkg("oncoPredict")` provides biomarker discovery from cell line
    screening
    data. `r pkg("tidyestimate")` infers tumor purity from expression data.
    `r pkg("driveR")` is a tool for personalized or batch analysis of genomic 
    data for cancer driver gene prioritization by combining genomic information 
    and prior biological knowledge. `r pkg("DRomics")` is dedicated to
    omics data obtained using a dose-response design, with a large number of 
    tested doses.    
-   *Nutrition*: `r pkg("BRINDA")` computes biomarkers reflecting inflammation
    and nutritional determinants of anemia.

### Links

-   [The Bioconductor project](https://www.bioconductor.org/)
