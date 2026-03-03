# Model-agnostic variable importance with pathway-level aggregation to support reproducible translational oncology: application to cervical cancer

## Context 

Reliable prognostic biomarkers remain an unmet need in advanced cervical cancer (CESC), where clinical factors alone provide limited risk stratification. Although multi-omics profiling has identified numerous candidate molecular markers, their translational impact has been hindered by limited reproducibility, high dimensionality, and model-dependent inference. Pathway-level aggregation offers a biologically informed strategy to reduce feature redundancy and improve interpretability, while model-agnostic variable importance measures (VIMs) enable robust assessment of incremental predictive value in survival settings. We leverage two independent and complementary CESC cohorts—Bio-RAIDs and TCGA-CESC—to evaluate the stability and reproducibility of pathway-level prognostic signals across studies. This work aims to provide methodological insights into robust biomarker discovery and to illustrate how reproducibility-driven analytical strategies can support translational oncology research.

## Prerequisites 

This analysis requires `R (≥ 4.2.0)` and several survival machine learning packages.

1. Model-agnostic variable importance framework 

`survML` package is required for:

* Global survival stacking
* Estimation of model-agnostic variable importance

You can install a stable version of `survML` from CRAN using the following code: 

```r
install.packages("survML")
```

You can alternatively install the development version of the package from Github using the `devtools` package as follows: 

```r
## install.packages("devtools") # run only if necessary
install_github(repo = "cwolock/survML")
```

2. Additional required packages 

```r
install.packages(c(
  "survival",
  "tidyverse",
  "ranger",
  "xgboost", 
  "glmnet",
  "SuperLearner"
))
```

If using TCGA downloads: 
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
```

## Overall analytical workflow


## Data access 

1. TCGA-CESC cohort

Data from *The Cancer Genome Atlas (TCGA)* project CESC (Cervical Squamous Carcinoma and Endocervical Adenocarcinoma) are publicly available through the *Genomic Data Commons*.

Data can be accessed via: 

* The GDC Data portal: https://portal.gdc.cancer.gov 
* The R package `TCGAbiolinks`

In this repository, TCGA data are downloaded programmatically using `TCGAbiolinks` and saved locally in: 

data/tcga/raw


2. Bio-RAIDs cohort (NCT02428842)

The Bio-RAIDs cohort consists of clinically annotated cervical cancer patients frim a multicenter European study. 

* These data are **not publicly available**.
* Access is subject to institutional agreements and data-sharing policies.
* The data are not stored locally and **not distributed with this repository**.

## References 

**Bio-RAIDs cohort**:
```
@article{Ngo2015,
  title = {From prospective biobanking to precision medicine: BIO-RAIDs – an EU study protocol in cervical cancer},
  volume = {15},
  ISSN = {1471-2407},
  url = {http://dx.doi.org/10.1186/s12885-015-1801-0},
  DOI = {10.1186/s12885-015-1801-0},
  number = {1},
  journal = {BMC Cancer},
  publisher = {Springer Science and Business Media LLC},
  author = {Ngo,  Charlotte and Samuels,  Sanne and Bagrintseva,  Ksenia and Slocker,  Andrea and Hupé,  Philippe and Kenter,  Gemma and Popovic,  Marina and Samet,  Nina and Tresca,  Patricia and von der Leyen,  Heiko and Deutsch,  Eric and Rouzier,  Roman and Belin,  Lisa and Kamal,  Maud and Scholl,  Suzy},
  year = {2015},
  month = nov 
}
```

**TCGA-CESC**: 
```
@misc{Lucchesi2016,
  title     = "The Cancer Genome Atlas Cervical Squamous Cell Carcinoma and
               Endocervical Adenocarcinoma collection ({TCGA-CESC})",
  author    = "Lucchesi, Fabiano R and Aredes, Nat{\'a}lia D",
  abstract  = "The Cancer Genome Atlas Cervical Squamous Cell Carcinoma and
               Endocervical Adenocarcinoma (TCGA-CESC) data collection is part
               of a larger effort to build a research community focused on
               connecting cancer phenotypes to genotypes by providing clinical
               images matched to subjects from The Cancer Genome Atlas (TCGA).
               Clinical, genetic, and pathological data resides in the Genomic
               Data Commons (GDC) Data Portal while the radiological data is
               stored on The Cancer Imaging Archive (TCIA). Matched TCGA
               patient identifiers allow researchers to explore the TCGA/TCIA
               databases for correlations between tissue genotype, radiological
               phenotype and patient outcomes. Tissues for TCGA were collected
               from many sites all over the world in order to reach their
               accrual targets, usually around 500 specimens per cancer type.
               For this reason the image data sets are also extremely
               heterogeneous in terms of scanner modalities, manufacturers and
               acquisition protocols. In most cases the images were acquired as
               part of routine care and not as part of a controlled research
               study or clinical trial.",
  publisher = "The Cancer Imaging Archive",
  year      =  2016
}
```

**Conditional survival estimation using global survival stacking:**
```
@article{wolock2024framework,
        title={A framework for leveraging machine learning tools to estimate personalized survival curves},
        author={Wolock, Charles J and Gilbert, Peter B and Simon, Noah and Carone, Marco},
        journal={Journal of Computational and Graphical Statistics},
        year={2024},
        volume = {33},
        number = {3},
        pages = {1098--1108},
        publisher={Taylor \& Francis},
        doi={10.1080/10618600.2024.2304070}
}
```

**Model-agnostic variable importance in survival analysis:**
```
@article{wolock2025assessing,
    author = {Wolock, Charles J and Gilbert, Peter B and Simon, Noah and Carone, Marco},
    title = {Assessing variable importance in survival analysis using machine learning},
    journal = {Biometrika},
    volume = {112},
    number = {2}
    pages = {asae061},
    year = {2025},
    doi = {10.1093/biomet/asae061},
    publisher={Oxford University Press}
}

```