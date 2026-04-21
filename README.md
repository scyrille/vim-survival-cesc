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

data/raw/tcga

2. Bio-RAIDs cohort (NCT02428842)

The Bio-RAIDs cohort consists of clinically annotated cervical cancer patients from a multicenter European study. 

* These data are **not publicly available**.
* Access is subject to institutional agreements and data-sharing policies.
* The data are **not distributed with this repository**.

## References 

**Bio-RAIDs cohort**:
```
@article{Ngo2015,
  title = {From prospective biobanking to precision medicine: BIO-RAIDs – an EU study protocol in cervical cancer},
  volume = {15},
  ISSN = {1471-2407},
  url = {http://dx.doi.org/10.1186/s12885-015-1801-0},
  doi = {10.1186/s12885-015-1801-0},
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
	title = {The {Cancer} {Genome} {Atlas} {Cervical} {Squamous} {Cell} {Carcinoma} and {Endocervical} {Adenocarcinoma} collection ({TCGA}-{CESC})},
	url = {https://www.cancerimagingarchive.net/collection/tcga-cesc/},
	doi = {https://doi.org/10.7937/K9/TCIA.2016.SQ4M8YP4},
	publisher = {The Cancer Imaging Archive},
	author = {Lucchesi, Fabiano R and Aredes, Natália D},
	year = {2016},
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