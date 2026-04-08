
# ----------------------------- Documentation ------------------------------#

#' This file provides helpers for downloading, importing and processing
#' TCGA-CESC data.

# ----------------------------- Dependencies -------------------------------#

library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(janitor)
library(readr)
library(TCGAbiolinks)
library(valr)

#-------------------------------  Functions --------------------------------#

#' Download and import TCGA-CESC data 
#' 
#' @param project TCGA project (default "TCGA-CESC")
#' @param directory Local directory used by TCGAbiolinks (default "data")
#' @param download Logical. If TRUE, runs GDCdownload for each query (default FALSE)
#' @param files_per_chunk Passed to GDCdownload (default 1)
#' @param clinical_info Vector of clinical.info blocks for GDCprepare_clinic
#'
#' @return A named list with: tcga_cesc_clin, tcga_cesc_CNVs, tcga_cesc_SNVs, tcga_cesc_RNA
#' @export
fetch_tcga_raw <- function(
    project = "TCGA-CESC",
    directory = "data/raw/tcga",
    download = FALSE,
    files_per_chunk = 1,
    clinical_info = c("patient", "stage_event","drug","admin","radiation")
) {
  stopifnot(dir.exists(directory) || dir.create(directory, recursive = TRUE))
  
  # Clinical 
  query_clin <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Clinical",
    data.format = "BCR XML"
  )
  
  if (isTRUE(download)) {
    TCGAbiolinks::GDCdownload(
      query = query_clin,
      directory = directory,
      method = "api",
      files.per.chunk = files_per_chunk
    )
  }
  
  clin <- lapply(
    clinical_info,
    function(x)
      TCGAbiolinks::GDCprepare_clinic(
        query = query_clin,
        directory = directory,
        clinical.info = x
      )
  )
  names(clin) <- clinical_info
  
  # CNVs 
  query_CNVs <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Copy Number Segment",
    experimental.strategy = "WGS",
    workflow.type = "GATK4 CNV",
    data.format = "TXT",
    platform = "Illumina",
    access = "open"
  )
  
  if (isTRUE(download)) {
    TCGAbiolinks::GDCdownload(
      query = query_CNVs,
      directory = directory,
      method = "api",
      files.per.chunk = files_per_chunk
    )
  }
  
  CNVs <- TCGAbiolinks::GDCprepare(
    query_CNVs,
    directory = directory
  )
  
  # SNVs (MAF)
  query_SNVs <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    experimental.strategy = "WXS",
    workflow.type =
      "Aliquot Ensemble Somatic Variant Merging and Masking",
    data.format = "MAF",
    platform = "Illumina",
    access = "open"
  )
  
  if (isTRUE(download)) {
    TCGAbiolinks::GDCdownload(
      query = query_SNVs,
      directory = directory,
      method = "api",
      files.per.chunk = files_per_chunk
    )
  }
  
  SNVs <- TCGAbiolinks::GDCprepare(
    query_SNVs,
    directory = directory
  )
  
  # RNA 
  query_RNA <- TCGAbiolinks::GDCquery(
    project = project,
    data.type = "Gene Expression Quantification",
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.format = "TSV",
    platform = "Illumina",
    access = "open",
    sample.type = "Primary Tumor"
  )
  
  if (isTRUE(download)) {
    TCGAbiolinks::GDCdownload(
      query = query_RNA,
      directory = directory,
      method = "api",
      files.per.chunk = files_per_chunk
    )
  }
  
  RNA <- TCGAbiolinks::GDCprepare(
    query_RNA,
    directory = directory,
    summarizedExperiment = FALSE
  )
  
  # Output 
  list(
    clin = clin,
    CNVs = CNVs,
    SNVs = SNVs,
    RNA = RNA
  )
  
}

#' Process TCGA CNV segment data and derive gene-level CNV status
#'
#' This function preprocesses TCGA copy-number variation (CNV) segment data,
#' optionally generates diagnostic histograms and sample-level CNV profile
#' plots, and derives gene-level CNV status by intersecting segment coordinates
#' with gene genomic locations.
#'
#' Segment-level data are cleaned and annotated with patient identifiers,
#' estimated copy number, segment length, and a coarse CNV classification
#' (`"DEL"`, `"AMP"`, or `"NORMAL"`). Gene-level CNV status is then assigned by
#' intersecting genomic segments with a BED file of gene coordinates and
#' retaining, for each sample-gene pair, the segment with the largest absolute
#' segment mean.
#'
#' @param cnv_df A data.frame containing TCGA CNV segment-level data. It must
#'   include at least columns corresponding to sample identifiers, chromosome,
#'   segment start and end positions, and segment mean values.
#'
#' @param gene_bed_path A character string specifying the path to a BED file
#'   containing gene genomic coordinates. Default is
#'   `"data/raw/gencode.v34.annotation_gene.bed"`.
#'
#' @param gene_set_path A character string specifying the path to a CSV file
#'   containing the list of genes to retain for gene-level annotation. Default is
#'   `"data/raw/raids/genomic_data/RAIDs_WES_DecisionAlgo_OncoKBgenes_Pathways.csv"`.
#'
#' @param out_dir_hist A character string specifying the output directory for
#'   histogram images of segment means. If `NULL`, histograms are not saved.
#'   Default is `"outputs/figures"`.
#'
#' @param out_dir_profiles A character string specifying the output directory
#'   for sample-level CNV profile plots. If `NULL`, profile plots are not
#'   saved. Default is `"outputs/figures"`.
#'
#' @param n_profiles Integer specifying the maximum number of sample CNV
#'   profiles to plot. Default is `5`.
#'
#' @param hist_xlim Numeric vector of length 2 specifying the x-axis limits for
#'   the first segment mean histogram. Default is `c(-6, 2)`.
#'
#' @param hist_breaks1 Integer specifying the number of histogram breaks for
#'   the first segment mean histogram. Default is `1500`.
#'
#' @param hist_breaks2 Integer specifying the number of histogram breaks for
#'   the second segment mean histogram. Default is `500`.
#'
#' @return A named list with two elements:
#' \itemize{
#'   \item `segments`: a data.frame of processed CNV segments including cleaned
#'   identifiers, estimated copy number, segment length, variant classification,
#'   and chromosome factor.
#'   \item `gene_status`: a data.frame of gene-level CNV calls containing gene
#'   names, tumor sample barcodes, patient identifiers, and variant
#'   classification.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Cleans and renames CNV segment columns.
#'   \item Derives patient identifiers, estimated copy number, segment length,
#'   and a coarse CNV class based on `segment_mean` and segment size.
#'   \item Optionally saves histograms of segment mean distributions.
#'   \item Optionally generates CNV profile plots for up to `n_profiles`
#'   samples.
#'   \item Reads gene genomic coordinates from a BED file.
#'   \item Optionally filters genes using the gene list provided in
#'   `gene_set_path`.
#'   \item Intersects CNV segments with gene coordinates using
#'   `valr::bed_intersect()`.
#'   \item For each sample-gene pair, keeps the intersecting segment with the
#'   largest absolute segment mean.
#'   \item Returns both processed segment-level and gene-level CNV data.
#' }
#'
#' CNV classification is currently defined as:
#' \itemize{
#'   \item `"DEL"` if `segment_mean <= -2`,
#'   \item `"AMP"` if segment length is at most 10 Mb and estimated copy number
#'   is greater than 7,
#'   \item `"NORMAL"` otherwise.
#' }
#'
#' @export
process_tcga_cnv <- function(
    cnv_df,
    gene_bed_path = "data/raw/gencode.v34.annotation_gene.bed",
    gene_set_path = "data/raw/dna_gene_sets.csv",
    out_dir_hist = "outputs/figures",
    out_dir_profiles = "outputs/figures",       
    n_profiles = 5,
    hist_xlim = c(-6, 2),
    hist_breaks1 = 1500,
    hist_breaks2 = 500
){
  # Packages used inside
  stopifnot(is.data.frame(cnv_df))
  stopifnot(file.exists(gene_bed_path))
  stopifnot(file.exists(gene_set_path))
  
  # 1) preprocess segments 
  segments <- cnv_df %>%
    janitor::clean_names() %>%
    dplyr::rename(
      tumor_sample_barcode = gdc_aliquot_id,
      chrom = chromosome
    ) %>%
    dplyr::mutate(
      bcr_patient_barcode = substring(sample, 1, 12),
      copy_number = 2 * (2^segment_mean),
      length_bp = end - start + 1,
      variant_classification = dplyr::case_when(
        segment_mean <= -2 ~ "DEL",
        length_bp <= 1e7 & copy_number > (2 + 5) ~ "AMP",
        TRUE ~ "NORMAL"
      ),
      chrom_f = factor(
        chrom,
        levels = paste0("chr", c(1:22, "M", "X", "Y")),
        labels = c(1:22, "M", "X", "Y")
      )
    )
  
  # 2) optional histograms 
  if (!is.null(out_dir_hist)) {
    dir.create(out_dir_hist, recursive = TRUE, showWarnings = FALSE)
    
    grDevices::png(file.path(out_dir_hist, "tcga_hist_segment_mean.png"),
                   width = 800, height = 800, units = "px", res = 120)
    hist(segments$segment_mean, breaks = hist_breaks1, xlim = hist_xlim)
    grDevices::dev.off()
    
    grDevices::png(file.path(out_dir_hist, "tcga_hist_segment_mean2.png"),
                   width = 800, height = 800, units = "px", res = 120)
    hist(segments$segment_mean, breaks = hist_breaks2)
    grDevices::dev.off()
  }
  
  # 3) optional CNV profile plots 
  plot_profile <- function(sample_barcode){
    ggplot2::ggplot(
      dplyr::filter(segments, tumor_sample_barcode == sample_barcode),
      ggplot2::aes(
        x = start, xend = end,
        y = segment_mean, yend = segment_mean,
        color = variant_classification
      )
    ) +
      ggplot2::geom_segment(linewidth = 2) +
      ggplot2::facet_grid(. ~ chrom_f, scales = "free") +
      ggplot2::scale_x_continuous(breaks = NULL) +
      ggplot2::scale_color_manual(values = c(
        "DEL" = "limegreen",
        "NORMAL" = "darkgoldenrod1",
        "AMP" = "firebrick"
      )) +
      ggplot2::theme(
        panel.spacing = grid::unit(0.3, "lines"),
        legend.position = "bottom",
        panel.background = ggplot2::element_rect(fill = NA, color = "gray"),
        panel.grid.major = ggplot2::element_line(color = "gray70", linetype = 2),
        panel.grid.minor = ggplot2::element_line(color = "gray80", linetype = 2)
      )
  }
  
  plotted_samples <- character(0)
  if (!is.null(out_dir_profiles)) {
    dir.create(out_dir_profiles, recursive = TRUE, showWarnings = FALSE)
    
    samples <- unique(segments$tumor_sample_barcode)
    samples <- samples[seq_len(min(n_profiles, length(samples)))]
    plotted_samples <- samples
    
    for (s in samples){
      grDevices::png(file.path(out_dir_profiles, paste0(s, ".png")),
                     width = 1800, height = 800, units = "px", res = 120)
      print(plot_profile(s))
      grDevices::dev.off()
    }
  }
  
  # 4) gene-level status via bed_intersect 
  gene_location <- valr::read_bed(here::here(gene_bed_path))
  
  # optional gene filtering (expects gene names in column 'name' from gencode bed)
  gene_set <- read.csv2(here::here(gene_set_path))%>%
    dplyr::select(gene)%>%
    distinct()%>%
    pull()
  
  if (!is.null(gene_set)) {
    gene_location <- dplyr::filter(gene_location, .data$name %in% gene_set)
  }
  
  seg_for_intersect <- segments %>%
    dplyr::mutate(variant_classification = as.character(variant_classification))
  
  seg_genes <- valr::bed_intersect(seg_for_intersect, gene_location) %>%
    dplyr::mutate(abs_mean = abs(segment_mean.x)) %>%
    dplyr::group_by(sample.x, name.y) %>%
    dplyr::slice_max(order_by = abs_mean, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  gene_status <- seg_genes %>%
    dplyr::select(
      gene_name = name.y,
      tumor_sample_barcode = tumor_sample_barcode.x,
      variant_classification = variant_classification.x
    ) %>%
    dplyr::mutate(
      bcr_patient_barcode = substr(tumor_sample_barcode, 1, 12)
    )
  
  list(
    segments = segments,
    gene_status = gene_status
  )
}

#' Process TCGA SNV data
#'
#' This function preprocesses TCGA single-nucleotide variant (SNV) data by
#' cleaning column names, standardizing identifiers, and selecting relevant
#' variables for downstream analyses.
#'
#' It derives patient-level identifiers from sample barcodes, harmonizes gene
#' names, and creates a binary indicator for COSMIC annotation.
#'
#' @param snv_df A data.frame containing TCGA SNV data. It must include at least
#'   the following columns:
#'   \itemize{
#'     \item `tumor_sample_barcode`: sample identifier,
#'     \item `hugo_symbol`: gene symbol,
#'     \item `variant_classification`: mutation type,
#'     \item `cosmic`: COSMIC annotation (may contain missing values).
#'   }
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item `bcr_patient_barcode`: patient identifier (first 12 characters of
#'   the sample barcode),
#'   \item `tumor_sample_barcode`: sample identifier,
#'   \item `gene_name`: gene symbol,
#'   \item `variant_classification`: mutation type,
#'   \item `is_cosmic`: character indicator (`"true"` or `"false"`) specifying
#'   whether the variant is annotated in COSMIC.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Cleans column names using `janitor::clean_names()`.
#'   \item Derives a patient identifier (`bcr_patient_barcode`) from the
#'   sample barcode.
#'   \item Creates a binary COSMIC annotation (`is_cosmic`), where missing
#'   values are considered `"false"`.
#'   \item Renames the gene column (`hugo_symbol` to `gene_name`).
#'   \item Selects and returns relevant columns for downstream analyses.
#' }
#'
#' @export
process_tcga_snv <- function(snv_df){
  snv_df %>%
    janitor::clean_names() %>%
    dplyr::mutate(
      bcr_patient_barcode = substring(tumor_sample_barcode, 1, 12),
      is_cosmic = ifelse(is.na(cosmic), "false", "true")
    ) %>%
    dplyr::rename(
      gene_name = hugo_symbol
    ) %>%
    dplyr::select(
      bcr_patient_barcode,
      tumor_sample_barcode,
      gene_name,
      variant_classification,
      is_cosmic
    )
}

#' Process TCGA DNA alterations and build gene- and pathway-level matrices
#'
#' This function integrates TCGA SNV and CNV data, annotates genes using the
#' OncoKB cancer gene list, filters genes to a predefined pathway gene set,
#' applies a rule-based alteration algorithm, and derives patient-level
#' gene- and pathway-level DNA alteration features.
#'
#' The function returns intermediate annotation tables, a gene-by-patient
#' alteration matrix, and a long-format pathway-level DNA table suitable for
#' downstream analyses.
#'
#' @param tcga_snv A data.frame containing processed TCGA SNV data. It is
#'   expected to include at least the following columns:
#'   \itemize{
#'     \item `tumor_sample_barcode`,
#'     \item `bcr_patient_barcode`,
#'     \item `gene_name`,
#'     \item `variant_classification`,
#'     \item `is_cosmic`.
#'   }
#'
#' @param tcga_cnv A data.frame containing processed TCGA CNV data. It is
#'   expected to include at least the following columns:
#'   \itemize{
#'     \item `tumor_sample_barcode`,
#'     \item `bcr_patient_barcode`,
#'     \item `gene_name`,
#'     \item `variant_classification`.
#'   }
#'
#' @param oncokb_path A character string specifying the path to the OncoKB gene
#'   annotation file. The file is expected to be a tab-separated table
#'   containing at least the columns `Hugo Symbol`, `Is Oncogene`,
#'   `Is Tumor Suppressor Gene`, and `GRCh38 Isoform`. Default is
#'   `"data/raw/cancerGeneList.tsv"`.
#'
#' @param gene_set_path A character string specifying the path to a CSV file
#'   containing the gene-to-pathway mapping used to filter genes and derive
#'   pathway-level features. The file is expected to contain at least the
#'   columns `gene` and `pathway`. Default is
#'   `"data/raw/raids/genomic_data/RAIDs_WES_DecisionAlgo_OncoKBgenes_Pathways.csv"`.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item `snv_cnv`: merged SNV and CNV data,
#'   \item `oncokb_genes`: processed OncoKB gene annotation table,
#'   \item `oncokb_joined`: merged SNV/CNV table annotated with OncoKB genes,
#'   \item `oncokb_joined_filtered`: annotated table filtered to pathway genes,
#'   \item `mat`: a gene-by-patient character matrix of retained alteration
#'   classes,
#'   \item `dna_pathway_long`: a long-format table containing gene-level SNV,
#'   CNV, altered-gene, and altered-pathway indicators.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Merges processed SNV and CNV data by sample, patient, gene, and
#'   alteration class.
#'   \item Reads and preprocesses the OncoKB cancer gene annotation table.
#'   \item Restricts the analysis to patients with both SNV and CNV data.
#'   \item Reads the pathway gene set file and filters the merged DNA table to
#'   genes present in that reference set.
#'   \item Applies a rule-based alteration algorithm:
#'   \enumerate{
#'     \item oncogenes are altered by amplification or selected COSMIC-annotated
#'     activating SNVs,
#'     \item tumor suppressor genes are altered by deletion, truncating/splicing
#'     variants, or selected COSMIC-annotated SNVs,
#'     \item genes annotated as both, or neither oncogene nor tumor suppressor,
#'     are altered by either amplification, deletion, truncating/splicing
#'     variants, or selected COSMIC-annotated SNVs.
#'   }
#'   \item Builds a gene-by-patient matrix of retained alteration classes.
#'   \item Adds missing patients as empty columns and enforces a consistent
#'   patient order.
#'   \item Derives binary CNV and SNV indicators in long format.
#'   \item Joins gene-to-pathway annotations and computes gene-level and
#'   pathway-level alteration indicators.
#' }
#'
#' In the output `dna_pathway_long` table:
#' \itemize{
#'   \item `cnv` indicates whether the retained alteration includes `DEL` or
#'   `AMP`,
#'   \item `snv` indicates whether the retained alteration includes one of the
#'   selected SNV classes,
#'   \item `altered_gene` indicates whether either `snv` or `cnv` is present,
#'   \item `altered_pathway` indicates whether at least one gene in the pathway
#'   is altered in the patient.
#' }
#'
#' @export
process_tcga_dna <- function(
    tcga_snv,
    tcga_cnv,
    oncokb_path = "data/raw/cancerGeneList.tsv",
    gene_set_path = "data/raw/dna_gene_sets.csv"
){
  # joins SNVs + CNVs 
  tcga_snv_cnv <- dplyr::full_join(
    tcga_snv, tcga_cnv,
    by = c("tumor_sample_barcode",
           "bcr_patient_barcode",
           "gene_name",
           "variant_classification")
  )
  
  # OncoKB gene list processing 
  oncokb_genes <- readr::read_tsv(here::here(oncokb_path), 
                                  show_col_types = FALSE) %>%
    dplyr::select(`Hugo Symbol`, `Is Oncogene`,
                  `Is Tumor Suppressor Gene`, `GRCh38 Isoform`) %>%
    janitor::clean_names() %>%
    dplyr::rename(
      gene_id   = gr_ch38_isoform,
      gene_name = hugo_symbol,
      is_tsg    = is_tumor_suppressor_gene
    ) %>%
    dplyr::mutate(
      is_both = ifelse(is_oncogene == "Yes" & is_tsg == "Yes", "Yes", "No")
    )
  
  # Attach OncoKB + filter to DNA cohort
  tcga_dna_id <- intersect(unique(tcga_cnv$bcr_patient_barcode), 
                           unique(tcga_snv$bcr_patient_barcode))
  
  tcga_snv_cnv_oncokb <- tcga_snv_cnv %>%
    dplyr::full_join(oncokb_genes, by = "gene_name") %>%
    dplyr::filter(.data$bcr_patient_barcode %in% tcga_dna_id)
  
  # Pathway filter
    gene_set <- read.csv2(here::here(gene_set_path))%>%
      dplyr::rename(gene_name = gene)%>%
      dplyr::select(gene_name, pathway)%>%
      distinct()
    
    gene_filter <- gene_set %>%
      dplyr::select(gene_name)%>%
      distinct()%>%
      pull()
    tcga_snv_cnv_oncokb_filter <- tcga_snv_cnv_oncokb %>%
      dplyr::filter(.data$gene_name %in% gene_filter)
  
  # Apply alteration algorithm and build matrix 
  mat <- tcga_snv_cnv_oncokb_filter %>%
    dplyr::filter(
      # oncogene (not tsg): AMP or COSMIC in-frame del / missense
      ((is_oncogene == "Yes" & is_tsg == "No") & (
        (variant_classification == "AMP") |
          (variant_classification == "In_Frame_Del" & is_cosmic == "true") |
          (variant_classification == "Missense_Mutation" & is_cosmic == "true")
      )) |
        # (variant_classification == "TERT") |
        
        # tsg (not oncogene): DEL or COSMIC in-frame del / missense or trunc/splice
        ((is_tsg == "Yes" & is_oncogene == "No") & (
          (variant_classification == "DEL") | #(gene_name == "POLE") | (gene_name == "POLD1") |
            (variant_classification == "In_Frame_Del" & is_cosmic == "true") |
            (variant_classification == "Missense_Mutation" & is_cosmic == "true") |
            (variant_classification %in% c(
              "Nonsense_Mutation",
              "Frame_Shift_Del",
              "Frame_Shift_Ins",
              "Splice_Region",
              "Splice_Site"
            ))
        )) |
        
        # both or neither: AMP/DEL/COSMIC in-frame del/missense or trunc/splice
        ((is_both == "Yes" | (is_oncogene == "No" & is_tsg == "No")) & (
          (variant_classification == "AMP") |
            (variant_classification == "DEL") |
            (variant_classification == "In_Frame_Del" & is_cosmic == "true") |
            (variant_classification == "Missense_Mutation" & is_cosmic == "true") |
            (variant_classification %in% c(
              "Nonsense_Mutation",
              "Frame_Shift_Del",
              "Frame_Shift_Ins",
              "Splice_Region",
              "Splice_Site"
            ))
        ))
    ) %>%
    dplyr::select(gene_name, bcr_patient_barcode, variant_classification) %>%
    dplyr::distinct() %>%
    dplyr::group_by(gene_name, bcr_patient_barcode) %>%
    dplyr::summarise(
      variant_classification = paste(unique(variant_classification), 
                                     collapse = ";"),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = bcr_patient_barcode,
      values_from = variant_classification
    ) %>%
    tibble::column_to_rownames("gene_name")
  
  # add missing samples (columns) with NA
  missing_samples <- setdiff(tcga_dna_id, colnames(mat))
  if (length(missing_samples) > 0) {
    mat <- dplyr::mutate(
      mat,
      !!!stats::setNames(rep(list(NA_character_), 
                             length(missing_samples)), missing_samples)
    )
  }
  
  # enforce column order = tcga_DNA_id
  mat <- mat[, tcga_dna_id, drop = FALSE]
  
  # CNV : yes/no
  tcga_cnv_long <- mat %>%
    as.data.frame()%>%
    mutate(across(everything(), ~ifelse(.x %in% c("DEL","AMP"), 1, 0)))%>%
    rownames_to_column("gene_name")%>%
    pivot_longer(-gene_name, 
                 values_to = "cnv", 
                 names_to = "bcr_patient_barcode")%>%
    left_join(gene_set, by = "gene_name", 
              relationship = "many-to-many")
  
  # SNV : yes/no
  tcga_snv_long <- mat %>%
    as.data.frame()%>%
    rownames_to_column("gene_name")%>%
    mutate(across(-gene_name, ~ifelse(.x %in% c("Nonsense_Mutation",
                                                "Frame_Shift_Del",
                                                "Frame_Shift_Ins",
                                                "Splice_Region",
                                                "Splice_Site",
                                                "In_Frame_Del",
                                                "Missense_Mutation"), 1, 0)))%>%
    pivot_longer(-gene_name, 
                 values_to = "snv", 
                 names_to = "bcr_patient_barcode")%>%
    left_join(gene_set, by = "gene_name", 
              relationship = "many-to-many")
  
  # Pathway-level aggregation 
  tcga_dna_pathway_long <- 
    dplyr::inner_join(tcga_snv_long, tcga_cnv_long, 
                      by = c("gene_name","bcr_patient_barcode","pathway"))%>%
    dplyr::mutate(altered_gene = ifelse(snv == 1 | cnv == 1, 1, 0)) %>%
    dplyr::group_by(pathway, bcr_patient_barcode) %>%
    dplyr::mutate(altered_pathway = max(altered_gene, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dplyr::across(c(altered_pathway, altered_gene, cnv, snv), as.integer)
    )%>%
    dplyr::relocate(bcr_patient_barcode, gene_name, pathway)
  
  # Output
  list(
    snv_cnv = tcga_snv_cnv,
    oncokb_genes = oncokb_genes,
    oncokb_joined = tcga_snv_cnv_oncokb,
    oncokb_joined_filtered = tcga_snv_cnv_oncokb_filter,
    mat = as.matrix(mat),
    dna_pathway_long = tcga_dna_pathway_long
  )
}

#' Prepare TCGA RNA-seq expression data
#'
#' This function preprocesses TCGA RNA-seq expression data by selecting TPM
#' expression columns, harmonizing gene identifiers, removing unnecessary
#' annotation variables, filtering out missing gene names, and standardizing
#' sample column names to patient-level TCGA barcodes.
#'
#' In particular, version suffixes are removed from Ensembl gene identifiers,
#' except for pseudoautosomal region Y (`PAR_Y`) genes, which are preserved and
#' distinguished by appending `"_PAR_Y"` to the gene symbol.
#'
#' @param rna_df A data.frame containing TCGA RNA-seq data. It must include at
#'   least the following columns:
#'   \itemize{
#'     \item `gene_id`,
#'     \item `gene_name`,
#'     \item `gene_type`,
#'     \item one or more TPM expression columns starting with `"tpm_"`.
#'   }
#'
#' @return A tibble containing:
#' \itemize{
#'   \item `gene_name`: harmonized gene symbol,
#'   \item `gene_id`: processed gene identifier,
#'   \item one column per patient containing TPM expression values, with column
#'   names truncated to the first 12 characters of the TCGA barcode.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Selects gene annotation columns from `gene_id` to `gene_type` and
#'   all TPM expression columns.
#'   \item Stores the original TCGA gene identifier in a temporary column.
#'   \item Removes Ensembl version suffixes from `gene_id` values, except for
#'   identifiers containing `"PAR_Y"`.
#'   \item Appends `"_PAR_Y"` to `gene_name` for genes located in the
#'   pseudoautosomal Y region.
#'   \item Removes `gene_type` and the temporary original gene identifier.
#'   \item Filters out rows with empty gene names.
#'   \item Renames TPM sample columns by removing the `"tpm_unstranded_"`
#'   prefix and truncating barcodes to the first 12 characters, corresponding
#'   to TCGA patient identifiers.
#' }
#'
#' @export
prep_tcga_rna <- function(rna_df){
  rna_df %>%
    dplyr::select(gene_id:gene_type, 
                  dplyr::starts_with("tpm_"))%>% # select TPM data
    dplyr::mutate(
      gene_id_tcga = gene_id,
      gene_id = ifelse(!grepl("PAR_Y", gene_id_tcga), 
                       gsub("\\..*", "", gene_id_tcga), # remove . and all number after .
                       gene_id_tcga),
      gene_name = ifelse(grepl("PAR_Y", gene_id_tcga), 
                         paste0(gene_name, "_PAR_Y"),
                         gene_name))%>%
    dplyr::select(-c(gene_type, gene_id_tcga))%>%
    dplyr::filter(gene_name!="")%>%
    dplyr::rename_with(~ substring(gsub("tpm_unstranded_","", .x), 1, 12), 
                       -c(gene_name, gene_id)) 
}

#' Prepare TCGA RNA-seq data at the gene level in long format
#'
#' This function preprocesses TCGA RNA-seq data using `prep_tcga_rna()`,
#' aggregates expression values across duplicated gene symbols by taking the
#' mean, and reshapes the result into a long-format table with one row per
#' patient-gene pair.
#'
#' It is useful when gene-level TPM values are needed in a tidy format for
#' downstream analyses.
#'
#' @param rna_df A data.frame containing TCGA RNA-seq expression data. It must
#'   include the columns required by `prep_tcga_rna()`, namely:
#'   \itemize{
#'     \item `gene_id`,
#'     \item `gene_name`,
#'     \item `gene_type`,
#'     \item one or more TPM expression columns starting with `"tpm_"`.
#'   }
#'
#' @return A tibble in long format with the following columns:
#' \itemize{
#'   \item `gene_name`: gene symbol,
#'   \item `bcr_patient_barcode`: TCGA patient identifier,
#'   \item `tpm`: mean TPM expression value for the gene in the patient.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Calls `prep_tcga_rna()` to preprocess TCGA RNA-seq expression data.
#'   \item Removes the `gene_id` column.
#'   \item Aggregates duplicated gene symbols by computing the mean expression
#'   across all samples.
#'   \item Reshapes the data from wide to long format using `pivot_longer()`.
#' }
#'
#' Aggregation by mean is used to handle multiple rows mapping to the same
#' gene symbol.
#'
#' @export
prep_tcga_rna_gene <- function(rna_df){
  
  prep_tcga_rna(rna_df)%>%
  # We need to calculate the gene mean for each tumor sample (or patient)
    dplyr::select(-gene_id)%>%
    group_by(gene_name)%>%
    summarise(across(everything(),mean))%>%
    pivot_longer(-gene_name, 
                 names_to = "bcr_patient_barcode",
                 values_to = "tpm")
}
