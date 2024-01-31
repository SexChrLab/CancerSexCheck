library(GenomicDataCommons)
library(tidyverse)
library(TCGAutils)

### TCGA Study Abbreviations
### - LAML - Acute Myeloid Leukemia
### - ACC  - Adrenocortical Carcinoma
### - BLCA - Bladder Urothelial Carcinoma
### - LGG  - Brain Lower Grade Glioma
### - BRCA - Breast Invasive Carcinoma
### - CESC - Cervical Squamous cell Carcinoma and Endocervical Adenocarcinoma
### - CHOL - Cholangiocarcinoma
### - COAD - Colon Adenocarcinoma
### - ESCA - Esophageal Adenocarcinoma
### - GBM  - Glioblastoma
### - HNSC - Head and Neck Squamous cell Carcinoma
### - KICH - Kidney Chromophobe
### - KIRC - Kidney Renal Clear cell Carcinoma
### - KIRP - Kidney Renal Papillary cell Carcinoma
### - LIHC - Liver Hepatocellular Carcinoma
### - LUAD - Lung Adenocarcinoma
### - LUSC - Lung Squamous cell Carcinoma
### - DLBC - Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
### - MESO - Mesothelioma
### - OV   - Ovarian serous cystadenocarcinoma
### - PAAD - Pancreatic Adenocarcinoma
### - PCPG - Pheochromocytoma and Paraganglioma
### - PRAD - Prostate Adenocarcinoma
### - READ - Rectum Adenocarcinoma
### - SARC - Sarcoma
### - SKCM - Skin Cutaneous Melanoma
### - STAD - Stomach Adenocarcinoma
### - TGCT - Testicular Germ Cell Tumor
### - THCA - Thyroid gland
### - THYM - Thymoma
### - UCS  - Uterine Carcinosarcoma
### - UCEC - Uterine Corpus Endometrial Carcinoma
### - UVM  - Uveal Melanoma
### From https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
###  - edited to match with actual sets available on the GDC portal
###  - specifically, removed LCML, CNTL, FPPP, and MISC; added THCA

### List the TCGA types to run
#tcga_studies_to_prep <- c("LGG", "GBM", "LIHC")

working_path <- "~/Desktop/2023-Fall-CIREN/"
setwd(working_path)
data_out_path <- "/data/CEM/shared/public_data/TCGA_RNAseq_counts/"



########################################
### LOOP OVER SPECIFIED TCGA STUDIES ###
########################################

#for (t in tcga_studies_to_prep) {
  t = "GBM"
  study <- paste("TCGA", t, sep = "-")

  ################
  ### METADATA ###
  ################

  # Identify and Extract Relevant Metadata:
  case_ids <- cases() %>%
    GenomicDataCommons::filter(~ project.project_id == study) %>%
    ids()
  clindat <- gdc_clinical(case_ids)

  # Pull out clinical data of interest into lists
  mf_list <- clindat[["demographic"]][["gender"]]
  status_list <- clindat[["demographic"]][["vital_status"]]
  surv_times  <- clindat[["demographic"]][["days_to_death"]]
  follow_list <- clindat[["diagnoses"]][["days_to_last_follow_up"]]

  # Combine these into a slimmed down dataframe of metadata
  meta_df <- data.frame(cbind(case_ids, mf_list, status_list,
                              surv_times, follow_list))

  # Save out metadata as a tsv file
  meta_fname <- paste(study, 'META.tsv', sep = "-")
  meta_path <- paste(data_out_path, meta_fname, sep = "")
  write_tsv(
    meta_df,
    meta_path,
    na = "NA",
    append = FALSE,
    col_names = TRUE,
    quote = "none",
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress()
  )

  ########################
  ### STAR COUNTS DATA ###
  ########################

  ### STEP 1: DOWNLOAD THE STAR COUNTS DATA

  # Create a manifest of files to download:
  ge_manifest <- files() %>%
    GenomicDataCommons::filter(cases.project.project_id == study) %>%
    GenomicDataCommons::filter(type == "gene_expression") %>%
    GenomicDataCommons::filter(analysis.workflow_type == "STAR - Counts")  %>%
    GenomicDataCommons::filter(access == "open") %>%
    manifest()

  # Set up directory for download
  # - this will be a subdirectory in the working directory
  download_pathname <- paste(study, "StarCounts", sep = "-")
  dir.create(download_pathname)
  download_path <- paste(download_pathname, "/", sep = "")
  gdc_set_cache(directory = download_pathname)

  # Download the counts data:
  fnames <- lapply(ge_manifest$id, gdcdata)
  # can download a test subset: fnames <- lapply(ge_manifest$id[1:15], gdcdata)

  ### STEP 2: REARRANGE THE STAR COUNTS DATA
  # Create a list of files that contain the counts:
  path <- paste(working_path, download_path, sep = "")

  files <- dir(path = path, pattern = "*.counts", recursive = TRUE) %>%
    paste0(path, .)

  # Initialize dataframes using the first star counts data file in the list

  # Dataframe to contain all of the combined counts data
  all_df <- read_tsv(files[1], skip = 1, comment = "N_",
                     col_types = list(gene_id = "c",
                                      gene_name = "c",
                                      gene_type = "c",
                                      unstranded = "i",
                                      stranded_first = "i",
                                      stranded_second = "i",
                                      tpm_unstranded = "d",
                                      fpkm_unstranded = "d",
                                      fpkm_uq_unstranded = "d"))

  # Dataframes to contain each type of count separated out
  # - TPM, FPKM, and FPKM_UQ
  tpm_df <- subset(all_df, select = -c(unstranded, stranded_first,
                                       stranded_second, fpkm_unstranded,
                                       fpkm_uq_unstranded))
  fpkm_df <- subset(all_df, select = -c(unstranded, stranded_first,
                                        stranded_second, tpm_unstranded,
                                        fpkm_uq_unstranded))
  fpkm_uq_df <- subset(all_df, select = -c(unstranded, stranded_first,
                                           stranded_second, tpm_unstranded,
                                           fpkm_unstranded))

  # Extract the uuid from the filename string & find corresponding case uuid;
  # then concatenate these to create a combination uuid 
  # - having both will help identify any repeated samples from same case
  # - as well as match up with metadata
  pre_fuuid <- gsub(paste(working_path, download_path, sep = ""), "", files[1])
  file_uuid <- strsplit(pre_fuuid, "[/]")[[1]][1]
  case_uuid <- UUIDtoUUID(file_uuid, to_type = "case_id")
  combo_uuid <- paste(file_uuid, case_uuid[1, 2], sep = "_")

  # Add these combined ids as row names for the dataframes
  tpm_df <- tpm_df %>% rename_at("tpm_unstranded", ~combo_uuid)
  fpkm_df <- fpkm_df %>% rename_at("fpkm_unstranded", ~combo_uuid)
  fpkm_uq_df <- fpkm_uq_df %>% rename_at("fpkm_uq_unstranded", ~combo_uuid)

  # Loop through remaining files in the list to add them to the dataframes
  for (filename in files[2:length(files)]) {

    # extract the uuid from the filename string & find corresponding case uuid
    pre_fuuid <- gsub(paste(working_path, download_path, sep = ""), "", filename)
    file_uuid <- strsplit(pre_fuuid, "[/]")[[1]][1]
    case_uuid <- UUIDtoUUID(file_uuid, to_type = "case_id")
    combo_uuid <- paste(file_uuid, case_uuid[1, 2], sep = "_")

    # read in the file
    df <- read_tsv(filename, skip = 1, comment = "N_",
                   col_types = list(gene_id = "c",
                                    gene_name = "c",
                                    gene_type = "c",
                                    unstranded = "i",
                                    stranded_first = "i",
                                    stranded_second = "i",
                                    tpm_unstranded = "d",
                                    fpkm_unstranded = "d",
                                    fpkm_uq_unstranded = "d"))


    ### -- TPM -- ###
    # subset to just the tpm column
    mini_df <- subset(df, select = -c(gene_id, gene_name, gene_type, unstranded,
                                      stranded_first, stranded_second,
                                      fpkm_unstranded, fpkm_uq_unstranded))

    # add the current file's tpm counts to the aggregate set:
    tpm_df <- bind_cols(tpm_df, mini_df)

    # rename column from uuid to the combination id:
    tpm_df <- tpm_df %>% rename_at("tpm_unstranded", ~combo_uuid)


    ### -- FPKM -- ###
    # subset to just the fpkm column
    mini_df <- subset(df, select = -c(gene_id, gene_name, gene_type, unstranded,
                                      stranded_first, stranded_second,
                                      tpm_unstranded, fpkm_uq_unstranded))

    # add the current file's tpm counts to the aggregate set:
    fpkm_df <- bind_cols(fpkm_df, mini_df)

    # rename column from uuid to the combination id:
    fpkm_df <- fpkm_df %>% rename_at("fpkm_unstranded", ~combo_uuid)


    ### -- FPKM -- ###
    # subset to just the fpkm column
    mini_df <- subset(df, select = -c(gene_id, gene_name, gene_type, unstranded,
                                      stranded_first, stranded_second,
                                      tpm_unstranded, fpkm_unstranded))

    # add the current file's tpm counts to the aggregate set:
    fpkm_uq_df <- bind_cols(fpkm_uq_df, mini_df)

    # rename column from uuid to the combination id:
    fpkm_uq_df <- fpkm_uq_df %>% rename_at("fpkm_uq_unstranded", ~combo_uuid)

  }

  # STEP 3: WRITE DATAFRAMES TO TSV FILES:

  write_tsv(
    tpm_df,
    paste(data_out_path, study, "-TPM.tsv", sep = ""),
    na = "NA",
    append = FALSE,
    col_names = TRUE,
    quote = "none",
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress()
  )

  write_tsv(
    fpkm_df,
    paste(data_out_path, study, "-FPKM.tsv", sep = ""),
    na = "NA",
    append = FALSE,
    col_names = TRUE,
    quote = "none",
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress()
  )

  write_tsv(
    fpkm_uq_df,
    paste(data_out_path, study, "-FPKM_UQ.tsv", sep = ""),
    na = "NA",
    append = FALSE,
    col_names = TRUE,
    quote = "none",
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress()
  )

#}
