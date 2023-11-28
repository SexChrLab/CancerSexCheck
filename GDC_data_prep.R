library(GenomicDataCommons)
library(tidyverse)
library(TCGAutils)


# Identify and Pull down meta-data:
case_ids <- cases() %>%
  GenomicDataCommons::filter(~ project.project_id == "TCGA-LIHC") %>%
  ids()
clindat <- gdc_clinical(case_ids)

# Pull out clinical data of interest into lists
# (for now; could streamline later when I know what I'm planning to do with it)
mf_list <- clindat[["demographic"]][["gender"]]
status_list <- clindat[["demographic"]][["vital_status"]]
surv_times  <- clindat[["demographic"]][["days_to_death"]]
follow_list <- clindat[["diagnoses"]][["days_to_last_follow_up"]]

meta_df <- data.frame(cbind(case_ids, mf_list, status_list,
                            surv_times, follow_list))

# write_tsv(
#   meta_df,
#   "TCGA_LIHC_META.tsv",
#   na = "NA",
#   append = FALSE,
#   col_names = TRUE,
#   quote = "none",
#   eol = "\n",
#   num_threads = readr_threads(),
#   progress = show_progress()
# )

### -----------------------------
# Create a manifest of files to download:

ge_manifest <- files() %>%
  GenomicDataCommons::filter(cases.project.project_id == "TCGA-LIHC") %>%
  GenomicDataCommons::filter(type == "gene_expression") %>%
  GenomicDataCommons::filter(analysis.workflow_type == "STAR - Counts")  %>%
  GenomicDataCommons::filter(access == "open") %>%
  manifest()

# Set up directory for download:
setwd("C:/Users/scmassey/Documents/CIREN/")
dir.create("TCGA-LIHC-StarCounts")
gdc_set_cache(directory = "TCGA-LIHC-StarCounts/")

# Download the counts data:
# can download just a subset to test: lapply(ge_manifest$id[1:100], gdcdata)
fnames <- lapply(ge_manifest, gdcdata)

# Create a list of files that contain the counts:
path <- "C:/Users/scmassey/Documents/CIREN/TCGA-LIHC-StarCounts/"

files <- dir(path = path, pattern = "*.counts", recursive = TRUE) %>%
  paste0(path, .)


### RECONFIGURING DATA FROM SEPARATE FILES NESTED IN INDIVIDUAL FOLDERS
### INTO THREE FILES, ONE FOR EACH TYPE OF COUNTS: TPM, FPKM, FPKM-UQ

# Initialize dataframes
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

# strip out the uuid from the filename string & find corresponding case uuid
file_uuid <- substr(files[1], 56, 91)
case_uuid <- UUIDtoUUID(file_uuid, to_type = "case_id")
combo_uuid <- paste(file_uuid, case_uuid[1, 2], sep = "_")

tpm_df <- subset(all_df, select = -c(unstranded, stranded_first,
                                     stranded_second, fpkm_unstranded,
                                     fpkm_uq_unstranded))
fpkm_df <- subset(all_df, select = -c(unstranded, stranded_first,
                                      stranded_second, tpm_unstranded,
                                      fpkm_uq_unstranded))
fpkm_uq_df <- subset(all_df, select = -c(unstranded, stranded_first,
                                         stranded_second, tpm_unstranded,
                                         fpkm_unstranded))

tpm_df <- tpm_df %>% rename_at("tpm_unstranded", ~combo_uuid)
fpkm_df <- fpkm_df %>% rename_at("fpkm_unstranded", ~combo_uuid)
fpkm_uq_df <- fpkm_uq_df %>% rename_at("fpkm_uq_unstrande", ~combo_uuid)


# LOOP THROUGH REMAINING FILES TO APPEND TPM COUNTS WTIH TCGA BARCODE COL LABELS

for (filename in files[2:length(files)]) {

  # strip out the uuid from the filename string & find corresponding case uuid
  file_uuid <- substr(filename, 56, 91)
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

  # below is a remaining bit from when tried to use joins instead of bind col -
  # if want to try this again, need to subset the mini_df differently to include
  # gene_name (or gene_id if want to use that instead)
  # tpm_df <- left_join(tpm_df, mini_df, by = "gene_name")
  # this was giving an error - "unexpected many-to-many"
}

# WRITE DATAFRAMEs TO TSV FILES:
write_tsv(
  meta_df,
  "TCGA_LIHC_META.tsv",
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

write_tsv(
  tpm_df,
  "TCGA_LIHC_TPM.tsv",
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
  "TCGA_LIHC_FPKM.tsv",
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
  "TCGA_LIHC_FPKM_UQ.tsv",
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)
