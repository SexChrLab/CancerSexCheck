library(GenomicDataCommons)

# Create a manifest of files to download:
ge_manifest <- files() %>%
  filter( cases.project.project_id == 'TCGA-LIHC') %>% 
  filter( type == 'gene_expression' ) %>%
  filter( analysis.workflow_type == 'STAR - Counts')  %>%
  filter( access=="open") %>%
  manifest()

### CAN ALSO SPECIFY WHICH FILTER FUNCTION TO USE:
# ge_manifest <- files() %>%
#   GenomicDataCommons::filter( cases.project.project_id == 'TCGA-LIHC') %>% 
#   GenomicDataCommons::filter( type == 'gene_expression' ) %>%
#   GenomicDataCommons::filter( analysis.workflow_type == 'STAR - Counts')  %>%
#   GenomicDataCommons::filter( access=="open") %>%
#   manifest()

# head(ge_manifest)

### NEED TO ALSO DO SOMETHING LIKE THIS TO PULL DOWN METADATA 
### TO SEE WHO IS LISTED AS "F"/"M" AND ALSO GET SURVIVAL DATA
### -- STILL WORKING ON THIS PART - WILL NEED TO MATCH UP DIFFERENT ID TYPES --
case_ids <-cases() %>%
  GenomicDataCommons::filter(~ project.project_id == 'TCGA-LIHC') %>%
  ids()
clindat = gdc_clinical(case_ids)
# head(clindat[["demographic"]]) 
### ----------------------------- 

# Set up directory for download:
setwd('C:/Users/scmassey/Documents/CIREN/')
dir.create('TCGA-LIHC-StarCounts')
gdc_set_cache(directory = 'TCGA-LIHC-StarCounts/')

# Download the counts data:
fnames <- lapply(ge_manifest$id[1:100], gdcdata)

# Create a list of files that contain the counts:
path = 'C:/Users/scmassey/Documents/CIREN/TCGA-LIHC-StarCounts/'

files <- dir(path=path, pattern="*.counts", recursive = TRUE) %>%
  paste0(path, .)

# head(files)


### RECONFIGURE DATA FROM SEPARATE FILES NESTED IN INDIVIDUAL FOLDERS
### INTO THREE FILES, ONE FOR EACH TYPE OF COUNTS: TPM, FPKM, FPKM-UQ

library(tidyverse)
library(TCGAutils)

# Initialize dataframes 
all_df <- read_tsv(files[1],skip=1,comment="N_", col_types=list(
  gene_id = "c", 
  gene_name = "c",
  gene_type = "c",
  unstranded = "i",
  stranded_first = "i",
  stranded_second = "i",
  tpm_unstranded = "d",
  fpkm_unstranded = "d",
  fpkm_uq_unstranded = "d")) 

uuid = substr(files[1],56,91) # uuid from folder name  (surprisingly, #uuid = substr(files[1],93,128) doesn't return a recognizable uuid
tcga_bar = UUIDtoBarcode(uuid, from_type = "file_id")

tpm_df <- subset(all_df, select = -c(unstranded,stranded_first,stranded_second,
                                     fpkm_unstranded,fpkm_uq_unstranded))
fpkm_df <- subset(all_df,select = -c(unstranded,stranded_first,stranded_second,
                                      tpm_unstranded,fpkm_uq_unstranded))
fpkm_uq_df<-subset(all_df,select= -c(unstranded,stranded_first,stranded_second,
                                     tpm_unstranded,fpkm_unstranded))

tpm_df <- tpm_df %>% rename_at('tpm_unstranded', ~tcga_bar[1,2])
fpkm_df <- fpkm_df %>% rename_at('fpkm_unstranded', ~tcga_bar[1,2])
fpkm_uq_df <- fpkm_uq_df %>% rename_at('fpkm_uq_unstranded', ~tcga_bar[1,2])


# LOOP THROUGH REMAINING FILES TO APPEND TPM COUNTS WTIH TCGA BARCODE COL LABELS

for (filename in files[2:length(files)]){  
  
  # strip out the uuid from the filename string and convert to barcode
  uuid = substr(filename,56,91)
  tcga_bar = UUIDtoBarcode(uuid, from_type = "file_id")
  
  # read in the file
  df <- read_tsv(filename,skip=1,comment="N_",col_types=list(
    gene_id = "c",
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
  mini_df <- subset(df, select = -c(gene_id,gene_name,gene_type,unstranded,
                                    stranded_first,stranded_second,
                                    fpkm_unstranded,fpkm_uq_unstranded))
  
  # add the current file's tpm counts to the aggregate set:
  tpm_df <- bind_cols(tpm_df, mini_df)
  
   # rename column from uuid to the tcga barcode:
  tpm_df <- tpm_df %>% rename_at('tpm_unstranded', ~tcga_bar[1,2]) 
  
  
  ### -- FPKM -- ###
  # subset to just the fpkm column
  mini_df <- subset(df, select = -c(gene_id,gene_name,gene_type,unstranded,
                                    stranded_first,stranded_second,
                                    tpm_unstranded,fpkm_uq_unstranded))
  
  # add the current file's tpm counts to the aggregate set:
  fpkm_df <- bind_cols(fpkm_df, mini_df)
  
  # rename column from uuid to the tcga barcode:
  fpkm_df <- fpkm_df %>% rename_at('fpkm_unstranded', ~tcga_bar[1,2])   
  
  
  ### -- FPKM -- ###
  # subset to just the fpkm column
  mini_df <- subset(df, select = -c(gene_id,gene_name,gene_type,unstranded,
                                    stranded_first,stranded_second,
                                    tpm_unstranded,fpkm_unstranded))
  
  # add the current file's tpm counts to the aggregate set:
  fpkm_uq_df <- bind_cols(fpkm_uq_df, mini_df)
  
  # rename column from uuid to the tcga barcode:
  fpkm_uq_df <- fpkm_uq_df %>% rename_at('fpkm_uq_unstranded', ~tcga_bar[1,2]) 
  
  # below is a remaining bit from when tried to use joins instead of bind col - 
  # if want to try this again, need to subset the mini_df differently to include 
  # gene_name (or gene_id if want to use that instead)
  # tpm_df <- left_join(tpm_df, mini_df, by = "gene_name") 
  # this was giving an error - "unexpected many-to-many"
  
}

# WRITE DATAFRAMEs TO TSV FILES: 
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
