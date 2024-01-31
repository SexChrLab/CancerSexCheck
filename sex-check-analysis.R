# WORKFLOW for CIREN TCGA SEX chromosomes vs Survival

# Import libraries to use
library(tidyverse)
library(ggplot2)
library(plotly)

# Specify which TCGA study to analyze
study = "LGG"

# Set directories
setwd("~/Desktop/2023-Fall-CIREN/")
filepath <- "/data/CEM/shared/public_data/TCGA_RNAseq_counts/"

counts_fname <- paste("TCGA", study, "TPM.tsv", sep = "-")
meta_fname <- paste("TCGA", study, "META.tsv", sep = "-")

counts_path <- paste(filepath, counts_fname, sep = "")
meta_path <- paste(filepath, meta_fname, sep = "")

# Read in data
counts <- read.delim(counts_path, row.names = 1)
metadf <- read.delim(meta_path, row.names = 1)

# Modify row id's in metadf to correspond to similar format in counts
meta_ids <- rownames(metadf)
meta_ids <- gsub("[-]", ".", meta_ids)
rownames(metadf) <- meta_ids

# Genes under consideration for inference:
#   chromosome Y
#     AMELY  ENSG00000099721
#     DDX3Y  ENSG00000067048
#     EIF1AY ENSG00000198692
#     KDM5D  ENSG00000012817
#     NLGN4Y ENSG00000165246
#     PRKY   ENSG00000099725
#     TMSB4Y ENSG00000154620
#     USP9Y  ENSG00000114374
#     UTY    ENSG00000183878
#     ZFY    ENSG00000067646
#     SRY    ENSG00000184895
#     TSPY   (may be challenging due to multicopy, there are 10+, left out)
#
#   chromosome X
#     XIST   ENSG00000229807
#     AR     ENSG00000169083

ychr_genes <- c("ENSG00000099721", "ENSG00000067048", "ENSG00000198692",
                "ENSG00000012817", "ENSG00000165246", "ENSG00000099725",
                "ENSG00000154620", "ENSG00000114374", "ENSG00000183878",
                "ENSG00000067646", "ENSG00000184895")
ychr_gnames <- c("AMELY", "DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "PRKY",
                 "TMSB4Y", "USP9Y", "UTY", "ZFY", "SRY")

xchr_genes  <- c("ENSG00000229807", "ENSG00000169083")
xchr_gnames <- c("XIST", "AR")

ychr_counts <- counts[ychr_genes, ]
xchr_counts <- counts[xchr_genes, ]

## may want a version of new_counts for x and y?:
new_counts1 <- data.frame(t(subset(xchr_counts,
                                   select = -c(gene_name, gene_type))))
colnames(new_counts1) <- xchr_counts[, 1]

new_counts2 <- data.frame(t(subset(ychr_counts,
                                   select = -c(gene_name, gene_type))))
colnames(new_counts2) <- ychr_counts[, 1]

new_counts <- data.frame(c(new_counts1, new_counts2))
row.names(new_counts) <- row.names(new_counts1)

################################################
###---CHECK FOR REPEATS AND LOOK AT VALUES---###
################################################

# Return list of row labels (ids):
id_list <- rownames(new_counts)

# Extract file and case uuid's from combo:
splt_id_list <- strsplit(id_list, "[_]")
c_uuid <- vector(mode = "character", length = length(id_list))
f_uuid <- vector(mode = "character", length = length(id_list))

for (j in seq_len(length(id_list))) {
  f_uuid[j] <- splt_id_list[[j]][1]
  c_uuid[j] <- splt_id_list[[j]][2]
}

# Organize ids by case uuid (for viewing)
id_sets <- split(id_list, c_uuid)

repeats_df <- data.frame()
rep_cases_list <- list()

c_uuid <- unique(c_uuid)

for (x in c_uuid) {
  if (length(id_sets[[x]]) > 1) {
    
    # extract labels needed for the next step
    rows_keep <- id_sets[[x]]
    
    # PUTTING THESE INTO A SEPARATE DF TO *LOOK* AT:
    repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[1], ])
    repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[2], ])
    
    if (length(id_sets[[x]]) > 2) {
      repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[3], ])
    }
    
    rep_cases_list <- c(rep_cases_list, x)
  }
}

# repeats_df dataframe to tsv
repeats_fname <- paste("TCGA", study, "repeat_samples_TPM_counts.tsv", sep = "_")
write_tsv(
  repeats_df %>% rownames_to_column(),
  repeats_fname,
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

###################################################
###---SEX CHECK OF THE SAMPLES---###
###################################################

sex_check <- data.frame(matrix(0, length(c_uuid), 5),
                        row.names = c_uuid)
colnames(sex_check) <- c("status_XIST", "status_Y", "annotated_sex",
                         "survival_status", "time_to_status")

dummy_df <- data.frame(matrix(0, length(id_list), 5))

counts_plus <- data.frame(c(new_counts, dummy_df))
row.names(counts_plus) <- row.names(new_counts)
colnames(counts_plus) <- c(colnames(new_counts), "status_XIST", "status_Y",
                        "annotated_sex", "survival", "followed")

### THE CODE BELOW HANDLES ONLY TWO SAMPLES PER CASE

for (i in c_uuid){
  sex_check[i, "annotated_sex"] <- metadf[i, "mf_list"]
  if (metadf[i, "status_list"] == "Alive") {
    sex_check[i, "survival_status"] <- 0
    sex_check[i, "time_to_status"] <- metadf[i, "follow_list"]
  } else {
    sex_check[i, "survival_status"] <- 1
    sex_check[i, "time_to_status"] <- metadf[i, "surv_times"]
  }
  # Check if there are replicates and if so, handle those first:
  if (length(id_sets[[i]]) > 1) {
    iid1 <- id_sets[[i]][1]
    iid2 <- id_sets[[i]][2]
    counts_plus[iid1, "annotated_sex"] <- metadf[i, "mf_list"]
    counts_plus[iid2, "annotated_sex"] <- metadf[i, "mf_list"]
    counts_plus[iid1, "survival"] <- metadf[i, "surv_times"]
    counts_plus[iid2, "survival"] <- metadf[i, "surv_times"]
    counts_plus[iid1, "followed"] <- metadf[i, "follow_list"]
    counts_plus[iid2, "followed"] <- metadf[i, "follow_list"]
    # Evaluate each replicate individually
    # -- Y chromosome:
    if (all(new_counts[iid1, "DDX3Y"] < 1.0,
            new_counts[iid1, "USP9Y"] < 1.0,
            new_counts[iid1, "UTY"] < 1.0,
            new_counts[iid1, "ZFY"] < 1.0)) {
      counts_plus[iid1, "status_Y"] <- "no"
    } else {
      counts_plus[iid1, "status_Y"] <- "yes"
    }
    if (all(new_counts[iid2, "DDX3Y"] < 1.0,
            new_counts[iid2, "USP9Y"] < 1.0,
            new_counts[iid2, "UTY"] < 1.0,
            new_counts[iid2, "ZFY"] < 1.0)) {
      counts_plus[iid2, "status_Y"] <- "no"
    } else {
      counts_plus[iid2, "status_Y"] <- "yes"
    }
    # -- XIST: 
    if (new_counts[iid1, "XIST"] > 1.0) {
      counts_plus[iid1, "status_XIST"] <- "yes"
    } else {
      counts_plus[iid1, "status_XIST"] <- "no"
    }
    if (new_counts[iid2, "XIST"] > 1.0) {
      counts_plus[iid2, "status_XIST"] <- "yes"
    } else {
      counts_plus[iid2, "status_XIST"] <- "no"
    }
    # Now evaluate the *pairs* of replicates (only written for 2 reps!)
    if ((counts_plus[iid1, "status_Y"] == counts_plus[iid2, "status_Y"]) &&
        (counts_plus[iid1, "status_XIST"] == counts_plus[iid2, "status_XIST"])) {
      # -- If the pairs match, take the first one as the *case* status
      sex_check[i, "status_Y"] <- counts_plus[iid1, "status_Y"]
      sex_check[i, "status_XIST"] <- counts_plus[iid1, "status_XIST"]
      # -- If the pairs DO NOT match, check the three possibilities:
      # -- 1. Do both look totally different?
      #       (i.e. "yes, no" vs "no, yes" - or - "yes yes" vs "no no")
    } else if ((counts_plus[iid1, "status_Y"] != counts_plus[iid2, "status_Y"]) &&
               (counts_plus[iid1, "status_XIST"] != counts_plus[iid2, "status_XIST"])) {
      # this is an unexpected edge case, so for now, just "flag" it as unusual
      sex_check[i, "status_Y"] <- "FLAG"
      sex_check[i, "status_XIST"] <- "FLAG"
      # -- 2. The first pair has "yes yes" or "no no":
    } else if (counts_plus[iid1, "status_XIST"] == counts_plus[iid1, "status_Y"]) {
      # Take "yes XIST, yes Y" or "no XIST, no Y" as the overall case status
      sex_check[i, "status_Y"] <- counts_plus[iid1, "status_Y"]
      sex_check[i, "status_XIST"] <- counts_plus[iid1, "status_XIST"]
      # -- 3. The second pair has "yes yes" or "no no":
    } else if (counts_plus[iid2, "status_XIST"] == counts_plus[iid2, "status_Y"]) {
      # Take "yes XIST, yes Y" or "no XIST, no Y" as the overall case status
      sex_check[i, "status_Y"] <- counts_plus[iid2, "status_Y"]
      sex_check[i, "status_XIST"] <- counts_plus[iid2, "status_XIST"]
    } else {
      # That should have handled all the cases, but just in case
      # something goes wrong with the evaluation of the conditionals
      print("we missed something")
    }
    
    # Now look at cases with only a single sample:
  } else {
    iid <- id_sets[[i]]
    counts_plus[iid, "annotated_sex"] <- metadf[i, "mf_list"]
    counts_plus[iid, "survival"] <- metadf[i, "surv_times"]
    counts_plus[iid, "followed"] <- metadf[i, "follow_list"]
    if (all(new_counts[iid, "DDX3Y"] < 1.0,
            new_counts[iid, "USP9Y"] < 1.0,
            new_counts[iid, "UTY"] < 1.0,
            new_counts[iid, "ZFY"] < 1.0)) {
      counts_plus[iid, "status_Y"] <- "no"
      sex_check[i, "status_Y"] <- "no"
    } else {
      counts_plus[iid, "status_Y"] <- "yes"
      sex_check[i, "status_Y"] <- "yes"
    }
    if (new_counts[iid, "XIST"] > 1.0) {
      counts_plus[iid, "status_XIST"] <- "yes"
      sex_check[i, "status_XIST"] <- "yes"
    } else {
      counts_plus[iid, "status_XIST"] <- "no"
      sex_check[i, "status_XIST"] <- "no"
    }
  }
}

## THE CODE ABOVE HANDLES ONLY UP TO TWO SAMPLES PER CASE, BUT
## CASE d6486001.240a.455a.980c.e06c25c61fa5 HAS THREE SAMPLES.
## THEY ARE: id_sets[["d6486001.240a.455a.980c.e06c25c61fa5"]]
## "X78a0f8f9.e010.4a10.978c.94c8bb9157cd_d6486001.240a.455a.980c.e06c25c61fa5"
## "X7b90b9fd.0015.47b9.9148.f040c1cfcb5a_d6486001.240a.455a.980c.e06c25c61fa5"
## "c7a64911.e1b0.4615.9521.98d4cd4a9882_d6486001.240a.455a.980c.e06c25c61fa5"

### USING A MANUAL APPROACH TO COMPLETE THE CORRESPONDING COUNTS_PLUS DATAFRAME
third_id <- "c7a64911.e1b0.4615.9521.98d4cd4a9882_d6486001.240a.455a.980c.e06c25c61fa5"
tcase_id <- "d6486001.240a.455a.980c.e06c25c61fa5"
other_id <- "X7b90b9fd.0015.47b9.9148.f040c1cfcb5a_d6486001.240a.455a.980c.e06c25c61fa5"
counts_plus[third_id, "status_XIST"] <- counts_plus[other_id, "status_XIST"]
counts_plus[third_id, "status_Y"] <- counts_plus[other_id, "status_Y"]
counts_plus[third_id, "annotated_sex"] <- metadf[tcase_id, "mf_list"]
counts_plus[third_id, "survival"] <- metadf[tcase_id, "surv_times"]
counts_plus[third_id, "followed"] <- metadf[tcase_id, "follow_list"]

###################################################
###---WRITE DATAFRAMES TO TSV FILES---###
###################################################

# sex_check dataframe to tsv
sex_check_fname <- paste("TCGA", study, "case_XIST-Y_outcomes.tsv", sep = "_")
write_tsv(
  sex_check %>% rownames_to_column(),
  sex_check_fname,
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

# counts_plus dataframe to tsv
counts_plus_fname <- paste("TCGA", study, "sample_XIST-Y_outcomes.tsv", sep = "_")

write_tsv(
  counts_plus %>% rownames_to_column(),
  counts_plus_fname,
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

###################################################
###---SURVIVAL ANALYSIS---###
###################################################
library(dplyr)
library(survival)
library(survminer)

sex_check_m <- sex_check %>% filter(annotated_sex == "male")
df_points2 <- sex_check_m %>%
  transmute(XIST_Y = paste(sex_check_m$status_XIST, "XIST_",
                           sex_check_m$status_Y, "Y", sep = ""))
sex_check_m <- cbind(sex_check_m, df_points2)

sex_check_f <- sex_check %>% filter(annotated_sex == "female")
df_points3 <- sex_check_f %>%
  transmute(XIST_Y = paste(sex_check_f$status_XIST, "XIST_",
                           sex_check_f$status_Y, "Y", sep = ""))
sex_check_f <- cbind(sex_check_f, df_points3)

km_m <- survfit(Surv(time_to_status, survival_status) ~ XIST_Y, data = sex_check_m)
km_f <- survfit(Surv(time_to_status, survival_status) ~ XIST_Y, data = sex_check_f)


km_m_plot_fname <- paste("KM", study, "Male.png", sep = "_")
png(km_m_plot_fname)

km_m %>%
  ggsurvplot(
    data = sex_check_m,
    fun = "pct",
    # linetype = "strata", # Change line type by groups
    # pval = TRUE, # Not sure if want
    # conf.int = TRUE, # Not sure if want
    risk.table = TRUE,
    fontsize = 3, # used in risk table
    surv.median.line = "hv", # median horizontal and vertical ref lines
    ggtheme = theme_light(),
    palette = c("goldenrod", "sienna", "tomato", "cadetblue"),
    title = "LIHC - Male - Kaplan-Meier Survival Function Estimate",
    legend.title = "",
    legend.labs = levels(sex_check_m$XIST_Y)
  )
dev.off()

km_f_plot_fname <- paste("KM", study, "Female.png", sep = "_")
png(km_f_plot_fname)

km_f %>%
  ggsurvplot(
    data = sex_check_f,
    fun = "pct",
    # linetype = "strata", # Change line type by groups
    # pval = TRUE, # Not sure if want
    # conf.int = TRUE, # Not sure if want
    risk.table = TRUE,
    fontsize = 3, # used in risk table
    surv.median.line = "hv", # median horizontal and vertical ref lines
    ggtheme = theme_light(),
    palette = c("goldenrod", "sienna", "tomato","cadetblue"),
    title = "LIHC - Female - Kaplan-Meier Survival Function Estimate",
    legend.title = "",
    legend.labs = levels(sex_check_f$XIST_Y)
  )
dev.off()

###################################################
###---COUNT GROUP MEMBERSHIP FOR TABLE---###
###################################################

# among those annotated "male":
n_male_nxny <- nrow(sex_check_m[sex_check_m$XIST_Y == "noXIST_noY", ])
n_male_yxyy <- nrow(sex_check_m[sex_check_m$XIST_Y == "yesXIST_yesY", ])
n_male_nxyy <- nrow(sex_check_m[sex_check_m$XIST_Y == "noXIST_yesY", ])
n_male_yxny <- nrow(sex_check_m[sex_check_m$XIST_Y == "yesXIST_noY", ])

# among those annotated "female":
n_female_nxny <- nrow(sex_check_f[sex_check_f$XIST_Y == "noXIST_noY", ])
n_female_yxyy <- nrow(sex_check_f[sex_check_f$XIST_Y == "yesXIST_yesY", ])
n_female_nxyy <- nrow(sex_check_f[sex_check_f$XIST_Y == "noXIST_yesY", ])
n_female_yxny <- nrow(sex_check_f[sex_check_f$XIST_Y == "yesXIST_noY", ])


#####################################################
###-VIOLIN PLOTS OF TPM COUNTS OF INTEREST BY M/F-###
###-TO SHOW IF THRESHOLD WAS CHOSEN APPROPRIATELY-###
#####################################################
library(patchwork)

violin1_fname <- paste("Log_Violin", study, "TPM_XIST.png", sep = "_")
png(violin1_fname)
pXIST <- ggplot(data = counts_plus, aes(factor(annotated_sex), XIST))
pXIST <- pXIST + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)
dev.off()

violin2_fname <- paste("Log_Violin", study, "TPM_DDX3Y.png", sep = "_")
png(violin2_fname)
pDDX3Y <- ggplot(data = counts_plus, aes(factor(annotated_sex), DDX3Y))
pDDX3Y <- pDDX3Y + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)
dev.off()

violin3_fname <- paste("Log_Violin", study, "TPM_USP9Y.png", sep = "_")
png(violin3_fname)
pUSP9Y <- ggplot(data = counts_plus, aes(factor(annotated_sex), USP9Y))
pUSP9Y <- pUSP9Y + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)
dev.off()

violin4_fname <- paste("Log_Violin", study, "TPM_UTY.png", sep = "_")
png(violin4_fname)
pUTY <- ggplot(data = counts_plus, aes(factor(annotated_sex), UTY))
pUTY <- pUTY + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)
dev.off()

violin5_fname <- paste("Log_Violin", study, "TPM_ZFY.png", sep = "_")
png(violin5_fname)
pZFY <- ggplot(data = counts_plus, aes(factor(annotated_sex), ZFY))
pZFY <- pZFY + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)
dev.off()

violin6_fname <- paste("Log_Violin", study, "TPM_All.png", sep = "_")
png(violin6_fname)
p_all <-pXIST | (pDDX3Y | pUSP9Y) / (pUTY | pZFY)
p_all
dev.off()

# Violins in plotly for interactive viewing
figDDX3Y <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~DDX3Y,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

figUSP9Y <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~USP9Y,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

figUTY <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~UTY,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

figZFY <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~ZFY,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

figXIST <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~XIST,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

fig <- subplot(figDDX3Y, figUSP9Y, figUTY, figZFY, figXIST,
               nrows = 3, shareY = TRUE) %>%
  layout(title = "Distribution of TPM counts for genes of interest by sex",
         plot_bgcolor = "#e5ecf6",
         showlegend = FALSE
  )

fig

#####################################################
###---LINE PLOTS OF TPM COUNTS FOR THE ODD COMBOS-###
###---OF LOW XIST & LOW Y OR HIGH XIST & HIGH Y---###
###---TO ASSESS FOR ANY PATTERNS WITHIN SAMPLES---###
#####################################################

# Two Groups of Interest:

# 1. no XIST and no Y
nxny_counts <- counts_plus %>% filter(status_XIST == "no" & status_Y == "no")

# 2. yes XIST and yes Y
yxyy_counts <- counts_plus %>% filter(status_XIST == "yes" & status_Y == "yes")

# Prepping list of gene names for plot functions below:
gene_names <- c("XIST", "DDX3Y", "USP9Y", "UTY", "ZFY")

# Reshape to work well with plot functions:
nxny_plots <- reshape(nxny_counts,
                      varying = gene_names,
                      drop = c("AR", "AMELY", "EIF1AY", "KDM5D",
                               "NLGN4Y", "PRKKY", "TMSB4Y", "SRY",
                               "status_XIST", "status_Y",
                               "survival", "followed"),
                      v.names = "TPM_counts",
                      timevar = "gene_names",
                      times = gene_names,
                      ids = row.names(nxny_counts),
                      direction = "long")

yxyy_plots <- reshape(yxyy_counts,
                      varying = gene_names,
                      drop = c("AR", "AMELY", "EIF1AY", "KDM5D",
                               "NLGN4Y", "PRKKY", "TMSB4Y", "SRY",
                               "status_XIST", "status_Y",
                               "survival", "followed"),
                      v.names = "TPM_counts",
                      timevar = "gene_names",
                      times = gene_names,
                      ids = row.names(yxyy_counts),
                      direction = "long")

# Plot the lines

# Linear y-axis
lineplot1_fname <- paste("TPM", study, "LowXIST-LowY_Color-by-Sex.png", sep = "_")
png(lineplot1_fname)
ggplot(data = nxny_plots, aes(x = gene_names, y = TPM_counts,
                              group = id, color = annotated_sex)) +  # color = id
  scale_color_discrete(guide = "none") + geom_point() + geom_line() +
  ggtitle("LIHC samples with No/Low XIST No/Low Y chr, inclusive of repeats") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

lineplot2_fname <- paste("TPM", study, "HighXIST-HighY_Color-by-Sex.png", sep = "_")
png(lineplot2_fname)
ggplot(data = yxyy_plots, aes(x = gene_names, y = TPM_counts,
                              group = id, color = annotated_sex)) +  # color = id
  scale_color_discrete(guide = "none") + geom_point() + geom_line() +
  ggtitle("LIHC samples with XIST & Y markers TPM>1.0, inclusive of repeats") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Log y-axis
loglineplot1_fname <- paste("TPM", study, "LowXIST-LowY_Color-by-Sex_Log-axis.png", sep = "_")
png(loglineplot1_fname)
ggplot(data = nxny_plots, aes(x = gene_names, y = TPM_counts,
                              group = id, color = annotated_sex)) +
  scale_color_discrete(guide = "none") +
  geom_point() + geom_line() + scale_y_continuous(trans = "log10")
dev.off()


loglineplot2_fname <- paste("TPM", study, "HighXIST-HighY_Color-by-Sex_Log-axis.png", sep = "_")
png(loglineplot2_fname)
ggplot(data = yxyy_plots, aes(x = gene_names, y = TPM_counts,
                              group = id, color = annotated_sex)) +
  scale_color_discrete(guide = "none") +
  geom_point() + geom_line() + scale_y_continuous(trans = "log10")
dev.off()
