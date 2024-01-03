# WORKFLOW for CIREN TCGA SEX chromosomes vs Survival

# Import libraries to use
library(tidyverse)
library(ggplot2)
library(plotly)

# Read in data files of interest/to use
setwd("~/Desktop/2023-Fall-CIREN/")
counts <- read.delim("/data/CEM/shared/public_data/TCGA_RNAseq_counts/TCGA_LIHC_TPM.tsv", row.names = 1)
metadf <- read.delim("/data/CEM/shared/public_data/TCGA_RNAseq_counts/TCGA_LIHC_META.tsv", row.names = 1)

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

# make one violin plot each for reported males and
# one for reported females (with/without log transformation)

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
    # whatdo
    rows_keep <- id_sets[[x]]

    # PUTTING THESE INTO A SEPARATE DF TO *LOOK* AT:
    repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[1], ])
    repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[2], ])

    rep_cases_list <- c(rep_cases_list, substr(rows_keep[1], 39, 74))
  }
}

###################################################
###--END - CHECK FOR REPEATS AND LOOK AT VALUES-###
###################################################

### FINALLY, SEX CHECK OF THE SAMPLES!

sex_check <- data.frame(matrix(0, length(c_uuid), 5),
                        row.names = c_uuid)
colnames(sex_check) <- c("status_XIST", "status_Y", "assigned",
                         "survival", "followed")

dummy_df <- data.frame(matrix(0, length(id_list), 5))

counts_plus <- data.frame(c(new_counts, dummy_df))
row.names(counts_plus) <- row.names(new_counts)
colnames(counts_plus) <- c(colnames(new_counts), "status_XIST", "status_Y",
                           "assigned", "survival", "followed")

for (i in c_uuid){
  sex_check[i, "assigned"] <- metadf[i, "mf_list"]
  sex_check[i, "survival"] <- metadf[i, "surv_times"]
  sex_check[i, "followed"] <- metadf[i, "follow_list"]
  if (length(id_sets[[i]]) > 1) {
    iid1 <- id_sets[[i]][1]
    iid2 <- id_sets[[i]][2]
    counts_plus[iid1, "assigned"] <- metadf[i, "mf_list"]
    counts_plus[iid2, "assigned"] <- metadf[i, "mf_list"]
    counts_plus[iid1, "survival"] <- metadf[i, "surv_times"]
    counts_plus[iid2, "survival"] <- metadf[i, "surv_times"]
    counts_plus[iid1, "followed"] <- metadf[i, "follow_list"]
    counts_plus[iid2, "followed"] <- metadf[i, "follow_list"]
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
    if (counts_plus[iid1, "status_Y"] == "yes" &&
        counts_plus[iid2, "status_Y"] == "yes") {
      sex_check[i, "status_Y"] <- "yes"
    } else {
      sex_check[i, "status_Y"] <- "no"
    }
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
    if (counts_plus[iid1, "status_XIST"] == "yes" &&
        counts_plus[iid2, "status_XIST"] == "yes"){
      sex_check[i ,"status_XIST"] == "yes"
    } else {
      sex_check[i, "status_XIST"] == "no"
    }
  } else {
    iid <- id_sets[[i]]
    counts_plus[iid, "assigned"] <- metadf[i, "mf_list"]
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

# NEED TO ADD ROW IDS IN A COLUMN TO USE
write_tsv(
  sex_check %>% rownames_to_column(),
  "TCGA_LIHC_inferences.tsv",
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

###################################################
###################################################
# PLOT CODE - NOT YET UPDATED SINCE CHANGE TO XIST-yes/no vs Y-chr-yes/no

# # First,
# # want to look at counts for each gene among the no_inference group
# # with a line to connect subjects across the different genes
# 
# not_inferred <- counts_plus %>% filter(inferred == "no_inference")

# Two Groups of Interest:

# 1. no XIST and no Y
no_Xist_no_Y <- sex_check %>% filter(status_XIST == "no" & status_Y == "no")
no_Xist_no_Y_counts <- counts_plus %>% filter(status_XIST == "no" & status_Y == "no")

# 2. yes XIST and yes Y
yes_Xist_yes_Y <- sex_check %>% filter(status_XIST == "yes" & status_Y == "yes")
yes_Xist_yes_Y_counts <- counts_plus %>% filter(status_XIST == "yes" & status_Y == "yes")

# Prepping list of gene names for plot functions below:
# gene_names <- c(xchr_gnames, ychr_gnames)
# # Drop SRY
# gene_names <- gene_names[-13]
# # Drop genes not used for inference
gene_names_short <- c("XIST", "DDX3Y", "USP9Y", "UTY", "ZFY")

# Reshape to work well with plot functions:
no_Xist_no_Y_plots <- reshape(no_Xist_no_Y_counts,
                              varying = gene_names_short,
                              drop = c("AR", "AMELY", "EIF1AY", "KDM5D",
                                       "NLGN4Y", "PRKKY", "TMSB4Y", "SRY",
                                       "status_XIST", "status_Y", "assigned",
                                       "survival", "followed"),
                              v.names = "TPM_counts",
                              timevar = "gene_names_short",
                              times = gene_names_short,
                              ids = row.names(no_Xist_no_Y_counts),
                              direction = "long")

yes_Xist_yes_Y_plots <- reshape(yes_Xist_yes_Y_counts,
                              varying = gene_names_short,
                              drop = c("AR", "AMELY", "EIF1AY", "KDM5D",
                                       "NLGN4Y", "PRKKY", "TMSB4Y", "SRY",
                                       "status_XIST", "status_Y", "assigned",
                                       "survival", "followed"),
                              v.names = "TPM_counts",
                              timevar = "gene_names_short",
                              times = gene_names_short,
                              ids = row.names(yes_Xist_yes_Y_counts),
                              direction = "long")

# Plots

# with Plotly - these don't look quite right yet

fig1 <- no_Xist_no_Y_plots %>%
  group_by(id) %>%
  plot_ly(
    x = ~gene_names_short,
    y = ~TPM_counts,
    type = "scatter",
    mode = "lines+markers"
  )

fig1 <- fig1 %>%
  layout(
    xaxis = list(
      title = "gene"
    ),
    yaxis = list(
      type = "log",
      title = "TPM",
      zeroline = FALSE
    )
  )

fig1

fig2 <- yes_Xist_yes_Y_plots %>%
  group_by(id) %>%
  plot_ly(
    x = ~gene_names_short,
    y = ~TPM_counts,
    type = "scatter",
    mode = "lines+markers"
  )

fig2 <- fig2 %>%
  layout(
    xaxis = list(
      title = "gene"
    ),
    yaxis = list(
      type = "log",
      title = "TPM",
      zeroline = FALSE
    )
  )

fig2

# with ggplot

ggplot(data = no_Xist_no_Y_plots, aes(x = gene_names_short, y = TPM_counts,
                                       group = id, color = id)) +
       scale_color_discrete(guide = "none") +
       geom_point() + geom_line() + scale_y_continuous(trans = "log10")

ggplot(data = yes_Xist_yes_Y_plots, aes(x = gene_names_short, y = TPM_counts,
                                        group = id, color = id)) +
       scale_color_discrete(guide = "none") +
       geom_point() + geom_line() + scale_y_continuous(trans = "log10")


# Consider retaining an additional column for "assigned" = M/F
not_inferred_plots <- reshape(not_inferred,
                              varying = gene_names,
                              drop = c("SRY", "inferred", "assigned",
                                       "survival", "followed"),
                              v.names = "TPM_counts",
                              timevar = "gene_names",
                              times = gene_names,
                              ids = row.names(not_inferred),
                              direction = "long")

fig1 <- not_inferred_plots %>%
  group_by(id) %>%
  plot_ly(
    x = ~gene_names,
    y = ~TPM_counts,
    type = "scatter",
    mode = "lines+markers"
  )

fig1 <- fig1 %>%
  layout(
    xaxis = list(
      title = "gene"
    ),
    yaxis = list(
      type = "log",
      title = "TPM",
      zeroline = FALSE
    )
  )

fig1

# different ways to specify the layout here
fig <- layout(fig, yaxis = list(type = "log"))

# Not interactive, but ggplot is giving something
# more like what I think we have in mind:
ggplot(data = not_inferred_plots, aes(x = gene_names, y = TPM_counts,
                                      group = id)) +
  geom_point() + geom_line() + scale_y_continuous(trans = "log10")


not_inferred_plots2 <- reshape(not_inferred,
                               varying = gene_names_short,
                               drop = c("AR", "AMELY", "EIF1AY", "KDM5D",
                                        "NLGN4Y", "PRKKY", "TMSB4Y", "SRY",
                                        "inferred", "assigned",
                                        "survival", "followed"),
                               v.names = "TPM_counts",
                               timevar = "gene_names_short",
                               times = gene_names_short,
                               ids = row.names(not_inferred),
                               direction = "long")

fig2 <- not_inferred_plots2 %>%
  group_by(id) %>%
  plot_ly(
    x = ~gene_names_short,
    y = ~TPM_counts,
    type = "scatter",
    mode = "lines+markers"
  )

fig2 <- fig2 %>%
  layout(
    xaxis = list(
      title = "gene"
    ),
    yaxis = list(
      type = "log",
      title = "TPM",
      zeroline = FALSE
    )
  )

fig2


write_tsv(
  not_inferred2 %>% rownames_to_column(),
  "TCGA_LIHC_inferences-counts.tsv",
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

ggplot(data = not_inferred_plots2, aes(x = gene_names_short, y = TPM_counts,
                                       group = id, color = id)) +
  scale_color_discrete(guide = "none") +
  geom_point() + geom_line() + scale_y_continuous(trans = "log10")
