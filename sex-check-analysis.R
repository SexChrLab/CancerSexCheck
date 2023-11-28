# WORKFLOW for CIREN TCGA SEX chromosomes vs Survival

# Import libraries to use
library(tidyverse)
library(ggplot2)
library(plotly)

# Read in data files of interest/to use
setwd("C:/Users/scmassey/Documents/CIREN/2023-Fall-MW/")
counts <- read.delim("data/TCGA_LIHC_TPM.tsv", row.names = 1)
metadf <- read.delim("data/TCGA_LIHC_META.tsv", row.names = 1)

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

for (j in seq_along(length(id_list))) {
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

sex_check <- data.frame(matrix(0, length(c_uuid), 2),
                        row.names = c_uuid)
colnames(sex_check) <- c("inferred", "assigned")

for (xid in c_uuid){
  kid <- id_sets[[xid]]
  if (all(new_counts[kid, "XIST"] > 1.0,
          new_counts[kid, "DDX3Y"] < 1.0,
          new_counts[kid, "USP9Y"] < 1.0,
          new_counts[kid, "UTY"] < 1.0,
          new_counts[kid, "ZFY"] < 1.0)) {
    sex_check[xid, "inferred"] <- "XX"
  } else if (new_counts[kid, "XIST"] < 1.0 && any(new_counts[kid, "DDX3Y"] > 1.0,
                                                  new_counts[kid, "USP9Y"] > 1.0,
                                                  new_counts[kid, "UTY"] > 1.0,
                                                  new_counts[kid, "ZFY"] > 1.0)) {
    sex_check[xid, "inferred"] <- "XY"
  } else {
    sex_check[xid, "inferred"] <- "no_inference"
  }
  sex_check[xid, "assigned"] <- metadf[xid, "mf_list"]
}

# note that using ids, we should be able to handle the double observations of
# repeated samples, esp with all/any, but might want to tread carefully there
# >>> NOPE! GOT WARNINGS ABOUT THIS..
# WILL NEED TO RETHINK HOW TO HANDLE THE MULTIPLE SAMPLE CASES HERE

# AND NEXT, look at survival info and/or disease stage.
# Q: add as additional columns to the sex_check df?
#    or keep separate as e.g., "outcomes" and pair up by row ids?
# A: look at what plot functions of interest expect as inputs and fit that


###################################################
###################################################
# OLD PLOT CODE USED - REVISE AFTER ISOLATING THE DATA WE WANT TO LOOK AT.
plotable_counts <- reshape(new_counts,
                           varying = xchr_gnames,
                           v.names = "counts",
                           timevar = "gene_names",
                           times = xchr_gnames,
                           direction = "long")

fig <- plotable_counts %>%
  plot_ly(
    x = ~gene_names,
    y = ~counts,
    type = "violin",
    box = list(
      visible = TRUE
    ),
    meanline = list(
      visible = TRUE
    ),
    x0 = xchr_gnames
  )

fig <- fig %>%
  layout(
    xaxis = list(
      title = "gene"
    ),
    yaxis = list(
      title = "TPM",
      zeroline = FALSE
    )
  )

fig

new_counts2 <- data.frame(t(subset(ychr_counts,
                                   select = -c(gene_name, gene_type))))
colnames(new_counts2) <- ychr_counts[, 1]

plotable_counts2 <- reshape(new_counts2,
                            varying = ychr_gnames,
                            v.names = "counts",
                            timevar = "gene_names",
                            times = ychr_gnames,
                            direction = "long")

fig <- plotable_counts2 %>%
  plot_ly(
    x = ~gene_names,
    y = ~counts,
    type = "violin",
    box = list(
      visible = TRUE
    ),
    meanline = list(
      visible = TRUE
    ),
    x0 = ychr_gnames
  )

fig <- fig %>%
  layout(
    xaxis = list(
      title = "gene"
    ),
    yaxis = list(
      title = "TPM",
      zeroline = FALSE
    )
  )

fig