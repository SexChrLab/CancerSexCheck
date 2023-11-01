# WORKFLOW for CIREN TCGA SEX chromosomes vs Survival

# Import libraries to use
library(ggplot2)
library(plotly)

# Read in data file(s) of interest/to use
counts = read.delim(filename, row.names = 1) # for use from saved file
# counts = data.frame(tpm_df, row.names = 1) # for testing - already in workspace

# plot violins (with jitter) of: 
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
#     TSPY   (may be challenging due to multicopy) # LEFT OUT FOR NOW -  
#                                                  # NOT SURE WHICH ONE/S TO USE 
#                                                  # ... AND THERE ARE 10+
# 
#   chromosome X
#     XIST   ENSG00000229807
#     AR     ENSG00000169083

# make one each for reported males and one for reported females (with/without log transformation)
# start with ALL together - curious if obvious bimodality

ychr_genes = c('ENSG00000099721','ENSG00000067048','ENSG00000198692','ENSG00000012817',
               'ENSG00000165246','ENSG00000099725','ENSG00000154620','ENSG00000114374',
               'ENSG00000183878','ENSG00000067646','ENSG00000184895')
ychr_gnames= c('AMELY','DDX3Y','EIF1AY','KDM5D','NLGN4Y','PRKY','TMSB4Y',
               'USP9Y','UTY','ZFY','SRY')

xchr_genes = c('ENSG00000229807','ENSG00000169083')
Xchr_gnames= c('XIST','AR')

ychr_counts <- counts[ychr_genes,]
xchr_counts <- counts[xchr_genes,]

## may want a version of new_counts for x and y?:
new_counts1 <- data.frame(t(subset(xchr_counts, select = -c(gene_name,gene_type))))
colnames(new_counts1) <- xchr_counts[, 1]

new_counts2 <- data.frame(t(subset(ychr_counts, select = -c(gene_name,gene_type))))
colnames(new_counts2) <- ychr_counts[, 1]

new_counts <- data.frame(c(new_counts1,new_counts2))
row.names(new_counts) <- row.names(new_counts1)

################################################
###---CHECK FOR REPEATS AND LOOK AT VALUES---###
################################################

# Return list of row labels (ids): 
id_list <- rownames(new_counts)

# Extract file and case uuid's from combo:
f_uuid <- substr(id_list,2,37)
c_uuid <- substr(id_list,39,74)

# Organize ids by case uuid (for viewing)
id_sets <- split(id_list,c_uuid)

repeats_df <- data.frame()
rep_cases_list <- list()

for (x in c_uuid){
    if (length(id_sets[[x]])>1){
      # whatdo
      rows_keep <- id_sets[[x]]
      
      # PUTTING THESE INTO A SEPARATE DF TO *LOOK* AT:
      repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[1],])
      repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[2],])
      
      rep_cases_list <- c(rep_cases_list, substr(rows_keep[1],39,74))
    }
}

# WRITE REPEATS DATAFRAME TO TSV FILE:
write_tsv(
  repeats_df,
  "TCGA_LIHC_TPM_Repeated_Counts.tsv",
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)
###################################################
###--END - CHECK FOR REPEATS AND LOOK AT VALUES-###
###################################################

### FINALLY, SEX CHECK OF THE SAMPLES!

# Pseudo code to get started:

# first, look at values of XIST and DDX3Y in each observation (by id)
# and assign a category - assumed XX, assumed XY, or neither

##  XIST > 1.0 && DDX3Y < 1.0 => infer XX -> F in metadata?
##  XIST < 1.0 && DDX3Y > 1.0 => infer XY -> M in metadata?
##  ALL ELSE.. => something interesting might be going on!
##  (e.g., both low might suggest loss of a Y from an XY individual)

# note that using ids, we should be able to handle the double observations
# of repeated samples, but might want to tread carefully there.

# put that categorization in a new df with same combo uuids for the rownames

# then, take the ids x categories df and add M/F from corresponding metadata

# lastly, can start doing fun things - plots,
# comparing with survival/disease stage, etc.


###################################################
###################################################
# OLD PLOT CODE - WILL REVISE AFTER ISOLATING THE DATA WE WANT TO LOOK AT.
plotable_counts <- reshape(new_counts,
               varying = Xchr_gnames,
               v.names = 'counts',
               timevar = "gene_names",
               times = Xchr_gnames,
               direction = 'long')

# #vpy <- ggplot(ychr_counts)
# #vpx <- ggplot(xchr_counts)
# 
# fig <- plot_ly(data = new_counts, x = ~XIST, y = ~AR)
# fig
# 
# ## THis WORKED YESTERDAY 9/26 BUT ISN'T WORKING NOW - NOT SURE WHY...
# fig <- new_counts %>%
#   plot_ly(
#     y = ~XIST,
#     type = 'violin',
#     box = list(
#       visible = T
#     ),
#     meanline = list(
#       visible = T
#     ),
#     x0 = 'XIST',
#   )
# 
# fig <- fig %>%
#   layout(
#     yaxis = list(
#       title = "",
#       zeroline = F
#     )
#   )
# 
# fig
# 
# fig <- new_counts %>%
#   plot_ly(
#     y = ~AR,
#     type = 'violin',
#     box = list(
#       visible = T
#     ),
#     meanline = list(
#       visible = T
#     ),
#     x0 = 'AR',
#   )
# 
# fig <- fig %>%
#   layout(
#     yaxis = list(
#       title = "",
#       zeroline = F
#     )
#   )
# 
# fig

## THIS DOES WORK THOUGH! 
fig <- plotable_counts %>%
  plot_ly(
    x = ~gene_names,
    y = ~counts,
    type = 'violin',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    x0 = Xchr_gnames
  )

fig <- fig %>%
  layout(
    xaxis = list(
      title = "gene"
    ),
    yaxis = list(
      title = "TPM",
      zeroline = F
    )
  )

fig


# right now we mostly want to see the distribution of these to decide what's useful; 
# may also want a reference gene (or two, three..) that we expect to be unaffected




new_counts2 <- data.frame(t(subset(ychr_counts, select = -c(gene_name,gene_type))))
colnames(new_counts2) <- ychr_counts[, 1]

plotable_counts2 <- reshape(new_counts2,
                           varying = ychr_gnames,
                           v.names = 'counts',
                           timevar = "gene_names",
                           times = ychr_gnames,
                           direction = 'long')

fig <- plotable_counts2 %>%
  plot_ly(
    x = ~gene_names,
    y = ~counts,
    type = 'violin',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
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
      zeroline = F
    )
  )

fig

