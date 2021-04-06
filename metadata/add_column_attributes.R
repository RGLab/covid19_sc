library(SeuratDisk)
library(hdf5r)
all_files <- list.files(file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed"), pattern = "_processed.h5seurat")
# all_files <- all_files[grepl("stephenson|grant|schrepping", all_files)]

makeLabel <- function(x) {
  x <- gsub("_standard$", "", x)
  c <- strsplit(x, "_")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}
ConnectToDataset <- function(dataset) {
  Connect(file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed", paste0(dataset, "_processed.h5seurat")), mode = "r+")
}

# Add attributes for standard fields
for (file in all_files) {
  message("----- ", file, " -----")
  filepath <- file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed", file)

  hfile <- Connect(filepath, mode = "r+")
  hmeta <- hfile[["meta.data"]]

  res <- try({
    # Add standard metadata fields
    for (column in names(hmeta)) {
      h5attr(hmeta[[column]], "standard") <- grepl("_standard$", column)
      h5attr(hmeta[[column]], "color_option") <- grepl("_standard$", column)
      h5attr(hmeta[[column]], "filter_option") <- grepl("_standard$", column)
      h5attr(hmeta[[column]], "label") <- ifelse(grepl("_standard$", column), makeLabel(column), column)
    }

    # Update analysis columns
    for (level in c("1", "2", "3")) {
      column <- hmeta[[paste0("predicted.celltype.l", level)]]
      h5attr(column, "standard") <- TRUE
      h5attr(column, "label") <- paste0("Predicted Cell Type (level ", level, ")")
      h5attr(column, "filter_option") <- TRUE
      h5attr(column, "color_option") <- TRUE
      column$close()

      column <- hmeta[[paste0("predicted.celltype.l", level, ".score")]]
      h5attr(column, "standard") <- TRUE
      h5attr(column, "label") <-  paste0("Predicted Cell Type (level ", level, ") score")
      h5attr(column, "filter_option") <- FALSE
      h5attr(column, "color_option") <- TRUE
      column$close()
    }

    column <- hmeta[["mapping.score"]]
    h5attr(column, "standard") <- TRUE
    h5attr(column, "label") <- "Seurat Mapping Score"
    h5attr(column, "filter_option") <- TRUE
    h5attr(column, "color_option") <- TRUE
    column$close()

    column <- hmeta[["rnaClusterID"]]
    h5attr(column, "standard") <- TRUE
    h5attr(column, "label") <- "Louvain Clusters (RNA)"
    h5attr(column, "filter_option") <- FALSE
    h5attr(column, "color_option") <- TRUE
    column$close()

    column <- hmeta[["patient"]]
    h5attr(column, "standard") <- TRUE
    h5attr(column, "label") <- "Donor ID"
    h5attr(column, "filter_option") <- FALSE
    h5attr(column, "color_option") <- TRUE
    column$close()

    column <- hmeta[["sample"]]
    h5attr(column, "standard") <- TRUE
    h5attr(column, "label") <- "Sample ID"
    h5attr(column, "filter_option") <- FALSE
    h5attr(column, "color_option") <- TRUE
    column$close()

    # Fix any non-default color options for standard fields
    h5attr(hmeta[["sample_type_note_standard"]], "color_option") <- FALSE
    h5attr(hmeta[["race_reported_standard"]], "color_option") <- length(unique(hmeta[["race_reported_standard"]][])) > 1
    h5attr(hmeta[["days_since_symptom_onset_estimated_standard"]], "color_option") <- FALSE
    h5attr(hmeta[["days_since_hospitalization_estimated_standard"]], "color_option") <- FALSE

    # Fix any non-default filter options for standard fields
    h5attr(hmeta[["race_reported_standard"]], "filter_option") <- length(unique(hmeta[["race_reported_standard"]][])) > 1

    # Fix non-default labels
    h5attr(hmeta[["days_since_symptom_onset_estimated_standard"]], "label") <- "Days Since Sympom Onset Imputed"
    h5attr(hmeta[["days_since_hospitalization_estimated_standard"]], "label") <- "Days Since Hospitalization Imputed"

  })
  if ("try-error" %in% class(res)) message(res)

  hfile$close_all()

}

# Now go through individual datasets
useLabel <- function(hfile, oldLabel, newLabel = NULL) {
  if (is.null(newLabel)) newLabel <- oldLabel
  column <- hfile[["meta.data"]][[oldLabel]]
  h5attr(column, "label") <- newLabel
  h5attr(column, "color_option") <- TRUE
  column$close()
}





# ----- wilk -----
hfile <- ConnectToDataset("wilk_2020")
useLabel(hfile, "Admission.level", "Admission Level")
useLabel(hfile, "Ventilated")
useLabel(hfile, "cell.type.coarse", "Cell Type (coarse)")
useLabel(hfile, "cell.type", "Cell Type (fine)")
hfile$close_all()
# ----- chua -----
hfile <- ConnectToDataset("chua_2020")
useLabel(hfile, "celltype", "Cell Type")
useLabel(hfile, "celltype_sub", "Cell Type (sub)")
useLabel(hfile, "severity", "Severity (reported)")
useLabel(hfile, "virus_pos", "SARS-CoV-2 Positive")
hfile$close_all()
# ----- silvin -----
hfile <- ConnectToDataset("silvin_2020")
useLabel(hfile, "Factor.Value.sampling.time.point.", "Sampling Time Point")
useLabel(hfile, "Factor.Value.clinical.history.", "Severity")
hfile$close_all()
# ----- wen ------

# ----- lee -----

# ----- liao -----

# ----- zhu -----
hfile <- ConnectToDataset("zhu_2020")
useLabel(hfile, "cell_type", "Cell Type")
hfile$close
# ----- arunachalam -----
hfile <- ConnectToDataset("arunachalam_2020")
# Standard for CITE-seq
column <- hfile[["meta.data"]][["adtClusterID"]]
h5attr(column, "standard") <- TRUE
h5attr(column, "label") <- "Louvain Clusters (Protein)"
h5attr(column, "filter_option") <- FALSE
h5attr(column, "color_option") <- TRUE
column$close()
column <- hfile[["meta.data"]][["wsnnClusterID"]]
h5attr(column, "standard") <- TRUE
h5attr(column, "label") <- "Louvain Clusters (WNN)"
h5attr(column, "filter_option") <- FALSE
h5attr(column, "color_option") <- TRUE
column$close()
hfile$close_all()

# ----- yu -----
hfile <- ConnectToDataset("yu_2020")
useLabel(hfile, "description", "Description")
hfile$close_all()

# ----- su -----
hfile <- ConnectToDataset("su_2020")
column <- hfile[["meta.data"]][["adtClusterID"]]
h5attr(column, "standard") <- TRUE
h5attr(column, "label") <- "Louvain Clusters (Protein)"
h5attr(column, "filter_option") <- FALSE
h5attr(column, "color_option") <- TRUE
column$close()
column <- hfile[["meta.data"]][["wsnnClusterID"]]
h5attr(column, "standard") <- TRUE
h5attr(column, "label") <- "Louvain Clusters (WNN)"
h5attr(column, "filter_option") <- FALSE
h5attr(column, "color_option") <- TRUE
column$close()
useLabel(hfile, "batch", "Batch")
hfile$close_all()

# ----- meckiff -----
# TCR, TRA, TRB, Clonotype size, Clonotype proportion can be visualized
hfile <- ConnectToDataset("meckiff_2020")
useLabel(hfile, "clonotype.tag", "Clonotype")
useLabel(hfile, "TRA.tag", "TRA")
useLabel(hfile, "TRB.tag", "TRB")
useLabel(hfile, "TCR.tag", "TCR")
useLabel(hfile, "clon.proportion.tag", "Clonotype proportion")
hfile$close_all()
# ----- kusnadi -----
hfile <- ConnectToDataset("kusnadi_2021")
useLabel(hfile, "orig.asthma", "Asthma")
useLabel(hfile, "orig.diabetes", "Diabetes")
useLabel(hfile, "orig.severity_score", "Severity Score")
useLabel(hfile, "orig.virus", "Virus")
useLabel(hfile, "clonotype.tag", "Clonotype")
useLabel(hfile, "TRA.tag", "TRA")
useLabel(hfile, "TRB.tag", "TRB")
useLabel(hfile, "TCR.tag", "TCR")
useLabel(hfile, "clon.proportion.tag", "Clonotype proportion")
hfile$close_all()

# ----- bacher ------
hfile <- ConnectToDataset("bacher_2020")
useLabel(hfile, "batch", "Batch")
useLabel(hfile, "new_cluster_names", "Cell Type")
useLabel(hfile, "raw_clonotype_id", "Clonotype")
hfile$close_all()

# ----- schulte_schrepping -----
hfile <- ConnectToDataset("schulte-schrepping_2020")
useLabel(hfile, "age", "Age Group")
useLabel(hfile, "blueprint.labels", "Blueprint Labels")
useLabel(hfile, "cluster_labels_res.0.4", "Cluster Labels (res 0.4)")
useLabel(hfile, "cluster_labels_res.0.8", "Cluster Labels (res 0.8)")
useLabel(hfile, "cohort", "Cohort")
useLabel(hfile, "hemato.labels", "Hemato Labels")
useLabel(hfile, "hpca.labels", "HPCA Labels")
useLabel(hfile, "id.celltype", "Cell Type")
useLabel(hfile, "immune.labels", "Immune Labels")
useLabel(hfile, "monaco.labels", "Monaco Labels")
useLabel(hfile, "outcome", "Outcome")
useLabel(hfile, "disease_stage", "Disease Stage")
useLabel(hfile, "who_per_sample", "WHO score")
hfile$close_all()

# ----- grant -----
hfile <- ConnectToDataset("grant_2021")
useLabel(hfile, "Day.after.intubation", "Day After Intubation")
hfile$close_all()

# ----- stephenson -----
hfile <- ConnectToDataset("stephenson_2021")
useLabel(hfile, "Collection_Day", "Collection Day")
useLabel(hfile, "Outcome")
useLabel(hfile, "Site")
useLabel(hfile, "Swab_result", "Swab Result")
useLabel(hfile, "initial_clustering", "Cell Type (initial)")
useLabel(hfile, "full_clustering", "Cell Type (full)")
# Standard for CITE-seq
column <- hfile[["meta.data"]][["adtClusterID"]]
h5attr(column, "standard") <- TRUE
h5attr(column, "label") <- "Louvain Clusters (Protein)"
h5attr(column, "filter_option") <- FALSE
h5attr(column, "color_option") <- TRUE
column$close()
column <- hfile[["meta.data"]][["wsnnClusterID"]]
h5attr(column, "standard") <- TRUE
h5attr(column, "label") <- "Louvain Clusters (WNN)"
h5attr(column, "filter_option") <- FALSE
h5attr(column, "color_option") <- TRUE
column$close()
hfile$close_all()

# ----- heming -----
hfile <- ConnectToDataset("heming_2021")
useLabel(hfile, "celltype", "Cell Type")
useLabel(hfile, "diagnosis", "Diagnosis")
useLabel(hfile, "diagnosis", "T-Cell Subtype")
hfile$close_all()

# ------ trump -----
hfile <- ConnectToDataset("trump_2020")
useLabel(hfile, "cardiac_condition", "Cardiac Condition")
useLabel(hfile, "medication", "Medication")
hfile$close_all()

# ----- yao -----
hfile <- ConnectToDataset("yao_2021")
useLabel(hfile, "celltype5", "Cell Type")
hfile$close_all()

# ----- combes -----

# ----- bost ------
hfile <- ConnectToDataset("bost_2021")
useLabel(hfile, "clinical_outcome", "Outcome")
hfile$close_all()


