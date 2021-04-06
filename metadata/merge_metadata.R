# Add standard names and values for standardized metadata fields
# Add severity_standardized and days_since_symptom_onset_standardized

library(SeuratDisk)
library(data.table)

all_meta_path <- file.path(Sys.getenv("DATA_DIR"), "metadata", "all_meta.tsv")
h5seuratDir <- file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed")

all_meta <- fread(all_meta_path)

merge_standard_fields <- function(datasetname, patient_source_field) {

  hfile <- Connect(file.path(h5seuratDir, paste0(datasetname, "_processed.h5seurat")), mode = "r+")
  meta <- all_meta[dataset == datasetname]
  sample_source_field <- unique(meta$sample_source_field)

  res <- try({
    metadataInfo <- hfile[["meta.data"]]$ls()
    metadataList <- lapply(metadataInfo$name, function(x) {
      d <- hfile[[paste0("meta.data/", x)]]
      if ("H5D" %in% class(d)) {
        d <- d[]
      } else {
        d <- d[["values"]][]
      }
      dt <- data.table(x = d)
      setnames(dt, "x", x)
      return(dt)
    })
    metadata <- Reduce(cbind, metadataList)
    # grant sample ids are read as integers, resulting in merge error
    metadata[[sample_source_field]] <- as.character(metadata[[sample_source_field]])

    # Heming_2021 does not come with patient/sample ids, so need to merge them in from map
    if (datasetname == "heming_2021" & !"patient" %in% names(metadata)) {
      map_heming <- fread(file.path(Sys.getenv("DATA_DIR"), "metadata", "supplemental", "annotation_patients.csv"), header = TRUE)
      metadata <- merge(metadata, map_heming, by.x = "barcode", by.y = "V1", all.x = TRUE, all.y = FALSE, sort = FALSE)
      setnames(metadata, "x", "patient_id")
      hfile[["meta.data"]][["patient_id"]] <- metadata$patient_id
    }

    merged <- merge(metadata, meta[, .(patient,
                                       sample,
                                       tissue,
                                       sample_type,
                                       sample_type_note,
                                       race_reported,
                                       sex,
                                       age,
                                       covid19,
                                       disease_status,
                                       disease_severity,
                                       days_since_symptom_onset,
                                       days_since_symptom_onset_estimated,
                                       days_since_hospitalization,
                                       days_since_hospitalization_estimated)],
                    by.x = sample_source_field,
                    by.y = "sample",
                    sort = FALSE,
                    all.x = TRUE,
                    all.y = FALSE)

    hmeta <- hfile[["meta.data"]]
    if (patient_source_field %in% names(hmeta)) {
      hmeta[["patient"]] <- hmeta[[patient_source_field]]
    } else if (!"patient" %in% names(hmeta)) {
      hmeta[["patient"]] <- merged$patient
    }
    hmeta[["sample"]] <- hmeta[[sample_source_field]]

    for (fieldname in c("tissue",
                        "sample_type",
                        "sample_type_note",
                        "race_reported",
                        "sex",
                        "age",
                        "covid19",
                        "disease_status",
                        "disease_severity",
                        "days_since_symptom_onset",
                        "days_since_symptom_onset_estimated",
                        "days_since_hospitalization",
                        "days_since_hospitalization_estimated")) {
      fieldname_standard <- paste0(fieldname, "_standard")
      if (fieldname_standard %in% names(hmeta)) {
        message("Overwriting ", fieldname_standard)
        hmeta$link_delete(fieldname_standard)
      }
      if (fieldname %in% names(hmeta)) {
        hmeta[[fieldname_standard]] <- merged[[paste0(fieldname, ".y")]]
      } else {
        hmeta[[fieldname_standard]] <- merged[[fieldname]]
      }

    }

    TRUE
  }, silent = TRUE)

  hfile$close_all()
  return(res)
}

datasets <- c(
  "wilk_2020",
  "chua_2020",
  "silvin_2020",
  "wen_2020",
  "lee_2020",
  "liao_2020",
  "zhu_2020",
  "arunachalam_2020",
  "yu_2020",
  "su_2020",
  "meckiff_2020",
  "schulte-schrepping_2020",
  "kusnadi_2021",
  "bacher_2020",
  "grant_2021",
  "stephenson_2021",
  "heming_2021",
  "trump_2020",
  "yao_2021",
  "combes_2021",
  "bost_2021")
patient_source_fields <- c(
  "Donor",
  "patient",
  "Characteristics.individual.",
  "patient",
  "Patient",
  "sample_new",
  "subject",
  "sample_name",
  "patient",
  "patient",
  "orig.donor",
  "donor",
  "orig.donor",
  "sample",
  "Patient",
  "Patient ID",
  "patient_id",
  "sample",
  "patient",
  "sample",
  "subject_id"
)

x <- mapply(merge_standard_fields,
            datasetname = datasets,
            patient_source_field = patient_source_fields,
            USE.NAMES = TRUE)

