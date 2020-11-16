# Add standard names and values for standardized metadata fields
# Add severity_standardized and days_since_symptom_onset_standardized

library(SeuratDisk)
library(data.table)

all_meta_path <- file.path(Sys.getenv("DATA_DIR"), "metadata", "all_meta.tsv")
h5seuratDir <- file.path(Sys.getenv("DATA_DIR"), "h5seurat")

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

    merged <- merge(metadata, meta[, .(patient,
                                       sample,
                                       sex,
                                       age,
                                       tissue,
                                       disease_status,
                                       disease_severity,
                                       days_since_symptom_onset,
                                       days_since_symptom_onset_estimated,
                                       days_since_hospitalization,
                                       days_since_hospitalization_estimated)],
                    by.x = sample_source_field,
                    by.y = "sample",
                    sort = FALSE)

    hmeta <- hfile[["meta.data"]]
    if (patient_source_field %in% names(hmeta)) {
      hmeta[["patient"]] <- hmeta[[patient_source_field]]
    } else if (!"patient" %in% names(hmeta)) {
      hmeta[["patient"]] <- merged$patient
    }
    hmeta[["sample"]] <- hmeta[[sample_source_field]]

    for (fieldname in c("tissue",
                        "sex",
                        "age",
                        "disease_status",
                        "disease_severity",
                        "days_since_symptom_onset",
                        "days_since_symptom_onset_estimated",
                        "days_since_hospitalization",
                        "days_since_hospitalization_estimated")) {
      fieldname_standard <- paste0(fieldname, "_standard")
      if (!fieldname_standard %in% names(hmeta)) {
        if (fieldname %in% names(hmeta)) {
          hmeta[[fieldname_standard]] <- merged[[paste0(fieldname, ".y")]]
        } else {
          hmeta[[fieldname_standard]] <- merged[[fieldname]]
        }
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
  "meckiff_2020")
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
  "orig.donor"
)

x <- mapply(merge_standard_fields,
       datasetname = datasets,
       patient_source_field = patient_source_fields,
       USE.NAMES = TRUE)
