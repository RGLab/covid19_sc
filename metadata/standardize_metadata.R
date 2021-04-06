# Standardize metadata. Reads in metadata provided with the dataset
# and standardizes fields for all datasets, writing the result to
# a file. Looks for metadata written to a tsv file in metadataDir.
# These tsv files can be written using write_metadata.R.

library(data.table)
setDTthreads(4)
library(readxl)

metadataDir <- file.path(Sys.getenv("DATA_DIR"), "metadata")
supplementalDir <- file.path(Sys.getenv("DATA_DIR"), "covid19_sc", "metadata", "supplemental")
files <- list.files(metadataDir, pattern = ".tsv")
files <- files[!grepl("all_meta", files)]
mList <- lapply(file.path(metadataDir, files), fread)
names(mList) <- gsub(".tsv", "", files)
names <- lapply(mList, names)
names(names) <- files
files


# tissue, sample_type, sample_type_sub
# Look for race/ethnicity where available

# dataset, sample_source, patient, sample, tissue, sex (male, female, other), age, disease_status (COVID-19, healthy, flu), disease_severity (mild, severe)

# ----- wilk_2020 -----
# severe = hospitalized
wilk <- unique(mList$wilk_2020[, .(Donor, Donor.full, DPS, Sex, Status, Admission.level)])[, .(
  dataset = "wilk_2020",
  sample_source_field = "Donor.full",
  patient = Donor,
  sample = Donor.full,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = NA,
  sex = ifelse(Sex == "M", "male", "female"),
  age = NA,
  race = NA,
  covid19 = (Status == "COVID"),
  disease_status = ifelse(Status == "COVID", "COVID-19", "healthy"),
  disease_severity = ifelse(Admission.level == "ICU", "severe",
                            ifelse(Admission.level == "Floor", "moderate", "healthy")),
  days_since_symptom_onset = DPS,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)]



# ----- chua -----
# TODO: add more informative timepoint info
# moderate vs severe per WHO guidelines (all were hospitalized)
# Metadata from https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-020-0602-4/MediaObjects/41587_2020_602_MOESM1_ESM.pdf
meta_chua <- list(
  "BIH-CoV-15" = list(age = 21, sex = "male", disease_severity = "moderate"),
  "BIH-CoV-14" = list(age = 45, sex = "male", disease_severity = "moderate"),
  "BIH-CoV-16" = list(age = 50, sex = "female", disease_severity = "moderate"),
  "BIH-CoV-13" = list(age = 52, sex = "male", disease_severity = "moderate"),
  "BIH-CoV-19" = list(age = 64, sex = "male", disease_severity = "moderate"),
  "BIH-CoV-12" = list(age = 71, sex = "male", disease_severity = "moderate"),
  "BIH-CoV-17" = list(age = 73, sex = "male", disease_severity = "moderate"),
  "BIH-CoV-18" = list(age = 76, sex = "male", disease_severity = "moderate"),

  "BIH-CoV-07" = list(age = 32, sex = "male", disease_severity = "severe"),
  "BIH-CoV-11" = list(age = 61, sex = "male", disease_severity = "severe"),
  "BIH-CoV-06" = list(age = 63, sex = "male", disease_severity = "severe"),
  "BIH-CoV-05" = list(age = 75, sex = "male", disease_severity = "severe"),
  "BIH-CoV-09" = list(age = 53, sex = "female", disease_severity = "severe"),
  "BIH-CoV-03" = list(age = 52, sex = "male", disease_severity = "severe"),
  "BIH-CoV-08" = list(age = 62, sex = "female", disease_severity = "severe"),
  "BIH-CoV-04" = list(age = 73, sex = "male", disease_severity = "severe"),
  "BIH-CoV-10" = list(age = 68, sex = "male", disease_severity = "severe"),
  "BIH-CoV-01" = list(age = 50, sex = "female", disease_severity = "severe"),
  "BIH-CoV-02" = list(age = 54, sex = "male", disease_severity = "severe"),

  "BIH-Con-01" = list(age = NA, sex = NA, disease_severity = NA),
  "BIH-Con-02" = list(age = NA, sex = NA, disease_severity = NA),
  "BIH-Con-03" = list(age = NA, sex = NA, disease_severity = NA),
  "BIH-Con-04" = list(age = NA, sex = NA, disease_severity = NA),
  "BIH-Con-05" = list(age = NA, sex = NA, disease_severity = NA)
)
names <- names(meta_chua)
meta_chua <- data.table(matrix(unlist(meta_chua), nrow=length(meta_chua), byrow=T))
meta_chua$patient <- names
setnames(meta_chua, c("V1", "V2", "V3"), c("age", "sex", "disease_severity"))




# All patients were hospitalized. Severe = ICU
chua <- unique(mList$chua_2020[, .(patient, sample, infection, location)])
chua <- merge(meta_chua, chua, by = "patient", all.y = TRUE, sort = FALSE)
# Add days since symptom onset
chua[, days_since_symptom_onset := as.numeric(NA)]
chua[sample == "BIH-CoV-15_NS_1", days_since_symptom_onset := 25]
chua[sample == "BIH-CoV-15_NS_2", days_since_symptom_onset := 27]
chua[sample == "BIH-CoV-15_NS_4", days_since_symptom_onset := 30]
chua[sample == "BIH-CoV-15_NS_3", days_since_symptom_onset := 33]
chua[sample == "BIH-CoV-14_NS_1", days_since_symptom_onset := 11]
chua[sample == "BIH-CoV-14_NS_2", days_since_symptom_onset := 16]
chua[sample == "BIH-CoV-16_NS_1", days_since_symptom_onset := 10]
chua[sample == "BIH-CoV-13_NS_1", days_since_symptom_onset := 10]
chua[sample == "BIH-CoV-19_NS_1", days_since_symptom_onset := 11]
chua[sample == "BIH-CoV-12_NS_1", days_since_symptom_onset := 16]
chua[sample == "BIH-CoV-12_NS_2", days_since_symptom_onset := 18]
chua[sample == "BIH-CoV-12_NS_3", days_since_symptom_onset := 21]
chua[sample == "BIH-CoV-17_NS_1", days_since_symptom_onset := 30]
chua[sample == "BIH-CoV-18_NS_1", days_since_symptom_onset := 3]
chua[sample == "BIH-CoV-07_NS_1", days_since_symptom_onset := 8]
chua[sample == "BIH-CoV-07_NS_2", days_since_symptom_onset := 11]
chua[sample == "BIH-CoV-11_NS_1", days_since_symptom_onset := 11]
chua[sample == "BIH-CoV-06_NS_1", days_since_symptom_onset := 4]
chua[sample == "BIH-CoV-06_NS_2", days_since_symptom_onset := 7]
chua[sample == "BIH-CoV-05_NS_1", days_since_symptom_onset := 7]
chua[sample == "BIH-CoV-09_NS_1", days_since_symptom_onset := 12]
chua[sample == "BIH-CoV-03_NS_1", days_since_symptom_onset := 10]
chua[sample == "BIH-CoV-08_NS_1", days_since_symptom_onset := 8]
chua[sample == "BIH-CoV-04_BL_1", days_since_symptom_onset := 13]
chua[sample == "BIH-CoV-04_PS_1", days_since_symptom_onset := 13]
chua[sample == "BIH-CoV-04_NS_1", days_since_symptom_onset := 13]
chua[sample == "BIH-CoV-10_NS_1", days_since_symptom_onset := 20]
chua[sample == "BIH-CoV-01_BL_1", days_since_symptom_onset := 17]
chua[sample == "BIH-CoV-01_PS_1", days_since_symptom_onset := 17]
chua[sample == "BIH-CoV-01_NS_1", days_since_symptom_onset := 17]
chua[sample == "BIH-CoV-02_NS_1", days_since_symptom_onset := 7]

chua <- chua[, .(
  dataset = "chua_2020",
  sample_source_field = "sample",
  sample = sample,
  patient = patient,
  tissue = sapply(location, switch,
                  "NS" = "Nasopharynx",
                  "BL" = "Bronchus",
                  "PSB" = "Bronchus"),
  sample_type = sapply(location, switch,
                       "NS" = "Nasopharyngeal or Pooled Nasopharyngeal/Pharyngeal Swabs",
                       "BL" = "Bronchial Lavage",
                       "PSB" = "Bronchial Protected Specimen Brush"),
  sample_type_note = NA,
  race = NA,
  sex = sex,
  age = age,
  covid19 = infection == "SARS-CoV-2",
  disease_status = ifelse(infection == "SARS-CoV-2", "COVID-19", "healthy"),
  disease_severity = disease_severity,
  days_since_symptom_onset = days_since_symptom_onset,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)]

# ----- silvin -----
silvin <- unique(mList$silvin_2020[, .(
  Characteristics.individual.,
  sample,
  Characteristics.sex.,
  Characteristics.age.,
  Characteristics.disease.,
  Characteristics.clinical.history.,
  Factor.Value.sampling.time.point.)])
# Get time since hospitalization from fig 2
offset <- data.table(
  Characteristics.individual. = c("SARS-CoV2 pos Severe #1",
              "SARS-CoV2 pos Mild",
              "SARS-CoV2 pos Severe #2",
              "SARS-CoV2 pos Severe #3"),
  day0_offset = c(14,
           0,
           10,
           0)
)
silvin <- merge(silvin, offset, all.x = TRUE)

silvin <- silvin[, .(
    dataset = "silvin_2020",
    sample_source_field = "sample",
    sample = sample,
    patient = Characteristics.individual.,
    tissue = "Blood",
    sample_type = "Whole Blood",
    sample_type_note = "Red blood cells depleted whole blood",
    race = NA,
    sex = ifelse(Characteristics.sex. == "not available", NA, Characteristics.sex.),
    age = ifelse(Characteristics.age. == "not available", NA, as.numeric(Characteristics.age.)),
    covid19 = Characteristics.disease. == "COVID-19",
    disease_status = ifelse(Characteristics.disease. == "COVID-19", "COVID-19", "healthy"),
    disease_severity = ifelse(Characteristics.clinical.history. == "severe COVID-19",
                              "severe",
                              ifelse(Characteristics.clinical.history. == "mild COVID-19", "mild", NA)
                              ),
    days_since_symptom_onset = NA,
    days_since_hospitalization = day0_offset + Factor.Value.sampling.time.point.,
    days_since_symptom_onset_estimated = TRUE,
    days_since_hospitalization_estimated = FALSE
  )]



# ----- wen -----
# data from supplementary figure https://static-content.springer.com/esm/art%3A10.1038%2Fs41421-020-0168-9/MediaObjects/41421_2020_168_MOESM1_ESM.pdf

wen_meta <- data.table(
  patient = c("ERS1", "ERS2", "ERS3", "ERS4", "ERS5", "LRS1", "LRS2", "LRS3", "LRS4", "LRS5"),
  sex = c("male", "male", "female", "female", "female", "male", "male", "male", "female", "female"),
  disease_severity = c("moderate", "moderate", "severe", "severe", "moderate", "severe", "moderate", "moderate", "severe", "severe")
)


wen <- unique(mList$wen_2020[, .(patient, sample)])
wen <- merge(wen, wen_meta, all.x = TRUE)
wen <- wen[, .(
  dataset = "wen_2020",
  sample_source_field = "sample",
  sample = sample,
  patient = patient,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = NA,
  race = NA,
  sex = sex,
  age = NA,
  covid19 = !grepl("Healthy", patient),
  disease_status = ifelse(grepl("Healthy", patient), "healthy", "COVID-19"),
  disease_severity = ifelse(grepl("Healthy", patient), "healthy", disease_severity),
  days_since_symptom_onset = ifelse(grepl("ERS", patient), 15,
                                    ifelse(grepl("LRS", patient), 29, NA)),
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = TRUE,
  days_since_hospitalization_estimated = TRUE
)]

# ----- lee -----

# NOTE: NA cells belong with the preceding patient. merged cells.
meta_lee <- data.table(read_xlsx(file.path(supplementalDir, "lee_2020_Table_S1.xlsx")))
for (i in 1:nrow(meta_lee)) {
  if (is.na(meta_lee[i, "Patient ID"])) {
    meta_lee[i, `:=`("Patient ID" = meta_lee[i-1, "Patient ID"],
                     Age = meta_lee[i-1, "Age"],
                     Sex = meta_lee[i-1, "Sex"])
             ]
  }

}
meta_lee[, SampleName := paste0(`Sample ID`, " scRNA-seq")]

meta_lee2 <- data.table(read_xlsx(file.path(supplementalDir, "lee_2020_Table_S3.xlsx")))
meta_lee2 <- meta_lee2[grepl("nCoV", `Sample ID`), .(`Sample ID`, `Hospital \nday`)]
meta_lee2[, SampleName := paste0(`Sample ID`, " scRNA-seq")]

lee <- unique(mList$lee_2020[, .(SampleName = Patient, SampleID = Sample)])
lee <- merge(lee, meta_lee, by = "SampleName", all.x = TRUE)
lee <- merge(lee, meta_lee2, by = "SampleName", all.x = TRUE)
lee <- lee[, .(
  dataset = "lee_2020",
  sample_source_field = "Patient",
  sample = SampleName,
  patient = `Patient ID`,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = NA,
  race = NA,
  sex = ifelse(Sex == "M", "male", "female"),
  age = Age,
  covid19 = grepl("COVID", `Disease group`),
  disease_status = ifelse(grepl("influenza", `Disease group`), "influenza",
                          ifelse(grepl("COVID", `Disease group`), "COVID-19", "healthy")),
  disease_severity = ifelse(grepl("COVID", `Disease group`),
                            ifelse(grepl("severe", `Disease group`), "severe", "moderate"),
                            ifelse(grepl("influenza", `Disease group`), "severe", "healthy")),
  days_since_symptom_onset = NA,
  days_since_hospitalization = `Hospital \nday`,
  days_since_symptom_onset_estimated = TRUE,
  days_since_hospitalization_estimated = FALSE

)]

# ----- liao -----
# Disease severity: moderate vs severe
# NOTE: sex and age missing for controls
if ( !file.exists(file.path(supplementalDir, "liao_2020_meta.txt")) ) {
  # Pull in new version of metadata
  meta_liao <- fread("https://raw.githubusercontent.com/zhangzlab/covid_balf/master/meta.txt")
  fwrite(meta_liao, file.path(supplementalDir, "liao_2020_meta.txt"))
}
meta_liao <- fread(file.path(supplementalDir, "liao_2020_meta.txt"))

meta2_liao <- data.table(
  patient = c("M1", "M2", "M3", "S1", "S2", "S3", "S4", "S5", "S6"),
  age = c(36, 37, 35, 62, 66, 63, 65, 57, 46),
  sex = c("male", "female", "male", "male", "male", "male", "female", "female", "male"),
  symptom_onset_date = c("2020-01-09", "2020-01-11", "2020-01-09", "2020-01-11", "2020-01-03", "2020-01-08", "2020-01-04", "2020-01-21", "2020-01-21"),
  hospitalization_date = c("2020-01-16", "2020-01-16", "2020-01-15", "2020-01-18", "2020-01-11", "2020-01-15", "2020-01-20", "2020-01-23", "2020-01-22"),
  sampling_date = c("2020-01-20", "2020-01-20", "2020-01-22", "2020-01-22", "2020-01-21", "2020-01-22", "2020-01-29", "2020-01-29", "2020-02-02")
)
meta2_liao[, `:=`(
  symptom_onset_date = strptime(symptom_onset_date, "%Y-%m-%d"),
  hospitalization_date = strptime(hospitalization_date, "%Y-%m-%d"),
  sampling_date = strptime(sampling_date, "%Y-%m-%d")
)]


liao <- unique(mList$liao_2020[, .(sample)])
liao <- merge(liao, meta_liao, by = "sample")
liao <- merge(liao, meta2_liao, by.x = "sample_new", by.y = "patient", all.x = TRUE)
liao <- liao[, .(
  dataset = "liao_2020",
  sample_source_field = "sample",
  sample = sample,
  patient = sample_new,
  tissue = "Bronchus",
  sample_type = "Bronchoalveolar Lavage Fluid",
  sample_type_note = NA,
  race = NA,
  sex = sex,
  age = age,
  covid19 = disease == "Y",
  disease_status = ifelse(disease == "Y", "COVID-19", "healthy"),
  disease_severity = ifelse(grepl("HC", sample_new_old), "healthy",
                            ifelse(grepl("C", sample_new_old), "severe", # labeled "critical"
                                   "moderate")),
  days_since_symptom_onset = as.numeric(difftime(sampling_date, symptom_onset_date, units = "days")),
  days_since_hospitalization = as.numeric(difftime(sampling_date, hospitalization_date, units = "days")),
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = FALSE
)]

# ----- zhu -----
# age and sex missing for control subjects
# severity: mild vs severe

meta_zhu <- data.table(read_xlsx(file.path(supplementalDir, "zhu_2020_mmc2.xlsx"), skip = 2, col_names = FALSE))
names <- unlist(meta_zhu[, 1])
meta_zhu <- data.table(t(meta_zhu[, 2:ncol(meta_zhu)]))
colnames(meta_zhu) <- names
setnames(meta_zhu, "V1", "subject")
meta_zhu[, subject := gsub("IAV", "Flu", subject)]

day1_offset <- data.table(
  subject = c("COV-1", "COV-2", "COV-3", "COV-4", "COV-5", "Flu-1", "Flu-2"),
  symptom_offset = c(7, 7, 13, 10, 3, NA, NA),
  hospitalization_offset = c(4, 4, 4, 2, 1, 1, 5)
)

zhu <- unique(mList$zhu_2020[, .(batch, Stage)])[, subject := gsub("-D\\d+", "", batch)]
zhu <- merge(zhu, meta_zhu, by = "subject", all.x = TRUE)
zhu <- merge(zhu, day1_offset, by = "subject", all.x = TRUE)
zhu[, timepoint := as.numeric(gsub(".+-D", "", batch))]
zhu <- zhu[, .(
  dataset = "zhu_2020",
  sample_source_field = "batch",
  sample = batch,
  patient = subject,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = NA,
  race = NA,
  sex = tolower(Gender),
  age = Age,
  covid19 = grepl("COV", subject),
  disease_status = ifelse(grepl("COV", subject), "COVID-19",
                          ifelse(grepl("Flu", subject), "influenza", "healthy")),
  disease_severity = ifelse(grepl("COV", subject), "moderate",
                            ifelse(grepl("Flu", subject), "severe", "healthy")),
  days_since_symptom_onset = symptom_offset + timepoint - 1,
  days_since_hospitalization = hospitalization_offset + timepoint -1,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = FALSE
)]




# ----- arunachalam -----
meta_arunachalam <- fread(file.path(supplementalDir, "arunachalam_table_2.tsv"))

# Severity: severe vs moderate
arunachalam <- unique(mList$arunachalam_2020[, .(sample_name,
                                                 sex,
                                                 disease_status,
                                                 disease_severity,
                                                 age,
                                                 days_since_symptom_onset)])[, .(
  dataset = "arunachalam_2020",
  sample_source_field = "sample_name",
  sample = sample_name,
  patient = sample_name,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = "DCs enriched",
  race = NA,
  sex = ifelse(sex == "M", "male", "female"),
  age = age,
  covid19 = disease_status == "COVID-19",
  disease_status = ifelse(disease_status == "COVID-19", "COVID-19", "healthy"),
  disease_severity = ifelse(disease_status == "COVID-19", "moderate", "healthy"),
  days_since_symptom_onset = days_since_symptom_onset,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)]

# update to severe based on supplementary table ICU
arunachalam[age == 60 & sex == "female", disease_severity := "severe"]

# ----- yu -----
# NOTE:  missing severe samples
# NOTE:  PM = post mild disease (hospitalized), CM = convalescence of mild disease (no longer hospitalized)
timepoints <- data.table(
  patient = c("CM1", "CM2", "CM3", "PM1", "PM2", "PM3"),
  days_since_hospitalization = c(14, 9, 14, 14, 7, 7)
)
yu <- unique(mList$yu_2020[, .(patient, sample, sex, description)])
yu <- merge(yu, timepoints, all.x = TRUE)
yu <- yu[, .(
  dataset = "yu_2020",
  sample_source_field = "sample",
  sample = sample,
  patient = patient,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = NA,
  race = NA,
  sex = sex,
  age = NA,
  covid19 = !grepl("HD", patient),
  disease_status = ifelse(grepl("HD", patient), "healthy", "COVID-19"),
  disease_severity = ifelse(grepl("mild", description), "moderate", "healthy"),
  days_since_symptom_onset = NA,
  days_since_hospitalization = days_since_hospitalization,
  days_since_symptom_onset_estimated = TRUE,
  days_since_hospitalization_estimated = FALSE
)]


# NA timepoine for PM3CM3-Combine
yu[is.na(days_since_hospitalization), `:=`(days_since_symptom_onset_estimated = NA,
                                           days_since_hospitalization_estimated = NA)]




# ----- su -----
# severity: ICU vs no ICU. Also have info on hospital vs home.
meta_su_cov <- data.table(read_xlsx(file.path(supplementalDir, "su_2020_Table_S1.xlsx"), sheet = 2))

meta_su_cov <- meta_su_cov[, .(
  sample = gsub("INCOV0*", "", `Sample ID`),
  patient = gsub("INCOV0*", "", `Study Subject ID`),
  age = Age,
  sex = tolower(Sex),
  location = `Patient Location`,
  WHO = `Who Ordinal Scale`,
  disease_status = "COVID-19",
  timepoint = `Blood draw time point`,
  race = Race
)]
meta_su_healthy <- data.table(read_xlsx(file.path(supplementalDir, "su_2020_Table_S1.xlsx"), sheet = 3))
meta_su_healthy <- meta_su_healthy[, .(
  sample = `Sample ID`,
  patient = `Sample ID`,
  sex = tolower(sex),
  age = age,
  disease_status = "healthy",
  location = "healthy"
)]
meta_su <- rbind(meta_su_cov, meta_su_healthy, fill = TRUE)
# Fill in missing severity info for 129-2
meta_su[sample == "129-2", location := "Home (mobile phlebotomy)"]

su <- unique(mList$su_2020[, .(patient, sample)])
su <- merge(su, meta_su, by = c("sample", "patient"), all.x = TRUE)
su <- su[, .(
  dataset = "su_2020",
  sample_source_field = "sample",
  patient = patient,
  sample = sample,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = NA,
  race = race,
  sex = sex,
  age = as.numeric(age),
  covid19 = disease_status == "COVID19",
  disease_status = disease_status,
  disease_severity = ifelse(location == "Hospital", "moderate",
                            ifelse(location == "ICU", "severe",
                                   ifelse(location == "healthy", "healthy", "mild"))),
  days_since_symptom_onset = ifelse(timepoint == "T1", 7,
                                    ifelse(timepoint == "T2", 14, NA)),
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = TRUE,
  days_since_hospitalization_estimated = TRUE
)]



# ----- meckiff -----
# Severity: ICU vs ward vs no hospital (severe, moderate, mild)
# no demographics for controls
meta_meckiff <- data.table(read_xlsx(file.path(supplementalDir, "meckiff_2020_1-s2.0-S0092867420313076-mmc1.xlsx"),
                                     sheet = 2,
                                     skip = 3))

meckiff <- unique(mList$meckiff_2020[, .(orig.donor, orig.sex, orig.hospital)])
meckiff <- merge(meckiff, meta_meckiff, by.x = "orig.donor", by.y = "PATIENT ID", all.x = TRUE)
meckiff <- meckiff[, .(
  dataset = "meckiff_2020",
  sample_source_field = "orig.donor",
  patient = orig.donor,
  sample = orig.donor,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = "CD4 ARTE",
  race = ETHNICITY,
  sex = tolower(orig.sex),
  age = as.numeric(AGE),
  covid19 = !is.na(HOSPITALIZATION),
  disease_status = ifelse(is.na(HOSPITALIZATION), "healthy", "COVID-19"),
  disease_severity = ifelse(is.na(HOSPITALIZATION), "healthy",
                            ifelse(grepl("ICU", HOSPITALIZATION), "severe",
                                   ifelse(grepl("Ward", HOSPITALIZATION), "moderate", "mild"))),
  days_since_symptom_onset = `INTERVAL BETWEEN SYMPTOM ONSET TO SAMPLE COLLECTION (DAYS)`,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)]

# ----- kusnadi_2021 -----
meta_kusnadi <- data.table(read_excel(file.path(supplementalDir, "kusnadi_2021_Table_S1.xlsx"),
                                     sheet = 1,
                                     skip = 3))
meta_kusnadi[, pid := gsub("00", "P", `PATIENT ID`)]
meta_kusnadi_control <- data.table(read_excel(file.path(supplementalDir, "kusnadi_2021_Table_S1.xlsx"),
                                              sheet = 2,
                                              skip = 3))

meta_kusnadi_control[, pid := gsub("-", "", `PATIENT ID`)]

meta_kusnadi <- rbind(meta_kusnadi, meta_kusnadi_control,
                            use.names = TRUE,
                            fill = TRUE,
                            idcol = TRUE)

kusnadi <- unique(mList$kusnadi_2021[, .(
  orig.donor,
  orig.hospital,
  orig.severity,
  orig.severity_score,
  orig.severity_x,
  orig.sex,
  orig.unit,
  orig.virus,
  orig.virus2,
  ethnicity
)])

kusnadi <- merge(kusnadi, meta_kusnadi, by.x = "orig.donor", by.y = "pid", all.x = TRUE)

kusnadi <- unique(kusnadi[, .(
  dataset = "kusnadi_2021",
  sample_source_field = "orig.donor",
  patient = orig.donor,
  sample = orig.donor,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = "CD8 ARTE",
  race = ETHNICITY,
  sex = tolower(ifelse(orig.sex == "void", NA, orig.sex)),
  age = as.numeric(AGE),
  covid19 = grepl("SARS-CoV-2", `VIRUS (PCR POSITIVE)`),
  disease_status = ifelse(grepl("SARS-CoV-2", `VIRUS (PCR POSITIVE)`), "COVID-19", "healthy"),
  disease_severity = ifelse(is.na(HOSPITALIZATION), "healthy",
                            ifelse(HOSPITALIZATION == "Not admitted", "mild",
                                   ifelse(HOSPITALIZATION == "Ward", "moderate",
                                          ifelse(HOSPITALIZATION == "ICU", "severe", NA)))),
  days_since_symptom_onset = `INTERVAL BETWEEN SYMPTOM ONSET TO SAMPLE COLLECTION (DAYS)`,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)])

# ----- bacher_2020 ----
# Additional demographics info in supp table, but no info on how to map to data.
# https://www.cell.com/cms/10.1016/j.immuni.2020.11.016/attachment/11145884-7256-4dea-a99b-004757319db3/mmc1.pdf

meta_bacher <- data.table(read_excel(file.path(supplementalDir, "bacher_2020_demographics_ID.xlsx")))
meta_bacher[, sample := gsub("_RNA", "", GEO_ID...1)]

bacher <- unique(mList$bacher_2020[, .(
  diagnosis,
  sample
)])
bacher <- merge(bacher, meta_bacher, by = "sample")


bacher <- bacher[, .(
  dataset = "bacher_2020",
  sample_source_field = "sample",
  patient = sample,
  sample = sample,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = "CD4 ARTE",
  race = NA,
  sex = Gender,
  age = Age,
  covid19 = diagnosis == "COVID19",
  disease_status = ifelse(diagnosis == "COVID19", "COVID-19", "healthy"),
  disease_severity = ifelse(Classification == "Non-hospitalized", "mild",
                            ifelse(Classification == "mild-moderate", "moderate",
                                   Classification)),
  days_since_symptom_onset = `Days since symptom onset`,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)]


# ----- schulte-schrepping -----
schulte_schrepping <- unique(mList$`schulte-schrepping_2020`[, .(donor,
                                                                 sampleid_unique,
                                                                 sex,
                                                                 disease_stage,
                                                                 group_per_sample,
                                                                 days_after_onset,
                                                                 age,
                                                                 who_per_sample,
                                                                 cells)])
# NOTE: Using min of age range for age
schulte_schrepping <- schulte_schrepping[, .(
  dataset = "schulte-schrepping_2020",
  sample_source_field = "sampleid_unique",
  patient = donor,
  sample = sampleid_unique,
  tissue = "Blood",
  sample_type = ifelse(cells == "Whole_blood",
                       "Whole Blood", "Peripheral Blood Mononuclear Cells"),
  sample_type_note = ifelse(cells == "CD45-sorted frozen PBMC",
                            cells, NA),
  race = NA,
  sex = ifelse(sex == "n/a", NA,
               ifelse(sex == "M", "male",
                      ifelse(sex == "F", "female", sex))),
  age = as.numeric(sub("_\\d+", "", age)),
  covid19 = who_per_sample > 0,
  disease_status = ifelse(who_per_sample > 0, "COVID-19", "healthy"),
  disease_severity = ifelse(who_per_sample == 0, "healthy",
                            ifelse(who_per_sample < 3, "mild",
                                   ifelse(who_per_sample <= 5, "moderate",
                                          ifelse(who_per_sample > 5, "severe", NA)))),
  days_since_symptom_onset = days_after_onset,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)]


# ----- grant_2021 -----
meta_grant2 <- fread(file.path(supplementalDir, "grant_2021_supplement-metadata.csv"))
meta_grant3 <- fread(file.path(supplementalDir, "grant_supplemental.csv"))
grant <- mList$grant_2021


setnames(meta_grant2, "V1", "barcode")
grant <- merge(mList$grant_2021, meta_grant2, by.x = "_index", by.y = "barcode")
grant <- unique(grant[, .(
  patient = Patient.y,
  sample = Sample.x,
  sample_type = `Sample type`,
  covid_status = `COVID-19`
)])
grant <- merge(grant, meta_grant3,
               by.x = "patient",
               by.y = "scRNAseq_id",
               all.x = TRUE,
               all.y = FALSE)

grant <- grant[, .(
  dataset = "grant_2021",
  sample_source_field = "Sample",
  patient = patient,
  sample = sample,
  tissue = "Bronchus",
  sample_type = "Bronchoalveolar Lavage Fluid",
  sample_type_note = "sorted T cells and macrophages",
  race = NA,
  sex = tolower(sex),
  age = trunc(age),
  covid19 = covid_status,
  disease_status = ifelse(diagnosis == "Non-Pneumonia Control", "hospitalized non-COVID", diagnosis),
  disease_severity = ifelse(covid_status, "severe", NA),
  days_since_symptom_onset = NA,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = NA,
  days_since_hospitalization_estimated = NA
)]

# ----- stephenson_2021 -----
meta_stephenson <- data.table(read_excel(file.path(supplementalDir, "stephenson_2021_supplemental.xlsx"),
                                         sheet = 1))

stephenson <- unique(mList$stephenson_2021[, .(
  Age,
  Collection_Day,
  Days_from_onset,
  Sex,
  Status,
  Status_on_day_collection,
  Status_on_day_collection_summary,
  orig.ident,
  sample_id
)])

stephenson <- merge(stephenson, meta_stephenson, by.y = c("Sample ID", "Status"), by.x = c("sample_id", "Status"))

# impute days_since_onset
stephenson[Status == "Covid",
           days_since_symptom_onset_imputed := round(median(as.numeric(Days_from_onset), na.rm = TRUE)),
           c("Status", "Status_on_day_collection_summary")]

# Correct days_since_onset For Collection_Day > 0
stephenson[, days_since_symptom_onset_corrected := as.numeric(Days_from_onset) + as.numeric(gsub("D", "", Collection_Day))]
stephenson[, days_since_symptom_onset_imputed_corrected := days_since_symptom_onset_imputed + as.numeric(gsub("D", "", Collection_Day))]

stephenson <- unique(stephenson[, .(
  dataset = "stephenson_2021",
  sample_source_field = "sample_id",
  patient = `Patient ID`,
  sample = sample_id,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = NA,
  race = NA,
  sex = tolower(Sex),
  age = Age.y,
  covid19 = Status == "Covid",
  disease_status = ifelse(Status == "Covid", "COVID-19",
                          ifelse(Status == "Healthy", "healthy",
                                 ifelse(Status == "Non_covid", "hospitalized non-COVID",
                                        ifelse(Status == "LPS", "IV-LPS", Status)))),
  disease_severity = ifelse(Status_on_day_collection_summary == "Healthy", "healthy",
                            ifelse(Status_on_day_collection_summary == "Asymptomatic", "mild",
                                   ifelse(Status_on_day_collection_summary %in% c("Mild", "Moderate"), "moderate",
                                          ifelse(Status_on_day_collection_summary %in% c("Severe", "Critical"), "severe", NA)))),
  days_since_symptom_onset = ifelse(is.na(days_since_symptom_onset_corrected), days_since_symptom_onset_imputed_corrected, days_since_symptom_onset_corrected),
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = is.na(as.numeric(Days_from_onset)),
  days_since_hospitalization_estimated = TRUE
)])

# ----- heming_2021 -----
# Single-cell atlas of cerebrospinal fluid in Neuro-COVID and controls
# Only includes diagnosis (IIH, VE, MS, or COVID) and no other metadata.
map_heming <- fread(file.path(supplementalDir, "annotation_patients.csv"), header = TRUE)
setnames(map_heming,
         c("V1", "x"),
         c("barcode", "patient_id"))
meta_heming <- data.table(read_excel(file.path(supplementalDir, "heming_2021_tableS1.xlsx")))
setnames(meta_heming,
         c("COVID-19\n  severity",
                        "original\npseudonym"),
         c("severity", "patient_id_orig"))

# annotation_patients uses slightly different IDs for some of the patients.
# > setdiff(meta_heming$patient_id, map_heming$patient_id)
# [1] "MSC79670" "MSC76177" "MSC77654" "MSC25719" "MSC90896" "PST32190" "PST41540" "PST45044"
# [9] "PST85037"
# > setdiff(map_heming$patient_id, meta_heming$patient_id)
# [1] "IIH32190" "IIH41540" "IIH45044" "IIH85037" "MS25719"  "MS76177"  "MS77654"  "MS79670"
# [9] "MS90896"
meta_heming[, patient_id := gsub("MSC", "MS", patient_id_orig)]
meta_heming[, patient_id := gsub("PST", "IIH", patient_id)]


all_meta_heming <- merge(map_heming,
                         meta_heming[, c(
                           "Patient",
                           "Sex",
                           "Age",
                           "severity",
                           "patient_id"),
                           with = FALSE],
                         by = "patient_id",
                         all.x = TRUE)

heming <- mList$heming_2021
heming <- merge(heming,
                all_meta_heming,
                by = "barcode",
                all.x = TRUE)

heming <- unique(heming[, .(
  dataset = "heming_2021",
  sample_source_field = "patient_id",
  patient = patient_id,
  sample = patient_id,
  tissue = "Central Nervous System",
  sample_type = "Cerebrospinal Fluid", # cerebrospinal fluid
  sample_type_note = NA,
  race = NA,
  sex = ifelse(Sex == "f", "female", "male"),
  age = Age,
  covid19 = diagnosis == "COVID",
  disease_status = ifelse(diagnosis == "COVID", "COVID-19", diagnosis), # TODO: Not really healthy. Oher conditions.
  disease_severity = ifelse(severity == "1", "mild",
                            ifelse(severity == "2", "moderate",
                                   ifelse(severity == "3", "severe", NA))),
  days_since_symptom_onset = NA,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = NA,
  days_since_hospitalization_estimated = NA
)])


# ----- trump_2020 -----
# NOTE: this dataset includes a re-analysis of 17 samples from the chua_2020 dataset
trump <- unique(mList$trump_2020[, .(
  sample,
  age,
  covid,
  severity,
  sex,
  dps
)])

trump <- unique(trump[, .(
  dataset = "trump_2020",
  sample_source_field = "sample",
  patient = sample,
  sample = sample,
  tissue = "Nasopharynx",
  sample_type = "Nasopharyngeal Swab", # nasopharyngeal
  sample_type_note = NA,
  race = NA,
  sex = ifelse(sex == "m", "male", "female"),
  age = age,
  covid19 = covid == "pos",
  disease_status = ifelse(covid == "pos", "COVID-19", "healthy"),
  disease_severity = ifelse(grepl("control", severity), "healthy", "severe"),
  days_since_symptom_onset = dps,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)])

# ----- yao_2021 -----
# Additional demographic and clinical info in supplemental table
# but no way to map to patient ids.
# https://www.cell.com/cms/10.1016/j.celrep.2020.108590/attachment/74f8aff2-daef-4a76-b1cd-cd7ca22f78dc/mmc2.xlsx
yao <- unique(mList$yao_2021[, .(
  group,
  orig.ident,
  patient
)])

# Get median days since onset and hospitalization
# for each severity group
meta_yao <- data.table(read_xlsx(file.path(supplementalDir, "yao_2021_Table_S1.xlsx"),
                                 n_max = 17))
meta_yao[, Group := c(rep("Moderate", 5), rep("Severe", 6), rep("Recovered", 6))]
yao_days <- meta_yao[, .(days_since_symptom_onset = round(median(`Symptoms           (Days prior to admit)` + `Sample Collection (Days after admission)`)),
                         days_since_hospitalization = median(`Sample Collection (Days after admission)`)),
                     Group]
yao <- merge(yao, yao_days, by.x = "group", by.y = "Group")

yao <- yao[, .(
  dataset = "yao_2021",
  sample_source_field = "patient",
  patient = patient,
  sample = patient,
  tissue = "Blood",
  sample_type = "Peripheral Blood Mononuclear Cells",
  sample_type_note = NA,
  race = NA,
  sex = NA,
  age = NA,
  covid19 = TRUE,
  disease_status = "COVID-19",
  disease_severity = ifelse(group == "Moderate", "moderate", "severe"),
  days_since_symptom_onset = days_since_symptom_onset,
  days_since_hospitalization = days_since_hospitalization,
  days_since_symptom_onset_estimated = TRUE,
  days_since_hospitalization_estimated = TRUE
)]

# ----- combes_2021 ----
meta_combes_path <- file.path(supplementalDir, "combes_2021_supplemental.tsv")
if ( !file.exists(meta_combes_path) ) {
  # Pull in new version of metadata
  meta_combes <- fread("https://raw.githubusercontent.com/UCSF-DSCOLAB/combes_et_al_COVID_2020/master/auxiliary_files/COMET_10X_CLINICAL_SCORES_PAPER.csv")
  fwrite(meta_combes, meta_combes_path)
} else {
  meta_combes <- fread(meta_combes_path)
}

meta_combes[, sample_id := {
  si = gsub("\\.", "-", SAMPLE.by.SNPs)
  si = gsub("-H", "-HS", si)
  si = gsub("^M1", "MVIR1", si)
  si = gsub("^X1", "XHLT1", si)
  si = gsub("-S1", "-SCG1", si)
  si = gsub("B", "BLD", si)
  si
}]

combes <- unique(mList$combes_2021[, .(
  covid_status,
  phenotype,
  pooled,
  sample,
  tissue,
  title
)])

combes <- merge(combes, meta_combes,
                by.x = "sample",
                by.y = "sample_id",
                all.x = TRUE,
                all.y = FALSE)

# impute days_since_onset
combes[covid_status.y == "POS", days_since_symptom_onset_imputed := median(Day_after_onset, na.rm = TRUE), Qualitative_score]

combes <- combes[, .(
  dataset = "combes_2021",
  sample_source_field = "sample",
  patient = sample,
  sample = sample,
  tissue = "Blood",
  sample_type = "Whole Blood",
  sample_type_note = NA,
  race = NA,
  sex = tolower(Gender),
  age = trunc(Age),
  covid19 = covid_status.y == "POS",
  disease_status = ifelse(covid_status.y == "POS", "COVID-19",
                          ifelse(covid_status.y == "NEG", "hospitalized non-COVID", "healthy")),
  disease_severity = ifelse(Qualitative_score == "SEVERE", "severe",
                            ifelse(Qualitative_score == "MILD", "moderate", "healthy")),
  days_since_symptom_onset = ifelse(is.na(Day_after_onset), days_since_symptom_onset_imputed, Day_after_onset),
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = is.na(Day_after_onset),
  days_since_hospitalization_estimated = TRUE
)]

# ----- bost_2021 -----
bost <- unique(mList$bost_2021[, .(
  age,
  clinic_status,
  sex,
  sofa_score,
  subject_id,
  sample,
  tissue
)])


bost <- unique(bost[, .(
  dataset = "bost_2021",
  sample_source_field = "sample",
  patient = subject_id,
  sample = sample,
  tissue = ifelse(tissue == "BAL", "Bronchus", "Blood"),
  sample_type = ifelse(tissue == "BAL",
                       "Bronchoalveolar Lavage Fluid",
                       "Whole Blood"),
  sample_type_note = NA,
  race = NA,
  sex = tolower(sex),
  age = age,
  covid19 = grepl("COVID", clinic_status),
  disease_status = ifelse(grepl("COVID", clinic_status),
                          "COVID-19",
                          "healthy"),
  disease_severity = sapply(clinic_status, switch,
                            "Severe COVID" = "severe",
                            "Mild COVID" = "moderate",
                            "Healthy control" = "healthy"),
  days_since_symptom_onset = sapply(clinic_status, switch,
                                    "Severe COVID" = 6,
                                    "Mild COVID" = 5,
                                    "Healthy control" = NA),
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = TRUE,
  days_since_hospitalization_estimated = TRUE
)])

# ----- merge all metadata -----

all_meta <- rbind(
  wilk,
  chua,
  silvin,
  wen,
  lee,
  liao,
  zhu,
  arunachalam,
  yu,
  su,
  meckiff,
  kusnadi,
  bacher,
  schulte_schrepping,
  grant,
  stephenson,
  heming,
  trump,
  yao,
  combes,
  bost
)
all_meta[, `:=`(days_since_symptom_onset = as.numeric(days_since_symptom_onset),
                days_since_hospitalization = as.numeric(days_since_hospitalization))]
all_meta[!disease_status %in% c("COVID-19", "healthy"), disease_status := "other"]

# estimate days_since_symptom_onset and days_since_hospitalization where missing
all_meta[days_since_symptom_onset_estimated == TRUE, days_since_symptom_onset := ifelse(is.na(days_since_symptom_onset), (days_since_hospitalization + 7), days_since_symptom_onset)]
all_meta[days_since_hospitalization_estimated == TRUE, days_since_hospitalization := ifelse(is.na(days_since_hospitalization), (days_since_symptom_onset - 7), days_since_hospitalization)]

# remove estimated days_since_hospitalization for non-hospitalized patients
all_meta[disease_severity %in% c("mild", "healthy"), `:=`(days_since_hospitalization = NA,
                                                          days_since_hospitalization_estimated = NA)]

# floor "days_since_hospitalization" at 0
all_meta[days_since_hospitalization < 0, days_since_hospitalization := 0]

all_meta[disease_status != "COVID-19", `:=`(
  days_since_symptom_onset = NA,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = NA,
  days_since_hospitalization_estimated = NA,
  disease_severity = NA)]

setnames(all_meta, "race", "race_reported")

message(">>> Writing outuput to ", file.path(metadataDir, "all_meta.tsv"))
fwrite(all_meta, file.path(metadataDir, "all_meta.tsv"), sep = "\t")


