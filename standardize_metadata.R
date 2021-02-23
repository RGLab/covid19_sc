# Standardize metadata. Reads in metadata provided with the dataset
# and standardizes fields for all datasets, writing the result to
# a file. Looks for metadata written to a tsv file in metadataDir.
# These tsv files can be written using write_metadata.R.

library(data.table)
setDTthreads(4)
library(readxl)

metadataDir <- file.path(Sys.getenv("DATA_DIR"), "metadata")
files <- list.files(metadataDir, pattern = ".tsv")
files <- files[!grepl("all_meta", files)]
mList <- lapply(file.path(metadataDir, files), fread)
names(mList) <- gsub(".tsv", "", files)
names <- lapply(mList, names)
names(names) <- files
files




# dataset, sample_source, patient, sample, tissue, sex (male, female, other), age, disease_status (COVID-19, healthy, flu), disease_severity (mild, severe)

# ----- wilk_2020 -----
# severe = hospitalized
wilk <- unique(mList$wilk_2020[, .(Donor, Donor.full, DPS, Sex, Status, Admission.level)])[, .(
  dataset = "wilk_2020",
  sample_source_field = "Donor.full",
  patient = Donor,
  sample = Donor.full,
  tissue = "PBMC",
  sex = ifelse(Sex == "M", "male", "female"),
  age = NA,
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
  tissue = location,
  sex = sex,
  age = age,
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
    tissue = "RBC",
    sex = ifelse(Characteristics.sex. == "not available", NA, Characteristics.sex.),
    age = ifelse(Characteristics.age. == "not available", NA, as.numeric(Characteristics.age.)),
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
  tissue = "PBMC",
  sex = sex,
  age = NA,
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
meta_lee <- data.table(read_xlsx(file.path(metadataDir, "supplemental", "lee_2020_Table_S1.xlsx")))
for (i in 1:nrow(meta_lee)) {
  if (is.na(meta_lee[i, "Patient ID"])) {
    meta_lee[i, `:=`("Patient ID" = meta_lee[i-1, "Patient ID"],
                     Age = meta_lee[i-1, "Age"],
                     Sex = meta_lee[i-1, "Sex"])
             ]
  }

}
meta_lee[, SampleName := paste0(`Sample ID`, " scRNA-seq")]

meta_lee2 <- data.table(read_xlsx(file.path(metadataDir, "supplemental", "lee_2020_Table_S3.xlsx")))
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
  tissue = "PBMC",
  sex = ifelse(Sex == "M", "male", "female"),
  age = Age,
  disease_status = ifelse(grepl("influenza", `Disease group`), "flu",
                          ifelse(grepl("COVID", `Disease group`), "COVID-19", "healthy")),
  disease_severity = ifelse(grepl("COVID", `Disease group`),
                            ifelse(grepl("severe", `Disease group`), "severe", "moderate"),
                            ifelse(grepl("influenza", `Disease group`), "flu", "healthy")),
  days_since_symptom_onset = NA,
  days_since_hospitalization = `Hospital \nday`,
  days_since_symptom_onset_estimated = TRUE,
  days_since_hospitalization_estimated = FALSE

)]

# ----- liao -----
# Disease severity: moderate vs severe
# NOTE: sex and age missing for controls
if ( !file.exists(file.path(metadataDir, "supplemental", "liao_2020_meta.txt")) ) {
  # Pull in new version of metadata
  meta_liao <- fread("https://raw.githubusercontent.com/zhangzlab/covid_balf/master/meta.txt")
  fwrite(meta_liao, file.path(metadataDir, "supplemental", "liao_2020_meta.txt"))
}
meta_liao <- fread(file.path(metadataDir, "supplemental", "liao_2020_meta.txt"))

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
  tissue = "BL", # bronchial lavage fluid
  sex = sex,
  age = age,
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

meta_zhu <- data.table(read_xlsx(file.path(metadataDir, "supplemental", "zhu_2020_mmc2.xlsx"), skip = 2, col_names = FALSE))
names <- unlist(meta_zhu[, 1])
meta_zhu <- data.table(t(meta_zhu[, 2:ncol(meta_zhu)]))
colnames(meta_zhu) <- names
setnames(meta_zhu, "V1", "subject")

day1_offset <- data.table(
  subject = c("COV-1", "COV-2", "COV-3", "COV-4", "COV-5", "IAV-1", "IAV-2"),
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
  tissue = "PBMC",
  sex = tolower(Gender),
  age = Age,
  disease_status = ifelse(grepl("COV", subject), "COVID-19",
                          ifelse(grepl("Flu", subject), "flu", "healthy")),
  disease_severity = ifelse(grepl("COV", subject), "moderate",
                            ifelse(grepl("Flu", subject), "flu", "healthy")),
  days_since_symptom_onset = symptom_offset + timepoint - 1,
  days_since_hospitalization = hospitalization_offset + timepoint -1,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = FALSE
)]




# ----- arunachalam -----
meta_arunachalam <- fread(file.path(metadataDir, "supplemental", "arunachalam_table_2.tsv"))

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
  tissue = "PBMC",
  sex = ifelse(sex == "M", "male", "female"),
  age = age,
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
  tissue = "PBMC",
  sex = sex,
  age = NA,
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
meta_su_cov <- data.table(read_xlsx(file.path(metadataDir, "supplemental", "su_2020_Table_S1.xlsx"), sheet = 2))

meta_su_cov <- meta_su_cov[, .(
  sample = gsub("INCOV0*", "", `Sample ID`),
  patient = gsub("INCOV0*", "", `Study Subject ID`),
  age = Age,
  sex = tolower(Sex),
  location = `Patient Location`,
  WHO = `Who Ordinal Scale`,
  disease_status = "COVID-19",
  timepoint = `Blood draw time point`
)]
meta_su_healthy <- data.table(read_xlsx(file.path(metadataDir, "supplemental", "su_2020_Table_S1.xlsx"), sheet = 3))
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
  tissue = "PBMC",
  sex = sex,
  age = age,
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
meta_meckiff <- data.table(read_xlsx(file.path(metadataDir, "supplemental", "meckiff_2020_1-s2.0-S0092867420313076-mmc1.xlsx"),
                                     sheet = 2,
                                     skip = 3))

meckiff <- unique(mList$meckiff_2020[, .(orig.donor, orig.sex, orig.hospital)])
meckiff <- merge(meckiff, meta_meckiff, by.x = "orig.donor", by.y = "PATIENT ID", all.x = TRUE)
meckiff <- meckiff[, .(
  dataset = "meckiff_2020",
  sample_source_field = "orig.donor",
  patient = orig.donor,
  sample = orig.donor,
  tissue = "PBMC",
  sex = tolower(orig.sex),
  age = as.numeric(AGE),
  disease_status = ifelse(is.na(HOSPITALIZATION), "healthy", "COVID-19"),
  disease_severity = ifelse(is.na(HOSPITALIZATION), "healthy",
                            ifelse(grepl("ICU", HOSPITALIZATION), "severe",
                                   ifelse(grepl("Ward", HOSPITALIZATION), "moderate", "mild"))),
  days_since_symptom_onset = `INTERVAL BETWEEN SYMPTOM ONSET TO SAMPLE COLLECTION (DAYS)`,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)]

# ----- schulte-schrepping -----
schulte_schrepping <- unique(mList$`schulte-schrepping_2020`[, .(donor,
                                                                 sampleID,
                                                                 sex,
                                                                 disease_stage,
                                                                 group_per_sample,
                                                                 days_after_onset,
                                                                 age,
                                                                 who_per_sample)])
# NOTE: Using min of age range for age
schulte_schrepping <- schulte_schrepping[, .(
  dataset = "schulte-schrepping_2020",
  sample_source_field = "sampleID",
  patient = donor,
  sample = sampleID,
  tissue = "PBMC",
  sex = ifelse(sex == "n/a", NA, sex),
  age = as.numeric(sub("_\\d+", "", age)),
  disease_status = ifelse(who_per_sample > 0, "COVID-19", "healthy"),
  disease_severity = ifelse(who_per_sample == 0, "healthy",
                            ifelse(who_per_sample < 3, "mild",
                                   ifelse(who_per_sample < 5, "moderate",
                                          ifelse(who_per_sample >= 5, "severe", NA)))),
  days_since_symptom_onset = days_after_onset,
  days_since_hospitalization = NA,
  days_since_symptom_onset_estimated = FALSE,
  days_since_hospitalization_estimated = TRUE
)]


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
  schulte_schrepping
)
all_meta <- all_meta[disease_status == "healthy", disease_severity := "healthy"]
all_meta[, `:=`(days_since_symptom_onset = as.numeric(days_since_symptom_onset),
                days_since_hospitalization = as.numeric(days_since_hospitalization))]

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
  days_since_hospitalization_estimated = NA)]

fwrite(all_meta, file.path(metadataDir, "all_meta.tsv"), sep = "\t")
