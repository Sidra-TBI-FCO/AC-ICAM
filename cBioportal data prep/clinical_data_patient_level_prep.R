

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 10                                        # SET cut-off
Group.of.interest = "ICR_HML"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
exclude_adjuvant = "" #"adjuvant_treated_excluded" or ""
exclude_stage_IV = ""
CMS = ""
ICR_cluster_with_metastasis = ""
exclude = c("Conpair_lower_90_percent", "non-epithelial")


# Read in the clinical data file
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")
# Add ICR as a variable and assign ICR cluster according to table cluster assignment

load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")

table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset[, Group.of.interest] = factor(Merged_dataset[, Group.of.interest], levels = c("ICR High", "ICR Medium", "ICR Low"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]
Merged_dataset$CMS[which(is.na(Merged_dataset$CMS))] = "mixed"

# Exclude adjuvant treated
if(exclude_adjuvant == "adjuvant_treated_excluded"){
  Merged_dataset = Merged_dataset[which(Merged_dataset$Adjuvant_treatment == "No"),]
}

# Exclude Stage IV
if(exclude_stage_IV == "exclude_stage_IV" ){
  Merged_dataset = Merged_dataset[-which(Merged_dataset$ajcc_pathologic_tumor_stage == "4"),]
}

if(CMS == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$CMS == CMS),]
}

if(exclude == ""){
}else{
  Merged_dataset = Merged_dataset[-which(Merged_dataset$Patient_ID %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]),]
}
# for cbioportal
cbio = Merged_dataset
cbio$Patient_ID = paste("SER-SILU-CC-P0", cbio$Patient_ID, sep = "")
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "adenocarcinoma nno")] = "adenocarcinoma"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "adenocarcinoom in villeus adenoom")] = "adenocarcinoma in villeus adenoma"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "adenocarcinoom met gemengde subtypes")] = "adenocarcinoma with mixed subtypes"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "adenocarcinoom, intestinaal type")] = "adenocarcinoma intestinal type"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "cribriform carcinoom")] = "cribriform carcinoma"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "grootcellig neuroendocrien carcinoom")] = "large cell neuroendocrine carcinoma"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "maligne lymfoom, groot b-cel, diffuus, nno")] = "malignant large b cell lymphoma"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "mucineus adenocarcinoom")] = "mucineus adenocarcinoma"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "neoplasma, maligne")] = "adenocarcinoma"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "neuro-endocrien carcinoom")] = "neuroendocrine carcinoma"
cbio$Tumor_morphology[which(cbio$Tumor_morphology == "zegelringcel carcinoom")] = "signet ring cell carcinoma"

cbio$cause_of_death[which(cbio$cause_of_death == "infectie")] = "infection"
cbio$cause_of_death[which(cbio$cause_of_death %in% c("onbekende oorzaak", "onbekende en niet gespecificeerde oorzaken van ziekte"))] = "unknown cause"
cbio$cause_of_death[which(cbio$cause_of_death == "maligne hersentumor")] = "malignant brain tumor"
cbio$cause_of_death[which(cbio$cause_of_death == "aandoening van zenuwstelsel/neurologische aandoening")] = "CNS-related illness"
cbio$cause_of_death[which(cbio$cause_of_death == "alvleesklierkanker")] = "pancreatic tumor"
cbio$cause_of_death[which(cbio$cause_of_death == "baarmoederkanker")] = "uterine cancer"
cbio$cause_of_death[which(cbio$cause_of_death == "blaaskanker")] = "bladder cancer"
cbio$cause_of_death[which(cbio$cause_of_death == "dikkedarmkanker")] = "colon cancer"
cbio$cause_of_death[which(cbio$cause_of_death == "eierstokkanker")] = "ovary cancer"
cbio$cause_of_death[which(cbio$cause_of_death == "endeldarmkanker")] = "rectal cancer"
cbio$cause_of_death[which(cbio$cause_of_death == "hart- en vaatziekten")] = "cardiovascular disease"
cbio$cause_of_death[which(cbio$cause_of_death == "longkanker")] = "lung cancer"
cbio$cause_of_death[which(cbio$cause_of_death == "meerdere maligniteiten")] = "multiple malignancies"
cbio$cause_of_death[which(cbio$cause_of_death == "mesothelioom")] = "mesothelioma"
cbio$cause_of_death[which(cbio$cause_of_death == "myelo\x95de leukemie")] = "myeloid leukemia"
cbio$cause_of_death[which(cbio$cause_of_death == "niet oncologisch")] = "not cancer-related"
cbio$cause_of_death[which(cbio$cause_of_death == "non-hodgkinlymfoom")] = "non-hodgkin lymphoma"
cbio$cause_of_death[which(cbio$cause_of_death == "onbekende en niet gespecificeerde oorzaken van ziekte")] = "not cancer-related"
cbio$cause_of_death[which(cbio$cause_of_death == "non-hodgkinlymfoom")] = "non-hodgkin lymphoma"
cbio$cause_of_death[which(cbio$cause_of_death == "ongevallen (door pati\x91nt zelf of door anderen)")] = "accident"
cbio$cause_of_death[which(cbio$cause_of_death == "schaamlipkanker")] = "labia cancer"
cbio$cause_of_death[which(cbio$cause_of_death == "slokdarmkanker")] = "esophagus cancer"
cbio$cause_of_death[which(cbio$cause_of_death == "ziekte van ademhalingsstelsel")] = "respiratory disease"
cbio$cause_of_death[which(cbio$cause_of_death == "ziekte van spijsverteringsstelsel")] = "disease of digestive system"
cbio$cause_of_death[which(cbio$cause_of_death == "ziekte van urogenitaal stelsel")] = "disease of genitourinary system"

table(cbio$cause_of_death)

# Dead to DECEASED and Alive to LIVING
cbio$OS.Status[which(cbio$OS.Status == "Dead")] = "DECEASED"
cbio$OS.Status[which(cbio$OS.Status == "Alive")] = "LIVING"

# DFS_STATUS expected values are DiseaseFree, Recurred/Progressed Recurred or Progressed
cbio$DFS.Status[which(cbio$DFS.Status == "Disease Free")] = "DiseaseFree"
cbio$DFS.Status[which(cbio$DFS.Status == "Event")] = "Recurred"

cbio = cbio[,c(1:26, 28:34, 38, 43:50, 53, 64:65)]

#added_patient = c("SER-SILU-CC-P0017", "FEMALE", "[Not Available]", "[Not Available]", "colon transversum",
 #                 "adenocarcinoma", "Yes",
  #                "No", "No", "3", "0", "1", "4", "[Not Applicable]", "0", "1", "4",
   #               "liver, peritoneum", "metastasectomy", "64", "233384", "2003", rep(NA, 23))
#cbio[349,] = added_patient

cbio$DFS_def2.Status = NULL
cbio$DFS_def2.Time = NULL
cbio$ethnicity = NULL
cbio$race = NULL

colnames(cbio)[which(nchar(colnames(cbio)) > 25)]
colnames(cbio)[which(colnames(cbio) == "tumour_anatomic_site")] = "tumor_anatomic_location"
colnames(cbio)[which(colnames(cbio) == "ajcc_tumor_pathologic_pt")] = "Path_tumor_stage"
colnames(cbio)[which(colnames(cbio) == "ajcc_nodes_pathologic_pn")] = "Path_nodes_stage"
colnames(cbio)[which(colnames(cbio) == "ajcc_metastasis_pathologic_pm")] = "Path_metastasis_stage"
colnames(cbio)[which(colnames(cbio) == "ajcc_pathologic_tumor_stage")] = "ajcc_path_stage"
colnames(cbio)[which(colnames(cbio) == "age_at_initial_pathologic_diagnosis")] = "age_at_dx"
colnames(cbio)[which(colnames(cbio) == "year_of_initial_diagnosis")] = "year_of_dx"
colnames(cbio)[which(colnames(cbio) == "location_metastasis_at_diagnosis")] = "location_metastasis_at_dx"
colnames(cbio)[which(colnames(cbio) == "history_other_malignancy")] = "history_of_cancer"
colnames(cbio)[which(colnames(cbio) == "therapy_metastasis_present_at_diagnosis")] = "therapy_metastasis_at_dx"
colnames(cbio)[which(colnames(cbio) == "history_neoadjuvant_treatment")] = "neoadjuvant_treated"
colnames(cbio)[which(colnames(cbio) == "new_primary_tumor_in_FU_dx_indicator")] = "new_primary_tumor_in_FU"
colnames(cbio) = gsub("\\.", "_", colnames(cbio))

# Add TCR data
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
cbio$TCR_productive_templates = TCR_Overview$productive_templates[match(substring(cbio$Patient_ID, 15, 17),
                                                                        TCR_Overview$Patient_ID)]

cbio$TCR_productive_rearrangements = TCR_Overview$productive_rearrangements[match(substring(cbio$Patient_ID, 15, 17),
                                                                             TCR_Overview$Patient_ID)]
  
cbio$TCR_productive_clonality = TCR_Overview$productive_clonality[match(substring(cbio$Patient_ID, 15, 17),
                                                                        TCR_Overview$Patient_ID)]

colnames(cbio)
colnames(cbio)[which(colnames(cbio) == "DFS_Status")] = "PFS_Status"
colnames(cbio)[which(colnames(cbio) == "DFS_Time")] = "PFS_Time"

# Add IES classification
load("./Analysis/WES/022f_IES_categories/022f_IES_df.Rdata")

cbio$IES = IES_df$IES[match(substring(cbio$Patient_ID, 15, 17), IES_df$Patient_ID)]

write.csv(cbio, file = "./Processed_Data/Shared_Data/For cbioportal/SER-SILU-CC-clinical_data_for_cbioportal_1609_2021.csv", row.names = FALSE)
write.table(cbio, file = "./Processed_Data/Shared_Data/For cbioportal/SER-SILU-CC-clinical_data_for_cbioportal_1609_2021.tsv", sep = "\t", row.names = FALSE)
