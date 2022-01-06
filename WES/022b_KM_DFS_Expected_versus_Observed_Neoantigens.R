

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(required.packages = c("survival", "plyr"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

# Load data
load("./Analysis/WES/020_Mut_Load_and_Neoantigen/Nonsilent_mutation_frequency_and_filtered_Neoantigen_count.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Set parameters
subset = ""  #"nonhypermutated" "hypermutated"
Group.of.interest = "Immunoedited"  #"Immunoedited_extremes"
x = 1
Surv.cutoff.years = 10
Stages = "all stages" # "all stages" # c(3, 4) # "all stages"  or c(1, 2) or c(3, 4)
Stage_name = "All" #"Stage I&II" #"All"
Treatment = "" # "No_treatment" or "Adjuvant_treated"

# Hypermutation status
frequency_df$Mutation_cat = NA
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"

# Ratio calculation
frequency_df$Ratio = frequency_df$Neoantigen_count / frequency_df$Non_silent_Mutation_frequency

fit = lm(Neoantigen_count ~ Non_silent_Mutation_frequency, data = frequency_df)
summary(fit)
print(fit)
# Formula: Neoantigen_count = -2.38770 + (0.09171  * Non_silent_Mutation_frequency)
# For nonhypermutated formula: Neoantigen_count = -0.61144 + (0.08724  * Non_silent_Mutation_frequency)
# For hypermutated formula: Neoantigen_count = -0.61144 + (0.08724  * Non_silent_Mutation_frequency)
if(subset == ""){
  frequency_df$Expected_neoantigen_count = -2.38770 + (0.09171 * frequency_df$Non_silent_Mutation_frequency)}
if(subset == "nonhypermutated"){
  frequency_df$Expected_neoantigen_count = -0.61144 + (0.08724 * frequency_df$Non_silent_Mutation_frequency)}
if(subset == "hypermutated"){
  frequency_df$Expected_neoantigen_count = -11.88801 + (0.09503 * frequency_df$Non_silent_Mutation_frequency)}

frequency_df$Immunoediting_score = (frequency_df$Neoantigen_count / frequency_df$Expected_neoantigen_count) #/ frequency_df$Non_silent_Mutation_frequency

frequency_df$Immunoedited = NA
frequency_df$Immunoedited[which(frequency_df$Immunoediting_score < 1)] = "immunoedited"
frequency_df$Immunoedited[which(frequency_df$Immunoediting_score >= 1)] = "less immunoedited"
table(frequency_df$Immunoedited)

frequency_df$Immunoedited_extremes = NA
frequency_df$Immunoedited_extremes[which(frequency_df$Immunoediting_score < 0.9)] = "immunoedited"
frequency_df$Immunoedited_extremes[which(frequency_df$Immunoediting_score >= 1.1)] = "less immunoedited"
frequency_df$Immunoedited_extremes[which(is.na(frequency_df$Immunoedited_extremes))] = "medium immunoediting"
table(frequency_df$Immunoedited_extremes)

frequency_df$ICRscore = table_cluster_assignment$ICRscore[match(frequency_df$Patient_ID,
                                                                substring(rownames(table_cluster_assignment), 1, 3))]

MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
frequency_df$MSI_status = MANTIS$MSI[match(frequency_df$Patient_ID,
                                           MANTIS$Patient_ID)]

Merged_dataset = clinical_data[which(clinical_data$Patient_ID %in% frequency_df$Patient_ID),]
Merged_dataset = merge(Merged_dataset, frequency_df, by = "Patient_ID")

#Merged_dataset = Merged_dataset[-which(Merged_dataset$ajcc_pathologic_tumor_stage == "4"),]

table(Merged_dataset$ajcc_pathologic_tumor_stage)
if(Stages == "all stages"){
}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ajcc_pathologic_tumor_stage %in% Stages),]
}
table(Merged_dataset$ajcc_pathologic_tumor_stage)

Merged_dataset$ajcc_pathologic_tumor_stage = as.numeric(Merged_dataset$ajcc_pathologic_tumor_stage)
Merged_dataset$age_at_initial_pathologic_diagnosis = as.numeric(Merged_dataset$age_at_initial_pathologic_diagnosis)
Merged_dataset$MSI_status = factor(Merged_dataset$MSI_status, levels = c("MSS", "MSI-H"))

table(Merged_dataset$Adjuvant_treatment, Merged_dataset$ajcc_pathologic_tumor_stage)
if(Treatment == "No_treatment"){
  Merged_dataset = Merged_dataset[which(Merged_dataset$Adjuvant_treatment == "No"),]
}
if(Treatment == "Adjuvant_treated"){
  Merged_dataset = Merged_dataset[-which(Merged_dataset$Adjuvant_treatment == "No"),]
}



### KM
# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$DFS.Status == "Disease Free", c("DFS.Status", "DFS.Time", Group.of.interest,
                                                                 "ICRscore", "ajcc_pathologic_tumor_stage", 
                                                                 "age_at_initial_pathologic_diagnosis", "MSI_status")]
colnames(TS.Alive) = c("Status","Time", Group.of.interest, "ICRscore", "Stage", "Age", "MSI")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$DFS.Status == "Event", c("DFS.Status", "DFS.Time", Group.of.interest,
                                                               "ICRscore", "ajcc_pathologic_tumor_stage",
                                                               "age_at_initial_pathologic_diagnosis", "MSI_status")]
colnames(TS.Dead) = c("Status","Time", Group.of.interest, "ICRscore", "Stage", "Age", "MSI")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Disease Free"
TS.Dead$Time[TS.Dead$Time > Y] = Y


TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Event"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

if(Group.of.interest == "Immunoedited_extremes"){
  TS.Surv[,Group.of.interest] = factor(TS.Surv[,Group.of.interest], levels = c("immunoedited", "medium immunoediting", "less immunoedited"))
}
# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv[,Group.of.interest],conf.type = "log-log")

# Calculations (Needs manual adaptation!)
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

Cal.Surv = TS.Surv
#Cal.Surv[,Group.of.interest] = as.factor(Cal.Surv[,Group.of.interest])
#Cal.Surv[,Group.of.interest] = relevel(Cal.Surv[,Group.of.interest], "ICR3")
Cal.Surv[,Group.of.interest] = as.factor(Cal.Surv[,Group.of.interest])
if(Group.of.interest == "Immunoedited_extremes"){
  Cal.Surv[,Group.of.interest] = factor(Cal.Surv[,Group.of.interest], levels = c("immunoedited", "medium immunoediting", "less immunoedited"))
}

mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("immunoedited", "less immunoedited"))
mHR.extract = extract.coxph(mHR, include.aic = TRUE,
                            include.rsquared = TRUE, include.maxrs=TRUE,
                            include.events = TRUE, include.nobs = TRUE,
                            include.missings = TRUE, include.zph = TRUE)
HRtxt = paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta = coef(mHR)
se   = sqrt(diag(mHR$var))
p    = 1 - pchisq((beta/se)^2, 1)
CI   = confint(mHR)
CI   = round(exp(CI),2)

if(x == 1){
  palette = c("#FF5800", "#8085E9")
}
if(x == 2){
  palette = c("#FF5800", "#FFCC00","#8085E9")
}
png(paste0("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/ggkm_DFS_", subset, "_by_immunoediting_status_", 
           Group.of.interest, "_", Stage_name, "_", Treatment,".png"),
    res=600,height=6,width=10,unit="in")
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,Group.of.interest]),
     ystrataname = "Legend",
     main= paste0("Survival curve across ", Group.of.interest),
     legend = FALSE,
     xlabs = "Time in months",
     palette = palette,
     ylabs = "DFS probability",
     PLOT_P = signif(p[x],3),
     PLOT_HR = round(signif(exp(mHR.extract@coef),3)[x], 3),
     PLOT_CI1 = CI[x,1],
     PLOT_CI2 = CI[x,2])
dev.off()

dir.create("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens", showWarnings = FALSE)
#save(frequency_df, Merged_dataset, file = "./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/frequency_df.Rdata")

TS.Surv$Immunoedited = factor(TS.Surv$Immunoedited, levels = c("less immunoedited", "immunoedited"))

#multivariate = coxph(formula = Surv(Time, Status) ~ ICRscore + MSI + Stage + Age + Immunoedited, data = TS.Surv)
univariate_immunoedited = coxph(formula = Surv(Time, Status) ~ Immunoedited, data = TS.Surv)
summary(univariate_immunoedited)

univariate_ICRscore = coxph(formula = Surv(Time, Status) ~ ICRscore, data = TS.Surv)
summary(univariate_ICRscore)

univariate_stage = coxph(formula = Surv(Time, Status) ~ Stage, data = TS.Surv)
summary(univariate_stage)

univariate_age = coxph(formula = Surv(Time, Status) ~ Age, data = TS.Surv)
summary(univariate_age)

multivariate = coxph(formula = Surv(Time, Status) ~ ICRscore + Immunoedited + Stage + Age , data = TS.Surv)
summary(multivariate)
