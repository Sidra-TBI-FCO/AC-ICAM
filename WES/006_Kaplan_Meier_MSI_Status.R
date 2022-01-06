
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
Group.of.interest = "MSI"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
exclude_adjuvant = "" #"adjuvant_treated_excluded" or ""
exclude_stage_IV = ""
CMS = ""
ICR = ""

# Create folders and log file
dir.create("./Figures/",showWarnings = FALSE)
dir.create("./Figures/WES/006_MSI_KM", showWarnings = FALSE)
dir.create("./Analysis/WES/006_Survival Analysis", showWarnings = FALSE)

# Read in the clinical data file
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

# Divide in groups
#frequency_df$Mutational_load_category = NA
#frequency_df$Mutational_load_category[which(frequency_df$Nonsilent_mutational_burden_per_Mb < 12)] = "Low"
#frequency_df$Mutational_load_category[which(frequency_df$Nonsilent_mutational_burden_per_Mb >= 12)] = "High"

clinical_data = clinical_data[which(clinical_data$Patient_ID %in% frequency_df$Patient_ID),]
Merged_dataset = clinical_data
Merged_dataset$MSI = MANTIS$MSI[match(Merged_dataset$Patient_ID, MANTIS$Patient_ID)]
  
Merged_dataset$MSI = factor(Merged_dataset$MSI, levels = c("MSI-H", "MSS"))
Merged_dataset$ICR_cluster = table_cluster_assignment$ICR_HML[match(Merged_dataset$Patient_ID,
                                                                    substring(rownames(table_cluster_assignment), 1, 3))]

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

if(ICR == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ICR_cluster == ICR),]
}

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$OS.Status == "Alive", c("OS.Status", "OS.Time", Group.of.interest)]
colnames(TS.Alive) = c("Status","Time", Group.of.interest)
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$OS.Status == "Dead", c("OS.Status", "OS.Time", Group.of.interest)]
colnames(TS.Dead) = c("Status","Time", Group.of.interest)
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y


TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

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
mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("MSI-H", "MSS"))
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


# plots
png(paste0("./Figures/WES/006_MSI_KM/Overall_Survival_", Group.of.interest, "_", exclude_adjuvant, "_",
           exclude_stage_IV, ICR, "_Surv_cutoff_years_",Surv.cutoff.years,".png"),res=600,height=6,width=8,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,Group.of.interest]),
     ystrataname = "Legend",
     main= paste0("Survival curve across ", Group.of.interest),
     xlabs = "Time in months",
     palette = c("#FF5800", "#8085E9"),
     ylabs = "OS probability",
     PLOT_P = signif(p[1],3),
     PLOT_HR = round(signif(exp(mHR.extract@coef),3)[1], 3),
     PLOT_CI1 = CI[1,1],
     PLOT_CI2 = CI[1,2])
dev.off()

# Cox regression (continuous)
TS.Surv$Mutational_load_per_Mb = log(TS.Surv$Mutational_load_per_Mb + 1, 10)
uni_variate_ICRscore = coxph(formula = Surv(Time, Status) ~ Mutational_load_per_Mb, data = TS.Surv)
summary(uni_variate_ICRscore)

# Multivariate
multivariate_ICR_and_stage = coxph(formula = Surv(Time, Status) ~ ICRscore + pathologic_stage, data = TS.Surv)
summary(multivariate_ICR_and_stage)


