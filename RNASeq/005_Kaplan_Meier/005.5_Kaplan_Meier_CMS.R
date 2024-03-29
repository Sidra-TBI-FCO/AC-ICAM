
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr", "survminer")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 20                                        # SET cut-off
Group.of.interest = "CMS"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
exclude_adjuvant = "" #adjuvant_treated_excluded" #"adjuvant_treated_excluded" or ""
exclude_stage_IV = "" 
ICR = "ICR High"
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Create folders and log file
dir.create("./Figures/Trimmed_p",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Kaplan Meier Plots", showWarnings = FALSE)
dir.create("./Logfiles/Kaplan Meier Plots", showWarnings = FALSE)
dir.create("./Analysis/Trimmed_p", showWarnings = FALSE)
dir.create("./Analysis/Trimmed_p/Survival Analysis", showWarnings = FALSE)

# Read in the clinical data file
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Add ICR as a variable and assign ICR cluster according to table cluster assignment
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")

Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]
Merged_dataset[, Group.of.interest] = factor(Merged_dataset[, Group.of.interest], levels = c("CMS4", "CMS3", "CMS2", "CMS1"))

Merged_dataset = Merged_dataset[-which(is.na(Merged_dataset$CMS)),]

# Exclude adjuvant treated
if(exclude_adjuvant == "adjuvant_treated_excluded"){
  Merged_dataset = Merged_dataset[which(Merged_dataset$Adjuvant_treatment == "No"),]
}

# Exclude Stage IV
if(exclude_stage_IV == "stage_IV_excluded" ){
  Merged_dataset = Merged_dataset[-which(Merged_dataset$ajcc_pathologic_tumor_stage == "IV"),]
}

if(exclude == ""){
}else{
  Merged_dataset = Merged_dataset[-which(Merged_dataset$Patient_ID %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]),]
}

if(ICR > 1){
  Merged_dataset = Merged_dataset[which(Merged_dataset$ICR_HML == ICR),]
}




# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$vital_status == "Alive", c("vital_status", "last_contact_days_to", Group.of.interest, "ICRscore", "ICR_HML")]
colnames(TS.Alive) = c("Status","Time", Group.of.interest, "ICRscore")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$vital_status == "Dead", c("vital_status", "death_days_to", Group.of.interest, "ICRscore", "ICR_HML")]
colnames(TS.Dead) = c("Status","Time", Group.of.interest, "ICRscore")
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
mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("CMS1", "CMS2", "CMS3", "CMS4"))
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

summary(mHR)
logrank_test(Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv[which(TS.Surv[, Group.of.interest] %in% c("CMS1", "CMS2", "CMS3", "CMS4")),])

# plots
png(paste0("./Figures/Trimmed_p/Kaplan Meier Plots/Sept_2020_", ICR, "_",
           str_c(exclude, collapse = "_"), "_", exclude_stage_IV, "_", exclude_adjuvant,"_Overall_Survival_", Group.of.interest,"_Surv_cutoff_years_",Surv.cutoff.years,".png"),
    res=600,height=6,width=8,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,Group.of.interest]),
     ystrataname = "Legend",
     main= paste0("Survival curve across ", Group.of.interest),
     xlabs = "Time in months",
     palette = c("#009E74", "#E97AA8", "#0074AF", "#FF9F21"),
     #palette = c("#FF9F21", "#0074AF", "#E97AA8", "#009E74"),
     PLOT_P = signif(p[3],3),
     PLOT_HR = round(signif(exp(mHR.extract@coef),3)[3], 3),
     PLOT_CI1 = CI[3,1],
     PLOT_CI2 = CI[3,2],
     ylabs = "OS probability")
dev.off()

TS.Surv[,Group.of.interest] = factor(TS.Surv[,Group.of.interest], levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
mfit = survfit(msurv~TS.Surv[,Group.of.interest],conf.type = "log-log")

ggsurv_plot = ggsurvplot(mfit,
                         xlim = c(0, 90),
                         break.time.by = 30,
                         data = TS.Surv,
                         censor = TRUE,
                         risk.table = TRUE,
                         tables.y.text.col = TRUE,
                         tables.y.text = FALSE,
                         tables.height = 0.3,
                         tables.theme = theme_cleantable(),
                         #tables.col = "strata",
                         risk.table.pos = "out",
                         legend = "none",
                         ylab = "",
                         xlab = "Time in months",
                         fontsize = 4.5,
                         font.x = 18,
                         font.tickslab = 18,
                         censor.shape = 3,
                         censor.size = 1.5,
                         #pval = TRUE,
                         palette = c("#FF9F21", "#0074AF", "#E97AA8", "#009E74")
)

dir.create("./Figures/Trimmed_p/005_ggsurv_plots", showWarnings = FALSE)
png(paste0("./Figures/Trimmed_p/005_ggsurv_plots/Dec_2020_", ICR, "_",
           str_c(exclude, collapse = "_"), "_", exclude_stage_IV, "_", exclude_adjuvant,"_Overall_Survival_", Group.of.interest,"_Surv_cutoff_years_",Surv.cutoff.years,".png"),
    res=600,height=3.8,width=4.2,unit="in")  # set filename
print(ggsurv_plot)
dev.off()


dev.new()
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,Group.of.interest]),
     #ystrataname = "Legend",
     main= paste0("Survival curve across ", Group.of.interest),
     xlabs = "Time in months",
     ylabs = "OS probability",
     palette = c("#FF9F21", "#0074AF", "#E97AA8", "#009E74"),
     PLOT_P = signif(p[2],3),
     PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3),
     PLOT_CI1 = CI[2,1],
     PLOT_CI2 = CI[2,2],
     legend = FALSE)

# Cox regression (continuous)
uni_variate_ICRscore = coxph(formula = Surv(Time, Status) ~ ICRscore, data = TS.Surv)
summary(uni_variate_ICRscore)

# Multivariate
multivariate_ICR_and_stage = coxph(formula = Surv(Time, Status) ~ ICRscore + pathologic_stage, data = TS.Surv)
summary(multivariate_ICR_and_stage)

