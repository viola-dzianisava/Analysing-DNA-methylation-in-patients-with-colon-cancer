library(TCGAbiolinks)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)

# Define download directory
download_dir <- "/Users/linxijun/Desktop/GDCdata"
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)
}


# Query clinical data
clin_coad <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical", save.csv = FALSE)

# Check structure
str(clin_coad)

head(clin_coad)

# Deal with days_to_death and vital_status column
# Switch vital_status into numberic（0 = alive, 1 = deceased）
clin_coad$status <- as.numeric(clin_coad$vital_status == "Dead")

clin_coad$days_to_death[is.na(clin_coad$days_to_death)] <- max(clin_coad$days_to_death, na.rm = TRUE)

library(survival)
library(survminer)

# Create Surv objective
surv_obj <- Surv(time = clin_coad$days_to_death, event = clin_coad$status)

# Using Cox model
cox_model <- coxph(surv_obj ~ 1, data = clin_coad)  
summary(cox_model)

# Create Kaplan-Meier survival cruve
km_fit <- survfit(surv_obj ~ 1, data = clin_coad)
ggsurvplot(km_fit, data = clin_coad, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Estimate for TCGA-COAD Patients",
           xlab = "Days to Death", ylab = "Survival Probability")



##======stage
library(TCGAbiolinks)
library(survival)
library(survminer)
library(RColorBrewer)
library(dplyr)

# load data
clin_coad <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical", save.csv = FALSE)


# Simplified into 4 stage of cancer
clin_coad <- clin_coad %>%
  mutate(ajcc_pathologic_stage_simple = case_when(
    grepl("Stage IV", ajcc_pathologic_stage) ~ "Stage IV",
    grepl("Stage III", ajcc_pathologic_stage) ~ "Stage III",
    grepl("Stage II", ajcc_pathologic_stage) ~ "Stage II",
    grepl("Stage I", ajcc_pathologic_stage) ~ "Stage I",
    TRUE ~ NA_character_
  ))

# Check
print(table(clin_coad$ajcc_pathologic_stage_simple))

# Clean NA values in ajcc_pathologic_stage_simple column
clin_coad <- clin_coad %>%
  filter(!is.na(ajcc_pathologic_stage_simple))

# Check vital_status and days_to_death columns are exist and correct
clin_coad$status <- as.numeric(clin_coad$vital_status == "Dead")
clin_coad$days_to_death[is.na(clin_coad$days_to_death)] <- max(clin_coad$days_to_death, na.rm = TRUE)

# Check
print(table(clin_coad$ajcc_pathologic_stage_simple))
print(table(clin_coad$status))

# Create sure obj
surv_obj <- Surv(time = clin_coad$days_to_death, event = clin_coad$status)

# 4 stages
km_fit <- survfit(surv_obj ~ ajcc_pathologic_stage_simple, data = clin_coad)

# Check the fit results
print(summary(km_fit))

# Color
num_stages <- length(unique(clin_coad$ajcc_pathologic_stage_simple))
colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(num_stages)

## Create Kaplan-Meier
ggsurvplot(km_fit, data = clin_coad, pval = TRUE, risk.table = TRUE,
           palette = colors,
           title = "Kaplan-Meier Estimate by Cancer Stage for TCGA-COAD Patients",
           xlab = "Days to Death", ylab = "Survival Probability",
           legend.labs = c("Stage I", "Stage II", "Stage III", "Stage IV"),
           legend.title = "Stage",
           risk.table.height = 0.25,   
           risk.table.y.text.col = TRUE, 
           risk.table.y.text = TRUE,  
           risk.table.y.text.angle = 0, 
           risk.table.fontsize = 3)  
#======stage and gender
library(TCGAbiolinks)
library(survival)
library(survminer)
library(RColorBrewer)
library(dplyr)


clin_coad <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical", save.csv = FALSE)


if (!("vital_status" %in% colnames(clin_coad))) {
  stop("The vital_status column is missing from the clinical data.")
}
if (!("days_to_death" %in% colnames(clin_coad))) {
  stop("The days_to_death column is missing from the clinical data.")
}
if (!("gender" %in% colnames(clin_coad))) {
  stop("The gender column is missing from the clinical data.")
}

clin_coad <- clin_coad %>%
  mutate(ajcc_pathologic_stage_simple = case_when(
    grepl("Stage IV", ajcc_pathologic_stage) ~ "Stage IV",
    grepl("Stage III", ajcc_pathologic_stage) ~ "Stage III",
    grepl("Stage II", ajcc_pathologic_stage) ~ "Stage II",
    grepl("Stage I", ajcc_pathologic_stage) ~ "Stage I",
    TRUE ~ NA_character_
  ))


clin_coad <- clin_coad %>%
  filter(!is.na(ajcc_pathologic_stage_simple))


clin_coad$status <- as.numeric(clin_coad$vital_status == "Dead")
clin_coad$days_to_death[is.na(clin_coad$days_to_death)] <- max(clin_coad$days_to_death, na.rm = TRUE)


print(table(clin_coad$ajcc_pathologic_stage_simple))
print(table(clin_coad$status))
print(table(clin_coad$gender))


surv_obj <- Surv(time = clin_coad$days_to_death, event = clin_coad$status)

# stage and gender
km_fit <- survfit(surv_obj ~ ajcc_pathologic_stage_simple + gender, data = clin_coad)


print(summary(km_fit))

num_stages_gender <- length(unique(clin_coad$ajcc_pathologic_stage_simple)) * 2

colors <- c(
  "#A1D99B", "#FCBBA1", # Stage I Male, Stage I Female
  "#74C476", "#FB6A4A", # Stage II Male, Stage II Female
  "#41AB5D", "#EF3B2C", # Stage III Male, Stage III Female
  "#238B45", "#CB181D"  # Stage IV Male, Stage IV Female
)

ggsurvplot(km_fit, data = clin_coad, pval = TRUE, risk.table = TRUE,
           palette = colors,
           title = "Kaplan-Meier Estimate by Cancer Stage and Gender for TCGA-COAD Patients",
           xlab = "Days to Death", ylab = "Survival Probability",
           legend.labs = c("Stage I - Male", "Stage I - Female", 
                           "Stage II - Male", "Stage II - Female", 
                           "Stage III - Male", "Stage III - Female", 
                           "Stage IV - Male", "Stage IV - Female"),
           legend.title = "Stage and Gender",
           risk.table.height = 0.4,   # 调整风险表格高度
           risk.table.y.text.col = TRUE,  # 确保风险表格中的文本颜色
           risk.table.y.text = TRUE,  # 确保风险表格中的文本可见
           risk.table.y.text.angle = 0,  # 文本角度调整为水平
           risk.table.fontsize = 3)  # 调整风险表格字体大小

#======gender only
library(TCGAbiolinks)
library(survival)
library(survminer)
library(RColorBrewer)
library(dplyr)


clin_coad <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical", save.csv = FALSE)


if (!("vital_status" %in% colnames(clin_coad))) {
  stop("The vital_status column is missing from the clinical data.")
}
if (!("days_to_death" %in% colnames(clin_coad))) {
  stop("The days_to_death column is missing from the clinical data.")
}
if (!("gender" %in% colnames(clin_coad))) {
  stop("The gender column is missing from the clinical data.")
}


clin_coad$status <- as.numeric(clin_coad$vital_status == "Dead")
clin_coad$days_to_death[is.na(clin_coad$days_to_death)] <- max(clin_coad$days_to_death, na.rm = TRUE)


print(table(clin_coad$status))
print(table(clin_coad$gender))


surv_obj <- Surv(time = clin_coad$days_to_death, event = clin_coad$status)

km_fit <- survfit(surv_obj ~ gender, data = clin_coad)


print(summary(km_fit))

colors <- c("#41AB5D", "#EF3B2C") # Male, Female

=
ggsurvplot(km_fit, data = clin_coad, pval = TRUE, risk.table = TRUE,
           palette = colors,
           title = "Kaplan-Meier Estimate by Gender for TCGA-COAD Patients",
           xlab = "Days to Death", ylab = "Survival Probability",
           ylim = c(0.6, 1), 
           pval.coord = c(500, 0.65),  
           legend.labs = c("Male", "Female"),
           legend.title = "Gender",
           risk.table.height = 0.2,   
           risk.table.y.text.col = TRUE, 
           risk.table.y.text = TRUE,  
           risk.table.y.text.angle = 0, 
           risk.table.fontsize = 3)  


#======race only
library(TCGAbiolinks)
library(survival)
library(survminer)
library(RColorBrewer)
library(dplyr)


clin_coad <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical", save.csv = FALSE)


if (!("vital_status" %in% colnames(clin_coad))) {
  stop("The vital_status column is missing from the clinical data.")
}
if (!("days_to_death" %in% colnames(clin_coad))) {
  stop("The days_to_death column is missing from the clinical data.")
}
if (!("race" %in% colnames(clin_coad))) {
  stop("The race column is missing from the clinical data.")
}


clin_coad <- clin_coad %>%
  filter(!is.na(race) & race != "not reported" & race != "american indian or alaska native")

clin_coad$status <- as.numeric(clin_coad$vital_status == "Dead")
clin_coad$days_to_death[is.na(clin_coad$days_to_death)] <- max(clin_coad$days_to_death, na.rm = TRUE)

print(table(clin_coad$status))
print(table(clin_coad$race))

surv_obj <- Surv(time = clin_coad$days_to_death, event = clin_coad$status)

km_fit <- survfit(surv_obj ~ race, data = clin_coad)

print(summary(km_fit))

num_races <- length(unique(clin_coad$race))
colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a") 

# label
legend_labs <- c("Asian", "Black or African American", "White")

ggsurvplot(km_fit, data = clin_coad, pval = TRUE, risk.table = TRUE,
           palette = colors,
           title = "Kaplan-Meier Estimate by Race for TCGA-COAD Patients",
           xlab = "Days to Death", ylab = "Survival Probability",
           ylim = c(0.6, 1), 
           pval.coord = c(500, 0.65), 
           legend.labs = legend_labs,
           legend.title = "Race",
           risk.table.height = 0.2, 
           risk.table.y.text.col = TRUE,  
           risk.table.y.text = TRUE,
           risk.table.y.text.angle = 0, 
           risk.table.fontsize = 3) 
