# Title: Graft and Recipient Survival by HLA Mismatch
# Author: Dave Merola
# Date Created: 12-20-2024
# Description: This study investigates the effect of HLA mismatch on graft and 
# recipient survival.

rm(list = ls())

# Constants ----
# Specify directories
DIR_RAW_DATA <- "~/OneDrive/Research/UNOS Data/STAR_SAS/SAS Dataset 202409/Kidney_ Pancreas_ Kidney-Pancreas/"
DIR_FORMATS <- "~/OneDrive/Research/UNOS Data/STAR_SAS/CODE DICTIONARY - FORMATS 202409/Kidney_ Pancreas_ Kidney-Pancreas/"
DIR_DATA <- "~/OneDrive/Research/Transplant HLA MM Kidney/Data/"
DIR_OUTPUT <- "~/OneDrive/Research/Transplant HLA MM Kidney/Output/"

# Load Libraries ----
library(survey)
library(doParallel)
library(doRNG)
library(WeightIt)
library(nnet)
library(adjustedCurves)
library(tableone)
library(mice)
library(haven)
library(tidyverse)

# Functions ----
AlleleMismatchLevel <- function(don_allele_1, don_allele_2, rec_allele_1, rec_allele_2) {
  
  # Description:
  # This function takes four vectors (two each from donor & recipient) and makes  
  # rowwise comparisons to determine whether there are 0, 2, or 4 mismatches
  
  # Arguments: 
  # don_allele_1: Vector containing each donor individual's locus for donor allele 1
  # don_allele_2: Vector containing each donor individual's locus for donor allele 2
  # rec_allele_1: Vector containing each donor individual's locus for receipient allele 1
  # rec_allele_2: Vector containing each donor individual's locus for recipient allele 2
  
  # Output:
  # A variable containing the mismatch level (i.e., 0, 2, or 4) for each patient
  
  don_allele_matrix <- rbind(t(don_allele_1),
                             t(don_allele_2))
  
  rec_allele_matrix <- rbind(t(rec_allele_1),
                             t(rec_allele_2))
  
  for (i in 1:ncol(don_allele_matrix)) {
    
    if (i == 1) {output <- rep(NA, ncol(don_allele_matrix))}
    
    don_match <- sum(don_allele_matrix[, i] %in% rec_allele_matrix[, i])
    
    rec_match <- sum(rec_allele_matrix[, i] %in% don_allele_matrix[, i])
    
    match_level <- case_when(don_match + rec_match == 0 ~ 4,
                             don_match + rec_match == 2 ~ 2,
                             don_match + rec_match == 3 ~ 2,
                             don_match + rec_match == 4 ~ 0,
                             TRUE ~ NA_real_)
    
    output[i] <- match_level
    
  }
  
  return(output)
  
}

# Load Data ----

# Convert files to lower-memory .RData format
# data_kidney <- read_sas(paste0(DIR_RAW_DATA, "kidpan_data.sas7bdat"), NULL)
# save(data_kidney, file = paste0(DIR_DATA, "kidpan_data.RData"))

# data_kidney_immuno <- read_sas(paste0(DIR_RAW_DATA, "Immunosuppression/kidpan_immuno_discharge_data.sas7bdat"), NULL)
# save(data_kidney_immuno, file = paste0(DIR_DATA, "data_kidney_immuno.RData"))

# data_kidney_pra <- read_sas(paste0(DIR_RAW_DATA, "PRA and Crossmatch/kidpan_pra_crossmatch_data.sas7bdat"), NULL)
# save(data_kidney_pra, file = paste0(DIR_DATA, "data_kidney_pra.RData"))

# data_kidney_hla <- read_sas(paste0(DIR_RAW_DATA, "HLA Additional/kidpan_addtl_hla.sas7bdat"), NULL)
# save(data_kidney_hla, file = paste0(DIR_DATA, "data_kidney_hla.RData"))


# Load .RData files (Runtime ~ 2 mins)
load(file = paste0(DIR_DATA, "kidpan_data.RData"))
load(file = paste0(DIR_DATA, "data_kidney_immuno.RData"))
load(file = paste0(DIR_DATA, "data_kidney_hla.RData"))

formats <- read.delim(paste0(DIR_FORMATS, "KIDPAN_FORMATS_FLATFILE.DAT"), 
                      header = FALSE, 
                      col.names = c("label", "fmtname", "type", "code"))

# Data Management ----

## Variable Selection ----
# Pull key variables from Immuno table
data_kidney_immuno_select <- select(data_kidney_immuno,
                                    TRR_ID_CODE, ATGAM_IND, THYMOGLOBULIN_IND, 
                                    OKT3_IND, CAMPATH_IND, SIMULECT_IND, 
                                    STEROIDS_IND, ENVARSUSXR_MAINT, PROGRAF_MAINT, 
                                    GENERICCYCLOSPORIN_MAINT, CELLCEPT_MAINT, 
                                    GENERICMYFORTIC_MAINT)

# Pull key variables from HLA table
data_kidney_hla_select <- select(data_kidney_hla, TRR_ID_CODE, DC1, DC2)

# Select key variables from KIDPAN table and assemble analytic cohort 
data_kidney_cohort <- data_kidney %>% 
  
  # Join kidney_immuno table
  left_join(data_kidney_immuno_select, by = 'TRR_ID_CODE') %>% 
  
  # Join HLA table
  left_join(data_kidney_hla_select, by = 'TRR_ID_CODE') %>% 
  
  # Select variables of interest and tag as outcome, eligibility, exposure, or covariate
  select(trr_id_code = "TRR_ID_CODE", 
         pt_code = "PT_CODE", 
         out_gstatus = "GSTATUS_KI", 
         out_gtime = "GTIME_KI", 
         out_pstatus = "PSTATUS", 
         out_ptime = "PTIME", 
         out_composite_death_date = "COMPOSITE_DEATH_DATE",
         elig_organ = "ORGAN",
         elig_tx_date = "TX_DATE", 
         elig_age = "AGE", 
         elig_val_dt_tcr = "VAL_DT_TCR", 
         elig_val_dt_trr = "VAL_DT_TRR", 
         elig_multiorg = "MULTIORG", 
         cov_a1 = "A1",
         cov_a2 = "A2",
         cov_age = "AGE",
         cov_age_don = "AGE_DON",
         exp_amis = "AMIS",
         cov_b1 = "B1",
         cov_b2 = "B2",
         cov_bmi_tcr = "BMI_TCR",
         exp_bmis = "BMIS",
         cov_c1 = "C1",
         cov_c2 = "C2",
         cov_cmv_don = "CMV_DON",
         cov_cmv_status = "CMV_STATUS",
         cov_cold_isch_ki = "COLD_ISCH_KI",
         cov_current_pra = "CURRENT_PRA",
         cov_da1 = "DA1",
         cov_da2 = "DA2",
         cov_dayswait_chron = "DAYSWAIT_CHRON",
         cov_db1 = "DB1",
         cov_db2 = "DB2",
         cov_ddr1 = "DDR1",
         cov_ddr2 = "DDR2",
         cov_diag = "DIAG_KI",
         cov_dialysis_date = "DIALYSIS_DATE",
         cov_don_ty = "DON_TY",
         cov_dr1 = "DR1",
         cov_dr2 = "DR2",
         exp_drmis = "DRMIS",
         cov_ebv_serostatus = "EBV_SEROSTATUS",
         cov_ebv_dna_don = "EBV_DNA_DON",
         cov_ebv_igg_cad_don = "EBV_IGG_CAD_DON",
         cov_ebv_igg_don = "EBV_IGG_DON",
         cov_ebv_igm_don = "EBV_IGM_DON",
         cov_ethcat = "ETHCAT",
         cov_func_stat_trr = "FUNC_STAT_TRR",
         cov_gender = "GENDER",
         cov_hgt_cm_tcr = "HGT_CM_TCR",
         cov_init_wgt_kg = "INIT_WGT_KG",
         cov_los = "LOS",
         cov_on_dialysis = "ON_DIALYSIS",
         cov_atgam_ind = "ATGAM_IND",
         cov_thymoglobulin_ind = "THYMOGLOBULIN_IND",
         cov_okt3_ind = "OKT3_IND",
         cov_campath_ind = "CAMPATH_IND",
         cov_simulect_ind = "SIMULECT_IND",
         cov_steroids_ind = "STEROIDS_IND",
         cov_envarsusxr_maint = "ENVARSUSXR_MAINT",
         cov_prograf_maint = "PROGRAF_MAINT",
         cov_genericcyclosporin_maint = "GENERICCYCLOSPORIN_MAINT",
         cov_cellcept_maint = "CELLCEPT_MAINT",
         cov_genericmyfortic_maint = "GENERICMYFORTIC_MAINT",
         cov_dc1 = "DC1",
         cov_dc2 = "DC2") %>% 
  
  ## Cohort Selection ----
  # Apply eligibility criteria (nrow = 1,208,213 non-unique patient observations)
  
  # Select kidney transplant records (n = 539,477 unique patients)
  filter(trr_id_code != "" & elig_organ == "KI") %>% 
  
  # Inclusion 1: Initial transplant between 1/1/2004 - 9/30/2024 (n = 365,280)
  group_by(pt_code) %>% 
  arrange(elig_tx_date) %>% 
  mutate(i = seq(1, n(), 1)) %>% 
  ungroup() %>% 
  filter(i == 1 & 
           elig_tx_date >= as.Date("2004-01-01") & 
           elig_tx_date <= as.Date("2024-09-30")) %>% 
  select(-i) %>% 
  
  # Inclusion 2: Recipient aged <18 years at time of transplant (n = 14,829)
  filter(elig_age < 18) %>% 
  
  # Exclusion 1: Multi-organ transplant recipient (n = 14,495)
  filter(elig_multiorg != "Y") %>% 
  
  # Exclusion 2: Missing HLA mismatch or allele data (n = 12,940)
  filter(complete.cases(exp_amis, exp_bmis, exp_drmis,
                        cov_a1, cov_b1, cov_dr1,
                        cov_a2, cov_b2, cov_dr2,
                        cov_da1, cov_db1, cov_ddr1,
                        cov_da2, cov_db2, cov_ddr2)) %>% 

  # Exclusion 3: Transplant and candidate registrations not validated by UNOS (n = 12,840)
  filter(is.na(elig_val_dt_tcr) == FALSE & 
           is.na(elig_val_dt_trr) == FALSE) 


## Variable Coding ---- 
data_kidney_clean <- data_kidney_cohort %>% 
  mutate(out_gstatus = as.logical(out_gstatus),
         exp_amis = as.factor(exp_amis),
         cov_bmi_tcr = factor(case_when(cov_bmi_tcr < 18.5 ~ "<18.5",
                                        cov_bmi_tcr >= 18.5 & cov_bmi_tcr <= 25 ~ "18.5-25",
                                        cov_bmi_tcr > 25 ~ ">25",
                                        TRUE ~ NA_character_),
                              levels = c("<18.5","18.5-25",">25")),
         exp_bmis = as.factor(exp_bmis),
         cov_cmv_don = as.factor(case_when(cov_cmv_don == "P" ~ "Positive",
                                            cov_cmv_don == "N" ~ "Negative",
                                            TRUE ~ NA)),
         cov_cmv_status = as.factor(case_when(cov_cmv_status == "P" ~ "Positive",
                                               cov_cmv_status == "N" ~ "Negative",
                                               TRUE ~ NA)),
         cov_current_pra = factor(case_when(cov_current_pra == 0 ~ "0",
                                            cov_current_pra %in% c(1:50) ~ "1 - 50",
                                            cov_current_pra %in% c(51:100) ~ "51 - 100",
                                            TRUE ~ NA_character_)),
         cov_diag = factor(case_when(cov_diag %in% c("3007", "3025", "3026", 
                                                     "3036", "3052") ~ "CAKUT",
                                     cov_diag %in% c("3000", "3001", "3002", 
                                                     "3003", "3004", "3005", 
                                                     "3006", "3024", "3033", 
                                                     "3034", "3035", "3040", 
                                                     "3041", "3042", "3043", 
                                                     "3051", "3054", "3055", 
                                                     "3068") ~ "Glomerular",
                                     cov_diag %in% c("3009", "3010", "3015",
                                                     "3031", "3032", "3064") ~ "Hereditary Nephritis",
                                     cov_diag %in% c("3008", "3028", "3029") ~ "Cystic Kidney Disease",
                                     cov_diag %in% c("999", "3013", "3014", 
                                                     "3016", "3018", "3019", 
                                                     "3020", "3021", "3027", 
                                                     "3030", "3037", "3044", 
                                                     "3045", "3046", "3047", 
                                                     "3048", "3050", "3053", 
                                                     "3057", "3058", "3059", 
                                                     "3060" ,"3062", "3063",
                                                     "3069", "3070", "3071",
                                                     "3072", "3073", "3074") ~ "Other",
                                     TRUE ~ NA_character_),
                           levels = c("CAKUT", "Hereditary Nephritis", 
                                      "Glomerular", "Cystic Kidney Disease",
                                      "Other")),
         cov_dialysis_need = as.logical(!is.na(cov_dialysis_date)),
         cov_dialysis_duration = as.numeric(elig_tx_date - cov_dialysis_date),
         cov_dialysis_duration_cat = as.factor(case_when(is.na(cov_dialysis_duration) ~ "N/A",
                                                         cov_dialysis_duration %in% (0:365) ~ "0 - 1 Year",
                                                         cov_dialysis_duration %in% (366:730) ~ "1 - 2 Years",
                                                         cov_dialysis_duration %in% (731:1095) ~ "2 - 3 Years",
                                                         cov_dialysis_duration %in% (1096:1460) ~ "3 - 4 Years",
                                                         cov_dialysis_duration > 1460 ~ "4+ Years")),
         cov_don_ty = as.factor(case_when(cov_don_ty == "C" ~ "Deceased",
                                          cov_don_ty == "L" ~ "Living")),
         exp_drmis = factor(exp_drmis),
         cov_ebv_serostatus = as.factor(case_when(cov_ebv_serostatus == "P" ~ "Positive",
                                                  cov_ebv_serostatus == "N" ~ "Negative",
                                                  TRUE ~ NA)),
         cov_ebv_serostatus_don = as.factor(case_when(cov_ebv_dna_don == "P" |
                                                        cov_ebv_igg_cad_don == "P" |
                                                        cov_ebv_igg_don == "P" | 
                                                        cov_ebv_igm_don == "P" ~ "Positive",
                                                      cov_ebv_dna_don == "N" |
                                                        cov_ebv_igg_cad_don == "N" |
                                                        cov_ebv_igg_don == "N" | 
                                                        cov_ebv_igm_don == "N" ~ "Negative",
                                                      TRUE ~ NA)),
         cov_ethcat = as.factor(case_when(cov_ethcat == 1 ~ "White, Non-Hispanic",
                                          cov_ethcat == 2 ~ "Black, Non-Hispanic",
                                          cov_ethcat == 4 ~ "Hispanic/Latino",
                                          cov_ethcat == 5 ~ "Asian, Non-Hispanic",
                                          cov_ethcat == 6 ~ "Amer Ind/Alaska Native, Non-Hispanic",
                                          cov_ethcat == 7 ~ "Native Hawaiian/other Pacific Islander, Non-Hispanic",
                                          cov_ethcat == 9 ~ "Multiracial, Non-Hispanic",
                                          TRUE ~ "Unknown")),
         cov_func_stat_trr = factor(case_when(cov_func_stat_trr == 996 | elig_age < 1 ~ "N/A (age < 1 year)",
                                              cov_func_stat_trr %in% c(3, 4010, 4020, 4030) ~ "10%-30%",
                                              cov_func_stat_trr %in% c(2, 4040, 4050) ~ "40%-60%",
                                              cov_func_stat_trr %in% c(1, 4060, 4070, 4080, 4090, 4100) ~ "70%-100%",
                                              TRUE ~ NA_character_),
                                    levels = c("10%-30%",
                                               "40%-60%",
                                               "70%-100%",
                                               "N/A (age < 1 year)")),
         cov_gender = factor(cov_gender),
         cov_c1 = case_when(cov_c1 %in% c(97, 99) ~ NA,
                            TRUE ~ cov_c1),
         cov_c2 = case_when(cov_c2 %in% c(97, 99) ~ NA,
                            TRUE ~ cov_c2),
         cov_dc1 = case_when(cov_dc1 %in% c(97, 99) ~ NA,
                            TRUE ~ cov_dc1),
         cov_dc2 = case_when(cov_dc2 %in% c(97, 99) ~ NA,
                            TRUE ~ cov_dc2),
         exp_cmis = as.factor(case_when(is.na(cov_c1) |  
                                           is.na(cov_c2) | 
                                           is.na(cov_dc1) | 
                                           is.na(cov_dc2) ~ NA_real_,
                                         TRUE ~ AlleleMismatchLevel(cov_c1, cov_c2,
                                                                    cov_dc1, cov_dc2) / 2))) %>% 
  
  # Remove variables that will not be used in analysis
  select(!starts_with("elig_"), 
         -c(out_pstatus, out_composite_death_date, out_ptime, cov_a1, cov_a2, 
            cov_b1, cov_b2,
            cov_c1, cov_c2, cov_da1, cov_da2, cov_db1, cov_db2, cov_ddr1, cov_ddr2,
            cov_dialysis_date, cov_dr1, cov_dr2, cov_ebv_dna_don, 
            cov_ebv_igg_cad_don, cov_ebv_igg_don, cov_ebv_igm_don, cov_on_dialysis,
            cov_dc1, cov_dc2)) %>% 
  
  # Arrange order of remaining variables in data table
  select(trr_id_code, pt_code,
         starts_with("out_"),
         starts_with("exp_"),
         starts_with("cov_"))

# Missing Data ----

## Exploration ----
# Assess percent missing values in each variable and determine variables with 
# a high amount (i.e., >30%) of missingness, which will be excluded from use 
# as predictors in imputation trees
pct_missing <- unlist(lapply(data_kidney_clean, function(x) sum(is.na(x))))/nrow(data_kidney_clean)
sort(pct_missing[pct_missing > 0], decreasing = TRUE)
VARS_EXCL_MI_PREDICTORS <- names(sort(pct_missing[pct_missing > 0.05], decreasing = TRUE))


## Imputation ----
### Parameters & Specification ---- 
# Run mice algorithm with no iterations
imp <- mice(data_kidney_clean, maxit = 0)

# Review for any problem variables that the package picked up.
# Note: mice() has appropriately identified several maintenance treatment  
# variables and an ID variable that are constants/non-informative.
imp$loggedEvents

# Extract predictorMatrix and methods of imputation 
predM <- imp$predictorMatrix
meth <- imp$method

# Set non-informative variables to zero in the predictor matrix
predM[, c("trr_id_code", "pt_code",   # Txplant & patient ID variables
          VARS_EXCL_MI_PREDICTORS,    # Variables with >5% missingness
          "cov_hgt_cm_tcr",           # Redundant with BMI
          "cov_init_wgt_kg",          # Redundant with BMI
          "cov_dialysis_need",        # Removed to facilitate convergence
          "cov_don_ty"                # Removed to facilitate convergence
          )] <- 0

# Specify random forest as the method for imputation of all variables
meth[!meth == ""] <- "rf"

# Specify variables that should not be imputed
meth[c("exp_cmis",
       "cov_hgt_cm_tcr",                 # Already adjusting for BMI
       "cov_init_wgt_kg",                # Already adjusting for BMI
       "cov_dialysis_duration",          # Missing is N/A; Already adjusting for categorical transformation
       "cov_atgam_ind",                  # Too much missingness to include as covariate
       "cov_thymoglobulin_ind",          # Too much missingness to include as covariate
       "cov_okt3_ind",                   # Too much missingness to include as covariate
       "cov_campath_ind",                # Too much missingness to include as covariate
       "cov_simulect_ind",               # Too much missingness to include as covariate
       "cov_steroids_ind",               # Too much missingness to include as covariate
       "cov_envarsusxr_maint",           # Too much missingness to include as covariate
       "cov_prograf_maint",              # Too much missingness to include as covariate
       "cov_genericcyclosporin_maint",   # Too much missingness to include as covariate
       "cov_cellcept_maint",             # Too much missingness to include as covariate
       "cov_genericmyfortic_maint"      # Too much missingness to include as covariate
       )] <- ""

### Run Imputation Algorithm -----
# data_imputed <- mice(data_kidney_clean,
#                      m = 10,
#                      maxit = 10,
#                      predictorMatrix = predM,
#                      method = meth,
#                      seed = 777,
#                      visitSequence = "monotone")

# Save/Load mids object
# save(data_imputed, file = paste0(DIR_DATA, "data_imputed.Rdata"))
load(file = paste0(DIR_DATA, "data_imputed.Rdata"))

### Algorithm Convergence Check -----
# Review for any problems that occurred during the imputation procedure
# Note: ethnicity variable appropriately dropped from PRA imputation model 
# due to sparsity
# View(data_imputed$loggedEvents)

### Inspection of Imputed Values ----
# QC - Evaluate convergence of algorithm by plotting:
# mean + SD of each imputed variable (y-axis) vs. iteration (x-axis)
# for every dataset (separate lines)
# 
# Inspection reveals no major trends in means and that the curves mix well.
# A lack of trend in the means suggests the algorithm converged.
# When inspecting the variability (SD) in the plots, a high degree of variation  
# between the datasets (lines) vs within the datasets indicates that there is a 
# relatively high increase in variance due to missing data. This is not evident.
plot(data_imputed, layout = c(5,5))

### Extract Imputed Datasets ----
data_imputed_all <- complete(data_imputed, "long")

# Analysis -----

## Kaplan Meier Estimator (IP Weighted) -----

### PS Model Specification ----
# Specify list of confounders for propensity score model
CONFOUNDERS <- data_kidney_clean %>% 
  select(starts_with("cov_"), 
         -c(cov_hgt_cm_tcr, cov_init_wgt_kg, cov_dialysis_duration)) %>% 
  select(-ends_with("_ind")) %>% 
  select(-ends_with("_maint")) %>% 
  names()
       
EXPOSURES <- data_kidney_clean %>% 
  select(starts_with("exp_")) %>% 
  select(-exp_cmis) %>% 
  names()

ALLELES <- c("a", "b", "dr")

# Specify formulas for propensity score models for each allele (HLA-A, -B, & -DR)
FORMULA_PS_AMIS <- as.formula(paste("exp_amis ~ ", paste(CONFOUNDERS, collapse = " + ")))
FORMULA_PS_BMIS <- as.formula(paste("exp_bmis ~ ", paste(CONFOUNDERS, collapse = " + ")))
FORMULA_PS_DRMIS <- as.formula(paste("exp_drmis ~ ", paste(CONFOUNDERS, collapse = " + ")))

# Estimate IP weighted KM curves for each allele (Run time ~ 1 hour)
### HLA-A ----
# adjusted_surv_amis <- adjustedsurv(data = data_imputed,
#                                    variable = "exp_amis",
#                                    ev_time = "out_gtime",
#                                    event = "out_gstatus",
#                                    method = "iptw_km",
#                                    stabilize = TRUE,
#                                    #trim_quantiles = c(0.01, 0.99),
#                                    treatment_model = FORMULA_PS_AMIS,
#                                    conf_int = TRUE, 
#                                    conf_level = 0.95,
#                                    times = NULL, 
#                                    bootstrap = TRUE,
#                                    n_boot = 500,
#                                    n_cores = 8)

# Save object 
# save(adjusted_surv_amis, file = paste0(DIR_DATA, "adjusted_surv_amis.Rdata"))
load(file = paste0(DIR_DATA, "adjusted_surv_amis.Rdata"))


### HLA-B ----
# adjusted_surv_bmis <- adjustedsurv(data = data_imputed,
#                                    variable = "exp_bmis",
#                                    ev_time = "out_gtime",
#                                    event = "out_gstatus",
#                                    method = "iptw_km",
#                                    stabilize = TRUE,
#                                    #trim_quantiles = c(0.01, 0.99),
#                                    treatment_model = FORMULA_PS_BMIS,
#                                    conf_int = TRUE, 
#                                    conf_level = 0.95,
#                                    times = NULL, 
#                                    bootstrap = TRUE,
#                                    n_boot = 500,
#                                    n_cores = 8)

# Save object
# save(adjusted_surv_bmis, file = paste0(DIR_DATA, "adjusted_surv_bmis.Rdata"))
load(file = paste0(DIR_DATA, "adjusted_surv_bmis.Rdata"))

### HLA-DR ----
# adjusted_surv_drmis <- adjustedsurv(data = data_imputed,
#                                    variable = "exp_drmis",
#                                    ev_time = "out_gtime",
#                                    event = "out_gstatus",
#                                    method = "iptw_km",
#                                    stabilize = TRUE,
#                                    #trim_quantiles = c(0.01, 0.99),
#                                    treatment_model = FORMULA_PS_DRMIS,
#                                    conf_int = TRUE, 
#                                    conf_level = 0.95,
#                                    times = NULL, 
#                                    bootstrap = TRUE,
#                                    n_boot = 500,
#                                    n_cores = 8)

# Save object
# save(adjusted_surv_drmis, file = paste0(DIR_DATA, "adjusted_surv_drmis.Rdata"))
load(file = paste0(DIR_DATA, "adjusted_surv_drmis.Rdata"))


### KM Plots ----
# Plot Pooled and Weighted Kaplan-Meier Estimates
plot_hla_a_adj_km <- plot(adjusted_surv_amis, 
                          risk_table = TRUE,
                          risk_table_stratify = TRUE,
                          use_boot = TRUE,
                          conf_int = TRUE,
                          xlab = "Time Since Transplant (Days)",
                          ylab = "Weighted Survival Probability (%)",
                          legend.title = "No. Mismatches",
                          legend.position = "bottom",
                          conf_int_alpha = 0.25,
                          x_breaks = c(0, 1825, 3650, 5475, 7300),
                          y_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                          risk_table_ylab = "",
                          gg_theme = theme_minimal())

plot_hla_b_adj_km <- plot(adjusted_surv_bmis, 
                          risk_table = TRUE,
                          risk_table_stratify = TRUE,
                          use_boot = TRUE,
                          conf_int = TRUE,
                          xlab = "Time Since Transplant (Days)",
                          ylab = "Weighted Survival Probability (%)",
                          legend.title = "No. Mismatches",
                          legend.position = "bottom",
                          conf_int_alpha = 0.25,
                          x_breaks = c(0, 1825, 3650, 5475, 7300),
                          y_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                          risk_table_ylab = "",
                          gg_theme = theme_minimal())

plot_hla_dr_adj_km <- plot(adjusted_surv_drmis, 
                          risk_table = TRUE,
                          risk_table_stratify = TRUE,
                          use_boot = TRUE,
                          conf_int = TRUE,
                          xlab = "Time Since Transplant (Days)",
                          ylab = "Weighted Survival Probability (%)",
                          legend.title = "No. Mismatches",
                          legend.position = "bottom",
                          conf_int_alpha = 0.25,
                          x_breaks = c(0, 1825, 3650, 5475, 7300),
                          y_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                          risk_table_ylab = "",
                          gg_theme = theme_minimal())

## IP Weight Inspection ----
# Collect object names in a vector to facilitate inspection of weights below
adjustedsurv_obj <- c("adjusted_surv_amis", "adjusted_surv_bmis", "adjusted_surv_drmis")

# Print weight distribution in each imputed dataset to check for extreme values, 
# suggestive of positivity violations, or means deviating from 1, suggestive of
# a mis-specification.
# Note: Weight means are all approximately 1 and there are no extreme values. 
for (o in 1:length(adjustedsurv_obj)) {
  for (m in 1:10) {
    print(paste("Object:", adjustedsurv_obj[o], "- Imputed Dataset No.", m ))
    print(summary(get(adjustedsurv_obj[o])[["mids_analyses"]][[m]][["weights"]]))
  }
}

## IP Weight Extraction -----
# Create data frame shell to store IP weights
weights <- as.data.frame(cbind(data_imputed_all$.id,
                               data_imputed_all$.imp,
                               rep(NA, length(data_imputed_all$trr_id_code)),
                               rep(NA, length(data_imputed_all$trr_id_code)),
                               rep(NA, length(data_imputed_all$trr_id_code)))) %>% 
  rename(.id = V1, 
         .imp = V2, 
         ipw_hla_a = V3, 
         ipw_hla_b = V4, 
         ipw_hla_dr = V5)

# Populate shell with IP weights from HLA-A, HLA-B, and HLA-DR analyses
for (o in 1:length(adjustedsurv_obj)) {
  for (m in 1:10) {
    if (m == 1){
      temp_weights <- as.data.frame(get(adjustedsurv_obj[o])[["mids_analyses"]][[m]][["weights"]])
    }
    if (m > 1){
      temp_weights <- rbind(temp_weights, 
                            as.data.frame(get(adjustedsurv_obj[o])[["mids_analyses"]][[m]][["weights"]]))
    }
    if (m == 10) {
      weights[, o + 2] <- temp_weights
      rm(temp_weights)
      }
  }
}

## Restricted Mean Survival Time (Weighted) ---- 
# Runtime ~ 2 mins

### HLA-A ----
# adjusted_rmst_amis <- adjusted_rmst(adjsurv = adjusted_surv_amis, 
#                                     from = 0,
#                                     to = c(365, 1825, 3650, 5475), 
#                                     conf_int = TRUE,
#                                     conf_level = 0.95)

# Save object
# save(adjusted_rmst_amis, file = paste0(DIR_DATA, "adjusted_rmst_amis.Rdata"))
load(file = paste0(DIR_DATA, "adjusted_rmst_amis.Rdata"))

### HLA-B ----
# adjusted_rmst_bmis <- adjusted_rmst(adjsurv = adjusted_surv_bmis, 
#                                     from = 0,
#                                     to = c(365, 1825, 3650, 5475), 
#                                     conf_int = TRUE,
#                                     conf_level = 0.95)

# Save object
# save(adjusted_rmst_bmis, file = paste0(DIR_DATA, "adjusted_rmst_bmis.Rdata"))
load(file = paste0(DIR_DATA, "adjusted_rmst_bmis.Rdata"))


### HLA-DR ----
# adjusted_rmst_drmis <- adjusted_rmst(adjsurv = adjusted_surv_drmis, 
#                                     from = 0,
#                                     to = c(365, 1825, 3650, 5475), 
#                                     conf_int = TRUE,
#                                     conf_level = 0.95)

# Save object
# save(adjusted_rmst_drmis, file = paste0(DIR_DATA, "adjusted_rmst_drmis.Rdata"))
load(file = paste0(DIR_DATA, "adjusted_rmst_drmis.Rdata"))

# Output ----

## Baseline Characteristics ----
# Define Table 1 elements for CreateTableOne() function
names(data_imputed_all)
vars <- c("cov_age",
          "cov_age_don",
          "cov_init_wgt_kg",
          "cov_hgt_cm_tcr",
          "cov_bmi_tcr",
          "cov_gender",
          "cov_ethcat",
          "cov_diag",
          "cov_dialysis_need",
          "cov_dialysis_duration",
          "cov_cmv_status",
          "cov_cmv_don",
          "cov_ebv_serostatus",
          "cov_ebv_serostatus_don",
          "cov_atgam_ind",
          "cov_thymoglobulin_ind",
          "cov_okt3_ind",
          "cov_campath_ind",
          "cov_simulect_ind",
          "cov_steroids_ind",
          "cov_envarsusxr_maint",
          "cov_prograf_maint",
          "cov_genericcyclosporin_maint",
          "cov_cellcept_maint",
          "cov_genericmyfortic_maint",
          "cov_current_pra",
          "cov_don_ty",
          "cov_cold_isch_ki",
          "cov_dayswait_chron",
          "cov_los",
          "cov_func_stat_trr")

### Cohort Before Imputation and Before Weighting ----
#### Overall & Stratified by HLA Mismatch ----
for (e in 1:length(EXPOSURES)) {

    # Generate baseline table
    temp_baseline_table <- CreateTableOne(vars = vars, 
                                          strata = EXPOSURES[e],
                                          data = data_kidney_clean,
                                          smd = TRUE,
                                          test = FALSE,
                                          addOverall = TRUE)
    
    # Extract SMDs for all pairwise comparisons across strata
    temp_smd <- ExtractSmd(temp_baseline_table)
    
    # Print baseline table
    temp_baseline_table_print <- print(temp_baseline_table,
                                       showAllLevels = TRUE,
                                       missing = TRUE,
                                       smd = TRUE)
    
    # Output baseline table
    write.csv(temp_baseline_table_print, 
              file = paste0(DIR_OUTPUT, "Baseline Table - Prior to Imputation & Weighting - ",
                            "Stratified by HLA-", toupper(ALLELES[e]), " Mismatch.csv"))
    
    # Output SMDs
    write.csv(temp_smd, 
              file = paste0(DIR_OUTPUT, "Baseline Table - Prior to Imputation & Weighting - ",
                            "Stratified by HLA-", toupper(ALLELES[e]), " Mismatch - SMD.csv"))
}



### Cohort After Imputation and Before Weighting ----

#### Overall & Stratified by HLA Mismatch ----
for (e in 1:length(EXPOSURES)) {
  for (m in 1:length(unique(data_imputed_all$.imp))){
    
    # Select subset of imputed dataset 'm'
    temp_data <- filter(data_imputed_all, .imp == m)
    
    # Generate baseline table
    temp_baseline_table <- CreateTableOne(vars = vars, 
                                          strata = EXPOSURES[e],
                                          data = temp_data,
                                          smd = TRUE,
                                          test = FALSE,
                                          addOverall = TRUE)
    
    # Extract SMDs for all pairwise comparisons across strata
    temp_smd <- ExtractSmd(temp_baseline_table)
    
    # Print baseline table
    temp_baseline_table_print <- print(temp_baseline_table,
                                       showAllLevels = TRUE,
                                       missing = TRUE,
                                       smd = TRUE)
    
    # Output baseline table
    write.csv(temp_baseline_table_print, 
              file = paste0(DIR_OUTPUT, "Baseline Table - Imputed Dataset ", m, " - Unweighted - ",
                            "Stratified by HLA-", toupper(ALLELES[e]), 
                            " Mismatch.csv"))
    
    # Output SMDs
    write.csv(temp_smd, 
              file = paste0(DIR_OUTPUT, "Baseline Table - Imputed Dataset ", m, " - Unweighted - ",
                            "Stratified by HLA-", toupper(ALLELES[e]), 
                            " Mismatch - SMD.csv"))
    
  }
}

### Cohort After Imputation and After Weighting ----

# Merge IP weights onto analytic file
data_imputed_all_weights <- left_join(data_imputed_all, weights, by = c(".imp", ".id"))

# Specify weight variables for each exposure to facilitate dynamic calling in loop
weight_varnames <- c("ipw_hla_a", "ipw_hla_b", "ipw_hla_dr")

# Drop variables that don't have enough non-missing values, which cause errors 
# in svyCreateTableOne()
vars_mod <- vars[!vars %in% c('cov_okt3_ind', 'cov_genericcyclosporin_maint')]


# Create survey design objects to be fed to svyCreateTableOne() function
for (e in 1:length(EXPOSURES)) {
  for (m in 1:length(unique(data_imputed_all$.imp))){
    
    # Select subset of imputed dataset 'm'
    temp_data <- filter(data_imputed_all_weights, .imp == m)
    
    # Create survey design object
    
    temp_survey_design <- svydesign(ids = ~ 1, 
                                    strata = ~ get(EXPOSURES[e]),
                                    weights = ~ temp_data[[weight_varnames[e]]],
                                    data = temp_data)
    
    # Generate baseline table
    temp_baseline_table <- svyCreateTableOne(vars = vars_mod, 
                                             strata = EXPOSURES[e],
                                             data = temp_survey_design,
                                             smd = TRUE,
                                             test = FALSE,
                                             addOverall = TRUE)
    
    # Extract SMDs for all pairwise comparisons across strata
    temp_smd <- ExtractSmd(temp_baseline_table)
    
    # Print baseline table
    temp_baseline_table_print <- print(temp_baseline_table,
                                       showAllLevels = TRUE,
                                       missing = TRUE,
                                       smd = TRUE)
    
    # Output baseline table
    write.csv(temp_baseline_table_print, 
              file = paste0(DIR_OUTPUT, "Baseline Table - Imputed Dataset ", m, " - Weighted - ",
                            "Stratified by HLA-", toupper(ALLELES[e]), 
                            " Mismatch.csv"))
    
    # Output SMDs
    write.csv(temp_smd, 
              file = paste0(DIR_OUTPUT, "Baseline Table - Imputed Dataset ", m, " - Weighted - ",
                            "Stratified by HLA-", toupper(ALLELES[e]), 
                            " Mismatch - SMD.csv"))
    
  }
}

## KM Plots ----
ggsave(paste0(DIR_OUTPUT, "Kaplan-Meier Plot - Weighted & Pooled - Stratified by HLA-A Mismatch.pdf"), plot = plot_hla_a_adj_km, width = 8, height = 6)

ggsave(paste0(DIR_OUTPUT, "Kaplan-Meier Plot - Weighted & Pooled - Stratified by HLA-B Mismatch.pdf"), plot = plot_hla_b_adj_km, width = 8, height = 6)

ggsave(paste0(DIR_OUTPUT, "Kaplan-Meier Plot - Weighted & Pooled - Stratified by HLA-DR Mismatch.pdf"), plot = plot_hla_dr_adj_km, width = 8, height = 6)

## RMST -----
adjustedrmst_obj <- c("adjusted_rmst_amis", "adjusted_rmst_bmis", "adjusted_rmst_drmis")

for (o in 1:length(adjustedrmst_obj)) {
  
  # Index 
  if (o == 1){
    i = 1
  }
  
  if (o > 1){
    i <- i + 1
  }
  
  # Clean up output
  table_clean <- get(adjustedrmst_obj[o]) %>% 
    mutate(`Time (Years)`= to / 365,
           `Mismatch No.` = group,
           `RMST - Days (95% CI)` = paste0(round(rmst), " (", round(ci_lower), 
                                           ", ", round(ci_upper), ")"),
           `RMST - Years (95% CI)` = paste0(round(rmst / 365, 2), 
                                            " (", round(ci_lower / 365, 2), ", ", 
                                            round(ci_upper / 365, 2), ")")) %>% 
    select(-c(to, group, rmst, se, ci_lower, ci_upper, n_boot))
  
  # Save file
  write.csv(table_clean, 
            file = paste0(DIR_OUTPUT, "RMST Graft Survival & Mortality - Stratified by HLA- ", toupper(ALLELES[i])," Mismatch.csv"),
            row.names = FALSE)
}

# Misc Code, Notes, and QC ----
# Table delivered for quantitative biologist:
kidney_table_grantham <- select(data_kidney_cohort, 
                                 trr_id_code, pt_code, 
                                 out_gstatus, out_gtime, 
                                 out_pstatus, out_ptime, out_composite_death_date,
                                 cov_amis, cov_bmis, cov_drmis,
                                 cov_a1, cov_a2,
                                 cov_da1, cov_da2,
                                 cov_b1, cov_b2, 
                                 cov_db1, cov_db2,
                                 cov_dr1, cov_dr2,
                                 cov_ddr1, cov_ddr2)

write.csv(kidney_table_grantham, 
          file = paste0(DIR_OUTPUT, "Kidney Cohort for Grantham Distance.csv"),
          row.names = FALSE)

write.csv(formats, 
          file = paste0(DIR_OUTPUT, "Formats (Variable Coding Definitions).csv"),
          row.names = FALSE)

# QC: Spot check output from AlleleMismatchLevel() function when applied to HLA-C
data_kidney_clean %>% 
  select(cov_c1, cov_c2,
         cov_dc1, cov_dc2,
         cov_cmis) %>% 
  View()


  
  
  
  
  
  
  
  
  

