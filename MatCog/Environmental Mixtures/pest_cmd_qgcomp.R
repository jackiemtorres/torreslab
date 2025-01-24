################################################################################################
# NAME:         Pesticides and cardiometabolic risk factors (pest_cmd_qgcomp.R)
# AUTHORS:      Lucia Calderon
# CREATED:      12/16/2024
# PURPOSE:      Script to run quantile G-computation to analyze associations of mixtures of pesticides, 
#               measured using Pesticide Use Reporting Data and urinary biomarkers, and cardio-metabolic risk 
#               factors in CHAMACOS moms
# COVARIATES:   MC1: Age, marital status, work status, poverty status, urinary DAPs
#               Baseline: Nativity/age arrived in US, education
#               (age_mc1 ageusa_cat_bl meducat_bl marcat_mc1 agwork_mc1 povcat_rc_mc1, log_daps)
# SENSITIVITY ANALYSES:
#     
# UPDATES: 
#  
################################################################################################

# Clear workspace
#rm(list=ls())

# Set working directory
#setwd("/Users/lcalderon1/Library/CloudStorage/Box-Box/CHAM Maternal Cognition Study/Analyses and Manuscripts/Lucia's analyses/CHAM PUR Cardiometabolic/Analysis 2")


library(tableone)
library(fastDummies)
library(R2jags)  
library(coda)
library(dplyr)
library(haven)
library(multcomp)
library(ggplot2)
library(qgcomp)
library(gfoRmula)
library(broom)
library(boot)

# Read in Data from Stata (this dataset only includes the analytic sample)
data<-read_dta("../Code/pest_cmd_analysis.dta") 

# Drop outliers for continuous outcomes
data$mbmi_mc1[data$mbmi_mc1_out==1]<-NA
data$mwaist_mc1[data$mwaist_mc1_out==1]<-NA
data$mbpsys_mc1[data$mbpsys_mc1_out==1]<-NA
data$mbpdia_mc1[data$mbpdia_mc1_out==1]<-NA
data$mbppul_mc1[data$mbppul_mc1_out==1]<-NA
data$mbpmap_mc1[data$mbpmap_mc1_out==1]<-NA
data$gluc_mc1[data$gluc_mc1_out==1]<-NA
data$hba1c_mc1[data$hba1c_mc1_out==1]<-NA
data$log_trig_mc1[data$log_trig_mc1_out==1]<-NA
data$hdl_mc1[data$hdl_mc1_out==1]<-NA
data$log_crp_mc1[data$log_crp_mc1_out==1]<-NA
data$log_il6_mc1[data$log_il6_mc1_out==1]<-NA

# Exclude those taking meds from relevant continuous outcomes
data$mbpsys_mc1[data$mbpsys_mc1_med==1]<-NA
data$mbpdia_mc1[data$mbpdia_mc1_med==1]<-NA
data$mbppul_mc1[data$mbppul_mc1_med==1]<-NA
data$mbpmap_mc1[data$mbpmap_mc1_med==1]<-NA
data$gluc_mc1[data$gluc_mc1_med==1]<-NA
data$hba1c_mc1[data$hba1c_mc1_med==1]<-NA
data$log_trig_mc1[data$log_trig_mc1_med==1]<-NA
data$hdl_mc1[data$hdl_mc1_med==1]<-NA

# Exclude non-fasting glucose
data$gluc_mc1[data$gluc_mc1_nf==1]<-NA

# Table 1 to check consistency with Stata analysis 
myvars<- c("age_mc1","ageusa_cat_bl","meducat_bl","marcat_mc1","agwork_mc1","povcat_rc_mc1", "tot",
           "mhtn_mc1","diabetes_mc1","mobese_mc1","waist_ind_mc1","bp_ind_mc1","gluc_ind_mc1","trig_ind_mc1","chdl_ind_mc1","metsyn_mc1","crp_cat_mc1",
           "mbpsys_mc1","mbpdia_mc1","mbppul_mc1","mbpmap_mc1","mbmi_mc1","mwaist_mc1","gluc_mc1","hba1c_mc1","log_trig_mc1","hdl_mc1","log_crp_mc1","log_il6_mc1")
mycatvars<-c("ageusa_cat_bl","meducat_bl","marcat_mc1","agwork_mc1","povcat_rc_mc1",
             "mhtn_mc1","diabetes_mc1","mobese_mc1","waist_ind_mc1","bp_ind_mc1","gluc_ind_mc1","trig_ind_mc1","chdl_ind_mc1","metsyn_mc1","crp_cat_mc1")
table1<-CreateTableOne(vars=myvars, data=data, factorVars=mycatvars)
tab1print<-print(table1, showAlllevels=T, catDigits=1, contDigits=1, pDigits=1, test=F, printToggle = FALSE)

#Subset out relevant covariates
data<-dummy_cols(data, select_columns = c("ageusa_cat_bl","meducat_bl","agwork_mc1","povcat_rc_mc1"), 
                 remove_first_dummy = TRUE, remove_selected_columns = FALSE, 
                 remove_most_frequent_dummy = FALSE, ignore_na = FALSE)
W <- subset(data, select = c("age_mc1","ageusa_cat_bl_2","ageusa_cat_bl_3","meducat_bl_2","meducat_bl_3",
                             "marcat_mc1","agwork_mc1_2","agwork_mc1_3","povcat_rc_mc1_2","povcat_rc_mc1_3", "log_tot"))
W<-as.matrix(W)

#Exposures are log-10 transformed 
X<- subset(data, select = c("log_ace","log_mala","log_bens",
                            "log_perm",
                            "log_meth",
                            "log_imid",
                            "log_mn",
                            "log_gly"))

X <- as.matrix(X) 

##Vectorize individual outcomes of interest from Y
Y_sys<-as.vector(data$mbpsys_mc1)
Y_dia<-as.vector(data$mbpdia_mc1)
Y_pul<-as.vector(data$mbppul_mc1)
Y_map<-as.vector(data$mbpmap_mc1)
Y_bmi<-as.vector(data$mbmi_mc1)
Y_wai<-as.vector(data$mwaist_mc1)
Y_glu<-as.vector(data$gluc_mc1)
Y_hba<-as.vector(data$hba1c_mc1)
Y_tri<-as.vector(data$log_trig_mc1)
Y_hdl<-as.vector(data$hdl_mc1)
Y_crp<-as.vector(data$log_crp_mc1)
Y_il6<-as.vector(data$log_il6_mc1)

Y_hyper<-as.vector(data$mhtn_mc1)
Y_diabe<-as.vector(data$diabetes_mc1)
Y_obese<-as.vector(data$mobese_mc1)
Y_waind<-as.vector(data$waist_ind_mc1)
Y_bpind<-as.vector(data$bp_ind_mc1)
Y_glind<-as.vector(data$gluc_ind_mc1)
Y_trind<-as.vector(data$trig_ind_mc1)
Y_chind<-as.vector(data$chdl_ind_mc1)
Y_metab<-as.vector(data$metsyn_mc1)
Y_crpca<-as.vector(data$crp_cat_mc1)

##Complete Case Flag for each Outcome
is.cc.sys <- complete.cases(cbind(W,X,Y_sys))
is.cc.dia <- complete.cases(cbind(W,X,Y_dia))
is.cc.pul <- complete.cases(cbind(W,X,Y_pul))
is.cc.map <- complete.cases(cbind(W,X,Y_map))
is.cc.bmi <- complete.cases(cbind(W,X,Y_bmi))
is.cc.wai <- complete.cases(cbind(W,X,Y_wai))
is.cc.glu <- complete.cases(cbind(W,X,Y_glu))
is.cc.hba <- complete.cases(cbind(W,X,Y_hba))
is.cc.tri <- complete.cases(cbind(W,X,Y_tri))
is.cc.hdl <- complete.cases(cbind(W,X,Y_hdl))
is.cc.crp <- complete.cases(cbind(W,X,Y_crp))
is.cc.il6 <- complete.cases(cbind(W,X,Y_il6))

is.cc.hyper <- complete.cases(cbind(W,X,Y_hyper))
is.cc.diabe <- complete.cases(cbind(W,X,Y_diabe))
is.cc.obese <- complete.cases(cbind(W,X,Y_obese))
is.cc.waind <- complete.cases(cbind(W,X,Y_waind))
is.cc.bpind <- complete.cases(cbind(W,X,Y_bpind))
is.cc.glind <- complete.cases(cbind(W,X,Y_glind))
is.cc.trind <- complete.cases(cbind(W,X,Y_trind))
is.cc.chind <- complete.cases(cbind(W,X,Y_chind))
is.cc.metab <- complete.cases(cbind(W,X,Y_metab))
is.cc.crpca <- complete.cases(cbind(W,X,Y_crpca))

##Subset outcomes as complete cases
Y.cc.sys <- subset(Y_sys,is.cc.sys) 
Y.cc.dia <- subset(Y_dia,is.cc.dia) 
Y.cc.pul <- subset(Y_pul,is.cc.pul) 
Y.cc.map <- subset(Y_map,is.cc.map) 
Y.cc.bmi <- subset(Y_bmi,is.cc.bmi)
Y.cc.wai <- subset(Y_wai,is.cc.wai)
Y.cc.glu <- subset(Y_glu,is.cc.glu)
Y.cc.hba <- subset(Y_hba,is.cc.hba)
Y.cc.tri <- subset(Y_tri,is.cc.tri)
Y.cc.hdl <- subset(Y_hdl,is.cc.hdl)
Y.cc.crp <- subset(Y_crp,is.cc.crp)
Y.cc.il6 <- subset(Y_il6,is.cc.il6)

Y.cc.hyper <- subset(Y_hyper,is.cc.hyper)
Y.cc.diabe <- subset(Y_diabe,is.cc.diabe)
Y.cc.obese <- subset(Y_obese,is.cc.obese)
Y.cc.waind <- subset(Y_waind,is.cc.waind)
Y.cc.bpind <- subset(Y_bpind,is.cc.bpind)
Y.cc.glind <- subset(Y_glind,is.cc.glind)
Y.cc.trind <- subset(Y_trind,is.cc.trind)
Y.cc.chind <- subset(Y_chind,is.cc.chind)
Y.cc.metab <- subset(Y_metab,is.cc.metab)
Y.cc.crpca <- subset(Y_crpca,is.cc.crpca)

# Set up data frames for regression models
data.sys <- as.data.frame(cbind(Y.cc.sys, X[is.cc.sys,], W[is.cc.sys,]))
data.dia <- as.data.frame(cbind(Y.cc.dia, X[is.cc.dia,], W[is.cc.dia,]))
data.pul <- as.data.frame(cbind(Y.cc.pul, X[is.cc.pul,], W[is.cc.pul,]))
data.map <- as.data.frame(cbind(Y.cc.map, X[is.cc.map,], W[is.cc.map,]))
data.bmi <- as.data.frame(cbind(Y.cc.bmi, X[is.cc.bmi,], W[is.cc.bmi,]))
data.wai <- as.data.frame(cbind(Y.cc.wai, X[is.cc.wai,], W[is.cc.wai,]))
data.glu <- as.data.frame(cbind(Y.cc.glu, X[is.cc.glu,], W[is.cc.glu,]))
data.hba <- as.data.frame(cbind(Y.cc.hba, X[is.cc.hba,], W[is.cc.hba,]))
data.tri <- as.data.frame(cbind(Y.cc.tri, X[is.cc.tri,], W[is.cc.tri,]))
data.hdl <- as.data.frame(cbind(Y.cc.hdl, X[is.cc.hdl,], W[is.cc.hdl,]))
data.crp <- as.data.frame(cbind(Y.cc.crp, X[is.cc.crp,], W[is.cc.crp,]))
data.il6 <- as.data.frame(cbind(Y.cc.il6, X[is.cc.il6,], W[is.cc.il6,]))

data.hyper <- as.data.frame(cbind(Y.cc.hyper, X[is.cc.hyper,], W[is.cc.hyper,]))
data.diabe <- as.data.frame(cbind(Y.cc.diabe, X[is.cc.diabe,], W[is.cc.diabe,]))
data.obese <- as.data.frame(cbind(Y.cc.obese, X[is.cc.obese,], W[is.cc.obese,]))
data.waind <- as.data.frame(cbind(Y.cc.waind, X[is.cc.waind,], W[is.cc.waind,]))
data.bpind <- as.data.frame(cbind(Y.cc.bpind, X[is.cc.bpind,], W[is.cc.bpind,]))
data.glind <- as.data.frame(cbind(Y.cc.glind, X[is.cc.glind,], W[is.cc.glind,]))
data.trind <- as.data.frame(cbind(Y.cc.trind, X[is.cc.trind,], W[is.cc.trind,]))
data.chind <- as.data.frame(cbind(Y.cc.chind, X[is.cc.chind,], W[is.cc.chind,]))
data.metab <- as.data.frame(cbind(Y.cc.metab, X[is.cc.metab,], W[is.cc.metab,]))
data.crpca <- as.data.frame(cbind(Y.cc.crpca, X[is.cc.crpca,], W[is.cc.crpca,]))

## Testing qgcomp package (From example 3 of A Keil vignette https://cran.r-project.org/web/packages/qgcomp/vignettes/qgcomp-vignette.html)

  # Mixture names
  Xnm.all <-c("log_ace","log_mala","log_bens",
            "log_perm",
            "log_meth",
            "log_imid",
            "log_mn",
            "log_gly") ##all exposures

#### Continuous outcomes (starting with glucose)

## qgcomp.glm.noboot (no bootstrapping)  
qc.fit.glu <- qgcomp.glm.noboot(Y.cc.glu ~ #outcome
                                  age_mc1 + #covars
                                  ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                  meducat_bl_2 + meducat_bl_3 +
                                  marcat_mc1 +
                                  agwork_mc1_2 + agwork_mc1_3 +
                                  povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                  log_tot +
                                  log_ace + log_mala + log_bens + #exposures
                                  log_perm +
                                  log_meth +
                                  log_imid +
                                  log_mn +
                                  log_gly,
                                expnms=Xnm.all, #exposure names
                                data.glu, #dataset
                                family=gaussian(), 
                                q=4)

qc.fit.glu

## KM's code:
# underlying model of qgcomp fit
qc.fit.glu$qx$Y.cc.glu <- qc.fit.glu$fit$data$Y.cc.glu # bring outcome back into quantized data

#add the covariates too
qc.fit.glu$qx$age_mc1 <- qc.fit.glu$fit$data$age_mc1 
qc.fit.glu$qx$ageusa_cat_bl_2 <- qc.fit.glu$fit$data$ageusa_cat_bl_2 
qc.fit.glu$qx$ageusa_cat_bl_3 <- qc.fit.glu$fit$data$ageusa_cat_bl_3 
qc.fit.glu$qx$meducat_bl_2 <- qc.fit.glu$fit$data$meducat_bl_2 
qc.fit.glu$qx$meducat_bl_3 <- qc.fit.glu$fit$data$meducat_bl_3 
qc.fit.glu$qx$marcat_mc1 <- qc.fit.glu$fit$data$marcat_mc1 
qc.fit.glu$qx$agwork_mc1_2 <- qc.fit.glu$fit$data$agwork_mc1_2 
qc.fit.glu$qx$agwork_mc1_3 <- qc.fit.glu$fit$data$agwork_mc1_3 
qc.fit.glu$qx$povcat_rc_mc1_2 <- qc.fit.glu$fit$data$povcat_rc_mc1_2 
qc.fit.glu$qx$povcat_rc_mc1_3 <- qc.fit.glu$fit$data$povcat_rc_mc1_3 
qc.fit.glu$qx$log_tot <- qc.fit.glu$fit$data$log_tot 

# adjusted coefficients from underlying model
# absolute change in outcome by increasing each exposure by 1
newfit.glu <- lm(Y.cc.glu ~ age_mc1 + #covars
                   ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                   meducat_bl_2 + meducat_bl_3 +
                   marcat_mc1 +
                   agwork_mc1_2 + agwork_mc1_3 +
                   povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                   log_tot +
                   log_ace_q + log_mala_q + log_bens_q + #exposures
                   log_perm_q +
                   log_meth_q +
                   log_imid_q +
                   log_mn_q +
                   log_gly_q, data = qc.fit.glu$qx)

newfit.glu ##this is where the individual coefficients are
tidy(newfit.glu, conf.int = T) # and confidence intervals

        ## QUESTIONS FOR KM: 
            # 1. Why don't these match with coefficients from the archived code below?
            #     (Results starting line 763)

sum(newfit.glu$coefficients[-c(1:12)]) # sums to psi_1 
#[1] 0.9605

coef(qc.fit.glu)
#(intercept)        psi1 
#    69.4000      0.9605 
  

## qgcomp.glm.boot (bootstrap CI)
qcboot.fit.glu <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                age_mc1 + #covars
                                ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                meducat_bl_2 + meducat_bl_3 +
                                marcat_mc1 +
                                agwork_mc1_2 + agwork_mc1_3 +
                                povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                log_tot +
                                  log_ace + log_mala + log_bens + #exposures
                                  log_perm +
                                  log_meth +
                                  log_imid +
                                  log_mn +
                                  log_gly,
                               expnms=Xnm.all, #exposure names
                               data.glu, #dataset
                               family=gaussian(), 
                               q=4, # quantiles
                               B=10, # Bootstrap iterations, should be 200-500+ in practice
                               seed=125)

qcboot.fit.glu ## This only gives the total exposure effect, no scaled effect size for each component of the mixture
#            Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
#(Intercept)    91.73       2.05    87.70    95.76   44.65   <2e-16
#psi1            0.96       1.78    -2.54     4.46    0.54     0.59

## KM's code:
# underlying model of qgcomp fit
qcboot.fit.glu$qx$Y.cc.glu <- qcboot.fit.glu$fit$data$Y.cc.glu # bring outcome back into quantized data

#add the covariates too
qcboot.fit.glu$qx$age_mc1 <- qcboot.fit.glu$fit$data$age_mc1 
qcboot.fit.glu$qx$ageusa_cat_bl_2 <- qcboot.fit.glu$fit$data$ageusa_cat_bl_2 
qcboot.fit.glu$qx$ageusa_cat_bl_3 <- qcboot.fit.glu$fit$data$ageusa_cat_bl_3 
qcboot.fit.glu$qx$meducat_bl_2 <- qcboot.fit.glu$fit$data$meducat_bl_2 
qcboot.fit.glu$qx$meducat_bl_3 <- qcboot.fit.glu$fit$data$meducat_bl_3 
qcboot.fit.glu$qx$marcat_mc1 <- qcboot.fit.glu$fit$data$marcat_mc1 
qcboot.fit.glu$qx$agwork_mc1_2 <- qcboot.fit.glu$fit$data$agwork_mc1_2 
qcboot.fit.glu$qx$agwork_mc1_3 <- qcboot.fit.glu$fit$data$agwork_mc1_3 
qcboot.fit.glu$qx$povcat_rc_mc1_2 <- qcboot.fit.glu$fit$data$povcat_rc_mc1_2 
qcboot.fit.glu$qx$povcat_rc_mc1_3 <- qcboot.fit.glu$fit$data$povcat_rc_mc1_3 
qcboot.fit.glu$qx$log_tot <- qcboot.fit.glu$fit$data$log_tot 

# adjusted coefficients from underlying model
# absolute change in outcome by increasing each exposure by 1
newboot.fit.glu <- lm(Y.cc.glu ~ age_mc1 + #covars
                   ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                   meducat_bl_2 + meducat_bl_3 +
                   marcat_mc1 +
                   agwork_mc1_2 + agwork_mc1_3 +
                   povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                   log_tot +
                   log_ace_q + log_mala_q + log_bens_q + #exposures
                   log_perm_q +
                   log_meth_q +
                   log_imid_q +
                   log_mn_q +
                   log_gly_q, data = qcboot.fit.glu$qx)

newboot.fit.glu ##this is where the individual coefficients are
tidy(newboot.fit.glu, conf.int = T) # and confidence intervals

        ## QUESTIONS FOR KM: 
            # 2. Why don't these match with coefficients from the archived code below?
            #     (Results starting line 763)
            # 3. Why are these exactly the same as the no-bootstrapping version above?

sum(newboot.fit.glu$coefficients[-c(1:12)]) # sums to psi_1 
#[1] 0.9604836

coef(qcboot.fit.glu)
#(intercept)        psi1 
#91.7310955   0.9604836



#### Binary outcomes (start with diabetes)

# conditional odds ratio
qc.fit.diabe <- qgcomp.glm.noboot(Y.cc.diabe ~ #outcome
                                    age_mc1 + #covars
                                    ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                    meducat_bl_2 + meducat_bl_3 +
                                    marcat_mc1 +
                                    agwork_mc1_2 + agwork_mc1_3 +
                                    povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                    log_tot +
                                    log_ace + log_mala + log_bens + #exposures
                                    log_perm +
                                    log_meth +
                                    log_imid +
                                    log_mn +
                                    log_gly,
                                  expnms=Xnm.all, #exposure names
                                  data.diabe, #dataset
                                  family=binomial(), q = 4)
qc.fit.diabe

# underlying model of qgcomp fit
qc.fit.diabe$qx$Y.cc.diabe <- qc.fit.diabe$fit$data$Y.cc.diabe # bring outcome back into quantized data

#add the covariates too
qc.fit.diabe$qx$age_mc1 <- qc.fit.diabe$fit$data$age_mc1 
qc.fit.diabe$qx$ageusa_cat_bl_2 <- qc.fit.diabe$fit$data$ageusa_cat_bl_2 
qc.fit.diabe$qx$ageusa_cat_bl_3 <- qc.fit.diabe$fit$data$ageusa_cat_bl_3 
qc.fit.diabe$qx$meducat_bl_2 <- qc.fit.diabe$fit$data$meducat_bl_2 
qc.fit.diabe$qx$meducat_bl_3 <- qc.fit.diabe$fit$data$meducat_bl_3 
qc.fit.diabe$qx$marcat_mc1 <- qc.fit.diabe$fit$data$marcat_mc1 
qc.fit.diabe$qx$agwork_mc1_2 <- qc.fit.diabe$fit$data$agwork_mc1_2 
qc.fit.diabe$qx$agwork_mc1_3 <- qc.fit.diabe$fit$data$agwork_mc1_3 
qc.fit.diabe$qx$povcat_rc_mc1_2 <- qc.fit.diabe$fit$data$povcat_rc_mc1_2 
qc.fit.diabe$qx$povcat_rc_mc1_3 <- qc.fit.diabe$fit$data$povcat_rc_mc1_3 
qc.fit.diabe$qx$log_tot <- qc.fit.diabe$fit$data$log_tot 

# adjusted coefficients from underlying model
# absolute change in outcome by increasing each exposure by 1
# can also use qc.fit2$fit
newfit.diabe <- glm(Y.cc.diabe ~ age_mc1 + #covars
                      ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                      meducat_bl_2 + meducat_bl_3 +
                      marcat_mc1 +
                      agwork_mc1_2 + agwork_mc1_3 +
                      povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                      log_tot +
                      log_ace_q + log_mala_q + log_bens_q + #exposures
                      log_perm_q +
                      log_meth_q +
                      log_imid_q +
                      log_mn_q +
                      log_gly_q, 
                    data = qc.fit.diabe$qx, family = binomial())

tidy(newfit.diabe, conf.int = T) #estimates and CIs for individual exposures

sum(newfit.diabe$coefficients[-c(1:12)]) # sums to psi_1 
#[1] 0.2106

coef(qc.fit.diabe)
#(intercept)  psi1 
#-4.9032      0.2106 




# marginal odds ratio

qcboot.fit.diabe <- qgcomp.glm.boot(Y.cc.diabe ~ #outcome
                                    age_mc1 + #covars
                                    ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                    meducat_bl_2 + meducat_bl_3 +
                                    marcat_mc1 +
                                    agwork_mc1_2 + agwork_mc1_3 +
                                    povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                    log_tot +
                                    log_ace + log_mala + log_bens + #exposures
                                    log_perm +
                                    log_meth +
                                    log_imid +
                                    log_mn +
                                    log_gly,
                                  expnms=Xnm.all, #exposure names
                                  data.diabe, #dataset
                                  family=binomial(), 
                                  q=4, # quantiles
                                  B=10, # Bootstrap iterations, should be 200-500+ in practice
                                  seed=125, rr = FALSE)

qcboot.fit.diabe # total scaled effect

# individual exposure beta coefficients and CIs
tidy(qcboot.fit.diabe$fit, conf.int = T) 

# underlying model of qgcomp fit with bootstrapping
qcboot.fit.diabe$qx$Y.cc.diabe <- qcboot.fit.diabe$fit$data$Y.cc.diabe # bring outcome back into quantized data

#add the covariates too
qcboot.fit.diabe$qx$age_mc1 <- qcboot.fit.diabe$fit$data$age_mc1 
qcboot.fit.diabe$qx$ageusa_cat_bl_2 <- qcboot.fit.diabe$fit$data$ageusa_cat_bl_2 
qcboot.fit.diabe$qx$ageusa_cat_bl_3 <- qcboot.fit.diabe$fit$data$ageusa_cat_bl_3 
qcboot.fit.diabe$qx$meducat_bl_2 <- qcboot.fit.diabe$fit$data$meducat_bl_2 
qcboot.fit.diabe$qx$meducat_bl_3 <- qcboot.fit.diabe$fit$data$meducat_bl_3 
qcboot.fit.diabe$qx$marcat_mc1 <- qcboot.fit.diabe$fit$data$marcat_mc1 
qcboot.fit.diabe$qx$agwork_mc1_2 <- qcboot.fit.diabe$fit$data$agwork_mc1_2 
qcboot.fit.diabe$qx$agwork_mc1_3 <- qcboot.fit.diabe$fit$data$agwork_mc1_3 
qcboot.fit.diabe$qx$povcat_rc_mc1_2 <- qcboot.fit.diabe$fit$data$povcat_rc_mc1_2 
qcboot.fit.diabe$qx$povcat_rc_mc1_3 <- qcboot.fit.diabe$fit$data$povcat_rc_mc1_3 
qcboot.fit.diabe$qx$log_tot <- qcboot.fit.diabe$fit$data$log_tot 

# adjusted coefficients from underlying model
# absolute change in outcome by increasing each exposure by 1
newbootfit.diabe <- glm(Y.cc.diabe ~ age_mc1 + #covars
                          ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                          meducat_bl_2 + meducat_bl_3 +
                          marcat_mc1 +
                          agwork_mc1_2 + agwork_mc1_3 +
                          povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                          log_tot +
                          log_ace_q + log_mala_q + log_bens_q + #exposures
                          log_perm_q +
                          log_meth_q +
                          log_imid_q +
                          log_mn_q +
                          log_gly_q, data = qcboot.fit.diabe$qx, family = binomial())

newbootfit.diabe

boot_func <- function(data, indices) {
  # Resample data
  boot_data <- data[indices, ]
  
  # Fit the glm model on the resampled data
  fit <- glm(Y.cc.diabe ~ age_mc1 + #covars
               ageusa_cat_bl_2 + ageusa_cat_bl_3 +
               meducat_bl_2 + meducat_bl_3 +
               marcat_mc1 +
               agwork_mc1_2 + agwork_mc1_3 +
               povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
               log_tot +
               log_ace_q + log_mala_q + log_bens_q + #exposures
               log_perm_q +
               log_meth_q +
               log_imid_q +
               log_mn_q +
               log_gly_q, data = boot_data, family = binomial())
  
  # Return the coefficients
  return(coef(fit))
}

set.seed(123)  # For reproducibility
boot_results <- boot(
  data = qcboot.fit.diabe$qx,  # Your dataset
  statistic = boot_func,  # The bootstrap function
  R = 1000                # Number of bootstrap replicates
)

coef_names <- names(coef(newbootfit.diabe))

# Calculate 95% bootstrap confidence intervals for each coefficient
ci <- t(apply(boot_results$t, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
}))

# Add coefficient names for clarity
rownames(ci) <- coef_names
colnames(ci) <- c("Lower 95% CI", "Upper 95% CI")

# Mean of bootstrap coefficients
boot_means <- colMeans(boot_results$t)

# Combine results
bootstrap_summary <- data.frame(
  Coefficient = coef_names,
  Estimate = boot_means,
  `Lower 95% CI` = ci[, 1],
  `Upper 95% CI` = ci[, 2]
)

print(bootstrap_summary)
sum(bootstrap_summary$Estimate[-c(1:12)])
#[1] 0.1568298

sum(newbootfit.diabe$coefficients[-c(1:12)]) # sums to psi_1 
#[1] 0.2105957

coef(qcboot.fit.diabe)
#(intercept)        psi1 
#-0.7952465   0.2002556

          ## QUESTIONS FOR KM: 
            # 4. These don't match (0.21 and 0.20). It looks like the newbootfit
            #     is estimating the conditional OR. What is going wrong?
      



### ARCHIVE:


## Jackie's suggestion for pulling our individual estimates

## First specifying all exposures (Xnm.all) 
qcboot.fit.glu <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.all, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)

p.glu = plot(qcboot.fit.glu) ##assess linearity of total exposure effect
#plot(qcboot.fit.glu, pointwiseref = 1) ##change referent category for pointwise comparisons
#pointwisebound.boot(qcboot.fit.glu, pointwiseref = 1) ##pointwise comparison confidence intervals
#qgcomp:::modelbound.boot(qcboot.fit.glu) ##confidence intervals for regression line

## Now trying to pull out estimates of individual components (Xnm.ace, Xnm.mala, etc.)

Xnm.ace <-c("log_ace") ##individual exposures
Xnm.mala <-c("log_mala")
Xnm.bens <-c("log_bens")
Xnm.perm <-c("log_perm")
Xnm.meth <-c("log_meth")
Xnm.imid <-c("log_imid")
Xnm.mn <-c("log_mn")
Xnm.gly <-c("log_gly")

### This is what we DONT want to do because it's not looking at the quantiles of the other exposures
qcboot.fit.glu.ace.original <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.ace, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)

test_ace_dat_original <- data.frame(qcboot.fit.glu.ace.original$qx)
test_ace_dat_original <- test_ace_dat_original %>% rename(log_ace_q = qcboot.fit.glu.ace.original.qx)
# underlying model of qgcomp fit
test_ace_dat_original$Y.cc.glu <- qcboot.fit.glu.ace.original$fit$data$Y.cc.glu # bring outcome back into quantized data

#add the covariates too
test_ace_dat_original$age_mc1 <- qcboot.fit.glu.ace.original$fit$data$age_mc1 
test_ace_dat_original$ageusa_cat_bl_2 <- qcboot.fit.glu.ace.original$fit$data$ageusa_cat_bl_2 
test_ace_dat_original$ageusa_cat_bl_3 <- qcboot.fit.glu.ace.original$fit$data$ageusa_cat_bl_3 
test_ace_dat_original$meducat_bl_2 <- qcboot.fit.glu.ace.original$fit$data$meducat_bl_2 
test_ace_dat_original$meducat_bl_3 <- qcboot.fit.glu.ace.original$fit$data$meducat_bl_3 
test_ace_dat_original$marcat_mc1 <- qcboot.fit.glu.ace.original$fit$data$marcat_mc1 
test_ace_dat_original$agwork_mc1_2 <- qcboot.fit.glu.ace.original$fit$data$agwork_mc1_2 
test_ace_dat_original$agwork_mc1_3 <- qcboot.fit.glu.ace.original$fit$data$agwork_mc1_3 
test_ace_dat_original$povcat_rc_mc1_2 <- qcboot.fit.glu.ace.original$fit$data$povcat_rc_mc1_2 
test_ace_dat_original$povcat_rc_mc1_3 <- qcboot.fit.glu.ace.original$fit$data$povcat_rc_mc1_3 
test_ace_dat_original$log_tot <- qcboot.fit.glu.ace.original$fit$data$log_tot 

# and other exposures
test_ace_dat_original$log_mala <- qcboot.fit.glu.ace.original$fit$data$log_mala
test_ace_dat_original$log_bens <- qcboot.fit.glu.ace.original$fit$data$log_bens
test_ace_dat_original$log_perm <- qcboot.fit.glu.ace.original$fit$data$log_perm
test_ace_dat_original$log_meth <- qcboot.fit.glu.ace.original$fit$data$log_meth
test_ace_dat_original$log_imid <- qcboot.fit.glu.ace.original$fit$data$log_imid
test_ace_dat_original$log_mn <- qcboot.fit.glu.ace.original$fit$data$log_mn
test_ace_dat_original$log_gly <- qcboot.fit.glu.ace.original$fit$data$log_gly

# adjusted coefficients from underlying model
# absolute change in outcome by increasing each exposure by 1
test_ace_fit_original <- lm(Y.cc.glu ~ #outcome
                     age_mc1 + #covars
                     ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                     meducat_bl_2 + meducat_bl_3 +
                     marcat_mc1 +
                     agwork_mc1_2 + agwork_mc1_3 +
                     povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                     log_tot +
                     log_ace_q + log_mala + log_bens + #exposures
                     log_perm +
                     log_meth +
                     log_imid +
                     log_mn +
                     log_gly, data = test_ace_dat_original)

test_ace_fit_original ##this is where the individual coefficients are
tidy(test_ace_fit_original, conf.int = T) # and confidence intervals


## This IS what we want to do to match the research question in line 247
# find out how increasing exposure (acetate) by 1 quantile impacts glucose outcome while holding **quantiles** of other covariates constant
qc.fit.glu$qx$log_ace <- data.glu$log_ace

qcboot.fit.glu.ace <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                    age_mc1 + #covars
                                    ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                    meducat_bl_2 + meducat_bl_3 +
                                    marcat_mc1 +
                                    agwork_mc1_2 + agwork_mc1_3 +
                                    povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                    log_tot +
                                    log_ace + log_mala_q + log_bens_q + #exposures
                                    log_perm_q +
                                    log_meth_q +
                                    log_imid_q +
                                    log_mn_q +
                                    log_gly_q,
                                  expnms=Xnm.ace, #exposure names
                                  qc.fit.glu$qx, #dataset
                                  family=gaussian(), 
                                  q=4, # quantiles
                                  B=10, # Bootstrap iterations, should be 200-500+ in practice
                                  seed=125)
qcboot.fit.glu.ace


test_ace_dat <- data.frame(qcboot.fit.glu.ace$qx)
test_ace_dat <- test_ace_dat %>% rename(log_ace_q = qcboot.fit.glu.ace.qx)
# underlying model of qgcomp fit
test_ace_dat$Y.cc.glu <- qcboot.fit.glu.ace$fit$data$Y.cc.glu # bring outcome back into quantized data

#add the covariates too
test_ace_dat$age_mc1 <- qcboot.fit.glu.ace$fit$data$age_mc1 
test_ace_dat$ageusa_cat_bl_2 <- qcboot.fit.glu.ace$fit$data$ageusa_cat_bl_2 
test_ace_dat$ageusa_cat_bl_3 <- qcboot.fit.glu.ace$fit$data$ageusa_cat_bl_3 
test_ace_dat$meducat_bl_2 <- qcboot.fit.glu.ace$fit$data$meducat_bl_2 
test_ace_dat$meducat_bl_3 <- qcboot.fit.glu.ace$fit$data$meducat_bl_3 
test_ace_dat$marcat_mc1 <- qcboot.fit.glu.ace$fit$data$marcat_mc1 
test_ace_dat$agwork_mc1_2 <- qcboot.fit.glu.ace$fit$data$agwork_mc1_2 
test_ace_dat$agwork_mc1_3 <- qcboot.fit.glu.ace$fit$data$agwork_mc1_3 
test_ace_dat$povcat_rc_mc1_2 <- qcboot.fit.glu.ace$fit$data$povcat_rc_mc1_2 
test_ace_dat$povcat_rc_mc1_3 <- qcboot.fit.glu.ace$fit$data$povcat_rc_mc1_3 
test_ace_dat$log_tot <- qcboot.fit.glu.ace$fit$data$log_tot 

# and other exposures
test_ace_dat$log_mala_q <- qcboot.fit.glu.ace$fit$data$log_mala_q
test_ace_dat$log_bens_q <- qcboot.fit.glu.ace$fit$data$log_bens_q
test_ace_dat$log_perm_q <- qcboot.fit.glu.ace$fit$data$log_perm_q
test_ace_dat$log_meth_q <- qcboot.fit.glu.ace$fit$data$log_meth_q
test_ace_dat$log_imid_q <- qcboot.fit.glu.ace$fit$data$log_imid_q
test_ace_dat$log_mn_q <- qcboot.fit.glu.ace$fit$data$log_mn_q
test_ace_dat$log_gly_q <- qcboot.fit.glu.ace$fit$data$log_gly_q

# adjusted coefficients from underlying model
# absolute change in outcome by increasing each exposure by 1
test_ace_fit <- lm(Y.cc.glu ~ #outcome
                   age_mc1 + #covars
                   ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                   meducat_bl_2 + meducat_bl_3 +
                   marcat_mc1 +
                   agwork_mc1_2 + agwork_mc1_3 +
                   povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                   log_tot +
                   log_ace_q + log_mala_q + log_bens_q + #exposures
                   log_perm_q +
                   log_meth_q +
                   log_imid_q +
                   log_mn_q +
                   log_gly_q, data = test_ace_dat)

test_ace_fit ##this is where the individual coefficients are
tidy(test_ace_fit, conf.int = T) # and confidence intervals


qcboot.fit.glu.mala <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.mala, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)

qcboot.fit.glu.bens <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.bens, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)

qcboot.fit.glu.perm <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot + 
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.perm, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)

qcboot.fit.glu.meth <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.meth, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)

qcboot.fit.glu.imid <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.imid, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)

qcboot.fit.glu.mn <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.mn, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)

qcboot.fit.glu.gly <- qgcomp.glm.boot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                                      expnms=Xnm.gly, #exposure names
                                      data.glu, #dataset
                                      family=gaussian(), 
                                      q=4, # quantiles
                                      B=10, # Bootstrap iterations, should be 200-500+ in practice
                                      seed=125)
qcboot.fit.glu.ace.original
#            Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
#(Intercept)   92.759      1.769    89.29    96.23   52.44   <2e-16
#psi1          -0.281      1.483    -3.19     2.63   -0.19     0.85

qcboot.fit.glu.mala
#            Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
#(Intercept)   91.099      1.975   87.229    94.97   46.14   <2e-16
#psi1           0.966      0.855   -0.709     2.64    1.13     0.26

qcboot.fit.glu.bens
#            Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
#(Intercept)   93.241      3.309    86.76    99.73    28.2   <2e-16
#psi1          -0.462      1.535    -3.47     2.55    -0.3     0.76

qcboot.fit.glu.perm
#            Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
#(Intercept)   93.595      1.939    89.80   97.395   48.28   <2e-16
#psi1          -0.698      0.828    -2.32    0.924   -0.84      0.4

qcboot.fit.glu.meth
#            Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
#(Intercept)    95.09       3.97    87.30   102.88   23.93   <2e-16
#psi1           -1.69       2.21    -6.03     2.64   -0.77     0.44

qcboot.fit.glu.imid
#            Estimate Std. Error Lower CI  Upper CI t value Pr(>|t|)
#(Intercept)    99.41       5.06    89.49  109.318   19.66   <2e-16
#psi1           -4.57       2.56    -9.59    0.445   -1.79    0.075

qcboot.fit.glu.mn
#            Estimate Std. Error Lower CI  Upper CI t value Pr(>|t|)
#(Intercept)    97.04       3.55    90.08  103.990   27.35   <2e-16
#psi1           -2.99       1.68    -6.28    0.298   -1.78    0.076

qcboot.fit.glu.gly
#            Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
#(Intercept)   91.108      1.340   88.482    93.73   68.00   <2e-16
#psi1           0.960      0.831   -0.669     2.59    1.15     0.25

## See how these estimates compare to the weights/"scaled effect size" from qgcomp.glm.noboot
qcnoboot.fit.glu <- qgcomp.glm.noboot(Y.cc.glu ~ #outcome
                                        age_mc1 + #covars
                                        ageusa_cat_bl_2 + ageusa_cat_bl_3 +
                                        meducat_bl_2 + meducat_bl_3 +
                                        marcat_mc1 +
                                        agwork_mc1_2 + agwork_mc1_3 +
                                        povcat_rc_mc1_2 + povcat_rc_mc1_3 + 
                                        log_tot +
                                        log_ace + log_mala + log_bens + #exposures
                                        log_perm +
                                        log_meth +
                                        log_imid +
                                        log_mn +
                                        log_gly,
                             expnms=Xnm.all,
                             data.glu, family=gaussian(), q=4)
qcnoboot.fit.glu
#Scaled effect size (positive direction, sum of positive coefficients = 4.9)
#log_mala log_bens  log_gly  log_ace log_perm log_meth 
#0.2923   0.2140   0.2015   0.1702   0.0809   0.0410 
#
#Scaled effect size (negative direction, sum of negative coefficients = -3.94)
#log_imid   log_mn 
#0.798    0.202 
#
#Mixture slope parameters (Delta method CI):
#  
#  Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
#(Intercept)    69.40      11.44     47.0    91.83    6.06  4.3e-09
#psi1            0.96       1.61     -2.2     4.12    0.60     0.55

plot(qcnoboot.fit.glu)

    ##The directionality of the estimates from pulling out individual exposures 
    ##from the qgcomp.glm.boot command only somewhat match the direction of 
    ##"scaled effect sizes" from the qgcomp.glm.noboot command.

## Testing gfoRmula package
test <- gformula_continuous_eof(
  covnames = c("log_ace","log_mala","log_bens",
                "log_perm",
                "log_meth",
                "log_imid",
                "log_mn",
                "log_gly"),  
  outcome_name = "Y.cc.glu",     
  covtypes = c("continuous", "continuous", "continuous", "continuous", 
               "continuous", "continuous", "continuous", "continuous"), 
  covparams <- list(NA, NA, NA, NA, NA, NA, NA, NA),
  outcome_type = "continuous",         
  data = data.glu,
  nsimul = nrow(data.glu)
)

## Getting errors. Kelsey is going to try to do this manually.
