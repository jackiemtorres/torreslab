rm(list=ls())

setwd("Library/CloudStorage/Box-Box/CHAM Maternal Cognition Study/Preliminary Results/Data for Jackie 12-19-22")
getwd()

# COVID stressors and cardiometabolic outcomes
library(mice)
library(dplyr)
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(msm)
library(readxl)
library(gridExtra)
library(tidyverse)

#Importing dataset created in STATA
MatCog_Stress <- read_excel("Library/CloudStorage/Box-Box/CHAM Maternal Cognition Study/Preliminary Results/Data for Jackie 12-19-22/MatCog_Stress_for_R.xls")

# --------------------------------------------------------------------------------------------------------------------------
# Multiple imputation procedures
# --------------------------------------------------------------------------------------------------------------------------

#Imputing with mice
md.pattern(MatCog_Stress)

MatCog <- mice(MatCog_Stress,m=5,maxit=50,meth='pmm',seed=500)
summary(MatCog)
densityplot(MatCog)

# --------------------------------------------------------------------------------------------------------------------------
# Correlation Matrices - stressor variables only (note: come back and do with imputed dataset)
# --------------------------------------------------------------------------------------------------------------------------

#Subsetting to stressors to get correlation matrix
MatCog_Stress_Only <- dplyr::select(MatCog_Stress, born_in_mexico, aces_binary_merged1, str_care_mcog1:high_worry_deportation)

MatCog_Stress_Only <- MatCog_Stress_Only %>%
  na.omit() #not needed if replacing this with an imputed dataset

# Corr plot = stressors
M1 = cor(MatCog_Stress_Only)

#Subsetting to those born in Mexico
MatCog_Stress_Mex <- dplyr::filter(MatCog_Stress_Only, born_in_mexico==1)
MatCog_Stress_Mex <- dplyr::select(MatCog_Stress_Mex, -born_in_mexico)

# Corr plot = stressors
M2 = cor(MatCog_Stress_Mex)

#Combine into two-paneled graph
par(mfrow=c(1,2))

corplot_all <- corrplot(M1, method = 'circle', tl.cex = .5, tl.col = "black")
corplot_mex <- corrplot(M2, method = 'circle', tl.cex = .5, tl.col = "black")

# --------------------------------------------------------------------------------------------------------------------------
# analyses: individual and combined
# --------------------------------------------------------------------------------------------------------------------------

stressor_measures <- exposures <- c("str_care_mcog1","str_cov_mcog1","str_else_mcog1","str_fin_mcog1",
        "str_hlthp_mcog1", "str_hlths_mcog1", "str_job_mcog1", "str_natdis_mcog1", "str_pare_mcog1", 
        "str_rel_mcog1", "high_worry_deportation")

# --------------------------------------------------------------------------------------------------------------------------
# estimate individual exposure propensity scores and record which people have propensity scores within the area of support
# --------------------------------------------------------------------------------------------------------------------------

pscore <- function(i, j) {
df_t <- mice::complete(MatCog, i)
  x <- predict(with(df_t, glm(get(stressor_measures[j]) ~  age_mcog1 + born_in_mexico + educat + qx_language_mcog1 + married + work_ag_3cat, family="binomial")), type="response")
 return(x)
}

include <- list()

for (j in 1:length(stressor_measures)) {
  
# calculate propensity scores for each imputed data set 
pscores_list <- lapply(1:5, function(x) pscore(x, j))
pscores <- do.call(cbind, pscores_list) 

# make matrix of exposure levels for each imputed data set 
observed_ls <- lapply(1:5, function(x) mice::complete(MatCog, x) %>% select(stressor_measures[j]))
observed_e <- do.call(cbind, observed_ls)

# find the minimum of propensity scores for people who were exposed
# and the maximum of propensity scores for people who were unexposed 
# to ensure there is adequate support 

pscores_support_ls <- lapply(1:5, function(x) cbind(min(pscores[,x][observed_e[,x]==1]), max(pscores[,x][observed_e[,x]==0])))
pscores_support <- do.call(rbind, pscores_support_ls)

# identify whether propensity scores are within the area of support for each imputed data set 
include_ls <- lapply(1:5, function(x) as.numeric(pscores[,x] >= pscores_support[x,1] & pscores[,x]<= pscores_support[x,2]))

include[[j]] <- do.call(cbind, include_ls)
}

names(include) <- stressor_measures

# --------------------------------------------------------------------------------------------------------------------------
# estimate association with birth weight for gestational age for individual exposures  
# --------------------------------------------------------------------------------------------------------------------------

gcomp_mi <- function(m, exposure) {

df_p <- mice::complete(MatCog, m)

# include those who fell within the area of support for each imputed data set 
df_p <- df_p %>% filter(include[[exposure]][,m]==1)

df_p$exposure <- df_p[[exposure]]

#COME BACK AND INCLUDE OUTCOME VAR
fit <- glm(cesd_score ~ exposure + age_mcog1 + born_in_mexico + educat + qx_language_mcog1 + married + work_ag_3cat, data=df_p)
nobs <- nobs(fit)

# calculate the ATE via just the beta coefficient -- this works because it's a linear model 
ate <- coef(fit)["exposure"] 

cov <- vcov(fit)
cov <- cov["exposure","exposure"]

results <- data.frame(cbind(ate = ate, var = cov, N=nobs)) 
return(results)
}

# --------------------------------------------------------------------------------------------------------------------------
# combine results across imputations using Rubin's combining rules and write results to file    
# --------------------------------------------------------------------------------------------------------------------------

all_results <- data.frame(matrix(ncol=5, nrow=length(exposures)))

for (i in 1:length(exposures)) {

# variance calculations from https://amstat.tandfonline.com/doi/pdf/10.1080/01621459.1996.10476908?needAccess=true
# var = average of complete data variances + variance of complete data means 
ate_M_ls <- lapply(1:5, function(x) gcomp_mi(x, exposures[i]))
ate_M_df <- data.frame(do.call(rbind, ate_M_ls))
names(ate_M_df) <- c("ate", "var", "N")

results <- data.frame(cbind(rd = mean(ate_M_df$ate), var = mean(ate_M_df$var) + ((30+1)/30)*var(ate_M_df$ate)))
results$lb <- results$rd - 1.96*sqrt(results$var)
results$ub <- results$rd + 1.96*sqrt(results$var)
results$N <- mean(ate_M_df$N)

all_results[i,] <- results 
}

rownames(all_results) <- exposures
colnames(all_results) <-  c("rd","var","lb","ub","N")

write.csv(all_results, file="stress_cesd_indiv.csv")

# --------------------------------------------------------------------------------------------------------------------------
# estimate mixture propensity scores -- propensity scores for each exposure while controlling for the others 
# --------------------------------------------------------------------------------------------------------------------------

# currently exposures are all stressors. decide if you want to focus on specific ones 

pscore_mix <- function(i, j) {
  df_t <- mice::complete(MatCog, i)
  exposure <- exposures[j]
  
  coexposures <- exposures[-(match(exposure, exposures))]
  
  # regress exposure of interest on other exposures plus covariates 
  formula <- paste(exposure, "~",   paste(coexposures, collapse="+ "), "+ age_mcog1 + born_in_mexico + educat + qx_language_mcog1 + married + work_ag_3cat")
  
  
  x <- predict(glm(formula=formula, data=df_t, family="binomial"), type="response")
  return(x)
}


include_mix <- list()

for (j in 1:length(exposures)) {
  
  # calculate propensity scores for each imputed data set 
  pscores_list <- lapply(1:5, function(x) pscore_mix(x, j))
  pscores_mix <- do.call(cbind, pscores_list) 
  
  # make matrix of exposure levels for each imputed data set 
  observed_ls <- lapply(1:5, function(x) mice::complete(MatCog, x) %>% select(stressor_measures[j]))
  observed_e <- do.call(cbind, observed_ls)
  
  # find the minimum of propensity scores for people who were exposed
  # and the maximum of propensity scores for people who were unexposed 
  # to ensure there is adequate support 
  
  pscores_support_ls <- lapply(1:5, function(x) cbind(min(pscores_mix[,x][observed_e[,x]==1]), max(pscores_mix[,x][observed_e[,x]==0])))
  pscores_support_mix <- do.call(rbind, pscores_support_ls)
  
  # identify whether propensity scores are within the area of support for each imputed data set 
  include_ls <- lapply(1:5, function(x) as.numeric(pscores_mix[,x] >= pscores_support_mix[x,1] & pscores_mix[,x]<= pscores_support_mix[x,2]))
  
  include_mix[[j]] <- do.call(cbind, include_ls)
}

names(include_mix) <- stressor_measures

# --------------------------------------------------------------------------------------------------------------------------
# estimate G-computation for CESD for individual exposures while controlling for the others 
# --------------------------------------------------------------------------------------------------------------------------

gcomp_mi_mix <- function(m, exposure) {
  
  df_p <- mice::complete(MatCog, m)
  
  df_p <- df_p %>% filter(include_mix[[exposure]][,m]==1)
  
  df_p$exposure <- df_p[[exposure]]
  
  coexposures <- exposures[-(match(exposure, exposures))]
  
  formula <- paste("cesd_score ~ exposure +",   paste(coexposures, collapse="+ "), "+ age_mcog1 + born_in_mexico + educat + qx_language_mcog1 + married + work_ag_3cat")
  
  fit <- glm(formula, data=df_p)
  nobs <- nobs(fit)
  
  # calculate the ATE via the beta coefficient -- this works because it's a linear model 
  ate <- coef(fit)["exposure"] 
  
  cov <- vcov(fit)
  cov <- cov["exposure","exposure"]
  
  results <- data.frame(cbind(ate = ate, var = cov, N=nobs)) 
  return(results)
}

# --------------------------------------------------------------------------------------------------------------------------
# combine results across imputations using Rubin's combining rules and write results to file  
# --------------------------------------------------------------------------------------------------------------------------

all_results_mix <- data.frame(matrix(ncol=5, nrow=length(exposures)))

for (i in 1:length(exposures)) {
  # variance calculations from https://amstat.tandfonline.com/doi/pdf/10.1080/01621459.1996.10476908?needAccess=true
  # var = average of complete data variances + variance of complete data means 
  ate_M_ls <- lapply(1:5, function(x) gcomp_mi_mix(x, exposures[i]))
  ate_M_df <- data.frame(do.call(rbind, ate_M_ls))
  names(ate_M_df) <- c("ate", "var", "N")
  
  results <- data.frame(cbind(rd = mean(ate_M_df$ate), var = mean(ate_M_df$var) + ((30+1)/30)*var(ate_M_df$ate)))
  results$lb <- results$rd - 1.96*sqrt(results$var)
  results$ub <- results$rd + 1.96*sqrt(results$var)
  results$N <- mean(ate_M_df$N)
  
  all_results_mix[i, ] <- results
}

rownames(all_results_mix) <- exposures 
colnames(all_results_mix) <- c("rd","var","lb","ub","N")

write.csv(all_results_mix, file="stress_cesd_mix.csv")


# --------------------------------------------------------------------------------------------------------------------------
# estimate G-computation for CESD for joint effect, all exposures at once
# --------------------------------------------------------------------------------------------------------------------------

gcomp_mi_joint <- function(m, ej) {
  
  df_p <- mice::complete(MatCog, m)
  
  df_p <- df_p %>% filter(include_mix[[exposures_joint[[ej]][1]]][,m]==1 & include_mix[[exposures_joint[[ej]][2]]][,m]==1)
  
  coexposures <- exposures[-(match(exposures_joint[[ej]], exposures))]
  
  formula <- paste("cesd_score ~ ", paste0(exposures_joint[[ej]], collapse="*"), " +", paste(coexposures, collapse=" + "), "+ age_mcog1 + born_in_mexico + educat + qx_language_mcog1 + married + work_ag_3cat")
  
  fit <- glm(formula, data=df_p)
  nobs <- nobs(fit)
  
  # calculate the ATE via sum of beta coefficients -- this works because it's a linear model 
  ate <- sum(coef(fit)[c(exposures_joint[[ej]], paste0(exposures_joint[[ej]], collapse=":"))])

  # calculate the SE via the delta method 
    
  cov <- vcov(fit)
  cov <- cov[c(exposures_joint[[ej]], paste0(exposures_joint[[ej]], collapse=":")),c(exposures_joint[[ej]], paste0(exposures_joint[[ej]], collapse=":"))]
  
  dm_var <- deltamethod(~ x1 + x2 + x3, 
              mean = coef(fit)[c(exposures_joint[[ej]], paste0(exposures_joint[[ej]], collapse=":"))], cov = cov, ses=T)^2
  
  results <- data.frame(cbind(ate = ate, var = dm_var, N=nobs)) 
  return(results)
}

# --------------------------------------------------------------------------------------------------------------------------
# combine results across imputations using Rubin's combining rules and write results to file  
# --------------------------------------------------------------------------------------------------------------------------
stressor_measures <- exposures <- c("str_care_mcog1", "str_fin_mcog1",
                                    "str_hlths_mcog1", "str_job_mcog1",  "high_worry_deportation")

exposures_joint <- list(c("str_care_mcog1","str_fin_mcog1"), 
                        c("str_care_mcog1","str_hlths_mcog1"), 
                        c("str_care_mcog1", "str_job_mcog1"), 
                        c("str_care_mcog1", "high_worry_deportation"), 
                        c("str_fin_mcog1", "str_care_mcog1"), 
                        c("str_fin_mcog1", "str_hlths_mcog1"), 
                        c("str_fin_mcog1", "str_job_mcog1"), 
                        c("str_fin_mcog1", "high_worry_deportation"), 
                        c("str_hlths_mcog1", "str_care_mcog1"), 
                        c("str_hlths_mcog1", "str_fin_mcog1"), 
                        c("str_hlths_mcog1","str_job_mcog1"),
                        c("str_hlths_mcog1","high_worry_deportation"),
                        c("str_job_mcog1","str_care_mcog1"),
                        c("str_job_mcog1","str_fin_mcog1"),
                        c("str_job_mcog1", "str_hlths_mcog1"), 
                        c("str_job_mcog1","high_worry_deportation"))

# variance calculations from https://amstat.tandfonline.com/doi/pdf/10.1080/01621459.1996.10476908?needAccess=true
# var = average of complete data variances + variance of complete data means 

all_results_joint_ls <- list()

for (ej in 1:length(exposures_joint)) {
  ate_M_ls <- lapply(1:5, function(x) gcomp_mi_joint(x, ej))
  ate_M_df <- data.frame(do.call(rbind, ate_M_ls))
  names(ate_M_df) <- c("ate", "var", "N")
  
  results <- data.frame(cbind(rd = mean(ate_M_df$ate), var = mean(ate_M_df$var) + ((30+1)/30)*var(ate_M_df$ate)))
  results$lb <- results$rd - 1.96*sqrt(results$var)
  results$ub <- results$rd + 1.96*sqrt(results$var)
  
  results$N <- mean(ate_M_df$N)

  all_results_joint_ls[[ej]] <- results 
}

all_results_joint <- data.frame(do.call(rbind, all_results_joint_ls)) 
rownames(all_results_joint) <- unlist(lapply(1:length(exposures_joint), function(x) paste0(exposures_joint[[x]], collapse=":")))

write.csv(all_results_joint, file="stress_cesd_joint.csv")

# --------------------------------------------------------------------------------------------------------------------------
# create a graphic to compare estimates
# --------------------------------------------------------------------------------------------------------------------------

#Individual
indiv <- read.csv("stress_cesd_indiv.csv")

ggplot(indiv, aes(x=X, y=rd, ymin=lb, ymax=ub, colour = X)) + 
  #specify position here
  geom_linerange(size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(size=3, shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Mean Difference", limits = c(-10, 10)) +
  #ggtitle ("Prevalence Ratio") +
  coord_flip() + theme_minimal() + theme(legend.title = element_blank())

#Mixed
mixed <- read.csv("stress_cesd_mix.csv")

ggplot(mixed, aes(x=X, y=rd, ymin=lb, ymax=ub, colour = X)) + 
  #specify position here
  geom_linerange(size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(size=3, shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Mean Difference", limits = c(-10, 10)) +
  #ggtitle ("Prevalence Ratio") +
  coord_flip() + theme_minimal() + theme(legend.title = element_blank())

#Combine plots
# add a group column to both datasets
indiv$model_type <- "Individual"
mixed$model_type <- "Mixed"
 
dotCOLS = c("#a6d8f0","#f9b282")
barCOLS = c("#008fd5","#de6b35")

# combine the two datasets                      
combined = rbind(indiv, mixed)

#Combined Graphic
ggplot(combined, aes(x=X, y=rd, ymin=lb, ymax=ub, col=model_type, fill=model_type)) + 
  #specify position here
  geom_linerange(size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(size=3, shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  scale_x_discrete(name="") +
  scale_y_continuous(name="Mean Difference", limits = c(-10, 10)) +
  #ggtitle ("Prevalence Ratio") +
  coord_flip() + theme_minimal() + theme(legend.title = element_blank())
