################################################################################################
# NAME:         Household Crowding and Neurocognitive Performance in a Cohort of Middle-Aged Latina Women: The CHAMACOS Maternal Cognition Study
# AUTHOR:       Kelsey MacCuish
# CREATED:      12/06/2024
# CONTENTS:     Table 1. Overall Descriptive Characteristics, CHAMACOS Maternal Cognition Study (lines 507-570)
#               Table 2. Mean Difference and 95% Confidence Intervals, Household Crowding and Neurocognitive Z-Scores (lines 1424-1837)
#               Table 3. Mean Difference and 95% Confidence Intervals, Household Crowding and Candidate Mediators (lines 1078-1420)
#               Table 4. Coefficients and 95% Confidence Intervals, Household Crowding and Neurocognitive Z-Scores with Sleep Duration (lines 2023-2290)
#               Figure 2. Coefficients and 95% Confidence Intervals for Natural Direct and Indirect Estimates of the Association between Household Crowding and Neurocognitive Z-Scores, Mediated via Nightly Hours of Sleep (lines 2294-2422)
#               eTable 1. Comparison of Baseline Characteristics Overall and by Maternal Cognition Study Participants (lines 3577-3714)
#               eTable 5. Descriptive Characteristics of the Analytic Sample Stratified by Household Crowding Status (People per Bedroom) (lines 576-597)
#               eTable 6. Descriptive Characteristics of the Analytic Sample Stratified by Household Crowding Status (People per Room) (lines 601-622)
#               eTable 7. Sensitivity Analysis with Path-Specific Effects (lines 2429-2800)
#               eTable 8. Sensitivity Analysis - Housing Density Calculated as Number of Rooms Minus 1 (lines 2804-2990)
#               eTable 9. Sensitivity Analysis - Excluding Outliers of Nightly Hours of Sleep (lines 2994-3378)
#               eTable 10. Sensitivity Analysis - Using Alternative Cognitive Scores (lines 3382-3573)
################################################################################################

## Outcomes: ## 
# primary:
# cognitive performance (memory, executive function, verbal fluency, global composition)

# secondary:
# sleep duration, quality, naps (slpwkend_mc1, slpwkday_mc1, slpqual_mc1, napfq_mc1)
# depression (depress_mc1)
# anxiety (gadscore_mc1)

## Exposures: ##
# primary:
# household crowding
# number of people per room (density) (density1_mc1)
# number of people per bedroom (ppb_mc1)
# <=2 vs. >2 people per bedroom (ppb2_mc1)

# secondary exposure: 
# total number of people in household (lvhome_n_mc1)
# number of children < 18 in household (lvhome_n18_mc1)

## Covariates ##
# age
# nativity
# language of assessment 
# educational attainment
# occupational status (mat cog)
# marital status (mat cog)
# income/poverty (mat cog)


#### Load libraries ####
library(haven)
library(dplyr)
library(tableone)
library(broom)
library(ggplot2)
library(tidyr)
library(lme4)
library(ggdag)
library(broom.mixed)
library(ipw)
library(geepack)
library(nnet)
library(multcomp)
library(RColorBrewer)
library(tidycensus)
library(usmap)
library(tigris)
library(boot)
library(marginaleffects)
library(forcats)
library(mediation)
library(devtools)
library(lmtp)
library(SuperLearner)
library(medoutcon)
#remotes::install_github("tlverse/tmle3")
#devtools::install_github("nt-williams/HDmediation/HDmediation")
#devtools::install_github("nt-williams/HDmediation@mlr3superlearner", subdir = "HDmediation")
#remotes::install_github("nt-williams/crumble")
library(tmle3)
library(sl3)
library(Rsolnp)
library(crumble)
library(mma)
library(mlr3extralearners)
library(HDmediation)
library(hal9001)
library(arm)
library(ggsci)
library(earth)
library(viridis)
library(xgboost)
library(psych)
library(corrr)
library(data.table)
library(S7)
set.seed(124)



#### Read in data ####
cham_all <- read_dta("torres_stressor_01i_nohsn.dta") # read in CHAMACOS data

head(cham_all) # first 6 rows of dataset
str(cham_all) # variables and variable types
dim(cham_all) # number of observations and variables
names(cham_all) # names of variables



# nativity
# calculate age participant came back to the US at matcog visit 1
# yrsusa_matcog = yrsusa_dl + (date_matcog - dlvrydt)/365.25 
# yrsusa_matcog is years in the US at maternal cognition study visit 1
cham <- cham_all %>% mutate(age_came_back_to_us = age_qx_mc1 - yrsusa_matcog)

# compare ageusa_bl with calculated age got back to the US
cham %>% dplyr::select(age_qx_mc1, yrsusa_matcog, ageusa_bl, age_came_back_to_us) %>% filter(!is.na(age_qx_mc1))

# if ageusa_bl is -9, participant was born in the US
cham <- cham %>% mutate(ageusa_bl = ifelse(ageusa_bl == -9, 0,
                                           ageusa_bl))

# if ageusa_bl is < 0, set to 0
cham <- cham %>% mutate(ageusa_bl = ifelse(ageusa_bl < 0, 0, ageusa_bl))

cham %>% group_by(ageusa_bl) %>% count()

# if ageusa_bl is missing, and age_qx_mc1 = yrsusa_matcog, replace ageusa_bl with age_came_back_to_us
# if age at survey = number of years lived in the US, person has lived in the US their whole life
# i.e., they did not leave and come back to the US
cham <- cham %>% mutate(ageusa_bl = ifelse(is.na(ageusa_bl) & 
                                           age_qx_mc1 - yrsusa_matcog < 1, age_came_back_to_us, 
                                           ageusa_bl))

# check instances where there are discrepancies between ageusa_bl and calculated age_came_back_to_us
cham %>% filter(abs(ageusa_bl - age_came_back_to_us) > 2) %>% 
  dplyr::select(newid, age_qx_mc1, yrsusa_matcog, ageusa_bl, age_came_back_to_us, country_mom)

# check distribution of age came to the US variable
ggplot(data = cham, aes(x = ageusa_bl)) + 
  geom_histogram()

cham %>% filter(is.na(ageusa_bl)) %>% dplyr::select(newid, ageusa_bl, age_came_back_to_us,
                                                    country_mom, age_qx_mc1, yrsusa_matcog)

# set ageusa_bl to 0 for individuals with country_mom = 1 
cham <- cham %>% mutate(ageusa_bl = ifelse(is.na(ageusa_bl) & country_mom == 1, 0, ageusa_bl))


# categorize age came to the US variable 
cham <- cham %>% mutate(ageusa_cat = 
                          ifelse(country_mom == 1, "0",
                                 ifelse(ageusa_bl > 0 & ageusa_bl < 18, "<1-17", 
                                        ifelse(ageusa_bl >= 18 & ageusa_bl < 25, "18-24",
                                               ifelse(ageusa_bl >= 25, "25+", NA)))))

# check age came to the us distribution among matcog participatns 
cham %>% filter(!is.na(age_qx_mc1)) %>% group_by(ageusa_cat) %>% count()


# check those missing age that came to the US variable
cham %>% filter(is.na(ageusa_cat)) %>% dplyr::select(newid, country_mom, ageusa_bl, ageusa_cat,
                                                     age_came_back_to_us, newid, yrsusa_matcog, 
                                                     age_qx_mc1)

# 7 individuals with ageusa_bl = 0 and country_mom = 2
# recategorize these into "<1-17" group
cham <- cham %>% mutate(ageusa_cat = ifelse(country_mom == 2 & is.na(ageusa_cat), "<1-17", ageusa_cat))

# replace the two missing their ageusa_cat values with their calculated age_came_back_to_us values
cham <- cham %>% mutate(ageusa_cat = ifelse(is.na(ageusa_cat), age_came_back_to_us, ageusa_cat))

cham <- cham %>% mutate(ageusa_cat = ifelse(ageusa_cat == 17.3915138244629, "<1-17", 
                                            ifelse(ageusa_cat == 23.0636558532715, "18-24", 
                                                   ageusa_cat)))

cham %>% group_by(ageusa_cat) %>% count()

# create adulthood arrival to the US variable: born in US, <18 and 18+
cham <- cham %>% mutate(ageusa_18 = ifelse(ageusa_cat == "0", "0",
                                           ifelse(ageusa_cat == "<1-17", "<18", "18+")))

# create categorized language variable
# 1: English, 2: Spanish
cham <- cham %>% mutate(lang_exam_mc1 = ifelse(lang_exam_mc1 == 1, "English", "Spanish"))

# create overall parental education variable
# if don't know both parents education -> Don't know/Don't remember
# if don't know mom & dad never attended or don't know dad & mom never attended or neither attended -> None, never attended school
# all others -> any formal education
cham <- cham %>% mutate(parent_educ = 
                          ifelse(medu_mc1 == 9 &
                                   fedu_mc1 == 9,
                                 "Dont know/dont remember",
                                 ifelse(medu_mc1 == 9 &
                                          fedu_mc1 == 1
                                        | medu_mc1 == 1 & 
                                          fedu_mc1 == 9 | 
                                          medu_mc1 == 1 & 
                                          fedu_mc1 == 1, 
                                        "No formal education",
                                        "any formal educ")))

# dichotomize marriage variable to married or not married
# 1 or 2: married, otherwise not married
cham <- cham %>% mutate(married = 
                          ifelse(married_mc1 == 1 | 
                                   married_mc1 == 2, "married", "not married"))


# categorize worked since last visit variable
# 1 = worked since last visit, 0 = did not work since last visit
cham <- cham %>% mutate(work_mc1 = 
                          ifelse(work_mc1 == 1, "Worked since last visit",
                                 ifelse(work_mc1 == 0, "Did not work since last visit", NA)))


# categorize working now variable
# 1: working now, 0: not working now
cham <- cham %>% mutate(workc_mc1 = ifelse(workc_mc1 == 1, "Working now",
                                           ifelse(workc_mc1 == 0, "Not working now",
                                                  NA)))


# categorize mom's own education
# 1: <= 6th grade, 2: 7-12th grade, 3: >= HS graduate
cham <- cham %>% mutate(educcat_mom = ifelse(educcat_mom == 1, "<=6th grade", 
                                             ifelse(educcat_mom == 2, "7-11th grade", 
                                                    ifelse(educcat_mom == 3, ">= HS graduate", NA))))


# poverty status: povcat_mc1 
cham <- cham %>% mutate(povcat_mc1 = ifelse(povcat_mc1 == 1, "At or below poverty",
                                            ifelse(povcat_mc1 == 2, "200% FPL",
                                                   ifelse(povcat_mc1 == 3, ">=200% FPL", "Missing"))))
cham %>% group_by(povcat_mc1) %>% count()


# more than 2 people per bedroom: ppb2_mc1 
cham <- cham %>% mutate(ppb2_mc1 = ifelse(ppb2_mc1 == 0, "<=2", ">2"))


# sleep quality: slpqual_mc1 
cham <- cham %>% mutate(slpqual_mc1 = ifelse(slpqual_mc1 == 0, "very sound/restful", 
                                             ifelse(slpqual_mc1 == 1, "sound/restful", 
                                                    ifelse(slpqual_mc1 == 2, "average", 
                                                           ifelse(slpqual_mc1 == 3, "restless", "very restless")))))

# nap frequency: napfq_mc1 
cham <- cham %>% mutate(napfq_mc1 = ifelse(napfq_mc1 == 0, "0", 
                                           ifelse(napfq_mc1 == 1, "1 or 2", 
                                                  ifelse(napfq_mc1 == 2, "3 or 4", "5+"))))

# more than 1 person per room: ppr_1
cham <- cham %>% mutate(ppr_1 = ifelse(dens1_cat_mc1 == 1 | dens1_cat_mc1 == 2, "<=1", 
                                       ifelse(dens1_cat_mc1 == 3 | dens1_cat_mc1 == 4, ">1", NA)))

cham <- cham %>% mutate(ppr_1_update = ifelse(dens2_cat_mc1 == 1 | dens2_cat_mc1 == 2, "<=1", 
                                              ifelse(dens2_cat_mc1 == 3 | dens2_cat_mc1 == 4, ">1", NA)))


# number of residents <18 years of age
cham %>% filter(!is.na(age_qx_mc1)) %>% group_by(lvhome_n18_mc1) %>% count()


# impute newid 763's weekend sleep with the average difference between weekday sleep and weekend sleep
# average difference is 1.60 more hours of sleep on weekends than weekdays
cham <- cham %>% filter(!is.na(age_qx_mc1)) %>% mutate(sleep_diff = slpwkend_mc1 - slpwkday_mc1)
cham %>% dplyr::select(newid, slpwkend_mc1, slpwkday_mc1, sleep_diff)

cham %>% filter(!is.na(sleep_diff)) %>% filter(!is.na(age_qx_mc1)) %>% summarize(med = median(sleep_diff), 
                                                                                 mean = mean(sleep_diff))

cham <- cham %>% mutate(slpwkend_mc1 = ifelse(newid == 763, slpwkday_mc1 + 1.595583, slpwkend_mc1))


# create average sleep per night (weighted average between weeknight and weekend sleep)
cham <- cham %>% mutate(avg_sleep = ((slpwkday_mc1*5) + (slpwkend_mc1*2))/7)


# create binary indicator for whether person has < 6 hrs sleep/night
cham <- cham %>% mutate(less_6_hours = ifelse(avg_sleep < 6, "<6", ">=6"))


# create binary poverty variable
cham <- cham %>% mutate(pov_binary = ifelse(is.na(povcat_mc1) | 
                                              povcat_mc1 == "At or below poverty", 
                                            "At or below poverty", "Above poverty"))

# create restless sleep variable
cham <- cham %>% mutate(restless_sleep = ifelse(slpqual_mc1 == "restless" | slpqual_mc1 == "very restless", "restless", "restful"))

cham %>% group_by(restless_sleep) %>% count()

# create people per household variable
cham <- cham %>% mutate(people_per_household = density1_mc1 * rooms_mc1)

# filter out woman missing most baseline data
cham <- cham %>% filter(!newid == 314)


# filter to only include matcog participants
cham_tab1 <- cham %>% filter(!is.na(age_qx_mc1))


# calculate summary score of self-rated doctor diagnosed health conditions
# diabetes
# high blood pressure 
# heart problems
# cancer
# asthma 
# thyroid problems

# re-code don't knows to 999
cham_tab1 <- cham_tab1 %>% 
  mutate(diab_mc1 = ifelse(diab_mc1 == 9, 999, diab_mc1),
         hbp_mc1 = ifelse(hbp_mc1 == 9, 999, hbp_mc1),
         # chol_mc1 = ifelse(chol_mc1 == 9, 999, chol_mc1),
         heart_mc1 = ifelse(heart_mc1 == 9, 999, heart_mc1), 
         cancer_mc1 = ifelse(cancer_mc1 == 9, 999, cancer_mc1), 
         asth_mc1 = ifelse(asth_mc1 == 9, 999, asth_mc1), 
         thyr_mc1 = ifelse(thyr_mc1 == 9, 999, thyr_mc1)) 

# check distribution of each self-rated doctor diagnosed health condition
cham_tab1 %>% group_by(diab_mc1) %>% count() # only 1 person missing diabetes score
cham_tab1 %>% group_by(hbp_mc1) %>% count() # 3 missing HBP score
cham_tab1 %>% group_by(chol_mc1) %>% count() # 6 missing cholesterol
cham_tab1 %>% group_by(heart_mc1) %>% count() # 4 missing heart issues
cham_tab1 %>% group_by(cancer_mc1) %>% count() # 4 missing cancer
cham_tab1 %>% group_by(asth_mc1) %>% count() # 2 missing asthma
cham_tab1 %>% group_by(thyr_mc1) %>% count() # 0 missing thyroid
cham_tab1 %>% group_by(dep_mc1) %>% count() # 7 missing depression 
cham_tab1 %>% group_by(anx_mc1) %>% count() # 9 missing anxiety


# For those who don't know high blood pressure: 
# Check bpcat_mc1 variable for blood pressure information
# Replace the don't knows with 1 if have high blood pressure, 0 if no HBP
# only don't have hypertension if 0 for hbp_mc1 and bpcat_mc1 is 1 or 2 
cham_tab1 %>% filter(hbp_mc1 == 999)
cham_tab1 %>% filter(hbp_mc1 == 0 & bpcat_mc1 %in% c(3, 4))

recat_bp <- cham_tab1 %>% filter(hbp_mc1 == 0 & bpcat_mc1 %in% c(3, 4))

cham_tab1 <- cham_tab1 %>% mutate(hbp_mc1 = ifelse(newid %in% c(recat_bp$newid), 1, 
                                                   hbp_mc1))

cham_tab1 <- cham_tab1 %>% mutate(hbp_mc1 = ifelse(newid == 99 & bpcat_mc1 == 3, 1,
                                                   ifelse(newid == 102 & bpcat_mc1 == 1, 0,
                                                          ifelse(newid == 617 & bpcat_mc1 == 1, 0,
                                                                 hbp_mc1))))


# check that there are no more "don't know" responses for high blood pressure
cham_tab1 %>% group_by(hbp_mc1) %>% count()


# For those who don't know diabetes:
# Make sure they had fasting blood taken and check gluccat_mc1 variable for diabetes information
# Replace the don't knows with 1 if diabetes and 0 if no diabetes/pre-diabetes (if fasting)
# if not fasting, check hba1c score 

cham_tab1 %>% filter(diab_mc1 == 999)

# check: among those who self report no diabetes, how many do not have a fasting blood
# sugar score and have a hba1ccat_mc1 score that show diabetes
cham_tab1 %>% filter(diab_mc1 == 0 & fast_bld_mc1 == 0 & hba1ccat_mc1 == 4)

# check: among those who report no diabetes, how many do have a fasting blood sugar
# and blood sugar is in diabetes category
recat_diab <- cham_tab1 %>% filter(diab_mc1 == 0 & fast_bld_mc1 == 1 & gluccat_mc1 == 3)

cham_tab1 <- cham_tab1 %>% mutate(diab_mc1 = ifelse(newid %in% c(recat_diab$newid), 1, 
                                                    diab_mc1))

cham_tab1 <- cham_tab1 %>% mutate(diab_mc1 = 
                                    ifelse(newid == 222 & gluccat_mc1 == 2, 0, diab_mc1))

cham_tab1 %>% group_by(diab_mc1) %>% count()

cham_tab1 %>% filter(diab_mc1 == 0 & is.na(fast_bld_mc1))

# 4 missing heart issues - can't adjust b/c other heart variables are from 18 year and show 
# no heart issues - can't be sure they didn't develop heart issues between then and now 


# For those who don't know cancer: 
# Check cancer_m18y variable
# Replace don't knows with 1 if ever had cancer 
cham_tab1 %>% filter(cancer_mc1 == 999)

cham_tab1 <- cham_tab1 %>% mutate(cancer_mc1 = 
                                    ifelse(newid == 607 & cancer_m18y == 1, 1, cancer_mc1))

cham_tab1 %>% group_by(cancer_mc1) %>% count()
# 3 left over missing cancer info

# 2 missing asthma - can't adjust b/c other asthma vars from 16y follow-up 
# show no asthma - can't be sure they didn't get asthma between then and now


# filter to only include mat cog participants to check distribution of self-rated health scores
cham_tab1 <- cham_tab1 %>% filter(!is.na(age_qx_mc1))

# calculate overall self-rated health score
cham_tab1 <- cham_tab1 %>% mutate(summary_health = diab_mc1 + hbp_mc1 +
                                    heart_mc1 + cancer_mc1 + 
                                    asth_mc1 + thyr_mc1)


# check distribution of scores
cham_tab1 %>% group_by(summary_health) %>% count()


# re-code overall self-rated doctor diagnosed health scores into categories
# 0: no major health issues
# 1: 1 health issue
# 2-5: 2+ health issues
# 999: 1 don't know, rest no
# 1000: 1 yes, 1 dont know, rest no
# 1001: 2 yes, 1 dont know, rest no
# 1002: 3 yes, 1 dont know, rest no
# 2000: 2 yes, 2 don't know, rest no
# NA: missing

# Group by yeses:
# No major health issues: 0
# 1 health issue: 1, 1000
# 2-5 health issues: 2-7, 1001, 1002, 2000
# Missing: NA
cham_tab1 <- cham_tab1 %>% mutate(summary_health_cat_yes = 
                                    ifelse(summary_health == 0 | 
                                             summary_health == 999, "no health issues", 
                                           ifelse(summary_health == 1 |
                                                    summary_health == 1000, "1 health issue", 
                                                  ifelse(summary_health > 1 & summary_health <= 7| 
                                                           summary_health %in% c(1001, 1002, 2000), "2+",
                                                         NA))))

# create alcohol variable 
# 0 if haven't used alcohol since last visit or have had 0 drinks in the last 30 days; 1 otherwise
cham_tab1 %>% group_by(alcslv_mc1) %>% count()
cham_tab1 %>% group_by(alc30_mc1) %>% count()

cham_tab1 <- cham_tab1 %>% mutate(alc = ifelse(is.na(alcslv_mc1), NA, 
                                               ifelse(alcslv_mc1 == 0 | alc30_mc1 == 0, 0, 1)))

cham_tab1 %>% group_by(alc) %>% count()


# menopause variables
cham_tab1 %>% dplyr::select(menoind_mc1, menocat_mc1, agmeno_mc1)

# menocat_mc1: 1 = premenopause, 2 = natural menopause, 3 = surg meno
cham_tab1 <- cham_tab1 %>% mutate(menoind_mc1 = ifelse(menoind_mc1 == 0, 0, 1))

# pregnancy since last visit variables 
cham_tab1 <- cham_tab1 %>% mutate(live_preg_slv = ifelse(pregslv_mc1 == 1 & preg1out_mc1 == 1, 1, 0))


# select relevant table 1 variables and outcomes for regression models
cham_household <- cham_tab1 %>% dplyr::select(newid, cham, age_qx_mc1, ageusa_bl, yrsusa_matcog,
                                              age_came_back_to_us, ageusa_cat, ageusa_18,
                                              lang_exam_mc1, parent_educ,
                                              educcat_mom, married,
                                              work_mc1, workc_mc1, povcat_mc1, pov_binary, 
                                              ficat_mc1,
                                              hbp_bl, hbpage_bl, diab_bl, diabage_bl,
                                              cancer_bl, cancerage_bl, dlvrydt, date_matcog,
                                              qx_matcog,
                                              memory_mc1, execfun_rc_mc1, verbal_rc_mc1,
                                              global_rc_mc1, execfun_tbex_mc1,
                                              global_tbex_mc1,
                                              density1_mc1, density2_mc1,
                                              ppb_mc1, ppb2_mc1, 
                                              lvhome_n18_mc1, 
                                              rooms_mc1, lvhome_n_mc1,
                                              ppr_1, ppr_1_update, slpwkend_mc1, slpwkday_mc1, slpqual_mc1,
                                              restless_sleep,
                                              avg_sleep, less_6_hours,
                                              napfq_mc1, gadscore_mc1, depress_mc1, 
                                              depresscat_mc1, depress1_mc1, depress2_mc1, depress3_mc1, 
                                              depress4_mc1, depress5_mc1, depress6_mc1, depress7_mc1, 
                                              depress8_mc1, depress9_mc1, depress10_mc1,
                                              gad4cat_mc1, gad1_mc1, gad2_mc1, gad3_mc1, gad4_mc1, 
                                              gad5_mc1, gad6_mc1, gad7_mc1, 
                                              marstat_bl, worksp_bl, workc_bl, 
                                              povcat_bl, ipovcat_bl, ipovcat_bl_met,
                                              summary_health_cat_yes, summary_health,
                                              alcslv_mc1, alc30_mc1, alc,
                                              menoind_mc1, menocat_mc1, agmeno_mc1, people_per_household, 
                                              pregslv_mc1, pregslv_n_mc1, preg1out_mc1, preg1mult_mc1,
                                              live_preg_slv, tmtb_fail_mc1)



#### Correlation between number of children <18 and household density ####

# calculate correlation between number of children <18 and number of people per bedroom
cor.test(cham_household$lvhome_n18_mc1, cham_household$ppb_mc1)

# calculate correlation between number of children <18 and number of people per room
cor.test(cham_household$lvhome_n18_mc1, cham_household$density1_mc1)



#### Table 1. Overall Descriptive Characteristics, CHAMACOS Maternal Cognition Study ####

# filter out individuals missing covariates for table 1
cham_tab1_analytic <- cham_tab1 %>% filter(!is.na(age_qx_mc1) & !is.na(ageusa_18) & !is.na(lang_exam_mc1) &
                                             !is.na(parent_educ) & !is.na(educcat_mom) & !is.na(married) &
                                             !is.na(work_mc1) & !is.na(pov_binary) & !is.na(avg_sleep) & 
                                             !is.na(depress_mc1) & !is.na(gadscore_mc1) & !is.na(lvhome_n18_mc1) &
                                             !is.na(alc) & !is.na(summary_health_cat_yes) & !is.na(menoind_mc1) &
                                             !is.na(ppb2_mc1) & !is.na(ppr_1))

# set up variables for CreateTableOne 
vars <- c("age_qx_mc1", "ageusa_cat", "ageusa_18", "lang_exam_mc1", 
          "parent_educ",
          "educcat_mom", "married", "work_mc1", "workc_mc1", 
          "pov_binary", "ficat_mc1", "depress_mc1", "depresscat_mc1", 
          "gadscore_mc1", "gad4cat_mc1", "people_per_household",
          "density1_mc1", "ppb_mc1", "ppb2_mc1", "lvhome_n_mc1", "rooms_mc1",
          "lvhome_n18_mc1",
          "ppr_1", "slpwkend_mc1", "slpwkday_mc1", "slpqual_mc1", "avg_sleep", 
          "less_6_hours", "restless_sleep", "summary_health_cat_yes",
          "menoind_mc1", "menocat_mc1", "agmeno_mc1", "alcslv_mc1", "alc", "live_preg_slv")

# make age came to the US a factor variable and set 0 as reference
cham_tab1_analytic$ageusa_cat <- factor(cham_tab1_analytic$ageusa_cat)
cham_tab1_analytic$ageusa_cat <- relevel(cham_tab1_analytic$ageusa_cat, ref = "0")

# make married variable a factor variable and set not married as reference
cham_tab1_analytic$married <- factor(cham_tab1_analytic$married)
cham_tab1_analytic$married <- relevel(cham_tab1_analytic$married, ref = "not married")

# make ppb2_mc1 a factor variable 
cham_tab1_analytic$ppb2_mc1 <- factor(cham_tab1_analytic$ppb2_mc1)
cham_tab1_analytic$ppb.2_mc1 <- relevel(cham_tab1_analytic$ppb2_mc1, ref = "<=2")

# make sleep quality a factor variable 
cham_tab1_analytic$slpqual_mc1 <- factor(cham_tab1_analytic$slpqual_mc1)
cham_tab1_analytic$slpqual_mc1 <- relevel(cham_tab1_analytic$slpqual_mc1, ref = "very sound/restful")

# make depresscat_mc1 a factor variable
cham_tab1_analytic$depresscat_mc1 <- factor(cham_tab1_analytic$depresscat_mc1)
cham_tab1_analytic$depresscat_mc1 <- relevel(cham_tab1_analytic$depresscat_mc1, ref = "0")

# make gad4cat_mc1 a factor variable
cham_tab1_analytic$gad4cat_mc1 <- factor(cham_tab1_analytic$gad4cat_mc1)
cham_tab1_analytic$gad4cat_mc1 <- relevel(cham_tab1_analytic$gad4cat_mc1, ref = "0")

# make menoind variable a factor variable
cham_tab1_analytic$menoind_mc1 <- factor(cham_tab1_analytic$menoind_mc1)
cham_tab1_analytic$menoind_mc1 <- relevel(cham_tab1_analytic$menoind_mc1, ref = "0")

# make alcohol variable a factor variable
cham_tab1_analytic$alc <- factor(cham_tab1_analytic$alc)
cham_tab1_analytic$alc <- relevel(cham_tab1_analytic$alc, ref = "0")

# make alcohol since last visit variable a factor variable
cham_tab1_analytic$alcslv_mc1 <- factor(cham_tab1_analytic$alcslv_mc1)
cham_tab1_analytic$alcslv_mc1 <- relevel(cham_tab1_analytic$alcslv_mc1, ref = "0")

# make pregnancy since last visit variable a factor variable
cham_tab1_analytic$live_preg_slv <- factor(cham_tab1_analytic$live_preg_slv)
cham_tab1_analytic$live_preg_slv <- relevel(cham_tab1_analytic$live_preg_slv, ref = "0")

# create table 1
CreateTableOne(vars = vars, data = cham_tab1_analytic)

cham_tab1_analytic %>% group_by(tmtb_fail_mc1) %>% count()



#### eTable 5. Descriptive Characteristics of the Analytic Sample Stratified by Household Crowding Status (People per Bedroom) ####
vars <- c("age_qx_mc1", "ageusa_cat", "ageusa_18", 
          "lang_exam_mc1", "parent_educ",
          "educcat_mom", "married", "work_mc1", "workc_mc1", 
          "pov_binary", "ficat_mc1", "depress_mc1", "depresscat_mc1", 
          "gadscore_mc1", "gad4cat_mc1", "people_per_household",
          "density1_mc1", "ppb_mc1", "lvhome_n_mc1", "rooms_mc1",
          "lvhome_n18_mc1",
          "ppr_1", "slpwkend_mc1", "slpwkday_mc1", "slpqual_mc1", "avg_sleep", 
          "less_6_hours", "restless_sleep", "summary_health_cat_yes",
          "menoind_mc1", "menocat_mc1", "agmeno_mc1", "alc", "alcslv_mc1",
          "live_preg_slv")

# less than or equal to 2 people per bedroom
less_eq_2_ppb <- cham_tab1_analytic %>% filter(ppb2_mc1 == "<=2")

CreateTableOne(vars = vars, data = less_eq_2_ppb)

# more than 2 people per bedroom
more_2_ppb <- cham_tab1_analytic %>% filter(ppb2_mc1 == ">2")

CreateTableOne(vars = vars, data = more_2_ppb)



#### eTable 6. Descriptive Characteristics of the Analytic Sample Stratified by Household Crowding Status (People per Room) ####
vars <- c("age_qx_mc1", "ageusa_cat", "ageusa_18", 
          "lang_exam_mc1", "parent_educ",
          "educcat_mom", "married", "work_mc1", "workc_mc1", 
          "pov_binary", "ficat_mc1", "depress_mc1", "depresscat_mc1", 
          "gadscore_mc1", "gad4cat_mc1", "people_per_household",
          "density1_mc1", "ppb_mc1", "ppb2_mc1", "lvhome_n_mc1", "rooms_mc1",
          "lvhome_n18_mc1",
          "slpwkend_mc1", "slpwkday_mc1", "slpqual_mc1", "avg_sleep", 
          "less_6_hours", "restless_sleep", "summary_health_cat_yes",
          "menoind_mc1", "menocat_mc1", "agmeno_mc1", "alc", "alcslv_mc1",
          "live_preg_slv")

# less than or equal to 2 people per bedroom
less_eq_1_ppr <- cham_tab1_analytic %>% filter(ppr_1 == "<=1")

CreateTableOne(vars = vars, data = less_eq_1_ppr)

# more than 2 people per bedroom
more_1_ppr <- cham_tab1_analytic %>% filter(ppr_1 == ">1")

CreateTableOne(vars = vars, data = more_1_ppr)



#### Make all variables numeric and calculate analytic sample ####

# select relevant variables
cham_household_analytic <- cham_household %>% dplyr::select(age_qx_mc1, ageusa_18, lang_exam_mc1, educcat_mom, 
                                                            married, pov_binary, work_mc1, ppb2_mc1, 
                                                            ppr_1, 
                                                            lvhome_n18_mc1,
                                                            depress_mc1, depresscat_mc1,
                                                            gadscore_mc1, gad4cat_mc1,
                                                            summary_health_cat_yes, menoind_mc1,
                                                            alcslv_mc1, alc30_mc1, alc,
                                                            avg_sleep, memory_mc1, execfun_rc_mc1, verbal_rc_mc1, 
                                                            global_rc_mc1, density1_mc1, ppb2_mc1,
                                                            ppb_mc1,
                                                            slpwkday_mc1, slpwkend_mc1,
                                                            restless_sleep,
                                                            ppr_1_update, execfun_tbex_mc1, global_tbex_mc1)

## calculate number of individuals missing covariates
cham_household_analytic <- cham_household_analytic %>% filter(!is.na(ppb2_mc1) & !is.na(age_qx_mc1) & !is.na(ageusa_18) & 
                                                                !is.na(lang_exam_mc1) & !is.na(educcat_mom) &
                                                                !is.na(married) & !is.na(pov_binary) & !is.na(work_mc1) &
                                                                !is.na(lvhome_n18_mc1) &
                                                                !is.na(avg_sleep) & !is.na(depress_mc1) & !is.na(gadscore_mc1) & 
                                                                !is.na(summary_health_cat_yes) & !is.na(menoind_mc1) & !is.na(alc)) 

# calculate number of individuals missing neurocog outcomes                                                       
cham_household_analytic %>% filter(!is.na(memory_mc1))  
cham_household_analytic %>% filter(!is.na(execfun_rc_mc1))
cham_household_analytic %>% filter(!is.na(verbal_rc_mc1))
cham_household_analytic %>% filter(!is.na(global_rc_mc1)) 



#### Household crowding and cognitive outcomes: linear models with >2 people per bedroom ####

# Memory composite
memory_lm <- lm(memory_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                  educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                data = cham_household_analytic)

tidy(memory_lm, conf.int = TRUE)
nobs(memory_lm)


# Executive function
exec_lm <- lm(execfun_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
              data = cham_household_analytic)

tidy(exec_lm, conf.int = TRUE)
nobs(exec_lm)


# Verbal fluency
verbal_lm <- lm(verbal_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                  educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                data = cham_household_analytic)

tidy(verbal_lm, conf.int = TRUE)
nobs(verbal_lm)


# Global function
global_lm <- lm(global_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                  educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                data = cham_household_analytic)

tidy(global_lm, conf.int = TRUE)
nobs(global_lm)



#### Household crowding and cognitive outcomes: linear models with >1 person per room ####

# Memory composite
memory_lm_room <- lm(memory_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                     data = cham_household_analytic)

tidy(memory_lm_room, conf.int = TRUE)
nobs(memory_lm)


# Executive function
exec_lm_room <- lm(execfun_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                     educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                   data = cham_household_analytic)

tidy(exec_lm_room, conf.int = TRUE)
nobs(exec_lm_room)


# Verbal fluency
verbal_lm_room <- lm(verbal_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                     data = cham_household_analytic)

tidy(verbal_lm_room, conf.int = TRUE)
nobs(verbal_lm_room)


# Global function
global_lm_room <- lm(global_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                     data = cham_household_analytic)

tidy(global_lm_room, conf.int = TRUE)
nobs(global_lm_room)



#### Household crowding and cognitive outcomes: linear models with >2 PPBR and controlling for depression and anxiety ####

# Memory composite
memory_bed_dep_anx <- lm(memory_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                           depress_mc1 + gadscore_mc1, 
                         data = cham_household_analytic)

tidy(memory_bed_dep_anx, conf.int = TRUE)
nobs(memory_bed_dep_anx)


# Executive function
exec_bed_dep_anx <- lm(execfun_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                         educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                         depress_mc1 + gadscore_mc1, 
                       data = cham_household_analytic)

tidy(exec_bed_dep_anx, conf.int = TRUE)
nobs(exec_bed_dep_anx)


# Verbal fluency
verbal_bed_dep_anx <- lm(verbal_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                           depress_mc1 + gadscore_mc1, 
                         data = cham_household_analytic)

tidy(verbal_bed_dep_anx, conf.int = TRUE)
nobs(verbal_bed_dep_anx)


# Global function
global_bed_dep_anx <- lm(global_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                           depress_mc1 + gadscore_mc1, 
                         data = cham_household_analytic)

tidy(global_bed_dep_anx, conf.int = TRUE)
nobs(global_bed_dep_anx)



#### Household crowding and cognitive outcomes: linear models with >1 PPR and controlling for depression and anxiety ####

# Memory composite
memory_room_dep_anx <- lm(memory_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                            educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                            depress_mc1 + gadscore_mc1, 
                          data = cham_household_analytic)

tidy(memory_room_dep_anx, conf.int = TRUE)
nobs(memory_room_dep_anx)


# Executive function
exec_room_dep_anx <- lm(execfun_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                          educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                          depress_mc1 + gadscore_mc1, 
                        data = cham_household_analytic)

tidy(exec_room_dep_anx, conf.int = TRUE)
nobs(exec_room_dep_anx)


# Verbal fluency
verbal_room_dep_anx <- lm(verbal_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                            educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                            depress_mc1 + gadscore_mc1, 
                          data = cham_household_analytic)

tidy(verbal_room_dep_anx, conf.int = TRUE)
nobs(verbal_room_dep_anx)


# Global function
global_room_dep_anx <- lm(global_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                            educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                            depress_mc1 + gadscore_mc1, 
                          data = cham_household_analytic)

tidy(global_room_dep_anx, conf.int = TRUE)
nobs(global_room_dep_anx)



#### Household crowding and cognitive outcomes: linear models with >2 PPBR and controlling for sleep ####

# Memory composite
memory_bed_sleep <- lm(memory_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                         educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                         avg_sleep, 
                       data = cham_household_analytic)

tidy(memory_bed_sleep, conf.int = TRUE)
nobs(memory_bed_sleep)


# Executive function
exec_bed_sleep <- lm(execfun_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                       avg_sleep, 
                     data = cham_household_analytic)

tidy(exec_bed_sleep, conf.int = TRUE)
nobs(exec_bed_sleep)


# Verbal fluency
verbal_bed_sleep <- lm(verbal_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                         educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                         avg_sleep, 
                       data = cham_household_analytic)

tidy(verbal_bed_sleep, conf.int = TRUE)
nobs(verbal_bed_sleep)


# Global function
global_bed_sleep <- lm(global_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                         educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                         avg_sleep, 
                       data = cham_household_analytic)

tidy(global_bed_sleep, conf.int = TRUE)
nobs(global_bed_sleep)



#### Household crowding and cognitive outcomes: linear models with >1 PPR and controlling for sleep ####

# Memory composite
memory_room_sleep <- lm(memory_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                          educcat_mom + married + pov_binary + work_mc1 +lvhome_n18_mc1 +
                          avg_sleep, 
                        data = cham_household_analytic)

tidy(memory_room_sleep, conf.int = TRUE)
nobs(memory_room_sleep)


# Executive function
exec_room_sleep <- lm(execfun_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                        educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                        avg_sleep, 
                      data = cham_household_analytic)

tidy(exec_room_sleep, conf.int = TRUE)
nobs(exec_room_sleep)


# Verbal fluency
verbal_room_sleep <- lm(verbal_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                          educcat_mom + married + pov_binary + work_mc1 +lvhome_n18_mc1 +
                          avg_sleep, 
                        data = cham_household_analytic)

tidy(verbal_room_sleep, conf.int = TRUE)
nobs(verbal_room_sleep)


# Global function
global_room_sleep <- lm(global_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                          educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                          avg_sleep, 
                        data = cham_household_analytic)

tidy(global_room_sleep, conf.int = TRUE)
nobs(global_room_sleep)



#### Household crowding and sleep outcomes: linear models with >2 PPBR ####

# Weekday sleep
weekday_sleep_bed <- lm(slpwkday_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                          educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                        data = cham_household_analytic)

tidy(weekday_sleep_bed, conf.int = TRUE)
nobs(weekday_sleep_bed)


# Weekend sleep
weekend_sleep_bed <- lm(slpwkend_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                          educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1,
                        data = cham_household_analytic)

tidy(weekend_sleep_bed, conf.int = TRUE)
nobs(weekend_sleep_bed)


# Average sleep
avg_sleep_bed <- lm(avg_sleep ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                      educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1,
                    data = cham_household_analytic)

tidy(avg_sleep_bed, conf.int = TRUE)
nobs(avg_sleep_bed)


# Probability of restless sleep

## make restless sleep numeric variable
cham_household <- cham_household %>% mutate(restless_num = ifelse(restless_sleep == "restful", 0, 1))

## only include individuals with all exposure and covariate values
cham_household_glm <- cham_household %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) |
                                                    is.na(ageusa_18) | is.na(lang_exam_mc1) |
                                                    is.na(educcat_mom) | is.na(married) |
                                                    is.na(pov_binary) | is.na(work_mc1) | is.na(lvhome_n18_mc1)))

restless_sleep_bed <- geeglm(restless_num ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                               educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                             family = poisson(link = "log"), data = cham_household_glm,
                             id = seq.int(nrow(cham_household_glm)),
                             corstr = "exchangeable")

tidy(restless_sleep_bed, conf.int = TRUE)
nobs(restless_sleep_bed)

## pull term column from tidied dataframe  
term <- tidy(restless_sleep_bed, conf.int = TRUE) %>% pull(term) 

## pull statistic column from tidied dataframe
stat <- tidy(restless_sleep_bed, conf.int = TRUE) %>% pull(statistic) 

## pull p-values from tidied dataframe
pval <- tidy(restless_sleep_bed, conf.int = TRUE) %>% pull(p.value) 

## exponentiate point est, SE, lower CI, upper CI and combine with term, statistic, and p-value into dataframe
exp_results <- cbind(term, exp(tidy(restless_sleep_bed, conf.int = TRUE)[, c(2, 3, 6, 7)]), 
                     stat, pval) 

exp_results



#### Household crowding and depression/anxiety outcomes: linear models with >2 PPBR ####

# Depression
depression_bed <- lm(depress_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                     data = cham_household_analytic)

tidy(depression_bed, conf.int = TRUE)
nobs(depression_bed)


# Anxiety
anx_bed <- lm(gadscore_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
              data = cham_household_analytic)

tidy(anx_bed, conf.int = TRUE)
nobs(anx_bed)



##### Household crowding and sleep outcomes: linear models with >1 PPR #####

# Weekday sleep
weekday_sleep_room <- lm(slpwkday_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                         data = cham_household_analytic)

tidy(weekday_sleep_room, conf.int = TRUE)
nobs(weekday_sleep_room)


# Weekend sleep
weekend_sleep_room <- lm(slpwkend_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1,
                         data = cham_household_analytic)

tidy(weekend_sleep_room, conf.int = TRUE)
nobs(weekend_sleep_room)


# Average sleep
avg_sleep_room <- lm(avg_sleep ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1,
                     data = cham_household_analytic)

tidy(avg_sleep_room, conf.int = TRUE)
nobs(avg_sleep_room)


# Probability of restless sleep

## only include individuals with all exposure and covariate values
cham_household_glm %>% group_by(ppr_1) %>% count()

restless_sleep_room <- geeglm(restless_num ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                                educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                              family = poisson(link = "log"), data = cham_household_glm,
                              id = seq.int(nrow(cham_household_glm)),
                              corstr = "exchangeable")

tidy(restless_sleep_room, conf.int = TRUE)
nobs(restless_sleep_room)

## pull term column from tidied dataframe  
term <- tidy(restless_sleep_room, conf.int = TRUE) %>% pull(term) 

## pull statistic column from tidied dataframe
stat <- tidy(restless_sleep_room, conf.int = TRUE) %>% pull(statistic) 

## pull p-values from tidied dataframe
pval <- tidy(restless_sleep_room, conf.int = TRUE) %>% pull(p.value) 

## exponentiate point est, SE, lower CI, upper CI and combine with term, statistic, and p-value into dataframe
exp_results <- cbind(term, exp(tidy(restless_sleep_room, conf.int = TRUE)[, c(2, 3, 6, 7)]), 
                     stat, pval) 

exp_results



#### Household crowding and depression/anxiety outcomes: linear models with >1 PPR ####

# Depression
depression_room <- lm(depress_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                        educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                      data = cham_household_analytic)

tidy(depression_room, conf.int = TRUE)
nobs(depression_room)


# Anxiety
anx_room <- lm(gadscore_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                 educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
               data = cham_household_analytic)

tidy(anx_room, conf.int = TRUE)
nobs(anx_room)



#### Table 3. Mean Difference and 95% Confidence Intervals, Household Crowding and Candidate Mediators ####


### Mean difference in candidate mediators under different household crowding scenarios (>2 people per bedroom) using tmle3 package ###

# Depression outcome

# get rid of observations missing exposure, outcome, covars
cham_household_depress_mean <- cham_household_analytic %>% filter(!is.na(depress_mc1))
cham_household_depress_mean <- cham_household_depress_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                          is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                          is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))

# make variables factor vars
cham_household_depress_mean$ppb2_mc1 <- as.factor(cham_household_depress_mean$ppb2_mc1)
cham_household_depress_mean <- cham_household_depress_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                                ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                           ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_depress_mean$menoind_mc1 <- as.factor(cham_household_depress_mean$menoind_mc1)

cham_household_depress_mean$alc <- as.factor(cham_household_depress_mean$alc)

# create node list
node_list <- list(A = "ppb2_mc1",
                  Y = "depress_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_depress_mean, node_list)
cham_household_depress_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">2",
                     control_level = "<=2")


# Check learners by data type
sl3_list_learners(properties = "binomial");
sl3_list_learners(properties = "continuous")

# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_glmnet <- make_learner(Lrnr_glmnet)
lrnr_earth <- make_learner(Lrnr_earth)
lrnr_xgboost <- make_learner(Lrnr_xgboost)


sl_Y <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_glm, lrnr_glmnet, lrnr_earth, lrnr_xgboost))


sl_A <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_glm))


learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_depress_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Anxiety outcome
cham_household_anx_mean <- cham_household_analytic %>% filter(!is.na(gadscore_mc1))
cham_household_anx_mean <- cham_household_anx_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                  is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                  is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_anx_mean$ppb2_mc1 <- as.factor(cham_household_anx_mean$ppb2_mc1)
cham_household_anx_mean <- cham_household_anx_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                        ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                   ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_anx_mean$menoind_mc1 <- as.factor(cham_household_anx_mean$menoind_mc1)

cham_household_anx_mean$alc <- as.factor(cham_household_anx_mean$alc)

# create node list
node_list <- list(A = "ppb2_mc1",
                  Y = "gadscore_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_anx_mean, node_list)
cham_household_anx_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">2",
                     control_level = "<=2")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_anx_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Sleep outcome
# get rid of observations missing exposure, outcome, covars
cham_household_sleep_mean <- cham_household_analytic %>% filter(!is.na(avg_sleep))
cham_household_sleep_mean <- cham_household_sleep_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                      is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                      is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))

# make variables factor vars
cham_household_sleep_mean$ppb2_mc1 <- as.factor(cham_household_sleep_mean$ppb2_mc1)
cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                            ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                       ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_sleep_mean$menoind_mc1 <- as.factor(cham_household_sleep_mean$menoind_mc1)

cham_household_sleep_mean$alc <- as.factor(cham_household_sleep_mean$alc)

# create node list
node_list <- list(A = "ppb2_mc1",
                  Y = "avg_sleep",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_sleep_mean, node_list)
cham_household_sleep_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">2",
                     control_level = "<=2")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_sleep_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)



### Mean difference in candidate mediators under different household crowding scenarios (>1 person per room) using tmle3 package ###

# Depression outcome 

# get rid of observations missing exposure, outcome, covars
cham_household_depress_mean <- cham_household_analytic %>% filter(!is.na(depress_mc1))
cham_household_depress_mean <- cham_household_depress_mean %>% filter(!(is.na(ppr_1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                          is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                          is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))

# make variables factor vars
cham_household_depress_mean$ppr_1 <- as.factor(cham_household_depress_mean$ppr_1)
cham_household_depress_mean <- cham_household_depress_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                                ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                           ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_depress_mean <- cham_household_depress_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_depress_mean$menoind_mc1 <- as.factor(cham_household_depress_mean$menoind_mc1)

cham_household_depress_mean$alc <- as.factor(cham_household_depress_mean$alc)

# create node list
node_list <- list(A = "ppr_1",
                  Y = "depress_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_depress_mean, node_list)
cham_household_depress_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">1",
                     control_level = "<=1")

sl_Y <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_glm, lrnr_glmnet, lrnr_earth, lrnr_xgboost))


sl_A <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_glm))


learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_depress_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Anxiety outcome
cham_household_anx_mean <- cham_household_analytic %>% filter(!is.na(gadscore_mc1))
cham_household_anx_mean <- cham_household_anx_mean %>% filter(!(is.na(ppr_1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                  is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                  is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_anx_mean$ppr_1 <- as.factor(cham_household_anx_mean$ppr_1)
cham_household_anx_mean <- cham_household_anx_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                        ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                   ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_anx_mean <- cham_household_anx_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_anx_mean$menoind_mc1 <- as.factor(cham_household_anx_mean$menoind_mc1)

cham_household_anx_mean$alc <- as.factor(cham_household_anx_mean$alc)

# create node list
node_list <- list(A = "ppr_1",
                  Y = "gadscore_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_anx_mean, node_list)
cham_household_anx_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">1",
                     control_level = "<=1")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_anx_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Sleep outcome
# get rid of observations missing exposure, outcome, covars
cham_household_sleep_mean <- cham_household_analytic %>% filter(!is.na(avg_sleep))
cham_household_sleep_mean <- cham_household_sleep_mean %>% filter(!(is.na(ppr_1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                      is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                      is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))

# make variables factor vars
cham_household_sleep_mean$ppr_1 <- as.factor(cham_household_sleep_mean$ppr_1)
cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                            ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                       ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_sleep_mean$menoind_mc1 <- as.factor(cham_household_sleep_mean$menoind_mc1)

cham_household_sleep_mean$alc <- as.factor(cham_household_sleep_mean$alc)

# create node list
node_list <- list(A = "ppr_1",
                  Y = "avg_sleep",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_sleep_mean, node_list)
cham_household_sleep_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">1",
                     control_level = "<=1")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_sleep_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)



#### Table 2. Mean Difference and 95% Confidence Intervals, Household Crowding and Neurocognitive Z-Scores ####


### Mean difference in cognitive scores under different household crowding scenarios (>2 people per bedroom) using tmle3 package ###

# Memory outcome
cham_household_mem_mean <- cham_household_analytic %>% filter(!is.na(memory_mc1))
cham_household_mem_mean <- cham_household_mem_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                  is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                  is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_mem_mean$ppb2_mc1 <- as.factor(cham_household_mem_mean$ppb2_mc1)
cham_household_mem_mean <- cham_household_mem_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                        ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                   ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_mem_mean$menoind_mc1 <- as.factor(cham_household_mem_mean$menoind_mc1)

cham_household_mem_mean$alc <- as.factor(cham_household_mem_mean$alc)

# create node list
node_list <- list(A = "ppb2_mc1",
                  Y = "memory_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_mem_mean, node_list)
cham_household_mem_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">2",
                     control_level = "<=2")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_mem_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Executive function outcome
cham_household_exec_mean <- cham_household_analytic %>% filter(!is.na(execfun_rc_mc1))
cham_household_exec_mean <- cham_household_exec_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                    is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                    is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_exec_mean$ppb2_mc1 <- as.factor(cham_household_exec_mean$ppb2_mc1)
cham_household_exec_mean <- cham_household_exec_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                          ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                     ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_exec_mean$menoind_mc1 <- as.factor(cham_household_exec_mean$menoind_mc1)

cham_household_exec_mean$alc <- as.factor(cham_household_exec_mean$alc)

# create node list
node_list <- list(A = "ppb2_mc1",
                  Y = "execfun_rc_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_exec_mean, node_list)
cham_household_exec_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">2",
                     control_level = "<=2")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_exec_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Verbal fluency outcome
cham_household_verbal_mean <- cham_household_analytic %>% filter(!is.na(verbal_rc_mc1))
cham_household_verbal_mean <- cham_household_verbal_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                        is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                        is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_verbal_mean$ppb2_mc1 <- as.factor(cham_household_verbal_mean$ppb2_mc1)
cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                              ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                         ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_verbal_mean$menoind_mc1 <- as.factor(cham_household_verbal_mean$menoind_mc1)

cham_household_verbal_mean$alc <- as.factor(cham_household_verbal_mean$alc)

# create node list
node_list <- list(A = "ppb2_mc1",
                  Y = "verbal_rc_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_verbal_mean, node_list)
cham_household_verbal_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">2",
                     control_level = "<=2")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_verbal_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Global function outcome
cham_household_global_mean <- cham_household_analytic %>% filter(!is.na(global_rc_mc1))
cham_household_global_mean <- cham_household_global_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                        is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                        is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_global_mean$ppb2_mc1 <- as.factor(cham_household_global_mean$ppb2_mc1)
cham_household_global_mean <- cham_household_global_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                              ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_global_mean <- cham_household_global_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_global_mean <- cham_household_global_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                         ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_global_mean <- cham_household_global_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_global_mean <- cham_household_global_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_global_mean <- cham_household_global_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_global_mean$menoind_mc1 <- as.factor(cham_household_global_mean$menoind_mc1)

cham_household_global_mean$alc <- as.factor(cham_household_global_mean$alc)

# create node list
node_list <- list(A = "ppb2_mc1",
                  Y = "global_rc_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_global_mean, node_list)
cham_household_global_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">2",
                     control_level = "<=2")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_global_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)



### Mean difference in cognitive scores under different household crowding scenarios (>1 person per room) using tmle3 package ###

# Memory outcome
cham_household_mem_mean <- cham_household_analytic %>% filter(!is.na(memory_mc1))
cham_household_mem_mean <- cham_household_mem_mean %>% filter(!(is.na(ppr_1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                  is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                  is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_mem_mean$ppr_1 <- as.factor(cham_household_mem_mean$ppr_1)
cham_household_mem_mean <- cham_household_mem_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                        ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                   ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_mem_mean <- cham_household_mem_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_mem_mean$menoind_mc1 <- as.factor(cham_household_mem_mean$menoind_mc1)

cham_household_mem_mean$alc <- as.factor(cham_household_mem_mean$alc)

# create node list
node_list <- list(A = "ppr_1",
                  Y = "memory_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_mem_mean, node_list)
cham_household_mem_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">1",
                     control_level = "<=1")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_mem_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Executive function outcome
cham_household_exec_mean <- cham_household_analytic %>% filter(!is.na(execfun_rc_mc1))
cham_household_exec_mean <- cham_household_exec_mean %>% filter(!(is.na(ppr_1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                    is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                    is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_exec_mean$ppr_1 <- as.factor(cham_household_exec_mean$ppr_1)
cham_household_exec_mean <- cham_household_exec_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                          ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                     ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_exec_mean <- cham_household_exec_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_exec_mean$menoind_mc1 <- as.factor(cham_household_exec_mean$menoind_mc1)

cham_household_exec_mean$alc <- as.factor(cham_household_exec_mean$alc)

# create node list
node_list <- list(A = "ppr_1",
                  Y = "execfun_rc_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_exec_mean, node_list)
cham_household_exec_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">1",
                     control_level = "<=1")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_exec_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Verbal fluency outcome
cham_household_verbal_mean <- cham_household_analytic %>% filter(!is.na(verbal_rc_mc1))
cham_household_verbal_mean <- cham_household_verbal_mean %>% filter(!(is.na(ppr_1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                        is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                        is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_verbal_mean$ppr_1 <- as.factor(cham_household_verbal_mean$ppr_1)
cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                              ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                         ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_verbal_mean <- cham_household_verbal_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_verbal_mean$menoind_mc1 <- as.factor(cham_household_verbal_mean$menoind_mc1)

cham_household_verbal_mean$alc <- as.factor(cham_household_verbal_mean$alc)

# create node list
node_list <- list(A = "ppr_1",
                  Y = "verbal_rc_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_verbal_mean, node_list)
cham_household_verbal_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">1",
                     control_level = "<=1")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_verbal_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)


# Global function outcome
cham_household_global_mean <- cham_household_analytic %>% filter(!is.na(global_rc_mc1))
cham_household_global_mean <- cham_household_global_mean %>% filter(!(is.na(ppr_1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                        is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                        is.na(work_mc1) | is.na(lvhome_n18_mc1) | is.na(menoind_mc1) | is.na(alc)))


# make variables factor vars
cham_household_global_mean$ppr_1 <- as.factor(cham_household_global_mean$ppr_1)
cham_household_global_mean <- cham_household_global_mean %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                                              ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_global_mean <- cham_household_global_mean %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_global_mean <- cham_household_global_mean %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                                         ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_global_mean <- cham_household_global_mean %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_global_mean <- cham_household_global_mean %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_global_mean <- cham_household_global_mean %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_global_mean$menoind_mc1 <- as.factor(cham_household_global_mean$menoind_mc1)

cham_household_global_mean$alc <- as.factor(cham_household_global_mean$alc)

# create node list
node_list <- list(A = "ppr_1",
                  Y = "global_rc_mc1",
                  W = c("age_qx_mc1", "ageusa18_cat_num", "lang_exam_mc1_num", "educcat_num", "married_num", "pov_num", "work_num", "lvhome_n18_mc1",
                        "menoind_mc1", "alc"))

# process any missings
processed <- process_missing(cham_household_global_mean, node_list)
cham_household_global_mean <- processed$data
node_list <- processed$node_list

# create a spec object
ate_spec <- tmle_ATE(treatment_level = ">1",
                     control_level = "<=1")

learner_list <- list(A = sl_A, Y = sl_Y)

tmle_fit <- tmle3(ate_spec, cham_household_global_mean, node_list, learner_list)
print(tmle_fit)

estimates <- tmle_fit$summary$psi_transformed
print(estimates)



#### Create numeric exposure variables ####

# make people per bedroom numeric 
cham_household <- cham_household %>% mutate(ppb2_num = ifelse(ppb2_mc1 == "<=2", 0, 
                                                              ifelse(ppb2_mc1 == ">2", 1, NA)))
cham_household %>% group_by(ppb2_num) %>% count()

# make people per room numeric 
cham_household <- cham_household %>% mutate(ppr1_num = ifelse(ppr_1 == "<=1", 0, 
                                                              ifelse(ppr_1 == ">1", 1, NA)))

cham_household %>% group_by(ppr1_num) %>% count()



#### Make all variables numeric and calculate analytic sample ####

# select relevant variables
cham_household_med <- cham_household %>% dplyr::select(age_qx_mc1, ageusa_18, lang_exam_mc1, educcat_mom, 
                                                       married, pov_binary, work_mc1, ppb2_num, 
                                                       ppr1_num, lvhome_n18_mc1,
                                                       depress_mc1, depresscat_mc1,
                                                       gadscore_mc1, gad4cat_mc1,
                                                       summary_health_cat_yes, menoind_mc1,
                                                       alcslv_mc1, alc30_mc1, alc,
                                                       avg_sleep, sleep_sq, memory_mc1, execfun_rc_mc1, verbal_rc_mc1, 
                                                       global_rc_mc1, density1_mc1, ppb2_mc1,
                                                       ppb_mc1,
                                                       ppr_1_update, execfun_tbex_mc1, global_tbex_mc1)

cham_household_med <- cham_household_med %>% mutate(ageusa18_cat_num = ifelse(ageusa_18 == "0", 0, 
                                                                              ifelse(ageusa_18 == "<18", 1, 2)))

cham_household_med <- cham_household_med %>% mutate(lang_exam_mc1_num = ifelse(lang_exam_mc1 == "English", 0, 1))

cham_household_med <- cham_household_med %>% mutate(educcat_num = ifelse(educcat_mom == "<=6th grade", 0, 
                                                                         ifelse(educcat_mom == "7-11th grade", 1, 2)))

cham_household_med <- cham_household_med %>% mutate(married_num = ifelse(married == "not married", 0, 1))

cham_household_med <- cham_household_med %>% mutate(pov_num = ifelse(pov_binary == "At or below poverty", 0, 1))

cham_household_med <- cham_household_med %>% mutate(work_num = ifelse(work_mc1 == "Did not work since last visit", 0, 1))

cham_household_med <- cham_household_med %>% mutate(summary_health_num = 
                                                      ifelse(summary_health_cat_yes == "no health issues", 0, 
                                                             ifelse(summary_health_cat_yes == "1 health issue", 1, 2)))

cham_household_med <- cham_household_med %>% mutate(ppr_1_update = ifelse(ppr_1_update == "<=1", 0, 1))

# create numeric variable for depression
cham_household_med <- cham_household_med %>% mutate(depresscat_mc1 = ifelse(depresscat_mc1 == 0, 0, 1))

# create binary variable for anxiety (minimal/mild anxiety vs. moderate/severe anxiety)
cham_household_med <- cham_household_med %>% mutate(anx_bin = ifelse(gad4cat_mc1 == 0 | gad4cat_mc1 == 1, 0, 1))

## calculate number of individuals missing covariates
cham_household_med <- cham_household_med %>% filter(!is.na(ppb2_num) & !is.na(age_qx_mc1) & !is.na(ageusa18_cat_num) & 
                                                      !is.na(lang_exam_mc1_num) & !is.na(educcat_num) &
                                                      !is.na(married_num) & !is.na(pov_num) & !is.na(work_num) &
                                                      !is.na(lvhome_n18_mc1) &
                                                      !is.na(avg_sleep) & !is.na(depress_mc1) & !is.na(gadscore_mc1) & 
                                                      !is.na(summary_health_num) & !is.na(menoind_mc1) & !is.na(alc)) 
cham_household_med

# calculate number of individuals missing neurocog outcomes                                                       
cham_household_med %>% filter(!is.na(memory_mc1)) # same number missing as verbal_rc_mcl   
cham_household_med %>% filter(!is.na(global_rc_mc1)) # same number misssing as execfun_rc_mc1



#### Rename variables to use in Crumble package (>2 PPBR) ####

# Memory composite outcome
cham_household_med_mem <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                               educcat_num, married_num, pov_num, work_num,
                                                               lvhome_n18_mc1,
                                                               ppb2_num, depress_mc1, gadscore_mc1,
                                                               summary_health_num, menoind_mc1, 
                                                               alc, avg_sleep, memory_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep,
         Y = memory_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_mem <- as.data.frame(cham_household_med_mem)


# Executive function outcome
cham_household_med_exec <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                educcat_num, married_num, pov_num, work_num,
                                                                lvhome_n18_mc1,
                                                                ppb2_num, depress_mc1, gadscore_mc1,
                                                                summary_health_num, menoind_mc1, 
                                                                alc, avg_sleep, execfun_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = execfun_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_exec <- as.data.frame(cham_household_med_exec)


# Verbal fluency outcome
cham_household_med_verbal <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                  educcat_num, married_num, pov_num, work_num,
                                                                  lvhome_n18_mc1,
                                                                  ppb2_num, depress_mc1, gadscore_mc1,
                                                                  summary_health_num, menoind_mc1, 
                                                                  alc, avg_sleep, verbal_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = verbal_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_verbal <- as.data.frame(cham_household_med_verbal)


# Global function outcome
cham_household_med_global <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                  educcat_num, married_num, pov_num, work_num,
                                                                  lvhome_n18_mc1, 
                                                                  ppb2_num, depress_mc1, gadscore_mc1,
                                                                  summary_health_num, menoind_mc1, 
                                                                  alc, avg_sleep, global_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num,
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = global_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_global <- as.data.frame(cham_household_med_global)



#### Table 4. Coefficients and 95% Confidence Intervals, Household Crowding and Neurocognitive Z-Scores with Sleep Duration ####

### Mediation analysis (natural effects) with average sleep and > 2 people per bedroom (Crumble package) ####

# Memory outcome
data(cham_household_med_mem, package = "mma")
crumble_memory <- crumble(
  data = cham_household_med_mem,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)
crumble_memory


# Executive function outcome
data(cham_household_med_exec, package = "mma")
crumble_exec <- crumble(
  data = cham_household_med_exec,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)
crumble_exec


# Verbal fluency outcome
data(cham_household_med_verbal, package = "mma")
crumble_verbal <- crumble(
  data = cham_household_med_verbal,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)
crumble_verbal


# Global function outcome
data(cham_household_med_global, package = "mma")
crumble_global <- crumble(
  data = cham_household_med_global,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)
crumble_global



#### Rename variables to use in Crumble package (>1 PPR) ####

# Memory outcome
cham_household_med_mem_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                    educcat_num, married_num, pov_num, work_num,
                                                                    lvhome_n18_mc1,
                                                                    ppr1_num, depress_mc1, gadscore_mc1,
                                                                    summary_health_num, menoind_mc1, 
                                                                    alc, avg_sleep, memory_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = memory_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_mem_1ppr <- as.data.frame(cham_household_med_mem_1ppr)


# Executive function outcome
cham_household_med_exec_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                     educcat_num, married_num, pov_num, work_num,
                                                                     lvhome_n18_mc1,
                                                                     ppr1_num, depress_mc1, gadscore_mc1,
                                                                     summary_health_num, menoind_mc1, 
                                                                     alc, avg_sleep, execfun_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = execfun_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_exec_1ppr <- as.data.frame(cham_household_med_exec_1ppr)


# Verbal fluency outcome
cham_household_med_verbal_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                       educcat_num, married_num, pov_num, work_num,
                                                                       lvhome_n18_mc1,
                                                                       ppr1_num, depress_mc1, gadscore_mc1,
                                                                       summary_health_num, menoind_mc1, 
                                                                       alc, avg_sleep, verbal_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num,
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = verbal_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_verbal_1ppr <- as.data.frame(cham_household_med_verbal_1ppr)


# Global function outcome
cham_household_med_global_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                       educcat_num, married_num, pov_num, work_num,
                                                                       lvhome_n18_mc1,
                                                                       ppr1_num, depress_mc1, gadscore_mc1,
                                                                       summary_health_num, menoind_mc1, 
                                                                       alc, avg_sleep, global_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = global_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_global_1ppr <- as.data.frame(cham_household_med_global_1ppr)



### Mediation analysis (natural effects) with average sleep and > 1 person per room (Crumble package) ###

# Memory outcome
data(cham_household_med_mem_1ppr, package = "mma")
crumble_memory_ppr1 <- crumble(
  data = cham_household_med_mem_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)
crumble_memory_ppr1


# Executive function outcome
data(cham_household_med_exec_1ppr, package = "mma")
crumble_exec_ppr1 <- crumble(
  data = cham_household_med_exec_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)
crumble_exec_ppr1


# Verbal fluency outcome
data(cham_household_med_verbal_1ppr, package = "mma")
crumble_verbal_ppr1 <- crumble(
  data = cham_household_med_verbal_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)
crumble_verbal_ppr1


# Global function outcome
data(cham_household_med_global_1ppr, package = "mma")
crumble_global_ppr1 <- crumble(
  data = cham_household_med_global_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)
crumble_global_ppr1



#### Figure 2. Coefficients and 95% Confidence Intervals for Natural Direct and Indirect Estimates of the Association between Household Crowding and Neurocognitive Z-Scores, Mediated via Nightly Hours of Sleep ####


### Plot direct, indirect, total effects for >2 PPBR and >1 PPR with sleep mediation ###

exposure <- rep(c(rep(">2 people per bedroom", 3), rep(">1 person per room", 3)), 4)

type <- rep(c("Direct Effect", "Indirect Effect", 
              "Total Effect"), 8)

outcome <- c(rep("Memory", 6), rep("Executive Function", 6), rep("Verbal Fluency", 6), rep("Global Function", 6))


coefs <- c(crumble_memory$estimates$direct@x, 
           crumble_memory$estimates$indirect@x,
           crumble_memory$estimates$ate@x,
           crumble_memory_ppr1$estimates$direct@x, 
           crumble_memory_ppr1$estimates$indirect@x, 
           crumble_memory_ppr1$estimates$ate@x,
           crumble_exec$estimates$direct@x, 
           crumble_exec$estimates$indirect@x,
           crumble_exec$estimates$ate@x,
           crumble_exec_ppr1$estimates$direct@x, 
           crumble_exec_ppr1$estimates$indirect@x, 
           crumble_exec_ppr1$estimates$ate@x,
           crumble_verbal$estimates$direct@x, 
           crumble_verbal$estimates$indirect@x,
           crumble_verbal$estimates$ate@x,
           crumble_verbal_ppr1$estimates$direct@x, 
           crumble_verbal_ppr1$estimates$indirect@x, 
           crumble_verbal_ppr1$estimates$ate@x,
           crumble_global$estimates$direct@x, 
           crumble_global$estimates$indirect@x,
           crumble_global$estimates$ate@x,
           crumble_global_ppr1$estimates$direct@x, 
           crumble_global_ppr1$estimates$indirect@x, 
           crumble_global_ppr1$estimates$ate@x)

lower <- c(crumble_memory$estimates$direct@x - 1.96*crumble_memory$estimates$direct@std_error,
           crumble_memory$estimates$indirect@x - 1.96*crumble_memory$estimates$indirect@std_error,
           crumble_memory$estimates$ate@x - 1.96*crumble_memory$estimates$ate@std_error,
           crumble_memory_ppr1$estimates$direct@x - 1.96*crumble_memory_ppr1$estimates$direct@std_error, 
           crumble_memory_ppr1$estimates$indirect@x - 1.96*crumble_memory_ppr1$estimates$indirect@std_error,
           crumble_memory_ppr1$estimates$ate@x - 1.96*crumble_memory_ppr1$estimates$ate@std_error,
           crumble_exec$estimates$direct@x - 1.96*crumble_exec$estimates$direct@std_error,
           crumble_exec$estimates$indirect@x - 1.96*crumble_exec$estimates$indirect@std_error,
           crumble_exec$estimates$ate@x - 1.96*crumble_exec$estimates$ate@std_error,
           crumble_exec_ppr1$estimates$direct@x - 1.96*crumble_exec_ppr1$estimates$direct@std_error, 
           crumble_exec_ppr1$estimates$indirect@x - 1.96*crumble_exec_ppr1$estimates$indirect@std_error,
           crumble_exec_ppr1$estimates$ate@x - 1.96*crumble_exec_ppr1$estimates$ate@std_error,
           crumble_verbal$estimates$direct@x - 1.96*crumble_verbal$estimates$direct@std_error,
           crumble_verbal$estimates$indirect@x - 1.96*crumble_verbal$estimates$indirect@std_error,
           crumble_verbal$estimates$ate@x - 1.96*crumble_verbal$estimates$ate@std_error,
           crumble_verbal_ppr1$estimates$direct@x - 1.96*crumble_verbal_ppr1$estimates$direct@std_error, 
           crumble_verbal_ppr1$estimates$indirect@x - 1.96*crumble_verbal_ppr1$estimates$indirect@std_error,
           crumble_verbal_ppr1$estimates$ate@x - 1.96*crumble_verbal_ppr1$estimates$ate@std_error,
           crumble_global$estimates$direct@x - 1.96*crumble_global$estimates$direct@std_error,
           crumble_global$estimates$indirect@x - 1.96*crumble_global$estimates$indirect@std_error,
           crumble_global$estimates$ate@x - 1.96*crumble_global$estimates$ate@std_error,
           crumble_global_ppr1$estimates$direct@x - 1.96*crumble_global_ppr1$estimates$direct@std_error, 
           crumble_global_ppr1$estimates$indirect@x - 1.96*crumble_global_ppr1$estimates$indirect@std_error,
           crumble_global_ppr1$estimates$ate@x - 1.96*crumble_global_ppr1$estimates$ate@std_error)

upper <- c(crumble_memory$estimates$direct@x + 1.96*crumble_memory$estimates$direct@std_error,
           crumble_memory$estimates$indirect@x + 1.96*crumble_memory$estimates$indirect@std_error,
           crumble_memory$estimates$ate@x + 1.96*crumble_memory$estimates$ate@std_error,
           crumble_memory_ppr1$estimates$direct@x + 1.96*crumble_memory_ppr1$estimates$direct@std_error, 
           crumble_memory_ppr1$estimates$indirect@x + 1.96*crumble_memory_ppr1$estimates$indirect@std_error,
           crumble_memory_ppr1$estimates$ate@x + 1.96*crumble_memory_ppr1$estimates$ate@std_error,
           crumble_exec$estimates$direct@x + 1.96*crumble_exec$estimates$direct@std_error,
           crumble_exec$estimates$indirect@x + 1.96*crumble_exec$estimates$indirect@std_error,
           crumble_exec$estimates$ate@x + 1.96*crumble_exec$estimates$ate@std_error,
           crumble_exec_ppr1$estimates$direct@x + 1.96*crumble_exec_ppr1$estimates$direct@std_error, 
           crumble_exec_ppr1$estimates$indirect@x + 1.96*crumble_exec_ppr1$estimates$indirect@std_error,
           crumble_exec_ppr1$estimates$ate@x + 1.96*crumble_exec_ppr1$estimates$ate@std_error,
           crumble_verbal$estimates$direct@x + 1.96*crumble_verbal$estimates$direct@std_error,
           crumble_verbal$estimates$indirect@x + 1.96*crumble_verbal$estimates$indirect@std_error,
           crumble_verbal$estimates$ate@x + 1.96*crumble_verbal$estimates$ate@std_error,
           crumble_verbal_ppr1$estimates$direct@x + 1.96*crumble_verbal_ppr1$estimates$direct@std_error, 
           crumble_verbal_ppr1$estimates$indirect@x + 1.96*crumble_verbal_ppr1$estimates$indirect@std_error,
           crumble_verbal_ppr1$estimates$ate@x + 1.96*crumble_verbal_ppr1$estimates$ate@std_error,
           crumble_global$estimates$direct@x + 1.96*crumble_global$estimates$direct@std_error,
           crumble_global$estimates$indirect@x + 1.96*crumble_global$estimates$indirect@std_error,
           crumble_global$estimates$ate@x + 1.96*crumble_global$estimates$ate@std_error,
           crumble_global_ppr1$estimates$direct@x + 1.96*crumble_global_ppr1$estimates$direct@std_error, 
           crumble_global_ppr1$estimates$indirect@x + 1.96*crumble_global_ppr1$estimates$indirect@std_error,
           crumble_global_ppr1$estimates$ate@x + 1.96*crumble_global_ppr1$estimates$ate@std_error)

exposure <- as.data.frame(exposure)
type <- as.data.frame(type)
outcome <- as.data.frame(outcome)
coefs <- as.data.frame(coefs)
lower <- as.data.frame(lower)
upper <- as.data.frame(upper)

sleep_mediation_results <- cbind(exposure, type, outcome, coefs, lower, upper)


sleep_mediation_results$exposure <- factor(sleep_mediation_results$exposure, levels = c(">2 people per bedroom", ">1 person per room"))
sleep_mediation_results$outcome <- factor(sleep_mediation_results$outcome, levels = c("Memory", "Executive Function", "Verbal Fluency", "Global Function"))


variable_names <- list("Memory" = "Memory", 
                       "Executive Function" = "Executive \n Function",
                       "Verbal Fluency" = "Verbal \n Fluency", 
                       "Global Function" = "Global \n Function")

exposure_names <- list(">2 people per bedroom" = ">2 people \n per bedroom",
                       ">1 person per bedroom" = ">1 person \n per room")

variable_labeller2 <- function(variable,value){
  if (variable=='outcome') {
    return(variable_names[value])
  } else {
    return(exposure_names)
  }
}


ggplot(data = sleep_mediation_results, aes(x = type, y = coefs, color = type)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) + 
  facet_grid(outcome~exposure, scales = "free", space = "free_x", labeller = variable_labeller2) +
  ylab("") + xlab("") +
  theme_minimal(base_size = 10) +
  theme(legend.position="none") + 
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'black') + 
  scale_x_discrete(labels = c("Direct \n Effect", "Indirect \n Effect", "Total \n Effect")) + 
  ylim(-0.25, 0.25)



#### Sensitivity Analyses ####


#### eTable 7. Sensitivity Analysis with Path-Specific Effects #### 


### Rename variables to use in Crumble package (>2 PPBR) ###

# Memory outcome
cham_household_med_mem <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                               educcat_num, married_num, pov_num, work_num,
                                                               lvhome_n18_mc1,
                                                               ppb2_num, depress_mc1, gadscore_mc1,
                                                               summary_health_num, menoind_mc1, 
                                                               alc, avg_sleep, memory_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep,
         Y = memory_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num,
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_mem <- as.data.frame(cham_household_med_mem)


# Executive function outcome
cham_household_med_exec <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                educcat_num, married_num, pov_num, work_num,
                                                                lvhome_n18_mc1,
                                                                ppb2_num, depress_mc1, gadscore_mc1,
                                                                summary_health_num, menoind_mc1, 
                                                                alc, avg_sleep, execfun_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = execfun_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num,
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_exec <- as.data.frame(cham_household_med_exec)


# Verbal fluency outcome
cham_household_med_verbal <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                  educcat_num, married_num, pov_num, work_num,
                                                                  lvhome_n18_mc1,
                                                                  ppb2_num, depress_mc1, gadscore_mc1,
                                                                  summary_health_num, menoind_mc1, 
                                                                  alc, avg_sleep, verbal_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = verbal_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num,
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_verbal <- as.data.frame(cham_household_med_verbal)


# Global function outcome
cham_household_med_global <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                  educcat_num, married_num, pov_num, work_num,
                                                                  lvhome_n18_mc1, 
                                                                  ppb2_num, depress_mc1, gadscore_mc1,
                                                                  summary_health_num, menoind_mc1, 
                                                                  alc, avg_sleep, global_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num,
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = global_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num,
         Z4 = menoind_mc1,
         Z5 = alc,) %>% drop_na()

cham_household_med_global <- as.data.frame(cham_household_med_global)



### Mediation analysis with recanting twins (>2 PPBR) ###

# Memory outcome
data(cham_household_med_mem, package = "mma")
crumble(
  data = cham_household_med_mem,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8"),
  mediators = c("M"),
  moc = c("Z1", "Z2", "Z3", "Z4", "Z5"),
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "RT",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Executive function outcome
data(cham_household_med_exec, package = "mma")
crumble(
  data = cham_household_med_exec,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8"),
  mediators = c("M"),
  moc = c("Z1", "Z2", "Z3", "Z4", "Z5"),
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "RT",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Verbal fluency outcome
data(cham_household_med_verbal, package = "mma")
crumble(
  data = cham_household_med_verbal,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8"),
  mediators = c("M"),
  moc = c("Z1", "Z2", "Z3", "Z4", "Z5"),
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "RT",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Global function outcome
data(cham_household_med_global, package = "mma")
crumble(
  data = cham_household_med_global,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8"),
  mediators = c("M"),
  moc = c("Z1", "Z2", "Z3", "Z4", "Z5"),
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "RT",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)



### Rename variables to use in Crumble package (>1 PPR) ###

# Memory outcome
cham_household_med_mem_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                    educcat_num, married_num, pov_num, work_num,
                                                                    lvhome_n18_mc1,
                                                                    ppr1_num, depress_mc1, gadscore_mc1,
                                                                    summary_health_num, menoind_mc1, 
                                                                    alc, avg_sleep, memory_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = memory_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num,
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_mem_1ppr <- as.data.frame(cham_household_med_mem_1ppr)


# Executive function outcome
cham_household_med_exec_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                     educcat_num, married_num, pov_num, work_num,
                                                                     lvhome_n18_mc1,
                                                                     ppr1_num, depress_mc1, gadscore_mc1,
                                                                     summary_health_num, menoind_mc1, 
                                                                     alc, avg_sleep, execfun_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = execfun_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_exec_1ppr <- as.data.frame(cham_household_med_exec_1ppr)


# Verbal fluency outcome
cham_household_med_verbal_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                       educcat_num, married_num, pov_num, work_num,
                                                                       lvhome_n18_mc1,
                                                                       ppr1_num, depress_mc1, gadscore_mc1,
                                                                       summary_health_num, menoind_mc1, 
                                                                       alc, avg_sleep, verbal_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num,
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = verbal_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_verbal_1ppr <- as.data.frame(cham_household_med_verbal_1ppr)


# Global function outcome
cham_household_med_global_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                       educcat_num, married_num, pov_num, work_num,
                                                                       lvhome_n18_mc1,
                                                                       ppr1_num, depress_mc1, gadscore_mc1,
                                                                       summary_health_num, menoind_mc1, 
                                                                       alc, avg_sleep, global_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         M = avg_sleep, 
         Y = global_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num, 
         Z4 = menoind_mc1,
         Z5 = alc) %>% drop_na()

cham_household_med_global_1ppr <- as.data.frame(cham_household_med_global_1ppr)



### Mediation analysis with recanting twins (>1 PPR) ###

# Memory outcome
data(cham_household_med_mem_1ppr, package = "mma")
crumble(
  data = cham_household_med_mem_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8"),
  mediators = c("M"),
  moc = c("Z1", "Z2", "Z3", "Z4", "Z5"),
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "RT",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Executive function outcome
data(cham_household_med_exec_1ppr, package = "mma")
crumble(
  data = cham_household_med_exec_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8"),
  mediators = c("M"),
  moc = c("Z1", "Z2", "Z3", "Z4", "Z5"),
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "RT",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Verbal fluency outcome
data(cham_household_med_verbal_1ppr, package = "mma")
crumble(
  data = cham_household_med_verbal_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8"),
  mediators = c("M"),
  moc = c("Z1", "Z2", "Z3", "Z4", "Z5"),
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "RT",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Global function outcome
data(cham_household_med_global_1ppr, package = "mma")
crumble(
  data = cham_household_med_global_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8"),
  mediators = c("M"),
  moc = c("Z1", "Z2", "Z3", "Z4", "Z5"),
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "RT",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)



#### eTable 8. Sensitivity Analysis - Housing Density Calculated as Number of Rooms Minus 1 ####


### Sensitivity analysis with average sleep and >1 person per room using dens2_cat_mc1 ###


### Rename variables to use in Crumble package (>1 PPR) ###

# Memory outcome
cham_household_med_mem_1ppr_sens <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                         educcat_num, married_num, pov_num, work_num,
                                                                         lvhome_n18_mc1,
                                                                         ppr_1_update, depress_mc1, gadscore_mc1,
                                                                         summary_health_num, menoind_mc1, 
                                                                         alc, avg_sleep, sleep_sq, memory_mc1) %>% 
  rename(A = ppr_1_update, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1, 
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = memory_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_mem_1ppr_sens <- as.data.frame(cham_household_med_mem_1ppr_sens)


# Executive function outcome
cham_household_med_exec_1ppr_sens <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                          educcat_num, married_num, pov_num, work_num,
                                                                          lvhome_n18_mc1,
                                                                          ppr_1_update, depress_mc1, gadscore_mc1,
                                                                          summary_health_num, menoind_mc1, 
                                                                          alc, avg_sleep, sleep_sq, execfun_rc_mc1) %>% 
  rename(A = ppr_1_update, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1, 
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = execfun_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()


# Verbal fluency outcome
cham_household_med_verbal_1ppr_sens <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                            educcat_num, married_num, pov_num, work_num,
                                                                            lvhome_n18_mc1,
                                                                            ppr_1_update, depress_mc1, gadscore_mc1,
                                                                            summary_health_num, menoind_mc1, 
                                                                            alc, avg_sleep, sleep_sq, verbal_rc_mc1) %>% 
  rename(A = ppr_1_update, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1, 
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = verbal_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_verbal_1ppr_sens <- as.data.frame(cham_household_med_verbal_1ppr_sens)


# Global function outcome
cham_household_med_global_1ppr_sens <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                            educcat_num, married_num, pov_num, work_num,
                                                                            lvhome_n18_mc1,
                                                                            ppr_1_update, depress_mc1, gadscore_mc1,
                                                                            summary_health_num, menoind_mc1, 
                                                                            alc, avg_sleep, sleep_sq, global_rc_mc1) %>% 
  rename(A = ppr_1_update, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1, 
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = global_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_global_1ppr_sens <- as.data.frame(cham_household_med_global_1ppr_sens)



### Mediation analysis with >1 person per room using dens2_cat_mc1 (>1 PPR) ###

# Memory composite outcome
data(cham_household_med_mem_1ppr_sens, package = "mma")
crumble(
  data = cham_household_med_mem_1ppr_sens,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, #for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Executive function outcome
data(cham_household_med_exec_1ppr_sens, package = "mma")
crumble(
  data = cham_household_med_exec_1ppr_sens,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Verbal fluency outcome
data(cham_household_med_verbal_1ppr_sens, package = "mma")
crumble(
  data = cham_household_med_verbal_1ppr_sens,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Global function outcome
data(cham_household_med_global_1ppr_sens, package = "mma")
crumble(
  data = cham_household_med_global_1ppr_sens,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)



#### eTable 9. Sensitivity Analysis - Excluding Outliers of Nightly Hours of Sleep ####


# get rid of outliers for average sleep
no_miss_sleep <- cham_household_med %>% filter(!is.na(avg_sleep))
Q1 <- quantile(no_miss_sleep$avg_sleep, 0.25)
Q3 <- quantile(no_miss_sleep$avg_sleep, 0.75)
IQR <- IQR(no_miss_sleep$avg_sleep)

# 9 outliers
outliers <- no_miss_sleep %>% filter(avg_sleep < Q1 - 1.5*IQR | avg_sleep > Q3 + 1.5*IQR)

# subset original dataset to get rid of outliers
cham_household_no_outliers <- subset(cham_household_med, cham_household_med$avg_sleep > (Q1- 1.5*IQR) & cham_household_med$avg_sleep < (Q3 + 1.5*IQR))


### Rename variables to use in Crumble package (>2 PPBR) ###

# Memory outcome
cham_household_med_mem_sens <- cham_household_no_outliers %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                            educcat_num, married_num, pov_num, work_num,
                                                                            lvhome_n18_mc1,
                                                                            ppb2_num, depress_mc1, gadscore_mc1,
                                                                            summary_health_num, menoind_mc1, 
                                                                            alc, avg_sleep, sleep_sq, memory_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = memory_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_mem_sens <- as.data.frame(cham_household_med_mem_sens)


# Executive function outcome
cham_household_med_exec_sens <- cham_household_no_outliers %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                             educcat_num, married_num, pov_num, work_num,
                                                                             lvhome_n18_mc1,
                                                                             ppb2_num, depress_mc1, gadscore_mc1,
                                                                             summary_health_num, menoind_mc1, 
                                                                             alc, avg_sleep, sleep_sq, execfun_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = execfun_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_exec_sens <- as.data.frame(cham_household_med_exec_sens)


# Verbal fluency outcome
cham_household_med_verbal_sens <- cham_household_no_outliers %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                               educcat_num, married_num, pov_num, work_num,
                                                                               lvhome_n18_mc1,
                                                                               ppb2_num, depress_mc1, gadscore_mc1,
                                                                               summary_health_num, menoind_mc1, 
                                                                               alc, avg_sleep, sleep_sq, verbal_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1, 
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = verbal_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_verbal_sens <- as.data.frame(cham_household_med_verbal_sens)


# Global function outcome
cham_household_med_global_sens <- cham_household_no_outliers %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                               educcat_num, married_num, pov_num, work_num,
                                                                               lvhome_n18_mc1,
                                                                               ppb2_num, depress_mc1, gadscore_mc1,
                                                                               summary_health_num, menoind_mc1, 
                                                                               alc, avg_sleep, sleep_sq, global_rc_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1, 
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = global_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_global_sens <- as.data.frame(cham_household_med_global_sens)



### Mediation analysis excluding outliers (>2 PPBR) ###

# Memory outcome
data(cham_household_med_mem_sens, package = "mma")
crumble(
  data = cham_household_med_mem_sens,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Executive function outcome
data(cham_household_med_exec_sens, package = "mma")
crumble(
  data = cham_household_med_exec_sens,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Verbal fluency outcome
data(cham_household_med_verbal_sens, package = "mma")
crumble(
  data = cham_household_med_verbal_sens,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Global function outcome
data(cham_household_med_global_sens, package = "mma")
crumble(
  data = cham_household_med_global_sens,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)



### Rename variables to use in Crumble package (>1 PPR) ###

# Memory outcome
cham_household_med_mem_sens_1ppr <- cham_household_no_outliers %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                                 educcat_num, married_num, pov_num, work_num,
                                                                                 lvhome_n18_mc1,
                                                                                 ppr1_num, depress_mc1, gadscore_mc1,
                                                                                 summary_health_num, menoind_mc1, 
                                                                                 alc, avg_sleep, sleep_sq, memory_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = memory_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_mem_sens_1ppr <- as.data.frame(cham_household_med_mem_sens_1ppr)


# Executive function outcome
cham_household_med_exec_sens_1ppr <- cham_household_no_outliers %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                                  educcat_num, married_num, pov_num, work_num,
                                                                                  lvhome_n18_mc1,
                                                                                  ppr1_num, depress_mc1, gadscore_mc1,
                                                                                  summary_health_num, menoind_mc1, 
                                                                                  alc, avg_sleep, sleep_sq, execfun_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = execfun_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_exec_sens_1ppr <- as.data.frame(cham_household_med_exec_sens_1ppr)


# Verbal fluency outcome
cham_household_med_verbal_sens_1ppr <- cham_household_no_outliers %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                                    educcat_num, married_num, pov_num, work_num,
                                                                                    lvhome_n18_mc1,
                                                                                    ppr1_num, depress_mc1, gadscore_mc1,
                                                                                    summary_health_num, menoind_mc1, 
                                                                                    alc, avg_sleep, sleep_sq, verbal_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = verbal_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_verbal_sens_1ppr <- as.data.frame(cham_household_med_verbal_sens_1ppr)


# Global function outcome
cham_household_med_global_sens_1ppr <- cham_household_no_outliers %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                                    educcat_num, married_num, pov_num, work_num,
                                                                                    lvhome_n18_mc1,
                                                                                    ppr1_num, depress_mc1, gadscore_mc1,
                                                                                    summary_health_num, menoind_mc1, 
                                                                                    alc, avg_sleep, sleep_sq, global_rc_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = global_rc_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_global_sens_1ppr <- as.data.frame(cham_household_med_global_sens_1ppr)



### Mediation analysis excluding outliers (>1 PPR) ###

# Memory composite outcome
data(cham_household_med_mem_sens_1ppr, package = "mma")
crumble(
  data = cham_household_med_mem_sens_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Executive function outcome
data(cham_household_med_exec_sens_1ppr, package = "mma")
crumble(
  data = cham_household_med_exec_sens_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Verbal fluency outcome
data(cham_household_med_verbal_sens_1ppr, package = "mma")
crumble(
  data = cham_household_med_verbal_sens_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Global function outcome
data(cham_household_med_global_sens_1ppr, package = "mma")
crumble(
  data = cham_household_med_global_sens_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)



#### eTable 10. Sensitivity Analysis - Using Alternative Cognitive Scores ####


### Rename variables to use in Crumble package (>2 PPBR) ###

# Executive function outcome
cham_household_med_exec_sens_alt <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                         educcat_num, married_num, pov_num, work_num,
                                                                         lvhome_n18_mc1,
                                                                         ppb2_num, depress_mc1, gadscore_mc1,
                                                                         summary_health_num, menoind_mc1, 
                                                                         alc, avg_sleep, sleep_sq, execfun_tbex_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = execfun_tbex_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_exec_sens_alt <- as.data.frame(cham_household_med_exec_sens_alt)


# Global function outcome
cham_household_med_global_sens_alt <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                           educcat_num, married_num, pov_num, work_num,
                                                                           lvhome_n18_mc1,
                                                                           ppb2_num, depress_mc1, gadscore_mc1,
                                                                           summary_health_num, menoind_mc1, 
                                                                           alc, avg_sleep, sleep_sq, global_tbex_mc1) %>% 
  rename(A = ppb2_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = global_tbex_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_global_sens_alt <- as.data.frame(cham_household_med_global_sens_alt)



### Mediation analysis using alternative cognitive scores (>2 PPBR) ###

# Executive function outcome
data(cham_household_med_exec_sens_alt, package = "mma")
crumble(
  data = cham_household_med_exec_sens_alt,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Global function outcome
data(cham_household_med_global_sens_alt, package = "mma")
crumble(
  data = cham_household_med_global_sens_alt,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)



### Rename variables to use in Crumble package (>1 PPR) ###

# Executive function outcome
cham_household_med_exec_sens_alt_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                              educcat_num, married_num, pov_num, work_num,
                                                                              lvhome_n18_mc1,
                                                                              ppr1_num, depress_mc1, gadscore_mc1,
                                                                              summary_health_num, menoind_mc1, 
                                                                              alc, avg_sleep, sleep_sq, execfun_tbex_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1,
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = execfun_tbex_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_exec_sens_alt_1ppr <- as.data.frame(cham_household_med_exec_sens_alt_1ppr)


# Global function outcome
cham_household_med_global_sens_alt_1ppr <- cham_household_med %>% dplyr::select(age_qx_mc1, ageusa18_cat_num, lang_exam_mc1_num, 
                                                                                educcat_num, married_num, pov_num, work_num,
                                                                                lvhome_n18_mc1,
                                                                                ppr1_num, depress_mc1, gadscore_mc1,
                                                                                summary_health_num, menoind_mc1, 
                                                                                alc, avg_sleep, sleep_sq, global_tbex_mc1) %>% 
  rename(A = ppr1_num, 
         W1 = age_qx_mc1, 
         W2 = ageusa18_cat_num, 
         W3 = lang_exam_mc1_num, 
         W4 = educcat_num, 
         W5 = married_num, 
         W6 = pov_num, 
         W7 = work_num, 
         W8 = lvhome_n18_mc1, 
         W9 = menoind_mc1,
         W10 = alc,
         M1 = avg_sleep,
         Y = global_tbex_mc1, 
         Z1 = depress_mc1, 
         Z2 = gadscore_mc1, 
         Z3 = summary_health_num) %>% drop_na()

cham_household_med_global_sens_alt_1ppr <- as.data.frame(cham_household_med_global_sens_alt_1ppr)



### Mediation analysis using alternative cognitive scores (>1 PPR) ###

# Executive function outcome
data(cham_household_med_exec_sens_alt_1ppr, package = "mma")
crumble(
  data = cham_household_med_exec_sens_alt_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)


# Global function outcome
data(cham_household_med_global_sens_alt_1ppr, package = "mma")
crumble(
  data = cham_household_med_global_sens_alt_1ppr,
  trt = "A", 
  outcome = "Y",
  covar = c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10"),
  mediators = c("M1"),
  moc = NULL, # for natural effects
  d0 = \(data, trt) rep(0, nrow(data)), 
  d1 = \(data, trt) rep(1, nrow(data)), 
  effect = "N",
  learners = c("mean", "glm", "glmnet", "earth", "xgboost"), 
  nn_module = sequential_module(),
  control = crumble_control(crossfit_folds = 10L, epochs = 20L)
)



#### eTable 1. Comparison of Baseline Characteristics Overall and by Maternal Cognition Study Participants ####


#### Baseline Characteristics ####

### Current characteristics of MatCog sample ###
cham_tab1 %>% group_by(educcat_mom) %>% count()

cham_tab1 %>% group_by(ageusa_cat) %>% count()

# categorize baseline race/ethnicity
cham_tab1 <- cham_tab1 %>% mutate(raceeth_mom = ifelse(raceeth_mom == 1, "Mexican", 
                                                       ifelse(raceeth_mom == 2, "Mexican Indian",
                                                              ifelse(raceeth_mom == 3, "Mexican-American/Chicana", 
                                                                     ifelse(raceeth_mom == 4, "Other Latina", 
                                                                            ifelse(raceeth_mom == 5, "Asian or Pacific Islander", 
                                                                                   ifelse(raceeth_mom == 6, "White non-Latina", 
                                                                                          ifelse(raceeth_mom == 7, "Black non-Latina",
                                                                                                 ifelse(raceeth_mom == 8, "Other", NA)))))))))

cham_tab1 %>% group_by(raceeth_mom) %>% count()

cham_tab1 %>% group_by(lang_exam_mc1) %>% count()

cham_tab1 %>% group_by(povcat_mc1) %>% count()

cham_tab1 %>% group_by(married_mc1) %>% count()

cham_tab1 %>% group_by(work_mc1) %>% count()
cham_tab1 %>% group_by(wkag_mc1) %>% count()

cham_tab1 %>% group_by(lvhome_n_mc1) %>% count()

cham_tab1 %>% group_by(country_mom) %>% count()



### Baseline variables for MatCog sample ###

# categorize baseline high blood pressure
cham_tab1 <- cham_tab1 %>% mutate(hbp_bl = 
                                    ifelse(hbp_bl == 0, "No", 
                                           ifelse(hbp_bl == 1, "Yes", 
                                                  ifelse(hbp_bl == 9, "Don't know", hbp_bl))))

# categorize baseline diabetes
cham_tab1 <- cham_tab1 %>% mutate(diab_bl = 
                                    ifelse(diab_bl == 0, "No", 
                                           ifelse(diab_bl == 1, "Yes", 
                                                  ifelse(diab_bl == 9, "Don't know", diab_bl))))

# categorize baseline cancer
cham_tab1 <- cham_tab1 %>% mutate(cancer_bl = 
                                    ifelse(cancer_bl == 0, "No", 
                                           ifelse(cancer_bl == 1, "Yes", 
                                                  ifelse(cancer_bl == 9, "Don't know", cancer_bl))))


# set missing health variables as "no" 
cham_tab1 <- cham_tab1 %>% mutate(hbp_bl = ifelse(is.na(hbp_bl), "No", hbp_bl))

cham_tab1 <- cham_tab1 %>% mutate(diab_bl = ifelse(is.na(diab_bl), "No", diab_bl))

cham_tab1 <- cham_tab1 %>% mutate(cancer_bl = ifelse(is.na(cancer_bl), "No", cancer_bl))


# create collapsed chronic health conditions variable
cham_tab1 <- cham_tab1 %>% mutate(health = ifelse(hbp_bl == "No" & diab_bl == "No" &
                                                    cancer_bl == "No", "None", 
                                                  ifelse(hbp_bl == "Yes" | diab_bl == "Yes" |
                                                           cancer_bl == "Yes", "1+", NA)))

cham_tab1 %>% group_by(health) %>% count()

# recategorize those missing collapsed health score into "None" category
cham_tab1 <- cham_tab1 %>% mutate(health = ifelse(is.na(health), "None", health))

cham_tab1 %>% group_by(dens2_cat_mc1) %>% count()

cham_tab1 %>% group_by(povcat_bl) %>% count()

cham_tab1 %>% group_by(marstat_bl) %>% count()

cham_tab1 %>% group_by(worksp_bl) %>% count()


# check number of participants who participated in MATCOG sample
cham_all %>% group_by(qx_matcog) %>% count() 

# check number of participants who participated in CHAM1 and CHAM2
cham_all %>% group_by(cham) %>% count()

# remove individual with 314 who is missing almost all covariates
cham_base <- cham_all %>% filter(!newid == 314)

# original 1999 CHAMACOS study participants 
cham_1999_cohort <- cham_base %>% filter(cham == 1)


# 2009 CHAMACOS study participants (348 CHAM1 and 287 CHAM2)
# include those who started at CHAM2
cham_2_only <- cham_base %>% filter(cham == 2) # 287 CHAM2 participants

# include those who started at CHAM1 and participated in the 9 year questionnaire 
cham_1_9yr <- cham_base %>% filter(cham == 1) %>% filter(qx_9y == 1)

# create one dataset with participants who participated in 2009 questionnaire
cham_2009_cohort <- rbind(cham_2_only, cham_1_9yr)

# still missing 14 observations who participated in CHAM1 and did the 9-yr questionnaire
# must have erroneously been labeled as 0 for qx_9y?
# if person participated in CHAM1 and matcog, must have done the 9 yr questionnaire? 
cham_1_missing_9yr <- cham_base %>% filter(cham == 1 & qx_9y == 0) %>% filter(!is.na(age_qx_mc1))

# combine the above 11 individuals who participated in matcog from CHAM1 who were erroneously labeled as not having participated in 9 year follow-up
cham_2009_cohort <- rbind(cham_2009_cohort, cham_1_missing_9yr)

# still missing 3 CHAM1 individuals who participated in 2009 cohort
cham_1_missing_9yr_have_16_18yr <- cham_base %>% filter(!newid %in% c(cham_2009_cohort$newid)) %>% filter(!is.na(date_m18y) | !is.na(married_m16y))

# combine to complete 2009 cohort
cham_2009_cohort <- rbind(cham_2009_cohort, cham_1_missing_9yr_have_16_18yr)


# did not participate in matcog 
no_matcog_cham1 <- cham_1999_cohort %>% filter(!newid %in% c(cham_tab1$newid))
no_matcog_cham2 <- cham_2009_cohort %>% filter(!newid %in% c(cham_tab1$newid))

no_matcog <- rbind(no_matcog_cham1, no_matcog_cham2)
no_matcog <- distinct(no_matcog, newid, .keep_all = T)

no_matcog %>% group_by(educcat_mom) %>% count()
no_matcog %>% group_by(ageusa_cat) %>% count()
no_matcog %>% group_by(raceeth_mom) %>% count()
no_matcog %>% group_by(ipovcat_bl) %>% count()
no_matcog %>% group_by(marstat_bl) %>% count()
no_matcog %>% group_by(worksp_bl) %>% count()
no_matcog %>% group_by(health) %>% count()



