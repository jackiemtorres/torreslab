######################################

# Household Crowding paper
# Written by: Kelsey MacCuish
# Last edited by: Kelsey MacCuish
# Date of last edit: 12/6/2024

######################################

## outcomes: 
# primary:
# cognitive performance (memory, executive function, verbal fluency, global composition)

# secondary:
# sleep duration, quality, naps (slpwkend_mc1, slpwkday_mc1, slpqual_mc1, napfq_mc1)
# depression (depress_mc1)
# anxiety (gadscore_mc1)

## exposures: 
# primary:
# household crowding
# number of people per room (density) (density1_mc1)
# number of people per bedroom (ppb_mc1)
# <=2 vs. >2 people per bedroom (ppb2_mc1)

# secondary exposure: 
# total number of people in household (lvhome_n_mc1)
# number of children < 18 in household (lvhome_n18_mc1)

## covariates
# age
# nativity
# language of assessment 
# educational attainment
# occupational status (mat cog)
# marital status (mat cog)
# income/poverty (mat cog)


##### Load libraries
library(haven)
library(dplyr)
library(tableone)
library(broom)
library(ggplot2)
library(tidyr)
library(lme4)
library(ggdag)
library(dagitty)
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
library(SuperLearner)
library(medoutcon)
#devtools::install_github("nt-williams/HDmediation/HDmediation")
#devtools::install_github("nt-williams/HDmediation@mlr3superlearner", subdir = "HDmediation")
#remotes::install_github("nt-williams/crumble")
library(HDmediation)
library(crumble)
library(mma)
library(mlr3extralearners)
library(hal9001)
library(arm)
library(ggsci)
library(earth)
library(viridis)
library(xgboost)
library(psych)
library(corrr)
set.seed(124)


##### Read data #####
cham_all <- read_dta("torres_stressor_01i_nohsn.dta") # read in CHAMACOS data

head(cham_all) # first 6 rows of dataset
str(cham_all) # variables and variable types
dim(cham_all) # number of observations and variables
names(cham_all) # names of variables


##### Data cleaning #####
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
# 1 = English, 2 = Spanish
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
# 1 or 2 = married, otherwise not married
cham <- cham %>% mutate(married = 
                          ifelse(married_mc1 == 1 | 
                                   married_mc1 == 2, "married", "not married"))


# categorize worked since last visit variable
# 1 = worked since last visit, 0 = did not work since last visit
cham <- cham %>% mutate(work_mc1 = 
                          ifelse(work_mc1 == 1, "Worked since last visit",
                                 ifelse(work_mc1 == 0, "Did not work since last visit", NA)))


# categorize working now variable
# 1 = working now, 0 = not working now
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

# more than 1 person per bedroom: ppr_1
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

# NOTE: removed cholesterol from summary score - effects of cholesterol on outcome 
# accounted for with 
# For those who don't know cholesterol: 
# Check hdlcat_mc1 variable for cholesterol information
# Replace the don't knows with 1 if have high cholesterol, 0 if no high cholesterol
#cham_add_var %>% filter(chol_mc1 == 999)

#cham_add_var <- cham_add_var %>% mutate(chol_mc1 = ifelse(newid == 458 & trigcat_mc1 == 2, 0,
#ifelse(newid == 635 & trigcat_mc1 == 1, 0,
#ifelse(newid == 675 & trigcat_mc1 == 1, 0,
#ifelse(newid == 695 & trigcat_mc1 == 3, 1,
#ifelse(newid == 809 & trigcat_mc1 == 2, 0, 
#ifelse(newid == 812 & trigcat_mc1 == 3, 1,
#chol_mc1)))))))

# check that there are no more "don't know" responses for cholesterol
#cham_add_var %>% group_by(chol_mc1) %>% count()

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

# For those missing depression: 
# Check depresscat_mc1, depresscat_m18y, and depresscat_m16y variables
# Replace don't knows with values based on depression category from CESD or 
# past depression categories
# cham_add_var %>% filter(dep_mc1 == 999)

# cham_add_var <- cham_add_var %>% mutate(dep_mc1 =
#ifelse(newid %in% c(194, 254, 467, 656) & 
# depresscat_mc1 == 1, 1,
#ifelse(newid %in% c(695, 726) &
# depresscat_m16y == 1, 1,
#ifelse(newid == 332 & depresscat_m16y == 0, 0,
#dep_mc1))))

#cham_add_var %>% group_by(dep_mc1) %>% count()


# 9 missing anxiety
#cham_add_var %>% filter(anx_mc1 == 999) %>% dplyr::select(newid, gad4cat_mc1, gad7_4cat_m18y, 
#gad7_4cat_m16y)

#cham_add_var <- cham_add_var %>% mutate(anx_mc1 = 
#ifelse(newid %in% c(300, 332) &
#    gad4cat_mc1 == 0, 0, 
# ifelse(newid %in% c(60, 194, 254, 261, 695, 
#726, 815), 1, anx_mc1)))
# one person listed as NA but has gad4cat score of 2
#cham_add_var <- cham_add_var %>% mutate(anx_mc1 = ifelse(newid == 498, 1, anx_mc1))

#cham_add_var %>% group_by(anx_mc1) %>% count()


# filter to only include mat cog participants to check distribution of self-rated health scores
cham_tab1 <- cham_tab1 %>% filter(!is.na(age_qx_mc1))

# calculate overall self-rated health score
cham_tab1 <- cham_tab1 %>% mutate(summary_health = diab_mc1 + hbp_mc1 +
                                    heart_mc1 + cancer_mc1 + 
                                    asth_mc1 + thyr_mc1)


# check distribution of scores
cham_tab1 %>% group_by(summary_health) %>% count()


# re-code overall self-rated doctor diagnosed health scores into categories
# 0 - no major health issues
# 1 - 1 health issue
# 2-5: 2+ health issues
# 999 - 1 don't know, rest no
# 1000 = 1 yes, 1 dont know, rest no
# 1001 = 2 yes, 1 dont know, rest no
# 1002 = 3 yes, 1 dont know, rest no
# 2000 = 2 yes, 2 don't know, rest no
# NA - missing

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


##### Correlation between number of children <18 and household density #####

# calculate correlation between number of children <18 and number of people per bedroom
cor.test(cham_household$lvhome_n18_mc1, cham_household$ppb_mc1)

# calculate correlation between number of children <18 and number of people per room
cor.test(cham_household$lvhome_n18_mc1, cham_household$density1_mc1)


##### Table 1 #####
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


##### Table 1 stratified by <= 2 people per bedroom #####
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


##### Table 1 stratified by <= 1 person per room #####
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


##### Make all variables numeric and calculate analytic sample #####
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


##### Household crowding and cognitive outcomes: linear models with >2 people per bedroom #####
# memory composite
memory_lm <- lm(memory_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                  educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                data = cham_household_analytic)

tidy(memory_lm, conf.int = TRUE)
nobs(memory_lm)

# executive function
exec_lm <- lm(execfun_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
              data = cham_household_analytic)

tidy(exec_lm, conf.int = TRUE)
nobs(exec_lm)

# verbal fluency
verbal_lm <- lm(verbal_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                  educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                data = cham_household_analytic)

tidy(verbal_lm, conf.int = TRUE)
nobs(verbal_lm)

# global function
global_lm <- lm(global_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                  educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                data = cham_household_analytic)

tidy(global_lm, conf.int = TRUE)
nobs(global_lm)


##### Household crowding and cognitive outcomes: linear models with >1 person per room #####
# memory composite
memory_lm_room <- lm(memory_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                     data = cham_household_analytic)

tidy(memory_lm_room, conf.int = TRUE)
nobs(memory_lm)

# executive function
exec_lm_room <- lm(execfun_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                     educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                   data = cham_household_analytic)

tidy(exec_lm_room, conf.int = TRUE)
nobs(exec_lm_room)

# verbal fluency
verbal_lm_room <- lm(verbal_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                     data = cham_household_analytic)

tidy(verbal_lm_room, conf.int = TRUE)
nobs(verbal_lm_room)

# global function
global_lm_room <- lm(global_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                     data = cham_household_analytic)

tidy(global_lm_room, conf.int = TRUE)
nobs(global_lm_room)


##### Household crowding and cognitive outcomes: linear models with >2 PPBR and controlling for depression and anxiety #####
# memory composite
memory_bed_dep_anx <- lm(memory_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                           depress_mc1 + gadscore_mc1, 
                         data = cham_household_analytic)

tidy(memory_bed_dep_anx, conf.int = TRUE)
nobs(memory_bed_dep_anx)

# executive function
exec_bed_dep_anx <- lm(execfun_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                         educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                         depress_mc1 + gadscore_mc1, 
                       data = cham_household_analytic)

tidy(exec_bed_dep_anx, conf.int = TRUE)
nobs(exec_bed_dep_anx)

# verbal fluency
verbal_bed_dep_anx <- lm(verbal_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                           depress_mc1 + gadscore_mc1, 
                         data = cham_household_analytic)

tidy(verbal_bed_dep_anx, conf.int = TRUE)
nobs(verbal_bed_dep_anx)

# global function
global_bed_dep_anx <- lm(global_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                           depress_mc1 + gadscore_mc1, 
                         data = cham_household_analytic)

tidy(global_bed_dep_anx, conf.int = TRUE)
nobs(global_bed_dep_anx)


##### Household crowding and cognitive outcomes: linear models with >1 PPR and controlling for depression and anxiety #####
# memory composite
memory_room_dep_anx <- lm(memory_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                            educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                            depress_mc1 + gadscore_mc1, 
                          data = cham_household_analytic)

tidy(memory_room_dep_anx, conf.int = TRUE)
nobs(memory_room_dep_anx)

# executive function
exec_room_dep_anx <- lm(execfun_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                          educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                          depress_mc1 + gadscore_mc1, 
                        data = cham_household_analytic)

tidy(exec_room_dep_anx, conf.int = TRUE)
nobs(exec_room_dep_anx)

# verbal fluency
verbal_room_dep_anx <- lm(verbal_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                            educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                            depress_mc1 + gadscore_mc1, 
                          data = cham_household_analytic)

tidy(verbal_room_dep_anx, conf.int = TRUE)
nobs(verbal_room_dep_anx)

# global function
global_room_dep_anx <- lm(global_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                            educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                            depress_mc1 + gadscore_mc1, 
                          data = cham_household_analytic)

tidy(global_room_dep_anx, conf.int = TRUE)
nobs(global_room_dep_anx)


##### Household crowding and cognitive outcomes: linear models with >2 PPBR and controlling for sleep #####
# memory composite
memory_bed_sleep <- lm(memory_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                         educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                         avg_sleep, 
                       data = cham_household_analytic)

tidy(memory_bed_sleep, conf.int = TRUE)
nobs(memory_bed_sleep)

# executive function
exec_bed_sleep <- lm(execfun_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                       avg_sleep, 
                     data = cham_household_analytic)

tidy(exec_bed_sleep, conf.int = TRUE)
nobs(exec_bed_sleep)

# verbal fluency
verbal_bed_sleep <- lm(verbal_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                         educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                         avg_sleep, 
                       data = cham_household_analytic)

tidy(verbal_bed_sleep, conf.int = TRUE)
nobs(verbal_bed_sleep)

# global function
global_bed_sleep <- lm(global_rc_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                         educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                         avg_sleep, 
                       data = cham_household_analytic)

tidy(global_bed_sleep, conf.int = TRUE)
nobs(global_bed_sleep)


##### Household crowding and cognitive outcomes: linear models with >1 PPR and controlling for sleep #####
# memory composite
memory_room_sleep <- lm(memory_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                          educcat_mom + married + pov_binary + work_mc1 +lvhome_n18_mc1 +
                          avg_sleep, 
                        data = cham_household_analytic)

tidy(memory_room_sleep, conf.int = TRUE)
nobs(memory_room_sleep)

# executive function
exec_room_sleep <- lm(execfun_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                        educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                        avg_sleep, 
                      data = cham_household_analytic)

tidy(exec_room_sleep, conf.int = TRUE)
nobs(exec_room_sleep)

# verbal fluency
verbal_room_sleep <- lm(verbal_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                          educcat_mom + married + pov_binary + work_mc1 +lvhome_n18_mc1 +
                          avg_sleep, 
                        data = cham_household_analytic)

tidy(verbal_room_sleep, conf.int = TRUE)
nobs(verbal_room_sleep)

# global function
global_room_sleep <- lm(global_rc_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 +
                          educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1 +
                          avg_sleep, 
                        data = cham_household_analytic)

tidy(global_room_sleep, conf.int = TRUE)
nobs(global_room_sleep)


##### Household crowding and sleep outcomes: linear models with >2 PPBR #####
# weekday sleep
weekday_sleep_bed <- lm(slpwkday_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                          educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                        data = cham_household_analytic)

tidy(weekday_sleep_bed, conf.int = TRUE)
nobs(weekday_sleep_bed)


# weekend sleep
weekend_sleep_bed <- lm(slpwkend_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                          educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1,
                        data = cham_household_analytic)

tidy(weekend_sleep_bed, conf.int = TRUE)
nobs(weekend_sleep_bed)


# average sleep
avg_sleep_bed <- lm(avg_sleep ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                      educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1,
                    data = cham_household_analytic)

tidy(avg_sleep_bed, conf.int = TRUE)
nobs(avg_sleep_bed)


# probability of restless sleep

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


##### Household crowding and depression/anxiety outcomes: linear models with >2 PPBR #####
# depression
depression_bed <- lm(depress_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                     data = cham_household_analytic)

tidy(depression_bed, conf.int = TRUE)
nobs(depression_bed)


# anxiety
anx_bed <- lm(gadscore_mc1 ~ ppb2_mc1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
              data = cham_household_analytic)

tidy(anx_bed, conf.int = TRUE)
nobs(anx_bed)


##### Household crowding and sleep outcomes: linear models with >1 PPR #####
# weekday sleep
weekday_sleep_room <- lm(slpwkday_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                         data = cham_household_analytic)

tidy(weekday_sleep_room, conf.int = TRUE)
nobs(weekday_sleep_room)


# weekend sleep
weekend_sleep_room <- lm(slpwkend_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                           educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1,
                         data = cham_household_analytic)

tidy(weekend_sleep_room, conf.int = TRUE)
nobs(weekend_sleep_room)


# average sleep
avg_sleep_room <- lm(avg_sleep ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                       educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1,
                     data = cham_household_analytic)

tidy(avg_sleep_room, conf.int = TRUE)
nobs(avg_sleep_room)


# probability of restless sleep

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


##### Household crowding and depression/anxiety outcomes: linear models with >1 PPR #####
# depression
depression_room <- lm(depress_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                        educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
                      data = cham_household_analytic)

tidy(depression_room, conf.int = TRUE)
nobs(depression_room)


# anxiety
anx_room <- lm(gadscore_mc1 ~ ppr_1 + age_qx_mc1 + ageusa_18 + lang_exam_mc1 + 
                 educcat_mom + married + pov_binary + work_mc1 + lvhome_n18_mc1, 
               data = cham_household_analytic)

tidy(anx_room, conf.int = TRUE)
nobs(anx_room)


##### Mean difference in candidate mediators with >2 PPBR using TMLE ##### 
### depression outcome
cham_household_depress_mean <- cham_household_analytic %>% filter(!is.na(depress_mc1))
cham_household_depress_mean <- cham_household_depress_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                          is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                          is.na(work_mc1) | is.na(lvhome_n18_mc1)))

## make ppb2_mc1 variable numeric
cham_household_depress_mean <- cham_household_depress_mean %>% mutate(A = ifelse(ppb2_mc1 == "<=2", 0, 1))

cham_household_depress_mean <- cham_household_depress_mean %>% dplyr::select(depress_mc1, A, age_qx_mc1, ageusa_18, lang_exam_mc1,
                                                                             educcat_mom, married, pov_binary, work_mc1, lvhome_n18_mc1)
ObsData <- cham_household_depress_mean
ObsData <- ObsData %>% rename(Y = depress_mc1)

# transform continuous outcome to be within range of [0, 1]
min.Y <- min(ObsData$Y)
max.Y <- max(ObsData$Y)
ObsData$Y.bounded <- (ObsData$Y - min.Y)/(max.Y - min.Y)
summary(ObsData$Y.bounded)

# initial G-comp estimate
set.seed(124)
ObsData.noY <- dplyr::select(ObsData, !c(Y, Y.bounded))
Y.fit.sl <- SuperLearner(Y = ObsData$Y.bounded, 
                         X = ObsData.noY,
                         cvControl = list(V = 10L),
                         SL.library = c("SL.mean", 
                                        "SL.glm",
                                        "SL.glmnet",
                                        "SL.earth",
                                        "SL.xgboost"), 
                         #"SL.earth", 
                         #"SL.ranger"), 
                         method = "method.CC_nloglik", 
                         family = "gaussian")

# get initial predictions
ObsData$init.Pred <- predict(Y.fit.sl, newdata = ObsData.noY, 
                             type = "response")$pred

summary(ObsData$init.Pred)

# get predictions under treatment A = 1
ObsData.noY$A <- 1
ObsData$Pred.Y1 <- predict(Y.fit.sl, newdata = ObsData.noY, 
                           type = "response")$pred
summary(ObsData$Pred.Y1)

# get predictions under treatment A = 0
ObsData.noY$A <- 0
ObsData$Pred.Y0 <- predict(Y.fit.sl, newdata = ObsData.noY, 
                           type = "response")$pred

# get initial treatment effect estimate
ObsData$Pred.TE <- ObsData$Pred.Y1 - ObsData$Pred.Y0
summary(ObsData$Pred.TE)

# perform targeted improvement - propensity score model
set.seed(124)
ObsData.noYA <- dplyr::select(ObsData, !c(Y, Y.bounded, A, init.Pred,
                                          Pred.Y1, Pred.Y0, Pred.TE))

PS.fit.SL <- SuperLearner(Y = ObsData$A, 
                          X = ObsData.noYA, 
                          cvControl = list(V = 10L), 
                          SL.library = c("SL.mean", 
                                         "SL.glm",
                                         "SL.glmnet",
                                         "SL.earth",
                                         "SL.xgboost"), 
                          #"SL.earth", 
                          #"SL.ranger"),
                          method = "method.CC_nloglik", 
                          family = "binomial")

# get propensity score predictions
all.pred <- predict(PS.fit.SL, type = "response")

ObsData$PS.SL <- all.pred$pred
summary(ObsData$PS.SL)

tapply(ObsData$PS.SL, ObsData$A, summary)


# plot propensities
plot(density(ObsData$PS.SL[ObsData$A==0]), 
     col = "red", main = "")
lines(density(ObsData$PS.SL[ObsData$A==1]), 
      col = "blue", lty = 2)
legend("topright", c("<=2 PPBR",">2 PPBR"), 
       col = c("red", "blue"), lty=1:2)

# Estiamte H
ObsData$H.A1L <- (ObsData$A) / ObsData$PS.SL 
ObsData$H.A0L <- (1-ObsData$A) / (1- ObsData$PS.SL)
ObsData$H.AL <- ObsData$H.A1L - ObsData$H.A0L
summary(ObsData$H.AL)

tapply(ObsData$H.AL, ObsData$A, summary)

t(apply(cbind(-ObsData$H.A0L,ObsData$H.A1L), 
        2, summary)) 

# estimate epsilon
eps_mod <- glm(Y.bounded ~ -1 + H.A1L + H.A0L +  
                 offset(qlogis(init.Pred)), 
               family = "binomial",
               data = ObsData)
epsilon <- coef(eps_mod)  
epsilon["H.A1L"]
epsilon["H.A0L"]

eps_mod1 <- glm(Y.bounded ~ -1 + H.AL +
                  offset(qlogis(init.Pred)),
                family = "binomial",
                data = ObsData)
epsilon1 <- coef(eps_mod1) 
epsilon1 

ObsData$Pred.Y1.update <- plogis(qlogis(ObsData$Pred.Y1) +  
                                   epsilon["H.A1L"]*ObsData$H.A1L)
ObsData$Pred.Y0.update <- plogis(qlogis(ObsData$Pred.Y0) + 
                                   epsilon["H.A0L"]*ObsData$H.A0L)
summary(ObsData$Pred.Y1.update)

summary(ObsData$Pred.Y0.update)  

# effect estimate
ATE.TMLE.bounded.vector <- ObsData$Pred.Y1.update -  
  ObsData$Pred.Y0.update
summary(ATE.TMLE.bounded.vector) 

ATE.TMLE.bounded <- mean(ATE.TMLE.bounded.vector, 
                         na.rm = TRUE) 
ATE.TMLE.bounded 

# rescale effect estimate
ATE.TMLE <- (max.Y-min.Y)*ATE.TMLE.bounded   
ATE.TMLE 

# confidence interval estimation
ci.estimate <- function(data = ObsData, H.AL.components = 1){
  min.Y <- min(data$Y)
  max.Y <- max(data$Y)
  # transform predicted outcomes back to original scale
  if (H.AL.components == 2){
    data$Pred.Y1.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y1.update + min.Y
    data$Pred.Y0.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y0.update + min.Y
  } 
  if (H.AL.components == 1) {
    data$Pred.Y1.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y1.update1 + min.Y
    data$Pred.Y0.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y0.update1 + min.Y
  }
  EY1_TMLE1 <- mean(data$Pred.Y1.update.rescaled, 
                    na.rm = TRUE)
  EY0_TMLE1 <- mean(data$Pred.Y0.update.rescaled, 
                    na.rm = TRUE)
  # ATE efficient influence curve
  D1 <- data$A/data$PS.SL*
    (data$Y - data$Pred.Y1.update.rescaled) + 
    data$Pred.Y1.update.rescaled - EY1_TMLE1
  D0 <- (1 - data$A)/(1 - data$PS.SL)*
    (data$Y - data$Pred.Y0.update.rescaled) + 
    data$Pred.Y0.update.rescaled - EY0_TMLE1
  EIC <- D1 - D0
  # ATE variance
  n <- nrow(data)
  varHat.IC <- var(EIC, na.rm = TRUE)/n
  # ATE 95% CI
  if (H.AL.components == 2) {
    ATE.TMLE.CI <- c(ATE.TMLE - 1.96*sqrt(varHat.IC), 
                     ATE.TMLE + 1.96*sqrt(varHat.IC))
  }
  if (H.AL.components == 1) {
    ATE.TMLE.CI <- c(ATE.TMLE1 - 1.96*sqrt(varHat.IC), 
                     ATE.TMLE1 + 1.96*sqrt(varHat.IC))
  }
  return(ATE.TMLE.CI) 
}

CI2 <- ci.estimate(data = ObsData, H.AL.components = 2) 
CI2

### anxiety outcome
cham_household_anx_mean <- cham_household_analytic %>% filter(!is.na(gadscore_mc1))
cham_household_anx_mean <- cham_household_anx_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                  is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                  is.na(work_mc1) | is.na(lvhome_n18_mc1)))

## make ppb2_mc1 variable numeric
cham_household_anx_mean <- cham_household_anx_mean %>% mutate(A = ifelse(ppb2_mc1 == "<=2", 0, 1))

cham_household_anx_mean <- cham_household_anx_mean %>% dplyr::select(gadscore_mc1, A, age_qx_mc1, ageusa_18, lang_exam_mc1,
                                                                     educcat_mom, married, pov_binary, work_mc1, lvhome_n18_mc1)
ObsData <- cham_household_anx_mean
ObsData <- ObsData %>% rename(Y = gadscore_mc1)

# transform continuous outcome to be within range of [0, 1]
min.Y <- min(ObsData$Y)
max.Y <- max(ObsData$Y)
ObsData$Y.bounded <- (ObsData$Y - min.Y)/(max.Y - min.Y)
summary(ObsData$Y.bounded)

# initial G-comp estimate
set.seed(124)
ObsData.noY <- dplyr::select(ObsData, !c(Y, Y.bounded))
Y.fit.sl <- SuperLearner(Y = ObsData$Y.bounded, 
                         X = ObsData.noY,
                         cvControl = list(V = 10L),
                         SL.library = c("SL.mean", 
                                        "SL.glm",
                                        "SL.glmnet",
                                        "SL.earth",
                                        "SL.xgboost"), 
                         #"SL.earth", 
                         #"SL.ranger"), 
                         method = "method.CC_nloglik", 
                         family = "gaussian")

# get initial predictions
ObsData$init.Pred <- predict(Y.fit.sl, newdata = ObsData.noY, 
                             type = "response")$pred

summary(ObsData$init.Pred)

# get predictions under treatment A = 1
ObsData.noY$A <- 1
ObsData$Pred.Y1 <- predict(Y.fit.sl, newdata = ObsData.noY, 
                           type = "response")$pred
summary(ObsData$Pred.Y1)

# get predictions under treatment A = 0
ObsData.noY$A <- 0
ObsData$Pred.Y0 <- predict(Y.fit.sl, newdata = ObsData.noY, 
                           type = "response")$pred

# get initial treatment effect estimate
ObsData$Pred.TE <- ObsData$Pred.Y1 - ObsData$Pred.Y0
summary(ObsData$Pred.TE)

# perform targeted improvement - propensity score model
set.seed(124)
ObsData.noYA <- dplyr::select(ObsData, !c(Y, Y.bounded, A, init.Pred,
                                          Pred.Y1, Pred.Y0, Pred.TE))

PS.fit.SL <- SuperLearner(Y = ObsData$A, 
                          X = ObsData.noYA, 
                          cvControl = list(V = 10L), 
                          SL.library = c("SL.mean", 
                                         "SL.glm",
                                         "SL.glmnet",
                                         "SL.earth",
                                         "SL.xgboost"), 
                          #"SL.earth", 
                          #"SL.ranger"),
                          method = "method.CC_nloglik", 
                          family = "binomial")

# get propensity score predictions
all.pred <- predict(PS.fit.SL, type = "response")

ObsData$PS.SL <- all.pred$pred
summary(ObsData$PS.SL)

tapply(ObsData$PS.SL, ObsData$A, summary)


# plot propensities
plot(density(ObsData$PS.SL[ObsData$A==0]), 
     col = "red", main = "")
lines(density(ObsData$PS.SL[ObsData$A==1]), 
      col = "blue", lty = 2)
legend("topright", c("<=2 PPBR",">2 PPBR"), 
       col = c("red", "blue"), lty=1:2)

# Estiamte H
ObsData$H.A1L <- (ObsData$A) / ObsData$PS.SL 
ObsData$H.A0L <- (1-ObsData$A) / (1- ObsData$PS.SL)
ObsData$H.AL <- ObsData$H.A1L - ObsData$H.A0L
summary(ObsData$H.AL)

tapply(ObsData$H.AL, ObsData$A, summary)

t(apply(cbind(-ObsData$H.A0L,ObsData$H.A1L), 
        2, summary)) 

# estimate epsilon
eps_mod <- glm(Y.bounded ~ -1 + H.A1L + H.A0L +  
                 offset(qlogis(init.Pred)), 
               family = "binomial",
               data = ObsData)
epsilon <- coef(eps_mod)  
epsilon["H.A1L"]
epsilon["H.A0L"]

eps_mod1 <- glm(Y.bounded ~ -1 + H.AL +
                  offset(qlogis(init.Pred)),
                family = "binomial",
                data = ObsData)
epsilon1 <- coef(eps_mod1) 
epsilon1 

ObsData$Pred.Y1.update <- plogis(qlogis(ObsData$Pred.Y1) +  
                                   epsilon["H.A1L"]*ObsData$H.A1L)
ObsData$Pred.Y0.update <- plogis(qlogis(ObsData$Pred.Y0) + 
                                   epsilon["H.A0L"]*ObsData$H.A0L)
summary(ObsData$Pred.Y1.update)

summary(ObsData$Pred.Y0.update)  

# effect estimate
ATE.TMLE.bounded.vector <- ObsData$Pred.Y1.update -  
  ObsData$Pred.Y0.update
summary(ATE.TMLE.bounded.vector) 

ATE.TMLE.bounded <- mean(ATE.TMLE.bounded.vector, 
                         na.rm = TRUE) 
ATE.TMLE.bounded 

# rescale effect estimate
ATE.TMLE <- (max.Y-min.Y)*ATE.TMLE.bounded   
ATE.TMLE 

# confidence interval estimation
ci.estimate <- function(data = ObsData, H.AL.components = 1){
  min.Y <- min(data$Y)
  max.Y <- max(data$Y)
  # transform predicted outcomes back to original scale
  if (H.AL.components == 2){
    data$Pred.Y1.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y1.update + min.Y
    data$Pred.Y0.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y0.update + min.Y
  } 
  if (H.AL.components == 1) {
    data$Pred.Y1.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y1.update1 + min.Y
    data$Pred.Y0.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y0.update1 + min.Y
  }
  EY1_TMLE1 <- mean(data$Pred.Y1.update.rescaled, 
                    na.rm = TRUE)
  EY0_TMLE1 <- mean(data$Pred.Y0.update.rescaled, 
                    na.rm = TRUE)
  # ATE efficient influence curve
  D1 <- data$A/data$PS.SL*
    (data$Y - data$Pred.Y1.update.rescaled) + 
    data$Pred.Y1.update.rescaled - EY1_TMLE1
  D0 <- (1 - data$A)/(1 - data$PS.SL)*
    (data$Y - data$Pred.Y0.update.rescaled) + 
    data$Pred.Y0.update.rescaled - EY0_TMLE1
  EIC <- D1 - D0
  # ATE variance
  n <- nrow(data)
  varHat.IC <- var(EIC, na.rm = TRUE)/n
  # ATE 95% CI
  if (H.AL.components == 2) {
    ATE.TMLE.CI <- c(ATE.TMLE - 1.96*sqrt(varHat.IC), 
                     ATE.TMLE + 1.96*sqrt(varHat.IC))
  }
  if (H.AL.components == 1) {
    ATE.TMLE.CI <- c(ATE.TMLE1 - 1.96*sqrt(varHat.IC), 
                     ATE.TMLE1 + 1.96*sqrt(varHat.IC))
  }
  return(ATE.TMLE.CI) 
}

CI2 <- ci.estimate(data = ObsData, H.AL.components = 2) 
CI2


### sleep outcome
cham_household_sleep_mean <- cham_household_analytic %>% filter(!is.na(avg_sleep))
cham_household_sleep_mean <- cham_household_sleep_mean %>% filter(!(is.na(ppb2_mc1) | is.na(age_qx_mc1) | is.na(ageusa_18) | 
                                                                      is.na(lang_exam_mc1) | is.na(educcat_mom) | is.na(married) | is.na(pov_binary) | 
                                                                      is.na(work_mc1) | is.na(lvhome_n18_mc1)))

## make ppb2_mc1 variable numeric
cham_household_sleep_mean <- cham_household_sleep_mean %>% mutate(A = ifelse(ppb2_mc1 == "<=2", 0, 1))

cham_household_sleep_mean <- cham_household_sleep_mean %>% dplyr::select(avg_sleep, A, age_qx_mc1, ageusa_18, lang_exam_mc1,
                                                                         educcat_mom, married, pov_binary, work_mc1, lvhome_n18_mc1)
ObsData <- cham_household_sleep_mean
ObsData <- ObsData %>% rename(Y = avg_sleep)

# transform continuous outcome to be within range of [0, 1]
min.Y <- min(ObsData$Y)
max.Y <- max(ObsData$Y)
ObsData$Y.bounded <- (ObsData$Y - min.Y)/(max.Y - min.Y)
summary(ObsData$Y.bounded)

# initial G-comp estimate
set.seed(124)
ObsData.noY <- dplyr::select(ObsData, !c(Y, Y.bounded))
Y.fit.sl <- SuperLearner(Y = ObsData$Y.bounded, 
                         X = ObsData.noY,
                         cvControl = list(V = 10L),
                         SL.library = c("SL.mean", 
                                        "SL.glm",
                                        "SL.glmnet",
                                        "SL.earth",
                                        "SL.xgboost"), 
                         #"SL.earth", 
                         #"SL.ranger"), 
                         method = "method.CC_nloglik", 
                         family = "gaussian")

# get initial predictions
ObsData$init.Pred <- predict(Y.fit.sl, newdata = ObsData.noY, 
                             type = "response")$pred

summary(ObsData$init.Pred)

# get predictions under treatment A = 1
ObsData.noY$A <- 1
ObsData$Pred.Y1 <- predict(Y.fit.sl, newdata = ObsData.noY, 
                           type = "response")$pred
summary(ObsData$Pred.Y1)

# get predictions under treatment A = 0
ObsData.noY$A <- 0
ObsData$Pred.Y0 <- predict(Y.fit.sl, newdata = ObsData.noY, 
                           type = "response")$pred

# get initial treatment effect estimate
ObsData$Pred.TE <- ObsData$Pred.Y1 - ObsData$Pred.Y0
summary(ObsData$Pred.TE)

# perform targeted improvement - propensity score model
set.seed(124)
ObsData.noYA <- dplyr::select(ObsData, !c(Y, Y.bounded, A, init.Pred,
                                          Pred.Y1, Pred.Y0, Pred.TE))

PS.fit.SL <- SuperLearner(Y = ObsData$A, 
                          X = ObsData.noYA, 
                          cvControl = list(V = 10L), 
                          SL.library = c("SL.mean", 
                                         "SL.glm",
                                         "SL.glmnet",
                                         "SL.earth",
                                         "SL.xgboost"), 
                          #"SL.earth", 
                          #"SL.ranger"),
                          method = "method.CC_nloglik", 
                          family = "binomial")

# get propensity score predictions
all.pred <- predict(PS.fit.SL, type = "response")

ObsData$PS.SL <- all.pred$pred
summary(ObsData$PS.SL)

tapply(ObsData$PS.SL, ObsData$A, summary)


# plot propensities
plot(density(ObsData$PS.SL[ObsData$A==0]), 
     col = "red", main = "")
lines(density(ObsData$PS.SL[ObsData$A==1]), 
      col = "blue", lty = 2)
legend("topright", c("<=2 PPBR",">2 PPBR"), 
       col = c("red", "blue"), lty=1:2)

# Estiamte H
ObsData$H.A1L <- (ObsData$A) / ObsData$PS.SL 
ObsData$H.A0L <- (1-ObsData$A) / (1- ObsData$PS.SL)
ObsData$H.AL <- ObsData$H.A1L - ObsData$H.A0L
summary(ObsData$H.AL)

tapply(ObsData$H.AL, ObsData$A, summary)

t(apply(cbind(-ObsData$H.A0L,ObsData$H.A1L), 
        2, summary)) 

# estimate epsilon
eps_mod <- glm(Y.bounded ~ -1 + H.A1L + H.A0L +  
                 offset(qlogis(init.Pred)), 
               family = "binomial",
               data = ObsData)
epsilon <- coef(eps_mod)  
epsilon["H.A1L"]
epsilon["H.A0L"]

eps_mod1 <- glm(Y.bounded ~ -1 + H.AL +
                  offset(qlogis(init.Pred)),
                family = "binomial",
                data = ObsData)
epsilon1 <- coef(eps_mod1) 
epsilon1 

ObsData$Pred.Y1.update <- plogis(qlogis(ObsData$Pred.Y1) +  
                                   epsilon["H.A1L"]*ObsData$H.A1L)
ObsData$Pred.Y0.update <- plogis(qlogis(ObsData$Pred.Y0) + 
                                   epsilon["H.A0L"]*ObsData$H.A0L)
summary(ObsData$Pred.Y1.update)

summary(ObsData$Pred.Y0.update)  

# effect estimate
ATE.TMLE.bounded.vector <- ObsData$Pred.Y1.update -  
  ObsData$Pred.Y0.update
summary(ATE.TMLE.bounded.vector) 

ATE.TMLE.bounded <- mean(ATE.TMLE.bounded.vector, 
                         na.rm = TRUE) 
ATE.TMLE.bounded 

# rescale effect estimate
ATE.TMLE <- (max.Y-min.Y)*ATE.TMLE.bounded   
ATE.TMLE 

# confidence interval estimation
ci.estimate <- function(data = ObsData, H.AL.components = 1){
  min.Y <- min(data$Y)
  max.Y <- max(data$Y)
  # transform predicted outcomes back to original scale
  if (H.AL.components == 2){
    data$Pred.Y1.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y1.update + min.Y
    data$Pred.Y0.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y0.update + min.Y
  } 
  if (H.AL.components == 1) {
    data$Pred.Y1.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y1.update1 + min.Y
    data$Pred.Y0.update.rescaled <- 
      (max.Y- min.Y)*data$Pred.Y0.update1 + min.Y
  }
  EY1_TMLE1 <- mean(data$Pred.Y1.update.rescaled, 
                    na.rm = TRUE)
  EY0_TMLE1 <- mean(data$Pred.Y0.update.rescaled, 
                    na.rm = TRUE)
  # ATE efficient influence curve
  D1 <- data$A/data$PS.SL*
    (data$Y - data$Pred.Y1.update.rescaled) + 
    data$Pred.Y1.update.rescaled - EY1_TMLE1
  D0 <- (1 - data$A)/(1 - data$PS.SL)*
    (data$Y - data$Pred.Y0.update.rescaled) + 
    data$Pred.Y0.update.rescaled - EY0_TMLE1
  EIC <- D1 - D0
  # ATE variance
  n <- nrow(data)
  varHat.IC <- var(EIC, na.rm = TRUE)/n
  # ATE 95% CI
  if (H.AL.components == 2) {
    ATE.TMLE.CI <- c(ATE.TMLE - 1.96*sqrt(varHat.IC), 
                     ATE.TMLE + 1.96*sqrt(varHat.IC))
  }
  if (H.AL.components == 1) {
    ATE.TMLE.CI <- c(ATE.TMLE1 - 1.96*sqrt(varHat.IC), 
                     ATE.TMLE1 + 1.96*sqrt(varHat.IC))
  }
  return(ATE.TMLE.CI) 
}

CI2 <- ci.estimate(data = ObsData, H.AL.components = 2) 
CI2
