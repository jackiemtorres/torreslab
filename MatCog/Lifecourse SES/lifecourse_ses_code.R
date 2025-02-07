################################################################################################
# NAME:         Pesticides and cardiometabolic risk factors (pest_cmd_qgcomp.R)
# AUTHORS:      Kelsey MacCuish
# CREATED:      11/06/2023
# PURPOSE:      Script to run quantile G-computation to analyze associations of mixtures of pesticides, 
#               measured using Pesticide Use Reporting Data and urinary biomarkers, and cardio-metabolic risk 
#               factors in CHAMACOS moms
#  
################################################################################################

#### Load libraries ####
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
library(car)
library(ggspatial)



#### Read in data ####
cham <- read_dta("../Code/torres_stressor_01i_nohsn.dta") # read in CHAMACOS data
head(cham) # first 6 rows of dataset
str(cham) # variables and variable types
dim(cham) # number of observations and variables
names(cham) # names of variables



#### Data cleaning ####

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

# create overall ACES variable using two existing ACES variables (at 18 yr follow up and current aces collected for those missing at 18 yr follow up)
# note: aces_tot_mc1 categorized as 0 = 0, 1 = 1, 2 = 2, 3 = 3, 4 = 4, 5 = 5+, 6 = NA
# if missing aces value at 18 yr follow up use aces_tot_mc1 aces value
cham %>% group_by()
cham <- cham %>% mutate(aces_tot_m18y_cat = ifelse(aces_tot_m18y == 0, 0, 
                                                   ifelse(aces_tot_m18y == 1, 1, 
                                                          ifelse(aces_tot_m18y == 2, 2, 
                                                                 ifelse(aces_tot_m18y == 3, 3, 
                                                                        ifelse(aces_tot_m18y == 4, 4, 
                                                                               ifelse(aces_tot_m18y >= 5, 5, NA)))))))

# if missing 
cham <- cham %>% mutate(aces = ifelse(is.na(aces_tot_mc1), aces_tot_m18y_cat, aces_tot_mc1))

cham <- cham %>% mutate(aces = ifelse(newid == 64, NA, aces))

cham %>% filter(!is.na(age_qx_mc1)) %>% group_by(aces) %>% count()

cham %>% filter(!is.na(age_qx_mc1)) %>% group_by(aces_tot_mc1) %>% count();
cham %>% filter(!is.na(age_qx_mc1)) %>% group_by(aces_tot_m18y) %>% count();
cham %>% filter(!is.na(age_qx_mc1)) %>% group_by(aces_tot_m18y_cat) %>% count()

# make ACES a factor variable
cham$aces <- as.factor(cham$aces)

# dichotomize aces to >= 1 or 0 (or missing)
cham <- cham %>% mutate(aces_cat = ifelse(aces == 0, "0", 
                                          ifelse(aces == "" | aces == "PREFER NOT TO SAY", NA,
                                                 ">=1")))

cham <- cham %>% mutate(aces_cat = ifelse(aces == 6 & newid == 64, NA, aces_cat))

cham %>% filter(!is.na(age_qx_mc1)) %>% group_by(aces_cat) %>% count()

# categorize mom's own education
# 1: <= 6th grade, 2: 7-12th grade, 3: >= HS graduate
cham <- cham %>% mutate(educcat_mom = ifelse(educcat_mom == 1, "<=6th grade", 
                                             ifelse(educcat_mom == 2, "7-11th grade", 
                                                    ifelse(educcat_mom == 3, ">= HS graduate", NA))))

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

# povcat_mc1 
cham <- cham %>% mutate(povcat_mc1 = ifelse(povcat_mc1 == 1, "At or below poverty",
                                            ifelse(povcat_mc1 == 2, "200% FPL",
                                                   ifelse(povcat_mc1 == 3, ">=200% FPL", "Missing"))))
cham %>% group_by(povcat_mc1) %>% count()

# ficat_mc1
cham %>% filter(is.na(ficat_mc1))
cham %>% group_by(ficat_mc1) %>% count()
cham <- cham %>% mutate(ficat_mc1 = ifelse(ficat_mc1 == 1, "High or marginal food security", 
                                           ifelse(ficat_mc1 == 2, "Low food security", 
                                                  ifelse(ficat_mc1 == 3, "Very low food security", NA))))


# calculate age participant came back to the US at matcog visit 1
# yrsusa_matcog = yrsusa_dl + (date_matcog - dlvrydt)/365.25 
# yrsusa_matcog is years in the US at maternal cognition study visit 1
cham <- cham %>% mutate(age_came_back_to_us = age_qx_mc1 - yrsusa_matcog) 

# compare ageusa_bl with calculated age got back to the US
cham %>% dplyr::select(age_qx_mc1, yrsusa_matcog, ageusa_bl, age_came_back_to_us, momdl_age) %>% filter(!is.na(age_qx_mc1))

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
cham %>% filter(is.na(ageusa_cat)) %>% dplyr::select(newid, country_mom, ageusa_bl, momdl_age, ageusa_cat,
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

# categorize number of people per household 
cham <- cham %>% mutate(num_ppl_house = ifelse(lvhome_n_mc1 <= 5, "0-5", 
                                               ifelse(lvhome_n_mc1 >= 6, "6+", NA)))

cham %>% group_by(num_ppl_house) %>% count()

# select relevant table 1 variables and outcomes for regression models
cham_tab1_col <- cham %>% dplyr::select(newid, cham, age_qx_mc1, ageusa_bl, yrsusa_matcog,
                                        age_came_back_to_us, ageusa_cat,
                                        lang_exam_mc1, parent_educ,
                                        aces_cat, aces, educcat_mom, married,
                                        work_mc1, workc_mc1, povcat_mc1, ficat_mc1,
                                        hbp_bl, hbpage_bl, diab_bl, diabage_bl,
                                        cancer_bl, cancerage_bl, dlvrydt, date_matcog,
                                        qx_matcog,
                                        resp_simp_bqx, resp_simp_18y,
                                        memory_mc1, execfun_rc_mc1, verbal_rc_mc1,
                                        global_rc_mc1, execfun_tbex_mc1,
                                        global_tbex_mc1,
                                        marstat_bl, worksp_bl, workc_bl,
                                        wkag_mc1, 
                                        #wkfield_mc1, wkshed_mc1, 
                                        #wknursery_mc1, wktruck_mc1, wkagoth_mc1, 
                                        povcat_bl, ipovcat_bl, ipovcat_bl_met, lvhome_n_mc1, 
                                        num_ppl_house)

# filter out woman missing most baseline data
cham <- cham %>% filter(!newid == 314)


# filter to only include matcog participants
cham_tab1 <- cham_tab1_col %>% filter(!is.na(age_qx_mc1))



#### Create table 1 ####

# set up variables for CreateTableOne 
vars <- c("age_qx_mc1", "ageusa_cat", "lang_exam_mc1", 
          "parent_educ", "aces_cat", 
          "educcat_mom", "married", "work_mc1", "workc_mc1", 
          "povcat_mc1", "ficat_mc1")

# make age came to the US a factor variable and set 0 as reference
cham_tab1$ageusa_cat <- factor(cham_tab1$ageusa_cat)
cham_tab1$ageusa_cat <- relevel(cham_tab1$ageusa_cat, ref = "0")

# make married variable a factor variable and set not married as reference
cham_tab1$married <- factor(cham_tab1$married)
cham_tab1$married <- relevel(cham_tab1$married, ref = "not married")

# make aces_cat variable a factor variable and set 0 as reference
cham_tab1$aces_cat <- factor(cham_tab1$aces_cat)
cham_tab1$aces_cat <- relevel(cham_tab1$aces_cat, ref = "0")

# create table 1
CreateTableOne(vars = vars, data = cham_tab1)



#### Create new variables for descriptive lifeourse tables ####

# dichotomize parental education to any education or no formal education/unknown
cham_desc <- cham_tab1 %>% mutate(parent_educ_dich = 
                                    ifelse(parent_educ == "any formal educ", 
                                           "any formal educ", 
                                           "no formal or unknown"))

# create dichotomized poverty variable
# At or below poverty level and above poverty level (>200% poverty and Poverty-200%)
cham_desc <- cham_desc %>% mutate(pov_dich = 
                                    ifelse(povcat_mc1 == "At or below poverty", 
                                           "at or below poverty level", 
                                           ifelse(povcat_mc1 == ">=200% FPL" |
                                                    povcat_mc1 == "200% FPL", 
                                                  "above poverty level", 
                                                  NA)))

# dichotomize education variable
# Primary education or below (<=6th grade) and above primary (7-12th, more than HS)
cham_desc <- cham_desc %>% mutate(educ_dich = 
                                    ifelse(educcat_mom == "<=6th grade", 
                                           "primary or below", "above primary"))

# check contingency table between dichotomized poverty level and dichotomized education
cham_desc_analytic <- cham_desc %>% filter(!(is.na(aces_cat)))
table(cham_desc_analytic$educ_dich, cham_desc_analytic$pov_dich, useNA = "always")

# create variable that combines poverty and education
cham_desc <- cham_desc %>% mutate(pov_educ = 
                                    ifelse(pov_dich == "at or below poverty level" &
                                             educ_dich == "primary or below", 
                                           "primary or below + at or below poverty level",
                                           ifelse(pov_dich == "above poverty level" &
                                                    educ_dich == "primary or below", 
                                                  "primary or below + above poverty level",
                                                  ifelse(pov_dich == "at or below poverty level" &
                                                           educ_dich == "above primary", 
                                                         "above primary + at or below poverty level", 
                                                         ifelse(pov_dich == "above poverty level" &
                                                                  educ_dich == "above primary",
                                                                "above primary + above poverty level", NA)))))

# check new pov_educ variable counts
cham_desc_analytic <- cham_desc %>% filter(!(is.na(aces_cat)))
cham_desc_analytic %>% group_by(pov_educ) %>% count()

# generate descriptive lifecourse table
table(cham_desc_analytic$pov_educ, cham_desc_analytic$parent_educ_dich, useNA = "always")



#### MEMORY LINENAR REGRESSION MODELS ####


### parent education and memory linear regression ###

# copy descriptive table dataset to include relevant variables for parent education model
cham_memory <- cham_desc

# factor parent education exposure variable and set no education as reference
cham_memory$parent_educ <- factor(cham_memory$parent_educ, 
                                  levels = c("No formal education", 
                                             "Dont know/dont remember", 
                                             "any formal educ"))

# run linear regression model with memory composite score as outcome and
# parent education as exposure
# adjust for age (continuous) and years in US (categorical)
parent_educ_memory <- lm(memory_mc1 ~ parent_educ + age_qx_mc1 + ageusa_cat, 
                         data = cham_memory)

# tidy model and include confidence intervals in output
tidy(parent_educ_memory, conf.int = TRUE)

# check number of observations used in linear model
nobs(parent_educ_memory)


### own education and memory linear regression ###

# factor own education exposure variable and set <=6th grade as reference
cham_memory$educcat_mom <- factor(cham_memory$educcat_mom, 
                                  levels = c("<=6th grade", 
                                             "7-11th grade", 
                                             ">= HS graduate"))

# run linear regression model with memory composite score as outcome and
# own education as main exposure
# adjust for parent education (categorical), age (continuous), years in US (categorical), ACES (categorical)
own_educ_memory <- lm(memory_mc1 ~ educcat_mom + parent_educ + 
                        age_qx_mc1 + factor(ageusa_cat) + factor(aces_cat), 
                      data = cham_memory)

# tidy model and include confidence intervals in output
tidy(own_educ_memory, conf.int = TRUE)

# check number of observations used in linear model
nobs(own_educ_memory)


### poverty level and memory linear regression ###

# Set NA values as "missing" in dichotomized poverty variable to create categorical var (want to calculate coef for those who are missing poverty exposure in regression model)
cham_memory <- cham_memory %>% 
  mutate(pov_cat = 
           ifelse(is.na(pov_dich), "Missing", pov_dich))

# factor poverty variable and set at or below poverty as reference
cham_memory$pov_cat <- factor(cham_memory$pov_cat, 
                              levels = c("at or below poverty level", 
                                         "above poverty level", 
                                         "Missing"))


poverty_memory <- lm(memory_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                       age_qx_mc1 + factor(ageusa_cat) + factor(aces_cat) + 
                       factor(married) + factor(work_mc1), 
                     data = cham_memory)

# tidy model and include confidence intervals in output
tidy(poverty_memory, conf.int = TRUE)

# check number of observations used in linear model 
nobs(poverty_memory)



#### EXECTUTIVE FUNCTION LINEAR REGRESSION MODELS ####


### parent education and executive function linear regression ###

# run linear regression model with executive function score as outcome and
# parent education as main exposure
# adjust for age (continuous), years in US (categorical)
parent_educ_exec <- lm(execfun_rc_mc1 ~ parent_educ + 
                         age_qx_mc1 + factor(ageusa_cat), 
                       data = cham_memory)

# tidy model and include confidence intervals in output
tidy(parent_educ_exec, conf.int = TRUE)

# check number of observations used in linear model 
nobs(parent_educ_exec)


### own education and executive function linear regression ###

# run linear regression model with executive function score as outcome and
# own education as main exposure
# adjust for parent education (categorical), age (continuous), years in US (categorical), ACES (categorical)
own_educ_exec <- lm(execfun_rc_mc1 ~ educcat_mom + parent_educ + 
                      age_qx_mc1 + factor(ageusa_cat) + factor(aces_cat), 
                    data = cham_memory)

# tidy model and include confidence intervals in output
tidy(own_educ_exec, conf.int = TRUE)

# check number of observations used in linear model 
nobs(own_educ_exec)


### poverty level and executive function linear regression ###

# run linear regression model with executive function score as outcome and
# poverty as main exposure
# adjust for own education (categorical), parent education (categorical), age (continuous), years in US (categorical), ACES (categorical), marital status (categorical), employment status (categorical)
poverty_exec <- lm(execfun_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                     age_qx_mc1 + factor(ageusa_cat) + factor(aces_cat) + 
                     factor(married) + factor(work_mc1), 
                   data = cham_memory)

# tidy model and include confidence intervals in output
tidy(poverty_exec, conf.int = TRUE)

# check number of observations used in linear model 
nobs(poverty_exec)



#### VERBAL FLUENCY LINEAR REGRESSION MODELS ####


### parental education and verbal fluency lienar regression ###

# run linear regression model with verbal fluency score as outcome and
# parent education as main exposure
# adjust for age (continuous), years in US (categorical)
parent_educ_verbal <- lm(verbal_rc_mc1 ~ parent_educ + 
                           age_qx_mc1 + factor(ageusa_cat), 
                         data = cham_memory)

# tidy model and include confidence intervals in output
tidy(parent_educ_verbal, conf.int = TRUE)

# check number of observations used in linear model 
nobs(parent_educ_verbal)


### own education and verbal fluency linear regression ###

# run linear regression model with verbal fluency score as outcome and
# own education as main exposure
# adjust for parent education (categorical), age (continuous), years in US (categorical), ACES (categorical)
own_educ_verbal <- lm(verbal_rc_mc1 ~ educcat_mom + parent_educ + 
                        age_qx_mc1 + factor(ageusa_cat) + factor(aces_cat), 
                      data = cham_memory)

# tidy model and include confidence intervals in output
tidy(own_educ_verbal, conf.int = TRUE)

# check number of observations used in linear model 
nobs(own_educ_verbal)


### poverty level and verbal fluency linear regression ###

# run linear regression model with verbal fluency score as outcome and
# poverty as main exposure
# adjust for own education (categorical), parent education (categorical), age (continuous), years in US (categorical), ACES (categorical), marital status (categorical), employment status (categorical)
poverty_verbal <- lm(verbal_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                       age_qx_mc1 + factor(ageusa_cat) + factor(aces_cat) + 
                       factor(married) + factor(work_mc1), 
                     data = cham_memory)

# tidy model and include confidence intervals in output
tidy(poverty_verbal, conf.int = TRUE)

# check number of observations used in linear model 
nobs(poverty_verbal)



#### GLOBAL FUNCTION LINEAR REGRESSION MODELS ####


### parental education and global function linear regression ###

# run linear regression model with global function score as outcome and
# parent education as main exposure
# adjust for age (continuous), years in US (categorical)
parent_educ_global <- lm(global_rc_mc1 ~ parent_educ + 
                           age_qx_mc1 + factor(ageusa_cat), 
                         data = cham_memory)

# tidy model and include confidence intervals in output
tidy(parent_educ_global, conf.int = TRUE)

# check number of observations used in linear model 
nobs(parent_educ_global)


### own education and global function linear regression ###

# run linear regression model with global function score as outcome and
# own education as main exposure
# adjust for parent education (categorical), age (continuous), years in US (categorical), ACES (categorical)
own_educ_global <- lm(global_rc_mc1 ~ educcat_mom + parent_educ + 
                        age_qx_mc1 + factor(ageusa_cat) + factor(aces_cat), 
                      data = cham_memory)

# tidy model and include confidence intervals in output
tidy(own_educ_global, conf.int = TRUE)

# check number of observations used in linear model 
nobs(own_educ_global)


### poverty level and global function linear regression ###

# run linear regression model with global function score as outcome and
# poverty as main exposure
# adjust for own education (categorical), parent education (categorical), age (continuous), years in US (categorical), ACES (categorical), marital status (categorical), employment status (categorical)
poverty_global <- lm(global_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                       age_qx_mc1 + factor(ageusa_cat) + factor(aces_cat) + 
                       factor(married) + factor(work_mc1), 
                     data = cham_memory)

# tidy model and include confidence intervals in output
tidy(poverty_global, conf.int = TRUE)

# check number of observations used in linear model 
nobs(poverty_global)



#### MARGINAL STRUCTURAL MODELS ####


### Create DAG ###
# set up DAG relationships using dagify
dag <- dagify(cog_perform ~ parent_educ + maternal_educ + poverty + 
                age_nativity + ACES + home,
              cog_perform~parent_educ,
              cog_perform~maternal_educ,
              maternal_educ ~ parent_educ, 
              poverty ~ maternal_educ,
              ACES ~ parent_educ,
              parent_educ ~ age_nativity, #+ maternal_educ,
              maternal_educ ~ age_nativity + ACES, # + poverty,
              poverty~ age_nativity + ACES + home,
              home ~ parent_educ + maternal_educ,
              labels = c(cog_perform = "cognitive\n performance", 
                         parent_educ = "parent\n education",
                         maternal_educ = "own\n education",
                         poverty = "poverty\n level",
                         age_nativity = "age \n nativity",
                         ACES = "ACES", 
                         home = "marital status\n employment status\n family size\n chronic health\n conditions"),
              coords = list(x = c(parent_educ = 1, maternal_educ = 2.5,
                                  poverty = 4, cog_perform = 8, 
                                  age_nativity = 3, ACES = 3.1, home = 4),
                            y = c(parent_educ = 2, 
                                  maternal_educ = 2.5, 
                                  poverty = 2.8, 
                                  cog_perform = 2, 
                                  age_nativity = 5, ACES = 4, home = 1)),
              exposure = "parent_educ",
              outcome = "cog_perform")
dag_colors <- dag %>% tidy_dagitty() %>%
  dplyr::mutate(color = ifelse(name == "cog_perform", "outcome", 
                               ifelse(name =="parent_educ" | name == "maternal_educ" |
                                        name == "poverty", "exposure", "covariate")))

# plot DAG using ggplot
ggplot(aes(x = x, y = y, xend = xend, yend = yend), data = dag_colors) +
  geom_dag_point(aes(color = color)) +
  geom_dag_text(color = "black", size = 2, labels = c("exposure", "outcome")) + 
  geom_dag_edges_diagonal() +
  theme_dag()

# find the variables that need to be adjusted for based on the DAG
adjustmentSets(dag)


### Clean up data for marginal structural models ###

# recode parental education numerically
cham_msm <- cham_memory %>% mutate(parent_educ_dich_num = 
                                     ifelse(parent_educ_dich == "no formal or unknown", 0, 1))

# recode own education numerically
cham_msm <- cham_msm %>% mutate(own_educ_dich_num = 
                                  ifelse(educ_dich == "primary or below", 0, 1))

# recode poverty level numerically
cham_msm <- cham_msm %>% mutate(poverty_cat_num = 
                                  ifelse(pov_cat == "at or below poverty level", 0 , 
                                         ifelse(pov_cat == "above poverty level", 1, 2)))

# add row number to dataset
cham_msm <- cham_msm %>% mutate(ID = row_number())

# covariates
# ACES (mediator)
# age (time varying confounder)
# nativity (age of arrival to US - time varying confounder)
# marital status (mediator)
# employment status (mediator)
# depression (CESD) (depress_mc1 or depresscat_mc1)
# anxiety (GAD) (gadscore_mc1 or gad4cat_mc1)
# summary score of of self-rated doctor diagnosed health conditions
# number of people in household (density1_mc1)


### Select additional covariates and summarize ###

# add additional covariates from original dataset
cham_add_var <- cham %>% dplyr::select(newid, age_qx_mc1, depress_mc1, depresscat_mc1, 
                                       gadscore_mc1, gad4cat_mc1,
                                       diab_mc1, hbp_mc1, chol_mc1, heart_mc1, cancer_mc1, 
                                       asth_mc1, thyr_mc1, dep_mc1, anx_mc1, density1_mc1, 
                                       bpcat_mc1, trigcat_mc1, trig_mc1, gluccat_mc1, fast_bld_mc1,
                                       hba1ccat_mc1, cancer_m18y, cancersp_m18y, masth_m16y, 
                                       masthmed_m16y, depresscat_m16y, depresscat_m18y, 
                                       gad7_4cat_m18y, gad7_4cat_m16y, tmtb_fail_mc1)

# calculate summary score of self-rated doctor diagnosed health conditions
# diabetes
# high blood pressure 
# heart problems
# cancer
# asthma 
# thyroid problems

# re-code don't knows to 999
cham_add_var <- cham_add_var %>% 
  mutate(diab_mc1 = ifelse(diab_mc1 == 9, 999, diab_mc1),
         hbp_mc1 = ifelse(hbp_mc1 == 9, 999, hbp_mc1),
         # chol_mc1 = ifelse(chol_mc1 == 9, 999, chol_mc1),
         heart_mc1 = ifelse(heart_mc1 == 9, 999, heart_mc1), 
         cancer_mc1 = ifelse(cancer_mc1 == 9, 999, cancer_mc1), 
         asth_mc1 = ifelse(asth_mc1 == 9, 999, asth_mc1), 
         thyr_mc1 = ifelse(thyr_mc1 == 9, 999, thyr_mc1)) 


# check distribution of each self-rated doctor diagnosed health condition
cham_add_var %>% group_by(diab_mc1) %>% count() # only 1 person missing diabetes score
cham_add_var %>% group_by(hbp_mc1) %>% count() # 3 missing HBP score
cham_add_var %>% group_by(chol_mc1) %>% count() # 6 missing cholesterol
cham_add_var %>% group_by(heart_mc1) %>% count() # 4 missing heart issues
cham_add_var %>% group_by(cancer_mc1) %>% count() # 4 missing cancer
cham_add_var %>% group_by(asth_mc1) %>% count() # 2 missing asthma
cham_add_var %>% group_by(thyr_mc1) %>% count() # 0 missing thyroid
cham_add_var %>% group_by(dep_mc1) %>% count() # 7 missing depression 
cham_add_var %>% group_by(anx_mc1) %>% count() # 9 missing anxiety


# For those who don't know high blood pressure: 
# Check bpcat_mc1 variable for blood pressure information
# Replace the don't knows with 1 if have high blood pressure, 0 if no HBP
# only don't have hypertension if 0 for hbp_mc1 and bpcat_mc1 is 1 or 2 
cham_add_var %>% filter(hbp_mc1 == 999)
cham_add_var %>% filter(hbp_mc1 == 0 & bpcat_mc1 %in% c(3, 4))

recat_bp <- cham_add_var %>% filter(hbp_mc1 == 0 & bpcat_mc1 %in% c(3, 4))

cham_add_var <- cham_add_var %>% mutate(hbp_mc1 = ifelse(newid %in% c(recat_bp$newid), 1, 
                                                         hbp_mc1))

cham_add_var <- cham_add_var %>% mutate(hbp_mc1 = ifelse(newid == 99 & bpcat_mc1 == 3, 1,
                                                         ifelse(newid == 102 & bpcat_mc1 == 1, 0,
                                                                ifelse(newid == 617 & bpcat_mc1 == 1, 0,
                                                                       hbp_mc1))))


# check that there are no more "don't know" responses for high blood pressure
cham_add_var %>% group_by(hbp_mc1) %>% count()
cham_add_var %>% filter(hbp_mc1 == 0 & is.na(bpcat_mc1))


# For those who don't know diabetes:
# Make sure they had fasting blood taken and check gluccat_mc1 variable for diabetes information
# Replace the don't knows with 1 if diabetes and 0 if no diabetes/pre-diabetes (if fasting)
# if not fasting, check hba1c score 

cham_add_var %>% filter(diab_mc1 == 999)

# check: among those who self report no diabetes, how many do not have a fasting blood
# sugar score and have a hba1ccat_mc1 score that show diabetes
cham_add_var %>% filter(diab_mc1 == 0 & fast_bld_mc1 == 0 & hba1ccat_mc1 == 4)

# check: among those who report no diabetes, how many do have a fasting blood sugar
# and blood sugar is in diabetes category
recat_diab <- cham_add_var %>% filter(diab_mc1 == 0 & fast_bld_mc1 == 1 & gluccat_mc1 == 3)

cham_add_var <- cham_add_var %>% mutate(diab_mc1 = ifelse(newid %in% c(recat_diab$newid), 1, 
                                                          diab_mc1))

cham_add_var <- cham_add_var %>% mutate(diab_mc1 = 
                                          ifelse(newid == 222 & gluccat_mc1 == 2, 0, diab_mc1))

cham_add_var %>% group_by(diab_mc1) %>% count()

cham_add_var %>% filter(diab_mc1 == 0 & is.na(fast_bld_mc1))

# 4 missing heart issues - can't adjust b/c other heart variables are from 18 year and show 
# no heart issues - can't be sure they didn't develop heart issues between then and now 


# For those who don't know cancer: 
# Check cancer_m18y variable
# Replace don't knows with 1 if ever had cancer 
cham_add_var %>% filter(cancer_mc1 == 999)

cham_add_var <- cham_add_var %>% mutate(cancer_mc1 = 
                                          ifelse(newid == 610 & cancer_m18y == 1, 1, cancer_mc1))

cham_add_var %>% group_by(cancer_mc1) %>% count()
# 3 left over missing cancer info

# 2 missing asthma - can't adjust b/c other asthma vars from 16y follow-up 
# show no asthma - can't be sure they didn't get asthma between then and now


# filter to only include mat cog participants to check distribution of self-rated health scores
cham_add_var <- cham_add_var %>% filter(!is.na(age_qx_mc1))

# calculate overall self-rated health score
cham_add_var <- cham_add_var %>% mutate(summary_health = diab_mc1 + hbp_mc1)


# check distribution of scores
cham_add_var %>% group_by(summary_health) %>% count()


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
cham_add_var <- cham_add_var %>% mutate(summary_health_cat_yes = 
                                          ifelse(summary_health == 0 | 
                                                   summary_health == 999, "no health issues", 
                                                 ifelse(summary_health == 1 |
                                                          summary_health == 1000, "1 health issue", 
                                                        ifelse(summary_health > 1 & summary_health <= 7| 
                                                                 summary_health %in% c(1001, 1002, 2000), "2+",
                                                               NA))))


# also try grouping by "don't knows"
# No major health issues: 0
# 1 health issue: 1
# 2-5 health issues: 2-5
# At least one don't know: 999, 1000, 1001, 1002, 2000
# NA - missing
cham_add_var <- cham_add_var %>% mutate(summary_health_cat_dk = 
                                          ifelse(summary_health == 0, "no health issues", 
                                                 ifelse(summary_health == 1, "1 health issue", 
                                                        ifelse(summary_health > 1 & summary_health <= 7, "2+",
                                                               ifelse(summary_health >=999, "at least 1 unknown", NA)))))


# Check distribution of categorized self-rated scores
cham_add_var %>% group_by(summary_health_cat_yes) %>% count();
cham_add_var %>% group_by(summary_health_cat_dk) %>% count()

# select variables to use in IPW calculations and filter to include only matcog participants
cham_other_vars <- cham_add_var %>% filter(!is.na(age_qx_mc1)) %>% 
  dplyr::select(depress_mc1, gadscore_mc1, 
                summary_health_cat_yes, summary_health_cat_dk,
                density1_mc1, depresscat_mc1, gad4cat_mc1, diab_mc1, hbp_mc1, heart_mc1,
                tmtb_fail_mc1) 

# combine original variables with additional variables in single dataset
cham_msm <- cbind(cham_msm, cham_other_vars)
cham_msm_full <- cham_msm

no_miss_cov <- cham_msm %>% filter(!(is.na(age_qx_mc1) | is.na(ageusa_cat) |
                                       is.na(aces_cat) | is.na(married) | is.na(work_mc1))) 


# create table of dataset not missing any covariates used in model
no_miss_cov_tab <- cham_msm %>% filter(!(is.na(age_qx_mc1) | is.na(ageusa_cat) |
                                           is.na(aces_cat) |
                                           is.na(married) | is.na(work_mc1))) 

cham_msm %>% filter(is.na(aces_cat)) 

cham_msm %>% group_by(aces) %>% count()


no_miss_cov_tab %>% summarize(mean = mean(age_qx_mc1), 
                              sd = sd(age_qx_mc1))

no_miss_cov_tab %>% group_by(ageusa_cat) %>% count() 

no_miss_cov_tab %>% group_by(lang_exam_mc1) %>% count() 

no_miss_cov_tab %>% group_by(aces_cat) %>% count()

no_miss_cov_tab %>% group_by(married) %>% count()

no_miss_cov_tab %>% group_by(work_mc1) %>% count()

no_miss_cov_tab %>% group_by(summary_health_cat_yes) %>% count()

no_miss_cov_tab %>% group_by(diab_mc1) %>% count()

no_miss_cov_tab %>% group_by(hbp_mc1) %>% count()

no_miss_cov_tab %>% summarize(mean = mean(density1_mc1), 
                              sd = sd(density1_mc1))

no_miss_cov_tab %>% group_by(parent_educ) %>% count()

no_miss_cov_tab %>% group_by(educcat_mom) %>% count()

no_miss_cov_tab %>% group_by(povcat_mc1) %>% count()

# calculate remaining descriptive statistics for table 1 
# including those missing covariate values
cham_msm_full %>% group_by(depresscat_mc1) %>% count()
cham_msm_full %>% group_by(gad4cat_mc1) %>% count()
cham_msm_full %>% group_by(summary_health_cat_yes) %>% count()
cham_msm_full %>% filter(!is.na(density1_mc1)) %>% summarize(mean = mean(density1_mc1), 
                                                             sd = sd(density1_mc1))


# Descriptive lifecourse SES table using analytic sample
# dichotomize parental education to any education or no formal education/unknown
cham_desc_analytic <- no_miss_cov %>% mutate(parent_educ_dich = 
                                               ifelse(parent_educ == "any formal educ", 
                                                      "any formal educ", 
                                                      "no formal or unknown"))

# create dichotomized poverty variable
# At or below poverty level and above poverty level (>200% poverty and Poverty-200%)
cham_desc_analytic <- cham_desc_analytic %>% mutate(pov_dich = 
                                                      ifelse(povcat_mc1 == "At or below poverty", 
                                                             "at or below poverty level", 
                                                             ifelse(povcat_mc1 == ">=200% FPL" |
                                                                      povcat_mc1 == "200% FPL", 
                                                                    "above poverty level", 
                                                                    NA)))

# dichotomize education variable
# Primary education or below (<=6th grade) and above primary (7-12th, more than HS)
cham_desc_analytic <- cham_desc_analytic %>% mutate(educ_dich = 
                                                      ifelse(educcat_mom == "<=6th grade", 
                                                             "primary or below", "above primary"))

# check contingency table between dichotomized poverty level and dichotomized education
table(cham_desc_analytic$educ_dich, cham_desc_analytic$pov_dich, useNA = "always")

# create variable that combines poverty and education
cham_desc_analytic <- cham_desc_analytic %>% mutate(pov_educ = 
                                                      ifelse(pov_dich == "at or below poverty level" &
                                                               educ_dich == "primary or below", 
                                                             "primary or below + at or below poverty level",
                                                             ifelse(pov_dich == "above poverty level" &
                                                                      educ_dich == "primary or below", 
                                                                    "primary or below + above poverty level",
                                                                    ifelse(pov_dich == "at or below poverty level" &
                                                                             educ_dich == "above primary", 
                                                                           "above primary + at or below poverty level", 
                                                                           ifelse(pov_dich == "above poverty level" &
                                                                                    educ_dich == "above primary",
                                                                                  "above primary + above poverty level", NA)))))

cham_desc_analytic %>% group_by(pov_educ) %>% count()

table(cham_desc_analytic$pov_educ, cham_desc_analytic$parent_educ_dich, useNA = "always")



#### CALCULATE INVERSE PROBABILITY WEIGHTS ####


### Weights for own education ###

# use previous exposure (parent education) to calculate numerator
model_num_own_educ <- 1

# denominator using parent education, age arrived to US, ACEs, and age
model_denom_own_educ <- glm(own_educ_dich_num ~ parent_educ_dich_num + 
                            ageusa_cat + aces + age_qx_mc1 + 
                            lang_exam_mc1,
                            data = cham_msm, 
                            family = binomial(link = "logit"))

# denominator using parent education, age arrived to US, ACEs, age, and squared age term
model_denom_own_educ_quad <- glm(own_educ_dich_num ~ parent_educ_dich_num + 
                                 ageusa_cat + aces + age_qx_mc1 + 
                                 lang_exam_mc1 + I(age_qx_mc1^2),
                                 data = cham_msm, 
                                 family = binomial(link = "logit"))
# not significantly different so don't need to include squared age term

# denominator using parent education, age arrived to US, ACEs, age, and 
# interactions between all
model_denom_own_educ_quad <- glm(own_educ_dich_num ~ parent_educ_dich_num + 
                                   ageusa_cat + aces + age_qx_mc1 +
                                   parent_educ_dich_num*ageusa_cat +
                                   parent_educ_dich_num*aces +
                                   parent_educ_dich_num*age_qx_mc1 +
                                   ageusa_cat*aces +
                                   ageusa_cat*age_qx_mc1 +
                                   aces*age_qx_mc1,
                                 data = cham_msm, 
                                 family = binomial(link = "logit"))


# predict outcomes using models
pred_denom_own_educ <- predict(model_denom_own_educ, type = "response")
pred_denom_own_educ_quad <- predict(model_denom_own_educ_quad, type = "response")

# plot residuals against model predictions to see if quadratic shape
res <- residuals(model_denom_own_educ, type = 'deviance')
plot(pred_denom_own_educ,res)


# calculate IPW for own education using predictions
# filter out records that have NAs for aces (only variable that has NAs) to calculate weights properly
cham_msm_own_educ <- cham_msm %>% filter(!is.na(aces)) %>%
  mutate(ipw_own_educ = ifelse(own_educ_dich_num == 0, 
                               1/(1-pred_denom_own_educ), 1/pred_denom_own_educ))


### Weights for poverty ###
model_num_poverty <- 1

# filter out observations that are missing any covariates (only covariates with missings are density1_mc1)
cham_msm_pov <- cham_msm_own_educ %>% filter(!is.na(lang_exam_mc1))  

# update working in agriculture variable
cham_msm_pov$wkag_mc1 <- as.numeric(cham_msm_pov$wkag_mc1)
cham_msm_pov <- cham_msm_pov %>% mutate(wkag_mc1 = ifelse(is.na(wkag_mc1), "Did not work", 
                                                          ifelse(wkag_mc1 == 0, "No", "Yes")))

cham_msm_pov %>% group_by(wkag_mc1) %>% count()


### Weights for poverty variable including missings as its own category ###

# denominator using own education, parent education, age arrived to US, ACEs, age, 
# married status, and working since last visit status
model_denom_poverty <- multinom(pov_cat ~ own_educ_dich_num + 
                                parent_educ_dich_num + age_qx_mc1 + ageusa_cat +
                                lang_exam_mc1 +
                                aces + married + 
                                work_mc1, 
                                data = cham_msm_pov)

cham_msm_pov %>% group_by(pov_cat) %>% count()

# predict outcomes using model
pred_denom_poverty <- predict(model_denom_poverty, type = "probs")

# create new dataframe with prediction probabilities
m1 <- pred_denom_poverty[,1]
m2 <- pred_denom_poverty[,2]
m3 <- pred_denom_poverty[,3]

pov_pred <- as.data.frame(cbind(m1, m2, m3))

# calculate inverse probabilities as 1/probability
pov_pred <- pov_pred %>% mutate(p1 = 1/m1, p2 = 1/m2, p3 = 1/m3)


# calculate IPW for poverty using predictions
cham_msm_pov$pov_weight <- ifelse(
  cham_msm_pov$pov_cat == "at or below poverty level", pov_pred$p1, 
  ifelse(cham_msm_pov$pov_cat == "above poverty level", pov_pred$p2, 
         ifelse(cham_msm_pov$pov_cat == "Missing", pov_pred$p3, NA)))

# check weights 
cbind(cham_msm_pov, pov_pred) %>% dplyr::select(pov_cat, pov_weight, m1, m2, m3, p1, p2, p3)

## weights for poverty variable re-categorizing missings

# denominator using own education, parent education, years in the US, ACEs, age, 
# married status, working since last visit status

# re-categorize those with missing poverty information as part of 
# at or below poverty level
cham_msm_no_missing_pov <- cham_msm_pov %>% 
  mutate(pov_bin = 
           ifelse(pov_cat == "Missing" |
                    pov_cat == "at or below poverty level", 0, 
                  1))

model_denom_poverty_no_miss <- glm(pov_bin ~ own_educ_dich_num + 
                                   parent_educ_dich_num + age_qx_mc1 + ageusa_cat +
                                   lang_exam_mc1 + 
                                   aces + married + 
                                   work_mc1, 
                                   data = cham_msm_no_missing_pov, 
                                   family = binomial(link = "logit"))

# predict outcomes using model
pred_denom_poverty_no_miss <- predict(model_denom_poverty_no_miss, type = "response")

# calculate IPW for poverty using predictions
cham_msm_no_missing_pov <- cham_msm_no_missing_pov %>% 
  mutate(ipw_poverty_no_miss = ifelse(pov_bin == 0, 
                                      1/(1-pred_denom_poverty_no_miss),
                                      1/pred_denom_poverty_no_miss))

# test whether weights are significantly different between poverty variable including missings
# as its own category and poverty variable categorizing missings as part of 
# at or below poverty level
t.test(cham_msm_no_missing_pov$ipw_poverty_no_miss, cham_msm_pov$pov_weight, paired = T)

# there is a significant difference - proceed with using model that includes poverty variable that includes missings



## weights for poverty variable using different coding of self-rated health score
# code based on "don't knows"

# Denominator using own education, parent education, years in the US, ACEs, age, 
# married status, working since last visit status

model_denom_poverty_dk <- multinom(pov_cat ~ own_educ_dich_num + 
                                   parent_educ_dich_num + age_qx_mc1 + ageusa_cat +
                                   lang_exam_mc1 +
                                   aces + married + 
                                   work_mc1, 
                                   data = cham_msm_own_educ)

# predict outcomes using model with different coding of self rated health score
pred_denom_poverty_dk <- predict(model_denom_poverty_dk, type = "probs")

# create new dataframe with prediction probabilities
m1_dk <- pred_denom_poverty_dk[,1]
m2_dk <- pred_denom_poverty_dk[,2]
m3_dk <- pred_denom_poverty_dk[,3]

pov_pred_dk <- as.data.frame(cbind(m1_dk, m2_dk, m3_dk))

# calculate inverse probabilities as 1/probability
pov_pred_dk <- pov_pred_dk %>% mutate(p1_dk = 1/m1_dk, p2_dk = 1/m2_dk, 
                                      p3_dk= 1/m3_dk)


# calculate IPW for poverty using predictions
cham_msm_pov_dk <- cham_msm_own_educ 

cham_msm_pov_dk$pov_weight_dk <- ifelse(
  cham_msm_pov_dk$pov_cat == "at or below poverty level", pov_pred_dk$p1_dk, 
  ifelse(cham_msm_pov_dk$pov_cat == "above poverty level", pov_pred_dk$p2_dk, 
         ifelse(cham_msm_pov_dk$pov_cat == "Missing", pov_pred_dk$p3_dk, NA)))

# plot residuals against model predictions to see if quadratic shape
res <- residuals(model_denom_poverty_dk, type = 'deviance')
plot(pred_denom_poverty_dk,res)

# ANOVA LRT to test if there is a difference in the models
#anova(model_denom_poverty, model_denom_poverty_dk)

# test whether weights are significantly different when using different scoring of 
# summary self-rated health score
#t.test(cham_msm_pov_dk$pov_weight_dk, cham_msm_pov$pov_weight, paired = T)

# significant difference based on weights - proceed with model that uses
# self-rated health score that codes based on yeses

# check how many individuals failed to complete Trails B because they did not attempt it (discontinued), they reached 5 errors, or they reached the time limit of 240 seconds.
cham_msm_pov %>% group_by(tmtb_fail_mc1) %>% count()

cham_msm_pov %>% filter(is.na(global_rc_mc1) | is.na(execfun_rc_mc1))


### Calculate multiplied weights ###

# multiply poverty weight and own education weight
cham_msm_pov <- cham_msm_pov %>% mutate(ipw_combined = ipw_own_educ * pov_weight)

cham_msm_pov_dk <- cham_msm_pov_dk %>% mutate(ipw_combined_dk = ipw_own_educ * pov_weight_dk)


# compare distributions of combined weights when using model that uses 
# health summary score grouped by "yeses" compared to health summary score grouped by "don't knows"
ggplot(data = cham_msm_pov, aes(x = ipw_combined)) + 
  geom_histogram()

ggplot(data = cham_msm_pov_dk, aes(x = ipw_combined_dk)) + 
  geom_histogram()

cham_msm_pov %>% summarize(min = min(ipw_combined), 
                           mean = mean(ipw_combined),
                           q1 = quantile(ipw_combined, .25),
                           median = median(ipw_combined), 
                           q3 = quantile(ipw_combined, .75), 
                           max = max(ipw_combined)); 
cham_msm_pov_dk %>% summarize(min = min(ipw_combined_dk), 
                              mean = mean(ipw_combined_dk),
                              q1 = quantile(ipw_combined_dk, .25),
                              median = median(ipw_combined_dk), 
                              q3 = quantile(ipw_combined_dk, .75), 
                              max = max(ipw_combined_dk))

# distribution of combined weights is similar - proceed with model that uses summary
# health score with fewer categories (categorized by yeses)


## trim weights at 99th percentile
cut_off <- quantile(cham_msm_pov$ipw_combined, .99)

cham_msm_pov <- cham_msm_pov %>% 
  mutate(ipw_combined = ifelse(ipw_combined > cut_off, cut_off, ipw_combined))

ggplot(data = cham_msm_pov, aes(x = ipw_combined)) + 
  geom_histogram()



#### MSM WITH MEMORY SCORE OUTCOME ####


### Memory: total effect estimate of parental education ###

# no weights
# run linear regression
msm_memory_parent <- lm(memory_mc1 ~ parent_educ + age_qx_mc1 +
                        ageusa_cat + lang_exam_mc1,
                        data = cham_msm_pov)

# summarize regression output
summary(msm_memory_parent)

# calculate 95% CIs
beta <- coef(msm_memory_parent)
SE <- coef(summary(msm_memory_parent))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
memory_parent_dat <- as.data.frame(cbind(beta, lcl, ucl))
memory_parent_dat <- memory_parent_dat %>% mutate(type = "Memory")
memory_parent_dat <- memory_parent_dat %>% mutate(model = "Model 1")

# calculate number of observations used in linear regression
nobs(msm_memory_parent)


### Memory: total effect estimate of own education and direct effect of parental education ### 

# use weights for own education

# run linear regression
msm_memory_own_educ <- lm(memory_mc1 ~ educcat_mom + parent_educ + 
                          age_qx_mc1 + ageusa_cat +
                          lang_exam_mc1,
                          data = cham_msm_pov,
                          weights = ipw_own_educ)


# summarize regression output
summary(msm_memory_own_educ)

# calculate 95% CIs
beta <- coef(msm_memory_own_educ)
SE <- coef(summary(msm_memory_own_educ))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
memory_educ_dat <- as.data.frame(cbind(beta, lcl, ucl))
memory_educ_dat <- memory_educ_dat %>% mutate(type = "Memory")
memory_educ_dat <- memory_educ_dat %>% mutate(model = "Model 2")

# calculate number of observations used in linear regression
nobs(msm_memory_own_educ)


### Memory: total effect estimate of poverty and direct effect estimate of parental and own education ###

# use multiplied weights

# run linear regression
msm_memory_poverty <- lm(memory_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                         age_qx_mc1 + ageusa_cat +
                         lang_exam_mc1, 
                         data = cham_msm_pov, 
                         weights = ipw_combined)

# summarize regression output
summary(msm_memory_poverty)

# calculate 95% CIs
beta <- coef(msm_memory_poverty)
SE <- coef(summary(msm_memory_poverty))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
memory_pov_dat <- as.data.frame(cbind(beta, lcl, ucl))
memory_pov_dat <- memory_pov_dat %>% mutate(type = "Memory")
memory_pov_dat <- memory_pov_dat %>% mutate(model = "Model 3")

# calculate number of observations used in linear regression model
nobs(msm_memory_poverty)



### MSM WITH EXECUTIVE FUNCTION SCORE OUTCOME #### 


### Executive function: total effect estimate of parental education ###

# no weights
# run linear regression
msm_exec_parent <- lm(execfun_rc_mc1 ~ parent_educ + age_qx_mc1 +
                      ageusa_cat + 
                      lang_exam_mc1,
                      data = cham_msm_pov)

# summarize regression output
summary(msm_exec_parent)

# calculate 95% CIs
beta <- coef(msm_exec_parent)
SE <- coef(summary(msm_exec_parent))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
exec_parent_dat <- as.data.frame(cbind(beta, lcl, ucl))
exec_parent_dat <- exec_parent_dat %>% mutate(type = "Executive Function")
exec_parent_dat <- exec_parent_dat %>% mutate(model = "Model 1")

# calculate number of observations used in linear regression model
nobs(msm_exec_parent)


### Executive function: total effect estimate of own education and direct effect of parental education ###

# use weights for own education

# run linear regression
# use dataset that removed NAs from ACES variable
msm_exec_own_educ <- lm(execfun_rc_mc1 ~ educcat_mom + parent_educ + 
                        age_qx_mc1 + ageusa_cat +
                        lang_exam_mc1, 
                        data = cham_msm_pov,
                        weights = ipw_own_educ)

# summarize regression output
summary(msm_exec_own_educ)

# calculate 95% CIs
beta <- coef(msm_exec_own_educ)
SE <- coef(summary(msm_exec_own_educ))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
exec_educ_dat <- as.data.frame(cbind(beta, lcl, ucl))
exec_educ_dat <- exec_educ_dat %>% mutate(type = "Executive Function")
exec_educ_dat <- exec_educ_dat %>% mutate(model = "Model 2")

# calculate number of observations used in linear regression model
nobs(msm_exec_own_educ)


### Executive function: total effect estimate of poverty and direct effect estimate of parental education and own education ###

# use multiplied weights
# run linear regression
msm_exec_poverty <- lm(execfun_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                       age_qx_mc1 + ageusa_cat + 
                       lang_exam_mc1, 
                       data = cham_msm_pov, 
                       weights = ipw_combined)

# summarize regression output
summary(msm_exec_poverty)

# calculate 95% CIs
beta <- coef(msm_exec_poverty)
SE <- coef(summary(msm_exec_poverty))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
exec_pov_dat <- as.data.frame(cbind(beta, lcl, ucl))
exec_pov_dat <- exec_pov_dat %>% mutate(type = "Executive Function")
exec_pov_dat <- exec_pov_dat %>% mutate(model = "Model 3")

# calculate number of observations used in linear regression model
nobs(msm_exec_poverty)



#### MSM WITH VERBAL FLUENCY OUTCOME ####


### Verbal fluency: total effect estimate of parental education ###

# no weights
# run linear regression
msm_verbal_parent <- lm(verbal_rc_mc1 ~ parent_educ + age_qx_mc1 +
                        ageusa_cat + 
                        lang_exam_mc1,
                        data = cham_msm_pov)

# summarize regression output
summary(msm_verbal_parent)

# calculate 95% CIs
beta <- coef(msm_verbal_parent)
SE <- coef(summary(msm_verbal_parent))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
verbal_parent_dat <- as.data.frame(cbind(beta, lcl, ucl))
verbal_parent_dat <- verbal_parent_dat %>% mutate(type = "Verbal Fluency")
verbal_parent_dat <- verbal_parent_dat %>% mutate(model = "Model 1")

# calculate number of observations used in model
nobs(msm_verbal_parent)


### Verbal fluency: total effect estimate of own education and direct effect estimate of parental education ###

# use weights for own education

# run linear regression
# use dataset that removed NAs from ACES variable
msm_verbal_own_educ <- lm(verbal_rc_mc1 ~ educcat_mom + parent_educ + 
                          age_qx_mc1 + ageusa_cat +
                          lang_exam_mc1, 
                          data = cham_msm_pov,
                          weights = ipw_own_educ)

# summarize regression output
summary(msm_verbal_own_educ)

# calculate 95% CIs
beta <- coef(msm_verbal_own_educ)
SE <- coef(summary(msm_verbal_own_educ))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
verbal_educ_dat <- as.data.frame(cbind(beta, lcl, ucl))
verbal_educ_dat <- verbal_educ_dat %>% mutate(type = "Verbal Fluency")
verbal_educ_dat <- verbal_educ_dat %>% mutate(model = "Model 2")

# calculate number of observations used in model
nobs(msm_verbal_own_educ)


### Verbal fluency: total effect estimate of poverty and direct effect estimate of parental education and own education ###

# use multiplied weights
# run linear regression
msm_verbal_poverty <- lm(verbal_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                         age_qx_mc1 + ageusa_cat +
                         lang_exam_mc1, 
                         data = cham_msm_pov, 
                         weights = ipw_combined)

# summarize regression output
summary(msm_verbal_poverty)

# calculate 95% CIs
beta <- coef(msm_verbal_poverty)
SE <- coef(summary(msm_verbal_poverty))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
verbal_pov_dat <- as.data.frame(cbind(beta, lcl, ucl))
verbal_pov_dat <- verbal_pov_dat %>% mutate(type = "Verbal Fluency")
verbal_pov_dat <- verbal_pov_dat %>% mutate(model = "Model 3")

# calculate number of observations used in model
nobs(msm_verbal_poverty)



#### MSM WITH GLOBAL FUNCTION OUTCOME ####


### Global function: total effect estimate of parental education ###

# no weights
# run linear regression
msm_global_parent <- lm(global_rc_mc1 ~ parent_educ + age_qx_mc1 +
                        ageusa_cat +
                        lang_exam_mc1,
                        data = cham_msm_pov)

# summarize regression output
summary(msm_global_parent)

# calculate 95% CIs
beta <- coef(msm_global_parent)
SE <- coef(summary(msm_global_parent))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
global_parent_dat <- as.data.frame(cbind(beta, lcl, ucl))
global_parent_dat <- global_parent_dat %>% mutate(type = "Global Function")
global_parent_dat <- global_parent_dat %>% mutate(model = "Model 1")

# calculate number of observations used in model
nobs(msm_global_parent)


### Global function: total effect estimate of own education and direct effect estimate of parental education ###

# use weights for own education

# run linear regression
# use dataset that removed NAs from ACES variable
msm_global_own_educ <- lm(global_rc_mc1 ~ educcat_mom + parent_educ + 
                          age_qx_mc1 + ageusa_cat + 
                          lang_exam_mc1,
                          data = cham_msm_pov,
                          weights = ipw_own_educ)

# summarize regression output
summary(msm_global_own_educ)

# calculate 95% CIs
beta <- coef(msm_global_own_educ)
SE <- coef(summary(msm_global_own_educ))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
global_educ_dat <- as.data.frame(cbind(beta, lcl, ucl))
global_educ_dat <- global_educ_dat %>% mutate(type = "Global Function")
global_educ_dat <- global_educ_dat %>% mutate(model = "Model 2")

# calculate number of observations used in linear regression model
nobs(msm_global_own_educ)


### Global function: total effect estimate of poverty and direct effect estimate of parental education and own education ####

# use multiplied weights
# run linear regression
msm_global_poverty <- lm(global_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                         age_qx_mc1 + ageusa_cat +
                         lang_exam_mc1, 
                         data = cham_msm_pov, 
                         weights = ipw_combined)

# summarize regression output
summary(msm_global_poverty)

# calculate 95% CIs
beta <- coef(msm_global_poverty)
SE <- coef(summary(msm_global_poverty))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
global_pov_dat <- as.data.frame(cbind(beta, lcl, ucl))
global_pov_dat <- global_pov_dat %>% mutate(type = "Global Function")
global_pov_dat <- global_pov_dat %>% mutate(model = "Model 3")

# calculate number of observations used in model
nobs(msm_global_poverty)


### Combine all results into errorbar plot ###

# combine all MSM results for models 1-3 and outcomes 1-4 into single dataframe
all_results <- rbind(memory_parent_dat, memory_educ_dat, memory_pov_dat, exec_parent_dat, 
                     exec_educ_dat, exec_pov_dat, verbal_parent_dat, verbal_educ_dat, 
                     verbal_pov_dat, global_parent_dat, global_educ_dat, global_pov_dat)

# create column name for the index column
all_results <- tibble::rownames_to_column(all_results, "exposure") 

# include relevant rows for errorbar plot
all_results <- all_results %>% filter(exposure %in% c("parent_educDont know/dont remember", 
                                                      "parent_educany formal educ", 
                                                      "educcat_mom7-11th grade", 
                                                      "educcat_mom>= HS graduate", 
                                                      "parent_educDont know/dont remember1", 
                                                      "parent_educany formal educ1", 
                                                      "pov_catabove poverty level", 
                                                      "pov_catMissing", 
                                                      "educcat_mom7-11th grade1", 
                                                      "educcat_mom>= HS graduate1", 
                                                      "parent_educDont know/dont remember2",
                                                      "parent_educany formal educ2",
                                                      "parent_educDont know/dont remember3",
                                                      "parent_educany formal educ3",
                                                      "educcat_mom7-11th grade2",
                                                      "educcat_mom>= HS graduate2",
                                                      "parent_educDont know/dont remember4",
                                                      "parent_educany formal educ4",
                                                      "pov_catabove poverty level1",
                                                      "pov_catMissing1", 
                                                      "educcat_mom7-11th grade3", 
                                                      "educcat_mom>= HS graduate3",
                                                      "parent_educDont know/dont remember5",
                                                      "parent_educany formal educ5",
                                                      "parent_educDont know/dont remember6",
                                                      "parent_educany formal educ6",
                                                      "educcat_mom7-11th grade4",
                                                      "educcat_mom>= HS graduate4",
                                                      "parent_educDont know/dont remember7",
                                                      "parent_educany formal educ7",
                                                      "pov_catabove poverty level2",
                                                      "pov_catMissing2",
                                                      "educcat_mom7-11th grade5",
                                                      "educcat_mom>= HS graduate5",
                                                      "parent_educDont know/dont remember8",
                                                      "parent_educany formal educ8",
                                                      "parent_educDont know/dont remember9",
                                                      "parent_educany formal educ9",
                                                      "educcat_mom7-11th grade6",
                                                      "educcat_mom>= HS graduate6",
                                                      "parent_educDont know/dont remember10",
                                                      "parent_educany formal educ10",
                                                      "pov_catabove poverty level3",
                                                      "pov_catMissing3",
                                                      "educcat_mom7-11th grade7",
                                                      "educcat_mom>= HS graduate7",
                                                      "parent_educDont know/dont remember11",
                                                      "parent_educany formal educ11"))

# rename exposures so they match between all models
all_results <- all_results %>% mutate(exposure = 
                                        ifelse(grepl("parent_educDont know/dont remember",
                                                     exposure),
                                               "parent_educDont know/dont remember",
                                               ifelse(grepl("parent_educany formal educ", exposure),
                                                      "parent_educany formal educ", 
                                                      ifelse(grepl("educcat_mom7-11th grade", exposure),
                                                             "educcat_mom7-11th grade", 
                                                             ifelse(grepl("educcat_mom>= HS graduate", exposure),
                                                                    "educcat_mom>= HS graduate", 
                                                                    ifelse(grepl("pov_catabove poverty level", exposure),
                                                                           "pov_catabove poverty level", 
                                                                           ifelse(grepl("pov_catMissing", exposure),
                                                                                  "pov_catMissing", NA)))))))

all_results <- all_results %>% 
  mutate(exposure = 
           ifelse(exposure == "parent_educDont know/dont remember",
                  "Parent Education: Don't know/don't remember: None", 
                  ifelse(exposure == "parent_educany formal educ",
                         "Parent Education: Any formal education: None",
                         ifelse(exposure == "educcat_mom7-11th grade", 
                                "Own Education: 7-11th grade: 6th grade or below", 
                                ifelse(exposure == "educcat_mom>= HS graduate", 
                                       "Own Education: >= HS: 6th grade or below",
                                       ifelse(exposure == "pov_catabove poverty level",
                                              "Poverty: Above FPL: At or below FPL", 
                                              ifelse(exposure == "pov_catMissing", 
                                                     "Poverty: Missing: At or below FPL", NA)))))))

# reorder exposure variable for plots
all_results <- all_results %>% mutate(exposure = 
                                        fct_relevel(exposure,
                                                    "Poverty: Above FPL: At or below FPL",
                                                    "Poverty: Missing: At or below FPL",
                                                    "Own Education: >= HS: 6th grade or below",
                                                    "Own Education: 7-11th grade: 6th grade or below",
                                                    "Parent Education: Don't know/don't remember: None",
                                                    "Parent Education: Any formal education: None"))

# model 3 only 
model3 <- all_results %>% filter(model == "Model 3")

ggplot(data = model3, aes(x = beta, y = exposure, color = type)) + 
  geom_point() +
  facet_wrap(~factor(type, levels = c("Memory", "Executive Function", 
                                      "Verbal Fluency", "Global Function"))) +
  geom_errorbar(aes(xmin = lcl, xmax = ucl)) +
  theme_minimal(base_size = 15) +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=0.5) +
  theme(axis.text = element_text(size = 5)) +
  ylab("") +
  xlab("Z-Score") +
  theme(axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 10)) +
  theme(legend.position = "none")


# excluding missing categories 
no_miss <- model3 %>% filter(exposure != "Poverty: Missing: At or below FPL" &
                               exposure != "Parent Education: Don't know/don't remember: None")

ggplot(data = no_miss, aes(x = beta, y = exposure, color = type)) + 
  geom_point() +
  facet_wrap(~factor(type, levels = c("Memory", "Executive Function", 
                                      "Verbal Fluency", "Global Function"))) +
  geom_errorbar(aes(xmin = lcl, xmax = ucl)) +
  theme_minimal(base_size = 15) +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=0.5) +
  theme(axis.text = element_text(size = 5)) +
  ylab("") +
  xlab("Z-Score") +
  theme(axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 10)) +
  theme(legend.position = "none")


# global function only
global_func <- all_results %>% filter(type == "Global Function")

ggplot(data = global_func, aes(x = beta, y = exposure, color = model)) + 
  geom_point() +
  facet_wrap(vars(model)) +
  geom_errorbar(aes(xmin = lcl, xmax = ucl)) +
  theme_minimal(base_size = 15) +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=0.5) +
  theme(axis.text = element_text(size = 5)) +
  ylab("") +
  xlab("Z-Score") + 
  theme(axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 10)) + 
  theme(legend.position = "none")


# memory, verbal fluency, executive function
mve <- all_results %>% filter(type != "Global Function")

ggplot(data = mve, aes(x = beta, y = exposure, color = type)) + 
  geom_point() +
  facet_wrap(vars(type, model)) +
  geom_errorbar(aes(xmin = lcl, xmax = ucl)) +
  theme_minimal(base_size = 15) +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=0.5) +
  theme(axis.text = element_text(size = 5)) +
  ylab("") +
  xlab("Z-Score") +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 10))


# all results
ggplot(data = all_results, aes(x = beta, y = exposure, color = type)) + 
  geom_point() +
  facet_wrap(vars(type, model), ncol = 3) +
  geom_errorbar(aes(xmin = lcl, xmax = ucl)) +
  theme_minimal(base_size = 15) +
  theme(axis.text = element_text(size = 5)) +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=0.5) +
  ylab("") +
  xlab("Z-Score") +
  theme(legend.position = "none")



#### SENSITIVITY ANALYSES ####


### Executive function using variable that excludes trails B ###

## Parental education
msm_exec_parent_sens <- lm(execfun_tbex_mc1 ~ parent_educ + age_qx_mc1 +
                           ageusa_cat + lang_exam_mc1,
                           data = cham_msm_pov)

# summarize regression output
summary(msm_exec_parent_sens)

# calculate 95% CIs
beta <- coef(msm_exec_parent_sens)
SE <- coef(summary(msm_exec_parent_sens))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_exec_parent_sens)


## Own education
msm_exec_own_educ_sens <- lm(execfun_tbex_mc1 ~ educcat_mom + parent_educ + 
                             age_qx_mc1 + ageusa_cat +
                             lang_exam_mc1, 
                             data = cham_msm_pov,
                             weights = ipw_own_educ)

# summarize regression output
summary(msm_exec_own_educ_sens)

# calculate 95% CIs
beta <- coef(msm_exec_own_educ_sens)
SE <- coef(summary(msm_exec_own_educ_sens))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_exec_own_educ_sens)


## Poverty
msm_exec_poverty_sens <- lm(execfun_tbex_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                            age_qx_mc1 + ageusa_cat +
                            lang_exam_mc1, 
                            data = cham_msm_pov, 
                            weights = ipw_combined)

# summarize regression output
summary(msm_exec_poverty_sens)

# calculate 95% CIs
beta_sens <- coef(msm_exec_poverty_sens)
SE_sens <- coef(summary(msm_exec_poverty_sens))[,2]
lcl_sens <- beta_sens-qnorm(0.975)*SE_sens
ucl_sens <- beta_sens+qnorm(0.975)*SE_sens
cbind(beta_sens, lcl_sens, ucl_sens)

# calculate number of observations used in model
nobs(msm_exec_poverty_sens)


### Global function using variable that excludes trail B ###

## Parental education
msm_global_parent_sens <- lm(global_tbex_mc1 ~ parent_educ + age_qx_mc1 +
                             ageusa_cat + lang_exam_mc1,
                             data = cham_msm_pov)

# summarize regression output
summary(msm_global_parent_sens)

# calculate 95% CIs
beta <- coef(msm_global_parent_sens)
SE <- coef(summary(msm_global_parent_sens))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in model
nobs(msm_global_parent_sens)


## Own education
msm_global_own_educ_sens <- lm(global_tbex_mc1 ~ educcat_mom + parent_educ + 
                               age_qx_mc1 + ageusa_cat +
                               lang_exam_mc1,
                               data = cham_msm_pov,
                               weights = ipw_own_educ)

# summarize regression output
summary(msm_global_own_educ_sens)

# calculate 95% CIs
beta <- coef(msm_global_own_educ_sens)
SE <- coef(summary(msm_global_own_educ_sens))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_global_own_educ_sens)


## Poverty
msm_global_poverty_sens <- lm(global_tbex_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                              age_qx_mc1 + ageusa_cat +
                              lang_exam_mc1, 
                              data = cham_msm_pov, 
                              weights = ipw_combined)

# summarize regression output
summary(msm_global_poverty_sens)

# calculate 95% CIs
beta <- coef(msm_global_poverty_sens)
SE <- coef(summary(msm_global_poverty_sens))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in model
nobs(msm_global_poverty_sens)


# test whether models are significantly different 
anova(msm_exec_poverty, msm_exec_poverty_sens)
anova(msm_global_poverty, msm_global_poverty_sens)



### Removing 7 individuals who completed the survey over the phone ###
phone_interview_participants <- cham %>% filter(newid %in% c(cham_msm_full$newid)) %>% 
  filter(adm_exam_mc1 == 3) %>% 
  dplyr::select(newid, execfun_rc_mc1, global_rc_mc1)

cham_no_phone <- cham_msm_pov %>% filter(!newid %in% c(phone_interview_participants$newid))


### Memory composite ###

## Parental education
msm_memory_parent_no_phone <- lm(memory_mc1 ~ parent_educ + age_qx_mc1 + ageusa_cat +
                                 lang_exam_mc1,
                                 data = cham_no_phone)

# calculate 95% CIs
beta <- coef(msm_memory_parent_no_phone)
SE <- coef(summary(msm_memory_parent_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

nobs(msm_memory_parent_no_phone)


## Own education
msm_memory_own_educ_no_phone <- lm(memory_mc1 ~ educcat_mom + parent_educ +
                                   age_qx_mc1 + ageusa_cat + lang_exam_mc1,                              
                                   data = cham_no_phone,                              
                                   weights = ipw_own_educ)

# calculate 95% CIs
beta <- coef(msm_memory_own_educ_no_phone)
SE <- coef(summary(msm_memory_own_educ_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

nobs(msm_memory_own_educ_no_phone)


## Poverty
msm_memory_poverty_no_phone <- lm(memory_mc1 ~ pov_cat + educcat_mom + parent_educ +
                                  age_qx_mc1 + ageusa_cat +
                                  lang_exam_mc1,
                                  data = cham_no_phone,
                                  weights = ipw_combined)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_memory_poverty_no_phone)
SE <- coef(summary(msm_memory_poverty_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_memory_poverty_no_phone)



### Executive function ###

## Parental education
msm_exec_parent_no_phone <- lm(execfun_rc_mc1 ~ parent_educ + age_qx_mc1 + ageusa_cat +
                               lang_exam_mc1,
                               data = cham_no_phone)

# calculate 95% CIs
beta <- coef(msm_exec_parent_no_phone)
SE <- coef(summary(msm_exec_parent_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

nobs(msm_exec_parent_no_phone)


## Own education
msm_exec_own_educ_no_phone <- lm(execfun_rc_mc1 ~ educcat_mom + parent_educ +
                                 age_qx_mc1 + ageusa_cat +
                                 lang_exam_mc1,                              
                                 data = cham_no_phone,                                 
                                 weights = ipw_own_educ)

# calculate 95% CIs
beta <- coef(msm_exec_own_educ_no_phone)
SE <- coef(summary(msm_exec_own_educ_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

nobs(msm_exec_own_educ_no_phone)


## Poverty
msm_exec_poverty_no_phone <- lm(execfun_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                                age_qx_mc1 + ageusa_cat +
                                lang_exam_mc1,
                                data = cham_no_phone,
                                weights = ipw_combined)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_exec_poverty_no_phone)
SE <- coef(summary(msm_exec_poverty_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_exec_poverty_no_phone)



### Verbal fluency ###

## Parental education
msm_verbal_parent_no_phone <- lm(verbal_rc_mc1 ~ parent_educ + age_qx_mc1 + ageusa_cat +
                                 lang_exam_mc1,
                                 data = cham_no_phone)

# calculate 95% CIs
beta <- coef(msm_verbal_parent_no_phone)
SE <- coef(summary(msm_verbal_parent_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

nobs(msm_verbal_parent_no_phone)


## Own education
msm_verbal_own_educ_no_phone <- lm(verbal_rc_mc1 ~ educcat_mom + parent_educ +
                                   age_qx_mc1 + ageusa_cat + 
                                   lang_exam_mc1,                              
                                   data = cham_no_phone,                             
                                   weights = ipw_own_educ)

# calculate 95% CIs
beta <- coef(msm_verbal_own_educ_no_phone)
SE <- coef(summary(msm_verbal_own_educ_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

nobs(msm_verbal_own_educ_no_phone)


## Poverty
msm_verbal_poverty_no_phone <- lm(verbal_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                                  age_qx_mc1 + ageusa_cat + 
                                  lang_exam_mc1,
                                  data = cham_no_phone,
                                  weights = ipw_combined)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_verbal_poverty_no_phone)
SE <- coef(summary(msm_verbal_poverty_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_verbal_poverty_no_phone)



### Global function ###

## Parental education
msm_global_parent_no_phone <- lm(global_rc_mc1 ~ parent_educ + age_qx_mc1 + ageusa_cat +
                                 lang_exam_mc1,
                                 data = cham_no_phone)

# calculate 95% CIs
beta <- coef(msm_global_parent_no_phone)
SE <- coef(summary(msm_global_parent_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

nobs(msm_global_parent_no_phone)


## Own education
msm_global_own_educ_no_phone <- lm(global_rc_mc1 ~ educcat_mom + parent_educ +
                                   age_qx_mc1 + ageusa_cat +
                                   lang_exam_mc1,                              
                                   data = cham_no_phone,                               
                                   weights = ipw_own_educ)

# calculate 95% CIs
beta <- coef(msm_global_own_educ_no_phone)
SE <- coef(summary(msm_global_own_educ_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)


nobs(msm_global_own_educ_no_phone)


## Poverty
msm_global_poverty_no_phone <- lm(global_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                                  age_qx_mc1 + ageusa_cat +
                                  lang_exam_mc1,
                                  data = cham_no_phone,
                                  weights = ipw_combined)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_global_poverty_no_phone)
SE <- coef(summary(msm_global_poverty_no_phone))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_global_poverty_no_phone)



### Binary poverty variable ###

## weights for poverty variable re-categorizing missings

# denominator using own education, parent education, years in the US, ACEs, age, 
# married status, working since last visit status, depression (CESD) score,
# anxiety (GAD) score, summary self-rated doctor diagnosed health score, and housing density

# re-categorize those with missing poverty information as part of 
# at or below poverty level
cham_msm_no_missing_pov <- cham_msm_pov %>% filter(!is.na(povcat_mc1))

cham_msm_no_missing_pov %>% group_by(pov_cat) %>% count()

model_denom_poverty_no_miss <- glm(pov_cat ~ own_educ_dich_num + 
                                   parent_educ_dich_num + age_qx_mc1 + ageusa_cat +
                                   aces + married + work_mc1 + 
                                   lang_exam_mc1, 
                                   data = cham_msm_no_missing_pov, 
                                   family = binomial(link = "logit"))

# predict outcomes using model
pred_denom_poverty_no_miss <- predict(model_denom_poverty_no_miss, type = "response")

# calculate IPW for poverty using predictions
cham_msm_no_missing_pov <- cham_msm_no_missing_pov %>% 
  mutate(ipw_poverty_no_miss = ifelse(pov_cat == "at or below poverty level", 
                                      1/(1-pred_denom_poverty_no_miss),
                                      1/pred_denom_poverty_no_miss))


# multiply poverty weight and own education weight
cham_msm_no_missing_pov <- cham_msm_no_missing_pov %>% mutate(ipw_combined_no_miss = ipw_own_educ * ipw_poverty_no_miss)



## trim weights at 99th percentile
trim_weight <- quantile(cham_msm_no_missing_pov$ipw_combined_no_miss, .99)

cham_msm_no_missing_pov <- cham_msm_no_missing_pov %>% 
  mutate(ipw_combined_no_miss = 
           ifelse(ipw_combined_no_miss > trim_weight, trim_weight, ipw_combined_no_miss))


### Memory composite ###
msm_memory_poverty <- lm(memory_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                         age_qx_mc1 + ageusa_cat +
                         lang_exam_mc1, 
                         data = cham_msm_no_missing_pov, 
                         weights = ipw_combined_no_miss)

# summarize regression output
summary(msm_memory_poverty)

# calculate 95% CIs
beta <- coef(msm_memory_poverty)
SE <- coef(summary(msm_memory_poverty))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_memory_poverty)


### Executive function ###
msm_exec_poverty <- lm(execfun_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                       age_qx_mc1 + ageusa_cat +
                       lang_exam_mc1, 
                       data = cham_msm_no_missing_pov, 
                       weights = ipw_combined_no_miss)

# summarize regression output
summary(msm_exec_poverty)

# calculate 95% CIs
beta <- coef(msm_exec_poverty)
SE <- coef(summary(msm_exec_poverty))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_exec_poverty)


### Verbal fluency ###
msm_verbal_poverty <- lm(verbal_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                         age_qx_mc1 + ageusa_cat +
                         lang_exam_mc1, 
                         data = cham_msm_no_missing_pov, 
                         weights = ipw_combined_no_miss)

# summarize regression output
summary(msm_verbal_poverty)

# calculate 95% CIs
beta <- coef(msm_verbal_poverty)
SE <- coef(summary(msm_verbal_poverty))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_verbal_poverty)


### Global function ###
msm_global_poverty <- lm(global_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                         age_qx_mc1 + ageusa_cat +
                         lang_exam_mc1, 
                         data = cham_msm_no_missing_pov, 
                         weights = ipw_combined_no_miss)

# summarize regression output
summary(msm_global_poverty)

# calculate 95% CIs
beta <- coef(msm_global_poverty)
SE <- coef(summary(msm_global_poverty))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_global_poverty)



### Language of assessment and nativity ###


# weights for own education

# use previous exposure (parent education) to calculate numerator
model_num_own_educ <- 1

# denominator using parent education, age arrived to US, ACEs, and age
model_denom_own_educ_lang <- glm(own_educ_dich_num ~ parent_educ_dich_num + 
                                 ageusa_cat + aces + age_qx_mc1 + lang_exam_mc1,
                                 data = cham_msm, 
                                 family = binomial(link = "logit"))

# check multicollinearity - does not show collinearity in model for own education
vif(model_denom_own_educ_lang)


# predict outcomes using models
pred_denom_own_educ_lang <- predict(model_denom_own_educ_lang, type = "response")


# plot residuals against model predictions to see if quadratic shape
res <- residuals(model_denom_own_educ_lang, type = 'deviance')
plot(pred_denom_own_educ_lang,res)


# calculate IPW for own education using predictions
# filter out records that have NAs for aces (only variable that has NAs) to calculate weights properly
cham_msm_own_educ_lang <- cham_msm %>% filter(!is.na(aces_cat)) %>%
  mutate(ipw_own_educ = ifelse(own_educ_dich_num == 0, 
                               1/(1-pred_denom_own_educ_lang), 1/pred_denom_own_educ_lang))



# weights for poverty 
model_num_poverty <- 1

## weights for poverty variable including missings as its own category

# denominator using own education, parent education, age arrived to US, ACEs, age, 
# married status, working since last visit status, depression (CESD) score,
# anxiety (GAD) score, summary self-rated doctor diagnosed health score, and housing density
model_denom_poverty_lang <- multinom(pov_cat ~ own_educ_dich_num + 
                                     parent_educ_dich_num + age_qx_mc1 + ageusa_cat +
                                     lang_exam_mc1 +
                                     aces + married + work_mc1 + 
                                     summary_health_cat_yes + density1_mc1, 
                                     data = cham_msm_own_educ_lang)

# check for multicollinearity - shows collinearity in poverty model
vif(model_denom_poverty_lang)

# predict outcomes using model
pred_denom_poverty_lang <- predict(model_denom_poverty_lang, type = "probs")

# create new dataframe with prediction probabilities
m1 <- pred_denom_poverty_lang[,1]
m2 <- pred_denom_poverty_lang[,2]
m3 <- pred_denom_poverty_lang[,3]

pov_pred_lang <- as.data.frame(cbind(m1, m2, m3))

# calculate inverse probabilities as 1/probability
pov_pred_lang <- pov_pred_lang %>% mutate(p1 = 1/m1, p2 = 1/m2, p3 = 1/m3)


# calculate IPW for poverty using predictions
cham_msm_pov_lang <- cham_msm_own_educ_lang %>% filter(
  !is.na(summary_health_cat_yes) & 
    !is.na(density1_mc1) & !is.na(lang_exam_mc1))

cham_msm_pov_lang$pov_weight <- ifelse(
  cham_msm_pov_lang$pov_cat == "at or below poverty level", pov_pred_lang$p1, 
  ifelse(cham_msm_pov_lang$pov_cat == "above poverty level", pov_pred_lang$p2, 
         ifelse(cham_msm_pov_lang$pov_cat == "Missing", pov_pred_lang$p3, NA)))

# check weights 
cbind(cham_msm_pov_lang, pov_pred_lang) %>% dplyr::select(pov_cat, pov_weight, m1, m2, m3, p1, p2, p3)



# multiply poverty weight and own education weight
cham_msm_pov_lang <- cham_msm_pov_lang %>% mutate(ipw_combined = ipw_own_educ * pov_weight)

## trim weights at 99th percentile
new_cut_off <- quantile(cham_msm_pov_lang$ipw_combined, .99)

cham_msm_pov_lang <- cham_msm_pov_lang %>% 
  mutate(ipw_combined = ifelse(ipw_combined > new_cut_off, new_cut_off, ipw_combined))

ggplot(data = cham_msm_pov_lang, aes(x = ipw_combined)) + 
  geom_histogram()


### Memory composite ###

## Poverty
msm_memory_poverty_lang <- lm(memory_mc1 ~ pov_cat + educcat_mom + parent_educ +
                              age_qx_mc1 + ageusa_cat + lang_exam_mc1,
                              data = cham_msm_pov_lang,
                              weights = ipw_combined)

vif(msm_memory_poverty_lang)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_memory_poverty_lang)
SE <- coef(summary(msm_memory_poverty_lang))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)


### Executive function ###

## Poverty
msm_exec_poverty_lang <- lm(execfun_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                            age_qx_mc1 + ageusa_cat + lang_exam_mc1,
                            data = cham_msm_pov_lang,
                            weights = ipw_combined)

vif(msm_exec_poverty_lang)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_exec_poverty_lang)
SE <- coef(summary(msm_exec_poverty_lang))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)


### Verbal fluency ###

## Poverty
msm_verbal_poverty_lang <- lm(verbal_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                              age_qx_mc1 + ageusa_cat + lang_exam_mc1,
                              data = cham_msm_pov_lang,
                              weights = ipw_combined)

vif(msm_verbal_poverty_lang)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_verbal_poverty_lang)
SE <- coef(summary(msm_verbal_poverty_lang))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)


### Global function ###

## Poverty
msm_global_poverty_lang <- lm(global_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                              age_qx_mc1 + ageusa_cat + lang_exam_mc1,
                              data = cham_msm_pov_lang,
                              weights = ipw_combined)

vif(msm_global_poverty_lang)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_global_poverty_lang)
SE <- coef(summary(msm_global_poverty_lang))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)



### Including langugae (not nativity) in models ###

# weights for own education

# use previous exposure (parent education) to calculate numerator
model_num_own_educ <- 1

# denominator using parent education, age arrived to US, ACEs, and age
model_denom_own_educ_lang <- glm(own_educ_dich_num ~ parent_educ_dich_num + 
                                 aces + age_qx_mc1 + lang_exam_mc1,
                                 data = cham_msm, 
                                 family = binomial(link = "logit"))

# check multicollinearity - does not show collinearity in model for own education
vif(model_denom_own_educ_lang)


# predict outcomes using models
pred_denom_own_educ_lang <- predict(model_denom_own_educ_lang, type = "response")


# plot residuals against model predictions to see if quadratic shape
res <- residuals(model_denom_own_educ_lang, type = 'deviance')
plot(pred_denom_own_educ_lang,res)


# calculate IPW for own education using predictions
# filter out records that have NAs for aces (only variable that has NAs) to calculate weights properly
cham_msm_own_educ_lang <- cham_msm %>% filter(!is.na(aces_cat)) %>%
  mutate(ipw_own_educ = ifelse(own_educ_dich_num == 0, 
                               1/(1-pred_denom_own_educ_lang), 1/pred_denom_own_educ_lang))

# weights for poverty 
model_num_poverty <- 1

## weights for poverty variable including missings as its own category

# denominator using own education, parent education, age arrived to US, ACEs, age, 
# married status, working since last visit status, 
# summary self-rated doctor diagnosed health score, and housing density
model_denom_poverty_lang <- multinom(pov_cat ~ own_educ_dich_num + 
                                     parent_educ_dich_num + age_qx_mc1 +
                                     lang_exam_mc1 +
                                     aces + married + work_mc1 + 
                                     summary_health_cat_yes + density1_mc1, 
                                     data = cham_msm_own_educ_lang)

# check for multicollinearity - shows collinearity in poverty model
vif(model_denom_poverty_lang)

# predict outcomes using model
pred_denom_poverty_lang <- predict(model_denom_poverty_lang, type = "probs")

# create new dataframe with prediction probabilities
m1 <- pred_denom_poverty_lang[,1]
m2 <- pred_denom_poverty_lang[,2]
m3 <- pred_denom_poverty_lang[,3]

pov_pred_lang <- as.data.frame(cbind(m1, m2, m3))

# calculate inverse probabilities as 1/probability
pov_pred_lang <- pov_pred_lang %>% mutate(p1 = 1/m1, p2 = 1/m2, p3 = 1/m3)


# calculate IPW for poverty using predictions
# filter out observations that are missing any covariates (only covariates with missings are summary_health_cat_yes, and density1_mc1)
cham_msm_pov_lang <- cham_msm_own_educ_lang %>% filter(
  !is.na(summary_health_cat_yes) & 
    !is.na(density1_mc1) & !is.na(lang_exam_mc1))

cham_msm_pov_lang$pov_weight <- ifelse(
  cham_msm_pov_lang$pov_cat == "at or below poverty level", pov_pred_lang$p1, 
  ifelse(cham_msm_pov_lang$pov_cat == "above poverty level", pov_pred_lang$p2, 
         ifelse(cham_msm_pov_lang$pov_cat == "Missing", pov_pred_lang$p3, NA)))

# check weights 
cbind(cham_msm_pov_lang, pov_pred_lang) %>% dplyr::select(pov_cat, pov_weight, m1, m2, m3, p1, p2, p3)



# multiply poverty weight and own education weight
cham_msm_pov_lang <- cham_msm_pov_lang %>% mutate(ipw_combined = ipw_own_educ * pov_weight)

## trim weights at 99th percentile
quantile(cham_msm_pov_lang$ipw_combined, .99)

cham_msm_pov_lang <- cham_msm_pov_lang %>% 
  mutate(ipw_combined = ifelse(ipw_combined > 73.0, 73.0, ipw_combined))

ggplot(data = cham_msm_pov_lang, aes(x = ipw_combined)) + 
  geom_histogram()



### Memory composite ###

## Poverty
msm_memory_poverty_lang <- lm(memory_mc1 ~ pov_cat + educcat_mom + parent_educ +
                              age_qx_mc1 + lang_exam_mc1,
                              data = cham_msm_pov_lang,
                              weights = ipw_combined)

vif(msm_memory_poverty_lang)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_memory_poverty_lang)
SE <- coef(summary(msm_memory_poverty_lang))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)



### Executive function ###

## Poverty
msm_exec_poverty_lang <- lm(execfun_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                            age_qx_mc1 + lang_exam_mc1,
                            data = cham_msm_pov_lang,
                            weights = ipw_combined)

vif(msm_exec_poverty_lang)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_exec_poverty_lang)
SE <- coef(summary(msm_exec_poverty_lang))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)


### Verbal fluency ###

## Poverty
msm_verbal_poverty_lang <- lm(verbal_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                              age_qx_mc1 + lang_exam_mc1,
                              data = cham_msm_pov_lang,
                              weights = ipw_combined)

vif(msm_verbal_poverty_lang)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_verbal_poverty_lang)
SE <- coef(summary(msm_verbal_poverty_lang))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)



### Global function ###

## Poverty
msm_global_poverty_lang <- lm(global_rc_mc1 ~ pov_cat + educcat_mom + parent_educ +
                              age_qx_mc1 + lang_exam_mc1,
                              data = cham_msm_pov_lang,
                              weights = ipw_combined)

vif(msm_global_poverty_lang)

# summarize regression output
# calculate 95% CIs
beta <- coef(msm_global_poverty_lang)
SE <- coef(summary(msm_global_poverty_lang))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)



### Correlation between nativity and language ###
cham_msm_pov_lang <- cham_msm_pov_lang %>% mutate(ageusa_cat_num = ifelse(ageusa_cat == 0, 0, 
                                                                          ifelse(ageusa_cat == "<1-17", 1, 
                                                                                 ifelse(ageusa_cat == "18-24", 2, 
                                                                                        3))))
cham_msm_pov_lang <- cham_msm_pov_lang %>% 
  mutate(lang_exam_num = ifelse(lang_exam_mc1 == "English", 1, 0))

model <- lm(lang_exam_num~ageusa_cat_num, data = cham_msm_pov_lang)
cor.test(cham_msm_pov_lang$lang_exam_num, cham_msm_pov_lang$ageusa_cat_num)
summary(model)



#### MEAN DIFFERENCE UNDER DIFFERENT LIFECOURSE EXPOSURE SCENARIOS ####


### MEMORY COMPOSITE ###

# create binary variable for mom's education
# 0 = <=6th grade 
# 1 = >=7th grade
cham_msm_mean <- cham_msm_pov %>% mutate(educcat_bin = ifelse(educcat_mom == "<=6th grade", 0, 1))

# create binary variable for poverty
# 0 = at or below poverty
# 1 = above poverty level
cham_msm_mean <- cham_msm_mean %>% mutate(pov_bin = ifelse(povcat_mc1 == "At or below poverty" |
                                                             is.na(povcat_mc1), 0, 1))


## memory outcome
memory_mean <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   lang_exam_mc1 +
                   aces + married + work_mc1, 
                   data = cham_msm_mean)


### HIGH EDUCATION, HIGH POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$memory_mc1 <- NA

# third copy: education and poverty set to 1, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 1
interv1$pov_bin <- 1
interv1$memory_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_memory <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_memory <- predict(base_memory, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_memory)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_memory)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_memory)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_memory) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_memory)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$memory_mc1 <- NA
  d1 <- d # 3rd copy: education and poverty set to 1, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 1
  d1$pov_bin <- 1
  d1$memory_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 +
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_memory <- bootstrap %>% mutate(outcome = "Memory", type = "High-High")


### HIGH EDUCATION, LOW POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$memory_mc1 <- NA

# third copy: education set to 1 and poverty set to 0, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 1
interv1$pov_bin <- 0
interv1$memory_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_memory <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 +
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_memory <- predict(base_memory, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_memory)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_memory)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_memory)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_memory) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_memory)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$memory_mc1 <- NA
  d1 <- d # 3rd copy: education set to 1 and poverty set to 0, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 1
  d1$pov_bin <- 0
  d1$memory_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_memory2 <- bootstrap %>% mutate(outcome = "Memory", type = "High-Low")
bootstrap_memory <- rbind(bootstrap_memory2, bootstrap_memory)


### LOW EDUCATION, HIGH POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$memory_mc1 <- NA

# third copy: education set to 0 and poverty set to 1, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 0
interv1$pov_bin <- 1
interv1$memory_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_memory <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 +
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_memory <- predict(base_memory, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_memory)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_memory)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_memory)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_memory) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_memory)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$memory_mc1 <- NA
  d1 <- d # 3rd copy: education set to 0 and poverty set to 1, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 0
  d1$pov_bin <- 1
  d1$memory_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1,  
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_memory3 <- bootstrap %>% mutate(outcome = "Memory", type = "Low-High")
bootstrap_memory <- rbind(bootstrap_memory3, bootstrap_memory)


### LOW EDUCATION, LOW POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$memory_mc1 <- NA

# third copy: education set to 0 and poverty set to 0, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 0
interv1$pov_bin <- 0
interv1$memory_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_memory <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 +
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_memory <- predict(base_memory, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_memory)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_memory)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_memory)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_memory) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_memory)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$memory_mc1 <- NA
  d1 <- d # 3rd copy: education and poverty set to 0, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 0
  d1$pov_bin <- 0
  d1$memory_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(memory_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_memory4 <- bootstrap %>% mutate(outcome = "Memory", type = "Low-Low")
bootstrap_memory <- rbind(bootstrap_memory4, bootstrap_memory)


### EXECUTIVE FUNCTION ###

exec_mean <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                 aces + married + work_mc1 +
                 lang_exam_mc1,
                 data = cham_msm_mean)


### HIGH EDUCATION, HIGH POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$execfun_rc_mc1 <- NA

# third copy: education and poverty set to 1, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 1
interv1$pov_bin <- 1
interv1$execfun_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_exec <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                 aces + married + work_mc1 + 
                 lang_exam_mc1, 
                 data = onesample)

onesample$predicted_mean_exec <- predict(base_exec, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_exec)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_exec)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_exec)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_exec) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_exec)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$execfun_rc_mc1 <- NA
  d1 <- d # 3rd copy: education and poverty set to 1, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 1
  d1$pov_bin <- 1
  d1$execfun_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_exec <- bootstrap %>% mutate(outcome = "Executive Function", type = "High-High")


### HIGH EDUCATION, LOW POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$execfun_rc_mc1 <- NA

# third copy: education set to 1 and poverty set to 0, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 1
interv1$pov_bin <- 0
interv1$execfun_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_exec <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                 aces + married + work_mc1 + 
                 lang_exam_mc1, 
                 data = onesample)

onesample$predicted_mean_exec <- predict(base_exec, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_exec)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_exec)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_exec)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_exec) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_exec)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$execfun_rc_mc1 <- NA
  d1 <- d # 3rd copy: education set to 1 and poverty set to 0, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 1
  d1$pov_bin <- 0
  d1$execfun_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 +
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_exec2 <- bootstrap %>% mutate(outcome = "Executive Function", type = "High-Low")
bootstrap_exec <- rbind(bootstrap_exec2, bootstrap_exec)


### LOW EDUCATION, HIGH POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$execfun_rc_mc1 <- NA

# third copy: education set to 0 and poverty set to 1, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 0
interv1$pov_bin <- 1
interv1$execfun_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_exec <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                 aces + married + work_mc1 + 
                 lang_exam_mc1, 
                 data = onesample)

onesample$predicted_mean_exec <- predict(base_exec, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_exec)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_exec)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_exec)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_exec) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_exec)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$execfun_rc_mc1 <- NA
  d1 <- d # 3rd copy: education set to 0 and poverty set to 1, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 0
  d1$pov_bin <- 1
  d1$execfun_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1,
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_exec3 <- bootstrap %>% mutate(outcome = "Executive Function", type = "Low-High")
bootstrap_exec <- rbind(bootstrap_exec3, bootstrap_exec)


### LOW EDUCATION, LOW POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$execfun_rc_mc1 <- NA

# third copy: education and poverty set to 0, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 0
interv1$pov_bin <- 0
interv1$execfun_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_exec <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                 aces + married + work_mc1 + 
                 lang_exam_mc1, 
                 data = onesample)

onesample$predicted_mean_exec <- predict(base_exec, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_exec)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_exec)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_exec)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_exec) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_exec)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$execfun_rc_mc1 <- NA
  d1 <- d # 3rd copy: education and poverty set to 0, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 0
  d1$pov_bin <- 0
  d1$execfun_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(execfun_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_exec4 <- bootstrap %>% mutate(outcome = "Executive Function", type = "Low-Low")
bootstrap_exec <- rbind(bootstrap_exec4, bootstrap_exec)



### VERBAL FLUENCY ###

verbal_mean <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1, 
                   data = cham_msm_mean)


### HIGH EDUCATION, HIGH POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$verbal_rc_mc1 <- NA

# third copy: education and poverty set to 1, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 1
interv1$pov_bin <- 1
interv1$verbal_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_verbal <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1,
                   data = onesample)

onesample$predicted_mean_verbal <- predict(base_verbal, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_verbal)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_verbal)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_verbal)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_verbal) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_verbal)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$verbal_rc_mc1 <- NA
  d1 <- d # 3rd copy: education and poverty set to 1, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 1
  d1$pov_bin <- 1
  d1$verbal_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_verbal <- bootstrap %>% mutate(outcome = "Verbal Fluency", type = "High-High")


### HIGH EDUCATION, LOW POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$verbal_rc_mc1 <- NA

# third copy: education set to 1 and poverty set to 0, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 1
interv1$pov_bin <- 0
interv1$verbal_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_verbal <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1,
                   data = onesample)

onesample$predicted_mean_verbal <- predict(base_verbal, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_verbal)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_verbal)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_verbal)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_verbal) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_verbal)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$verbal_rc_mc1 <- NA
  d1 <- d # 3rd copy: education set to 1 and poverty set to 0, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 1
  d1$pov_bin <- 0
  d1$verbal_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_verbal2 <- bootstrap %>% mutate(outcome = "Verbal Fluency", type = "High-Low")
bootstrap_verbal <- rbind(bootstrap_verbal2, bootstrap_verbal)


### LOW EDUCATION, HIGH POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$verbal_rc_mc1 <- NA

# third copy: education set to 0 and poverty set to 1, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 0
interv1$pov_bin <- 1
interv1$verbal_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_verbal <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_verbal <- predict(base_verbal, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_verbal)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_verbal)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_verbal)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_verbal) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_verbal)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$verbal_rc_mc1 <- NA
  d1 <- d # 3rd copy: education set to 0 and poverty set to 1, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 0
  d1$pov_bin <- 1
  d1$verbal_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_verbal3 <- bootstrap %>% mutate(outcome = "Verbal Fluency", type = "Low-High")
bootstrap_verbal <- rbind(bootstrap_verbal3, bootstrap_verbal)


### LOW EDUCATION, LOW POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$verbal_rc_mc1 <- NA

# third copy: education set to 0 and poverty set to 0, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 0
interv1$pov_bin <- 0
interv1$verbal_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_verbal <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_verbal <- predict(base_verbal, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_verbal)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_verbal)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_verbal)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_verbal) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_verbal)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$verbal_rc_mc1 <- NA
  d1 <- d # 3rd copy: education set to 0 and poverty set to 0, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 0
  d1$pov_bin <- 0
  d1$verbal_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(verbal_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1,  
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_verbal4 <- bootstrap %>% mutate(outcome = "Verbal Fluency", type = "Low-Low")
bootstrap_verbal <- rbind(bootstrap_verbal4, bootstrap_verbal)



#### GLOBAL FUNCTION ####

global_mean <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1, 
                   data = cham_msm_mean)

tidy(global_mean, conf.int = T)

### HIGH EDUCATION, HIGH POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$global_rc_mc1 <- NA

# third copy: education and poverty set to 1, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 1
interv1$pov_bin <- 1
interv1$global_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_global <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 +
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_global <- predict(base_global, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_global)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_global)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_global)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_global) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_global)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$global_rc_mc1 <- NA
  d1 <- d # 3rd copy: education and poverty set to 1, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 1
  d1$pov_bin <- 1
  d1$global_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1,  
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_global <- bootstrap %>% mutate(outcome = "Global Function", type = "High-High")


### HIGH EDUCATION, LOW POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$global_rc_mc1 <- NA

# third copy: education set to 1 and poverty set to 0, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 1
interv1$pov_bin <- 0
interv1$global_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_global <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_global <- predict(base_global, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_global)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_global)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_global)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_global) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_global)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$global_rc_mc1 <- NA
  d1 <- d # 3rd copy: education set to 1 and poverty set to 0, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 1
  d1$pov_bin <- 0
  d1$global_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 +
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_global2 <- bootstrap %>% mutate(outcome = "Global Function", type = "High-Low")
bootstrap_global <- rbind(bootstrap_global2, bootstrap_global)




### LOW EDUCATION, HIGH POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$global_rc_mc1 <- NA

# third copy: education set to 0 and poverty set to 1, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 0
interv1$pov_bin <- 1
interv1$global_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_global <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_global <- predict(base_global, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_global)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_global)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_global)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_global) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_global)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$global_rc_mc1 <- NA
  d1 <- d # 3rd copy: education set to 0 and poverty set to 1, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 0
  d1$pov_bin <- 1
  d1$global_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_global3 <- bootstrap %>% mutate(outcome = "Global Function", type = "Low-High")
bootstrap_global <- rbind(bootstrap_global3, bootstrap_global)


### LOW EDUCATION, LOW POVERTY STATUS vs. LOW EDUCATION, LOW POVERTY STATUS
# make copy of dataset with counterfactual exposure values 
# first copy equal to original one
cham_msm_mean$interv <- -1 

# second copy: education and poverty set to 0, outcome set to missing
interv0 <- cham_msm_mean
interv0$interv <- 0
interv0$educcat_bin <- 0
interv0$pov_bin <- 0
interv0$global_rc_mc1 <- NA

# third copy: education set to 0 and poverty set to 0, outcome set to missing
interv1 <- cham_msm_mean
interv1$interv <- 1
interv1$educcat_bin <- 0
interv1$pov_bin <- 0
interv1$global_rc_mc1 <- NA

# combine datasets
onesample <- rbind(cham_msm_mean, interv0, interv1)

# linear model to estimate mean outcome conditional on exposure and confounders
# parameters are estimated using original observations only (cham_msm_mean)
# parameter estimates are used to predict mean outcome for observations with 
## treatment set to 0 (interv=0) and to 1 (interv=1)

base_global <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 +
                   lang_exam_mc1, 
                   data = onesample)

onesample$predicted_mean_global <- predict(base_global, onesample)


## estimate mean outcome in each of the groups interv=0, and interv=1
## this mean outcome is a weighted average of the mean outcomes in each combination 
## of values of treatment and confounders, that is, the standardized outcome
mean(onesample[which(onesample$interv == -1), ]$predicted_mean_global)

mean(onesample[which(onesample$interv == 0), ]$predicted_mean_global)

mean(onesample[which(onesample$interv == 1), ]$predicted_mean_global)


mean(onesample[which(onesample$interv == 1), ]$predicted_mean_global) - mean(onesample[which(onesample$interv == 0), ]$predicted_mean_global)


## function to calculate difference in means
standardization <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices, ] # 1st copy: equal to original one
  d$interv <- -1
  d0 <- d # 2nd copy: education and poverty set to 0, outcome set to missing
  d0$interv <- 0
  d0$educcat_bin <- 0
  d0$pov_bin <- 0
  d0$global_rc_mc1 <- NA
  d1 <- d # 3rd copy: education set to 0 and poverty set to 0, outcome set to missing
  d1$interv <- 1
  d1$educcat_bin <- 0
  d1$pov_bin <- 0
  d1$global_rc_mc1 <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set
  # treatment (interv=0 and interv=1)
  fit <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
             aces + married + work_mc1 + 
             lang_exam_mc1, 
             data = d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
           mean(d.onesample$predicted_meanY[d.onesample$interv == 0]))
}

## bootstrap
results <- boot(data = cham_msm_mean, statistic = standardization, R = 1000)

## generating confidence intervals
se <- sd(results$t[, 1])
meant0 <- results$t0
ll <- meant0 - qnorm(0.975) * se
ul <- meant0 + qnorm(0.975) * se

bootstrap <- data.frame(
  " " = "Exposure - No Exposure",
  estimate = meant0,
  std.error = se,
  conf.low = ll,
  conf.high = ul,
  check.names = FALSE)
bootstrap

bootstrap_global4 <- bootstrap %>% mutate(outcome = "Global Function", type = "Low-Low")
bootstrap_global <- rbind(bootstrap_global4, bootstrap_global)



### each year of additional age corresponds to 0.018 lower cognitive score (see output for age variable in regression below)
### compared to those who experience relatively low SES, those who experience relatively high SES had estimated cognitive scores equivalent to 
### a difference of about 36 years 
### mean cog score output for relatively high SES = 0.656 from bootstrap_global table
### 0.656/0.018 = 36.4 years in age

global_mean <- glm(global_rc_mc1 ~ pov_bin + educcat_bin + age_qx_mc1 + ageusa_cat +
                   aces + married + work_mc1 + 
                   lang_exam_mc1,
                   data = cham_msm_mean)

tidy(global_mean, conf.int = T)



### Combine bootstrapped results ###

bootstrap_results <- rbind(bootstrap_memory, bootstrap_exec, bootstrap_verbal, bootstrap_global)
bootstrap_results

bootstrap_results <- bootstrap_results %>% mutate(type_long = 
                                                    ifelse(type == "Low-Low", 
                                                           "Relatively low lifecourse SES", 
                                                           ifelse(type == "High-Low", 
                                                                  "Change from relatively high to low SES",
                                                                  ifelse(type == "Low-High",
                                                                         "Change from relatively low to high SES", 
                                                                         "Relatively high lifecourse SES"))))



# reorder exposure variable for plots
bootstrap_results <- bootstrap_results %>% mutate(type_long = 
                                                    fct_relevel(type_long,
                                                                "Relatively high lifecourse SES",
                                                                "Change from relatively high to low SES",
                                                                "Change from relatively low to high SES",
                                                                "Relatively low lifecourse SES"
                                                    ))

# plot results
pdf(file = "../../Desktop/mean_diff.pdf", 
    width = 8, 
    height = 4)

ggplot(data = bootstrap_results, aes(x = estimate, y = type_long, color = outcome)) + 
  geom_point() +
  facet_wrap(~factor(outcome, levels = c("Memory", "Executive Function", 
                                         "Verbal Fluency", "Global Function"))) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high)) +
  theme_minimal(base_size = 15) +
  theme(axis.text = element_text(size = 10)) +
  ylab("") +
  xlab("Mean difference") + 
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10), 
        strip.text.x = element_text(size = 10)) + 
  scale_color_brewer(palette = "Set2", name = "outcome")

dev.off()



#### BASELINE SAMPLE CHARACTERISTICS ####


# select variables to identify baseline sample
cham_base <- cham %>% dplyr::select(newid, cham, ageusa_cat, ageusa_bl, age_bl, 
                                    parent_educ, educcat_mom,
                                    lang_exam_mc1,
                                    marstat_bl, worksp_bl, workc_bl, povcat_bl,
                                    ipovcat_bl, ipovcat_bl_met, hbp_bl, hbpage_bl,
                                    diab_bl, diabage_bl, cancer_bl, 
                                    cancerage_bl, qx_matcog, qx_9y, date_9y, 
                                    child_under18_9y, preg_9y, yrsusa_9y, age_qx_mc1, 
                                    date_m18y, married_m16y, lvhome_n_mc1, country_mom,
                                    racecat_mom, momdl_age, momdl_age2, age_qx_mc1, 
                                    date_9y, date_16y, date_m18y, age_exam_mc1, age_bld_mc1, 
                                    raceeth_mom, marcat_pg)
cham_base %>% dplyr::select(ageusa_bl, age_bl, momdl_age, momdl_age2)

# categorize baseline married status
cham_base <- cham_base %>% mutate(marstat_bl = 
                                    ifelse(marstat_bl == 1 | marstat_bl == 2, "Married", 
                                           ifelse(marstat_bl == 3, "Separated", 
                                                  ifelse(marstat_bl == 4, "Divorced", 
                                                         ifelse(marstat_bl == 5, "Widowed", 
                                                                ifelse(marstat_bl == 6, "Single", marstat_bl))))))

# categorize baseline race/ethnicity
cham_base <- cham_base %>% mutate(raceeth_mom = ifelse(raceeth_mom == 1, "Mexican", 
                                                       ifelse(raceeth_mom == 2, "Mexican Indian",
                                                              ifelse(raceeth_mom == 3, "Mexican-American/Chicana", 
                                                                     ifelse(raceeth_mom == 4, "Other Latina", 
                                                                            ifelse(raceeth_mom == 5, "Asian or Pacific Islander", 
                                                                                   ifelse(raceeth_mom == 6, "White non-Latina", 
                                                                                          ifelse(raceeth_mom == 7, "Black non-Latina",
                                                                                                 ifelse(raceeth_mom == 8, "Other", NA)))))))))

# categorize baseline work during index pregnancy (retrospective for CHAM2)
cham_base <- cham_base %>% mutate(worksp_bl =
                                    ifelse(worksp_bl == -9, "Not working", 
                                           ifelse(worksp_bl == 1, "Fieldwork", 
                                                  ifelse(worksp_bl == 2, "Agwork", 
                                                         ifelse(worksp_bl == 3, "Other work", worksp_bl)))))

# categorize work at baseline (9y for CHAM2)
cham_base <- cham_base %>% mutate(workc_bl =
                                    ifelse(workc_bl == -9, "Not working", 
                                           ifelse(workc_bl == 1, "Fieldwork", 
                                                  ifelse(workc_bl == 2, "Agwork", 
                                                         ifelse(workc_bl == 3, "Other work", workc_bl)))))

# categorize baseline poverty 
cham_base <- cham_base %>% mutate(povcat_bl = 
                                    ifelse(povcat_bl == 1, "At or below poverty", 
                                           ifelse(povcat_bl == 2, "Poverty-200%",
                                                  ifelse(povcat_bl == 3, ">=200% Poverty", povcat_bl))))

# categorize imputed baseline poverty 
cham_base <- cham_base %>% mutate(ipovcat_bl = 
                                    ifelse(ipovcat_bl == 1, "At or below poverty", 
                                           ifelse(ipovcat_bl == 2, "Poverty-200%",
                                                  ifelse(ipovcat_bl == 3, ">=200% Poverty", povcat_bl))))

# categorize baseline high blood pressure
cham_base <- cham_base %>% mutate(hbp_bl = 
                                    ifelse(hbp_bl == 0, "No", 
                                           ifelse(hbp_bl == 1, "Yes", 
                                                  ifelse(hbp_bl == 9, "Don't know", hbp_bl))))

# categorize baseline diabetes
cham_base <- cham_base %>% mutate(diab_bl = 
                                    ifelse(diab_bl == 0, "No", 
                                           ifelse(diab_bl == 1, "Yes", 
                                                  ifelse(diab_bl == 9, "Don't know", diab_bl))))

# categorize baseline cancer
cham_base <- cham_base %>% mutate(cancer_bl = 
                                    ifelse(cancer_bl == 0, "No", 
                                           ifelse(cancer_bl == 1, "Yes", 
                                                  ifelse(cancer_bl == 9, "Don't know", cancer_bl))))


# set missing health variables as "no" 
cham_base <- cham_base %>% mutate(hbp_bl = ifelse(is.na(hbp_bl), "No", hbp_bl))

cham_base <- cham_base %>% mutate(diab_bl = ifelse(is.na(diab_bl), "No", diab_bl))

cham_base <- cham_base %>% mutate(cancer_bl = ifelse(is.na(cancer_bl), "No", cancer_bl))


# create collapsed chronic health conditions variable
cham_base <- cham_base %>% mutate(health = ifelse(hbp_bl == "No" & diab_bl == "No" &
                                                    cancer_bl == "No", "None", 
                                                  ifelse(hbp_bl == "Yes" | diab_bl == "Yes" |
                                                           cancer_bl == "Yes", "1+", NA)))

cham_base %>% group_by(health) %>% count()

cham_base %>% filter(hbp_bl == "Don't know" | diab_bl == "Don't know" | cancer_bl == "Don't know")

# recategorize those missing collapsed health score into "None" category
cham_base <- cham_base %>% mutate(health = ifelse(is.na(health), "None", health))


# check number of participants who participated in MATCOG sample
cham_base %>% group_by(qx_matcog) %>% count() 

# check number of participants who participated in CHAM1 and CHAM2
cham_base %>% group_by(cham) %>% count()

# remove individual with 314 who is missing almost all covariates
cham_base <- cham_base %>% filter(!newid == 314)

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


### 1999 cohort ###

# create table of baseline characteristics for 1999 cohort

# MatCog sample - put values of 0 for HBP, diabetes, cancer if missing/don't know
# none of the matcog participants missing HBP 
cham_1999_cohort %>% filter(qx_matcog == 1) %>% group_by(hbp_bl) %>% count()

# 5 matcog baseline people missing diabetes, 1 doesn't know
cham_1999_cohort %>% filter(qx_matcog == 1) %>% group_by(diab_bl) %>% count()
cham_1999_cohort <- cham_1999_cohort %>% 
  mutate(diab_bl = ifelse(qx_matcog == 1 & (is.na(diab_bl) | diab_bl == "Don't know"), "No", diab_bl))

# No matcog baseline people missing cancer
cham_1999_cohort %>% filter(qx_matcog == 1) %>% group_by(cancer_bl) %>% count()


# select variables for baseline table 1
cham_1999_cohort_tab1 <- cham_1999_cohort %>%
  dplyr::select(newid, cham, ageusa_cat, age_bl, parent_educ, educcat_mom, 
                marstat_bl, country_mom, worksp_bl, worksp_bl,
                povcat_bl, ipovcat_bl, ipovcat_bl_met, hbp_bl,
                hbpage_bl, diab_bl, diabage_bl, cancer_bl,
                cancerage_bl, racecat_mom, momdl_age, 
                momdl_age2, age_qx_mc1, date_9y, date_16y, 
                date_m18y, age_exam_mc1, age_bld_mc1, raceeth_mom, 
                health, marcat_pg)

# calculate age in 2023 (start of MatCog)
cham_1999_cohort_tab1 <- cham_1999_cohort_tab1 %>% mutate(mom_age_today = momdl_age + 23.5)

vars <- c("ageusa_cat", "mom_age_today", "raceeth_mom",
          "parent_educ", "educcat_mom", "country_mom",
          "marstat_bl", "worksp_bl", "worksp_bl", "povcat_bl", 
          "ipovcat_bl", "ipovcat_bl_met", "health")

# update marital status for person missing marstat_bl with marcat_pg status
cham_1999_cohort_tab1 <- cham_1999_cohort_tab1 %>% 
  mutate(marstat_bl = ifelse(is.na(marstat_bl) & marcat_pg == 1, 
                             "Married", marstat_bl))
cham_1999_cohort_tab1 %>% group_by(marstat_bl) %>% count()

# make age came to the US a factor variable and set 0 as reference
cham_1999_cohort_tab1$ageusa_cat <- factor(cham_1999_cohort_tab1$ageusa_cat)
cham_1999_cohort_tab1$ageusa_cat <- relevel(cham_1999_cohort_tab1$ageusa_cat, ref = "0")

# make married variable a factor variable and set not married as reference
cham_1999_cohort_tab1$marstat_bl <- factor(cham_1999_cohort_tab1$marstat_bl)
cham_1999_cohort_tab1$marstat_bl <- relevel(cham_1999_cohort_tab1$marstat_bl, ref = "Married")

# make nativity variable a factor variable and set USA as reference
cham_1999_cohort_tab1$country_mom <- factor(cham_1999_cohort_tab1$country_mom)
cham_1999_cohort_tab1$country_mom <- relevel(cham_1999_cohort_tab1$country_mom, ref = 1)

# make race variable a factor variable and set 1 as reference
cham_1999_cohort_tab1$racecat_mom <- factor(cham_1999_cohort_tab1$racecat_mom)
cham_1999_cohort_tab1$racecat_mom <- relevel(cham_1999_cohort_tab1$racecat_mom, ref = 1)

# make work variable a factor variable and set 1 as reference
cham_1999_cohort_tab1$worksp_bl <- factor(cham_1999_cohort_tab1$worksp_bl)
cham_1999_cohort_tab1$worksp_bl <- relevel(cham_1999_cohort_tab1$worksp_bl, ref = "Agwork")

# create table 1
CreateTableOne(vars = vars, data = cham_1999_cohort_tab1)

# create table of dataset not missing any covariates
cham_1999_no_miss_cov <- cham_1999_cohort_tab1 %>% 
  filter(!(is.na(mom_age_today) | is.na(ageusa_cat) |
             is.na(marstat_bl) | is.na(worksp_bl))) 

cham_1999_no_miss_cov %>% summarize(mean = mean(mom_age_today), 
                                    sd = sd(mom_age_today))

cham_1999_no_miss_cov %>% group_by(ageusa_cat) %>% count() 

cham_1999_no_miss_cov %>% group_by(marstat_bl) %>% count()

cham_1999_no_miss_cov %>% group_by(worksp_bl) %>% count()

cham_1999_no_miss_cov %>% group_by(diab_bl) %>% count()

cham_1999_no_miss_cov %>% group_by(hbp_bl) %>% count()

cham_1999_no_miss_cov %>% group_by(parent_educ) %>% count()

cham_1999_no_miss_cov %>% group_by(educcat_mom) %>% count()

cham_1999_no_miss_cov %>% group_by(ipovcat_bl) %>% count()


### 2009 cohort ###

# create table of baseline characteristics for 2009 cohort

# MatCog sample - put values of 0 for HBP, diabetes, cancer if missing/don't know
# 2 matcog missing HBP, 2 don't know
cham_2009_cohort %>% filter(qx_matcog == 1) %>% group_by(hbp_bl) %>% count()
cham_2009_cohort <- cham_2009_cohort %>% 
  mutate(diab_bl = ifelse(qx_matcog == 1 & (is.na(diab_bl) | diab_bl == "Don't know"), "No", diab_bl))

# None of matcog baseline people missing diabetes
cham_2009_cohort %>% filter(qx_matcog == 1) %>% group_by(diab_bl) %>% count()

# 2 matcog baseline people missing cancer
cham_2009_cohort %>% filter(qx_matcog == 1) %>% group_by(cancer_bl) %>% count()
cham_2009_cohort <- cham_2009_cohort %>% 
  mutate(cancer_bl = ifelse(qx_matcog == 1 & (is.na(cancer_bl) | cancer_bl == "Don't know"), "No", cancer_bl))

# select variables for baseline table 1
cham_2009_cohort_tab1 <- cham_2009_cohort %>%
  dplyr::select(newid, cham, ageusa_cat, age_bl,
                lang_exam_mc1,
                parent_educ, educcat_mom, marstat_bl,
                country_mom, worksp_bl, worksp_bl, 
                povcat_bl, ipovcat_bl, ipovcat_bl_met, hbp_bl,
                hbpage_bl,diab_bl, diabage_bl, cancer_bl,
                cancerage_bl, racecat_mom, momdl_age, 
                momdl_age2, age_qx_mc1, date_9y, date_16y, 
                date_m18y, age_exam_mc1, age_bld_mc1, raceeth_mom,
                health, marcat_pg)

# calculate age in 2023 (start of MatCog)
cham_2009_cohort_tab1 <- cham_2009_cohort_tab1 %>% 
  mutate(mom_age_today = ifelse(cham == 1, momdl_age + 23.5, momdl_age + 13.5))

vars <- c("ageusa_cat", "mom_age_today", "raceeth_mom",
          "parent_educ", "educcat_mom", "lang_exam_mc1", 
          "country_mom",
          "marstat_bl", "worksp_bl", "povcat_bl", 
          "ipovcat_bl", "ipovcat_bl_met", "health")

# update marital status for person missing marstat_bl with marcat_pg status
cham_2009_cohort_tab1 <- cham_2009_cohort_tab1 %>% 
  mutate(marstat_bl = ifelse(is.na(marstat_bl) & marcat_pg == 1, 
                             "Married", marstat_bl))
cham_2009_cohort_tab1 %>% group_by(marstat_bl) %>% count()

# make age came to the US a factor variable and set 0 as reference
cham_2009_cohort_tab1$ageusa_cat <- factor(cham_2009_cohort_tab1$ageusa_cat)
cham_2009_cohort_tab1$ageusa_cat <- relevel(cham_2009_cohort_tab1$ageusa_cat, ref = "0")

# make married variable a factor variable and set not married as reference
cham_2009_cohort_tab1$marstat_bl <- factor(cham_2009_cohort_tab1$marstat_bl)
cham_2009_cohort_tab1$marstat_bl <- relevel(cham_2009_cohort_tab1$marstat_bl, ref = "Married")

# make nativity variable a factor variable and set USA as reference
cham_2009_cohort_tab1$country_mom <- factor(cham_2009_cohort_tab1$country_mom)
cham_2009_cohort_tab1$country_mom <- relevel(cham_2009_cohort_tab1$country_mom, ref = 1)

# create table 1
CreateTableOne(vars = vars, data = cham_2009_cohort_tab1)

# create table of dataset not missing any covariates
cham_2009_no_miss_cov <- cham_2009_cohort_tab1 %>% 
  filter(!(is.na(mom_age_today) | is.na(ageusa_cat) |
             is.na(marstat_bl) | is.na(worksp_bl))) 

cham_2009_no_miss_cov %>% summarize(mean = mean(mom_age_today), 
                                    sd = sd(mom_age_today))

cham_2009_no_miss_cov %>% group_by(ageusa_cat) %>% count() 

cham_2009_no_miss_cov %>% group_by(marstat_bl) %>% count()

cham_2009_no_miss_cov %>% group_by(worksp_bl) %>% count()

cham_2009_no_miss_cov %>% group_by(diab_bl) %>% count()

cham_2009_no_miss_cov %>% group_by(hbp_bl) %>% count()

cham_2009_no_miss_cov %>% group_by(health) %>% count()

cham_2009_no_miss_cov %>% group_by(parent_educ) %>% count()

cham_2009_no_miss_cov %>% group_by(educcat_mom) %>% count()

cham_2009_no_miss_cov %>% group_by(ipovcat_bl) %>% count()

cham_2009_no_miss_cov %>% group_by(cancer_bl) %>% count()



#### IPAW FOR 2009 COHORT ####


# weights for poverty 

# isolate the combined weights from the 515 women who were included from the matcog sample 
pov_weights_mat_cog <- cham_msm_pov %>% dplyr::select(newid, pov_weight, ipw_combined)

# create variable to indicate whether woman from 2009 sample participated in matcog
cham_2009_cohort_tab1 <- cham_2009_cohort_tab1 %>% mutate(participate_mat_cog = 
                                                            ifelse(newid %in% c(pov_weights_mat_cog$newid), 1, 0))

cham_2009_cohort_tab1 %>% group_by(participate_mat_cog) %>% count()

# add column of combined weights from matcog sample to 2009 baseline sample
# if person is not in matcog sample their weight should be 0
cham_2009_cohort_tab1 <- merge(cham_2009_cohort_tab1, 
                               pov_weights_mat_cog, by = "newid", all = TRUE) 
cham_2009_cohort_tab1 <- cham_2009_cohort_tab1 %>% 
  mutate(ipw_combined = ifelse(is.na(ipw_combined), 0, ipw_combined))

cham_2009_cohort_tab1 %>% filter(ipw_combined == 0)
# denominator using current age, age arrived to US, own education, 
# marital status, working at baseline/pregnancy status, poverty status, and composite health score

# 3 individuals missing work status (worksp_bl)
cham_2009_no_miss_cov <- cham_2009_cohort_tab1 %>% 
  filter(!(is.na(mom_age_today) | is.na(ageusa_cat) |
             is.na(educcat_mom) | is.na(marstat_bl) | 
             is.na(worksp_bl) | is.na(ipovcat_bl) |
             is.na(diab_bl) | is.na(hbp_bl))) 


model_denom_poverty_2009 <- glm(participate_mat_cog ~ mom_age_today +
                                ageusa_cat + educcat_mom + marstat_bl +
                                lang_exam_mc1 +
                                worksp_bl + ipovcat_bl + diab_bl + 
                                hbp_bl,
                                data = cham_2009_no_miss_cov, 
                                family = binomial(link = "logit"))

# predict outcomes using model
pred_denom_poverty_2009 <- predict(model_denom_poverty_2009, type = "response")

# calculate probability of attrition weights using predictions
cham_2009_no_miss_cov <- cham_2009_no_miss_cov %>%
  mutate(ipaw = ifelse(participate_mat_cog == 0, 
                       1/(1-pred_denom_poverty_2009), 1/pred_denom_poverty_2009))

# calculate full superweight using IPAW and combined IPWs from original model
cham_2009_no_miss_cov <- cham_2009_no_miss_cov %>% mutate(ipw_super = ipw_combined*ipaw)

# note that the 3 individuals missing baseline info are also in matcog - what to do about their weights? leave as combined (not including ipaw weights) or remove from sample?
matcog_ipaw <- cham_2009_no_miss_cov %>% filter(participate_mat_cog == 1)


# include 3 women missing baseline data - use original IPW weights 
remaining_3_missing_baseline <- cham_2009_cohort_tab1 %>% 
  filter(is.na(worksp_bl) | is.na(marstat_bl))

remaining_3_missing_baseline <- remaining_3_missing_baseline %>% mutate(ipaw = ipw_combined, 
                                                                        ipw_super = ipw_combined)

matcog_full_ipaw_ipw <- rbind(matcog_ipaw, remaining_3_missing_baseline)


# select only relevant columns for easier merging
matcog_ipaw <- matcog_ipaw %>% dplyr::select(newid, ipw_combined, ipaw, ipw_super)
matcog_full_ipaw_ipw <- matcog_full_ipaw_ipw %>% dplyr::select(newid, ipw_combined, ipaw, ipw_super)



### Memory composite ###

### Only including the 512 women who have true IPAW weights

# use super weights
# filter the IDs in the cham_msm_pov dataset to only include those who are in the 2009 dataset
cham_msm_pov_filtered <- cham_msm_pov %>% filter(newid %in% c(matcog_ipaw$newid))

# add IPAW weights to cham_msm_pov dataset
cham_msm_pov_filtered <- merge(cham_msm_pov_filtered, matcog_ipaw, by = "newid")

# run linear regression
msm_memory_poverty_ipaw <- lm(memory_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                              age_qx_mc1 + ageusa_cat +
                              lang_exam_mc1, 
                              data = cham_msm_pov_filtered, 
                              weights = ipw_super)

# summarize regression output
summary(msm_memory_poverty_ipaw)

# calculate 95% CIs
beta <- coef(msm_memory_poverty_ipaw)
SE <- coef(summary(msm_memory_poverty_ipaw))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_memory_poverty_ipaw)



### Including the 515 women - 3 of whom only have IPW weights instead of IPAW weights
# use super weights
# filter the IDs in the cham_msm_pov dataset to only include those who are in the 2009 dataset
cham_msm_pov_filtered_total <- cham_msm_pov %>% filter(newid %in% c(matcog_full_ipaw_ipw$newid))

# add IPAW weights to cham_msm_pov dataset
cham_msm_pov_filtered_total <- merge(cham_msm_pov_filtered_total, matcog_full_ipaw_ipw, by = "newid")

# run linear regression
msm_memory_poverty_ipaw_full <- lm(memory_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                                   age_qx_mc1 + ageusa_cat +
                                   lang_exam_mc1, 
                                   data = cham_msm_pov_filtered_total, 
                                   weights = ipw_super)

# summarize regression output
summary(msm_memory_poverty_ipaw_full)

# calculate 95% CIs
beta_full <- coef(msm_memory_poverty_ipaw_full)
SE <- coef(summary(msm_memory_poverty_ipaw_full))[,2]
lcl_full <- beta_full-qnorm(0.975)*SE 
ucl_full <- beta_full+qnorm(0.975)*SE
cbind(beta_full, lcl_full, ucl_full)

# calculate number of observations used in linear regression model
nobs(msm_memory_poverty_ipaw_full)


### Executive function ###

### Only including the women who have true IPAW weights

# run linear regression
msm_exec_poverty_ipaw <- lm(execfun_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                            age_qx_mc1 + ageusa_cat + 
                            lang_exam_mc1, 
                            data = cham_msm_pov_filtered, 
                            weights = ipw_super)

# summarize regression output
summary(msm_exec_poverty_ipaw)

# calculate 95% CIs
beta <- coef(msm_exec_poverty_ipaw)
SE <- coef(summary(msm_exec_poverty_ipaw))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_exec_poverty_ipaw)



### Including the 3 women who only have IPW weights instead of IPAW weights

# run linear regression
msm_exec_poverty_ipaw_full <- lm(execfun_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                                 age_qx_mc1 + ageusa_cat +
                                 lang_exam_mc1, 
                                 data = cham_msm_pov_filtered_total, 
                                 weights = ipw_super)

# summarize regression output
summary(msm_exec_poverty_ipaw_full)

# calculate 95% CIs
beta_full <- coef(msm_exec_poverty_ipaw_full)
SE <- coef(summary(msm_exec_poverty_ipaw_full))[,2]
lcl_full <- beta_full-qnorm(0.975)*SE 
ucl_full <- beta_full+qnorm(0.975)*SE
cbind(beta_full, lcl_full, ucl_full)

# calculate number of observations used in linear regression model
nobs(msm_exec_poverty_ipaw_full)


### Verbal fluency ###

### Only including the women who have true IPAW weights

# run linear regression
msm_verbal_poverty_ipaw <- lm(verbal_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                              age_qx_mc1 + ageusa_cat +
                              lang_exam_mc1, 
                              data = cham_msm_pov_filtered, 
                              weights = ipw_super)

# summarize regression output
summary(msm_verbal_poverty_ipaw)

# calculate 95% CIs
beta <- coef(msm_verbal_poverty_ipaw)
SE <- coef(summary(msm_verbal_poverty_ipaw))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_verbal_poverty_ipaw)



### Including the 3 women who only have IPW weights instead of IPAW weights

# run linear regression
msm_verbal_poverty_ipaw_full <- lm(verbal_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                                   age_qx_mc1 + ageusa_cat + 
                                   lang_exam_mc1, 
                                   data = cham_msm_pov_filtered_total, 
                                   weights = ipw_super)

# summarize regression output
summary(msm_verbal_poverty_ipaw_full)

# calculate 95% CIs
beta_full <- coef(msm_verbal_poverty_ipaw_full)
SE <- coef(summary(msm_verbal_poverty_ipaw_full))[,2]
lcl_full <- beta_full-qnorm(0.975)*SE 
ucl_full <- beta_full+qnorm(0.975)*SE
cbind(beta_full, lcl_full, ucl_full)

# calculate number of observations used in linear regression model
nobs(msm_verbal_poverty_ipaw_full)


### Global function ###

### Only including the women who have true IPAW weights

# run linear regression
msm_global_poverty_ipaw <- lm(global_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                              age_qx_mc1 + ageusa_cat +
                              lang_exam_mc1, 
                              data = cham_msm_pov_filtered, 
                              weights = ipw_super)

# summarize regression output
summary(msm_global_poverty_ipaw)

# calculate 95% CIs
beta <- coef(msm_global_poverty_ipaw)
SE <- coef(summary(msm_global_poverty_ipaw))[,2]
lcl <- beta-qnorm(0.975)*SE 
ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)

# calculate number of observations used in linear regression model
nobs(msm_global_poverty_ipaw)



### Including the 3 women who only have IPW weights instead of IPAW weights

# run linear regression
msm_global_poverty_ipaw_full <- lm(global_rc_mc1 ~ pov_cat + educcat_mom + parent_educ + 
                                   age_qx_mc1 + ageusa_cat +
                                   lang_exam_mc1, 
                                   data = cham_msm_pov_filtered_total, 
                                   weights = ipw_super)

# summarize regression output
summary(msm_global_poverty_ipaw_full)

# calculate 95% CIs
beta_full <- coef(msm_global_poverty_ipaw_full)
SE <- coef(summary(msm_global_poverty_ipaw_full))[,2]
lcl_full <- beta_full-qnorm(0.975)*SE 
ucl_full <- beta_full+qnorm(0.975)*SE
cbind(beta_full, lcl_full, ucl_full)

# calculate number of observations used in linear regression model
nobs(msm_global_poverty_ipaw_full)



#### BASELINE CHARACTERISTICS FOR MATCOG, 1999 COHORT, 2009 COHORT ####


# calculate baseline characteristics for those in matcog sample
matcog_baseline <- cham_2009_cohort_tab1 %>% filter(newid %in% c(cham_desc$newid))

matcog_baseline %>% group_by(educcat_mom) %>% count()
matcog_baseline %>% group_by(ageusa_cat) %>% count()
matcog_baseline %>% group_by(raceeth_mom) %>% count()
matcog_baseline %>% group_by(ipovcat_bl) %>% count()
matcog_baseline %>% group_by(marstat_bl) %>% count()
matcog_baseline %>% group_by(worksp_bl) %>% count()
matcog_baseline %>% group_by(diab_bl) %>% count()
matcog_baseline %>% group_by(hbp_bl) %>% count()


# did not participate in matcog 
no_matcog_cham1 <- cham_1999_cohort_tab1 %>% filter(!newid %in% c(cham_desc$newid))
no_matcog_cham2 <- cham_2009_cohort_tab1 %>% filter(!newid %in% c(cham_desc$newid))
no_matcog_cham2 <- no_matcog_cham2 %>% 
  dplyr::select(-participate_mat_cog, -pov_weight, -ipw_combined)

no_matcog <- rbind(no_matcog_cham1, no_matcog_cham2)
no_matcog <- distinct(no_matcog, newid, .keep_all = T)

no_matcog %>% group_by(educcat_mom) %>% count()
no_matcog %>% group_by(ageusa_cat) %>% count()
no_matcog %>% group_by(raceeth_mom) %>% count()
no_matcog %>% group_by(ipovcat_bl) %>% count()
no_matcog %>% group_by(marstat_bl) %>% count()
no_matcog %>% group_by(worksp_bl) %>% count()
no_matcog %>% group_by(health) %>% count()
no_matcog %>% group_by(diab_bl) %>% count()
no_matcog %>% group_by(hbp_bl) %>% count()


cham_1999_cohort %>% group_by(diab_bl) %>% count()
cham_1999_cohort %>% group_by(hbp_bl) %>% count()

cham_2009_cohort %>% group_by(diab_bl) %>% count()
cham_2009_cohort %>% group_by(hbp_bl) %>% count()



#### GENERATE MAPS OF STUDY SAMPLE ####


# generate county data
CA_county <- counties(state = "CA")

# create colors for two counties of interest
CA_county <- CA_county %>% mutate(fill = 
                                    ifelse(COUNTYFP == "087", "Santa Cruz", 
                                           ifelse(COUNTYFP == "053", "Monterey", "Other CA county")))


x11(type = "cairo")
ggplot(CA_county, aes(fill = fill)) + geom_sf() +
  scale_fill_manual(values = c(`Santa Cruz` = "red", `Monterey` = "blue", 
                               `Other CA county` = "grey")) +
  guides(fill = guide_legend(title = "Region")) + 
  theme_minimal(base_size = 10) + 
  annotation_scale() +
  annotation_north_arrow(which_north = "true", 
                         pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"), location = "tr")

# generate census tract data for Santa Cruz and Monterey
CA_tract <- tracts(state = "CA", county = c("Santa Cruz", "Monterey"))

CA_tract <- CA_tract %>% 
  mutate(fill = ifelse(NAME %in% 
                         c("1.01", "1.03", "1.04", "1.05", "1.06", "2", 
                           "3","4", "5.01", "5.02", "6", "7.01", "7.02", 
                           "8", "9", "12", "13", "14", "15", "16", "17",
                           "18.01", "18.02", "105.05", "105.06", "106.03",
                           "106.04", "106.05", "106.06", "106.07", "106.08",
                           "145", "9800"), "Salinas",
                       ifelse(NAME %in% c("112.02", "112.03", "112.04"), 
                              "Greenfield", 
                              ifelse(NAME %in% c("111.01", "111.03", "111.04",
                                                 "111.05", "111.06"), "Soledad", 
                                     ifelse(NAME %in% c("113.02", "113.03", "113.05",
                                                        "113.06"), "King City", 
                                            ifelse(NAME %in% c("108.04", "109", "148"), "Gonzales",
                                                   ifelse(NAME %in% c("104"), "Castroville", 
                                                          ifelse(NAME %in% c("101.01", "101.02", "102.02",
                                                                             "103.05", "103.06", "105.01",
                                                                             "105.04", "107.02", "114", "135",
                                                                             "136", "137", "138", "139",
                                                                             "141.02", "141.04", "141.05",
                                                                             "141.08", "141.09", "141.10",
                                                                             "142.01", "142.02", "143.01", 
                                                                             "143.02", "146.01", "147"), 
                                                                 "Other Salinas Valley Cities (Monterey)", 
                                                                 ifelse(NAME %in% c("1101.01", "1101.02", "1102.01",
                                                                                    "1102.02", "1103.01", "1103.02",
                                                                                    "1104.01", "1104.02", "1105.03",
                                                                                    "1105.04", "1105.05", "1105.06",
                                                                                    "1106.01", "1106.02", "1107", 
                                                                                    "1223", "1231", "1233"), 
                                                                        "Other Salinas Valley Cities (Santa Cruz)", "Other")))))))))


# create palette with 9 colors for plot
nb.cols <- 9
mycolors <- colorRampPalette(brewer.pal(5, "Set2"))(nb.cols)

ggplot(CA_tract, aes(fill = fill)) + geom_sf() + 
  scale_fill_manual(values = mycolors) +
  guides(fill = guide_legend(title = "Region")) + 
  theme_minimal(base_size = 10) + 
  annotation_scale() +
  annotation_north_arrow(which_north = "true", 
                         pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"), location = "tr")






