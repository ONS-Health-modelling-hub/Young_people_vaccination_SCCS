library(dplyr)
library(sparklyr)
library(tidyverse)
library(lubridate)
library(splines)
library(ggplot2)
library(sandwich)
library(lmtest)
library(stringr)
library(data.table)
library(survival)
library(gnm)

dir ="cen_dth_gps/Vaccination analysis/Deaths_post_vaccination/Young_people/"
results_folder <- paste0(dir, "Final_results/revisions_vaccinations/")
dir.create(results_folder)
plots_folder <- paste0(results_folder, "/figures_and_tables/")
dir.create(plots_folder)
data_folder <- paste0(dir, "Data_revisions/")


source(paste0(dir, "_Functions.R"))

#--------------------------
# Registrations data read in and variable transformations
#--------------------------


sccs_data <- read.csv(paste0(data_folder, "dataset_for_sccs_by_day_with_unvacc_2022-12-03_registrations.csv"))


# create numbers doses and age groups
sccs_data <- sccs_data %>%
  rename(id = nhsno) %>%
  mutate(current_dose = case_when(vacc_dose == "Unvaccinated" ~ 0,
                                 vacc_dose == "First" ~ 1,
                                 vacc_dose == "Second" ~ 2,
                                 vacc_dose == "Third" ~ 3),
       age_group = case_when(ageinyrs>=12 & ageinyrs<=17 ~ '12-17',
                             ageinyrs>=18 & ageinyrs<=24 ~ '18-24',
                             ageinyrs>=25 & ageinyrs<=29 ~ '25-29'))




# make factors
sccs_data$sex <-as.factor(sccs_data$sex)
sccs_data$age_group <-as.factor(sccs_data$age_group)
sccs_data$vacc_status <- relevel(as.factor(sccs_data$vacc_status), ref="Unvaccinated")


# create the interval value, the length of the week in days as is sometimes not 7
sccs_data_week <- sccs_data %>%
  group_by(id, week, vacc_dose) %>%
  mutate(interval = n()) %>%
  ungroup()

# drop daily columns to make weekly data
sccs_data_week <- sccs_data_week %>%
  select(-day, -vacc_day, -death_day) %>%
  distinct() %>%
  rename(death = death_week) %>%
  data.frame()

# cal_day_fortnight variable
sccs_data_week <- sccs_data_week %>%
  mutate(cal_day_fortnight = as.factor(floor(cal_day_week/14)))

# make 12 week sirks period variables
sccs_data_week <- make_risk_period_variables(sccs_data_week, 12)


# make cardiac deaths only dataset
sccs_data_week_cardiac <- sccs_data_week %>%
  filter(cardiac_cause_of_death == 1)

#--------------------------
# Quick checks
#--------------------------

# number of people in dataset 
n_people <- sccs_data_week %>%
  select(id) %>%
  distinct() %>%
  count()
n_people

# number of people in dataset 
n_people_unvacc <- sccs_data_week %>%
  filter(vacc_status == "Unvaccinated") %>%
  select(id) %>%
  distinct() %>%
  count()
n_people_unvacc



#--------------------------
# Numbers in risk periods for absolute rate calculations
#--------------------------

# number of women with cardiac death in 12 weeks after first dose vaccination
# for absolute rates
n_women_cardiac_vacc <- sccs_data_week_cardiac %>%
  filter(vacc_status == "First" & death == 1 & cardiac_cause_of_death == 1 & sex == 2 & risk_period == "risk") %>%
  count()
n_women_cardiac_vacc

# number of women with cardiac death in 12 weeks after first dose non mRNA vaccination
# for absolute rates
n_women_cardiac_vacc <- sccs_data_week %>%
  filter(vacc_status == "First" & death == 1 & cardiac_cause_of_death == 1 & sex == 2 & risk_period == "risk") %>%
  filter(first_vacc_vector == "not_mRNA_or_unknown") %>%
  count()
n_women_cardiac_vacc

# number of women with cardiac death in 12 weeks after first dose non mRNA vaccination, age 18-24
# for absolute rates
n_women_cardiac_vacc <- sccs_data_week %>%
  filter(vacc_status == "First" & death == 1 & cardiac_cause_of_death == 1 & sex == 2
         & risk_period == "risk") %>%
  filter(first_vacc_vector == "not_mRNA_or_unknown") %>%
  filter(age_group == '18-24') %>%
  count()
n_women_cardiac_vacc

# number of women with cardiac death in 12 weeks after second dose non mRNA vaccination, age 18-24
# for absolute rates
n_women_cardiac_vacc <- sccs_data_week %>%
  filter(vacc_status == "Second" & death == 1 & cardiac_cause_of_death == 1 & sex == 2
         & risk_period == "risk") %>%
  filter(second_vacc_vector == "not_mRNA_or_unknown") %>%
  filter(age_group == '18-24') %>%
  count()
n_women_cardiac_vacc


# number of men with cardiac death in 12 weeks after second dose vaccination
# for absolute rates
n_men_vacc <- sccs_data_week %>%
  filter(vacc_status == "Second" & death == 1 & sex == 1 & cardiac_cause_of_death & risk_period == "risk") %>%
  filter(second_vacc_vector == "mRNA") %>%
  count()
n_men_vacc



#--------------------------
# Vaccine vector tables
#--------------------------


# vaccination types
vacc_vec <- sccs_data_week %>%
  select(id, first_vacc_vector, second_vacc_vector, third_vacc_vector) %>%
  distinct() 

# first and second only
vacc_vec_table <- vacc_vec %>%
  group_by(first_vacc_vector, second_vacc_vector) %>%
  count()
vacc_vec_table

write.csv(vacc_vec_table, paste0(plots_folder, "vacc_vector_table_reg_1_and_2.csv"), row.names = FALSE)

# third dose
third_dose <- vacc_vec %>% group_by(third_vacc_vector) %>% count() 
third_dose
write.csv(third_dose, paste0(plots_folder, "vacc_vector_third_dose.csv"), row.names = FALSE)



# median dosing intervals

dosing <- sccs_data %>% 
  filter(vacc_status != vacc_dose) %>%
  group_by(id, vacc_dose) %>%
  summarise(day_of_next_vacc = max(day) + 1) %>%
  mutate(vacc_dose = as.character(vacc_dose)) %>%
  mutate(vacc = case_when(vacc_dose == "Unvaccinated" ~ "First",
                         vacc_dose == "First" ~ "Second",
                         vacc_dose == "Second" ~ "Third")) %>%
  select(-vacc_dose) %>%
  pivot_wider(names_from = vacc, values_from = day_of_next_vacc) %>%
  mutate(week_1_to_2 = floor((Second - First)/7)) %>%
  mutate(week_2_to_3 = floor((Third - Second)/7))

# median dosing intervals
dose_1_2 <- median(dosing$week_1_to_2, na.rm = TRUE)
dose_2_3 <- median(dosing$week_2_to_3, na.rm = TRUE)

dosing_intervals <- data.frame(dose_1_2 = dose_1_2,
                              dose_2_3 = dose_2_3)
write.csv(dosing_intervals, paste0(plots_folder, "dosing_intervals.csv"), row.names = FALSE)

# number of AZ before and after 7 April 2021 (day 120)
first_dose_az <- sccs_data_week %>% 
  filter(risk_week == "week_01" & current_dose == 1) %>%
  mutate(timing = ifelse(cal_day_week <120, "pre_April_7", "post_April_7")) %>%
  group_by(timing, first_vacc_vector) %>%
  count() %>%
  rename(vacc_vector = first_vacc_vector)
first_dose_az$dose <- "first"

second_dose_az <- sccs_data_week %>% 
  filter(risk_week == "week_01" & current_dose == 2) %>%
  mutate(timing = ifelse(cal_day_week <120, "pre_April_7", "post_April_7")) %>%
  group_by(timing, second_vacc_vector) %>%
  count() %>%
  rename(vacc_vector = second_vacc_vector)
second_dose_az$dose <- "second"



all <- rbind(first_dose_az, second_dose_az)

all

write.csv(all, paste0(plots_folder, "vacc_vector_pre_post_april_7.csv"), row.names = FALSE)


#--------------------------
# HES data read in and variable transformations
#--------------------------


sccs_data_hes <- read.csv(paste0(data_folder, "dataset_for_sccs_by_day_with_unvacc_2022-12-03_hes.csv"))


sccs_data_hes <- sccs_data_hes %>%
  rename(id = nhsno) %>%
  mutate(current_dose = case_when(vacc_dose == "Unvaccinated" ~ 0,
                                 vacc_dose == "First" ~ 1,
                                 vacc_dose == "Second" ~ 2,
                                 vacc_dose == "Third" ~ 3),
       age_group = case_when(ageinyrs>=12 & ageinyrs<=17 ~ '12-17',
                             ageinyrs>=18 & ageinyrs<=24 ~ '18-24',
                             ageinyrs>=25 & ageinyrs<=29 ~ '25-29'))


# make factors
sccs_data_hes$sex <-as.factor(sccs_data_hes$sex)
sccs_data_hes$age_group <-as.factor(sccs_data_hes$age_group)
sccs_data$vacc_status <- relevel(as.factor(sccs_data$vacc_status), ref="Unvaccinated")

# create the interval value, the lenght of the week in days as is sometimes not 7
sccs_data_week_hes <- sccs_data_hes %>%
  group_by(id, week, vacc_dose) %>%
  mutate(interval = n()) %>%
  ungroup()

# drop daily columns to make weekly data
sccs_data_week_hes <- sccs_data_week_hes %>%
  select(-day, -vacc_day, -death_day) %>%
  distinct() %>%
  rename(death = death_week) %>%
  data.frame()

# cal_day_fortnight variable
sccs_data_week_hes <- sccs_data_week_hes %>%
  mutate(cal_day_fortnight = as.factor(floor(cal_day_week/14)))

sccs_data_week_hes <- make_risk_period_variables(sccs_data_week_hes, 12)

#--------------------------
# Quick Checks
#--------------------------

# number of people in dataset 
n_people <- sccs_data_week_hes %>%
  select(id) %>%
  distinct() %>%
  count()
n_people

# number of people in dataset 
n_people_unvacc <- sccs_data_week_hes %>%
  filter(vacc_status == "Unvaccinated") %>%
  select(id) %>%
  distinct() %>%
  count()
n_people_unvacc


#--------------------------
# Numbers in risk periods for absolute rate calculations
#--------------------------

# number of women with hospital death in 12 weeks after first dose non mRNA vaccination, age 18-24
# for absolute rates
n_women_hospital_vacc <- sccs_data_week_hes %>%
  filter(vacc_status == "First" & death == 1 & sex == 2
         & risk_period == "risk") %>%
  filter(first_vacc_vector == "not_mRNA_or_unknown") %>%
  filter(age_group == '18-24') %>%
  count()
n_women_hospital_vacc

# number of women with hospital death in 12 weeks after second dose non mRNA vaccination, age 18-24
# for absolute rates
n_women_hospital_vacc <- sccs_data_week_hes %>%
  filter(vacc_status == "Second" & death == 1 & sex == 2
         & risk_period == "risk") %>%
  filter(second_vacc_vector == "not_mRNA_or_unknown") %>%
  filter(age_group == '18-24') %>%
  count()
n_women_hospital_vacc

# number of men with cardiac death in 12 weeks after second dose vaccination, 18-24
# for absolute rates
n_men_hospital_vacc <- sccs_data_week %>%
  filter(vacc_status == "Second" & death == 1 & sex == 1 & cardiac_cause_of_death & risk_period == "risk") %>%
  filter(second_vacc_vector == "mRNA") %>%
  filter(age_group == '18-24') %>%
  count()
n_men_hospital_vacc

#--------------------------
# Vaccine vector tables
#--------------------------


vacc_vec <- sccs_data_week_hes %>%
  select(id, first_vacc_vector, second_vacc_vector) %>%
  distinct() 

# first and second doses
vacc_vec_table <- vacc_vec %>%
  group_by(first_vacc_vector, second_vacc_vector) %>%
  count()
vacc_vec_table

write.csv(vacc_vec_table, paste0(plots_folder, "vacc_vector_table_hes_1_and_2.csv"), row.names = FALSE)

vacc_vec <- sccs_data_week_hes %>%
  select(id, first_vacc_vector, second_vacc_vector, third_vacc_vector) %>%
  distinct() 

# third dose
third_dose <- vacc_vec %>% group_by(third_vacc_vector) %>% count() 
third_dose
write.csv(third_dose, paste0(plots_folder, "vacc_vector_third_dose_hes.csv"), row.names = FALSE)


# number infected on epistart

n_infect_epistart <- sccs_data_week_hes %>%
  select(id, infect_on_epistart) %>%
  distinct() %>%
  group_by(infect_on_epistart) %>%
  count()
n_infect_epistart

write.csv(n_infect_epistart, paste0(plots_folder, "infect_on_epistart_hes.csv"), row.names = FALSE)

#--------------------------
# Table 1
#--------------------------



sccs_data_week <- sccs_data_week %>%
  mutate(risk_period_for_table = case_when(risk_period == "risk" ~ "12 weeks or less",
                                          risk_period == "baseline" & vacc_status == "Unvaccinated" ~ "Unvaccinated",
                                          TRUE ~ "13+ weeks"))

sccs_data_week_cardiac <- sccs_data_week_cardiac %>%
  mutate(risk_period_for_table = case_when(risk_period == "risk" ~ "12 weeks or less",
                                          risk_period == "baseline" & vacc_status == "Unvaccinated" ~ "Unvaccinated",
                                          TRUE ~ "13+ weeks"))

sccs_data_week_hes <- sccs_data_week_hes %>%
  mutate(risk_period_for_table = case_when(risk_period == "risk" ~ "12 weeks or less",
                                          risk_period == "baseline" & vacc_status == "Unvaccinated" ~ "Unvaccinated",
                                          TRUE ~ "13+ weeks"))

result_all <- get_cat_vars(filter(sccs_data_week, death==1), c('sex', 'age_group', 'vacc_status', 
                                                                  'risk_period_for_table'), total = TRUE)


result_all <-  result_all %>%
  rename('All-cause registered deaths' = Count..n.)
  

result_cardiac <- get_cat_vars(filter(sccs_data_week_cardiac, death==1), c('sex', 
                                              'age_group', 'vacc_status', 'risk_period_for_table'), total = TRUE)
result_cardiac <-  result_cardiac %>%
  rename('Cardiac registered deaths' = Count..n.)

result_hes <- get_cat_vars(filter(sccs_data_week_hes, death==1), c('sex', 
                                              'age_group', 'vacc_status', 'risk_period_for_table'), total = TRUE)
result_hes <-  result_hes %>%
  rename('All-cause hospital deaths' = Count..n.)

summary <- merge(result_all, result_cardiac)
summary <- merge(summary, result_hes)

summary$variable <- as.factor(summary$variable <- relevel(as.factor(summary$variable), ref="Total"))
summary <- summary[order(summary$variable),]
summary$variable <- as.character(summary$variable)

summary <- summary %>%
  mutate(variable = case_when(variable == 'age_group' ~ 'Age group',
                             variable == 'sex' ~ 'Sex',
                             variable == 'vacc_status' ~ 'Most recent vaccination',
                             variable == 'risk_period_for_table' ~ 'Risk period',
                             TRUE ~ variable),
        breakdown = case_when(breakdown == '1' ~ 'Male',
                             breakdown == '2' ~ 'Female',
                           #  breakdown == 'baseline' ~ 'Pre vaccination or 13+ weeks',
                           #  breakdown == 'risk' ~ '12 weeks or less',
                             TRUE ~ breakdown))

summary <- summary[,c('variable', 'breakdown', 'All-cause registered deaths', 'Cardiac registered deaths', 
                      'All-cause hospital deaths')]

summary <- summary %>%
  rename(Characteristic = variable,
        Level = breakdown)

summary

write.csv(summary, paste0(plots_folder, "Table1.csv"), row.names = FALSE)


#--------------------------
# Deaths by calendar week
#--------------------------

sccs_data_with_week <- sccs_data %>%
  mutate(week_from_dec_8 = floor(day/7) +1)

sccs_data_with_week_hes <- sccs_data_hes %>%
  mutate(week_from_dec_8 = floor(day/7) +1)


# check deaths by day to see where to set eos
deaths_by_week <- sccs_data_with_week %>%
  group_by(week_from_dec_8, vacc_status) %>%
  summarise(n_deaths = sum(death_day)) %>%
  mutate(cause = 'All-cause deaths')

deaths_by_week_cardiac <- sccs_data_with_week %>%
  filter(cardiac_cause_of_death == 1) %>%
  group_by(week_from_dec_8, vacc_status) %>%
  summarise(n_deaths = sum(death_day)) %>%
  mutate(cause = 'Cardiac deaths')

deaths_by_week_hes <- sccs_data_with_week_hes %>%
  group_by(week_from_dec_8, vacc_status) %>%
  summarise(n_deaths = sum(death_day)) %>%
  mutate(cause = 'All-cause deaths (HES)')

all_deaths_by_week <- rbind(deaths_by_week, deaths_by_week_cardiac)
all_deaths_by_week <- rbind(all_deaths_by_week, deaths_by_week_hes)



write.csv(all_deaths_by_week, paste0(plots_folder, "SupFig3_deaths_all_weeks.csv"), row.names=FALSE)

all_deaths_by_week$vacc_status <- as.character(all_deaths_by_week$vacc_status)

all_deaths_by_week <- all_deaths_by_week %>%
  mutate(vacc_status = ifelse(cause=="Cardiac deaths" & week_from_dec_8 >= 60, "Suppressed", vacc_status)) %>%
  group_by(week_from_dec_8 , vacc_status, cause) %>%
  mutate(n_deaths = sum(n_deaths)) %>%
  ungroup() %>%
  distinct()



write.csv(all_deaths_by_week, paste0(plots_folder, "SupFig3_deaths_all_weeks_suppressed.csv"), row.names=FALSE)

all_deaths_by_week$vacc_status <- as.factor(all_deaths_by_week$vacc_status)

all_deaths_by_week$vacc_status <- factor(all_deaths_by_week$vacc_status,
                                         levels = c("Suppressed", "Third", "Second", "First", "Unvaccinated"))
levels(all_deaths_by_week$vacc_status) <- c("Suppressed", "Third dose", "Second dose", "First dose", "Unvaccinated")

all_deaths_by_week$cause <- factor(all_deaths_by_week$cause,
                     levels = c("All-cause deaths", "Cardiac deaths", "All-cause deaths (HES)"))
levels(all_deaths_by_week$cause) <- c("All-cause registered deaths", "Cardiac registered deaths", 
                                         "All-cause hospital deaths")

plt <- ggplot(all_deaths_by_week, aes(x=week_from_dec_8, y=n_deaths, fill = vacc_status)) + 
  geom_bar(stat='identity', position = "stack") +
  theme_light() +
  facet_wrap(~cause) +
  xlab("Week since 8 December 2020") +
  ylab("Count of deaths") +
  scale_fill_discrete(name = "Last vaccination received") +
  theme(text=element_text(size=20)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  xlim(0, 70)
plt

ggsave(paste0(plots_folder, "SupFig3_deaths_all_weeks.png"), width = 15, height = 5)





#--------------------------
# Deaths since dose
#--------------------------

weeks_since_dose <- function(data){

  up_to_death <- data %>%
    group_by(id) %>%
    arrange(cal_day_week) %>%
    mutate(has_died = cumsum(death)) %>%
    filter(has_died == 0 | death == 1)

  days_since_dose_1 <- up_to_death %>%
    group_by(id, current_dose, vacc_status) %>%
    summarise(days_since_dose = sum(interval)) %>%
    filter(current_dose != 0) %>%
    group_by(id) %>%
    arrange(current_dose) %>%
    mutate(days_since_dose_1 = cumsum(days_since_dose)) %>%
    mutate(weeks_since_dose = floor(days_since_dose_1/7)) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(max = max(current_dose)) %>%
    ungroup() %>%
    filter(max == current_dose) %>%
    select(-days_since_dose, -max) %>%
    group_by(vacc_status, weeks_since_dose) %>%
    summarise(n = n()) %>%
    mutate(series = "Weeks since first dose")

  days_since_dose_2 <- up_to_death %>%
    group_by(id, current_dose, vacc_status) %>%
    summarise(days_since_dose = sum(interval)) %>%
    filter(current_dose != 0 & current_dose != 1) %>%
    group_by(id) %>%
    arrange(current_dose) %>%
    mutate(days_since_dose_2 = cumsum(days_since_dose)) %>%
    mutate(weeks_since_dose = floor(days_since_dose_2/7)) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(max = max(current_dose)) %>%
    ungroup() %>%
    filter(max == current_dose) %>%
    select(-days_since_dose, -max) %>%
    group_by(vacc_status, weeks_since_dose) %>%
    summarise(n = n()) %>%
    mutate(series = "Weeks since second dose")

  days_since_dose_3 <- up_to_death %>%
    group_by(id, current_dose, vacc_status) %>%
    summarise(days_since_dose = sum(interval)) %>%
    filter(current_dose == 3) %>%
    mutate(days_since_dose_3 = days_since_dose) %>%
    mutate(weeks_since_dose = floor(days_since_dose_3/7)) %>%
    select(-days_since_dose) %>%
    group_by(vacc_status, weeks_since_dose) %>%
    summarise(n = n()) %>%
    mutate(series = "Weeks since third dose")

  weeks_since <- rbind(days_since_dose_1, days_since_dose_2, days_since_dose_3)
  
  return(weeks_since)
}


week_since_reg <- weeks_since_dose(sccs_data_week)
week_since_reg$cause <- "All-cause registered deaths"

week_since_hes <- weeks_since_dose(sccs_data_week_hes)
week_since_hes$cause <- "All-cause hospital deaths"

week_since_cardiac <- weeks_since_dose(sccs_data_week_cardiac)
week_since_cardiac$cause <- "Cardiac registered deaths"

week_since <- rbind(week_since_reg, week_since_hes, week_since_cardiac)

write.csv(week_since, paste0(plots_folder, "SupFig2_deaths_since_dose_stacked.csv"), row.names=FALSE)

week_since$vacc_status <- as.character(week_since$vacc_status)

week_since <- week_since %>%
  mutate(vacc_status = ifelse(cause == "Cardiac registered deaths", "Suppressed", vacc_status)) %>%
  group_by(vacc_status, weeks_since_dose, series, cause) %>%
  mutate(n = sum(n)) %>%
  ungroup() %>%
  distinct()

week_since$vacc_status <- as.factor(week_since$vacc_status)

week_since$cause <- factor(week_since$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", "All-cause hospital deaths"))

week_since$vacc_status <- factor(week_since$vacc_status,
                                         levels = c("Suppressed", "Third", "Second", "First", "Unvaccinated"))
levels(week_since$vacc_status) <- c("Suppressed", "Third dose", "Second dose", "First dose", "Unvaccinated")

write.csv(week_since, paste0(plots_folder, "SupFig2_deaths_since_dose_stacked_suppressed.csv"), row.names=FALSE)


plt <- ggplot(week_since, aes(x=weeks_since_dose, y=n, fill = vacc_status)) + 
  geom_bar(stat='identity', position = "stack") +
  theme_light() +
  facet_grid(series~cause) +
  xlab("Week since dose") +
  ylab("Count of deaths") +
  scale_fill_discrete(name = "Last vaccination received") +
  theme(text=element_text(size=20)) +
  guides(fill = guide_legend(reverse=TRUE))
plt

ggsave(paste0(plots_folder, "deaths_since_dose.png"), width = 16, height = 10)



#--------------------------
# Model 
#--------------------------

df_clean <- create_dataset(sccs_data_week, 12)
df_clean_hes <- create_dataset(sccs_data_week_hes, 12)
df_clean_cardiac <- create_dataset(sccs_data_week_cardiac, 12)



# run model for all doses combined and for individual doses and svae out results
run_model <- function(dataset, suffix, calendar_adjustment){
  
  # registrations dataset
  res_period_dose <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period_dose",
                 x= calendar_adjustment,
                 data = dataset,
                 model =F,
                 save_name = paste0("Risk_dose_results_last_vacc", suffix))
  res_period_dose$period = "period"

  res_week_dose <- eventdepexp_model(outcome = "death",
                 exposure = "risk_week_dose",
                 x= calendar_adjustment,
                 data = dataset,
                 model =F,
                 save_name = paste0("Risk_week_dose_results", suffix))
  res_week_dose$period = "weeks"

  res_week <- eventdepexp_model(outcome = "death",
                 exposure = "risk_week",
                 x= calendar_adjustment,
                 data = dataset,
                 model =F,
                 save_name = paste0("Risk_week_results", suffix))
  res_week$period = "weeks"

  res_period <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period",
                 x= calendar_adjustment,
                 data = dataset,
                 model =F,
                 save_name = paste0("Risk_period_results", suffix))
  res_period$period = "period"
  
  all_data <- do.call(rbind, list(res_period_dose, res_week_dose, res_week, res_period))
  return(all_data)
}



# spline cal adjustment
res_reg <- run_model(df_clean, "_spline_cal_reg", "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))")
res_reg$cause <- "All-cause registered deaths"

res_hes <- run_model(df_clean_hes, "_spline_cal_hes", "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))")
res_hes$cause <- "All-cause hospital deaths"

res_cardiac <- run_model(df_clean_cardiac, "_spline_cal_cardiac", "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))")
res_cardiac$cause <- "Cardiac registered deaths"

all_res <- do.call(rbind, list(res_reg, res_hes, res_cardiac))

all_res <- all_res %>%
    mutate(dose = case_when(grepl("D1",term) ~ "First",
                          grepl("D2",term) ~ "Second",
                          grepl("D3",term) ~ "Third",
                         TRUE ~ 'All doses'),
           parameter = case_when(grepl("week_01",term) ~ "Week 1",
                          grepl("week_02",term) ~ "Week 2",
                          grepl("week_03",term) ~ "Week 3",
                          grepl("week_04",term) ~ "Week 4",
                          grepl("week_05",term) ~ "Week 5",
                          grepl("week_06",term) ~ "Week 6",
                          grepl("week_07",term) ~ "Week 7",
                          grepl("week_08",term) ~ "Week 8",
                          grepl("week_09",term) ~ "Week 9",
                          grepl("week_10",term) ~ "Week 10",
                          grepl("week_11",term) ~ "Week 11",
                          grepl("week_12",term) ~ "Week 12",
                          TRUE ~ 'Weeks 1-12'))

all_res_parameters <- all_res[!grepl("cal_day", all_res$term), ]

write.csv(all_res_parameters, paste0(plots_folder, "/Main_results_risk_period_and_all_weeks_spline_cal.csv"),
         row.names = FALSE)



# fortnight cal adjustment
res_reg <- run_model(df_clean, "_fortnight_cal_reg", "cal_day_fortnight")
res_reg$cause <- "All-cause registered deaths"

res_hes <- run_model(df_clean_hes, "_fortnight_cal_hes", "cal_day_fortnight")
res_hes$cause <- "All-cause hospital deaths"

res_cardiac <- run_model(df_clean_cardiac, "_fortnight_cal_cardiac", "cal_day_fortnight")
res_cardiac$cause <- "Cardiac registered deaths"

all_res <- do.call(rbind, list(res_reg, res_hes, res_cardiac))

all_res <- all_res %>%
    mutate(dose = case_when(grepl("D1",term) ~ "First",
                          grepl("D2",term) ~ "Second",
                          grepl("D3",term) ~ "Third",
                         TRUE ~ 'All doses'),
           parameter = case_when(grepl("week_01",term) ~ "Week 1",
                          grepl("week_02",term) ~ "Week 2",
                          grepl("week_03",term) ~ "Week 3",
                          grepl("week_04",term) ~ "Week 4",
                          grepl("week_05",term) ~ "Week 5",
                          grepl("week_06",term) ~ "Week 6",
                          grepl("week_07",term) ~ "Week 7",
                          grepl("week_08",term) ~ "Week 8",
                          grepl("week_09",term) ~ "Week 9",
                          grepl("week_10",term) ~ "Week 10",
                          grepl("week_11",term) ~ "Week 11",
                          grepl("week_12",term) ~ "Week 12",
                          TRUE ~ 'Weeks 1-12'))

all_res_parameters <- all_res[!grepl("cal_day", all_res$term), ]

write.csv(all_res_parameters, paste0(plots_folder, "/Main_results_risk_period_and_all_weeks_fortnight_cal.csv"),
         row.names = FALSE)



#--------------------------
# Plot IRR
#--------------------------

type <- "_spline_cal"
plot_data <- read.csv(paste0(plots_folder, "/Main_results_risk_period_and_all_weeks", type, ".csv"))

plot_data <- plot_data[plot_data$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))

plot_data$dose <- factor(plot_data$dose,levels = c("All doses", "First", "Second", "Third"))
levels(plot_data$dose) <- c("All doses", "First dose", "Second dose", "Third dose")

plot_data$parameter <- factor(plot_data$parameter,
                     levels = c("Week 12", "Week 11", "Week 10", "Week 9", "Week 8", "Week 7", 
                               "Week 6", "Week 5", "Week 4", "Week 3", "Week 2", "Week 1", "Weeks 1-12"))




write.csv(plot_data, paste0(plots_folder, "Fig1_results_with_spline", type, ".csv"),
         row.names = FALSE)

ggplot(filter(plot_data, dose=="All doses"), aes(x=irr, y=parameter)) + 
    geom_point(aes(colour = period), size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour=period), width=.2) +
    facet_grid(dose~cause) +
    theme_light() +
    xlab("Relative incidence") +
    scale_x_continuous(trans='log2', breaks = c(0.0625/4, 0.0625, 0.25, 1, 4, 16, 64), 
                       labels = c(expression(1/64), expression(1/16), expression(1/4), 
                                  expression(1), expression(4), expression(16), expression(64)), 
                                                                         limits = c(0.0625/8, 128)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=14),
       legend.position = 'none',
       panel.grid.major.y = element_blank())

ggsave(paste0(plots_folder, "Fig1_results_with_spline_all_doses", type, ".png"), width = 13, height = 7)

ggplot(filter(plot_data, dose!="All doses"), aes(x=irr, y=parameter)) + 
    geom_point(aes(colour = period), size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour=period), width=.2) +
    facet_grid(dose~cause) +
    theme_light() +
    xlab("Relative incidence") +
    scale_x_continuous(trans='log2', breaks = c(0.0625/4, 0.0625, 0.25, 1, 4, 16, 64), 
                       labels = c(expression(1/64), expression(1/16), expression(1/4), 
                                  expression(1), expression(4), expression(16), expression(64)), 
                                                                         limits = c(0.0625/8, 128)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=14),
       legend.position = 'none',
       panel.grid.major.y = element_blank())

ggsave(paste0(plots_folder, "Fig1_results_with_spline_separate_doses", type, ".png"), width = 13, height = 18)




#--------------------------
# Different risk periods
#--------------------------

get_all_risk_data <- function(risk_period_length){
 
  df_clean <- create_dataset(sccs_data_week, risk_period_length)
  
  # all_deaths
  all <- eventdepexp_model(outcome = "death",
                   exposure = "risk_period",
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean,
                   model =F)

  all$cause <- "All-cause deaths"

  df_clean_cardiac <- create_dataset(sccs_data_week_cardiac, risk_period_length)
  # cardiac related
  cardiac <- eventdepexp_model(outcome = "death",
                   exposure = "risk_period",
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean_cardiac,
                   model =F)

  cardiac$cause <- "Cardiac deaths"

  
  df_clean_hes <- create_dataset(sccs_data_week_hes, risk_period_length)
  # cardiac related
  hes <- eventdepexp_model(outcome = "death",
                   exposure = "risk_period",
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean_hes,
                   model =F)

  hes$cause <- "All-cause deaths (HES)"
  

  all_data <- rbind(all, cardiac)
  all_data <- rbind(all_data, hes)

  all_data$risk_period_length <- risk_period_length
  
  return(all_data)
}

risk_period_length <- 2
all_data <- get_all_risk_data(risk_period_length)

risk_period_lengths = 3:40
for (i in risk_period_lengths){
  
  print(i)

  combined <- get_all_risk_data(i)

  all_data <- rbind(all_data, combined)
}


all_data

write.csv(all_data, 
          paste0(plots_folder, "SupFig6_different_risk_period_length.csv"),
         row.names = FALSE)

all_data <- read.csv(paste0(plots_folder, "SupFig6_different_risk_period_length.csv"))



all_data <- all_data %>%
  filter(term == 'risk_periodrisk')


all_data$cause <- factor(all_data$cause,
                     levels = c("All-cause deaths", "Cardiac deaths", "All-cause deaths (HES)"))
levels(all_data$cause) <- c("All-cause registered deaths", "Cardiac registered deaths", 
                                         "All-cause hospital deaths")

all_data_8 <- filter(all_data, risk_period_length <=12)

plt <- ggplot(all_data_8, aes(x=irr, y=risk_period_length)) + 
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
    facet_wrap(~cause) +
    theme_light() +
    xlab("Relative incidence") +
    ylab("Risk period length (weeks)") +
    scale_x_continuous(trans='log2', breaks = c(1/4, 1/2, 1, 2, 4), 
                       labels = c("1/4", "1/2", "1", "2", "4"), limits = c(1/8, 8)) +
    scale_y_continuous(limits = c(1.5,12.5), breaks = c(2,4,6,8,10,12), 
                     labels = c(2,4,6,8,10,12)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
       legend.position = 'none',
       panel.grid.major.y = element_blank(),
       panel.spacing.x = unit(1.5, "lines"),
       plot.margin = margin(10, 10, 10, 10))
#  xlim(0.125, 2)

plt

ggsave(paste0(plots_folder, "SupFig6_different_risk_period_length.png"), height=8, width=12)



all_data_long <- all_data

plt <- ggplot(all_data_long, aes(x=irr, y=risk_period_length)) + 
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
    facet_grid(~cause) +
    theme_light() +
    xlab("Relative incidence") +
    ylab("Risk period length (weeks)") +
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
    scale_x_continuous(trans='log2', breaks = c(1/4, 1/2, 1, 2, 4), 
                       labels = c("1/4", "1/2", "1", "2", "4"), limits = c(1/8, 8)) +
  theme(text=element_text(size=20),
       legend.position = 'none',
       panel.grid.major.y = element_blank(),
       panel.spacing.x = unit(1.5, "lines"),
       plot.margin = margin(10, 10, 10, 10))
#  xlim(0.125, 2)

plt

ggsave(paste0(plots_folder, "SupFig6b_different_risk_period_length_long_all.png"), height=8, width=12)

all_data_long <- filter(all_data, risk_period_length <=24)
all_data_long <- all_data_long %>%
  mutate(is_12 = as.factor(ifelse(risk_period_length != 12, 0, 1)))


plt <- ggplot(all_data_long, aes(x=irr, y=risk_period_length, color = is_12)) + 
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
    facet_grid(~cause) +
    theme_light() +
    xlab("Relative incidence") +
    ylab("Risk period length (weeks)") +
    scale_x_continuous(trans='log2', breaks = c(1/4, 1/2, 1, 2, 4), 
                       labels = c("1/4", "1/2", "1", "2", "4"), limits = c(1/8, 8)) +
    scale_y_reverse() +
    scale_color_manual(values=c("black", "orange")) +
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
       legend.position = 'none',
       panel.grid.major.y = element_blank(),
       panel.spacing.x = unit(1.5, "lines"),
       plot.margin = margin(10, 10, 10, 10))
#  xlim(0.125, 2)

plt

ggsave(paste0(plots_folder, "SupFig6b_different_risk_period_length_long_all_24.png"), height=8, width=12)


#--------------------------
# Different risk periods by dose
#--------------------------

get_all_risk_data_dose <- function(risk_period_length, data, cause_name){
  

  df_clean <- create_dataset(data, risk_period_length)
  
  
  # all_deaths
  all <- eventdepexp_model(outcome = "death",
                   exposure = "risk_period_dose",
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean,
                   model =F)

  all$cause <- cause_name
  
  all$risk_period_length <- risk_period_length
  
  return(all)
  
}



get_all_risk_period_data_dose <- function(data, cause){

  risk_period_length <- 2
  all_data <- get_all_risk_data_dose(risk_period_length, data, cause)


  risk_period_lengths = 3:(max(filter(data, current_dose==3)$week) - 1)
  for (i in risk_period_lengths){
  
    print(i)

    combined <- get_all_risk_data_dose(i, data, cause)

    all_data <- rbind(all_data, combined)
  }
  return(all_data)
}

data_all <- get_all_risk_period_data_dose(sccs_data_week, "All-cause deaths")

data_cardiac <- get_all_risk_period_data_dose(sccs_data_week_cardiac, "Cardiac deaths")

data_hes <- get_all_risk_period_data_dose(sccs_data_week_hes, "All-cause deaths (HES)")

all_data <- rbind(data_all, data_cardiac)
all_data <- rbind(all_data, data_hes)


all_data <- all_data %>%
  mutate(dose = case_when(grepl("D1",term) ~ "First",
                          grepl("D2",term) ~ "Second",
                          grepl("D3",term) ~ "Third",
                         TRUE ~ 'All doses'))

all_data$cause <- factor(all_data$cause,
                     levels = c("All-cause deaths", "Cardiac deaths", "All-cause deaths (HES)"))
levels(all_data$cause) <- c("All-cause registered deaths", "Cardiac registered deaths", 
                                         "All-cause hospital deaths")

write.csv(all_data, paste0(plots_folder, "/SupFig6b_different_risk_period_length_by_doses.csv"), row.names = FALSE)

all_data_dose <- all_data %>%
  filter(!grepl("cal_day_week", term))

all_data_dose$dose <- factor(all_data_dose$dose, levels = c("All doses", "First", "Second", "Third"))
levels(all_data_dose$dose) <- c("All doses", "First dose", "Second dose", "Third dose")

all_data_long <- all_data_dose[all_data_dose$irr_low_ci >= 1E-11, ]

plt <- ggplot(all_data_long, aes(x=irr, y=risk_period_length)) + 
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
    facet_grid(dose~cause, scales="free_x") +
    theme_light() +
    xlab("Relative incidence") +
    ylab("Risk period length (weeks)") +
    scale_x_continuous(trans='log2', breaks = c(0.0625, 0.25, 1, 4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(0.0625/4, 64)) +
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
       legend.position = 'none',
       panel.grid.major.y = element_blank(),
       panel.spacing.x = unit(1.5, "lines"),
       plot.margin = margin(10, 10, 10, 10))
#  xlim(0.125, 2)

plt

ggsave(paste0(plots_folder, "SupFig6b_different_risk_period_length_by_doses.png"), height=18, width=12)



all_data_short <- filter(all_data_long, risk_period_length <=24)
all_data_short <- all_data_short %>%
  mutate(is_12 = as.factor(ifelse(risk_period_length != 12, 0, 1)))

plt <- ggplot(all_data_short, aes(x=irr, y=risk_period_length, color = is_12)) + 
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
    facet_grid(dose~cause, scales="free_x") +
    theme_light() +
    xlab("Relative incidence") +
    ylab("Risk period length (weeks)") +
    scale_x_continuous(trans='log2', breaks = c(0.0625, 0.25, 1, 4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(0.0625/4, 64)) +
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
    scale_y_reverse() +
    scale_color_manual(values=c("black", "orange")) +
  theme(text=element_text(size=20),
       legend.position = 'none',
       panel.grid.major.y = element_blank(),
       panel.spacing.x = unit(1.5, "lines"),
       plot.margin = margin(10, 10, 10, 10))
#  xlim(0.125, 2)

plt

ggsave(paste0(plots_folder, "SupFig6b_different_risk_period_length_by_doses_24.png"), height=18, width=12)



#--------------------------
# Age sex breakdown
#--------------------------

exposure <- "risk_period"

sex_parameters_all <- data_for_risk_period_breakdown(df_clean, 'sex', 'All-cause deaths', exposure)
sex_parameters_cardiac <- data_for_risk_period_breakdown(df_clean_cardiac, 'sex', 'Cardiac deaths', exposure)
sex_parameters_hes <- data_for_risk_period_breakdown(df_clean_hes, 'sex', 'All-cause deaths (HES)', exposure)
sex_parameters <- rbind(sex_parameters_all, sex_parameters_cardiac, sex_parameters_hes)

age_group_parameters_all <- data_for_risk_period_breakdown(df_clean, 'age_group', 'All-cause deaths', exposure)
age_group_parameters_cardiac <- data_for_risk_period_breakdown(df_clean_cardiac, 'age_group', 'Cardiac deaths', exposure)
age_group_parameters_hes <- data_for_risk_period_breakdown(df_clean_hes, 'age_group', 'All-cause deaths (HES)', exposure)
age_group_parameters <- rbind(age_group_parameters_all, age_group_parameters_cardiac, age_group_parameters_hes)

# combine and save breakdown data
all_breakdowns1 <- rbind(sex_parameters, age_group_parameters)

exposure <- "risk_period_dose"

sex_parameters_all <- data_for_risk_period_breakdown(df_clean, 'sex', 'All-cause deaths', exposure)
sex_parameters_cardiac <- data_for_risk_period_breakdown(df_clean_cardiac, 'sex', 'Cardiac deaths', exposure)
sex_parameters_hes <- data_for_risk_period_breakdown(df_clean_hes, 'sex', 'All-cause deaths (HES)', exposure)
sex_parameters <- rbind(sex_parameters_all, sex_parameters_cardiac, sex_parameters_hes)

age_group_parameters_all <- data_for_risk_period_breakdown(df_clean, 'age_group', 'All-cause deaths', exposure)
age_group_parameters_cardiac <- data_for_risk_period_breakdown(df_clean_cardiac, 'age_group', 'Cardiac deaths', exposure)
age_group_parameters_hes <- data_for_risk_period_breakdown(df_clean_hes, 'age_group', 'All-cause deaths (HES)', exposure)
age_group_parameters <- rbind(age_group_parameters_all, age_group_parameters_cardiac, age_group_parameters_hes)

# combine and save breakdown data
all_breakdowns2 <- rbind(sex_parameters, age_group_parameters)

all_breakdowns <- rbind(all_breakdowns1, all_breakdowns2)


write.csv(all_breakdowns, 
          paste0(results_folder, "breakdown_results_together.csv"),
         row.names = FALSE)

all_breakdowns$breakdown_value <- recode_factor(all_breakdowns$breakdown_value, '1' = "Male", "2" = "Female")

all_breakdowns <- all_breakdowns[!grepl("cal_day_week", all_breakdowns$term), ]


all_breakdowns$cause <- factor(all_breakdowns$cause,
                     levels = c("All-cause deaths", "Cardiac deaths", "All-cause deaths (HES)"))
levels(all_breakdowns$cause) <- c("All-cause registered deaths", "Cardiac registered deaths", 
                                         "All-cause hospital deaths")

all_breakdowns$dose <- factor(all_breakdowns$dose,
                     levels = c("All doses", "First", "Second", "Third"))
levels(all_breakdowns$dose) <- c("All doses", "First dose", "Second dose", "Third dose")


write.csv(all_breakdowns, paste0(plots_folder, "/Fig2_risk_period_with_breakdowns.csv"), row.names = FALSE)


all_breakdowns <- all_breakdowns[all_breakdowns$irr_low_ci >= 1E-12, ]



plt <- ggplot(all_breakdowns, aes(x=irr, y=breakdown_value)) + 
geom_point(aes(colour = breakdown), size=3) +
geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour = breakdown), width = .2) +
facet_grid(dose~cause) +
theme_light() +
xlab("Relative incidence") +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
geom_vline(xintercept=1, linetype='dashed', colour='grey') +
scale_x_continuous(trans='log2', breaks = c(0.0625, 0.25, 1, 4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(0.0625/2, 32)) +
theme(text=element_text(size=20),
        axis.title.y=element_blank(),
 #      strip.text.y = element_blank(),
     legend.position = 'none',
     panel.grid.major.y = element_blank())
  
plt



ggsave(paste0(plots_folder, "FIG2_risk_period_with_breakdowns.png"), height=12, width=10)





#--------------------------
# Fortnightly calendar adjustment
#--------------------------

type <- "_fortnight_cal"
plot_data <- read.csv(paste0(plots_folder, "/Main_results_risk_period_and_all_weeks", type, ".csv"))

plot_data <- plot_data[plot_data$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))

plot_data$dose <- factor(plot_data$dose,
                     levels = c("All doses", "First", "Second", "Third"))
levels(plot_data$dose) <- c("All doses", "First dose", "Second dose", "Third dose")

plot_data$parameter <- factor(plot_data$parameter,
                     levels = c("Week 12", "Week 11", "Week 10", "Week 9", "Week 8", "Week 7", 
                               "Week 6", "Week 5", "Week 4", "Week 3", "Week 2", "Week 1", "Weeks 1-12"))


write.csv(plot_data, paste0(plots_folder, "Fig1_results_with_fortnight", type, ".csv"),
         row.names = FALSE)

ggplot(filter(plot_data, dose=="All doses"), aes(x=irr, y=parameter)) + 
    geom_point(aes(colour = period), size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour=period), width=.2) +
    facet_grid(dose~cause) +
    theme_light() +
    xlab("Relative incidence") +
    scale_x_continuous(trans='log2', breaks = c(0.0625/4, 0.0625, 0.25, 1, 4, 16, 64), 
                       labels = c("1/64", "1/16", "1/4", 
                                  "1", "4", "16", "64"), 
                                limits = c(0.0625/8, 128)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=14),
       legend.position = 'none',
       panel.grid.major.y = element_blank())

ggsave(paste0(plots_folder, "Fig1_results_with_fortnight_all_doses", type, ".png"), width = 13, height = 7)

ggplot(filter(plot_data, dose!="All doses"), aes(x=irr, y=parameter)) + 
    geom_point(aes(colour = period), size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour=period), width=.2) +
    facet_grid(dose~cause) +
    theme_light() +
    xlab("Relative incidence") +
    scale_x_continuous(trans='log2', breaks = c(0.0625/4, 0.0625, 0.25, 1, 4, 16, 64), 
                       labels = c("1/64", "1/16", "1/4", 
                                  "1", "4", "16", "64"), 
                                limits = c(0.0625/8, 128)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=14),
       legend.position = 'none',
       panel.grid.major.y = element_blank())

ggsave(paste0(plots_folder, "Fig1_results_with_fortnight_separate_doses", type, ".png"), width = 13, height = 18)




#--------------------------
# Exclude infected on epistart
#--------------------------


sccs_data_week_hes_not_infect_on_epistart <- filter(sccs_data_week_hes, infect_on_epistart == 0)

df_clean_hes <- create_dataset(sccs_data_week_hes_not_infect_on_epistart, 12)


res_hes <- run_model(df_clean_hes, "_spline_cal_hes", "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))")
res_hes$cause <- "All-cause hospital deaths"


all_res <- res_hes

all_res <- all_res %>%
    mutate(dose = case_when(grepl("D1",term) ~ "First",
                          grepl("D2",term) ~ "Second",
                          grepl("D3",term) ~ "Third",
                         TRUE ~ 'All doses'),
           parameter = case_when(grepl("week_01",term) ~ "Week 1",
                          grepl("week_02",term) ~ "Week 2",
                          grepl("week_03",term) ~ "Week 3",
                          grepl("week_04",term) ~ "Week 4",
                          grepl("week_05",term) ~ "Week 5",
                          grepl("week_06",term) ~ "Week 6",
                          grepl("week_07",term) ~ "Week 7",
                          grepl("week_08",term) ~ "Week 8",
                          grepl("week_09",term) ~ "Week 9",
                          grepl("week_10",term) ~ "Week 10",
                          grepl("week_11",term) ~ "Week 11",
                          grepl("week_12",term) ~ "Week 12",
                          TRUE ~ 'Weeks 1-12'))

all_res_parameters <- all_res[!grepl("cal_day", all_res$term), ]

write.csv(all_res_parameters, paste0(plots_folder, "/Main_results_risk_period_and_all_weeks_spline_cal_exclude_infect_on_epistart.csv"),
         row.names = FALSE)


type <- "_spline_cal"
plot_data <- read.csv(paste0(plots_folder, "/Main_results_risk_period_and_all_weeks_spline_cal_exclude_infect_on_epistart.csv"))

plot_data <- plot_data[plot_data$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))

plot_data$dose <- factor(plot_data$dose,
                     levels = c("All doses", "First", "Second", "Third"))
levels(plot_data$dose) <- c("All doses", "First dose", "Second dose", "Third dose")

plot_data$parameter <- factor(plot_data$parameter,
                     levels = c("Week 12", "Week 11", "Week 10", "Week 9", "Week 8", "Week 7", 
                               "Week 6", "Week 5", "Week 4", "Week 3", "Week 2", "Week 1", "Weeks 1-12"))


write.csv(plot_data, paste0(plots_folder, "Fig1_results_with_spline_exclude_infect_on_epistart", type, ".csv"),
         row.names = FALSE)

ggplot(plot_data, aes(x=irr, y=parameter)) + 
    geom_point(aes(colour = period), size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour=period), width=.2) +
    facet_grid(~dose) +
    theme_light() +
    xlab("Relative incidence") +
    scale_x_continuous(trans='log2', breaks = c(0.0625/4, 0.0625, 0.25, 1, 4, 16, 64), 
                       labels = c("1/64", "1/16", "1/4", 
                                  "1", "4", "16", "64"), 
                                limits = c(0.0625/8, 128)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=14),
       legend.position = 'none',
       panel.grid.major.y = element_blank())

ggsave(paste0(plots_folder, "Fig1_results_with_spline_all_doses_exclude_infect_on_epistart", type, ".png"), width = 16, height = 7)


#--------------------------
# Vaccine type brakdown
#--------------------------

counts_vacc_vecc <- function(data){
  data <- data %>%
    filter(death == 1) %>%
    mutate(mRNA_any = ifelse(!is.na(first_vacc_vector) & first_vacc_vector == "mRNA" | 
                           !is.na(second_vacc_vector) & second_vacc_vector == "mRNA" |
                           !is.na(third_vacc_vector) & third_vacc_vector == "mRNA", 1, 0),
         not_mRNA_any = ifelse(!is.na(first_vacc_vector) & first_vacc_vector == "not_mRNA_or_unknown" | 
                           !is.na(second_vacc_vector) & second_vacc_vector == "not_mRNA_or_unknown" |
                           !is.na(third_vacc_vector) & third_vacc_vector == "not_mRNA_or_unknown", 1, 0))

  counts_mRNA <- data %>% group_by(sex, mRNA_any) %>% count()
  counts_mRNA <- counts_mRNA %>% 
    filter(mRNA_any ==1) %>%
    rename(vacc_vecc = mRNA_any) %>%
    mutate(vacc_vecc = ifelse(vacc_vecc==1, "mRNA"))
  
  counts_non_mRNA <- data %>% group_by(sex, not_mRNA_any) %>% count()
  counts_non_mRNA <- counts_non_mRNA %>% 
    filter(not_mRNA_any ==1) %>%
    rename(vacc_vecc = not_mRNA_any) %>%
    mutate(vacc_vecc = ifelse(vacc_vecc==1, "not_mRNA"))
  
  counts <- rbind(counts_mRNA, counts_non_mRNA)
  return(counts)

}

counts_reg <- counts_vacc_vecc(sccs_data_week)
counts_reg$type <- "All-cause registered deaths"
counts_cardiac <- counts_vacc_vecc(sccs_data_week_cardiac)
counts_cardiac$type <- "Cardiac registered deaths"
counts_hes <- counts_vacc_vecc(sccs_data_week_hes)
counts_hes$type <- "All-cause hospital deaths"

all_counts <- rbind(counts_reg, counts_cardiac, counts_hes)
data.frame(all_counts)
write.csv(all_counts, paste0(plots_folder, "/vacc_vecc_counts.csv"),row.names = FALSE)

recount_doses <- function(data, vacc_vecc_type){
  # filter to people who have at least one of that type of vaccine, plus unvaccinated

  data <- data %>%
    mutate(mRNA_any = ifelse(!is.na(first_vacc_vector) & first_vacc_vector == "mRNA" | 
                           !is.na(second_vacc_vector) & second_vacc_vector == "mRNA" |
                           !is.na(third_vacc_vector) & third_vacc_vector == "mRNA", 1, 0),
         not_mRNA_any = ifelse(!is.na(first_vacc_vector) & first_vacc_vector == "not_mRNA_or_unknown" | 
                           !is.na(second_vacc_vector) & second_vacc_vector == "not_mRNA_or_unknown" |
                           !is.na(third_vacc_vector) & third_vacc_vector == "not_mRNA_or_unknown", 1, 0))


  # select people who have at least one dose of that vector or are unvaccinated
  if (vacc_vecc_type == "mRNA"){
  sccs_data_week_vacc <- filter(data, mRNA_any == 1 | vacc_status == "Unvaccinated")
    vacc_vecc_other <- "not_mRNA_or_unknown"
  }else if(vacc_vecc_type == "not_mRNA_or_unknown"){
    sccs_data_week_vacc <- filter(data, not_mRNA_any == 1 | vacc_status == "Unvaccinated")
      vacc_vecc_other <- "mRNA"
  }

# relabel doses in the order for the vector of interest only

  sccs_data_week_vacc <- sccs_data_week_vacc %>%
    mutate(current_dose_recoded = case_when(# two combinations of 1 vaccine
                                          current_dose == 1 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_type ~ 1, 
                                           current_dose == 1 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_other ~ 0,
                                           # 4 combinations of 2 vaccines
                                           current_dose == 2 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_type & 
                                          !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_type ~ 2,
                                        current_dose == 2 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_other & 
                                        !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_type ~ 1,
                                         current_dose == 2 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_type & 
                                          !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_other ~ 1,
                                         current_dose == 2 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_other & 
                                          !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_other ~ 0,
                              # 8 combinations of 3 vaccines
                                        current_dose == 3 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_type & 
                                      !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_type &
                                      !is.na(third_vacc_vector) & third_vacc_vector == vacc_vecc_type ~ 3,
                                        current_dose == 3 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_other & 
                                      !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_type &
                                      !is.na(third_vacc_vector) & third_vacc_vector == vacc_vecc_type ~ 2,
                                        current_dose == 3 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_type & 
                                      !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_other &
                                      !is.na(third_vacc_vector) & third_vacc_vector == vacc_vecc_type ~ 2,
                                        current_dose == 3 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_other & 
                                      !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_other &
                                      !is.na(third_vacc_vector) & third_vacc_vector == vacc_vecc_type ~ 1,
                                        current_dose == 3 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_type & 
                                      !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_type &
                                      !is.na(third_vacc_vector) & third_vacc_vector == vacc_vecc_other ~ 2,
                                        current_dose == 3 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_type & 
                                      !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_other &
                                      !is.na(third_vacc_vector) & third_vacc_vector == vacc_vecc_other ~ 1,
                                        current_dose == 3 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_other & 
                                      !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_type &
                                      !is.na(third_vacc_vector) & third_vacc_vector == vacc_vecc_other ~ 1,
                                        current_dose == 3 & !is.na(first_vacc_vector) & first_vacc_vector == vacc_vecc_other & 
                                      !is.na(second_vacc_vector) & second_vacc_vector == vacc_vecc_other &
                                      !is.na(third_vacc_vector) & third_vacc_vector == vacc_vecc_other ~ 0,
                                      TRUE ~ 0))
    
  sccs_data_week_vacc <- sccs_data_week_vacc %>%
      select(-current_dose) %>%
      rename(current_dose = current_dose_recoded)
    
    
    return(sccs_data_week_vacc)
    }



  
# run model for all doses combined and for individual doses and svae out results
run_model_period_dose <- function(vacc_vecc_type, all_data, hes_data, cardiac_data, risk_period=12){
  
  sccs_data_week_vacc <- recount_doses(all_data, vacc_vecc_type)    
  df_clean <- create_dataset(sccs_data_week_vacc, risk_period)

  sccs_data_week_hes_vacc <- recount_doses(hes_data, vacc_vecc_type)  
  df_clean_hes <- create_dataset(sccs_data_week_hes_vacc, risk_period)

  sccs_data_week_cardiac_vacc <- recount_doses(cardiac_data, vacc_vecc_type)  
  df_clean_cardiac <- create_dataset(sccs_data_week_cardiac_vacc, risk_period)
  
  # registrations dataset
  res_reg <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period_dose",
                 x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                 data = df_clean,
                 model =F)
  res_reg$cause <- "All-cause registered deaths"
  
  res_hes <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period_dose",
                 x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                 data = df_clean_hes,
                 model =F)
  res_hes$cause <- "All-cause hospital deaths"
  
  res_cardiac <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period_dose",
                 x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                 data = df_clean_cardiac,
                 model =F)
  res_cardiac$cause <- "Cardiac registered deaths"
  
  res_reg_period <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period",
                 x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                 data = df_clean,
                 model =F)
  res_reg_period$cause <- "All-cause registered deaths"
  
  res_hes_period <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period",
                 x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                 data = df_clean_hes,
                 model =F)
  res_hes_period$cause <- "All-cause hospital deaths"
  
  res_cardiac_period <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period",
                 x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                 data = df_clean_cardiac,
                 model =F)
  res_cardiac_period$cause <- "Cardiac registered deaths"

  all_res <- do.call(rbind, list(res_reg, res_hes, res_cardiac, 
                                 res_reg_period, res_hes_period, res_cardiac_period))

  
  all_res$vacc_vecc_type <- vacc_vecc_type

  return(all_res)
}


# all 
  
mRNA <- run_model_period_dose("mRNA", sccs_data_week, sccs_data_week_hes, sccs_data_week_cardiac)
non_mRNA <- run_model_period_dose("not_mRNA_or_unknown", sccs_data_week, sccs_data_week_hes, sccs_data_week_cardiac)
  
all_res <- rbind(mRNA, non_mRNA)

all_res <- all_res %>%
    mutate(dose = case_when(grepl("D1",term) ~ "First",
                          grepl("D2",term) ~ "Second",
                          grepl("D3",term) ~ "Third",
                         TRUE ~ 'All doses'))

all_res_parameters <- all_res[!grepl("cal_day", all_res$term), ]

write.csv(all_res_parameters, paste0(plots_folder, "/vacc_vecc_results_all.csv"),
         row.names = FALSE)

plot_data <- all_res_parameters[all_res_parameters$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))
  
plot_data$dose <- factor(plot_data$dose,
                                         levels = c("All doses", "First", "Second", "Third"))
levels(plot_data$dose) <- c("All doses", "First dose", "Second dose", "Third dose")
  
plot_data <- filter(plot_data, dose != "Third dose")
  
plot_data$vacc_vecc_type <- factor(plot_data$vacc_vecc_type,
                     levels = c("not_mRNA_or_unknown", "mRNA"))
levels(plot_data$vacc_vecc_type) <- c("Not mRNA or unknown", "mRNA")
  
plt <- ggplot(plot_data, aes(x=irr, y=vacc_vecc_type)) + 
geom_point(aes(colour = vacc_vecc_type), size=3) +
geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour = vacc_vecc_type), width = .2) +
facet_grid(dose~cause) +
theme_light() +
xlab("Relative incidence") +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
geom_vline(xintercept=1, linetype='dashed', colour='grey') +
scale_x_continuous(trans='log2', breaks = c(1/16, 0.25, 1, 4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(1/16, 16)) +
theme(text=element_text(size=20),
        axis.title.y=element_blank(),
 #      strip.text.y = element_blank(),
     legend.position = 'none',
     panel.grid.major.y = element_blank())
  
plt

ggsave(paste0(plots_folder, "vacc_vecc.png"), width = 13, height = 7)
  
  
# males 
sccs_data_week_male <- filter(sccs_data_week, sex == 1)
sccs_data_week_hes_male <- filter(sccs_data_week_hes, sex == 1)
sccs_data_week_cardiac_male <- filter(sccs_data_week_cardiac, sex == 1)

mRNA <- run_model_period_dose("mRNA", sccs_data_week_male, sccs_data_week_hes_male, sccs_data_week_cardiac_male)
non_mRNA <- run_model_period_dose("not_mRNA_or_unknown", sccs_data_week_male, sccs_data_week_hes_male, sccs_data_week_cardiac_male)
  
all_res <- rbind(mRNA, non_mRNA)

all_res <- all_res %>%
    mutate(dose = case_when(grepl("D1",term) ~ "First",
                          grepl("D2",term) ~ "Second",
                          grepl("D3",term) ~ "Third",
                         TRUE ~ 'All doses'))

all_res_parameters <- all_res[!grepl("cal_day", all_res$term), ]

write.csv(all_res_parameters, paste0(plots_folder, "/vacc_vecc_results_male.csv"),
         row.names = FALSE)

plot_data <- all_res_parameters[all_res_parameters$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))

plot_data$dose <- factor(plot_data$dose,
                                         levels = c("All doses", "First", "Second", "Third"))
levels(plot_data$dose) <- c("All doses", "First dose", "Second dose", "Third dose")
  
plot_data <- filter(plot_data, dose != "Third dose")
  
plot_data$vacc_vecc_type <- factor(plot_data$vacc_vecc_type,
                     levels = c("not_mRNA_or_unknown", "mRNA"))
levels(plot_data$vacc_vecc_type) <- c("Not mRNA or unknown", "mRNA")


plt <- ggplot(plot_data, aes(x=irr, y=vacc_vecc_type)) + 
geom_point(aes(colour = vacc_vecc_type), size=3) +
geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour = vacc_vecc_type), width = .2) +
facet_grid(dose~cause) +
theme_light() +
xlab("Relative incidence") +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
geom_vline(xintercept=1, linetype='dashed', colour='grey') +
scale_x_continuous(trans='log2', breaks = c(1/16, 0.25, 1, 4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(1/16, 16)) +
theme(text=element_text(size=20),
        axis.title.y=element_blank(),
 #      strip.text.y = element_blank(),
     legend.position = 'none',
     panel.grid.major.y = element_blank())
  
plt

ggsave(paste0(plots_folder, "vacc_vecc_male.png"), width = 13, height = 7)
  
  
# females
  
sccs_data_week_female <- filter(sccs_data_week, sex == 2)
sccs_data_week_hes_female <- filter(sccs_data_week_hes, sex == 2)
sccs_data_week_cardiac_female <- filter(sccs_data_week_cardiac, sex == 2)

mRNA <- run_model_period_dose("mRNA", sccs_data_week_female, sccs_data_week_hes_female, sccs_data_week_cardiac_female)
non_mRNA <- run_model_period_dose("not_mRNA_or_unknown", sccs_data_week_female, sccs_data_week_hes_female, sccs_data_week_cardiac_female)
  
  
all_res <- rbind(mRNA, non_mRNA)

all_res <- all_res %>%
    mutate(dose = case_when(grepl("D1",term) ~ "First",
                          grepl("D2",term) ~ "Second",
                          grepl("D3",term) ~ "Third",
                         TRUE ~ 'All doses'))

all_res_parameters <- all_res[!grepl("cal_day", all_res$term), ]

write.csv(all_res_parameters, paste0(plots_folder, "/vacc_vecc_results_female.csv"),
         row.names = FALSE)

plot_data <- all_res_parameters[all_res_parameters$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))

plot_data$dose <- factor(plot_data$dose,levels = c("All doses", "First", "Second", "Third"))
levels(plot_data$dose) <- c("All doses", "First dose", "Second dose", "Third dose")
  
plot_data <- filter(plot_data, dose != "Third dose")

plot_data$vacc_vecc_type <- factor(plot_data$vacc_vecc_type,
                     levels = c("not_mRNA_or_unknown", "mRNA"))
levels(plot_data$vacc_vecc_type) <- c("Not mRNA or unknown", "mRNA")


plt <- ggplot(plot_data, aes(x=irr, y=vacc_vecc_type)) + 
geom_point(aes(colour = vacc_vecc_type), size=3) +
geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour = vacc_vecc_type), width = .2) +
facet_grid(dose~cause) +
theme_light() +
xlab("Relative incidence") +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
geom_vline(xintercept=1, linetype='dashed', colour='grey') +
scale_x_continuous(trans='log2', breaks = c(1/16, 0.25, 1, 4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(1/16, 16)) +
theme(text=element_text(size=20),
        axis.title.y=element_blank(),
 #      strip.text.y = element_blank(),
     legend.position = 'none',
     panel.grid.major.y = element_blank())
  
plt

ggsave(paste0(plots_folder, "vacc_vecc_female.png"), width = 13, height = 7)

  
  ## with different risk periods

mRNA_male_risk_period <- lapply(seq(6,20, 2), function(i){
  
  mRNA <- run_model_period_dose("mRNA", 
                                sccs_data_week_male, 
                                sccs_data_week_hes_male, 
                                sccs_data_week_cardiac_male,
                               risk_period = i)%>%
  mutate(risk_period = i )
  return(mRNA)
})

  mRNA_male_risk_period_df <- do.call(rbind, mRNA_male_risk_period)

  d2_male_cardiac <-  mRNA_male_risk_period_df %>% 
                  filter(cause == "Cardiac registered deaths",
                        term == "risk_period_doseD2")
  
  # save out
write.csv(d2_male_cardiac, paste0(plots_folder, "/dose2_male_cardiac.csv"),
         row.names = FALSE)

  d2_male_cardiac <- d2_male_cardiac %>%
  mutate(is_12 = as.factor(ifelse(risk_period != 12, 0, 1)))
  
  ggplot( data = d2_male_cardiac, aes( y = risk_period, x = irr, colour=is_12))+
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
  geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme_light() +
  xlab("Relative incidence") +
  ylab("Risk period length (weeks)")+
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+

scale_x_continuous(trans='log2', breaks = c(0.25, 0.5, 1, 2, 4), 
                       labels = c("1/4", "1/2", "1", "2", "4"), limits = c(1/4, 4)) +
      scale_y_reverse() +
      scale_color_manual(values=c("black", "orange")) +
theme(text=element_text(size=20),
     #   axis.title.y=element_blank(),
      # strip.text.y = element_blank(),
     legend.position = 'none',
     panel.grid.major.y = element_blank())

ggsave(paste0(plots_folder, "cardiac_male_d2_risk_period_length.png"), width = 3, height = 6)

  
 