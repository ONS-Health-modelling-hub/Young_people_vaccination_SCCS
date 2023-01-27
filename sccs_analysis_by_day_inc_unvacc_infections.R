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
data_folder <- paste0(dir, "Data_revisions/")
source(paste0(dir, "_Functions.R"))

#infection_cat <- "vaccinated_at_infection"
infection_cat <- "unvaccinated_at_infection"
#infection_cat <- "all"



if (infection_cat == "vaccinated_at_infection"){
  results_folder <- paste0(dir, "Final_results/revisions_infections_vaccinated/")
  dir.create(results_folder)
}else if (infection_cat == "unvaccinated_at_infection"){
  results_folder <- paste0(dir, "Final_results/revisions_infections_unvaccinated/")
  dir.create(results_folder)
}else if (infection_cat == "all"){
  results_folder <- paste0(dir, "Final_results/revisions_infections_all/")
  dir.create(results_folder)
}

plots_folder <- paste0(results_folder, "/figures_and_tables/")
dir.create(plots_folder)


#--------------------------
# Deaths data
#--------------------------
sccs_data <- read.csv(paste0(data_folder, "dataset_for_sccs_by_day_with_unvacc_infections_2022-12-03_registrations.csv"))


# make age group and infection variables
sccs_data <- sccs_data %>%
  rename(id = nhsno) %>%
  mutate(age_group = case_when(ageinyrs>=12 & ageinyrs<=17 ~ '12-17',
                             ageinyrs>=18 & ageinyrs<=24 ~ '18-24',
                             ageinyrs>=25 & ageinyrs<=29 ~ '25-29'),
        had_infection = ifelse(!is.na(max_infection_date), 1, 0),
        infection_number = ifelse(infected == "Infected", 1, 0),
        infect_status = case_when(had_infection==0 ~ "No positive test",
                                 vacc_at_infection == 1 ~ "Positive test (vaccinated)",
                                vacc_at_infection == 0 ~ "Positive test (unvaccinated)"))

# make factors
sccs_data$sex <-as.factor(sccs_data$sex)
sccs_data$age_group <-as.factor(sccs_data$age_group)

# create the interval value, the lenght of the week in days as is sometimes not 7
sccs_data_week <- sccs_data %>%
  group_by(id, week, infected) %>%
  mutate(interval = n()) %>%
  ungroup()

# drop daily columns to make weekly data
sccs_data_week <- sccs_data_week %>%
  select(-day, -infect_day, -death_day) %>%
  distinct() %>%
  rename(death = death_week) %>%
  data.frame()

# cal_day_fortnight variable
sccs_data_week <- sccs_data_week %>%
  mutate(cal_day_fortnight = as.factor(floor(cal_day_week/14)))

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




#--------------------------
# HES data
#--------------------------

# hes data
sccs_data_hes <- read.csv(paste0(data_folder, "dataset_for_sccs_by_day_with_unvacc_infections_2022-12-03_hes.csv"))

# make age group and infection variables
sccs_data_hes <- sccs_data_hes %>%
  rename(id = nhsno) %>%
  mutate(age_group = case_when(ageinyrs>=12 & ageinyrs<=17 ~ '12-17',
                             ageinyrs>=18 & ageinyrs<=24 ~ '18-24',
                             ageinyrs>=25 & ageinyrs<=29 ~ '25-29'),
        had_infection = ifelse(!is.na(max_infection_date), 1, 0),
        infection_number = ifelse(infected == "Infected", 1, 0),
        infect_status = case_when(had_infection==0 ~ "No positive test",
                                 vacc_at_infection == 1 ~ "Positive test (vaccinated)",
                                vacc_at_infection == 0 ~ "Positive test (unvaccinated)"))


# make factors
sccs_data_hes$sex <-as.factor(sccs_data_hes$sex)
sccs_data_hes$age_group <-as.factor(sccs_data_hes$age_group)

# create the interval value, the lenght of the week in days as is sometimes not 7
sccs_data_week_hes <- sccs_data_hes %>%
  group_by(id, week, infected) %>%
  mutate(interval = n()) %>%
  ungroup()

# drop daily columns to make weekly data
sccs_data_week_hes <- sccs_data_week_hes %>%
  select(-day, -infect_day, -death_day) %>%
  distinct() %>%
  rename(death = death_week) %>%
  data.frame()

# cal_day_fortnight variable
sccs_data_week_hes <- sccs_data_week_hes %>%
  mutate(cal_day_fortnight = as.factor(floor(cal_day_week/14)))

#--------------------------
# Filter to not infected on epistart
#--------------------------

# number infected on epistart
n_infect_epistart <- sccs_data_week_hes %>%
  select(id, infect_on_epistart) %>%
  distinct() %>%
  group_by(infect_on_epistart) %>%
  count()
n_infect_epistart

write.csv(n_infect_epistart, paste0(plots_folder, "infect_on_epistart_hes.csv"), row.names = FALSE)

# filter to not infected on day 0
sccs_data_week_hes <- filter(sccs_data_week_hes, infect_on_epistart == 0)


#--------------------------
# Quick checks
#--------------------------

# number of people in dataset 
n_people <- sccs_data_week_hes %>%
  select(id) %>%
  distinct() %>%
  count()
n_people


#--------------------------
# Table 1
#--------------------------
sccs_data_week <- make_risk_period_variables_infections(sccs_data_week, 12)
sccs_data_week_hes <- make_risk_period_variables_infections(sccs_data_week_hes, 12)
sccs_data_week_cardiac <- make_risk_period_variables_infections(sccs_data_week_cardiac, 12)

sccs_data_week <- sccs_data_week %>%
  mutate(risk_period_for_table = case_when(risk_period == "risk" & infect_status == "Positive test (vaccinated)" ~ "12 weeks or less (vaccinated at infection)",
                                          risk_period == "risk" & infect_status == "Positive test (unvaccinated)" ~ "12 weeks or less (unvaccinated at infection)",
                                          risk_period == "baseline" & infect_status == "Positive test (unvaccinated)" & infection_number == 1 ~ "13+ weeks (unvaccinated at infection)",
                                          risk_period == "baseline" & infect_status == "Positive test (vaccinated)" & infection_number == 1 ~ "13+ weeks (vaccinated at infection)",
                                          TRUE ~ "No positive test"))
sccs_data_week_cardiac <- sccs_data_week_cardiac %>%
  mutate(risk_period_for_table = case_when(risk_period == "risk" & infect_status == "Positive test (vaccinated)" ~ "12 weeks or less (vaccinated at infection)",
                                          risk_period == "risk" & infect_status == "Positive test (unvaccinated)" ~ "12 weeks or less (unvaccinated at infection)",
                                          risk_period == "baseline" & infect_status == "Positive test (unvaccinated)" & infection_number == 1 ~ "13+ weeks (unvaccinated at infection)",
                                          risk_period == "baseline" & infect_status == "Positive test (vaccinated)" & infection_number == 1 ~ "13+ weeks (vaccinated at infection)",
                                          TRUE ~ "No positive test"))
sccs_data_week_hes <- sccs_data_week_hes %>%
  mutate(risk_period_for_table = case_when(risk_period == "risk" & infect_status == "Positive test (vaccinated)" ~ "12 weeks or less (vaccinated at infection)",
                                          risk_period == "risk" & infect_status == "Positive test (unvaccinated)" ~ "12 weeks or less (unvaccinated at infection)",
                                          risk_period == "baseline" & infect_status == "Positive test (unvaccinated)" & infection_number == 1 ~ "13+ weeks (unvaccinated at infection)",
                                          risk_period == "baseline" & infect_status == "Positive test (vaccinated)" & infection_number == 1 ~ "13+ weeks (vaccinated at infection)",
                                          TRUE ~ "No positive test"))


sccs_data_week %>% filter(death == 1 & risk_period_for_table=="No positive test") %>% group_by(risk_period, infect_status, infection_number) %>% count()

result_all <- get_cat_vars(filter(sccs_data_week, death==1), c('sex', 'age_group', 
                                                                  'risk_period_for_table', 'infect_status'), total = TRUE)

result_all <-  result_all %>%
  rename('All-cause registered deaths' = Count..n.)
  

result_cardiac <- get_cat_vars(filter(sccs_data_week_cardiac, death==1), c('sex', 
                                              'age_group', 'risk_period_for_table', 'infect_status'), total = TRUE)
result_cardiac <-  result_cardiac %>%
  rename('Cardiac registered deaths' = Count..n.)

result_hes <- get_cat_vars(filter(sccs_data_week_hes, death==1), c('sex', 
                                              'age_group', 'risk_period_for_table', 'infect_status'), total = TRUE)
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
                             variable == 'infect_status' ~ 'Positive SARS-CoV-2 test',
                             variable == 'risk_period_for_table' ~ 'Risk period',
                             TRUE ~ variable),
        breakdown = case_when(breakdown == '1' & variable == 'Sex' ~ 'Male',
                             breakdown == '2' & variable == 'Sex' ~ 'Female',
                             TRUE ~ breakdown))


summary <- summary[,c('variable', 'breakdown', 'All-cause registered deaths', 
                      'Cardiac registered deaths', 'All-cause hospital deaths')]


summary <- summary %>%
  rename(Characteristic = variable,
        Level = breakdown)

summary

write.csv(summary, paste0(plots_folder, "Table1_infections.csv"), row.names = FALSE)




#--------------------------
# Deaths by week
#--------------------------

sccs_data_with_week <- sccs_data %>%
  mutate(week_from_dec_8 = floor(day/7) +1)

sccs_data_with_week_hes <- sccs_data_hes %>%
  mutate(week_from_dec_8 = floor(day/7) +1)


deaths_by_week <- sccs_data_with_week %>%
  group_by(week_from_dec_8, infected) %>%
  summarise(n_deaths = sum(death_day)) %>%
  mutate(cause = 'All-cause deaths')

deaths_by_week_cardiac <- sccs_data_with_week %>%
  filter(cardiac_cause_of_death == 1) %>%
  group_by(week_from_dec_8, infected) %>%
  summarise(n_deaths = sum(death_day)) %>%
  mutate(cause = 'Cardiac-related deaths')

deaths_by_week_hes <- sccs_data_with_week_hes %>%
  group_by(week_from_dec_8, infected) %>%
  summarise(n_deaths = sum(death_day)) %>%
  mutate(cause = 'All-cause deaths (HES)')

all_deaths_by_week <- rbind(deaths_by_week, deaths_by_week_cardiac)
all_deaths_by_week <- rbind(all_deaths_by_week, deaths_by_week_hes)

write.csv(all_deaths_by_week, paste0(plots_folder, "SupFig4_all_weeks_infections.csv"), row.names=FALSE)

all_deaths_by_week$infected <- as.character(all_deaths_by_week$infected)

all_deaths_by_week <- all_deaths_by_week %>%
  mutate(infected = ifelse(cause=="Cardiac-related deaths", "Suppressed", infected)) %>%
  group_by(week_from_dec_8 , infected, cause) %>%
  mutate(n_deaths = sum(n_deaths)) %>%
  ungroup() %>%
  distinct()


write.csv(all_deaths_by_week, paste0(plots_folder, "SupFig4_all_weeks_infections_suppressed.csv"), row.names=FALSE)

all_deaths_by_week$cause <- factor(all_deaths_by_week$cause,
                     levels = c("All-cause deaths", "Cardiac-related deaths", "All-cause deaths (HES)"))
levels(all_deaths_by_week$cause) <- c("All-cause registered deaths", "Cardiac registered deaths", 
                                         "All-cause hospital deaths")

all_deaths_by_week$infected <- factor(all_deaths_by_week$infected,
                     levels = c("Suppressed", "Infected", "Not infected"))
levels(all_deaths_by_week$infected) <- c("Suppressed", "Positive SARS-CoV-2 test", "No positive SARS-CoV-2 test")

plt <- ggplot(all_deaths_by_week, aes(x=week_from_dec_8, y=n_deaths, fill = infected)) + 
  geom_bar(stat='identity', position = "stack") +
  theme_light() +
  facet_wrap(~cause) +
  xlab("Week since 8 December 2020") +
  ylab("Count of deaths") +
  scale_fill_discrete() +
  theme(text=element_text(size=20),
       legend.title = element_blank()) +
  guides(fill = guide_legend(reverse=TRUE)) +
  xlim(0, 70)
plt

ggsave(paste0(plots_folder, "SupFig4_all_weeks_infections.png"), width = 15, height = 5)



#--------------------------
# Deaths since dose
#--------------------------


# events vs. week number
events_per_week_all_dose <- sccs_data_week %>%
  filter(infected == "Infected") %>%
  group_by(week) %>%
  summarise(n_deaths = sum(death), n_events_per_week = (sum(death)/n())*100000) %>%
  mutate(cause = "All-cause deaths")


events_per_week_cardiac_dose <- sccs_data_week_cardiac %>%
  filter(infected == "Infected") %>%
  group_by(week) %>%
  summarise(n_deaths = sum(death), n_events_per_week = (sum(death)/n())*100000) %>%
  mutate(cause = "Cardiac-related deaths")


events_per_week_hes_dose <- sccs_data_week_hes %>%
  filter(infected == "Infected") %>%
  group_by(week) %>%
  summarise(n_deaths = sum(death), n_events_per_week = (sum(death)/n())*100000) %>%
  mutate(cause = "All-cause deaths (HES)")

events_per_week <- rbind(events_per_week_all_dose, events_per_week_cardiac_dose, events_per_week_hes_dose)

write.csv(events_per_week, paste0(plots_folder, "SupFig1b_deaths_since_dose_infections.csv"), row.names=FALSE)

events_per_week$cause <- factor(events_per_week$cause,
                     levels = c("All-cause deaths", "Cardiac-related deaths", "All-cause deaths (HES)"))
levels(events_per_week$cause) <- c("All-cause registered deaths", "Cardiac registered deaths", 
                                         "All-cause hospital deaths")

plt <- ggplot(data = events_per_week, aes(x=week, y=n_deaths)) +
  geom_bar(stat = 'identity') + facet_grid(~cause) +
  theme_light() +
  xlab("Week since infection") +
  ylab("Number of deaths") +
  theme(text=element_text(size=16)) +
  xlim(0, 70)
plt

ggsave(paste0(plots_folder, "SupFig1b_deaths_since_dose_infections.png"), width = 8, height = 3)


#--------------------------
# Create sample to use
#--------------------------

# filter to vaccinated/not vaccinated if required
if (infection_cat == "vaccinated_at_infection"){
  #filter to exclude people unvaccinated at infection
  sccs_data_week_hes <- sccs_data_week_hes %>%
    filter((vacc_at_infection == 1 & had_infection == 1) | had_infection == 0)
  sccs_data_week <- sccs_data_week %>%
    filter((vacc_at_infection == 1 & had_infection == 1) | had_infection == 0)
}else if (infection_cat == "unvaccinated_at_infection"){
  # filter to exclude people vaccinated at infection
  sccs_data_week_hes <- sccs_data_week_hes %>%
    filter(vacc_at_infection == 0)
  sccs_data_week <- sccs_data_week %>%
    filter(vacc_at_infection == 0)
}

sccs_data_week_cardiac <- sccs_data_week %>%
  filter(cardiac_cause_of_death == 1)



#--------------------------
# Create dataset
#--------------------------


# data for all (to estimate seasonal effects)
create_dataset_non_repeating <-function(df_week, risk_period_length) {
  
# create variables to indicate whether in risk period or not, and the week in the risk period
  df_week <- df_week %>%
    mutate(risk_period = case_when(week <= risk_period_length & infection_number != 0 ~ 'risk', 
                                   week > risk_period_length & infection_number != 0 ~ 'baseline',
                                  infection_number == 0 ~ 'baseline')) %>%
    mutate(risk_week = case_when(week <= risk_period_length & infection_number != 0 ~ paste0("week_", str_pad(week, 2, pad="0")), 
                                   week > risk_period_length & infection_number != 0 ~ 'baseline',
                                  infection_number == 0 ~ 'baseline'))
  

# create weights for the period and week by doses analysis, with a different weight for each period-dose / week-dose
  df_all  <- df_week %>% 
    mutate(eventweight_period_combined = ifelse(death == 1 & risk_period == "risk", infection_number, 0),
           eventweight_period_combined = ifelse(death == 1 & risk_period == "risk", week, 0),
          risk_period = "baseline",
          risk_week = "baseline",
          id2 = id) 
  
# data after dose 3
# take into account all doses
# set weights to zero
df_d1 <-  df_week %>% 
  filter(infection_number == 1)%>%
  mutate(eventweight_period_combined = 0 ,
         eventweight_period_combined = 0,
         id2 = id + max(df_all$id)*3)


  # combine 
  df_clean <-rbind(df_all, df_d1)
  


  df_clean$risk_period <- relevel(as.factor(df_clean$risk_period), ref="baseline")
  df_clean$risk_week <- relevel(as.factor(df_clean$risk_week), ref="baseline")

  
  return(df_clean)
}





df_clean <- create_dataset_non_repeating(sccs_data_week, 12)
df_clean_hes <- create_dataset_non_repeating(sccs_data_week_hes, 12)
df_clean_cardiac <- create_dataset_non_repeating(sccs_data_week_cardiac, 12)


#--------------------------
# Model
#--------------------------


run_model <- function(dataset, suffix, calendar_adjustment){

  res_week <- eventdepexp_model(outcome = "death",
                 exposure = "risk_week",
                 x= calendar_adjustment,
                 data = dataset,
                 model =F,
                 save_name = paste0("Risk_week_results_infections", suffix))
  res_week$period = "week"

  res_period <- eventdepexp_model(outcome = "death",
                 exposure = "risk_period",
                 x= calendar_adjustment,
                 data = dataset,
                 model =F,
                 save_name = paste0("Risk_period_results_infections", suffix))
  res_period$period = "period"
  
  all_data <- rbind(res_week, res_period)
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

write.csv(all_res_parameters, paste0(plots_folder, "/Main_results_risk_period_and_all_weeks_spline_cal_infections.csv"),
         row.names = FALSE)



# fortnight cal adjustment
res_reg <- run_model(df_clean, "_fortnight_cal_reg", "cal_day_fortnight")
res_reg$cause <- "All-cause registered deaths"

res_hes <- run_model(df_clean_hes, "_spline_cal_hes", "cal_day_fortnight")
res_hes$cause <- "All-cause hospital deaths"

res_cardiac <- run_model(df_clean_cardiac, "_spline_cal_cardiac", "cal_day_fortnight")
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

write.csv(all_res_parameters, paste0(plots_folder, "/Main_results_risk_period_and_all_weeks_fortnight_cal_infections.csv"),
         row.names = FALSE)



#--------------------------
# Plot IRR
#--------------------------


if (infection_cat == "vaccinated_at_infection"){
  width <- 8
  }else{width<-12}


type <- "_spline_cal"
plot_data <- read.csv(paste0(plots_folder, "Main_results_risk_period_and_all_weeks_spline_cal_infections.csv"))


plot_data <- plot_data[plot_data$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))

plot_data$parameter <- factor(plot_data$parameter,
                     levels = c("Week 12", "Week 11", "Week 10", "Week 9", "Week 8", "Week 7", 
                               "Week 6", "Week 5", "Week 4", "Week 3", "Week 2", "Week 1", "Weeks 1-12"))

if (infection_cat == "vaccinated_at_infection"){
  plot_data <- filter(plot_data, cause != "Cardiac registered deaths")
}

write.csv(plot_data, paste0(plots_folder, "Fig3_infections", type, ".csv"),
         row.names = FALSE)

plt <- ggplot(plot_data, aes(x=irr, y=parameter)) + 
    geom_point(aes(colour = period), size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour=period), width=.2) +
    facet_wrap(~cause) +
    theme_light() +
    xlab("Relative incidence") +
    scale_x_continuous(trans='log2', breaks = c(1/64, 1/8, 1, 8, 64), 
                       labels = c("1/64", "1/8", "1", "8", "64"), limits = c(1/128, 128)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
        axis.title.y=element_blank(),
       legend.position = 'none',
       panel.grid.major.y = element_blank())

  print(plt)


ggsave(paste0(plots_folder, "Fig3_infections", type, ".png"), width = width, height = 7)



#--------------------------
# Different risk periods
#--------------------------


get_all_risk_data <- function(risk_period_length){

  df_clean <- create_dataset_non_repeating(sccs_data_week, risk_period_length)
  # all_deaths
  all <- eventdepexp_model(outcome = "death",
                   exposure = "risk_period",
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean,
                   model =F)

  all$cause <- "All-cause deaths"
  

  
  if (infection_cat != "vaccinated_at_infection"){
    df_clean_cardiac <- create_dataset_non_repeating(sccs_data_week_cardiac, risk_period_length)
    # cardiac related
    cardiac <- eventdepexp_model(outcome = "death",
                   exposure = "risk_period",
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean_cardiac,
                   model =F)

    cardiac$cause <- "Cardiac-related deaths"
  }
  

  
  df_clean_hes <- create_dataset_non_repeating(sccs_data_week_hes, risk_period_length)

  hes <- eventdepexp_model(outcome = "death",
                   exposure = "risk_period",
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean_hes,
                   model =F)

  hes$cause <- "All-cause deaths (HES)"

  if (infection_cat == "vaccinated_at_infection"){
    all_data <- rbind(all, hes)
  }else{
  all_data <- rbind(all, cardiac)
  all_data <- rbind(all_data, hes)
  }
  
  all_data$risk_period_length <- risk_period_length
  
  return(all_data)
}

risk_period_length <- 2
all_data <- get_all_risk_data(risk_period_length)

risk_period_lengths = 3:40
for (i in risk_period_lengths){

  combined <- get_all_risk_data(i)

  all_data <- rbind(all_data, combined)
}

all_data


write.csv(all_data, 
          paste0(plots_folder, "SupFig6_different_risk_period_length_infections.csv"),
         row.names = FALSE)

all_data <- read.csv(paste0(plots_folder, "SupFig6_different_risk_period_length_infections.csv"))

all_data <- all_data %>%
  filter(term == 'risk_periodrisk')



all_data$cause <- factor(all_data$cause,
                     levels = c("All-cause deaths", "Cardiac-related deaths", "All-cause deaths (HES)"))
levels(all_data$cause) <- c("All-cause registered deaths", "Cardiac registered deaths", 
                                         "All-cause hospital deaths")


if (infection_cat == "vaccinated_at_infection"){
  width <- 8
  }else{width<-12}


all_data_12 <- filter(all_data, risk_period_length <=12)

plt <- ggplot(all_data_12, aes(x=irr, y=risk_period_length)) + 
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
    facet_wrap(~cause) +
    theme_light() +
    xlab("Relative incidence") +
    ylab("Risk period length (weeks)") +
    scale_x_continuous(trans='log2', breaks = c(1/16, 1/4,  1,  4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(1/32, 32)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.7, vjust = -1)+
    scale_y_continuous(limits = c(1.5,12.5), breaks = c(2,4,6,8,10,12), 
                     labels = c(2,4,6,8,10,12)) +
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
       legend.position = 'none',
       panel.grid.major.y = element_blank(),
       panel.spacing.x = unit(1.5, "lines"),
       plot.margin = margin(10, 10, 10, 10))
#  xlim(0.125, 2)

plt

ggsave(paste0(plots_folder, "SupFig6_different_risk_period_length_infections.png"), height=6, width=width)


all_data_long <- all_data

plt <- ggplot(all_data_long, aes(x=irr, y=risk_period_length)) + 
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
    facet_grid(~cause) +
    theme_light() +
    xlab("Relative incidence") +
    ylab("Risk period length (weeks)") +
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
    scale_x_continuous(trans='log2', breaks = c(1/16, 1/4, 1, 4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(1/32, 32)) +
  theme(text=element_text(size=20),
       legend.position = 'none',
       panel.grid.major.y = element_blank(),
       panel.spacing.x = unit(1.5, "lines"),
       plot.margin = margin(10, 10, 10, 10))
#  xlim(0.125, 2)

plt

ggsave(paste0(plots_folder, "SupFig6b_different_risk_period_length_long_all.png"), height=8, width=width)

all_data_long_24 <- filter(all_data, risk_period_length <=24)
all_data_long_24 <- filter(all_data_long_24, irr <=100)
all_data_long_24 <- all_data_long_24 %>%
  mutate(is_12 = as.factor(ifelse(risk_period_length != 12, 0, 1)))

if (infection_cat == "unvaccinated_at_infection"){
  breaks <- c(1/16, 1/4, 1,4,16)
  labels <- c("1/16", "1/4", "1","4","16")
  limits <- c(1/32, 32)
  }else{breaks = c(1/64, 1/8, 1, 8, 64)
       labels = c("1/64", "1/8", "1", "8", "64")
       limits = c(1/128, 128)}

plt <- ggplot(all_data_long_24, aes(x=irr, y=risk_period_length, color = is_12)) + 
    geom_point(size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci), width = .2) +
    facet_grid(~cause) +
    theme_light() +
    xlab("Relative incidence") +
    ylab("Risk period length (weeks)") +
    scale_x_continuous(trans='log2', breaks = breaks, 
                       labels = labels, limits = limits) +
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

ggsave(paste0(plots_folder, "SupFig6b_different_risk_period_length_long_all_24.png"), height=8, width=width)



#--------------------------
# Age sex breakdown
#--------------------------

df_clean <- create_dataset_non_repeating(sccs_data_week, 12)
df_clean_hes <- create_dataset_non_repeating(sccs_data_week_hes, 12)
df_clean_cardiac <- create_dataset_non_repeating(sccs_data_week_cardiac, 7)

exposure <- "risk_period"

sex_parameters_all <- data_for_risk_period_breakdown(df_clean, 'sex', 'All-cause deaths', exposure)
sex_parameters_hes <- data_for_risk_period_breakdown(df_clean_hes, 'sex', 'All-cause deaths (HES)', exposure)
if (infection_cat != "vaccinated_at_infection"){
  sex_parameters_cardiac <- data_for_risk_period_breakdown(df_clean_cardiac, 'sex', 'Cardiac-related deaths', exposure)
  sex_parameters <- rbind(sex_parameters_all, sex_parameters_cardiac, sex_parameters_hes)
}else{
  sex_parameters <- rbind(sex_parameters_all, sex_parameters_hes)
}

age_group_parameters_all <- data_for_risk_period_breakdown(df_clean, 'age_group', 'All-cause deaths', exposure)
age_group_parameters_hes <- data_for_risk_period_breakdown(df_clean_hes, 'age_group', 'All-cause deaths (HES)', exposure)

if (infection_cat != "vaccinated_at_infection"){
  age_group_parameters_cardiac <- data_for_risk_period_breakdown(df_clean_cardiac, 'age_group', 'Cardiac-related deaths', exposure)
  age_group_parameters <- rbind(age_group_parameters_all, age_group_parameters_cardiac, age_group_parameters_hes)
}else{
  age_group_parameters <- rbind(age_group_parameters_all, age_group_parameters_hes)
}


# combine and save breakdown data
all_breakdowns <- rbind(sex_parameters, age_group_parameters)


write.csv(all_breakdowns, 
          paste0(results_folder, "breakdown_results_together_infections.csv"),
         row.names = FALSE)

all_breakdowns$breakdown_value <- recode_factor(all_breakdowns$breakdown_value, '1' = "Male", "2" = "Female")

all_breakdowns <- all_breakdowns[!grepl("cal_day_week", all_breakdowns$term), ]


all_breakdowns$cause <- factor(all_breakdowns$cause,
                     levels = c("All-cause deaths", "Cardiac-related deaths", "All-cause deaths (HES)"))

all_breakdowns$dose <- factor(all_breakdowns$dose,
                     levels = c("First", "Second", "Third", "All doses"))


write.csv(all_breakdowns, paste0(plots_folder, "Fig4_risk_period_all_breakdowns_infections.csv"), row.names = FALSE)


all_breakdowns <- all_breakdowns[all_breakdowns$irr_low_ci >= 1E-12, ]


all_breakdowns$cause <- factor(all_breakdowns$cause,
                     levels = c("All-cause deaths", "Cardiac-related deaths", "All-cause deaths (HES)"))
levels(all_breakdowns$cause) <- c("All-cause registered deaths", "Cardiac registered deaths", 
                                         "All-cause hospital deaths")


if (infection_cat == "vaccinated_at_infection"){
  width <- 8
  }else{width<-12}


plt <- ggplot(all_breakdowns, aes(x=irr, y=breakdown_value)) + 
geom_point(aes(colour = breakdown), size=3) +
geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour = breakdown), width = .2) +
facet_grid(~cause) +
theme_light() +
xlab("Relative incidence") +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
geom_vline(xintercept=1, linetype='dashed', colour='grey') +
    scale_x_continuous(trans='log2', breaks = c(1/16, 1/4, 1, 4, 16), 
                       labels = c("1/16", "1/4", "1", "4", "16"), limits = c(1/64, 64)) +
theme(text=element_text(size=20),
        axis.title.y=element_blank(),
 #      strip.text.y = element_blank(),
     legend.position = 'none',
     panel.grid.major.y = element_blank())
  
plt



ggsave(paste0(plots_folder, "Fig4_risk_period_all_breakdowns_infections.png"), height=4, width=width)





#--------------------------
# Fortnightly calendar adjustment
#--------------------------



type <- "_fortnight_cal"
plot_data <- read.csv(paste0(plots_folder, "Main_results_risk_period_and_all_weeks_fortnight_cal_infections.csv"))


plot_data <- plot_data[plot_data$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))

plot_data$parameter <- factor(plot_data$parameter,
                     levels = c("Week 12", "Week 11", "Week 10", "Week 9", "Week 8", "Week 7", 
                               "Week 6", "Week 5", "Week 4", "Week 3", "Week 2", "Week 1", "Weeks 1-12"))

if (infection_cat == "vaccinated_at_infection"){
  plot_data <- filter(plot_data, cause != "Cardiac registered deaths")
}

write.csv(plot_data, paste0(plots_folder, "Fig3_infections", type, ".csv"),
         row.names = FALSE)

plt <- ggplot(plot_data, aes(x=irr, y=parameter)) + 
    geom_point(aes(colour = period), size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour=period), width=.2) +
    facet_wrap(~cause) +
    theme_light() +
    xlab("Relative incidence") +
    scale_x_continuous(trans='log2', breaks = c(1/64, 1/8, 1, 8, 64), 
                       labels = c("1/64", "1/8", "1", "8", "64"), limits = c(1/128, 128)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
        axis.title.y=element_blank(),
       legend.position = 'none',
       panel.grid.major.y = element_blank())

  print(plt)


ggsave(paste0(plots_folder, "Fig3_infections", type, ".png"), width = width, height = 7)


#--------------------------
# Include epistart = infect day
#--------------------------


# create the interval value, the lenght of the week in days as is sometimes not 7
sccs_data_week_hes <- sccs_data_hes %>%
  group_by(id, week, infected) %>%
  mutate(interval = n()) %>%
  ungroup()

sccs_data_week_hes <- make_risk_period_variables(sccs_data_week_hes, 12)

# drop daily columns to make weekly data
sccs_data_week_hes <- sccs_data_week_hes %>%
  select(-day, -infect_day, -death_day) %>%
  distinct() %>%
  rename(death = death_week) %>%
  data.frame()

# cal_day_fortnight variable
sccs_data_week_hes <- sccs_data_week_hes %>%
  mutate(cal_day_fortnight = as.factor(floor(cal_day_week/14)))

# number of people in dataset 
n_people <- sccs_data_week_hes %>%
  select(id) %>%
  distinct() %>%
  count()
n_people

# number infected on epistart
n_infect_epistart <- sccs_data_week_hes %>%
  select(id, infect_on_epistart) %>%
  distinct() %>%
  group_by(infect_on_epistart) %>%
  count()
n_infect_epistart

# filter to vaccinated/not vaccinated if required
if (infection_cat == "vaccinated_at_infection"){
  #filter to exclude people unvaccinated at infection
  sccs_data_week_hes <- sccs_data_week_hes %>%
    filter((vacc_at_infection == 1 & had_infection == 1) | had_infection == 0)
  sccs_data_week_hes <- sccs_data_week_hes %>%
    filter((vacc_at_infection == 1 & had_infection == 1) | had_infection == 0)
}else if (infection_cat == "unvaccinated_at_infection"){
  # filter to exclude people vaccinated at infection
  sccs_data_week_hes <- sccs_data_week_hes %>%
    filter(vacc_at_infection == 0)
  sccs_data_week_hes <- sccs_data_week_hes %>%
    filter(vacc_at_infection == 0)
}


df_clean_hes <- create_dataset_non_repeating(sccs_data_week_hes, 12)

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

write.csv(all_res_parameters, paste0(plots_folder, "/Main_results_risk_period_and_all_weeks_spline_cal_infections_include_infect_on_epistart.csv"),
         row.names = FALSE)


plot_data <- read.csv(paste0(plots_folder, "/Main_results_risk_period_and_all_weeks_spline_cal_infections_include_infect_on_epistart.csv"))

plot_data <- plot_data[plot_data$irr_low_ci >= 1E-12, ]

plot_data$cause <- factor(plot_data$cause,
                     levels = c("All-cause registered deaths", "Cardiac registered deaths", 
                                "All-cause hospital deaths"))

plot_data$parameter <- factor(plot_data$parameter,
                     levels = c("Week 12", "Week 11", "Week 10", "Week 9", "Week 8", "Week 7", 
                               "Week 6", "Week 5", "Week 4", "Week 3", "Week 2", "Week 1", "Weeks 1-12"))


write.csv(plot_data, paste0(plots_folder, "Fig3_infections_include_infect_on_epistart", ".csv"),
         row.names = FALSE)

plt <- ggplot(plot_data, aes(x=irr, y=parameter)) + 
    geom_point(aes(colour = period), size=3) +
    geom_errorbar(aes(xmin = irr_low_ci, xmax=irr_up_ci, colour=period), width=.2) +
    theme_light() +
    xlab("Relative incidence") +
    scale_x_continuous(trans='log2', breaks = c(1/64, 1/8, 1, 8, 64), 
                       labels = c("1/64", "1/8", "1", "8", "64"), limits = c(1/128, 128)) +
    geom_text(aes(label = paste0(sprintf("%.2f", round(irr, 2)), " (", sprintf("%.2f", round(irr_low_ci,2)),
                                 " - ", sprintf("%.2f", round(irr_up_ci,2)), ")")),  hjust = 0.5, vjust = -1)+
    geom_vline(xintercept=1, linetype='dashed', colour='grey') +
  theme(text=element_text(size=20),
        axis.title.y=element_blank(),
       legend.position = 'none',
       panel.grid.major.y = element_blank())

  print(plt)


ggsave(paste0(plots_folder, "Fig3_infections_include_infect_epistart.png"), width = 4, height = 7)