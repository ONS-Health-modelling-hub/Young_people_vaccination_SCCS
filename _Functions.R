## descriptive tables

get_cat_vars <- function(df, vars, total = FALSE){
  
  res <- list()
  
  if (total==TRUE){

      df_summary <- df %>%
        tally() %>%
        collect() %>%
        as.data.frame()%>%
        arrange()
    

      df_summary$breakdown <- rep(NA, nrow(df_summary))
      df_summary$n <- format(df_summary$n, big.mark=",", scientific=FALSE)
      df_summary[['Count (n)']] <- paste0(df_summary$n)
      df_summary <- select(df_summary, -n)
      df_summary$variable <- rep('Total', nrow(df_summary))


    res <- append(res, list(df_summary))

    rm(df_summary)
    gc()
  }

  for (v in vars){
    
    v <- as.symbol(v)
    v_str <- as.character(v)
    
    df_var <- df %>%
      select(v)
    
      df_summary <- df_var %>%
        group_by(!!v, .drop=FALSE) %>%
        tally() %>%
        collect() %>%
        as.data.frame()%>%
        arrange()


      df_summary$percentage <- df_summary$n / sum(df_summary$n) * 100
  
      df_summary$percentage <- round(df_summary$percentage, 2)

      df_summary$n <- format(df_summary$n, big.mark=",", scientific=FALSE)

      df_summary[[v_str]] <- as.character(df_summary[[v_str]])
      

      df_summary[['Count (n)']] <- paste0(df_summary$n, ' (', format(df_summary$percentage, nsmall = 2), ')')

      df_summary <- select(df_summary, -n, -percentage)
  
  
      df_summary$variable <- rep(v_str, nrow(df_summary))
    
      df_summary <- df_summary %>%
        rename(breakdown = v_str)
    
    res <- append(res, list(df_summary))
    
    rm(df_var, df_summary)
    gc()
    
  }
  result <- data.frame(rbindlist(res))
  return(result)
}


# function to format the tables
format_table<- function(data, name='Count (n)'){
  ind <- c(sapply(continuous_vars, function(x) which(columns == x)[[1]]))
  cat <- data[-ind]
  con <- data[ind]


  for (i in 1:length(cat)){
    print(colnames(cat[[i]]))
    colnames(cat[[i]]) <- c('value', name, 'variable')
    cat[[i]]<- arrange(cat[[i]],value)
  }

  cat <- do.call('rbind', cat) 
  cat <- cat %>%
    left_join(var_lookup)%>%
    group_by(variable)%>%
    filter(is.na(type)|type == ""|type=="binary"&as.numeric(value) ==1)%>%
    as.data.frame()
  
    con <- do.call('rbind', con)%>%
  mutate(value="", labels="", type="")
names(con) <-c(name, 'value', "variable" , "labels" ,   "type" )

rbind(cat, con)
  
}
                  

# sccs model                  
eventdepexp_model <- function(outcome, exposure, x, data, tolerance = 0.000001, display=T, 
                              model = F, save_name = NULL, daily = FALSE){

  
  # number of vaccine doses (3 for by dose analysis, 1 if combined) * number of risk periods (1 for by period, 6 for by week)
  k = nlevels(as.factor(data[[exposure]])) - 1
  print(k)

  beta = rep( 0, k )
  betadiff = 1

  # initialise the event number column to be 1 where an event occurs in that week, 0 otherwise
  data$nevents = data[[outcome]] 
  d <- data

  if (exposure == "risk_week" | exposure == "risk_week_vector"){
       d$eventweight = d$eventweight_week_combined
  }else if(exposure == "risk_week_dose" | exposure == "risk_week_dose_vector"){
       d$eventweight = d$eventweight_week
  }else if(exposure == "risk_period" | exposure == "risk_period_vector"){
       d$eventweight = d$eventweight_period_combined
  }else if(exposure == "risk_period_dose" | exposure == "risk_period_dose_vector"){
       d$eventweight = d$eventweight_period
     }
  
#  d$eventweight <- as.factor(d$eventweight)

  # counter for runs through the loop
  i = 0

  while(betadiff > tolerance){

    beta_in = beta 

    # update the event number column:
    # if the eventweight is not 0, get the value of beta at that eventweight position in the vector
    # and set the event numebr column to be 1/exp(this value of beta). Otherwise, leave it as it was.
    d$nevents <- ifelse(d$eventweight != 0 ,  1/exp(beta[d$eventweight]), d$nevents)
  

    X = paste(exposure, x, sep="+")

    formula = as.formula(paste0("nevents ~ ", X))

    if (daily == TRUE){
    m <- gnm(formula, eliminate=as.factor(id2), 
             data = d ,
             family = poisson)
    }else{
     m <- gnm(formula, eliminate=as.factor(id2), 
             data = d ,
             offset = log(interval), 
             family = poisson)
    }


    # set new beta to be the coefficients from the model
    beta <- m$coefficients[1:k]

    # calculate the sum of the absolute differences between the old and new betas and stop if this is less than the tolerance
    betadiff = sum(abs(beta_in - beta))

    i = i + 1
    if (display == TRUE){
    print(paste(i, betadiff))
    }
  }
  
  if(model == TRUE){
  return(m)
  }else{
    b = coef(m)
    se <-  sqrt(diag(vcov(m)))
    
    res <- data.frame(term = names(b),
                      irr = exp(b),
                      irr_low_ci = exp(b - 1.97 * se),
                      irr_up_ci = exp(b + 1.97 * se),
                      p_value = summary(m)$coefficients[,4])
    
    if (!is.null(save_name)){
      write.csv(res, paste0(results_folder, save_name, ".csv"), row.names = FALSE)
    }
    return(res)
  }
  
}

#--------------------------
# Prepare data for model
#--------------------------


# data for all (to estimate seasonal effects)
create_dataset <-function(df_week, risk_period_length) {
  
  df_week <- make_risk_period_variables(df_week, risk_period_length)


# create weights for the period and week by doses analysis, with a different weight for each period-dose / week-dose
  df_all  <- df_week %>% 
    mutate(eventweight_period = ifelse(death == 1 & risk_period == "risk", current_dose, 0),
           eventweight_week = ifelse(death == 1 & risk_period == "risk", week + (current_dose-1)*risk_period_length, 0),
          risk_period = "baseline",
          risk_period_dose = "baseline",
          risk_week_dose = "baseline",
          risk_week = "baseline",
          id2 = id) 
  
# data after dose 3
# take into account all doses
# set weights to zero
df_d3 <-  df_week %>% 
  filter(current_dose == 3)%>%
  mutate(eventweight_period = 0 ,
         eventweight_week = 0,
         id2 = id + max(df_all$id)*3)

# data after dose 2
# take into account doses 1 and 2, setting data for dose 3 to baseline
# if the current dose is 3,set the weights in the week a death occurred
df_d2 <-  df_week %>% 
  filter(current_dose >= 2)%>%
  mutate(eventweight_period = ifelse(death == 1 & risk_period == "risk" & current_dose == 3, current_dose, 0),
         eventweight_week = ifelse(death == 1 & risk_period == "risk" & current_dose == 3,  week + (current_dose-1)*risk_period_length , 0), 
         risk_period_dose = ifelse(risk_period == "risk" & current_dose == 3, "baseline", risk_period_dose),
         risk_week = ifelse(risk_period == "risk" & current_dose == 3, "baseline", risk_week),
         risk_week_dose = ifelse(risk_period == "risk" & current_dose == 3, "baseline", risk_week_dose),
         risk_period = ifelse(risk_period == "risk" & current_dose == 3, "baseline", risk_period),
         id2 = id + max(df_all$id)*2)
  

# data after dose 1
# take into account dose 1, setting data for dose 2 and 3 to baseline
# if the current dose is 2 or 3,set the weights in the week a death occurred
df_d1 <-  df_week %>% 
  filter(current_dose >= 1)%>%
  mutate(eventweight_period = ifelse(death == 1 & risk_period == "risk" & current_dose >= 2, current_dose, 0),
         eventweight_week = ifelse(death == 1 & risk_period == "risk" & current_dose >= 2,  week + (current_dose-1)*risk_period_length , 0),
         risk_period_dose = ifelse(risk_period == "risk" & current_dose >= 2, "baseline", risk_period_dose),
         risk_week = ifelse(risk_period == "risk" & current_dose >= 2, "baseline", risk_week),
         risk_week_dose = ifelse(risk_period == "risk" & current_dose >= 2, "baseline", risk_week_dose),
         risk_period = ifelse(risk_period == "risk" & current_dose >= 2, "baseline", risk_period),
         id2 = id + max(df_all$id))

  # combine 
  df_clean <-rbind(df_all, df_d1, df_d2, df_d3)
  
# chaneg the weights when combining doses to be just dependent on the week, not the dose
  df_clean <-df_clean %>%
      mutate(eventweight_period_combined = ifelse(eventweight_period > 0, 1,0),
            eventweight_week_combined = ifelse(eventweight_week > 0, week, 0))


  df_clean$risk_period_dose <- relevel(as.factor(df_clean$risk_period_dose), ref="baseline")
  df_clean$risk_week_dose <- relevel(as.factor(df_clean$risk_week_dose), ref="baseline")
  df_clean$risk_period <- relevel(as.factor(df_clean$risk_period), ref="baseline")
  df_clean$risk_week <- relevel(as.factor(df_clean$risk_week), ref="baseline")
 
  return(df_clean)
}
    
    

# data for all (to estimate seasonal effects)
create_dataset_vacc_vecc <-function(df_week, risk_period_length, vector) {
  
  df_week <- make_risk_period_variables(df_week, risk_period_length)
 
  df_week <- df_week %>%
   mutate(risk_week_dose_vector = ifelse(risk_week_dose != "baseline", paste0(risk_week_dose, "_", vaccine_vector), "baseline"),
          risk_week_vector = ifelse(risk_week != "baseline", paste0(risk_week, "_", vaccine_vector), "baseline"),
              risk_period_dose_vector = ifelse(risk_period_dose != "baseline", paste0(risk_period_dose, "_", vaccine_vector), "baseline"),
          risk_period_vector = ifelse(risk_period != "baseline", paste0(risk_period, "_", vaccine_vector), "baseline"))
      
    # create weights for the period and week by doses analysis, with a different weight for each period-dose / week-dose
    df_all  <- df_week %>% 
      mutate(eventweight_period = ifelse(death == 1 & risk_period == "risk", current_dose, 0),
           eventweight_week = ifelse(death == 1 & risk_period == "risk", week + (current_dose-1)*risk_period_length, 0))

    
    df_all  <- df_all %>% 
      mutate(eventweight_period = case_when(eventweight_period== 0 ~ 0,
                                            eventweight_period>0 & vaccine_vector == "mRNA" ~ eventweight_period*2 - 1,
                                            eventweight_period>0 & vaccine_vector == "not_mRNA_or_unknown" ~ eventweight_period*2),
             eventweight_week = case_when(eventweight_week== 0 ~ 0,
                                            eventweight_week>0 & vaccine_vector == "mRNA" ~ eventweight_week*2 - 1,
                                            eventweight_week>0 & vaccine_vector == "not_mRNA_or_unknown" ~ eventweight_week*2))
 
    df_all  <- df_all %>%        
       mutate(risk_period = "baseline",
        risk_period_dose = "baseline",
        risk_week_dose = "baseline",
        risk_week = "baseline",
        risk_week_dose_vector = "baseline",
        risk_week_vector = "baseline",
        risk_period_vector = "baseline",
        risk_period_dose_vector = "baseline",
        id2 = id)

  
# data after dose 3
# take into account all doses
# set weights to zero
df_d3 <-  df_week %>% 
  filter(current_dose == 3)%>%
  mutate(eventweight_period = 0 ,
         eventweight_week = 0,
         id2 = id + max(df_all$id)*3)

# data after dose 2
# take into account doses 1 and 2, setting data for dose 3 to baseline
# if the current dose is 3,set the weights in the week a death occurred
df_d2 <-  df_week %>% 
  filter(current_dose >= 2)%>%
  mutate(eventweight_period = ifelse(death == 1 & risk_period == "risk" & current_dose == 3, current_dose, 0),
         eventweight_week = ifelse(death == 1 & risk_period == "risk" & current_dose == 3,  week + (current_dose-1)*risk_period_length , 0)) %>%
      mutate(eventweight_period = case_when(eventweight_period== 0 ~ 0,
                                            eventweight_period>0 & vaccine_vector == "mRNA" ~ eventweight_period*2 - 1,
                                            eventweight_period>0 & vaccine_vector == "not_mRNA_or_unknown" ~ eventweight_period*2),
             eventweight_week = case_when(eventweight_week== 0 ~ 0,
                                            eventweight_week>0 & vaccine_vector == "mRNA" ~ eventweight_week*2 - 1,
                                            eventweight_week>0 & vaccine_vector == "not_mRNA_or_unknown" ~ eventweight_week*2)) %>%
  mutate(risk_period_dose = ifelse(risk_period == "risk" & current_dose == 3, "baseline", risk_period_dose),
         risk_week = ifelse(risk_period == "risk" & current_dose == 3, "baseline", risk_week),
         risk_week_dose = ifelse(risk_period == "risk" & current_dose == 3, "baseline", risk_week_dose),
         risk_period = ifelse(risk_period == "risk" & current_dose == 3, "baseline", risk_period),
         risk_week_dose_vector = ifelse(risk_week_dose_vector == "risk" & current_dose == 3, "baseline", risk_week_dose_vector),
         risk_period_dose_vector = ifelse(risk_period_dose_vector == "risk" & current_dose == 3, "baseline", risk_period_dose_vector),
         risk_week_vector = ifelse(risk_week_vector == "risk" & current_dose == 3, "baseline", risk_week_vector),
         risk_period_vector = ifelse(risk_period_vector == "risk" & current_dose == 3, "baseline", risk_period_vector),
         id2 = id + max(df_all$id)*2)
  

# data after dose 1
# take into account dose 1, setting data for dose 2 and 3 to baseline
# if the current dose is 2 or 3,set the weights in the week a death occurred
df_d1 <-  df_week %>% 
  filter(current_dose >= 1)%>%
  mutate(eventweight_period = ifelse(death == 1 & risk_period == "risk" & current_dose >= 2, current_dose, 0),
         eventweight_week = ifelse(death == 1 & risk_period == "risk" & current_dose >= 2,  week + (current_dose-1)*risk_period_length , 0)) %>%
      mutate(eventweight_period = case_when(eventweight_period== 0 ~ 0,
                                            eventweight_period>0 & vaccine_vector == "mRNA" ~ eventweight_period*2 - 1,
                                            eventweight_period>0 & vaccine_vector == "not_mRNA_or_unknown" ~ eventweight_period*2),
             eventweight_week = case_when(eventweight_week== 0 ~ 0,
                                            eventweight_week>0 & vaccine_vector == "mRNA" ~ eventweight_week*2 - 1,
                                            eventweight_week>0 & vaccine_vector == "not_mRNA_or_unknown" ~ eventweight_week*2)) %>%
    mutate(risk_period_dose = ifelse(risk_period == "risk" & current_dose >= 2, "baseline", risk_period_dose),
         risk_week = ifelse(risk_period == "risk" & current_dose >= 2, "baseline", risk_week),
         risk_week_dose = ifelse(risk_period == "risk" & current_dose >= 2, "baseline", risk_week_dose),
         risk_period = ifelse(risk_period == "risk" & current_dose >= 2, "baseline", risk_period),
          risk_week_dose_vector = ifelse(risk_week_dose_vector == "risk" & current_dose >=2, "baseline", risk_week_dose_vector),
         risk_period_dose_vector = ifelse(risk_period_dose_vector == "risk" & current_dose >=2, "baseline", risk_period_dose_vector),
         risk_week_vector = ifelse(risk_week_vector == "risk" & current_dose >=2, "baseline", risk_week_vector),
         risk_period_vector = ifelse(risk_period_vector == "risk" & current_dose >=2, "baseline", risk_period_vector),
         id2 = id + max(df_all$id))

  # combine 
  df_clean <-rbind(df_all, df_d1, df_d2, df_d3)
  
  df_clean <- df_clean %>%
    mutate(eventweight_period = ifelse(vaccine_vector != vector, 0, eventweight_period),
          eventweight_week = ifelse(vaccine_vector != vector, 0, eventweight_week),
          risk_week_dose_vector = ifelse(vaccine_vector != vector, "baseline", risk_week_dose_vector),                     
          risk_period_dose_vector = ifelse(vaccine_vector != vector, "baseline", risk_period_dose_vector),
          risk_week_vector = ifelse(vaccine_vector != vector, "baseline", risk_week_vector),                     
          risk_period_vector = ifelse(vaccine_vector != vector, "baseline", risk_period_vector))
  
# chaneg the weights when combining doses to be just dependent on the week, not the dose
  df_clean <-df_clean %>%
      mutate(eventweight_period_combined = case_when(eventweight_period > 0 & vaccine_vector == "mRNA" ~ 1,
                                                    eventweight_period > 0 & vaccine_vector == "not_mRNA_or_unknown" ~ 2,
                                                    TRUE ~ 0),
            eventweight_week_combined = case_when(eventweight_week > 0  & vaccine_vector == "mRNA" ~ week*2-1,
                                                 eventweight_week > 0  & vaccine_vector == "not_mRNA_or_unknown" ~ week*2,
                                                 TRUE ~ 0))

  

  df_clean$risk_period_dose <- relevel(as.factor(df_clean$risk_period_dose), ref="baseline")
  df_clean$risk_week_dose <- relevel(as.factor(df_clean$risk_week_dose), ref="baseline")
  df_clean$risk_period <- relevel(as.factor(df_clean$risk_period), ref="baseline")
  df_clean$risk_week <- relevel(as.factor(df_clean$risk_week), ref="baseline")
  df_clean$risk_week_dose_vector <- relevel(as.factor(df_clean$risk_week_dose_vector), ref="baseline")
  df_clean$risk_week_vector <- relevel(as.factor(df_clean$risk_week_vector), ref="baseline")
  df_clean$risk_period_vector <- relevel(as.factor(df_clean$risk_period_vector), ref="baseline")
  df_clean$risk_period_dose_vector <- relevel(as.factor(df_clean$risk_period_dose_vector), ref="baseline")
  
  return(df_clean)
}




# breakdowns
data_for_risk_period_breakdown <- function(df, breakdown, cause_name, exposure){
  
  breakdown_values <- unique(df[,breakdown])
  
  df_clean_breakdown <- df[df[,breakdown]==breakdown_values[1],]

  # all_deaths
  all <- eventdepexp_model(outcome = "death",
                   exposure = exposure,
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean_breakdown,
                   model =F)

  all$breakdown_value <- breakdown_values[1]
  
  for (i in breakdown_values[-1]){

    df_clean_breakdown <- df[df[,breakdown]==i,]

    # all_deaths
    all2 <- eventdepexp_model(outcome = "death",
                   exposure = exposure,
                   x= "ns(cal_day_week, df=3, Boundary.knots = quantile(cal_day_week, c(.10,.90)))",
                   data = df_clean_breakdown,
                   model =F)

    all2$breakdown_value <- i
    
    all <- rbind(all, all2)
    
  }
  
  all$cause <- cause_name
  all$breakdown <- breakdown
  
  all <- all %>% 
    mutate(dose = case_when(grepl("D1",term) ~ "First",
                          grepl("D2",term) ~ "Second",
                          grepl("D3",term) ~ "Third",
                         TRUE ~ 'All doses'))

  return(all)
}
    

make_risk_period_variables <- function(df_week, risk_period_length){
# create variables to indicate whether in risk period or not, and the week in the risk period
  df_week <- df_week %>%
    mutate(risk_period = case_when(week <= risk_period_length & current_dose != 0 ~ 'risk', 
                                   week > risk_period_length & current_dose != 0 ~ 'baseline',
                                  current_dose == 0 ~ 'baseline')) %>%
    mutate(risk_week = case_when(week <= risk_period_length & current_dose != 0 ~ paste0("week_", str_pad(week, 2, pad="0")), 
                                   week > risk_period_length & current_dose != 0 ~ 'baseline',
                                  current_dose == 0 ~ 'baseline'))
  
# create variables with dose information added to the risk period and risk weeks
  df_week <- df_week %>%
    mutate(risk_period_dose = case_when(risk_period == "risk" ~ paste0("D", current_dose),
                              TRUE ~ "baseline"),
        risk_week_dose = case_when(risk_period == "risk" ~ paste0("D", current_dose, "_week_", str_pad(week, 2, pad="0")),
                              TRUE ~ "baseline"))
  return(df_week)
}

make_risk_period_variables_infections <- function(df_week, risk_period_length){
# create variables to indicate whether in risk period or not, and the week in the risk period
  df_week <- df_week %>%
    mutate(risk_period = case_when(week <= risk_period_length & infection_number != 0 ~ 'risk', 
                                   week > risk_period_length & infection_number != 0 ~ 'baseline',
                                  infection_number == 0 ~ 'baseline')) %>%
    mutate(risk_week = case_when(week <= risk_period_length & infection_number != 0 ~ paste0("week_", str_pad(week, 2, pad="0")), 
                                   week > risk_period_length & infection_number != 0 ~ 'baseline',
                                  infection_number == 0 ~ 'baseline'))
  
  return(df_week)
}

