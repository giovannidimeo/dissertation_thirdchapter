library(haven)
library(did)
library(tidyverse)
library(xtable)
library(foreign)
library(extrafont)
library(Cairo)
library(cowplot)
library(purrr)
library(here)
library(plm)
library(lfe)
library(readxl)
library(future)
library(furrr)
library(tictoc)
library(data.table)
library(purrr)
library(dplyr)
options(datatable.rbindlist.check = "none") # To silence the rbindlist() related warning. The warning has no effect on the estimation, and occurs only because "did" uses data.table in the latest version.

# Define directories
results_simul <- ("~/results")
datasets_shp <- ("~")

#---------------------------------------------------------------------------------#
# Load simulation data

setwd(results_simul)

data_simul <- read_csv("dataready.csv", show_col_types = FALSE)

data_simul <- data_simul %>% 
  mutate(shock = ifelse(counterfactual == 0, 0, SH)) %>%
  filter(time <= 40) %>% 
  mutate(YT = YT/10)

#---------------------------------------------------------------------------------#
# Load SHP data 

setwd(datasets_shp)

data_shp <- read_dta("shp_allwaves_ready_hh.dta")

data_shp <- data_shp %>% 
  mutate(age_shock = ifelse(year == shock, max_age, 0))  %>%
  mutate(age_shock = ifelse(shock == 0, 0, age_shock))  %>%
  group_by(idhous) %>% 
  mutate(age_shock = max(age_shock)) %>%
  ungroup() %>%
  filter(age_shock <= 65) %>%
  filter(max_age <= 65) %>%
  mutate(DT = owner) %>% 
  mutate(time = year) %>% 
  mutate(id = idhous) %>% 
  mutate(npers = ifelse(single == 1, 1, 2)) %>% 
  mutate(YT = income_hh/npers/100000)

data_shp_nopreshocks  <- data_shp %>%  # to get rid of people returning to previous income
  filter(shock == 0 | (shock !=0 & increase_in_neg1 == 0)) 

#---------------------------------------------------------------------------------#
# Define functions
span_top <- 5
span_bottom <- 5
levels <- c(0.01, 0.05, 0.1)

cs_fun <- function(y, df_data, xformla_input, alpha) {
  
  set.seed(1215) # Setting seed in the function. This should prevent us generating results w/o explicitly setting the right seed. 
  
  # Estimate the ATT(g,t)
  result <- att_gt(yname = y,
                   gname = "shock",
                   idname = "id",
                   tname = "time",
                   est_method = "ipw",
                   control_group = "nevertreated",
                   panel = TRUE,
                   allow_unbalanced_panel = TRUE,
                   pl = FALSE,
                   xformla = xformla_input,
                   cores = 12,
                   alp = alpha,
                   base_period = "universal",
                   anticipation = 0,
                   data = df_data)
  
  # Aggregate to ES coefficients
  es <- aggte(result, type = "dynamic", na.rm = TRUE, min_e = -span_bottom, max_e = span_top)
  
  return(list(result = result, es = es))
}

xformla_shp <- ~ avg_age + single
xformla_simul <- ~ 1

varnames <- c("YT", "DT")

res <- cs_fun("DT", data_shp_nopreshocks, xformla_shp, 0.05)

#---------------------------------------------------------------------------------#
# Run functions
analysis_settings <- list(
  results_shp = list(
    data = data_shp_nopreshocks,  # DEFINE DATA SET
    varnames = varnames,  # DEFINE VARIABLE LIST
    xformla = xformla_shp  # DEFINE CONTROL VARIABLE
  ),
  
  results_simul = list(
    data = data_simul,  # DEFINE DATA SET
    varnames = varnames,  # DEFINE VARIABLE LIST
    xformla = xformla_simul  # DEFINE CONTROL VARIABLE
  )
)


# Initialize an empty list to store results for all variables and significance levels
MainResults_Megaloop <- list()

# Use future_map to iterate over analysis types in parallel
results_list <- future_map(names(analysis_settings), ~{
  type <- .x
  current_setting <- analysis_settings[[type]]
  current_data <- current_setting$data
  current_xformla <- current_setting$xformla
  current_varnames <- current_setting$varnames
  
  # Initialize a list to store results for the current analysis type
  analysis_results <- list()
  
  # Iterate over variables (this could also be parallelized with future_map if desired)
  for (var in current_varnames) {
    es_01 <- NULL
    es_05 <- NULL
    es_10 <- NULL
    
    # Run the analysis for each significance level
    for (l in levels) {
      result <- cs_fun(var, current_data, current_xformla, l)
      es <- result$es
      
      # Store results based on significance level
      if (l == 0.01) {
        es_01 <- es
      } else if (l == 0.05) {
        es_05 <- es
      } else if (l == 0.1) {
        es_10 <- es
      }
    }
    
    # Create the table for a given variable
    table_name <- paste0(type, "_", var)
    table_data <- data.frame(
      group = es_05$egt,
      effect = es_05$att.egt,
      se = es_05$se.egt,
      lower = es_05$att.egt - es_05$crit.val.egt * es_05$se.egt,
      upper = es_05$att.egt + es_05$crit.val.egt * es_05$se.egt,
      critval99 = es_01$crit.val.egt,
      critval95 = es_05$crit.val.egt,
      critval90 = es_10$crit.val.egt,
      att = es_05$overall.att,
      attse = es_05$overall.se,
      unique_pids = es_05$DIDparams$n
    )
    # Store the table in the analysis_results list
    analysis_results[[table_name]] <- table_data
  }
  analysis_results
}, .progress = TRUE)  # Optional: show progress

# Combine results from all analyses into a single list
MainResults_Megaloop <- do.call(c, results_list)

toc()


setwd(results_simul)
objects_to_save <- ls(pattern = "^(MainResults_Megaloop)")
save(list = objects_to_save, file = paste0("AllResults_new", ".RData"))


#---------------------------------------------------------------------------------#
# Define function for stars
myround <- function(x) {
  ifelse(x < 1 & x > -1, round(x, 4), round(x, 2) )
}

i <- 1L

stars <- function(df) {
  df$tstat[i] <- 0
  df$tstatatt[i] <- 0
  for (i in 1:11) {
    print(i)
    df$tstat[i] <-   abs(df$effect[i]/df$se[i])
    if (!is.na(df$tstat[i]) &  df$tstat[i] >= df$critval99[i]) {
      df$stars[i] <- "***"
    } else if ( !is.na(df$tstat[i])  &  df$tstat[i] < df$critval99[i] & df$tstat[i] >= df$critval95[i]) {
      df$stars[i] <- "**"
    } else if (!is.na(df$tstat[i])  &  df$tstat[i] < df$critval95[i] & df$tstat[i] >= df$critval90[i]) {
      df$stars[i] <- "*"
    } else { 
      df$stars[i] <- "" }
    
    df$tstatatt[i] <-abs(df$att[i]/df$attse[i])
    
    if (!is.na(df$tstatatt[i]) &  df$tstatatt[i] >= qnorm(1-0.01/2)) {
      df$starsatt[i] <- "***"
    } else if ( !is.na(df$tstatatt[i])  &  df$tstatatt[i] < qnorm(1-0.01/2) & df$tstatatt[i] >= qnorm(1-0.05/2)) {
      df$starsatt[i] <- "**"
    } else if (!is.na(df$tstatatt[i])  &  df$tstatatt[i] < qnorm(1-0.05/2) & df$tstatatt[i] >= qnorm(1-0.1/2)) {
      df$starsatt[i] <- "*"
    } else { 
      df$starsatt[i] <- "" }
    
  }
  
  for (i in 1:11) {
    df$effect[i] <-  paste0(toString(myround(as.numeric(df$effect[i]))),  df$stars[i], sep = "")
    df$att[i] <-  paste0(toString(myround(as.numeric(df$att[i]))),  df$starsatt[i], sep = "")
    
  }
  return(df)
}

# Adjusted create_info_text function to eliminate spaces around parentheses
create_info_text <- function(df) {
  df_processed <- stars(df)
  
  # Extract 'att' value and its stars separately
  att_value <- as.numeric(gsub("[^0-9.-]", "", df_processed$att[1])) # Keeps numeric value, including negative
  att_stars <- gsub("[0-9.-]", "", df_processed$att[1]) # Keeps stars
  
  # Apply myround to 'att' numeric value and combine with stars
  att_rounded_with_stars <- paste0(myround(att_value), att_stars)
  
  # Apply myround to 'attse' directly
  attse_rounded <- myround(as.numeric(df_processed$attse[1]))
  
  # Extract number of unique individuals
  unique <- df_processed$unique_pids[1]
  
  info_text <- paste0("\n", "ATT: ", att_rounded_with_stars, " (", attse_rounded, ")")
  
  return(info_text)
}

# Apply this adjusted function to each table in MainResults_Megaloop and collect the results
info_texts <- lapply(MainResults_Megaloop, create_info_text)

# Name the elements of info_texts for easy reference
names(info_texts) <- names(MainResults_Megaloop)

# _________________________________________________________________________________#
# CREATE TABLES
#---------------------------------------------------------------------------------#

myround <- function(x) {
  ifelse(x < 1 & x > -1, round(x, 4), round(x, 2) )
}

i <- 1L

stars <- function(df) {
  df$tstat[i] <- 0
  df$tstatatt[i] <- 0
  for (i in 1:11) {
    print(i)
    df$tstat[i] <-   abs(df$effect[i]/df$se[i])
    if (!is.na(df$tstat[i]) &  df$tstat[i] >= df$critval99[i]) {
      df$stars[i] <- "***"
    } else if ( !is.na(df$tstat[i])  &  df$tstat[i] < df$critval99[i] & df$tstat[i] >= df$critval95[i]) {
      df$stars[i] <- "**"
    } else if (!is.na(df$tstat[i])  &  df$tstat[i] < df$critval95[i] & df$tstat[i] >= df$critval90[i]) {
      df$stars[i] <- "*"
    } else { 
      df$stars[i] <- "" }
    
    df$tstatatt[i] <-abs(df$att[i]/df$attse[i])
    
    if (!is.na(df$tstatatt[i]) &  df$tstatatt[i] >= qnorm(1-0.01/2)) {
      df$starsatt[i] <- "***"
    } else if ( !is.na(df$tstatatt[i])  &  df$tstatatt[i] < qnorm(1-0.01/2) & df$tstatatt[i] >= qnorm(1-0.05/2)) {
      df$starsatt[i] <- "**"
    } else if (!is.na(df$tstatatt[i])  &  df$tstatatt[i] < qnorm(1-0.05/2) & df$tstatatt[i] >= qnorm(1-0.1/2)) {
      df$starsatt[i] <- "*"
    } else { 
      df$starsatt[i] <- "" }
    
  }
  
  for (i in 1:11) {
    #df$effect[i] <-  myround(df$effect[i])
    df$effect[i] <-  paste0(toString(myround(as.numeric(df$effect[i]))),  df$stars[i], sep = "")
    df$att[i] <-  paste0(toString(myround(as.numeric(df$att[i]))),  df$starsatt[i], sep = "")
    
  }
  return(df)
}

tab_function <- function(df, name){  
  df <- stars(df)
  
  # Initialize the first row (coefficients)
  table_tex1 <- data.frame("Dependent variable" = as.character(name),
                           check.names = FALSE)
  
  # Initialize the second row (standard errors)
  table_tex2 <- data.frame("Dependent variable" = "",
                           check.names = FALSE)
  
  # Add columns 0 to 30
  for (i in 1:11) {
    table_tex1[[as.character(i)]] <- df$effect[i]
    table_tex2[[as.character(i)]] <- paste0("(", myround(df$se[i]), ")")
  }
  
  # Add final columns
  table_tex1$Overall <- df$att[1]
  table_tex1$`Unique obs.` <- df$unique_pids[1]
  
  table_tex2$Overall <- paste0("(", myround(df$attse[1]), ")")
  table_tex2$`Unique obs.` <- ""
  
  tab <- rbind(table_tex1, table_tex2)
  return(tab)
}

generate_table <- function(list) {
  
  tab_temp <- data.frame(matrix(nrow = 25, ncol = length(list))) 
  
  for (v in 1:length(list)) {
    name <- paste0("(", as.character(v), ")")
    tab_or <- tab_function(list[[v]], name)
    for (m in 1:12) {
      k = m*2-1
      j = m+1
      tab_temp[[v]][k] <- tab_or[[j]][1]
      tab_temp[[v]][k+1] <- tab_or[[j]][2]
    }
    tab_temp[[v]][25] <-  tab_or[["Unique obs."]][1]
    colnames(tab_temp)[v] = name
  }
  
  tab_temp["Event time"] <- ""
  tab_temp  <- tab_temp %>% relocate("Event time")
  
  for (p in 1:12) {
    q = p*2 - 1
    tab_temp[["Event time"]][q] <- as.character(p-6)
    tab_temp[["Event time"]][q+1] = " "
  }
  tab_temp[["Event time"]][23] = "Overall"
  tab_temp[["Event time"]][25] = "Unique obs."
  
  print(xtable(tab_temp, type = "latex"), include.rownames=FALSE)
  print(tab_temp)
}



#Generate the table
var_list <- list( "(1)" = MainResults_Megaloop$results_shp_YT,
                  "(2)" = MainResults_Megaloop$results_simul_YT,
                  "(3)" = MainResults_Megaloop$results_shp_DT,
                  "(4)" = MainResults_Megaloop$results_simul_DT)




generate_table(var_list)


