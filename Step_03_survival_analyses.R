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
library(survival)
library(survminer)
library(misty)

# Define directories
results_simul <- ("~/results")

#---------------------------------------------------------------------------------#
# Load simulation data

setwd(results_simul)

data_simul <- read_csv("dataready.csv", show_col_types = FALSE)

data_simul2 <- data_simul %>% 
  mutate(flagdt = ifelse(event == 0 & DT == 1, 1, 0)) %>% 
  group_by(id)  %>% 
  mutate(flagdt = max(flagdt)) %>% 
  ungroup() %>%
  filter(flagdt == 1)  %>%
  filter(event >= 0)  %>%
  mutate(flagds = ifelse(DS == 1, event, 30))  %>%
  group_by(id)  %>% 
  mutate(flagds = min(flagds)) %>% 
  ungroup() %>% 
  mutate(censor = ifelse(flagds < 30, 1, 0))  %>% 
  filter(event == 0)


fit <-survfit(Surv(flagds, censor) ~ counterfactual,
              data = data_simul2)

ggsurvplot(fit, data = data_simul2)


data_export <- data.frame(time = fit$time,
                          surv = fit$surv,
                          std = fit$std.err,
                          lower = fit$lower, 
                          upper = fit$upper,
                          group = 1)

data_export$group[1:fit$strata[1]] <- 0

write.csv(data_export, "survival.csv")

#################################################
# Heterogeneity analysis - age 

data_simul2_younger <- data_simul2 %>%
  filter(SH <= 20)

fit <-survfit(Surv(flagds, censor) ~ counterfactual,
              data = data_simul2_younger)

ggsurvplot(fit, data = data_simul2_younger)


data_export <- data.frame(time = fit$time,
                          surv = fit$surv,
                          std = fit$std.err,
                          lower = fit$lower, 
                          upper = fit$upper,
                          group = 1)

data_export$group[1:fit$strata[1]] <- 0

write.csv(data_export, "survival_younger.csv")


data_simul2_older <- data_simul2 %>%
  filter(SH > 20)

fit <-survfit(Surv(flagds, censor) ~ counterfactual,
              data = data_simul2_older)

ggsurvplot(fit, data = data_simul2_older)


data_export <- data.frame(time = fit$time,
                          surv = fit$surv,
                          std = fit$std.err,
                          lower = fit$lower, 
                          upper = fit$upper,
                          group = 1)

data_export$group[1:fit$strata[1]] <- 0

write.csv(data_export, "survival_older.csv")

#################################################
# Heterogeneity analysis - income 
median_val <- median(filter(data_simul, event == -1)$YT)

data_simul <- data_simul %>%
  mutate(income_flag = ifelse(event == -1 & YT <= median_val, 0, 1)) %>%
  mutate(income_flag = ifelse(event != -1, 99, income_flag)) %>%
  group_by(id) %>% 
  mutate(income_flag = min(income_flag)) %>% 
  ungroup()

data_simul2 <- data_simul %>% 
  mutate(flagdt = ifelse(event == 0 & DT == 1, 1, 0)) %>% 
  group_by(id)  %>% 
  mutate(flagdt = max(flagdt)) %>% 
  ungroup() %>%
  filter(flagdt == 1)  %>%
  filter(event >= 0)  %>%
  mutate(flagds = ifelse(DS == 1, event, 30))  %>%
  group_by(id)  %>% 
  mutate(flagds = min(flagds)) %>% 
  ungroup() %>% 
  mutate(censor = ifelse(flagds < 30, 1, 0))  %>% 
  filter(event == 0)

data_simul2_poorer <- data_simul2 %>%
  filter(income_flag == 0)

fit <-survfit(Surv(flagds, censor) ~ counterfactual,
              data = data_simul2_poorer)

ggsurvplot(fit, data = data_simul2_poorer)


data_export <- data.frame(time = fit$time,
                          surv = fit$surv,
                          std = fit$std.err,
                          lower = fit$lower, 
                          upper = fit$upper,
                          group = 1)

data_export$group[1:fit$strata[1]] <- 0

write.csv(data_export, "survival_poorer.csv")


data_simul2_richer <- data_simul2 %>%
  filter(income_flag == 1)

fit <-survfit(Surv(flagds, censor) ~ counterfactual,
              data = data_simul2_richer)

ggsurvplot(fit, data = data_simul2_richer)


data_export <- data.frame(time = fit$time,
                          surv = fit$surv,
                          std = fit$std.err,
                          lower = fit$lower, 
                          upper = fit$upper,
                          group = 1)

data_export$group[1:fit$strata[1]] <- 0

write.csv(data_export, "survival_richer.csv")


#################################################
# No transaction costs

setwd(results_simul)

data_simul <- read_csv("dataready_app.csv", show_col_types = FALSE)

data_simul2 <- data_simul %>% 
  mutate(flagdt = ifelse(event == 0 & DT == 1, 1, 0)) %>% 
  group_by(id)  %>% 
  mutate(flagdt = max(flagdt)) %>% 
  ungroup() %>%
  filter(flagdt == 1)  %>%
  filter(event >= 0)  %>%
  mutate(flagds = ifelse(DS == 1, event, 30))  %>%
  group_by(id)  %>% 
  mutate(flagds = min(flagds)) %>% 
  ungroup() %>% 
  mutate(censor = ifelse(flagds < 30, 1, 0))  %>% 
  filter(event == 0)


fit <-survfit(Surv(flagds, censor) ~ counterfactual,
              data = data_simul2)


data_export <- data.frame(time = fit$time,
                          surv = fit$surv,
                          std = fit$std.err,
                          lower = fit$lower, 
                          upper = fit$upper,
                          group = 1)

data_export$group[1:fit$strata[1]] <- 0

write.csv(data_export, "survival_app.csv")