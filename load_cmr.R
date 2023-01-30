rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
raw.cmr <- read.csv("SES_PM25_CMR_data/County_annual_PM25_CMR.csv")
pm.25.avg <- raw.cmr  %>% group_by(by=FIPS) %>% summarise(avg = mean(PM2.5))
colnames(pm.25.avg) <- c("FIPS", "PM2.5")

cmr.2010 <- raw.cmr[raw.cmr$Year == 2010, c("FIPS", "CMR")]
colnames(cmr.2010) <- c("FIPS", "CMR.2010")
cmr.1990 <- raw.cmr[raw.cmr$Year == 1990, c("FIPS", "CMR")]
colnames(cmr.1990) <- c("FIPS", "CMR.1990")
cmr.data <- merge(cmr.2010, cmr.1990, by="FIPS")
cmr.data$CMR.diff <- cmr.data$CMR.2010 / cmr.data$CMR.1990
cmr.data <- cmr.data[, c("FIPS", "CMR.diff")]
colnames(cmr.data) <- c("FIPS", "CMR")

raw.var <- read.csv("SES_PM25_CMR_data/County_RAW_variables.csv")
cov.2000 <- raw.var[, c("FIPS", "civil_unemploy_2000", "median_HH_inc_2000", 
                        "femaleHH_ns_pct_2000", "vacant_HHunit_2000", 
                        "owner_occ_pct_2000", "eduattain_HS_2000",
                        "healthfac_2005",
                        "pctfam_pover_2000", "healthfac_1999", "population_2000")]
cov.1990 <- raw.var[, c("FIPS", "civil_unemploy_1990", "median_HH_inc_1990", 
                        "femaleHH_ns_pct_1990", "vacant_HHunit_1990", 
                        "owner_occ_pct_1990", "eduattain_HS_1990",
                        "pctfam_pover_1990")]

raw.diff <- raw.var %>% group_by(FIPS) %>% 
  summarise(civil_unemploy_2010 = civil_unemploy_2010-civil_unemploy_1990,
            median_HH_inc_diff = median_HH_inc_2010-median_HH_inc_1990,
            femaleHH_ns_pct_diff = femaleHH_ns_pct_2010-femaleHH_ns_pct_1990,
            vacant_HHunit_diff = vacant_HHunit_2010 - vacant_HHunit_1990,
            owner_occ_pct_diff = owner_occ_pct_2010 - owner_occ_pct_1990,
            eduattain_HS_diff = eduattain_HS_2010 - eduattain_HS_1990, 
            pctfam_pover_diff = pctfam_pover_2010 - pctfam_pover_1990)
raw.cov <- raw.var[, c("FIPS", "healthfac_2005_1999", "population_2000")]

df_list <- list(pm.25.avg, cmr.data, cov.2000, cov.1990, cmr.1990)
df_list <- df_list %>% reduce(full_join, by='FIPS')
df_list <- data.frame(df_list)
write.csv(df_list, "data/CMR_1990_2010_y.csv")
