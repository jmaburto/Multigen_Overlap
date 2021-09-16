# =================================================================================
# This is the code for the multigenerational overlapping
# project with Christiaan Monden
# 2021/09/08
# Calculate the expected years of multigenerational overlapping from a cohort perspective
# for all countries available at the HFD and HMD
# =================================================================================

rm(list=ls())
library(tidyverse)
Mycol <- c("#08306B", "#238B45", "#FD8D3C", "#D4B9DA", "#FFEDA0")

source("Function/Func_MxtoLT.R")
source("Function/Func_KGBbc.R")
source("Function/Func_KGBper.R")

#### ---- Cohort ----
DNK_bc <- Fun_KGBbc(CountryName = "DNK", agelim = 90)
FIN_bc <- Fun_KGBbc(CountryName = "FIN", agelim = 90)
SWE_bc <- Fun_KGBbc(CountryName = "SWE", agelim = 90)
CHE_bc <- Fun_KGBbc(CountryName = "CHE", agelim = 90)
GBR_bc <- Fun_KGBbc(CountryName = "GBRTENW", agelim = 90)

result_cohort <- DNK_bc %>% 
  rbind(FIN_bc) %>% 
  rbind(SWE_bc) %>% 
  rbind(CHE_bc) %>%
  rbind(GBR_bc) %>% 
  mutate(Country = case_when(Country == "DNK" ~ "Denmark",
                             Country == "FIN" ~ "Finland",
                             Country == "SWE" ~ "Sweden",
                             Country == "CHE" ~ "Switzerland",
                             Country == "GBRTENW" ~ "England & Wales"))
write.csv(result_cohort, "Results/Cohort_Multigen.csv", row.names = F)

result_cohort %>% 
  gather(key = index, value = duration, c("Cohort_EP1", "diff")) %>% 
  mutate(Country_index = paste0(Country, index),
         index = ifelse(index == "Cohort_EP1", "Expected duration", "Simple difference")) %>% 
  ggplot(aes(x = bc_gm, y = duration, group = Country_index, colour = Country, linetype = index)) +
  geom_line(size = 1.3) +
  scale_colour_manual(values = Mycol) +
  labs(x = "Birth cohort of grandmother", y = "Overlapping duration between grandmother and granddaughter") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box = "vertical")
ggsave("Graph/CohortEP1_Diff_all.pdf", width = 8, height = 6.5)

result_cohort %>% 
  ggplot(aes(x = bc_gm, y = Cohort_ER1, colour = Country)) +
  geom_line(size = 1.3) +
  scale_colour_manual(values = Mycol) +
  labs(x = "Birth cohort of grandmother", y = "Expected number of years that granddaughter aged 10-19 spend life with maternal grandmother") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box = "vertical")
ggsave("Graph/CohortER1_all.pdf", width = 8, height = 6.5)

#### ---- Period ----
DNK_per <- Fun_KGBper(CountryName = "DNK")
FIN_per <- Fun_KGBper(CountryName = "FIN")
SWE_per <- Fun_KGBper(CountryName = "SWE")
CHE_per <- Fun_KGBper(CountryName = "CHE")
GBR_per <- Fun_KGBper(CountryName = "GBRTENW")

result_period <- DNK_per %>% 
  rbind(FIN_per) %>% 
  rbind(SWE_per) %>% 
  rbind(CHE_per) %>%
  rbind(GBR_per) %>% 
  mutate(Country = case_when(Country == "DNK" ~ "Denmark",
                             Country == "FIN" ~ "Finland",
                             Country == "SWE" ~ "Sweden",
                             Country == "CHE" ~ "Switzerlad",
                             Country == "GBRTENW" ~ "England & Wales"))
write.csv(result_cohort, "Results/Period_Multigen.csv", row.names = F)

result_period %>% 
  gather(key = index, value = duration, c("Period_EP1", "diff")) %>% 
  mutate(Country_index = paste0(Country, index),
         index = ifelse(index == "Period_EP1", "Expected duration", "Simple difference")) %>% 
  ggplot(aes(x = Year_gm, y = duration, group = Country_index, colour = Country, linetype = index)) +
  geom_line(size = 1.3) +
  scale_colour_manual(values = Mycol) +
  labs(x = "Birth cohort of grandmother", y = "Overlapping duration between grandmother and granddaughter") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box = "vertical")

#### ---- Period vs Cohort ----
result_cohort2 <- result_cohort %>% 
  select(Country, bc_gm, Cohort_EP1, Cohort_ER1)

result_period2 <- result_period %>% 
  select(Country, Year_gm, Period_EP1)

result_cohort2 %>% 
  left_join(result_period2, by = c("Country", "bc_gm" = "Year_gm")) %>% 
  filter(Period_EP1 > 0) %>% 
  gather(key = Per_Coh, value = EP1, c("Period_EP1", "Cohort_EP1")) %>% 
  ggplot(aes(x = bc_gm, y = EP1, group = Per_Coh, colour = Per_Coh)) +
  geom_line(size = 1.3) +
  facet_wrap(~ Country) +
  scale_colour_manual(values = c(Mycol[3], Mycol[1]),
                      labels = c("Cohort", "Period")) +
  labs(x = "Time period / birth cohort", y = "Overlapping duration between grandmother and granddaughter") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("Graph/Comp_perbc_all.pdf", width = 8, height = 6.5)

#### ---- Cohort EP1 vs Cohort ER1 ----
result_cohort %>% 
  gather(key = index, value = value, c("Cohort_EP1", "Cohort_ER1")) %>% 
  ggplot(aes(x = bc_gm, y = value, colour = index)) +
  geom_line(size = 1.3) +
  facet_wrap(~ Country) +
  scale_colour_manual(values = c(Mycol[3], Mycol[1]),
                      labels = c("From grandmother", "From teenager")) +
  labs(x = "Birth cohort of grandmother", y = "Overlapping duration between grandmother and granddaughter") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box = "vertical")
ggsave("Graph/CohortEP1ER1_all.pdf", width = 8, height = 6.5)
