# VDPV2-n analyses
rm(list=ls())

#### Load Libraries ####
library(readr)
library(ggridges)
library(tidyverse)
library(lubridate)
library(PolisAPI)
library(ggrepel)
library(GGally)
library(ggh4x)
library(ggbreak)
library(gganimate)
library(patchwork)

#### Load workspace ####
load("VDPV2n_analyses.RData")

#### Helper functions ####
#Create region from ADM_0
Func_region = function(ADM_0){
  EMRO <- c("AFGHANISTAN", "EGYPT", "IRAQ", "PAKISTAN", "SOMALIA", "SUDAN", 
            "DJIBOUTI", "ERITREA", "JORDAN", "TÜRKIYE",
            "IRAN (ISLAMIC REPUBLIC OF)", "LIBYA", "SYRIAN ARAB REPUBLIC", "TUNISIA",
            "YEMEN","LEBANON", "OCCUPIED PALESTINIAN TERRITORY, INCLUDING EAST JERUSALEM")
  WPRO <- c("PHILIPPINES", "MALAYSIA","LAO PEOPLE'S DEMOCRATIC REPUBLIC", "CHINA", "PAPUA NEW GUINEA")
  SEARO <- c("INDONESIA", "MYANMAR", "INDIA", "NEPAL")
  EURO <- c("TAJIKISTAN", "GEORGIA","GERMANY", "RUSSIAN FEDERATION", "UKRAINE", "ISRAEL", "POLAND", "THE UNITED KINGDOM")
  AMRO <- c("UNITED STATES OF AMERICA", "CANADA", "ARGENTINA", "COLOMBIA", "GUATEMALA", "PERU")
  ifelse(ADM_0 %in% EMRO, "EMRO",
         ifelse(ADM_0 %in% WPRO, "WPRO",
                ifelse(ADM_0 %in% SEARO, "SEARO",
                       ifelse(ADM_0 %in% EURO, "EURO",
                              ifelse(ADM_0 %in% AMRO, "AMRO",
                              "AFRO")))))
}

Func_africa = function(region, ADM_0){
  ifelse(region == "AFRO" | ADM_0 %in% c("EGYPT", "SOMALIA", "SUDAN", "DJIBOUTI", "ERITREA", "LIBYA", "TUNISIA"), "Africa", "Elsewhere")
}

Func_period_to_quarter = function(period){
  decimal = period - floor(period)
  if_else(decimal < 0.333, 0,
          if_else(decimal < 0.5, 0.25,
                  if_else(decimal < 0.750, 0.5, 0.75)))
}

#### Load Target Pops from disaggregate_polis_pop ####
# After dvc pull latest polio immunity mapping results, run disaggregate_polis_pop.R
# Cross-check numbers against nOPV2 campaign tracker

polis_pops <- readRDS("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/sia_polis_target_pop.rds")
names(polis_pops)[13] <- "GUID"

# Add fields
polis_pops$period <- year(polis_pops$start_date) + (month(polis_pops$start_date)-1)/12
polis_pops$quarter <- year(polis_pops$start_date) + (quarter(polis_pops$start_date)-1)/4
polis_pops$source <- ifelse(polis_pops$vaccinetype == "nOPV2", "nOPV2", "Sabin2")

# unique(polis_pops$adm0_name)
polis_pops$region <- Func_region(polis_pops$adm0_name)
polis_pops$africa <- Func_africa(polis_pops$region, polis_pops$adm0_name)
polis_pops %>% group_by(adm0_name) %>% summarize(region = unique(region), africa = unique(africa)) %>%View()

# Mark campaigns that did happen:
polis_pops_raw <- polis_pops
# done_campaigns <- c("KEN-2023-003","MRT-2022-003", "TZA-2023-002") # Are marked as done now
# polis_pops_raw %>% filter(parentactivitycode %in% done_campaigns) %>% View()
# polis_pops[polis_pops$parentactivitycode %in% done_campaigns, "status"] <- "Done"
polis_pops <- polis_pops %>% filter(status == "Done")

# Manually remove campaigns that didn't happen
# not_campaigns <- c( "DZA-2024-001", "NER-2023-006", "NER-2024-003")  

polis_pops %>% filter(parentactivitycode %in% not_campaigns) %>% View()
polis_pops[polis_pops$parentactivitycode %in% not_campaigns, "target_pop"] <- 0

# Count vaccine used
polis_pops %>% ungroup() %>%
  filter(vaccinetype == "nOPV2")  %>%
  summarize(pop= sum(target_pop, na.rm = T))

#### Compare polis to nOPV2 vaccine tracker (deprecated: no longer being manually tracked) ####
# Import nOPV2 tracker
nOPV2_Tracker <- read_csv("nOPV2_Tracker.csv")
nOPV2_Tracker[nOPV2_Tracker$adm0_name == "COTE D'IVOIRE", "adm0_name"] <- "CÔTE D’IVOIRE"
head(nOPV2_Tracker, 10)
nOPV2_Tracker %>% summarize(doses = sum(doses))
nOPV2_Tracker %>% filter(region == "AFRO") %>% summarize(doses = sum(doses))

# check target pop compared with tracker
polis_pops %>% group_by(adm0_name, region) %>%
  filter(vaccinetype == "nOPV2") %>%
  summarize(target_pop = sum(target_pop, na.rm=T)) %>%
  left_join(nOPV2_Tracker) %>%
  mutate(diff = target_pop - doses,
         diff_perc = 100 * round(diff / target_pop, 2)) %>%
  View()

# Visually compare sias to detailed tracker
polis_pops %>%
  filter(region == "AFRO") %>%
  filter(vaccinetype == "nOPV2") %>%
  group_by(adm0_name, parentactivitycode, vaccinetype) %>%
  summarize(target_pop = sum(target_pop, na.rm=T),
            start_date = min(start_date)) %>%
  arrange(adm0_name, start_date) %>%
  View()

polis_pops %>%
  filter(region == "AFRO") %>%
  filter(vaccinetype == "nOPV2") %>%
  mutate(year = year(start_date)) %>%
  group_by(adm0_name, year, vaccinetype) %>%
  summarize(target_pop = sum(target_pop, na.rm=T),
            start_date = min(start_date)) %>%
  arrange(adm0_name, start_date) %>%
  View()

polis_pops %>% 
  filter(adm0_name == "ETHIOPIA") %>%
  filter(vaccinetype == "nOPV2") %>%
  select(parentactivitycode) %>% unique()

#### Add pop to missing campaigns ####
# eth_rows <- polis_pops %>% filter(parentactivitycode == "ETH-2021-003") %>% nrow()
# polis_pops[polis_pops$parentactivitycode == "ETH-2021-003", "target_pop"] <- 16712725/eth_rows

# !!!!!!! MANUAL !!!!!! Distribute 12m to Indonesia equally across districts
indo_rows <- polis_pops %>% filter(vaccinetype=="nOPV2", adm0_name == "INDONESIA") %>% nrow()
polis_pops[polis_pops$adm0_name == "INDONESIA" & polis_pops$vaccinetype=="nOPV2", "target_pop"] <- 12415310/indo_rows

# !!!!!!! MANUAL !!!!!! Distribute 14.5m to Philippines equally across districts (based on total post-switch use)
phl_rows <- polis_pops %>% filter(vaccinetype=="mOPV2", adm0_name == "PHILIPPINES") %>% nrow()
polis_pops[polis_pops$adm0_name == "PHILIPPINES" & polis_pops$vaccinetype=="mOPV2", "target_pop"] <- 14535328/phl_rows

# !!!!!!! MANUAL !!!!!! Distribute 1.5m to Malaysia equally across districts (based on total post-switch use)
mys_rows <- polis_pops %>% filter(vaccinetype=="mOPV2", adm0_name == "MALAYSIA") %>% nrow()
polis_pops[polis_pops$adm0_name == "MALAYSIA" & polis_pops$vaccinetype=="mOPV2", "target_pop"] <- 1549306/mys_rows

#!!!!!!!! MANUAL !!!!!! KEN-2024-001 missing target population. Confirmed it was an R3 in Garissa, Mandera, Wajir. Use pops from KEN-2023-002
polis_pops[polis_pops$parentactivitycode == "KEN-2024-001", "target_pop"] <- 
  polis_pops[polis_pops$parentactivitycode == "KEN-2023-002" & polis_pops$adm1_name %in% c("GARISSA", "MANDERA", "WAJIR"), "target_pop"]


# Focus on just EUL era
# polis_pops <- polis_pops %>% 
#   # filter(start_date < today()) %>%
#   filter(start_date <= "2023-12-31") 
# View(polis_pops %>% filter(vaccinetype == "nOPV2", status == "Planned"))

#### Gut check nOPV2 usage ####
polis_pops %>% 
  filter(vaccinetype == "nOPV2") %>% 
  summarize(pop = sum(target_pop, na.rm=T))
polis_pops %>% 
  filter(vaccinetype == "nOPV2") %>%
  filter(start_date < "2024-01-01") %>%
  summarize(pop = sum(target_pop, na.rm=T))

# vaccine usage post-switch
polis_pops %>% 
  filter(start_date > "2016-05-01") %>%
  group_by(source) %>%
  summarize(pop = sum(target_pop, na.rm=T))

polis_pops %>%
  filter(start_date > "2016-05-01") %>%
  # filter(start_date > "2021-08-01") %>%
  group_by(source, region) %>%
  summarize(pop = sum(target_pop, na.rm=T))

polis_pops %>%
  filter(start_date > "2016-05-01") %>%
  group_by(vaccinetype) %>%
  # group_by(region) %>%
  summarize(pop = sum(target_pop, na.rm=T))

polis_pops %>%
  filter(start_date > "2016-05-01", adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO", "PAKISTAN", "EGYPT")) %>%
  group_by(source, region, adm0_name) %>%
  summarize(pop = sum(target_pop, na.rm=T))

# Compare campaign size by vaccine type (For paper)
polis_pops %>% 
  group_by(source, parentactivitycode) %>%
  # group_by(source, vaccinetype) %>%
  # group_by(source, vaccinetype, region) %>% 
  # filter(region %in% c("AFRO", "EMRO", "EURO")) %>%
  # filter(region %in% c("AFRO")) %>%
  filter(africa %in% c("Africa")) %>%
  filter(start_date >= "2016-05-01") %>%
  summarize(sum = sum(target_pop, na.rm=T)) %>%
  ungroup() %>%
  group_by(source) %>%
  # group_by(vaccinetype, region) %>%
  # group_by(vaccinetype) %>%
  summarize(count = n(),
            min = min(sum),
            q1 = quantile(sum, 0.25),
            median = median(sum),
            mean = mean(sum),
            q3 = quantile(sum, 0.75),
            max = max(sum))

temp <- 
  polis_pops %>%
  # filter(region %in% c("AFRO")) %>%
  filter(africa %in% c("Africa")) %>%
  filter(start_date >= "2016-05-01") %>%
  group_by(source, parentactivitycode) %>%
  summarize(sum = sum(target_pop, na.rm=T)) 
wilcox.test(temp$sum ~ temp$source, paired = FALSE,
            exact = FALSE)

# Test for normality
temp %>%
  group_by(source) %>%
  summarise(`W Stat` = shapiro.test(sum)$statistic,
            p.value = shapiro.test(sum)$p.value)

sia_target_box <- polis_pops %>% 
  group_by(source, parentactivitycode) %>% 
  # filter(region %in% c("AFRO", "EMRO", "EURO")) %>%
  filter(africa %in% c("Africa")) %>%
  filter(start_date >= "2016-05-01") %>%
  summarize(sum = sum(target_pop, na.rm=T)) %>%
  filter(sum > 0) %>%
  ungroup() %>%
  ggplot() +
    geom_boxplot(aes(x = source, y = sum/1e6, color = source)) +
    scale_y_log10(name = "SIA Target Population (Million)", 
                  breaks = 10^(2-seq(1,5)),
                  minor_breaks = .5*10^(2-seq(1,5))) +
  xlab("") +
  theme_bw() +
  guides(color="none")
ggsave(plot = sia_target_box, "figures/target pop.png", device = "png", units = "in", width = 3, height = 3)

# Simplify fields
polis_pops_full <- polis_pops
polis_pops <- polis_pops %>% 
  filter(start_date >= "2016-05-01") %>%
  select(adm0_name, adm1_name, adm2_name, parentactivitycode,  childactivitycode, 
                                    start_date, vaccinetype, GUID,
                                    target_pop, period, quarter, source, region, africa)

sias_figure <- polis_pops %>% group_by(quarter, source) %>% summarize(target_pop = sum(target_pop, na.rm=T))
ggplot(sias_figure, aes(x = quarter, y = target_pop, color = source)) +  geom_line()

# Plot cumulative number of doses (For paper)
fig_cum_doses <- 
  polis_pops %>%
  filter(!is.na(target_pop)) %>%
  # filter(region %in% c("AFRO")) %>%
  filter(africa %in% c("Africa")) %>%
  ungroup() %>%
  group_by(source, period) %>%
  reframe(period = unique(period),
            doses = sum(target_pop)) %>%
  group_by(source) %>%
  arrange(period) %>%
  reframe(period = unique(period),
            doses_cumsum = cumsum(doses)) %>%
  ggplot() +
    geom_vline(xintercept = 2024.50, alpha = 0.25, size = 1) +
    geom_vline(xintercept = 2021.167, alpha = 0.25, color = "red", size = 1) +
    geom_line(aes(x = period, y = doses_cumsum/1e6, color = source), size = 1) +
    theme_bw() +
    scale_x_continuous(limits = c(2016.0, 2027), breaks = seq(2016, 2027, 2), name = "") +
    scale_y_continuous(name = "Cumulative Doses\nin Africa (Million)") +
    scale_color_discrete(name = "Vaccine Type") +
    force_panelsizes(rows = unit(2, "in"),
                     cols = unit(5, "in"))
ggsave(plot = fig_cum_doses, "figures/Cumulative Doses Africa.png", device = "png", units = "in", width = 7, height = 3)

anim <- polis_pops %>%
  filter(!is.na(target_pop)) %>%
  # filter(region %in% c("AFRO")) %>%
  filter(africa %in% c("Africa")) %>%
  ungroup() %>%
  group_by(source, period) %>%
  reframe(period = unique(period),
          doses = sum(target_pop)) %>%
  group_by(source) %>%
  arrange(period) %>%
  reframe(period = unique(period),
          doses_cumsum = cumsum(doses)) %>%
  ggplot() +
  geom_vline(xintercept = 2024.50, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, alpha = 0.25, color = "red", size = 1) +
  geom_line(aes(x = period, y = doses_cumsum/1e6, color = source), size = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(2016.0, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_y_continuous(name = "Cumulative Doses\nin Africa (Million)") +
  scale_color_discrete(name = "Vaccine Type") +
  # force_panelsizes(rows = unit(2, "in"),
  #                  cols = unit(5, "in")) +
  transition_reveal(period)
anim_save("figures/Cumulative Doses Africa.gif", anim)

# Plot quarterly number of doses (for paper)
polis_pops %>%
  filter(!is.na(target_pop)) %>%
  ungroup() %>%
  # mutate(year = floor(period)) %>%
  group_by(source, quarter, africa) %>%
  reframe(africa = unique(africa),
          quarter = unique(quarter),
          doses = sum(target_pop)) %>%
  ggplot() +
  geom_col(aes(x = quarter, y = doses/1e6, fill = africa), size = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(2016.0, 2024.5), breaks = seq(2016, 2024, 2), name = "") +
  scale_y_continuous(name = "Quarterly Doses\n(Million)") +
  scale_fill_discrete(name = "") +
  facet_grid(source~.) 
ggsave("figures/Quarterly Doses.png", device = "png", units = "in", width = 7, height = 3)

#### POLIS Virus data ####
token = readLines("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/data_local/token.txt")[1]
viruses_raw = get_polis_virus(token = token, min_date = "2016-05-01", virus_id = 4)

viruses = viruses_raw %>%
  select(id, epid, virus_date, surveillance_type_name, region_who_code,
         admin0name, country_iso3code, admin1name, admin1guid, admin2name, admin2guid, 
         virus_type_name, vdpv_classification_name, vdpv_emergence_group_name, 
         vdpv_nt_changes_from_sabin, surveillance_type_name, po_ns_seq_date) %>%
  mutate(virus_date = ymd(as.Date(virus_date)),
         id = as.numeric(id),
         dot_name = paste(admin0name, admin1name, sep = ":"),
         dot_year_month = tolower(paste(year(virus_date), month(virus_date), sep = ":")))

viruses$africa <- Func_africa(viruses$region_who_code, viruses$admin0name)
viruses %>% group_by(admin0name) %>% summarize(africa = unique(africa)) %>%View()

# Check emergence group names
novel_emergences <- c("RDC-SKV-1", "RDC-TAN-2", "RDC-KOR-1",
                      "CAF-KEM-1", "NIE-KBS-1", "RDC-HKA-2",
                      "CAF-BNG-3", "BOT-FRA-1", "EGY-NOR-1",
                      "CAE-EXT-1", "ZIM-HRE-1", "NIE-KTS-1",
                      "MOZ-MAN-1", "RSS-WEQ-1", "ANG-LNO-3",
                      "ETH-TIG-1", "RSS-JON-1", "RDC-TSH-2" )
novel_emergences %in% c(viruses$vdpv_emergence_group_name %>% unique())

viruses %>% filter(vdpv_emergence_group_name %in% novel_emergences) %>% View()

# cVDPV2 emergences
viruses = viruses %>% 
  filter(virus_type_name == "cVDPV2") %>%
  group_by(vdpv_emergence_group_name) %>%
  mutate(index_date = min(virus_date, na.rm=T),
         index_isolate = (virus_date == index_date))

# Identify repeats
repeats = viruses %>%
  group_by(vdpv_emergence_group_name) %>% 
  filter(sum(index_isolate)>1) %>% 
  select(vdpv_emergence_group_name) %>% 
  summarise(unique(vdpv_emergence_group_name))

for (i in unlist(repeats[,1])){
  cat(i, sep = "\n")
  temp <- NA
  temp <- viruses[viruses$vdpv_emergence_group_name %in% i & viruses$index_isolate == TRUE,"index_isolate"]
  viruses[viruses$vdpv_emergence_group_name %in% i & viruses$index_isolate == TRUE,"index_isolate"] <- c(TRUE, rep(FALSE, nrow(temp)-1)) #Arbitrarily set first to True
}

# Estimate date of seeding 9 changes per year +2 free
viruses$vdpv_nt_changes_from_sabin <- as.numeric(viruses$vdpv_nt_changes_from_sabin)
viruses <- viruses %>% mutate(seeding_date = virus_date - (vdpv_nt_changes_from_sabin-2)/9 * 365.25)

# Assign nOPV2-related viruses
viruses$source = "Sabin2"
sort(unique(viruses$vdpv_emergence_group_name))
viruses[viruses$vdpv_emergence_group_name %in% novel_emergences, "source"] <- "nOPV2"

# Summary of emergences post-switch
viruses %>% 
  filter(index_isolate == TRUE, seeding_date > "2016-05-01") %>%
  mutate(admin0name = str_replace(admin0name, "DEMOCRATIC REPUBLIC OF THE CONGO", "DRC")) %>%
  mutate(admin0name = str_replace(admin0name, "CENTRAL AFRICAN REPUBLIC", "CAR")) %>%
  mutate(admin0name = str_replace(admin0name, "THE UNITED KINGDOM", "UK")) %>%
  mutate(region_who_code_AFRO = if_else(region_who_code %in% c("AFRO"), "AFRO", "Other")) %>%
  ungroup() %>% 
  group_by(admin0name, source, region_who_code_AFRO) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = reorder(admin0name, -count), y = count, fill = source)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) + coord_flip() +
    facet_grid(region_who_code_AFRO~., scales = "free_y") +
    theme_bw() +
    xlab("") + ylab("Post-Switch cVDPV2 Emergences")
ggsave("figures/emergences by country.png", device = "png", units = "in", width = 5, height = 7)

# Calculate period and quarter
viruses$week <- year(viruses$virus_date) + (week(viruses$virus_date)-1)/52
viruses$period <- round(year(viruses$virus_date) + (month(viruses$virus_date)-1)/12, 3)
viruses$quarter <- round(year(viruses$virus_date) + (quarter(viruses$virus_date)-1)/4, 2)

viruses_count_week <- viruses %>% 
  filter(seeding_date > "2016-04-01") %>%
  group_by(week, source) %>% 
  summarize(emergences = sum(index_isolate==TRUE))

viruses_count_period <- viruses %>% 
  filter(seeding_date > "2016-04-01") %>%
  group_by(period, source) %>% 
  summarize(emergences = sum(index_isolate==TRUE))

viruses_count_quarter <- viruses %>% 
  group_by(quarter, source) %>%
  summarize(emergences = sum(index_isolate==TRUE))

viruses_count_period_region <- viruses %>% 
  filter(seeding_date > "2016-04-01") %>%
  group_by(period, source, region_who_code) %>% 
  summarize(emergences = sum(index_isolate==TRUE))

viruses_count_period_africa <- viruses %>% 
  filter(seeding_date > "2016-04-01") %>%
  group_by(period, source, africa) %>% 
  summarize(emergences = sum(index_isolate==TRUE))

# Epi curve
viruses %>% filter(virus_type_name == "cVDPV2", period > 2021,
                   surveillance_type_name == "AFP") %>% 
  mutate(country_group = ifelse(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO", "DRC", 
                                ifelse(admin0name == "NIGERIA", "Nigeria", "Other"))) %>%
  group_by(country_group, period, ) %>% summarize(count = n()) %>%
  ggplot(aes(x = period, y = count, fill = country_group)) + geom_col() +
  ylab("Monthly cVDPV2 Cases")

sias_figure <- left_join(sias_figure, viruses_count_quarter, by = c("quarter", "source"))

#### Immunity mapper data U5 ####
# Run immunity_calc_eag.R from polio-immunity-mapping project to produce new U5 immunity estimates to use here
immunity <- readRDS("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/results/immunity_age_groups_0.8coverage.rds")
immunity$period <- as.numeric(immunity$period)
immunity$immunity_u5 <- immunity$immunity_6_59
immunity_u5_data <- immunity %>% filter(serotype == 'p2', period >= 2014) %>%
  select(guid, serotype, period, immunity_u5)
names(immunity_u5_data)[1] <- "GUID"

immunity_u5_data$period <- as.numeric(immunity_u5_data$period)
immunity_u5_data$year <- floor(immunity_u5_data$period)
immunity_u5_data$month <- (immunity_u5_data$period - immunity_u5_data$year)*12
immunity_u5_data$week <-  immunity_u5_data$year + (week(ymd(paste(immunity_u5_data$year, round(immunity_u5_data$month)+1, "01"))))/52

immunity_u5_data <- immunity_u5_data %>% filter(period >= 2016, period < 2025)

# Add admin1 and admin2 to immunity_u5_data
immunity_u5_data <- immunity_u5_data %>% 
  left_join(polis_pops %>% select(GUID, adm0_name, adm1_name,adm2_name, target_pop) %>%
              group_by(GUID) %>% filter(row_number() == 1), 
            by = "GUID")

write.csv(immunity_u5_data, "immunity_6_59mo.csv")

#### Combine SIA and immunity datasets ####
# Join data on GUID, and period, noting that immunity from campaign doesn't impact that period yet
data <- polis_pops %>% 
  left_join(immunity_u5_data %>% select(GUID, period, week, immunity_u5), by = c("GUID", "period"))

data$population_total <- data$target_pop / 0.17

data %>% group_by(vaccinetype) %>% summarize(target_pop = sum(target_pop, na.rm=T))

# add region # seems unnecessary
# unique(data$adm0_name)
# data$Region <- sapply(data$adm0_name, Func_region)

# See missingness
# data %>% filter(is.na(immunity_u5)) %>% View()
# missing immunity data from Indonesia, Phillipines, and Malaysia. Exclude from analysis

polis_pops %>% filter(vaccinetype == "nOPV2") %>% summarize(sum = sum(target_pop, na.rm = T))
data %>% group_by(adm1_name, period) %>% filter(vaccinetype == "nOPV2") %>% 
  summarize(pop = sum(target_pop, na.rm = T)) %>% ungroup() %>% summarize(sum = sum(pop))

# Convert in province-level campaigns (sum across Admin2)
# When the period and the admin1 name are the same, then sum the population and create a pop-weighted-immunity estimate for the province
data_province <- data %>% group_by(adm1_name, period, vaccinetype) %>%
  reframe(target_pop_sum = sum(target_pop, na.rm = T), 
            population_total_sum = sum(population_total, na.rm = T),
            immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
            vaccinetype = unique(vaccinetype),
            # parentactivitycode = unique(parentactivitycode), # Causing an issue with the count
            region = unique(region),
            africa = unique(africa),
            adm0_name = unique(adm0_name),
            source = unique(source),
            start_date = min(start_date),
            quarter = unique(quarter),
            week = unique(week)) %>% ungroup()
sum(data_province[data_province$vaccinetype == "nOPV2", "target_pop_sum"], na.rm=T)
sum(data_province[data_province$vaccinetype == "nOPV2" & data_province$region == "AFRO", "target_pop_sum"], na.rm=T)
data_province <- data_province %>% filter(is.na(adm1_name) == F)

data_national <- data %>% group_by(adm0_name, period) %>%
  reframe(target_pop_sum = sum(target_pop), 
            population_total_sum = sum(population_total),
            immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
            vaccinetype = unique(vaccinetype),
            region = unique(region),
            africa = unique(africa),
            source = unique(source),
            quarter = unique(quarter))
sum(data_national[data_national$vaccinetype == "nOPV2", "target_pop_sum"], na.rm=T)
sum(data_national[data_national$vaccinetype == "mOPV2", "target_pop_sum"], na.rm=T)
sum(data_national[data_national$vaccinetype == "tOPV", "target_pop_sum"], na.rm=T)

data_quarter <- data %>% group_by(adm0_name, quarter) %>%
  reframe(target_pop_sum = sum(target_pop), 
            population_total_sum = sum(population_total),
            immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
            immunity_weighted_summary = weighted.mean(immunity_weighted, population_total_sum, na.rm=TRUE),
            vaccinetype = unique(vaccinetype),
            region = unique(region),            
            africa = unique(africa),
            source = unique(source),
            quarter = unique(quarter))
sias_figure <- left_join(sias_figure, data_quarter, by = c("quarter", "source"))

# gut check figures
data_province %>% 
  filter(adm0_name == "NIGERIA") %>%
  ggplot(aes(x = period, y = immunity_weighted, color = adm1_name)) +
  geom_point() +
  geom_line()

#### Visualize Pre-campaign immunity ####
# quarterly weighted average of pre-campaign immunity
ggplot(sias_figure, aes(x = quarter, fill = source, color = source)) + 
  geom_line(aes(y = immunity_weighted_summary), size = 1)

# quarterly weighted average of pre-campaign immunity, and each provincial SIA by size
ggplot(data_province, aes(x = period, color = source, size = target_pop_sum)) + 
  geom_point(aes(y = immunity_weighted)) +
  geom_line(data = sias_figure, aes(x = quarter, color = source, y = immunity_weighted), size = 1) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  scale_size(name = "U5 Provincial Target Pop") +
  theme_bw() +
  ylab("Estimated Pre-Campaign Type-2 Immunity")

ggplot(data_province %>% filter(!(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO"))), 
       aes(x = period, color = source)) + 
  geom_point(aes(y = immunity_weighted)) +
  geom_line(data = sias_figure, aes(x = quarter, color = source,y = immunity_weighted_summary), size = 1) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  scale_size(name = "U5 Provincial Target Pop") +
  theme_bw() +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Estimated Pre-Campaign Type-2 Immunity")

# Spotlight on provinces with <0.5 immunity
ggplot(data_province %>% filter(immunity_weighted < 0.5), 
       aes(x = period,y = immunity_weighted, color = source)) + 
  geom_point() + geom_text(aes(label = adm1_name)) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  scale_size(name = "U5 Provincial Target Pop") +
  theme_bw() +
  ylab("Estimated Pre-Campaign Type-2 Immunity")

# summary of pre-campaign immunity by vaccine type for provinces
data_province <- data_province %>% group_by(source) %>% mutate(pop_weight = target_pop_sum / sum(target_pop_sum, na.rm=T))

ggplot(data_province %>% 
         # filter(region %in% c("AFRO", "EMRO", "EURO")), 
         # filter(region %in% c("AFRO")),
         filter(africa %in% c("Africa")),
         aes(x = source, color = source, y = immunity_weighted)) +
  geom_violin(aes(weights = pop_weight), size = 0.75) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  theme(legend.position = "none") +
  # scale_x_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")+
  ylab("Estimated Pre-Campaign\nType-2 Immunity")

data_province %>% 
  filter(region %in% c("AFRO")) %>%
  ggplot(aes(x = source, color = source, y = weighted.mean(immunity_weighted, w = pop_weight))) +
    geom_boxplot() +
    theme(legend.position = "none") +
    ylab("Estimated Pre-Campaign\nType-2 Immunity")

ggplot(data_province %>% filter(!(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO"))), 
       aes(x = source, color = source, y = immunity_weighted)) +
  geom_violin(aes(weights = pop_weight), size = 0.75) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  theme(legend.position = "none") +
  # scale_x_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")+
  scale_y_continuous(limits = c(0,1)) +
  ylab("Estimated Pre-Campaign\nType-2 Immunity")

# summary statistics on pre-campaign immunity by vaccine type at province level
data_province %>%
  filter(is.na(immunity_weighted)==F) %>%
  # filter(!(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO"))) %>%
  group_by(source) %>%
  summarize(min = min(immunity_weighted),
            q1 = quantile(immunity_weighted, 0.25),
            median = median(immunity_weighted),
            mean = weighted.mean(immunity_weighted, target_pop_sum, na.rm=T),
            q3 = quantile(immunity_weighted, 0.75),
            max = max(immunity_weighted))

# Summary of pre-campaign immunity by vaccine type at the SIA level (For paper)
data_sia <- data %>% group_by(parentactivitycode, period) %>%
  reframe(target_pop_sum = sum(target_pop), 
          population_total_sum = sum(population_total),
          immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
          vaccinetype = unique(vaccinetype),
          region = unique(region),
          africa = unique(africa),
          adm0_name = unique(adm0_name),
          source = unique(source),
          start_date = min(start_date),
          quarter = unique(quarter),
          week = unique(week))
sum(data_sia[data_sia$vaccinetype == "nOPV2", "target_pop_sum"], na.rm=T)

temp <- data_sia %>%
  filter(is.na(immunity_weighted)==F) %>%
  # filter(region %in% c("AFRO")) %>%
  filter(africa %in% c("Africa")) %>%
  filter(start_date >= "2016-05-01")

temp %>% 
  group_by(source) %>%
  summarize(min = min(immunity_weighted),
            q1 = quantile(immunity_weighted, 0.25),
            median = median(immunity_weighted),
            mean = weighted.mean(immunity_weighted, target_pop_sum, na.rm=T),
            q3 = quantile(immunity_weighted, 0.75),
            max = max(immunity_weighted))

wilcox.test(temp$immunity_weighted ~ temp$source,
            exact = FALSE)

# Immunity boxplot (for paper)
sia_immunity_box <- temp %>%
  ggplot(aes(x = source, y = immunity_weighted, color = source)) +
  # geom_violin() +
  geom_boxplot() +
  # geom_jitter(alpha = 0.5) +
  theme_bw() +
  guides(color = "none") +
  ylab("Pre-Campaign Type-2 Immunity") +
  xlab("")
ggsave(plot = sia_immunity_box, "figures/immunity.png", device = "png", units = "in", width = 3, height = 3)

# SIA target and immunity boxplots (For paper)
plot_layout(ncol = 1, nrow = 2, sia_target_box / sia_immunity_box)
ggsave("figures/SIA target and immunity.png", device = "png", units = "in", width = 3, height = 6)

#### Define which provinces have ES ####
# make a list of provinces with an ES result within past year. Manually cross-checked with ES.world
viruses_es = get_polis_virus(token = token, min_date = today()-365)
es_provinces = viruses_es %>% filter(surveillance_type_name == "Environmental") %>%
  mutate(dot.name = paste(admin0name, admin1name, sep = ".")) %>%
  select(dot.name) %>% unique() %>%
  mutate(adm0_name = NA,
         adm1_name = NA)

for (i in 1:nrow(es_provinces)){
  es_provinces$adm0_name[i] <- strsplit(es_provinces$dot.name[i], "[.]")[[1]][1]
  es_provinces$adm1_name[i] <- strsplit(es_provinces$dot.name[i], "[.]")[[1]][2]
}

es_provinces = es_provinces %>%
  mutate(ES = 1) %>%
  select(adm0_name, adm1_name, ES)

# make a field in data_province for ES site TRUE/FALSE
data_province <- left_join(data_province, es_provinces, by = c("adm0_name", "adm1_name"))
data_province <- data_province %>% mutate(ES = !is.na(ES))

#### Compare size of campaigns by nOPV2 and Sabin2 ####
data_province %>% 
  filter(africa %in% c("Africa")) %>%
  ggplot(aes(x = period, y = target_pop_sum, color = source)) + 
  geom_point() +
  geom_smooth() +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  # facet_grid(Region~.) +
  scale_y_continuous(limits = c(0, 7e6), name = "Provincial Target Population")

data_province %>%
  filter(is.na(target_pop_sum)==F) %>%
  filter(africa %in% c("Africa")) %>%
  group_by(source) %>%
  summarize(doses = sum(target_pop_sum),
            min = min(target_pop_sum),
            q1 = quantile(target_pop_sum, 0.25),
            median = median(target_pop_sum),
            mean = mean(target_pop_sum),
            q3 = quantile(target_pop_sum, 0.75),
            max = max(target_pop_sum))

#### Compare locations of nOPV2 and Sabin2 use ####
temp2 <- polis_pops %>%
  filter(target_pop > 0) %>% # Note that philippines, malaysia perhaps others have missing population sizes
  filter(period >= 2016.417) %>%
  mutate(region_AFRO = if_else(region %in% c("AFRO"), "AFRO", "Other")) %>%
  mutate(adm0_name = str_replace(adm0_name, "DEMOCRATIC REPUBLIC OF THE CONGO", "DRC")) %>%
  mutate(adm0_name = str_replace(adm0_name, "CENTRAL AFRICAN REPUBLIC", "CAR")) %>%
  mutate(adm0_name = str_replace(adm0_name, "UNITED REPUBLIC OF TANZANIA", "TANZANIA")) %>%
  mutate(adm0_name = str_replace(adm0_name, "IRAN (ISLAMIC REPUBLIC OF)", "IRAN")) %>%
  mutate(adm0_name = str_replace(adm0_name, "SYRIAN ARAB REPUBIC", "SYRIA")) %>%
  group_by(adm0_name, source) %>%
  summarize(target_pop = sum(target_pop, na.rm=T),
            africa = unique(africa),
            region_AFRO = unique(region_AFRO)) 

ggplot(temp2, aes(fill = source, y = reorder(adm0_name, target_pop), x = target_pop/1e6)) +
  geom_col(position = position_dodge2(preserve = "single"), orientation = "y") +
  ylab(element_blank()) +
  xlab("Target U5 Population (Million)") +
  facet_grid(region_AFRO~., scales = "free") +
  scale_fill_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  theme_bw()

# National dosage usage (for paper)
ggplot(temp2 %>% filter(africa == "Africa"), aes(color = source, y = reorder(adm0_name, target_pop), x = target_pop/1e6)) +
  geom_point() +
  ylab(element_blank()) +
  scale_x_log10(name = "Target U5 Population (Million)", 
                breaks = 10^(seq(0,7,1)),
                minor_breaks = .5*10^(seq(0,7,1))) +
  # facet_grid(africa~., scales = "free") +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  theme_bw()
ggsave("figures/doses by country_africa.png", device = "png", units = "in", width = 6, height = 6)

#### Crude per-dose emergence expectation analysis ####

# Regional dosage usage
polis_pops %>%
  group_by(region) %>%
  filter(target_pop > 0) %>%
  filter(period >= 2016.417) %>%
  # filter(period < 2024) %>%
  # filter(period < 2024.417) %>%
  group_by(region, source) %>%
  summarize(target_pop = sum(target_pop, na.rm=T)) 

# Africa/Elsewhere dosage usage (for paper table 1)
polis_pops %>%
  group_by(africa) %>%
  filter(target_pop > 0) %>%
  filter(period >= 2016.417) %>%
  # filter(period < 2024) %>%
  # filter(period < 2024.417) %>%
  group_by(africa, source) %>%
  summarize(target_pop = sum(target_pop, na.rm=T)) 

# View(polis_pops %>% filter(region %in% c("WPRO"), start_date >= "2016-05-01"))
View(polis_pops %>% filter(region %in% c("SEARO"), start_date >= "2016-05-01"))

polis_pops %>%
  filter(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO",
                          "PAKISTAN", "EGYPT")) %>%
  group_by(adm0_name) %>%
  filter(target_pop > 0) %>%
  filter(period >= 2016.417) %>%
  group_by(adm0_name, source) %>%
  summarize(target_pop = sum(target_pop, na.rm=T)) 

# Description of emergences (for table 1)
viruses %>%
  filter(index_isolate == TRUE, virus_date > "2016-05-01") %>% 
  mutate(seeded_post_switch = if_else(seeding_date > "2016-05-01","post", "pre")) %>%
  group_by(seeded_post_switch, source) %>%
  reframe(count = n())
viruses %>%
  filter(index_isolate == TRUE, seeding_date > "2016-05-01") %>% 
  group_by(region_who_code, source) %>%
  reframe(count = n())
viruses %>%
  filter(index_isolate == TRUE, seeding_date > "2016-05-01") %>% 
  group_by(africa, source) %>%
  reframe(count = n())
viruses %>%
  filter(index_isolate == TRUE, seeding_date > "2016-05-01") %>% 
  group_by(region_who_code, source, admin0name) %>%
  reframe(count = n()) %>% View()

# Calculate crude regional per-dose emergence rate
for (i in c("AFRO", "EMRO", "EURO", "SEARO", "WPRO")){
  sabin2_doses <- polis_pops %>% 
    filter(region == i) %>%
    filter(source == "Sabin2") %>%
    filter(period >= 2016.417) %>%
    summarize(target_pop = sum(target_pop, na.rm=T)) %>% 
    as.numeric()
  
  nOPV2_doses <- polis_pops %>% 
    filter(region == i) %>%
    filter(source == "nOPV2") %>%
    filter(period >= 2016.417) %>%
    summarize(target_pop = sum(target_pop, na.rm=T)) %>% 
    as.numeric()
  
  sabin2_emergences <- viruses %>%
    filter(source == "Sabin2") %>%
    filter(region_who_code == i) %>%
    filter(index_isolate == TRUE, seeding_date > "2016-05-01") %>% 
    summarize(count = n()) %>%
    summarize(sum = sum(count)) %>% 
    as.numeric()
  
  nOPV2_emergences <- viruses %>%
    filter(source == "nOPV2") %>%
    filter(region_who_code == i) %>%
    filter(index_isolate == TRUE, seeding_date > "2016-05-01") %>% 
    summarize(count = n()) %>%
    summarize(sum = sum(count)) %>% 
    as.numeric()
  
  cat("In ", i, "\n",
      sabin2_doses, " Sabin2 doses\n", 
      nOPV2_doses, " nOPV2 doses\n",
      sabin2_emergences, " Sabin2 emergences\n",
      nOPV2_emergences, " nOPV2 emergences\n",
      sabin2_emergences/sabin2_doses*10000000, " Sabin2 emergences per million doses\n",
      nOPV2_emergences/nOPV2_doses*10000000, " nOPV2 emergences per million doses\n\n")
}

# Manual check of data (SEARO)
View(polis_pops %>% filter(region %in% c("SEARO"), start_date >= "2016-05-01"))
viruses %>% filter(vdpv_emergence_group_name == "INO-ACE-1") %>% View()
# Post-switch SIAs in Indonesia were in response to the emergence, not the source.

# Manual check of data (EURO)
View(polis_pops %>% filter(region %in% c("EURO"), start_date >= "2016-05-01"))
viruses %>% filter(vdpv_emergence_group_name == "IUUC-2022") %>% View()
# Unsure if originally UK or Israel, and unsure originating SIA (heard Pakistan at some point?)

# Manual check of data (WPRO)
View(polis_pops %>% filter(region %in% c("WPRO"), start_date >= "2016-05-01"))
viruses %>% filter(vdpv_emergence_group_name == "CHN-SIC-1 ") %>% View()
# SIAs in Malaysia and Philippines in 2019-2020 inconsistent with emergence in China in 2018 with 13 nt changes. Possible pre-switch? 

# Estimate global crude lambda 
sabin2_emergences = viruses %>%
  filter(source == "Sabin2") %>%
  filter(index_isolate == TRUE, seeding_date > "2016-05-01") %>% 
  summarize(count = n()) %>%
  summarize(sum = sum(count)) %>% 
  as.numeric()

sabin2_doses <- polis_pops %>% 
  filter(source == "Sabin2") %>%
  filter(period >= 2016.417) %>%
  summarize(target_pop = sum(target_pop, na.rm=T)) %>% 
  as.numeric()

nOPV2_doses <- polis_pops %>% 
  filter(source == "nOPV2") %>%
  filter(period >= 2016.417) %>%
  summarize(target_pop = sum(target_pop, na.rm=T)) %>% 
  as.numeric()

# Estimated emergences based on nOPV2 doses 
crude_lambda = sabin2_emergences/sabin2_doses
qpois(c(0.025, 0.5, 0.975), crude_lambda * nOPV2_doses)


#### AFR: Function to calculate expected count of emergence from campaign ####
theta_size_mean  <- 0.795 
theta_size_lower <- 0.503
theta_size_upper <- 1.10

theta_u5_mean  <- -1.98
theta_u5_lower <- -3.31
theta_u5_upper <- -0.647

Func_u = function(p,  #p is the total population (not U-5)
                  q,  #q is the immunity in children under 5-years
                  alpha = 2.025*10^-6, #alpha is the scaling factor for total number of observed emergences (see Supplemental section 4 in Gray 2022) 
                  theta_size = theta_size_mean, 
                  theta_u5 = theta_u5_mean,
                  crude = FALSE,
                  lambda = crude_lambda){ 
  if (crude == FALSE){
    alpha * exp((theta_size-1)*log(p) + theta_u5*q) * p
  } else (p*0.17 * lambda)
}
alpha = 2.16*10^-6 # Initial

# testing function
Func_u(1e6, c(.05, .5, .95), alpha) # fixed pop size, increasing immunity, should have decreasing emergences
Func_u(c(1e5, 1e6, 1e7), .5, alpha) # increasing pop size, fixed immunity, should have increasing emergences
Func_u(3e7, c(.2, .4, .6, .8)) # compare with figure 4 in Gray 2023 paper
round(Func_u(1e6, 0.5, alpha) * 1.1, 6) == round(Func_u(1e6, 0.5, alpha*1.1), 6) # Expected number of emergences scales linearly with Alpha
Func_u(1e6/0.17, crude = TRUE)

# testing function
data_province <- data_province %>% 
  mutate(U_mOPV2 = Func_u(p=population_total_sum,
                          q = immunity_weighted,
                          alpha = alpha))
data_national <- data_national %>% mutate(U_mOPV2 = Func_u(p=population_total_sum, q = immunity_weighted, alpha = alpha))

data_province_quarter <- data_province %>% group_by(quarter, source) %>%
  summarize(U_mOPV2_sum = sum(U_mOPV2, na.rm=TRUE),
            target_pop_sum = sum(target_pop_sum, na.rm=TRUE),
            immunity_weighted_summary = weighted.mean(immunity_weighted, population_total_sum, na.rm=TRUE)) %>%
  select(quarter, source, U_mOPV2_sum, immunity_weighted_summary, target_pop_sum)

# note, immunity_weighted_summary is used above. Consider moving up
sias_figure <- left_join(sias_figure, data_province_quarter %>% select("quarter", "source", "U_mOPV2_sum"), by = c("quarter", "source"))

# Assuming poisson distribution of #emergences, calculate total number of expected emergences from probability, based on initial alpha
data_province %>% group_by(vaccinetype) %>% summarize(pop = sum(target_pop_sum, na.rm=T), U_mOPV2 = sum(U_mOPV2, na.rm=T))
data_province %>% group_by(source) %>% summarize(pop = sum(target_pop_sum, na.rm=T), U_mOPV2 = sum(U_mOPV2, na.rm=T))

# Plot number of doses by vaccine type and number of emergences observed and expected
ggplot(sias_figure, aes(x = quarter, fill = source, color = source)) + 
  geom_line(aes(y = target_pop_sum), size = 1) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  theme_bw() +
  geom_col(aes(y = U_mOPV2_sum*1e6), position = position_dodge2(preserve = "single")) +
  geom_point(aes(y = emergences*1e6), color = "black")

# Observed VDPV2 emergences by quarter
# ggplot(sias_figure, aes(x = quarter)) + 
#   geom_col(aes(y = emergences, fill = source)) +
#   ylab("VDPV2 Emergences") +
#   theme_bw()

#### Surveillance/sequencing lag ####
# Import sequence lab fit (derived from IDM RA Model)
# surveillance_sequence_lag.R
seq_lag_fit <- read_csv("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/results/sequence_lag_fit.csv")
max_lag = 365

seq_lag_fit$region <- Func_region(seq_lag_fit$admin0name)

period_midpoints <- round(365/24 * seq(1, 11, 2))

# Convert from days into period using mean during that period (aka month)
seq_lag_fit_period <- seq_lag_fit %>% 
  mutate(date = as.Date("2022-01-01") + day) %>%
  mutate(period = month(date)) %>%
  ungroup %>%
  group_by(admin0name, period) %>%
  mutate(dseq_mean = mean(dseq),
         pseq_mean = mean(pseq)) %>%
  filter(day %in% period_midpoints) %>%
  select(-c("date", "dseq", "pseq")) %>%
  ungroup()

# seq_lag_fit_period %>% filter(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>% View()
# seq_lag_fit_period %>% filter(admin0name == "CENTRAL AFRICAN REPUBLIC") %>% View()

ggplot() +
  geom_line(aes(x=day,y=pseq),colour='black',linetype=2,data=seq_lag_fit %>% filter(region %in% c("AFRO"))) +
  geom_point(aes(x=day,y=pseq_mean),colour='red',data=seq_lag_fit_period %>% filter(region %in% c("AFRO"))) +
  facet_wrap(vars(admin0name))+
  theme(legend.position = 'none') +
  scale_y_continuous('Proportion of Sequences Available',breaks = c(0,0.5,1))+
  scale_x_continuous('Days after collection', limits=c(0,max_lag),oob=scales::squish)
ggsave("figures/seq lag afro.png", device = "png", units = "in", width = 8, height = 6) #(for paper)

#### Estimate time to linkage (ttl) lag from index isolate to confirmatory isolate ####

# Identify confirmatory isolate for each emergence group
viruses = viruses %>% ungroup() %>% group_by(vdpv_emergence_group_name) %>%
  mutate(ttl = as.numeric(virus_date - min(virus_date))) %>% ungroup()
viruses$conf_isolate <- FALSE

for (i in unique(viruses$vdpv_emergence_group_name)){
  if (nrow(viruses[viruses$vdpv_emergence_group_name == i,"ttl"]) > 1){
    min_ttl <- min(as.numeric(viruses[viruses$vdpv_emergence_group_name == i & 
                                        viruses$index_isolate == FALSE, "ttl"][[1]]))
    viruses[viruses$vdpv_emergence_group_name == i & 
              viruses$ttl == min_ttl, "conf_isolate"] <- TRUE
  }
}  

# Identify repeats
repeats_ttl = viruses %>%
  group_by(vdpv_emergence_group_name) %>% 
  filter(sum(conf_isolate)>1) %>% 
  select(vdpv_emergence_group_name) %>% 
  summarise(unique(vdpv_emergence_group_name))

for (i in unlist(repeats_ttl[,1])){
  cat(i, sep = "\n")
  temp <- NA
  temp <- viruses[viruses$vdpv_emergence_group_name %in% i & viruses$conf_isolate == TRUE,"conf_isolate"]
  viruses[viruses$vdpv_emergence_group_name %in% i & viruses$conf_isolate == TRUE,"conf_isolate"] <- c(TRUE, rep(FALSE, nrow(temp)-1)) #Arbitrarily set first to True
}

# Distribution for time-to-linkage
viruses %>% 
  filter(conf_isolate == TRUE) %>%
  # group_by(region_who_code) %>%
  summarize(count = n(),
            min = min(ttl),
            q1 = quantile(ttl, 0.25),
            median = median(ttl),
            mean = mean(ttl),
            q3 = quantile(ttl, 0.75),
            max = max(ttl))

ggplot(viruses %>% filter(conf_isolate == TRUE, region_who_code %in% c("AFRO", "EMRO")), aes(x = ttl, fill = region_who_code)) + 
  geom_histogram() +
  facet_wrap(region_who_code~.) +
  xlab("Days") +
  ylab("Emergence Groups")+
  ggtitle("Time to Linkage (Index to Confirmatory Isolate)")

#estimate distribution with empirical bayes---- (based on methods for sequence lag estimation)
# get initial value from pooled data:
max_lag = 365
xx = viruses %>% filter(conf_isolate == TRUE, ttl < max_lag) %>% pull(ttl)
est = c(log(mean(xx)), log(1))
est = optim(est, function(x){ -sum(dgamma(xx, rate = exp(x[2] - x[1]),shape = exp(x[2]),log=TRUE)) })$par
est

curve(dgamma(x,rate = exp(est[2]-est[1]), shape = exp(est[2])),from=0,to=150)

# fit separately to each country:
fit_eb = viruses %>% filter(conf_isolate == TRUE, ttl < max_lag) %>% 
  group_by(admin0name) %>% 
  mutate(n = n()) %>% ungroup %>% 
  # filter(n >= 2) %>% 
  split(.$admin0name) %>%
  map_df(~{
    par = optim(est, function(x){ -sum(dgamma(.x$ttl, rate = exp(x[2] - x[1]),shape = exp(x[2]),log=TRUE)) })$par
    tibble(lmu = par[1], lshape = par[2])
  }, .id = 'admin0name')

# point estimate of hyper parameters:
lmu_mean = mean(fit_eb$lmu)
lmu_sd = sd(fit_eb$lmu)

lshape_mean = mean(fit_eb$lshape)
lshape_sd = sd(fit_eb$lshape)

#iteratively get posterior mode and update hyper parameter point estimates:
for(i in 1:10){
  fit_eb = viruses %>% filter(conf_isolate == TRUE, ttl < max_lag) %>%
    group_by(admin0name) %>% mutate(n = n()) %>% ungroup %>% 
    split(.$admin0name) %>%
    map_df(~{
      par = optim(est, function(x){ 
        -sum(dgamma(.x$ttl, rate = exp(x[2] - x[1]), shape = exp(x[2]),log=TRUE)) + 
          0.5*((x[1]-lmu_mean)/lmu_sd)^2 + 0.5*((x[2]-lshape_mean)/lshape_sd)^2})$par
      tibble(lmu = par[1], lshape = par[2])
    }, .id = 'admin0name')
  
  lmu_mean = mean(fit_eb$lmu)
  lmu_sd = sd(fit_eb$lmu)
  
  lshape_mean = mean(fit_eb$lshape)
  lshape_sd = sd(fit_eb$lshape)
} 

fit_eb = fit_eb %>% 
  add_row(admin0name = 'other', 'lmu' = lmu_mean, 'lshape'=lshape_mean) 

df_fit_eb = expand_grid(admin0name =  c('other',unique(fit_eb$admin0name)), day=0:200)

df_fit_eb = df_fit_eb %>% 
  left_join(fit_eb) %>%
  mutate(dseq = dgamma(day,rate = exp(lshape)/exp(lmu), shape=exp(lshape)),
         pseq = pgamma(day,rate = exp(lshape)/exp(lmu), shape=exp(lshape)))

#plot
ggplot() + 
  geom_histogram(aes(y=after_stat(density),fill=admin0name,x=ttl),data=viruses %>% filter(conf_isolate == TRUE))+
  geom_line(aes(x=day,y=dseq),colour='black',linetype=2,data=df_fit_eb) + 
  facet_wrap(vars(admin0name),scales = 'free_y')+
  theme(legend.position = 'none') + 
  scale_x_continuous('Days Index to Confirmatory Isolate',limits=c(0,max_lag),oob=scales::squish)
# ggsave(glue::glue('results/seq_lag_fit_type{serotype_use}.png'),width=12,height=10,units='in')

ggplot() + 
  geom_step(aes(colour=admin0name,x=ttl),stat='ecdf',data=viruses %>% filter(conf_isolate, ttl < 365))+
  geom_line(aes(x=day,y=pseq),colour='black',linetype=2,data=df_fit_eb) + 
  facet_wrap(vars(admin0name))+
  theme(legend.position = 'none') + 
  scale_y_continuous('Proportion of cVDPV2 linkages confirmed',breaks = c(0,0.5,1))+
  scale_x_continuous('Days Index to Confirmatory Isolate', limits=c(0,max_lag),oob=scales::squish)
# ggsave(glue::glue('results/seq_lag_fit_ecdf_type{serotype_use}.png'),width=10,height=8,units='in')

# REPEAT, BUT FOR REGIONS INSTEAD OF COUNTRIES
# get initial value from pooled data:
max_lag = 365
xx = viruses %>% filter(conf_isolate == TRUE, ttl < max_lag) %>% pull(ttl)
est = c(log(mean(xx)), log(1))
est = optim(est, function(x){ -sum(dgamma(xx, rate = exp(x[2] - x[1]),shape = exp(x[2]),log=TRUE)) })$par
est

curve(dgamma(x,rate = exp(est[2]-est[1]), shape = exp(est[2])),from=0,to=150)

# fit separately to each region:
fit_eb = viruses %>% filter(conf_isolate == TRUE, ttl < max_lag) %>% 
  group_by(region_who_code) %>% 
  mutate(n = n()) %>% ungroup %>% 
  # filter(n >= 2) %>% 
  split(.$region_who_code) %>%
  map_df(~{
    par = optim(est, function(x){ -sum(dgamma(.x$ttl, rate = exp(x[2] - x[1]),shape = exp(x[2]),log=TRUE)) })$par
    tibble(lmu = par[1], lshape = par[2])
  }, .id = 'region_who_code')

# point estimate of hyper parameters:
lmu_mean = mean(fit_eb$lmu)
lmu_sd = sd(fit_eb$lmu)

lshape_mean = mean(fit_eb$lshape)
lshape_sd = sd(fit_eb$lshape)

#iteratively get posterior mode and update hyper parameter point estimates:
for(i in 1:10){
  fit_eb = viruses %>% filter(conf_isolate == TRUE, ttl < max_lag) %>%
    group_by(region_who_code) %>% mutate(n = n()) %>% ungroup %>% 
    split(.$region_who_code) %>%
    map_df(~{
      par = optim(est, function(x){ 
        -sum(dgamma(.x$ttl, rate = exp(x[2] - x[1]), shape = exp(x[2]),log=TRUE)) + 
          0.5*((x[1]-lmu_mean)/lmu_sd)^2 + 0.5*((x[2]-lshape_mean)/lshape_sd)^2})$par
      tibble(lmu = par[1], lshape = par[2])
    }, .id = 'region_who_code')
  
  lmu_mean = mean(fit_eb$lmu)
  lmu_sd = sd(fit_eb$lmu)
  
  lshape_mean = mean(fit_eb$lshape)
  lshape_sd = sd(fit_eb$lshape)
} 

fit_eb = fit_eb %>% 
  add_row(region_who_code = 'other', 'lmu' = lmu_mean, 'lshape'=lshape_mean) 

df_fit_eb = expand_grid(region_who_code =  c('other',unique(fit_eb$region_who_code)), day=0:200)

df_fit_eb = df_fit_eb %>% 
  left_join(fit_eb) %>%
  mutate(dseq = dgamma(day,rate = exp(lshape)/exp(lmu), shape=exp(lshape)),
         pseq = pgamma(day,rate = exp(lshape)/exp(lmu), shape=exp(lshape)))

#plot
ggplot() + 
  geom_histogram(aes(y=after_stat(density),fill=region_who_code,x=ttl),
                 data=viruses %>% filter(conf_isolate == TRUE, region_who_code %in% c("AFRO", "EMRO")))+
  geom_line(aes(x=day,y=dseq),colour='black',linetype=2,data=df_fit_eb %>% filter(region_who_code %in% c("AFRO", "EMRO"))) + 
  facet_wrap(vars(region_who_code),scales = 'free_y')+
  theme(legend.position = 'none') + 
  scale_x_continuous('Days Index to Confirmatory Isolate',limits=c(0,max_lag),oob=scales::squish)
# ggsave(glue::glue('results/seq_lag_fit_type{serotype_use}.png'),width=12,height=10,units='in')

ggplot() + 
  geom_step(aes(colour=region_who_code,x=ttl),stat='ecdf',
            data=viruses %>% filter(conf_isolate, ttl < 365, region_who_code %in% c("AFRO", "EMRO")))+
  geom_line(aes(x=day,y=pseq),colour='black',linetype=2,data=df_fit_eb %>% filter(region_who_code %in% c("AFRO", "EMRO"))) + 
  facet_wrap(vars(region_who_code))+
  theme(legend.position = 'none') + 
  scale_y_continuous('Proportion of cVDPV2 linkages confirmed',breaks = c(0,0.5,1))+
  scale_x_continuous('Days Index to Confirmatory Isolate', limits=c(0,max_lag),oob=scales::squish)
# ggsave(glue::glue('results/seq_lag_fit_ecdf_type{serotype_use}.png'),width=10,height=8,units='in')

# Convert to months (but they don't fit well...)
df_fit_eb_period <- df_fit_eb %>%
  mutate(period = (month(ymd(20230101) + day)-1)/12) %>%
  group_by(period, region_who_code) %>%
  mutate(dseq = sum(dseq),
         pseq = max(pseq)) %>%
  filter(row_number()==1)%>% 
  ungroup()

# Convert to weeks (months don't fit well...)
df_fit_eb_period <- df_fit_eb %>%
  mutate(period = (week(ymd(20230101) + day)-1)/52) %>%
  group_by(period, region_who_code) %>%
  mutate(dseq = sum(dseq),
         pseq = max(pseq)) %>%
  filter(row_number()==1) %>% 
  ungroup()

ggplot() + 
  geom_step(aes(colour=region_who_code,x=ttl),stat='ecdf',
            data=viruses %>% filter(conf_isolate, ttl < 365, region_who_code %in% c("AFRO", "EMRO")))+
  geom_line(aes(x=day,y=pseq),colour='black',linetype=2,data=df_fit_eb_period %>% filter(region_who_code %in% c("AFRO", "EMRO"))) + 
  facet_wrap(vars(region_who_code))+
  theme(legend.position = 'none') + 
  scale_y_continuous('Proportion of cVDPV2 linkages confirmed',breaks = c(0,0.5,1))+
  scale_x_continuous('Days Index to Confirmatory Isolate', limits=c(0,max_lag),oob=scales::squish)

#### Convolution of TTD, seq_lag ####

# Inputs from time-to-detection distribution: logmean, logsd
m <- 2.460459 #mean log AFP only
s <- 0.2689168 #sd log AFP only
m_es <- 2.16862 #mean log AFP and ES
s_es <- 0.2274121 #sd log AFP and ES

func_afp_only <- function(x) dlnorm(x, meanlog=m, sdlog = s)
func_afp_es <- function(x) dlnorm(x, meanlog=m_es, sdlog = s_es)

#visualize parameters
curve(dlnorm(x, meanlog = m, sdlog = s), from=0, to=36)
curve(dlnorm(x, meanlog = m_es, sdlog = s_es), from=0, to=36)
curve(dlnorm(x, meanlog = m*1.25, sdlog = s), from=0, to=36)

# ggplot
afp_only <- dlnorm(1:48, m, s)
afp_es <- dlnorm(1:48, m_es, s_es)

sum(afp_only)
cumsum(afp_only) # 18 month cumsum is 95.6%
cumsum(afp_es) # 13 month cumsum is 97.2%

ttd_df <- data.frame("Surveillance" = "AFP Only", "Month" = 1:length(afp_only), "Density" = afp_only)
ttd_df <- ttd_df %>% rbind(data.frame("Surveillance" = "AFP & ES", "Month" = 1:length(afp_es), "Density" = afp_es))

# for paper
ggplot() +
  stat_function(fun = func_afp_only, aes(color = "AFP Only")) +
  stat_function(fun = func_afp_es, aes(color = "AFP & ES")) +
  scale_x_continuous(name = "Months", limits = c(0, 36), breaks = seq(0,36,6)) +
  scale_y_continuous(name = "Density") +
  labs(color = "Surveillance") +
  scale_color_manual(values = c("blue", "black")) +
  theme_bw() +
  theme(legend.position = c(.7,.7))
ggsave("figures/time to detection.png", device = "png", units = "in", width = 3, height = 3)

# placeholder for 7 years of post-switch data, plus 4 years of emergence risk forward
max(data_province$period) - min(data_province$period)
months = (7+4) * 12
period = seq(as.Date("2016-04-01"), length = months+1, by = "1 month")
period <- round(year(period) + (month(period)-1)/12, digits=3)

# Function to extend vector with zeros
Func_zeros <- function(vector, max=60){
  if (length(vector) >= max){
    vector[1:max]
  } else {
    c(vector, rep(0, max-length(vector)))
  }
}

# Convolution function
convolve_discrete = function(max_month = 60L, 
                             probs_ttd,
                             probs_ttl, 
                             probs_seq,
                             stdize = T,
                             ttl_include){
  
  
  if(length(probs_ttd) != max_month |
     length(probs_ttl) != max_month |
     length(probs_seq) != max_month){
    stop("Input probabilities don't match output times")
  }
  
  if(stdize){   
    probs_ttd = probs_ttd / sum(probs_ttd)
    probs_ttl = probs_ttl / sum(probs_ttl) 
    probs_seq = probs_seq / sum(probs_seq) # may return Inf if all zero
  }
  
  # if(sum(probs_ttd == 0)) return(tibble(month = 1:max_month, prob = probs_ttd))
  
  if (ttl_include){
    tibble(month_ttd = 1:max_month, probs_ttd = probs_ttd) %>%
      expand_grid(
        tibble(month_ttl = 1:max_month, probs_ttl = probs_ttl)) %>%
      expand_grid(
        tibble(month_seq = 1:max_month, probs_seq = probs_seq)) %>%
      mutate(month = month_ttd + month_ttl + month_seq) %>%
      group_by(month) %>%
      summarize(prob = sum(probs_ttd*probs_ttl*probs_seq)) %>%
      ungroup() %>%
      select("month", "prob")
  } else{
    tibble(month_ttd = 1:max_month, probs_ttd = probs_ttd) %>%
      expand_grid(
        tibble(month_seq = 1:max_month, probs_seq = probs_seq)) %>%
      mutate(month = month_ttd + month_seq) %>%
      group_by(month) %>%
      summarize(prob = sum(probs_ttd*probs_seq)) %>%
      ungroup() %>%
      select("month", "prob")
  }
}

# time-to-convolution
Func_tt_conv <- function(U_mOPV2,
                         adm0_name, adm1_name, region,
                         ES,
                         mean_log = m, sd_log = s, 
                         mean_log_es = m_es, sd_log_es = s_es,
                         max = 60, 
                         output = "U_d",
                         ttl_include = T){
  
  if (ES){
    ttd <- dlnorm(1:48, mean_log_es, sd_log_es) %>% Func_zeros()
  } else{
    ttd <- dlnorm(1:48, mean_log, sd_log) %>% Func_zeros()
  }
  ttd = ttd/sum(ttd) # standardize
  
  if (region %in% df_fit_eb_period$region_who_code){
    ttl <- df_fit_eb_period %>% filter(region_who_code == region) %>% select(dseq)
    ttl <- ttl[[1]] %>% Func_zeros()
  } else {
    cat("Error 1: region_who_code = ", region_who_code)
    stop()
  }
  ttl = ttl/sum(ttl) # standardize
  
  if (adm0_name %in% seq_lag_fit_period$admin0name){
    seq <- seq_lag_fit_period %>% filter(admin0name == adm0_name) %>% select(dseq_mean)
    seq <- seq[[1]] %>% Func_zeros()
  } else {
    seq <- seq_lag_fit_period %>% filter(admin0name == "other") %>% select(dseq_mean)
    seq <- seq[[1]] %>% Func_zeros()
  }
  seq = seq/sum(seq) # standardize
  
  tt_conv = convolve_discrete(60, probs_ttd = ttd, probs_ttl = ttl, probs_seq = seq, ttl_include = ttl_include)
  
  min = min(tt_conv$month)
  if (min > 1){
    temp <- tibble(month = 1:(min-1), prob = 0)
    tt_conv <- rbind(temp, tt_conv)
  }
  
  U_d <- U_mOPV2 * tt_conv$prob
  
  if (output == "U_d"){
    U_d
  } else {
    list("ttd" = ttd,
         "ttl" = ttl,
         "seq" = seq,
         "tt_conv" = tt_conv,
         "U_d" = U_d)
  }
}

# Test plot for a sample province
tt_conv_test <- Func_tt_conv(U_mOPV2 = 0.1, adm0_name = "NIGERIA", adm1_name = "BORNO", region = "AFRO", ES = T, output = "all", ttl_include = T)
sum(tt_conv_test$tt_conv$prob)

ggplot() + 
  geom_line(data = tt_conv_test$tt_conv, aes(x = month, y = prob)) +
  geom_line(data = data.frame(ttd = tt_conv_test$ttd, month = 1:60), aes(x = month, y = ttd), color = "red") +
  geom_line(data = data.frame(ttl = tt_conv_test$ttl, month = 1:60), aes(x = month, y = ttl), color = "blue") +
  geom_line(data = data.frame(seq = tt_conv_test$seq/sum(tt_conv_test$seq), month = 1:60), aes(x = month, y = seq), color = "green") +
  scale_x_continuous(limits = c(0,30))

#### Functions to estimate U_d_i for each day d and campaign i ####
Func_daily_U_conv <- function(data_province, period_input = period){
  period = period_input
  tt_conv_data = data_province %>% 
    filter(is.na(U_mOPV2) == F) %>%
    select(adm0_name, adm1_name, region, period, vaccinetype, source, U_mOPV2, ES)
  period_mat <- matrix(nrow = nrow(tt_conv_data), ncol = length(period),  0)
  tt_conv_data <- data.frame(tt_conv_data, period_mat)
  names(tt_conv_data) <- c("adm0_name","adm1_name", "region", "period","vaccinetype","source","U_mOPV2", "ES",
                           period)
  
  # Output: Estimate U_d_i for each day d and each campaign i for the probability of detecting an emergence on that day from that campaign
  for (i in 1:nrow(tt_conv_data)){
    start = which(names(tt_conv_data) == round(tt_conv_data[i, "period"],3)) 
    end = start + 59
    tt_conv_data[i, start:end] <- Func_tt_conv(U_mOPV2 = tt_conv_data[i, "U_mOPV2"],
                                               adm0_name = tt_conv_data[i, "adm0_name"],
                                               adm1_name = tt_conv_data[i, "adm1_name"],
                                               region = tt_conv_data[i, "region"], 
                                               ES = tt_conv_data[i, "ES"],
                                               ttl_include = F)
  }
  
  # calculate sum of U_d for each day d across all campaigns i of each vaccine type. 
  daily_U_conv <- data.frame(period = rep(period, each=2))
  daily_U_conv$source <- c("Sabin2", "nOPV2")
  daily_U_conv$U_mOPV2 <- 0
  
  for (i in 9:ncol(tt_conv_data)){
    period_i = names(tt_conv_data)[i]
    U_sabin2_sum <- sum(tt_conv_data[tt_conv_data$source == "Sabin2", i])
    U_nOPV2_sum <- sum(tt_conv_data[tt_conv_data$source == "nOPV2", i])
    
    daily_U_conv[daily_U_conv$period == period_i & daily_U_conv$source == "Sabin2", "U_mOPV2"] <- U_sabin2_sum
    daily_U_conv[daily_U_conv$period == period_i & daily_U_conv$source == "nOPV2", "U_mOPV2"] <- U_nOPV2_sum
  }
  return(daily_U_conv)
}
# Test
daily_U_conv <-  data_province %>%
  mutate(U_mOPV2 = Func_u(p=population_total_sum,
                          q = immunity_weighted,
                          alpha = alpha)) %>%
  # filter(region %in% c("AFRO")) %>%
  filter(africa %in% c("Africa")) %>%
  Func_daily_U_conv()


# Function to calculate cdf before defined dates
Func_cdf_by_date <- function(daily_U_conv, end_date = "2033-01-01", source_input = c("Sabin2", "nOPV2")){
  end_period <- round(year(as.Date(end_date)) +(month(as.Date(end_date))-1)/12, 3)
  daily_U_conv %>%
    filter(source %in% source_input) %>%
    filter(period <= end_period) %>%
    group_by(source) %>%
    summarize(sum(U_mOPV2))
}
# Test
Func_cdf_by_date(daily_U_conv, "2023-11-17")
Func_cdf_by_date(daily_U_conv)
Func_cdf_by_date(daily_U_conv, source_input = c("nOPV2"))

#### Fit alpha value to current sabin2 situation ####
# Count post-switch Sabin-2 emergences
viruses %>% filter(source == "Sabin2", index_isolate == "TRUE", seeding_date > "2016-04-01") %>% 
  # filter(region_who_code %in% c("AFRO")) %>%
  filter(africa %in% c("Africa")) %>%
  # group_by(region_who_code) %>%
  # group_by(year(virus_date)) %>%
  summarize(count = n())
sabin_emerge <- 59 # 56 post-switch Sabin-2 emergences in AFRO. Need to update for Africa (59)

# Begin with default values
theta_size_input = theta_size_mean
theta_u5_input = theta_u5_mean
alpha = 2.305*10^-6 # Updated 8/1/2024 for Africa only
end_date = today()

# Wrapper function
Func_wrapper <- function(
    data_province,
    alpha,
    theta_size_input,
    theta_u5_input,
    end_date = today(),
    sabin_emerge_input = sabin_emerge,
    fast_input = T,
    crude = FALSE,
    lambda = crude_lambda,
    AFRO = F,
    africa = T
){
  if (AFRO == T) {
    data_province <- data_province %>% filter(region %in% c("AFRO"))
  }
  if (africa == T) {
    data_province <- data_province %>% filter(africa %in% c("Africa"))
  }
  
  # calculate total emergences due to sabin2 by date
  Sabin2_Expected <-
    data_province %>%
    filter(source %in% c("Sabin2")) %>%
    mutate(U_mOPV2 = Func_u(p=population_total_sum,
                            q = immunity_weighted,
                            alpha = alpha,
                            theta_size = theta_size_input,
                            theta_u5 = theta_u5_input,
                            crude = crude,
                            lambda = lambda)) %>% 
    Func_daily_U_conv() 
  
  Sabin2_Expected$alpha <- alpha
  Sabin2_Expected$theta_size <- theta_size_input
  Sabin2_Expected$theta_u5 <- theta_u5_input

  Sabin2_Expected_cdf <- Sabin2_Expected %>%
    Func_cdf_by_date(end_date, source = c("Sabin2")) %>% .[,2]
  
  # calculate difference between observed and expected
  delta = sabin_emerge_input / Sabin2_Expected_cdf
  
  # adjust alpha
  alpha_new = as.numeric(alpha*delta)
  
  # Run with new alpha to calculate nOPV2 expectation
  if (fast_input == T){
    daily_U_conv <-
      data_province %>%
      filter(source %in% c("nOPV2")) %>% # to speed it up, only run for nOPV2
      mutate(U_mOPV2 = Func_u(p=population_total_sum, # calculate total emergences due to sabin2
                              q = immunity_weighted,
                              alpha = alpha_new,
                              theta_size = theta_size_input,
                              theta_u5 = theta_u5_input,
                              crude = crude,
                              lambda = lambda)) %>% 
      Func_daily_U_conv() %>% filter(source %in% c("nOPV2"))
  } else {
    daily_U_conv <-
      data_province %>%
      mutate(U_mOPV2 = Func_u(p=population_total_sum, # calculate total emergences due to sabin2
                              q = immunity_weighted,
                              alpha = alpha_new,
                              theta_size = theta_size_input,
                              theta_u5 = theta_u5_input,
                              crude = crude,
                              lambda = lambda)) %>% 
      Func_daily_U_conv()

  }

  daily_U_conv$alpha <- alpha_new
  daily_U_conv$theta_size <- theta_size_input
  daily_U_conv$theta_u5 <- theta_u5_input
  
  # daily_U_conv <- daily_U_conv %>%
  #   bind_rows(Sabin2_Expected) # bind together. note alpha is new for nOPV2 but orig for Sabin2

  return(daily_U_conv)
}

daily_U_conv <- Func_wrapper(data_province,
                             alpha,
                             theta_size_mean,
                             theta_u5_mean,
                             end_date = end_date,
                             fast_input = F)

head(daily_U_conv)

# Estimate expectation
daily_U_conv %>%
  Func_cdf_by_date(end_date, source_input = c("nOPV2")) %>% .[,2] %>% as.numeric()
daily_U_conv %>%
  Func_cdf_by_date(end_date, source_input = c("Sabin2")) %>% .[,2] %>% as.numeric()

#### Crude analysis

daily_U_conv <- Func_wrapper(data_province,
                             alpha,
                             theta_size_mean,
                             theta_u5_mean,
                             end_date = end_date,
                             fast_input = F,
                             crude = T, 
                             africa = F)
daily_U_conv %>%
  Func_cdf_by_date(end_date, source_input = c("Sabin2")) %>% .[,2] %>% as.numeric()

daily_U_conv %>%
  Func_cdf_by_date(end_date, source_input = c("nOPV2")) %>% .[,2] %>% as.numeric()

daily_U_conv %>%
  Func_cdf_by_date(source_input = c("nOPV2")) %>% .[,2] %>% as.numeric()

(crude_expectation <- daily_U_conv %>%
  Func_cdf_by_date(end_date, source_input = c("nOPV2")) %>% .[,2] %>% as.numeric())

1-14/
  qpois(c(0.025, 0.5, 0.975), crude_expectation)

#### Repeat using theta posterior uncertainty ####
load("theta_chain_mcmc_output.rda")
head(thetachain); dim(thetachain)

sample_size = 100
sample_thetas <- thetachain[sample(1:nrow(thetachain), size = sample_size),]

daily_U_conv <- Func_wrapper(data_province,
                             alpha,
                             theta_size_mean,
                             theta_u5_mean,
                             end_date = today(),
                             fast_input = F,
                             crude = F, 
                             africa = T)

daily_U_conv_list <- list(daily_U_conv)

for (i in 1:nrow(sample_thetas)){
  theta_size_input <- sample_thetas[i, 2]
  theta_u5_input <- sample_thetas[i, 1]
  daily_U_conv_temp <- Func_wrapper(data_province,
                               alpha,
                               theta_size_input,
                               theta_u5_input,
                               fast_input = F,
                               crude = F,
                               africa = T,
                               end_date = end_date)
  
  daily_U_conv_temp$iteration = i
  
  daily_U_conv_list <- append(daily_U_conv_list, list(daily_U_conv_temp))
  cat(i,"\n")
}

# Estimate nOPV2 by date
daily_U_conv_list %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  filter(is.na(iteration) == T) %>% #to get the best fitted Theta values
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  # filter(period == 2024.250) %>%
  filter(period == 2027.25) %>%
  summarize(mean = mean(U_mOPV2_cumsum))

daily_U_conv_list %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  filter(period == 2024.50) %>%
  # filter(period == 2027.25) %>%
  ungroup() %>%
  summarize(mean = mean(U_mOPV2_cumsum),
            median = median(U_mOPV2_cumsum),
            mean = mean(U_mOPV2_cumsum),
            upper = quantile(U_mOPV2_cumsum, 0.975),
            lower = quantile(U_mOPV2_cumsum, 0.025))

# Create theta uncertainty bounds for novel expectation
bounds_nOPV2 <-
  daily_U_conv_list %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize(lower = quantile(U_mOPV2_cumsum, 0.025),
            upper = quantile(U_mOPV2_cumsum, 0.975))

bounds_nOPV2[bounds_nOPV2$period == 2024.250,]
tail(bounds_nOPV2)

# Create theta uncertainty bounds for Sabin2 expectation
bounds_Sabin2 <-
  daily_U_conv_list %>% 
  bind_rows() %>%
  filter(source == c("Sabin2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize(lower = quantile(U_mOPV2_cumsum, 0.025),
            upper = quantile(U_mOPV2_cumsum, 0.975))

bounds_Sabin2[bounds_Sabin2$period == 2024.250,]
tail(bounds_Sabin2)

#### Plot daily_U_conv data ####
# Plot points and lines for U_d by period
ggplot() +
  geom_vline(xintercept = 2024.5, size = 1) +
  geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) +
  geom_line(data = daily_U_conv, aes(x = period, y = U_mOPV2, color = source), size = 1) +
  # geom_point(data = daily_U_conv, aes(x = period, y = U_mOPV2, color = source), size = 1) +
  geom_point(data = viruses_count_period %>% filter(emergences > 0), aes(x = period, y = emergences, color = source, shape = source), size = 2) +
  theme_bw() +
  ylab("Expected cVDPV2 Emergences\n(Assuming Seeding at mOPV2 Rate)") +
  xlab("Month") +
  # scale_y_continuous(limits = c(0, 8)) +
  scale_shape_manual(values = c(19,1)) +
  theme(legend.position = c(0.85, 0.8)) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")

# Plot cumulative counts and expectation
daily_U_conv <- daily_U_conv %>%
  ungroup() %>%  group_by(source) %>%
  arrange(period) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) 
viruses_count_period <- viruses_count_period %>%
  ungroup() %>% group_by(source) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))
viruses_count_period_region <- viruses_count_period_region %>%
  ungroup() %>% group_by(source, region_who_code) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))
viruses_count_period_africa <- viruses_count_period_africa %>%
  ungroup() %>% group_by(source, africa) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))

# Plot cumsum points and lines for U_d by period
temp <- left_join(daily_U_conv, viruses_count_period_africa %>%
                    filter(africa %in% c("Africa")), by = c("period", "source"))
temp[is.na(temp$emergences), "emergences"] <- 0
temp <- temp %>% ungroup() %>% group_by(source) %>% arrange(period) %>% mutate(emergences_cumsum = cumsum(emergences))
temp[temp$period > 2024.5, c("emergences", "emergences_cumsum")] <- NA
temp <- temp %>% 
  select(c("period", "source", "U_mOPV2_cumsum", "emergences_cumsum")) %>%
  pivot_longer(cols = !c("period", "source"))
temp <- temp %>% filter(!(source %in%  c("Sabin2") & name %in% c("U_mOPV2_cumsum")))

# for paper
fig_cum_emergences <- 
  ggplot() +
    # geom_ribbon(data = bounds, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
      geom_ribbon(data = bounds_nOPV2, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
      # geom_ribbon(data = bounds_Sabin2, aes(x = period, ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5, size = 1, linetype = "dashed") +
    geom_line(data = temp, aes(x = period, y = value, color = source, linetype = name), size = 1) +
    geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
    geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) + theme_bw() +
    ylab("Cumulative cVDPV2\nEmergences in Africa") +
    scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
    scale_linetype_discrete(name = "Emergences",
                            labels = c("Observed", "Expected based\non Sabin 2 Rate")) +
    # theme(legend.position = c(0.85, 0.3)) +
    force_panelsizes(rows = unit(4, "in"),
                     cols = unit(5, "in")) +
    scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")
ggsave(plot = fig_cum_emergences,"figures/Cumulative Emergences Africa.png", device = "png", units = "in", width = 7, height = 5)

anim <- ggplot() +
  # geom_ribbon(data = bounds, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
  geom_ribbon(data = bounds_nOPV2, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
  # geom_ribbon(data = bounds_Sabin2, aes(x = period, ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5, size = 1, linetype = "dashed") +
  geom_line(data = temp, aes(x = period, y = value, color = source, linetype = name), size = 1) +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) + theme_bw() +
  ylab("Cumulative cVDPV2 Emergence\nin Africa") +
  scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_linetype_discrete(name = "Emergences",
                          labels = c("Observed", "Expected based\non Sabin 2 Rate")) +
  # theme(legend.position = c(0.85, 0.3)) +
  # force_panelsizes(rows = unit(4, "in"),
  #                  cols = unit(5, "in")) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") + 
  transition_reveal(period)
anim_save("figures/Cumulative Emergences Africa.gif", anim)


# Combined figure
plot_layout(fig_cum_doses / fig_cum_emergences)
ggsave("figures/Cumulative Africa.png", device = "png", units = "in", width = 7, height = 10)

#### Global crude analysis ####


#### Sensitivity analysis: DRC Only ####
polis_pops %>% filter(adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>% group_by(source)%>% filter(start_date < today()) %>% summarize(sum = sum(target_pop))

viruses_count_period_DRC <- viruses %>% 
  filter(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  filter(seeding_date > "2016-03-01") %>%
  group_by(period, source) %>% 
  summarize(emergences = sum(index_isolate==TRUE))

viruses %>% 
  filter(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  filter(source == "Sabin2", index_isolate == "TRUE", seeding_date > "2016-03-01") %>% 
  View() # 17 Sabin2-emergences
sabin_emerge_DRC <- 17
alpha_DRC <- 6.6e-06 # based on 17 DRC emergences from Sabin2 use

daily_U_conv_DRC <- Func_wrapper(data_province %>% filter(adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO"),
                                 alpha = alpha_DRC,
                                 theta_size_mean,
                                 theta_u5_mean,
                                 sabin_emerge_input = sabin_emerge_DRC,
                                 end_date = end_date,
                                 fast_input = F)
# Check new alpha versus old
daily_U_conv_DRC$alpha[1] / alpha_DRC

# Estimate nOPV2 by date
daily_U_conv_DRC %>%
  Func_cdf_by_date(end_date)

# Estimate nOPV2 total
daily_U_conv_DRC %>%
  Func_cdf_by_date()

# Plot points and lines for U_d
ggplot() +
  geom_vline(xintercept = 2024.5, size = 1) +
  # geom_vline(xintercept = 2021.667, color = "red", alpha = 0.25, size = 1) +
  geom_line(data = daily_U_conv_DRC, aes(x = period, y = U_mOPV2, color = source), size = 1) +
  geom_point(data = viruses_count_period_DRC %>% filter(emergences > 0), aes(x = period, y = emergences, color = source, shape = source), size = 2) +
  theme_bw() +
  ylab("Expected cVDPV2 Emergences\n(Assuming Seeding at mOPV2 Rate)") +
  xlab("Month") +
  scale_y_continuous(limits = c(0, 2.2)) +
  scale_shape_manual(values = c(19,1)) +
  theme(legend.position = c(0.85, 0.85)) +
  ggtitle("DRC Only") +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")

# Plot cumulative counts and expectation
daily_U_conv_DRC <- daily_U_conv_DRC %>%
  ungroup() %>%  group_by(source) %>%
  arrange(period) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) 
viruses_count_period_DRC <- viruses_count_period_DRC %>%
  ungroup() %>% group_by(source) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))

# Plot cumsum points and lines for U_d by period
temp_DRC <- left_join(daily_U_conv_DRC, viruses_count_period_DRC, by = c("period", "source"))
temp_DRC[is.na(temp_DRC$emergences), "emergences"] <- 0
temp_DRC <- temp_DRC %>% ungroup() %>% group_by(source) %>% arrange(period) %>% mutate(emergences_cumsum = cumsum(emergences))
temp_DRC[temp_DRC$period > 2024.5, c("emergences", "emergences_cumsum")] <- NA
temp_DRC <- temp_DRC %>% 
  select(c("period", "source", "U_mOPV2_cumsum", "emergences_cumsum")) %>%
  pivot_longer(cols = !c("period", "source"))
temp_DRC <- temp_DRC %>% filter(!(source %in%  c("Sabin2") & name %in% c("U_mOPV2_cumsum")))

fig_cum_emergences_DRC <- 
  ggplot() +
  # geom_ribbon(data = bounds, aes(x = period, ymin = lower, ymax = upper), fill = "pink", size = 1, linetype = "dashed") +
  geom_line(data = temp_DRC, aes(x = period, y = value, color = source, linetype = name), size = 1) +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) +theme_bw() +
  ylab("Cumulative cVDPV2\nEmergences in DRC") +
  scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_linetype_discrete(name = "Emergences",
                          labels = c("Observed", "Expected based\non Sabin 2 Rate")) +
  # theme(legend.position = c(0.85, 0.3)) +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(5, "in")) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")
ggsave(plot = fig_cum_emergences_DRC,"figures/Cumulative Emergences DRC.png", device = "png", units = "in", width = 7, height = 5)

fig_cum_doses_DRC <- 
  polis_pops %>%
  filter(!is.na(target_pop)) %>%
  filter(adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  ungroup() %>%
  group_by(source, period) %>%
  reframe(period = unique(period),
          doses = sum(target_pop)) %>%
  group_by(source) %>%
  arrange(period) %>%
  reframe(period = unique(period),
          doses_cumsum = cumsum(doses)) %>%
  ggplot() +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, alpha = 0.25, color = "red", size = 1) +
  geom_line(aes(x = period, y = doses_cumsum/1e6, color = source), size = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(2016.0, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_y_continuous(name = "Cumulative Doses in\nDRC (Million)") +
  scale_color_discrete(name = "Vaccine Type") +
  force_panelsizes(rows = unit(2, "in"),
                   cols = unit(5, "in"))
ggsave(plot = fig_cum_doses_DRC, "figures/Cumulative Doses DRC.png", device = "png", units = "in", width = 7, height = 3)

# Combined figure
plot_layout(fig_cum_doses_DRC / fig_cum_emergences_DRC)
ggsave("figures/Cumulative DRC.png", device = "png", units = "in", width = 7, height = 10)

# Repeat with uncertainty analysis
daily_U_conv_list_DRC <- list(daily_U_conv_DRC)

for (i in 1:nrow(sample_thetas)){
  theta_size_input <- sample_thetas[i, 2]
  theta_u5_input <- sample_thetas[i, 1]
  daily_U_conv_temp_DRC <- Func_wrapper(data_province %>% filter(adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO"),
                                    alpha_DRC,
                                    theta_size_input,
                                    theta_u5_input,
                                    fast_input = F,
                                    sabin_emerge_input = sabin_emerge_DRC,
                                    crude = F,
                                    end_date = end_date)
  
  daily_U_conv_temp_DRC$iteration = i
  
  daily_U_conv_list_DRC <- append(daily_U_conv_list_DRC, list(daily_U_conv_temp_DRC))
  cat(i,"\n")
}

# Estimate nOPV2 by date
daily_U_conv_list_DRC %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  filter(is.na(iteration) == T) %>% #to get the best fitted Theta values
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  # filter(period == 2024.50) %>%
  filter(period == 2027.25) %>%
  summarize(mean = mean(U_mOPV2_cumsum))

daily_U_conv_list_DRC %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  filter(period == 2024.50) %>%
  # filter(period == 2027.25) %>%
  ungroup() %>%
  summarize(mean = mean(U_mOPV2_cumsum),
            median = median(U_mOPV2_cumsum),
            mean = mean(U_mOPV2_cumsum),
            upper = quantile(U_mOPV2_cumsum, 0.975),
            lower = quantile(U_mOPV2_cumsum, 0.025))

# Create theta uncertainty bounds for novel expectation
bounds_nOPV2_DRC <-
  daily_U_conv_list_DRC %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize(lower = quantile(U_mOPV2_cumsum, 0.025),
            upper = quantile(U_mOPV2_cumsum, 0.975))

bounds_nOPV2_DRC[bounds_nOPV2_DRC$period == 2024.50,]
tail(bounds_nOPV2_DRC)

# Create theta uncertainty bounds for Sabin2 expectation
bounds_Sabin2_DRC <-
  daily_U_conv_list_DRC %>% 
  bind_rows() %>%
  filter(source == c("Sabin2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize(lower = quantile(U_mOPV2_cumsum, 0.025),
            upper = quantile(U_mOPV2_cumsum, 0.975))

bounds_Sabin2_DRC[bounds_Sabin2_DRC$period == 2024.50,]
tail(bounds_Sabin2_DRC)

# Plot cumulative counts and expectation
daily_U_conv_DRC <- daily_U_conv_DRC %>%
  ungroup() %>%  group_by(source) %>%
  arrange(period) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) 
viruses_count_period_DRC <- viruses_count_period_DRC %>%
  ungroup() %>% group_by(source) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))

# Plot cumsum points and lines for U_d by period
temp_DRC <- left_join(daily_U_conv_DRC, viruses_count_period_DRC, by = c("period", "source"))
temp_DRC[is.na(temp_DRC$emergences), "emergences"] <- 0
temp_DRC <- temp_DRC %>% ungroup() %>% group_by(source) %>% arrange(period) %>% mutate(emergences_cumsum = cumsum(emergences))
temp_DRC[temp_DRC$period > 2024.5, c("emergences", "emergences_cumsum")] <- NA
temp_DRC <- temp_DRC %>% 
  select(c("period", "source", "U_mOPV2_cumsum", "emergences_cumsum")) %>%
  pivot_longer(cols = !c("period", "source"))
temp_DRC <- temp_DRC %>% filter(!(source %in%  c("Sabin2") & name %in% c("U_mOPV2_cumsum")))

# for paper
fig_cum_emergences_DRC <- 
  ggplot() +
  # geom_ribbon(data = bounds_DRC, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
  geom_ribbon(data = bounds_nOPV2_DRC, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
  # geom_ribbon(data = bounds_Sabin2_DRC, aes(x = period, ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5, size = 1, linetype = "dashed") +
  geom_line(data = temp_DRC, aes(x = period, y = value, color = source, linetype = name), size = 1) +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) + theme_bw() +
  ylab("Cumulative cVDPV2\nEmergences in DRC") +
  scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_linetype_discrete(name = "Emergences",
                          labels = c("Observed", "Expected based\non Sabin 2 Rate")) +
  # theme(legend.position = c(0.85, 0.3)) +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(5, "in")) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")
ggsave(plot = fig_cum_emergences_DRC,"figures/Cumulative Emergences DRC.png", device = "png", units = "in", width = 7, height = 5)

# Combined figure
plot_layout(fig_cum_doses_DRC / fig_cum_emergences_DRC)
ggsave("figures/Cumulative DRC.png", device = "png", units = "in", width = 7, height = 10)

# Boxplot of SIA size by vaccine type (For paper)
sia_target_box_DRC <- polis_pops %>%
  filter(adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  group_by(source, parentactivitycode) %>% 
  filter(start_date >= "2016-05-01") %>%
  summarize(sum = sum(target_pop, na.rm=T)) %>%
  ungroup() %>%
  ggplot() +
  geom_boxplot(aes(x = source, y = sum/1e6, color = source)) +
  scale_y_log10(name = "SIA Target Population (Million)", 
                breaks = 10^(2-seq(1,5)),
                minor_breaks = .5*10^(2-seq(1,5))) +
  xlab("") +
  theme_bw() +
  guides(color="none")
ggsave(plot = sia_target_box_DRC, "figures/target pop_DRC.png", device = "png", units = "in", width = 3, height = 3)

# Summary of pre-campaign immunity by vaccine type at the SIA level (For paper)
data_sia_DRC <- data %>% 
  filter(adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  group_by(parentactivitycode, period) %>%
  reframe(target_pop_sum = sum(target_pop), 
          population_total_sum = sum(population_total),
          immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
          vaccinetype = unique(vaccinetype),
          Region = unique(region),
          adm0_name = unique(adm0_name),
          source = unique(source),
          start_date = min(start_date),
          quarter = unique(quarter),
          week = unique(week))
sum(data_sia_DRC[data_sia_DRC$vaccinetype == "nOPV2", "target_pop_sum"], na.rm=T)

temp_DRC <- data_sia_DRC %>%
  filter(is.na(immunity_weighted)==F) %>%
  filter(start_date >= "2016-05-01")

temp_DRC %>% 
  group_by(source) %>%
  summarize(min = min(immunity_weighted),
            q1 = quantile(immunity_weighted, 0.25),
            median = median(immunity_weighted),
            mean = weighted.mean(immunity_weighted, target_pop_sum, na.rm=T),
            q3 = quantile(immunity_weighted, 0.75),
            max = max(immunity_weighted))

wilcox.test(temp_DRC$immunity_weighted ~ temp_DRC$source,
            exact = FALSE)

sia_immunity_box_DRC <- temp_DRC %>%
  ggplot(aes(x = source, y = immunity_weighted, color = source)) +
  # geom_violin() +
  geom_boxplot() +
  # geom_jitter(alpha = 0.5) +
  theme_bw() +
  guides(color = "none") +
  ylab("Pre-Campaign Type-2 Immunity") +
  xlab("")
ggsave(plot = sia_immunity_box_DRC, "figures/immunity_DRC.png", device = "png", units = "in", width = 3, height = 3)

# SIA target and immunity boxplots (For paper)
plot_layout(sia_target_box_DRC / sia_immunity_box_DRC)
ggsave("figures/SIA target and immunity_DRC.png", device = "png", units = "in", width = 3, height = 6)

#### Sensitivity analysis: Nigeria Only ####
polis_pops %>% filter(adm0_name == "NIGERIA") %>% group_by(source)%>% filter(start_date < today()) %>% summarize(sum = sum(target_pop, na.rm=T))

viruses_count_period_NIE <- viruses %>% 
  filter(admin0name == "NIGERIA") %>%
  filter(seeding_date > "2016-03-01") %>%
  group_by(period, source) %>% 
  summarize(emergences = sum(index_isolate==TRUE))

viruses %>% 
  filter(admin0name == "NIGERIA") %>%
  filter(source == "Sabin2", index_isolate == "TRUE", seeding_date > "2016-03-01") %>%
  View() # 10 Sabin2-emergences
sabin_emerge_NIE <- 10
alpha_Nigeria <- 1.58e-06 # based on 10 NIGERIA emergences from Sabin2 use

daily_U_conv_NIE <- Func_wrapper(data_province %>% filter(adm0_name == "NIGERIA"),
                             alpha = alpha_Nigeria,
                             theta_size_mean,
                             theta_u5_mean,
                             sabin_emerge_input = sabin_emerge_NIE,
                             end_date = end_date,
                             fast_input = F)
# Check new alpha versus old
daily_U_conv_NIE$alpha[1] / alpha_Nigeria

# Estimate nOPV2 by date
daily_U_conv_NIE %>%
  Func_cdf_by_date(end_date)

# Estimate nOPV2 total
daily_U_conv_NIE %>%
  Func_cdf_by_date()

# Plot cumulative counts and expectation
daily_U_conv_NIE <- daily_U_conv_NIE %>%
  ungroup() %>%  group_by(source) %>%
  arrange(period) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) 
viruses_count_period_NIE <- viruses_count_period_NIE %>%
  ungroup() %>% group_by(source) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))

# Plot cumsum points and lines for U_d by period
temp_NIE <- left_join(daily_U_conv_NIE, viruses_count_period_NIE, by = c("period", "source"))
temp_NIE[is.na(temp_NIE$emergences), "emergences"] <- 0
temp_NIE <- temp_NIE %>% ungroup() %>% group_by(source) %>% arrange(period) %>% mutate(emergences_cumsum = cumsum(emergences))
temp_NIE[temp_NIE$period > 2024.5, c("emergences", "emergences_cumsum")] <- NA
temp_NIE <- temp_NIE %>% 
  select(c("period", "source", "U_mOPV2_cumsum", "emergences_cumsum")) %>%
  pivot_longer(cols = !c("period", "source"))
temp_NIE <- temp_NIE %>% filter(!(source %in%  c("Sabin2") & name %in% c("U_mOPV2_cumsum")))

fig_cum_emergences_NIE <- 
  ggplot() +
  # geom_ribbon(data = bounds, aes(x = period, ymin = lower, ymax = upper), fill = "pink", size = 1, linetype = "dashed") +
  geom_line(data = temp_NIE, aes(x = period, y = value, color = source, linetype = name), size = 1) +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) +theme_bw() +
  ylab("Cumulative cVDPV2 Emergences in\nNigeria") +
  scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_linetype_discrete(name = "Emergences",
                          labels = c("Observed", "Expected based\non Sabin 2 Rate")) +
  # theme(legend.position = c(0.85, 0.3)) +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(5, "in")) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")
ggsave(plot = fig_cum_emergences_NIE,"figures/Cumulative Emergences NIE.png", device = "png", units = "in", width = 7, height = 5)

fig_cum_doses_NIE <- 
  polis_pops %>%
  filter(!is.na(target_pop)) %>%
  filter(adm0_name == "NIGERIA") %>%
  ungroup() %>%
  group_by(source, period) %>%
  reframe(period = unique(period),
          doses = sum(target_pop)) %>%
  group_by(source) %>%
  arrange(period) %>%
  reframe(period = unique(period),
          doses_cumsum = cumsum(doses)) %>%
  ggplot() +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, alpha = 0.25, color = "red", size = 1) +
  geom_line(aes(x = period, y = doses_cumsum/1e6, color = source), size = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(2016.0, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_y_continuous(name = "Cumulative Doses in\nNigeria (Million)") +
  scale_color_discrete(name = "Vaccine Type") +
  force_panelsizes(rows = unit(2, "in"),
                   cols = unit(5, "in"))
ggsave(plot = fig_cum_doses_NIE, "figures/Cumulative Doses NIE.png", device = "png", units = "in", width = 7, height = 3)

# Combined figure
plot_layout(fig_cum_doses_NIE / fig_cum_emergences_NIE)
ggsave("figures/Cumulative NIE.png", device = "png", units = "in", width = 7, height = 10)

# Repeat with uncertainty analysis
daily_U_conv_list_NIE <- list(daily_U_conv_NIE)

for (i in 1:nrow(sample_thetas)){
  theta_size_input <- sample_thetas[i, 2]
  theta_u5_input <- sample_thetas[i, 1]
  daily_U_conv_temp_NIE <- Func_wrapper(data_province %>% filter(adm0_name == "NIGERIA"),
                                        alpha = alpha_Nigeria,
                                        theta_size_input,
                                        theta_u5_input,
                                        fast_input = F,
                                        sabin_emerge_input = sabin_emerge_NIE,
                                        crude = F,
                                        end_date = end_date)
  
  daily_U_conv_temp_NIE$iteration = i
  
  daily_U_conv_list_NIE <- append(daily_U_conv_list_NIE, list(daily_U_conv_temp_NIE))
  cat(i,"\n")
}

# Estimate nOPV2 by date
daily_U_conv_list_NIE %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  filter(is.na(iteration) == T) %>% #to get the best fitted Theta values
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  filter(period == 2024.250) %>%
  # filter(period == 2027.25) %>%
  summarize(mean = mean(U_mOPV2_cumsum))

daily_U_conv_list_NIE %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  filter(period == 2024.250) %>%
  # filter(period == 2027.25) %>%
  ungroup() %>%
  summarize(mean = mean(U_mOPV2_cumsum),
            median = median(U_mOPV2_cumsum),
            mean = mean(U_mOPV2_cumsum),
            upper = quantile(U_mOPV2_cumsum, 0.975),
            lower = quantile(U_mOPV2_cumsum, 0.025))

# Create theta uncertainty bounds for novel expectation
bounds_nOPV2_NIE <-
  daily_U_conv_list_NIE %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize(lower = quantile(U_mOPV2_cumsum, 0.025),
            upper = quantile(U_mOPV2_cumsum, 0.975))

bounds_nOPV2_NIE[bounds_nOPV2_NIE$period == 2024.250,]
tail(bounds_nOPV2_NIE)

# Create theta uncertainty bounds for Sabin2 expectation
bounds_Sabin2_NIE <-
  daily_U_conv_list_NIE %>% 
  bind_rows() %>%
  filter(source == c("Sabin2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize(lower = quantile(U_mOPV2_cumsum, 0.025),
            upper = quantile(U_mOPV2_cumsum, 0.975))

bounds_Sabin2_NIE[bounds_Sabin2_NIE$period == 2024.250,]
tail(bounds_Sabin2_NIE)

# Plot cumulative counts and expectation
daily_U_conv_NIE <- daily_U_conv_NIE %>%
  ungroup() %>%  group_by(source) %>%
  arrange(period) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) 
viruses_count_period_NIE <- viruses_count_period_NIE %>%
  ungroup() %>% group_by(source) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))

# Plot cumsum points and lines for U_d by period
temp_NIE <- left_join(daily_U_conv_NIE, viruses_count_period_NIE, by = c("period", "source"))
temp_NIE[is.na(temp_NIE$emergences), "emergences"] <- 0
temp_NIE <- temp_NIE %>% ungroup() %>% group_by(source) %>% arrange(period) %>% mutate(emergences_cumsum = cumsum(emergences))
temp_NIE[temp_NIE$period > 2024.333, c("emergences", "emergences_cumsum")] <- NA
temp_NIE <- temp_NIE %>% 
  select(c("period", "source", "U_mOPV2_cumsum", "emergences_cumsum")) %>%
  pivot_longer(cols = !c("period", "source"))
temp_NIE <- temp_NIE %>% filter(!(source %in%  c("Sabin2") & name %in% c("U_mOPV2_cumsum")))

# for paper
fig_cum_emergences_NIE <- 
  ggplot() +
  # geom_ribbon(data = bounds_NIE, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
  geom_ribbon(data = bounds_nOPV2_NIE, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
  # geom_ribbon(data = bounds_Sabin2_NIE, aes(x = period, ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5, size = 1, linetype = "dashed") +
  geom_line(data = temp_NIE, aes(x = period, y = value, color = source, linetype = name), size = 1) +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) + theme_bw() +
  ylab("Cumulative cVDPV2\nEmergences in Nigeria") +
  scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_linetype_discrete(name = "Emergences",
                          labels = c("Observed", "Expected based\non Sabin 2 Rate")) +
  # theme(legend.position = c(0.85, 0.3)) +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(5, "in")) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")
ggsave(plot = fig_cum_emergences_NIE,"figures/Cumulative Emergences NIE.png", device = "png", units = "in", width = 7, height = 5)

# Combined figure
plot_layout(fig_cum_doses_NIE / fig_cum_emergences_NIE)
ggsave("figures/Cumulative NIE.png", device = "png", units = "in", width = 7, height = 10)

# Boxplot of SIA size by vaccine type (For paper)
sia_target_box_NIE <- polis_pops %>%
  filter(adm0_name == "NIGERIA") %>%
  group_by(source, parentactivitycode) %>% 
  filter(start_date >= "2016-05-01") %>%
  summarize(sum = sum(target_pop, na.rm=T)) %>%
  ungroup() %>%
  ggplot() +
  geom_boxplot(aes(x = source, y = sum/1e6, color = source)) +
  scale_y_log10(name = "SIA Target Population (Million)", 
                breaks = 10^(2-seq(1,5)),
                minor_breaks = .5*10^(2-seq(1,5))) +
  xlab("") +
  theme_bw() +
  guides(color="none")
ggsave(plot = sia_target_box_NIE, "figures/target pop_NIE.png", device = "png", units = "in", width = 3, height = 3)

# Summary of pre-campaign immunity by vaccine type at the SIA level (For paper)
data_sia_NIE <- data %>% 
  filter(adm0_name == "NIGERIA") %>%
  group_by(parentactivitycode, period) %>%
  reframe(target_pop_sum = sum(target_pop), 
          population_total_sum = sum(population_total),
          immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
          vaccinetype = unique(vaccinetype),
          Region = unique(Region),
          adm0_name = unique(adm0_name),
          source = unique(source),
          start_date = min(start_date),
          quarter = unique(quarter),
          week = unique(week))
sum(data_sia_NIE[data_sia_NIE$vaccinetype == "nOPV2", "target_pop_sum"], na.rm=T)

temp_NIE <- data_sia_NIE %>%
  filter(is.na(immunity_weighted)==F) %>%
  filter(start_date >= "2016-05-01")

temp_NIE %>% 
  group_by(source) %>%
  summarize(min = min(immunity_weighted),
            q1 = quantile(immunity_weighted, 0.25),
            median = median(immunity_weighted),
            mean = weighted.mean(immunity_weighted, target_pop_sum, na.rm=T),
            q3 = quantile(immunity_weighted, 0.75),
            max = max(immunity_weighted))

wilcox.test(temp_NIE$immunity_weighted ~ temp_NIE$source,
            exact = FALSE)

sia_immunity_box_NIE <- temp_NIE %>%
  ggplot(aes(x = source, y = immunity_weighted, color = source)) +
  # geom_violin() +
  geom_boxplot() +
  # geom_jitter(alpha = 0.5) +
  theme_bw() +
  guides(color = "none") +
  ylab("Pre-Campaign Type-2 Immunity") +
  xlab("")
ggsave(plot = sia_immunity_box_NIE, "figures/immunity_NIE.png", device = "png", units = "in", width = 3, height = 3)

# SIA target and immunity boxplots (For paper)
plot_layout(sia_target_box_NIE / sia_immunity_box_NIE)
ggsave("figures/SIA target and immunity_NIE.png", device = "png", units = "in", width = 3, height = 6)

# Compare timing of nOPV2 and Sabin2 SIAs. Hypothesis that variable NPEV prevalance may influence emergence risk
polis_pops %>% 
  filter(adm0_name == "NIGERIA") %>%
  filter(start_date > "2016-05-01") %>%
  mutate(season = quarter - floor(quarter)) %>% 
  group_by(season, source) %>% 
  summarize(count = n(), sum = sum(target_pop)) %>%
  ggplot(aes(x = season, y = sum/1e6, fill = source)) +
    geom_col(position = "dodge", stat = "identity") +
    theme_bw() +
    scale_y_continuous(name = "Doses (million)") +
    scale_x_continuous(name = "Quarter")

#### Triple boxplot ####
plot_layout(ncol = 3, nrow = 2, sia_target_box + sia_immunity_box + 
              sia_target_box_NIE + sia_immunity_box_NIE + 
              sia_target_box_DRC + sia_immunity_box_DRC,
            byrow = T)
ggsave("figures/SIA target and immunity_AFR NIE DRC.png", device = "png", units = "in", width = 6, height = 6)

#### Sensitivity analysis: clusters of emergences ####
viruses %>%
  filter(index_isolate == TRUE, admin0name %in% c("ANGOLA", "CENTRAL AFRICAN REPUBLIC", "DEMOCRATIC REPUBLIC OF THE CONGO")) %>%
  select(vdpv_emergence_group_name, admin0name, admin1name, virus_date, vdpv_nt_changes_from_sabin) %>%
  arrange(admin0name, virus_date) %>%
  View()

cluster_emergences <- c("ANG-HUI-1", "ANG-LNO-2", "ANG-LUA-1",
                        #"ANG_LNO-2", #keep one from each cluster
                        #"CAF-BAM-1", #keep one from each cluster
                        "CAF-BIM-3", "CAF-BAM-2", "CAF-BIM-2", "CAF-BIM-1")

viruses_count_period_cluster <- viruses %>% 
  filter(!(vdpv_emergence_group_name %in% cluster_emergences)) %>%
  filter(region_who_code %in% c("AFRO")) %>%
  filter(seeding_date > "2016-03-01") %>%
  group_by(period, source) %>% 
  summarize(emergences = sum(index_isolate==TRUE))
viruses_count_period_cluster %>% ungroup() %>%
  filter(source == "Sabin2") %>% summarize(sum = sum(emergences)) # 48 non cluster AFRO emergences
sabin_emerge_cluster <- 49

alpha_cluster <- 2.41e-06

daily_U_conv_cluster <- Func_wrapper(data_province,
                                 alpha = alpha_cluster,
                                 theta_size_mean,
                                 theta_u5_mean,
                                 sabin_emerge_input = sabin_emerge_cluster,
                                 end_date = end_date,
                                 AFRO = T,
                                 fast_input = F)
# Check new alpha versus old
daily_U_conv_cluster$alpha[1] / alpha_cluster

# Estimate nOPV2 by date
daily_U_conv_cluster %>%
  Func_cdf_by_date(end_date)

# Estimate nOPV2 total
daily_U_conv_cluster %>%
  Func_cdf_by_date()

# Plot cumulative counts and expectation
daily_U_conv_cluster <- daily_U_conv_cluster %>%
  ungroup() %>%  group_by(source) %>%
  arrange(period) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) 
viruses_count_period_cluster <- viruses_count_period_cluster %>%
  ungroup() %>% group_by(source) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))

# Plot cumsum points and lines for U_d by period
temp_cluster <- left_join(daily_U_conv_cluster, viruses_count_period_cluster, by = c("period", "source"))
temp_cluster[is.na(temp_cluster$emergences), "emergences"] <- 0
temp_cluster <- temp_cluster %>% ungroup() %>% group_by(source) %>% arrange(period) %>% mutate(emergences_cumsum = cumsum(emergences))
temp_cluster[temp_cluster$period > 2024.25, c("emergences", "emergences_cumsum")] <- NA
temp_cluster <- temp_cluster %>% 
  select(c("period", "source", "U_mOPV2_cumsum", "emergences_cumsum")) %>%
  pivot_longer(cols = !c("period", "source"))
temp_cluster <- temp_cluster %>% filter(!(source %in%  c("Sabin2") & name %in% c("U_mOPV2_cumsum")))

fig_cum_emergences_cluster <- 
  ggplot() +
  # geom_ribbon(data = bounds, aes(x = period, ymin = lower, ymax = upper), fill = "pink", size = 1, linetype = "dashed") +
  geom_line(data = temp_cluster, aes(x = period, y = value, color = source, linetype = name), size = 1) +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) +theme_bw() +
  ylab("Cumulative cVDPV2 Emergences in\nCollapsing Clusters") +
  scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_linetype_discrete(name = "Emergences",
                          labels = c("Observed", "Expected based\non Sabin 2 Rate")) +
  # theme(legend.position = c(0.85, 0.3)) +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(5, "in")) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")
ggsave(plot = fig_cum_emergences_cluster,"figures/Cumulative Emergences cluster.png", device = "png", units = "in", width = 7, height = 5)

# Repeat with uncertainty analysis
daily_U_conv_list_cluster <- list(daily_U_conv_cluster)

for (i in 1:nrow(sample_thetas)){
  theta_size_input <- sample_thetas[i, 2]
  theta_u5_input <- sample_thetas[i, 1]
  daily_U_conv_temp_cluster <- Func_wrapper(data_province,
                                        alpha_cluster,
                                        theta_size_input,
                                        theta_u5_input,
                                        fast_input = F,
                                        sabin_emerge_input = sabin_emerge_cluster,
                                        crude = F,
                                        AFRO = T,
                                        end_date = end_date)
  
  daily_U_conv_temp_cluster$iteration = i
  
  daily_U_conv_list_cluster <- append(daily_U_conv_list_cluster, list(daily_U_conv_temp_cluster))
  cat(i,"\n")
}

# Estimate nOPV2 by date
daily_U_conv_list_cluster %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  filter(is.na(iteration) == T) %>% #to get the best fitted Theta values
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  filter(period == 2024.250) %>%
  # filter(period == 2027.25) %>%
  summarize(mean = mean(U_mOPV2_cumsum))

daily_U_conv_list_cluster %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  filter(period == 2024.250) %>%
  # filter(period == 2027.25) %>%
  ungroup() %>%
  summarize(mean = mean(U_mOPV2_cumsum),
            median = median(U_mOPV2_cumsum),
            mean = mean(U_mOPV2_cumsum),
            upper = quantile(U_mOPV2_cumsum, 0.975),
            lower = quantile(U_mOPV2_cumsum, 0.025))

# Create theta uncertainty bounds for novel expectation
bounds_nOPV2_cluster <-
  daily_U_conv_list_cluster %>% 
  bind_rows() %>%
  filter(source == c("nOPV2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize(lower = quantile(U_mOPV2_cumsum, 0.025),
            upper = quantile(U_mOPV2_cumsum, 0.975))

bounds_nOPV2_cluster[bounds_nOPV2_cluster$period == 2024.250,]
tail(bounds_nOPV2_cluster)

# Create theta uncertainty bounds for Sabin2 expectation
bounds_Sabin2_cluster <-
  daily_U_conv_list_cluster %>% 
  bind_rows() %>%
  filter(source == c("Sabin2")) %>%
  group_by(iteration) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) %>%
  ungroup() %>%
  group_by(period) %>%
  summarize(lower = quantile(U_mOPV2_cumsum, 0.025),
            upper = quantile(U_mOPV2_cumsum, 0.975))

bounds_Sabin2_cluster[bounds_Sabin2_cluster$period == 2024.250,]
tail(bounds_Sabin2_cluster)

# Plot cumulative counts and expectation
daily_U_conv_cluster <- daily_U_conv_cluster %>%
  ungroup() %>%  group_by(source) %>%
  arrange(period) %>%
  mutate(U_mOPV2_cumsum = cumsum(U_mOPV2)) 
viruses_count_period_cluster <- viruses_count_period_cluster %>%
  ungroup() %>% group_by(source) %>%
  arrange(period) %>%
  mutate(emergences_cumsum = cumsum(emergences))

# Plot cumsum points and lines for U_d by period
temp_cluster <- left_join(daily_U_conv_cluster, viruses_count_period_cluster, by = c("period", "source"))
temp_cluster[is.na(temp_cluster$emergences), "emergences"] <- 0
temp_cluster <- temp_cluster %>% ungroup() %>% group_by(source) %>% arrange(period) %>% mutate(emergences_cumsum = cumsum(emergences))
temp_cluster[temp_cluster$period > 2024.333, c("emergences", "emergences_cumsum")] <- NA
temp_cluster <- temp_cluster %>% 
  select(c("period", "source", "U_mOPV2_cumsum", "emergences_cumsum")) %>%
  pivot_longer(cols = !c("period", "source"))
temp_cluster <- temp_cluster %>% filter(!(source %in%  c("Sabin2") & name %in% c("U_mOPV2_cumsum")))

# for paper
fig_cum_emergences_cluster <- 
  ggplot() +
  # geom_ribbon(data = bounds_cluster, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
  geom_ribbon(data = bounds_nOPV2_cluster, aes(x = period, ymin = lower, ymax = upper), fill = "pink", alpha = 0.5, size = 1, linetype = "dashed") +
  # geom_ribbon(data = bounds_Sabin2_cluster, aes(x = period, ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5, size = 1, linetype = "dashed") +
  geom_line(data = temp_cluster, aes(x = period, y = value, color = source, linetype = name), size = 1) +
  geom_vline(xintercept = 2024.5, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.167, color = "red", alpha = 0.25, size = 1) + theme_bw() +
  ylab("Cumulative cVDPV2\nEmergences in AFR") +
  scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_linetype_discrete(name = "Emergences",
                          labels = c("Observed", "Expected based\non Sabin 2 Rate")) +
  # theme(legend.position = c(0.85, 0.3)) +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(5, "in")) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")
ggsave(plot = fig_cum_emergences_cluster,"figures/Cumulative Emergences cluster.png", device = "png", units = "in", width = 7, height = 5)

# Combined figure
plot_layout(fig_cum_doses / fig_cum_emergences_cluster)
ggsave("figures/Cumulative cluster.png", device = "png", units = "in", width = 7, height = 10)

#### Apply time-to-linkage to any known aVDPV2-n's ####
viruses %>% 
  filter(conf_isolate == TRUE) %>%
  group_by(region_who_code) %>%
  # group_by(conf_isolate == TRUE) %>%
  summarize(count = n(),
            min = min(ttl),
            q1 = quantile(ttl, 0.25),
            median = median(ttl),
            mean = mean(ttl),
            q3 = quantile(ttl, 0.75),
            q95 = quantile(ttl, 0.95),
            max = max(ttl))

# CAR
days_since_isolate <- as.numeric(today()-as.Date("2022-12-24"))
  
CAR_ttl <- df_fit_eb %>% 
  # filter(admin0name == "CENTRAL AFRICAN REPUBLIC",
  filter(region_who_code == "AFRO", 
         day <= days_since_isolate) 

CAR_lag <- seq_lag_fit %>% 
  filter(admin0name == "CENTRAL AFRICAN REPUBLIC",
                                  day <= days_since_isolate) %>%
  pull(pseq) # want CDF, so use pseq

CAR_ttl <- CAR_ttl %>%
  mutate(dseq_lag = dseq * rev(CAR_lag), #multiple prob each day against the cumulative proportion detectable by surv lag
         pseq_lag = cumsum(dseq_lag))

tail(CAR_ttl)

# Uganda
# Using AFRO data, 50% at 33 days, 75% at 54 days, 95% at 106
days_since_isolate <- as.numeric(today()-as.Date("2022-02-15"))

UGA_ttl <- df_fit_eb %>% 
  filter(region_who_code == "AFRO", 
         day <= days_since_isolate) 

UGA_lag <- seq_lag_fit %>% 
  filter(admin0name == "UGANDA",
         day <= days_since_isolate) %>%
  pull(pseq) # want CDF, so use pseq

UGA_ttl <- UGA_ttl %>%
  mutate(dseq_lag = dseq * rev(UGA_lag), #multiple prob each day against the cumulative proportion detectable by surv lag
         pseq_lag = cumsum(dseq_lag))

tail(UGA_ttl)

# Apply time-to-linkage to unobserved cVDPV2 seeding risk (dragging unknown further right)

#### Seeding date for novel VDPV2s ####
# DRC SKV-1
viruses %>% filter(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO", index_isolate == TRUE) %>%
  ggplot(aes(x = vdpv_nt_changes_from_sabin)) + geom_histogram() +
  xlab("Index Isolate NT Changes from Sabin") + ggtitle("DRC")

expand_grid(free = c(0,2),t = c(68, 152)/365) %>% 
  mutate(lower = qpois(0.025,free + 9*t), upper = qpois(0.975,free + 9*t))

index <- as.Date("2022-09-27")
SKV1_emerge <- data.frame(date = rev(seq.Date(from = index - 364, to = index, by = "day")), prob = NA)
SKV1_emerge[,"prob"] <- dpois(1:365, 365*4/9)

ggplot(SKV1_emerge, aes(x = date, y = prob)) +
  geom_line() + 
  geom_vline(xintercept = as.Date("2022-04-28")) +
  geom_vline(xintercept = as.Date("2022-07-21")) +
  ggtitle("SKV-1 Estimated Seeding Date\nBased on 6nt changes and Sabin-2 mutation rate")

viruses %>% filter(vdpv_emergence_group_name == "RDC-SKV-1")

# DRC TAN-2
index <- as.Date("2022-11-10")
TAN2_emerge <- data.frame(date = rev(seq.Date(from = index - 364, to = index, by = "day")), prob = NA)
TAN2_emerge[,"prob"] <- dpois(1:365, 365*4/9)

ggplot(TAN2_emerge, aes(x = date, y = prob)) +
  geom_line() + 
  geom_vline(xintercept = as.Date("2022-04-28")) +
  geom_vline(xintercept = as.Date("2022-07-21")) +
  ggtitle("TAN-2 Estimated Seeding Date\nBased on 6nt changes and Sabin-2 mutation rate")

ggplot() +
  geom_line(data = SKV1_emerge, aes(x = date, y = prob), color = "red") + 
  geom_line(data = TAN2_emerge, aes(x = date, y = prob), color = "blue") + 
  geom_vline(xintercept = as.Date("2022-04-28")) +
  geom_vline(xintercept = as.Date("2022-07-21")) +
  ggtitle("SKV-1 (red) and TAN-2 (blue) Estimated Seeding Date")


# CAR aVDPV2
index <- as.Date("2022-12-24")
CAR_emerge <- data.frame(date = rev(seq.Date(from = index - 364, to = index, by = "day")), prob = NA)
CAR_emerge[,"prob"] <- dpois(1:365, 365*5/9)

ggplot(CAR_emerge, aes(x = date, y = prob)) +
  geom_line() + 
  geom_vline(xintercept = as.Date("2022-06-05")) +
  geom_vline(xintercept = as.Date("2022-08-04")) +
  ggtitle("CAR aVDPV2 Estimated Seeding Date\nBased on 7nt changes and Sabin-2 mutation rate")

# DRC Campaigns and Emergences: Seeding Prediction Interval Spread (adapted from Hil)
R1 = ymd(20220428)
R2 = ymd(20220721)

df = tibble(start_date = R1,
            date = seq.Date(R1, today(), by = '1 day')) %>%
  expand_grid(round = c(1,2)) %>%
  mutate(ddiff_R1 = decimal_date(date) - decimal_date(start_date),
         ddiff_R2 = if_else(date > R2, decimal_date(date) - decimal_date(R2), 0),
         nt_diff_R1_025 = 1 + qpois(0.025, 9*ddiff_R1),
         nt_diff_R1_975 = 1 + qpois(0.975, 9*ddiff_R1),
         nt_diff_R2_025 = 1 + qpois(0.025, 9*ddiff_R2),
         nt_diff_R2_975 = 1 + qpois(0.975, 9*ddiff_R2)) 

# DRC Emergences Plot
ggplot() +
  geom_smooth(data = df %>% filter(round == 1), aes(x=date, y=nt_diff_R1_025),linetype = 1, color = "black", method='lm', se=FALSE) +
  geom_smooth(data = df %>% filter(round == 1), aes(x=date, y=nt_diff_R1_975),linetype = 1, color = "black", method='lm', se=FALSE) +
  geom_smooth(data = df %>% filter(round == 2, date > R2), aes(x=date, y=nt_diff_R2_025),linetype = 2, color = "black", method='lm', se=FALSE) +
  geom_smooth(data = df %>% filter(round == 2, date > R2), aes(x=date, y=nt_diff_R2_975),linetype = 2, color = "black", method='lm', se=FALSE) +
  geom_point(data = viruses %>% filter(vdpv_emergence_group_name == "RDC-SKV-1"),
             aes(x = virus_date, y = vdpv_nt_changes_from_sabin), shape = 1) +
  geom_point(data = viruses %>% filter(vdpv_emergence_group_name == "RDC-TAN-2"),
             aes(x = virus_date, y = vdpv_nt_changes_from_sabin), shape = 17) +
  geom_point(data = viruses %>% filter(vdpv_emergence_group_name == "RDC-SKV-1",
                                       epid == "RDC-SKV-UVI-22-010"),
             aes(x = virus_date, y = vdpv_nt_changes_from_sabin), shape = 1, color = "red") +
  geom_point(data = viruses %>% filter(vdpv_emergence_group_name == "RDC-TAN-2",
                                       epid == "RDC-TAN-KBL-22-069"),
             aes(x = virus_date, y = vdpv_nt_changes_from_sabin), shape = 17, color = "red") +
  theme_bw() +
  ylab("NT Changes from Sabin2 VP1") +
  ggtitle("RDC-SKV-1 and RDC-TAN-2")

# CAR Campaigns and VDPV2: Seeding Prediction Interval Spread
R1 = ymd(20220605)
R2 = ymd(20220804)

df = tibble(start_date = R1,
            date = seq.Date(R1, today(), by = '1 day')) %>%
  expand_grid(round = c(1,2)) %>%
  mutate(ddiff_R1 = decimal_date(date) - decimal_date(start_date),
         ddiff_R2 = if_else(date > R2, decimal_date(date) - decimal_date(R2), 0),
         nt_diff_R1_025 = 1 + qpois(0.025, 9*ddiff_R1),
         nt_diff_R1_975 = 1 + qpois(0.975, 9*ddiff_R1),
         nt_diff_R2_025 = 1 + qpois(0.025, 9*ddiff_R2),
         nt_diff_R2_975 = 1 + qpois(0.975, 9*ddiff_R2)) 

ggplot() +
#   geom_step(data = df %>% filter(round == 1), aes(x=date, y=nt_diff_R1_025), col = "grey", method='lm', se=FALSE) +
#   geom_step(data = df %>% filter(round == 1), aes(x=date, y=nt_diff_R1_975), col = "grey", method='lm', se=FALSE) +
#   geom_step(data = df %>% filter(round == 2, date > R2), aes(x=date, y=nt_diff_R2_025), col = "grey", method='lm', se=FALSE) +
#   geom_step(data = df %>% filter(round == 2, date > R2), aes(x=date, y=nt_diff_R2_975), col = "grey", method='lm', se=FALSE) +
  geom_smooth(data = df %>% filter(round == 1), aes(x=date, y=nt_diff_R1_025),linetype=1, col = "black", method='lm', se=FALSE) +
  geom_smooth(data = df %>% filter(round == 1), aes(x=date, y=nt_diff_R1_975),linetype=1, col = "black", method='lm', se=FALSE) +
  geom_smooth(data = df %>% filter(round == 2, date > R2), aes(x=date, y=nt_diff_R2_025),linetype=2, col = "black", method='lm', se=FALSE) +
  geom_smooth(data = df %>% filter(round == 2, date > R2), aes(x=date, y=nt_diff_R2_975),linetype=2, col = "black", method='lm', se=FALSE) +
  geom_point(aes(x = ymd(20221224), y = 7)) +
  ylab("NT Changes from Sabin2 VP1") +
  theme_bw() +
  ggtitle("CAR aVDPV2")

#### Growth of nOPV2 Emergences ####
names(viruses)

# Custom country groupings
viruses$emergence_country <- NA
for (group in viruses$vdpv_emergence_group_name){
  emergence_country = viruses[viruses$vdpv_emergence_group_name %in% c(group) &
                                viruses$index_isolate == TRUE, "admin0name"]
  viruses[viruses$vdpv_emergence_group_name %in% c(group), "emergence_country"] <- emergence_country
}

# Identify active outbreaks
viruses = viruses %>% 
  filter(virus_type_name == "cVDPV2") %>%
  group_by(vdpv_emergence_group_name) %>%
  mutate(most_recent = max(virus_date,na.rm=T),
         active = most_recent >= today()-months(6))

viruses$time_since_index = viruses$virus_date - viruses$index_date

active_groups <- viruses %>% filter(active == TRUE) %>%
  select(vdpv_emergence_group_name) %>%
  unique()
active_groups <- active_groups$vdpv_emergence_group_name

# Dataset for labeling ends
viruses_ends <- viruses %>%
  filter(time_since_index <= 365) %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(desc(virus_date)) %>%
  slice(1) %>%
  ungroup()

# Add a fake ES detection today for all active strains, so that the surveillance lag shows up
viruses_supplemented <- viruses
for (i in active_groups){
  example <- c(viruses_ends[viruses_ends$vdpv_emergence_group_name == i,])
  viruses_supplemented <- rbind(viruses_supplemented, example)
  viruses_supplemented[nrow(viruses_supplemented),]$epid <- "FAKE"
  viruses_supplemented[nrow(viruses_supplemented),]$virus_date <- ymd(today())
  viruses_supplemented[nrow(viruses_supplemented),]$time_since_index <- viruses_supplemented[nrow(viruses_supplemented),]$virus_date - viruses_supplemented[nrow(viruses_supplemented),]$index_date
  viruses_supplemented[nrow(viruses_supplemented),]$surveillance_type_name <- "Environmental"
}

# Calculate cumsum of AFP cases
viruses_supplemented$AFP_count <- 0
viruses_supplemented[viruses_supplemented$surveillance_type_name == "AFP", "AFP_count"] <- 1
viruses_supplemented <- viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(virus_date) %>%
  mutate(AFP_cumsum = cumsum(AFP_count)) %>%
  ungroup()

viruses_supplemented %>% 
  filter(active == TRUE) %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(cases = max(AFP_cumsum),
            source = unique(source),
            most_recent = max(most_recent)) %>%
  arrange(source, -cases)

# Calculate cumsum of ES
viruses_supplemented$ES_count <- 0
viruses_supplemented[viruses_supplemented$surveillance_type_name == "Environmental" &
                       viruses_supplemented$epid != "FAKE", "ES_count"] <- 1
viruses_supplemented <- viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(virus_date) %>%
  mutate(ES_cumsum = cumsum(ES_count)) %>%
  ungroup()

viruses_supplemented %>% 
  filter(active == TRUE) %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(ES = max(ES_cumsum),
            source = unique(source),
            most_recent = max(most_recent)) %>%
  arrange(source, -ES)

# Dataset for labeling ends
viruses_ends <- viruses_supplemented %>%
  filter(time_since_index <= 365,
         epid != "FAKE") %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(desc(virus_date)) %>%
  slice(1) %>%
  ungroup()

facet_names <- c('DEMOCRATIC REPUBLIC OF THE CONGO' = "DRC",
                 'CENTRAL AFRICAN REPUBLIC' = "CAR",
                 'NIGERIA' = "NIGERIA")
# Plot
plot <- 
  ggplot(data = viruses_supplemented %>% 
           filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                                  "CENTRAL AFRICAN REPUBLIC",
                                                                  "NIGERIA"),
                  epid != "FAKE") %>%
           ungroup(),
    aes(x = time_since_index, y = AFP_cumsum, 
                      group = vdpv_emergence_group_name, color = source,
                      alpha = active)) +
    geom_line(size = 1) +
  geom_line(data = viruses_supplemented %>%
              filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                              "CENTRAL AFRICAN REPUBLIC",
                                              "NIGERIA")) %>%
              filter(virus_date > today()-120),
            aes(x = time_since_index, y = AFP_cumsum,
                group = vdpv_emergence_group_name),
            color = "red",
            linetype = "dotted",
            size = 1) +
    scale_x_continuous(limits = c(0, 365), name = "Days Since Index Isolate") +
    scale_y_log10(name = "Cumulative AFP Cases Reported") +
    theme_bw() +
    scale_alpha_discrete(range = c(0.3, 0.9)) +
    scale_color_discrete(labels = c("n-derived", "Sabin-derived"), name = "Source") +
    facet_wrap(.~emergence_country, labeller = as_labeller(facet_names)) +
    theme(legend.position = c(.2, .8))

plot +
  geom_text_repel(aes(label = vdpv_emergence_group_name, color = source),
                  show.legend=FALSE, fontface = "bold",
                  data = viruses_ends  %>%
                    mutate(AFP_cumsum = AFP_cumsum*1.1) %>%
                    filter(source == "nOPV2") %>%
                    filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                    "CENTRAL AFRICAN REPUBLIC",
                                                    "NIGERIA"))) +
  geom_text_repel(aes(label = vdpv_emergence_group_name, color = source),
    show.legend=FALSE,
    data = viruses_ends  %>%
      mutate(AFP_cumsum = AFP_cumsum*1.1) %>%
      filter(active == TRUE) %>%
      filter(source == "Sabin2") %>%
      filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                    "CENTRAL AFRICAN REPUBLIC",
                                                    "NIGERIA"))) +
  guides(alpha = "none")

#### Analyze spread of emergences over time ####
# Measure cumulative number of unique districts affected by emergence group over time
viruses_supplemented$district_count <- 0
viruses_supplemented$province_count <- 0
viruses_supplemented$country_count <- 0

viruses_supplemented <- viruses_supplemented %>%
  arrange(vdpv_emergence_group_name, virus_date)

group <- viruses_supplemented[1, "vdpv_emergence_group_name"]
district_list <- viruses_supplemented[1, "admin2name"]
province_list <- viruses_supplemented[1, "admin1name"]
country_list <- viruses_supplemented[1, "admin0name"]
district_count <- 1
province_count <- 1
country_count <- 1

for(i in 1:nrow(viruses_supplemented)){
  if (viruses_supplemented[i, "vdpv_emergence_group_name"] != group){ # starting a new emergence group
    group <- viruses_supplemented[i, "vdpv_emergence_group_name"]
    district_list <- viruses_supplemented[i, "admin2name"]
    province_list <- viruses_supplemented[i, "admin1name"]
    country_list <- viruses_supplemented[i, "admin0name"]
    district_count <- 1
    province_count <- 1
    country_count <- 1
  } else { #in the same group
    if (!(viruses_supplemented[i, "admin2name"] %in% district_list)){
      district_list <- c(district_list, viruses_supplemented[i, "admin2name"])
      district_count <- district_count + 1
    }
    if (!(viruses_supplemented[i, "admin1name"] %in% province_list)){
      province_list <- c(province_list, viruses_supplemented[i, "admin1name"])
      province_count<- province_count + 1
    }
    if (!(viruses_supplemented[i, "admin0name"] %in% country_list)){
      country_list <- c(country_list, viruses_supplemented[i, "admin0name"])
      country_count<- country_count + 1
    }
  }
  viruses_supplemented[i, "district_count"] <- district_count
  viruses_supplemented[i, "province_count"] <- province_count
  viruses_supplemented[i, "country_count"] <- country_count
  
  if (i %% 100 == 0){cat(".")}
}

# Dataset for labeling ends
viruses_ends <- viruses_supplemented %>%
  filter(time_since_index <= 365,
         epid != "FAKE") %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(desc(virus_date)) %>%
  slice(1) %>%
  ungroup()

# Plot for districts
plot <- 
  ggplot(data = viruses_supplemented %>% 
           filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                           "CENTRAL AFRICAN REPUBLIC",
                                           "NIGERIA"),
                  epid != "FAKE") %>%
           ungroup(),
         aes(x = time_since_index, y = district_count, 
             group = vdpv_emergence_group_name, color = source,
             alpha = active)) +
  geom_line(size = 1) +
  geom_line(data = viruses_supplemented %>%
              filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                              "CENTRAL AFRICAN REPUBLIC",
                                              "NIGERIA")) %>%
              filter(virus_date > today()-120),
            aes(x = time_since_index, y = district_count,
                group = vdpv_emergence_group_name),
            color = "red",
            linetype = "dotted",
            size = 1) +
  scale_x_continuous(limits = c(0, 365), name = "Days Since Index Isolate") +
  scale_y_log10(limits = c(1, 100), name = "Cumulative Districts Affected") +
  theme_bw() +
  scale_alpha_discrete(range = c(0.3, 0.9)) +
  scale_color_discrete(labels = c("n-derived", "Sabin-derived"), name = "Source") +
  facet_wrap(.~emergence_country, labeller = as_labeller(facet_names)) +
  theme(legend.position = c(.2, .8))

plot +
  geom_text_repel(aes(label = vdpv_emergence_group_name, color = source),
                  show.legend=FALSE, fontface = "bold",
                  data = viruses_ends  %>% 
                    filter(source == "nOPV2") %>%
                    filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                    "CENTRAL AFRICAN REPUBLIC",
                                                    "NIGERIA"))) +
  geom_text_repel(aes(label = vdpv_emergence_group_name, color = source),
                  show.legend=FALSE,
                  data = viruses_ends  %>%
                    filter(active == TRUE) %>%
                    filter(source == "Sabin2") %>%
                    filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                    "CENTRAL AFRICAN REPUBLIC",
                                                    "NIGERIA"))) +
  guides(alpha = "none")

# Plot for provinces
plot <- 
  ggplot(data = viruses_supplemented %>% 
           filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                           "CENTRAL AFRICAN REPUBLIC",
                                           "NIGERIA"),
                  epid != "FAKE") %>%
           ungroup(),
         aes(x = time_since_index, y = province_count, 
             group = vdpv_emergence_group_name, color = source,
             alpha = active)) +
  geom_line(size = 1) +
  geom_line(data = viruses_supplemented %>%
              filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                              "CENTRAL AFRICAN REPUBLIC",
                                              "NIGERIA")) %>%
              filter(virus_date > today()-120),
            aes(x = time_since_index, y = province_count,
                group = vdpv_emergence_group_name),
            color = "red",
            linetype = "dotted",
            size = 1) +
  scale_x_continuous(limits = c(0, 365), name = "Days Since Index Isolate") +
  scale_y_continuous(limits = c(1, 20), name = "Cumulative Provinces Affected") +
  theme_bw() +
  scale_alpha_discrete(range = c(0.3, 0.9)) +
  scale_color_discrete(labels = c("n-derived", "Sabin-derived"), name = "Source") +
  facet_wrap(.~emergence_country, labeller = as_labeller(facet_names)) +
  theme(legend.position = c(.2, .8))

plot +
  geom_text_repel(aes(label = vdpv_emergence_group_name, color = source),
                  show.legend=FALSE, fontface = "bold",
                  data = viruses_ends  %>% 
                    filter(source == "nOPV2") %>%
                    filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                    "CENTRAL AFRICAN REPUBLIC",
                                                    "NIGERIA"))) +
  geom_text_repel(aes(label = vdpv_emergence_group_name, color = source),
                  show.legend=FALSE,
                  data = viruses_ends  %>%
                    filter(active == TRUE) %>%
                    filter(source == "Sabin2") %>%
                    filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                    "CENTRAL AFRICAN REPUBLIC",
                                                    "NIGERIA"))) +
  guides(alpha = "none")

# Plot for countries
plot <- 
  ggplot(data = viruses_supplemented %>% 
           filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                           "CENTRAL AFRICAN REPUBLIC",
                                           "NIGERIA"),
                  epid != "FAKE") %>%
           ungroup(),
         aes(x = time_since_index, y = country_count, 
             group = vdpv_emergence_group_name, color = source,
             alpha = active)) +
  geom_line(size = 1) +
  geom_line(data = viruses_supplemented %>%
              filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                              "CENTRAL AFRICAN REPUBLIC",
                                              "NIGERIA")) %>%
              filter(virus_date > today()-120),
            aes(x = time_since_index, y = country_count,
                group = vdpv_emergence_group_name),
            color = "red",
            linetype = "dotted",
            size = 1) +
  scale_x_continuous(limits = c(0, 365), name = "Days Since Index Isolate") +
  scale_y_continuous(limits = c(1, 5), name = "Cumulative Countries Affected") +
  theme_bw() +
  scale_alpha_discrete(range = c(0.3, 0.9)) +
  scale_color_discrete(labels = c("n-derived", "Sabin-derived"), name = "Source") +
  facet_wrap(.~emergence_country, labeller = as_labeller(facet_names)) +
  theme(legend.position = c(.2, .8))

plot +
  geom_text_repel(aes(label = vdpv_emergence_group_name, color = source),
                  show.legend=FALSE, fontface = "bold",
                  data = viruses_ends  %>% 
                    filter(source == "nOPV2") %>%
                    filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                    "CENTRAL AFRICAN REPUBLIC",
                                                    "NIGERIA"))) +
  geom_text_repel(aes(label = vdpv_emergence_group_name, color = source),
                  show.legend=FALSE,
                  data = viruses_ends  %>%
                    filter(active == TRUE) %>%
                    filter(source == "Sabin2") %>%
                    filter(emergence_country %in% c("DEMOCRATIC REPUBLIC OF THE CONGO",
                                                    "CENTRAL AFRICAN REPUBLIC",
                                                    "NIGERIA"))) +
  guides(alpha = "none")

# Histogram of number of affected districts by emergence group
viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(districts = max(district_count),
            source = unique(source),
            active = unique(active)) %>%
  mutate(source.active = ifelse(source == "Sabin2" & active == "TRUE", "Sabin2_Active",
                                ifelse(source == "Sabin2" & active == "FALSE", "Sabin2_Inactive",
                                       ifelse(source == "nOPV2" & active == "TRUE", "nOPV2_Active",
                                              "nOPV2_Inactive")))) %>%
  ggplot(aes(x = districts, fill = source.active)) +
  geom_histogram( bins = 5) +
  facet_grid(source~.) +
  scale_x_log10(name = "Affected Districts") +
  scale_y_continuous(name = "Number of Emergence Groups") +
  scale_fill_manual(values = c("#F8766D50", "#F8766D", "#00BFC450", "#00BFC4"),
                    name = "") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Histogram of number of affected provinces by emergence group
viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(provinces = max(province_count),
            source = unique(source),
            active = unique(active)) %>%
  mutate(source.active = ifelse(source == "Sabin2" & active == "TRUE", "Sabin2_Active",
                                ifelse(source == "Sabin2" & active == "FALSE", "Sabin2_Inactive",
                                       ifelse(source == "nOPV2" & active == "TRUE", "nOPV2_Active",
                                              "nOPV2_Inactive")))) %>%
  ggplot(aes(x = provinces, fill = source.active)) +
  geom_histogram( bins = 5) +
  facet_grid(source~.) +
  scale_x_log10(name = "Affected Provinces") +
  scale_y_continuous(name = "Number of Emergence Groups") +
  scale_fill_manual(values = c("#F8766D50", "#F8766D", "#00BFC450", "#00BFC4"),
                    name = "") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Histogram of number of affected countries by emergence group
viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(
    countries = max(country_count),
    source = unique(source),
    active = unique(active)
  ) %>%
  mutate(source.active = ifelse(source == "Sabin2" & active == "TRUE", "Sabin2_Active",
                                ifelse(source == "Sabin2" & active == "FALSE", "Sabin2_Inactive",
                                       ifelse(source == "nOPV2" & active == "TRUE", "nOPV2_Active",
                                              "nOPV2_Inactive")))) %>%
  ggplot(aes(x = countries, fill = source.active)) +
  geom_histogram( bins = 5) +
  facet_grid(source~.) +
  scale_x_log10(name = "Affected Countries") +
  scale_y_continuous(name = "Number of Emergence Groups") +
  scale_fill_manual(values = c("#F8766D50", "#F8766D", "#00BFC450", "#00BFC4"),
                    name = "") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Summarize for countries
viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(countries = max(country_count),
            source = unique(source),
            active = unique(active)) %>%
  mutate(source.active = ifelse(source == "Sabin2" & active == "TRUE", "Sabin2_Active",
                                ifelse(source == "Sabin2" & active == "FALSE", "Sabin2_Inactive",
                                       ifelse(source == "nOPV2" & active == "TRUE", "nOPV2_Active",
                                              "nOPV2_Inactive")))) %>%
  ungroup() %>%
  group_by(source, countries == 1) %>%
  group_by(source.active, countries == 1) %>%
  summarize(n = n())

#### Assess predictors of initial outbreak growth ####
# Plot distribution of final size over/under 10 AFP cases
viruses_final <- viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(desc(AFP_cumsum)) %>%
  slice(1) %>%
  ungroup()

viruses_final %>%
  group_by(AFP_cumsum >= 10, source) %>%
  # filter(active == FALSE) %>%
  summarize(count = n())

viruses_final %>%
  mutate(cases = ifelse(AFP_cumsum>=10, "10+ Cases", "<10 Cases")) %>%
  ggplot(aes(x = cases)) +
  geom_bar(aes(fill = source), stat="count", width = 1, position = "dodge") +
  # geom_density() +
  xlab("") +
  ylab("Number of Emergence groups") 
  # scale_x_binned(v = c("<10 AFP cases", "10 or more AFP cases")) +
  # facet_grid(.~active) +
  # scale_fill_discrete(labels = c("n-derived", "Sabin-derived"), name = "Source")

# Plot distribution of time-to-first-10 AFP cases
viruses_supplemented %>%
  filter(AFP_cumsum == 10,
         surveillance_type_name == "AFP") %>%
  ggplot(aes(x = time_since_index)) +
  geom_histogram(aes(fill = source) ) +
  # geom_density() +
  xlab("Days since Index Isolate") +
  ylab("Emergence groups reaching 10 AFP cases by day X")

# Join with immunity data at time of index
names(immunity_u5_data)
immunity_u5_data_temp <- immunity_u5_data %>% 
  group_by(adm0_name, adm1_name, period) %>%
  summarize(immunity_u5_weighted = weighted.mean(immunity_u5, w = target_pop, na.rm=T))
names(immunity_u5_data_temp)[c(1, 2, 3)] <- c("admin0name", "admin1name", "period_index")

viruses_supplemented <- viruses_supplemented %>% ungroup() %>%
  mutate(period_index = year(viruses_supplemented$index_date) + (month(viruses_supplemented$index_date)-1)/12) %>%
  left_join(immunity_u5_data_temp, by = c("period_index", "admin1name", "admin0name")) %>%
  mutate(time_since_index_numeric = as.numeric(time_since_index))

viruses_supplemented %>% 
  filter(AFP_cumsum == 10) %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(desc(AFP_cumsum)) %>%
  slice(1) %>%
  ungroup() %>%
  # filter(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
ggplot(aes(x = immunity_u5_weighted, y = time_since_index_numeric)) +
  geom_smooth(method = "lm", color = "black") +
    geom_label(aes(color = source, label = vdpv_emergence_group_name), nudge_y = -10, alpha = 0.5) +
  geom_point(aes(color = source), size = 2) +
  scale_color_discrete(labels = c("n-derived", "Sabin-derived"), name = "Source") +
  xlab("Immunity in province of emergence") +
  ylab("Days until 10 AFP cases") +
  theme_bw()

# Scatterplot matrix
viruses_supplemented %>% 
  filter(AFP_cumsum == 10) %>%
  filter(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  select(time_since_index_numeric, immunity_u5_weighted) %>%
  ggpairs()

# Histogram cases per emergence, noting which are active
viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(cases = max(AFP_cumsum),
            source = unique(source),
            active = unique(active)) %>%
  ggplot(aes(x = cases, fill = source)) +
  geom_histogram(position = "dodge", bins = 5) +
  scale_x_log10(name = "Cases") +
  scale_y_continuous(name = "Number of Emergence Groups") +
  theme(panel.grid = element_blank())

# Histogram cases per emergence, noting which are active (option 2)
viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(cases = max(AFP_cumsum),
            source = unique(source),
            active = unique(active)) %>%
  mutate(source.active = ifelse(source == "Sabin2" & active == "TRUE", "Sabin2_Active",
                                ifelse(source == "Sabin2" & active == "FALSE", "Sabin2_Inactive",
                                       ifelse(source == "nOPV2" & active == "TRUE", "nOPV2_Active",
                                              "nOPV2_Inactive")))) %>%
  ggplot(aes(x = cases, fill = source.active)) +
  geom_histogram( bins = 5) +
  facet_grid(source~.) +
  scale_x_log10(name = "Cases") +
  scale_y_continuous(name = "Number of Emergence Groups") +
  scale_fill_manual(values = c("#F8766D50", "#F8766D", "#00BFC450", "#00BFC4"),
                    name = "") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Histogram ES per emergence, noting which are active (option 2)
viruses_supplemented %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(ES = max(ES_cumsum),
            source = unique(source),
            active = unique(active)) %>%
  mutate(source.active = ifelse(source == "Sabin2" & active == "TRUE", "Sabin2_Active",
                                ifelse(source == "Sabin2" & active == "FALSE", "Sabin2_Inactive",
                                       ifelse(source == "nOPV2" & active == "TRUE", "nOPV2_Active",
                                              "nOPV2_Inactive")))) %>%
  ggplot(aes(x = ES, fill = source.active)) +
  geom_histogram( bins = 5) +
  facet_grid(source~.) +
  scale_x_log10(name = "ES Detections") +
  scale_y_continuous(name = "Number of Emergence Groups") +
  scale_fill_manual(values = c("#F8766D50", "#F8766D", "#00BFC450", "#00BFC4"),
                    name = "") +
  theme_bw() +
  theme(panel.grid = element_blank())

#### Look into inter-campaign interval ####
# View(polis_pops_full)
sias_provinces <- polis_pops_full %>%
  group_by(adm0_name, adm1_name, start_date) %>%
  summarize(adm0_name = unique(adm0_name),
            adm1_name = unique(adm1_name),
            start_date = unique(start_date),
            vaccinetype = unique(vaccinetype),
            target_pop_total = unique(target_pop_total),
            region = unique(region),
            parentactivitycode = list(unique(parentactivitycode)),
            source = unique(source))
# View(sias_provinces)

# Calculate months since previous SIA
sias_provinces$prev_sia <- as.Date("9999-01-01")
sias_provinces <- sias_provinces %>% arrange(adm0_name, adm1_name, start_date)

for (i in 2:nrow(sias_provinces)){ #Manual, slow.
  if (sias_provinces[i-1,"adm1_name"] == sias_provinces[i,"adm1_name"]){
    sias_provinces[i,"prev_sia"] <- sias_provinces[i-1,"start_date"]
  }
}

# Identify SIAs where there is more than 12 months since the previous. Call these "R1"
sias_provinces[sias_provinces$prev_sia == "9999-01-01", "prev_sia"] <- NA
sias_provinces$interval <- days(sias_provinces$start_date - sias_provinces$prev_sia)$day
sias_provinces$round_1 <- sias_provinces$interval > 365

# Identify the next SIA after an "R1"
rows_round_1 <- which(sias_provinces$round_1 == TRUE)
rows_round_1 <- rows_round_1[1:length(rows_round_1)-1] # remove the last row, since it can't have a second round
sias_provinces$round_2 <- FALSE
for (i in rows_round_1){ #Manual, kinda slow.
  if (sias_provinces[i,"adm1_name"] == sias_provinces[i+1,"adm1_name"]){  # same province
    if (sias_provinces[i+1, "round_1"] == FALSE){  # not a first round of its own
        sias_provinces[i+1, "round_2"] <- TRUE
    }
  }
}
sias_provinces %>% filter(adm0_name == "NIGERIA", start_date > "2016-05-01") %>% View() 

# priority countries
sias_provinces$adm0_group <- "Others - non AFRO"
sias_provinces[sias_provinces$region == "AFRO", "adm0_group"] <- "Others - AFRO"
sias_provinces[sias_provinces$adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO", "adm0_group"] <- "DRC"
sias_provinces[sias_provinces$adm0_name == "CENTRAL AFRICAN REPUBLIC", "adm0_group"] <- "CAR"
sias_provinces[sias_provinces$adm0_name == "NIGERIA", "adm0_group"] <- "NIGERIA"
sias_provinces[sias_provinces$adm0_name == "ETHIOPIA", "adm0_group"] <- "ETHIOPIA"

# Plot distribution of interval to R2 over time and by geo and by vaccinetype
sias_provinces %>% filter(round_2 == TRUE,
                          start_date > "2016-05-01") %>%
  ggplot(aes(x = start_date, y = interval/30.4, color = vaccinetype, shape = vaccinetype)) +
    geom_hline(yintercept = 1, color = "darkgray") +
    geom_smooth(method = "lm", se = FALSE) +
    geom_jitter(size = 2, width = .25, height = .25) +
    facet_wrap(adm0_group~.) +
    ggtitle("Post-switch Type-2 response inter-campaign intervals") +
    xlab("R2 Start Date") +
    scale_y_continuous(name = "Months between R1 and R2", breaks = seq(1, 12, 2))

sias_provinces %>% filter(round_2 == TRUE,
                          start_date > "2016-05-01") %>%
  group_by(adm0_group, vaccinetype) %>%
  summarize(count = n(),
            median = median(interval/30.4),
            low = quantile(interval/30.4, .25),
            high = quantile(interval/30.4, .75))

# Compare suspected seeding rounds to rest
sias_provinces %>% filter(round_2 == TRUE,
                          start_date != "2022-08-04", # CAR rounds
                          start_date != "2021-04-24", # Nigeria rounds
                          start_date != "2022-07-21", # DRC rounds
                          vaccinetype == "nOPV2") %>%
  group_by(region) %>%
  summarize(count = n(),
            median = median(interval/7),
            low = quantile(interval/7, .25),
            high = quantile(interval/7, .75))

#### Summary of emergence groups over time ####
# Find order of vdpv emergence groups
group_order <- viruses %>% 
  group_by(vdpv_emergence_group_name) %>%
  summarize(max_date = max(virus_date)) %>%
  arrange(max_date) %>%
  pull(vdpv_emergence_group_name) %>%
  unique()

plot <-
viruses %>%
  arrange(virus_date) %>%
  filter(region_who_code == "AFRO") %>%
  filter(surveillance_type_name %in% c("AFP", "Environmental")) %>%
  ggplot(aes(x = virus_date, 
    y = factor(vdpv_emergence_group_name, levels = rev(group_order)), 
    color = surveillance_type_name)) + 
  geom_point(shape = "O", size = 3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", minor_breaks = waiver()) +
  ylab("Emergence Group") +
  xlab("Virus Date") 
plot
ggsave(filename = "C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio/figures/figure_emergence_groups_over_time.png", plot, width = 10, height = 7, units = "in", dpi = 300)

# number of emergences
viruses_supplemented %>% 
  filter(index_isolate == TRUE) %>% 
  filter(virus_date > "2016-05-01") %>%
  filter(!(admin0name %in% c("AFGHANISTAN", "PAKISTAN"))) %>% View()

# Plot by country
viruses$vdpv_emergence_group_cluster <- "Other"
viruses[viruses$vdpv_emergence_group_name == "NIE-JIS-1", "vdpv_emergence_group_cluster"] <- "NIE-JIS-1"
viruses[viruses$vdpv_emergence_group_name == "NIE-ZAS-1", "vdpv_emergence_group_cluster"] <- "NIE-ZAS-1"
viruses[viruses$vdpv_emergence_group_name == "CHA-NDJ-1", "vdpv_emergence_group_cluster"] <- "CHA-NDJ-1"

country_order <- viruses %>% 
  group_by(admin0name) %>%
  summarize(max_date = max(virus_date)) %>%
  arrange(max_date) %>%
  pull(admin0name) %>%
  unique()

viruses %>%
  arrange(admin0name) %>%
  # filter(region_who_code == "AFRO") %>%
  filter(surveillance_type_name %in% c("AFP", "Environmental")) %>%
  ggplot(aes(x = virus_date, 
             y = factor(admin0name, levels = country_order), 
             color = vdpv_emergence_group_cluster)) + 
  geom_point( size = 3, alpha = .5) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", minor_breaks = waiver()) +
  ylab("Country") +
  xlab("Virus Date")

#### Add mOPV2-derived emergences to analysis ####
# Read data provided by Elizabeth
data_sabin_emergence <- read.csv("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio/emergence_probs_by_SIA.csv")
names(data_sabin_emergence)
names(data_sabin_emergence)[3] <- "parentactivitycode"

# Generate cumulative culpability probability for each campaign
data_sabin_emergence_cumulative <- data_sabin_emergence %>%
  group_by(parentactivitycode) %>%
  summarize(emergence_prob = sum(summed_probabilities))

sias_included <- data_sabin_emergence_cumulative$parentactivitycode[which(unique(data_sabin_emergence_cumulative$parentactivitycode) %in% unlist(sias_provinces$parentactivitycode))]

sias_provinces$emerge_data <- FALSE
sias_provinces$emerge_data_sia <- "NA"
for (i in sias_included){
  row_match <- grep(i, sias_provinces$parentactivitycode)
  sias_provinces[row_match, "emerge_data"] <- TRUE
  sias_provinces[row_match, "emerge_data_sia"] <- i
}

sias_provinces <- sias_provinces %>%
  left_join(data_sabin_emergence_cumulative, by = c("emerge_data_sia" = "parentactivitycode"))

# Plot distribution of interval to R2 over time and by geo and by vaccinetype
sias_provinces %>% filter(round_2 == TRUE,
                          emerge_data == TRUE) %>%
    mutate(adm0_group_new = case_when(adm0_name == "NIGERIA" ~ "Nigeria",
                                      adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO" ~ "DRC",
                                      adm0_group %in% c("Others - AFRO", "CAR", "ETHIOPIA", "Others - non AFRO") ~ "Others")) %>%
    group_by(emerge_data_sia) %>%
    summarize(emerge_data_sia = unique(emerge_data_sia),
              emergence_prob = unique(emergence_prob),
              interval = unique(interval),
              adm0_group_new = unique(adm0_group_new)) %>%
  ggplot(aes(x = interval/7, y = emergence_prob)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm") +
    facet_wrap(adm0_group_new~.) +
    ggtitle("Longer R1-to-R2 intervals correlated with higher probability of seeding") +
    labs(subtitle = "Sabin-2 campaigns 2016-2019",
         caption = "Emergence probability estimates provided by Gray et al 2023") +
    scale_x_continuous(name = "Weeks between R1 and R2", limits = c(0,NA)) +
    scale_y_continuous(name = "Estimated number of emergences by campaign")
ggsave(filename = "C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio/figures/figure_sabin_emergences_by_interval.png", width = 7, height = 4, units = "in", dpi = 300)

temp <- 
  sias_provinces %>% 
    mutate(adm0_group_new = case_when(adm0_name == "NIGERIA" ~ "Nigeria",
                                                       adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO" ~ "DRC",
                                                       adm0_group %in% c("Others - AFRO", "CAR", "ETHIOPIA", "Others - non AFRO") ~ "Others")) %>%
    filter(
           round_2 == TRUE,
           adm0_group_new == "Nigeria",
           emerge_data == TRUE) %>%
    group_by(emerge_data_sia) %>%
    summarize(emerge_data_sia = unique(emerge_data_sia),
              emergence_prob = unique(emergence_prob),
              interval = unique(interval),
              adm0_group_new = unique(adm0_group_new)) %>%
    group_by(adm0_group_new)
cor.test(temp$interval, temp$emergence_prob, method = "spearman")

sias_provinces %>% 
  mutate(adm0_group_new = case_when(adm0_name == "NIGERIA" ~ "Nigeria",
                                    adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO" ~ "DRC",
                                    adm0_group %in% c("Others - AFRO", "CAR", "ETHIOPIA", "Others - non AFRO") ~ "Others")) %>%
  filter(round_2 == TRUE, 
         emerge_data == TRUE) %>%
  group_by(emerge_data_sia) %>%
  summarize(emerge_data_sia = unique(emerge_data_sia),
            emergence_prob = unique(emergence_prob),
            interval = unique(interval),
            adm0_group_new = unique(adm0_group_new)) %>%
  group_by(adm0_group_new) %>%
  summarize(cor(interval, emergence_prob, method = "spearman"))

#### Summarize scope of response ####
sias_provinces %>% # Not working yet
  filter(round_2 == TRUE) %>% #focus on just R2, since that would avoid mixing R0 and R1
  group_by(parentactivitycode) %>%
  mutate(target_pop_total = sum(target_pop_total, na.rm=T)) %>%
  group_by(adm0_group, start_date > "2021-01-01") %>%
  summarize(count = n(),
            median = median(target_pop_total),
            low = quantile(target_pop_total, .25),
            high = quantile(target_pop_total, .75))



#### Save workspace locally ####
save.image(file = "VDPV2n_analyses.RData")

  