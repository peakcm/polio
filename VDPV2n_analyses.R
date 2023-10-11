# VDPV2-n analyses
rm(list=ls())

#### Load Libraries ####
library(ggridges)
library(tidyverse)
library(lubridate)
library(PolisAPI)
library(ggrepel)
library(GGally)
library(ggh4x)

#### Load workspace ####
load("VDPV2n_analyses.RData")

#### Helper functions ####
#Create region from ADM_0
Func_region = function(ADM_0){
  EMRO <- c("AFGHANISTAN", "EGYPT", "IRAQ", "PAKISTAN", "SOMALIA", "SUDAN", 
            "DJIBOUTI", "ERITREA", "JORDAN",
            "IRAN (ISLAMIC REPUBLIC OF)", "LIBYA", "SYRIAN ARAB REPUBLIC", 
            "YEMEN","LEBANON", "LIBYA", "OCCUPIED PALESTINIAN TERRITORY, INCLUDING EAST JERUSALEM")
  WPRO <- c("PHILIPPINES", "MALAYSIA","LAO PEOPLE'S DEMOCRATIC REPUBLIC")
  SEARO <- c("INDONESIA", "MYANMAR", "INDIA", "NEPAL")
  EURO <- c("TAJIKISTAN", "GEORGIA", "RUSSIAN FEDERATION", "UKRAINE")
  ifelse(ADM_0 %in% EMRO, "EMRO",
         ifelse(ADM_0 %in% WPRO, "WPRO",
                ifelse(ADM_0 %in% SEARO, "SEARO",
                       ifelse(ADM_0 %in% EURO, "EURO",
                              "AFRO"))))
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
polis_pops %>% group_by(adm0_name) %>% summarize(region = unique(region)) %>%View()

# Clean data
polis_pops <- polis_pops %>% filter(start_date < today())
# View(polis_pops %>% filter(vaccinetype == "nOPV2", status == "Planned"))

# Mark campaigns that did happen:
# done_campaigns <- c("MLI-2023-002", "BDI-2023-002")
# polis_pops[polis_pops$parentactivitycode %in% done_campaigns, "status"] <- "Done"
polis_pops <- polis_pops %>% filter(status == "Done")

# Manually remove campaigns that didn't happen (yet)
# not_campaigns <- c( "COD-2023-003") 
# polis_pops %>% filter(parentactivitycode %in% not_campaigns) %>% View()
# polis_pops[polis_pops$parentactivitycode %in% not_campaigns, "target_pop"] <- 0

# Count vaccine used
polis_pops %>% ungroup() %>% filter(vaccinetype == "nOPV2")  %>% summarize(pop= sum(target_pop, na.rm = T))

# Add pop to campaigns missing it
# polis_pops %>% filter(region == "AFRO", is.na(target_pop)) %>% View()
# eth_rows <- polis_pops %>% filter(parentactivitycode == "ETH-2021-003") %>% nrow()
# polis_pops[polis_pops$parentactivitycode == "ETH-2021-003", "target_pop"] <- 16712725/eth_rows

# check target pop compared with tracker
polis_pops %>% group_by(adm0_name) %>% 
  filter(vaccinetype == "nOPV2") %>%
  summarize(target_pop = sum(target_pop, na.rm=T)) %>% View()
# Lacking in Indonesia. Too many in Nigeria (440349303 vs 397787446)

# !!!!!!! MANUAL !!!!!! Scale down Nigeria by 0.903 (397787446 / 440349303)
polis_pops[polis_pops$adm0_name == "NIGERIA" & polis_pops$vaccinetype == "nOPV2", "target_pop"] %>% sum()
polis_pops[polis_pops$adm0_name == "NIGERIA" & polis_pops$vaccinetype == "nOPV2", "target_pop"] <- 
  polis_pops[polis_pops$adm0_name == "NIGERIA" & polis_pops$vaccinetype == "nOPV2", "target_pop"] * 0.903

# !!!!!!! MANUAL !!!!!! Distribute 12m to Indonesia equally across districts
indo_rows <- polis_pops %>% filter(vaccinetype=="nOPV2", adm0_name == "INDONESIA") %>% nrow()
polis_pops[polis_pops$adm0_name == "INDONESIA" & polis_pops$vaccinetype=="nOPV2", "target_pop"] <- 12415310/indo_rows

# Gut check nOPV2 usage
polis_pops %>% filter(vaccinetype == "nOPV2") %>% summarize(pop = sum(target_pop, na.rm=T))

# vaccine usage post-switch
polis_pops %>% filter(region %in% c("AFRO", "EMRO", "EURO"), start_date > "2016-05-01") %>% group_by(source) %>%
  summarize(pop = sum(target_pop, na.rm=T))

polis_pops %>%
  # filter(region %in% c("AFRO", "EMRO", "EURO")) %>%
  filter(start_date > "2016-05-01") %>%
  group_by(vaccinetype) %>%
  # group_by(region) %>%
  summarize(pop = sum(target_pop, na.rm=T))

polis_pops %>%
  filter(start_date > "2023-04-01", start_date < "2023-10-01") %>%
  group_by(vaccinetype) %>%
  summarize(pop = sum(target_pop, na.rm=T))

polis_pops %>%
  filter(start_date > "2016-05-01") %>%
  group_by(vaccinetype, region) %>%
  summarize(pop = sum(target_pop, na.rm=T))

polis_pops %>%
  filter(start_date > "2016-05-01", adm0_name == "NIGERIA") %>%
  group_by(vaccinetype, region) %>%
  summarize(pop = sum(target_pop, na.rm=T))

# Compare campaign size by vaccine type
polis_pops %>% 
  group_by(source, parentactivitycode) %>% 
  # filter(region %in% c("AFRO", "EMRO", "EURO")) %>%
  filter(region %in% c("AFRO")) %>%
  filter(start_date > "2016-05-01") %>%
  summarize(sum = sum(target_pop, na.rm=T)) %>%
  ungroup() %>%
  group_by(source) %>%
  summarize(count = n(),
            min = min(sum),
            q1 = quantile(sum, 0.25),
            median = median(sum),
            mean = mean(sum),
            q3 = quantile(sum, 0.75),
            max = max(sum))

temp <- 
  polis_pops %>%
  filter(region %in% c("AFRO")) %>%
  filter(start_date > "2016-05-01") %>%
  group_by(source, parentactivitycode) %>%
  summarize(sum = sum(target_pop, na.rm=T)) 
wilcox.test(temp$sum ~ temp$source,
            exact = FALSE)


# Simplify fields
polis_pops_full <- polis_pops
polis_pops <- polis_pops %>% 
  filter(start_date > "2016-05-01") %>%
  select(adm0_name, adm1_name, adm2_name, parentactivitycode,  childactivitycode, 
                                    start_date, vaccinetype, GUID,
                                    target_pop, period, quarter, source, region)

sias_figure <- polis_pops %>% group_by(quarter, source) %>% summarize(target_pop = sum(target_pop, na.rm=T))
ggplot(sias_figure, aes(x = quarter, y = target_pop, color = source)) +  geom_line()

# Plot number of doses
polis_pops %>%
  filter(!is.na(target_pop)) %>%
  ungroup() %>%
  group_by(source, period) %>%
  reframe(period = unique(period),
            doses = sum(target_pop)) %>%
  group_by(source) %>%
  arrange(period) %>%
  reframe(period = unique(period),
            doses_cumsum = cumsum(doses)) %>%
  ggplot() +
    geom_vline(xintercept = 2023.75, alpha = 0.25, size = 1) +
    geom_vline(xintercept = 2021.167, alpha = 0.25, color = "red", size = 1) +
    geom_line(aes(x = period, y = doses_cumsum/1e6, color = source), size = 1) +
    theme_bw() +
    scale_x_continuous(limits = c(2016.0, 2027), breaks = seq(2016, 2027, 2)) +
    scale_y_continuous(name = "Cumulative Doses\n(Million)") +
    scale_color_discrete(name = "Vaccine Type") +
    force_panelsizes(rows = unit(2, "in"),
                     cols = unit(5, "in"))
ggsave("figures/Cumulative Doses.png", device = "png", units = "in", width = 7, height = 3)

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

# Check emergence group names
novel_emergences <- c("RDC-SKV-1", "RDC-TAN-2", "RDC-KOR-1", "CAF-KEM-1", "NIE-KBS-1", "RDC-HKA-2","CAF-BNG-3", "BOT-FRA-1", "EGY-NOR-1")
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
  ungroup() %>% 
  group_by(admin0name, source) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = reorder(admin0name, -count), y = count, fill = source)) +
    # facet_grid(source~.) +
    geom_col() + coord_flip() + 
    xlab("") + ylab("Post-Switch cVDPV2 Emergences")

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
# Run immunity_calc_eag.R to produce new U5 immunity estimates to use here
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

immunity_u5_data <- immunity_u5_data %>% filter(period >= 2016, period < 2023.5)

# Add admin1 and admin2 to immunity_u5_data
immunity_u5_data <- immunity_u5_data %>% 
  left_join(polis_pops %>% select(GUID, adm0_name, adm1_name,adm2_name, target_pop) %>%
              group_by(GUID) %>% filter(row_number() == 1), 
            by = "GUID")

#### Combine SIA and immunity datasets ####
# Join data on GUID, and period, noting that immunity from campaign doesn't impact that period yet
data <- polis_pops %>% 
  left_join(immunity_u5_data %>% select(GUID, period, week, immunity_u5), by = c("GUID", "period"))

data$population_total <- data$target_pop / 0.17

data %>% group_by(vaccinetype) %>% summarize(target_pop = sum(target_pop, na.rm=T))

# add region
unique(data$adm0_name)
data$Region <- sapply(data$adm0_name, Func_region)

# See missingness
# data %>% filter(is.na(immunity_u5)) %>% View()
# missing immunity data from Indonesia, Phillipines, and Malaysia. Exclude from analysis


# Convert in province-level campaigns (sum across Admin2)
# When the period and the admin1 name are the same, then sum the population and create a pop-weighted-immunity estimate for the province
data_province <- data %>% group_by(adm1_name, period) %>%
  summarize(target_pop_sum = sum(target_pop), 
            population_total_sum = sum(population_total),
            immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
            vaccinetype = unique(vaccinetype),
            Region = unique(Region),
            adm0_name = unique(adm0_name),
            source = unique(source),
            start_date = min(start_date),
            quarter = unique(quarter),
            week = unique(week))
sum(data_province[data_province$vaccinetype == "nOPV2", "target_pop_sum"], na.rm=T)
data_province <- data_province %>% filter(is.na(adm1_name) == F)

data_national <- data %>% group_by(adm0_name, period) %>%
  summarize(target_pop_sum = sum(target_pop), 
            population_total_sum = sum(population_total),
            immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
            vaccinetype = unique(vaccinetype),
            Region = unique(Region),
            source = unique(source),
            quarter = unique(quarter))
sum(data_national[data_national$vaccinetype == "nOPV2", "target_pop_sum"], na.rm=T)
sum(data_national[data_national$vaccinetype == "mOPV2", "target_pop_sum"], na.rm=T)
sum(data_national[data_national$vaccinetype == "tOPV", "target_pop_sum"], na.rm=T)

data_quarter <- data %>% group_by(adm0_name, quarter) %>%
  summarize(target_pop_sum = sum(target_pop), 
            population_total_sum = sum(population_total),
            immunity_weighted = weighted.mean(immunity_u5, target_pop, na.rm=T),
            immunity_weighted_summary = weighted.mean(immunity_weighted, population_total_sum, na.rm=TRUE),
            vaccinetype = unique(vaccinetype),
            Region = unique(Region),
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

# summary plot of pre-campaign immunity by vaccine type
data_province <- data_province %>% group_by(source) %>% mutate(pop_weight = target_pop_sum / sum(target_pop_sum, na.rm=T))

ggplot(data_province %>% 
         # filter(Region %in% c("AFRO", "EMRO", "EURO")), 
         filter(Region %in% c("AFRO")), 
         aes(x = source, color = source, y = immunity_weighted)) +
  geom_violin(aes(weights = pop_weight), size = 0.75) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  theme(legend.position = "none") +
  # scale_x_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")+
  ylab("Estimated Pre-Campaign\nType-2 Immunity")

data_province %>% 
  group_by(source, parentactivitycode) %>% 
  # filter(region %in% c("AFRO", "EMRO", "EURO")) %>%
  filter(region %in% c("AFRO")) %>%
  filter(start_date > "2016-05-01") %>%
  summarize(sum = sum(target_pop, na.rm=T)) %>%
  ungroup() %>%
  group_by(source) %>%
  summarize(count = n(),
            min = min(sum),
            q1 = quantile(sum, 0.25),
            median = median(sum),
            mean = mean(sum),
            q3 = quantile(sum, 0.75),
            max = max(sum))

temp <- 
  polis_pops %>%
  filter(region %in% c("AFRO")) %>%
  filter(start_date > "2016-05-01") %>%
  group_by(source, parentactivitycode) %>%
  summarize(sum = sum(target_pop, na.rm=T)) 
wilcox.test(temp$sum ~ temp$source,
            exact = FALSE)


ggplot(data_province %>% filter(!(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO"))), 
       aes(x = source, color = source, y = immunity_weighted)) +
  geom_violin(aes(weights = pop_weight), size = 0.75) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  theme(legend.position = "none") +
  # scale_x_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")+
  scale_y_continuous(limits = c(0,1)) +
  ylab("Estimated Pre-Campaign\nType-2 Immunity")

# summary statistics on pre-campaign immunity by vaccine type
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
data_province %>% filter(Region %in% c("AFRO", "EMRO")) %>%
  ggplot(aes(x = period,y = target_pop_sum, color = source)) + 
  geom_point() +
  geom_smooth() +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  facet_grid(Region~.) +
  scale_y_continuous(limits = c(0, 7e6), name = "Provincial Target Population")

data_province %>%
  filter(is.na(target_pop_sum)==F) %>%
  group_by(source) %>%
  summarize(doses = sum(target_pop_sum),
            min = min(target_pop_sum),
            q1 = quantile(target_pop_sum, 0.25),
            median = median(target_pop_sum),
            mean = mean(target_pop_sum),
            q3 = quantile(target_pop_sum, 0.75),
            max = max(target_pop_sum))

#### Compare locations of nOPV2 and Sabin2 use ####
data_national %>%
  filter(target_pop_sum > 0) %>%
  group_by(adm0_name, source) %>%
  summarize(target_pop_sum = sum(target_pop_sum, na.rm=T)) %>%
  ggplot(aes(fill = source, x = sort(adm0_name), y = target_pop_sum)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  coord_flip() +
  xlab(element_blank()) +
  ylab("Target U5 Population") +
  scale_fill_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")

#### Function to calculate expected count of emergence from campaign ####
#p is the total population (not U-5)
#q is the immunity in children under 5-years
#alpha is the scaling factor for total number of observed emergences (see Supplemental section 4 in Gray 2022) 
Func_u = function(p, q, alpha = 2.025*10^-6){ # alpha in march 2023 was 2.33*10^-6, 2.025*10^-6 in June 2023
  alpha * exp(-0.205*log(p) - 1.98*q) * p
}
alpha = 2.16*10^-6 # Updated 8/17

# testing function
Func_u(1e6, c(.05, .5, .95), alpha = 6.1*10^-6)
Func_u(3e7, c(.05, .5, .95))

# Function to calculate total expected number of emergences using initial value for alpha, just to set things up.
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

##### note, immunity_weighted_summary is used above. Consider moving up
sias_figure <- left_join(sias_figure, data_province_quarter %>% select("quarter", "source", "U_mOPV2_sum"), by = c("quarter", "source"))

# Assuming poisson distribution of #emergences, calculate total number of expected emergences from probability
data_province %>% group_by(vaccinetype) %>% summarize(pop = sum(target_pop_sum, na.rm=T), U_mOPV2 = sum(U_mOPV2, na.rm=T))
data_province %>% group_by(source) %>% summarize(pop = sum(target_pop_sum, na.rm=T), U_mOPV2 = sum(U_mOPV2, na.rm=T))

# Plot number of doses by vaccine type and number of emergences observed and expected
ggplot(sias_figure, aes(x = quarter, fill = source, color = source)) + 
  geom_line(aes(y = target_pop_sum), size = 1) +
  # scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type") +
  theme_bw() +
  geom_col(aes(y = U_mOPV2_sum*1e6), position = position_dodge2(preserve = "single")) +
  geom_point(aes(y = emergences*1e6), color = "black")

# Observed VDPV2 emergences
ggplot(sias_figure, aes(x = quarter)) + 
  geom_col(aes(y = emergences, fill = source)) +
  ylab("VDPV2 Emergences") +
  theme_bw()

#### Surveillance/sequencing lag ####

# Import sequence lab fit (derived from IDM RA Model)
seq_lag_fit <- read_csv("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/results/sequence_lag_fit.csv")
max_lag = 365

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
  geom_line(aes(x=day,y=pseq),colour='black',linetype=2,data=seq_lag_fit) +
  geom_point(aes(x=day,y=pseq_mean),colour='red',data=seq_lag_fit_period) +
  facet_wrap(vars(admin0name))+
  theme(legend.position = 'none') +
  scale_y_continuous('Proportion of Sequences Available',breaks = c(0,0.5,1))+
  scale_x_continuous('Days after collection', limits=c(0,max_lag),oob=scales::squish)

#### Estimate (ttl) lag from index isolate to confirmatory isolate ####

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
  group_by(region_who_code) %>%
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

#visualize parameters
curve(dlnorm(x, meanlog = m, sdlog = s), from=0, to=36)
curve(dlnorm(x, meanlog = m_es, sdlog = s_es), from=0, to=36)
curve(dlnorm(x, meanlog = m*1.25, sdlog = s), from=0, to=36)

sum(dlnorm(1:48, m, s))
cumsum(dlnorm(1:48, m, s)) # 18 month cumsum is 95.6%
cumsum(dlnorm(1:48, m_es, s_es)) # 13 month cumsum is 97.2%

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

#### Fit alpha value to current situation ####

# Count post-switch Sabin-2 emergences
viruses %>% filter(source == "Sabin2", index_isolate == "TRUE", seeding_date > "2016-03-01") %>% 
  filter(region_who_code %in% c("AFRO", "EMRO", "EURO")) %>% 
  # group_by(region_who_code) %>% 
  # group_by(year(virus_date)) %>%
  summarize(count = n())
sabin_emerge <- 74 # 74 post-switch Sabin-2 emergences in AFRO/EMRO/EURO.

# Tune alpha to match the number of observed Sabin-2 emergences observed.
# alpha = 2.15*10^-6 # Updated 8/17
alpha = 2.117*10^-6 # Updated 9/27

data_province <- data_province %>% 
  mutate(U_mOPV2 = Func_u(p=population_total_sum,
                          q = immunity_weighted,
                          alpha = alpha))

# Run new alpha values, update data_province, pass through tt_conv_output chunk below to converge on sabin_emerge

#### tt_conv output: Estimate U_d_i for each day d and campaign i ####

# data_province %>%
#   filter(is.na(U_mOPV2) == F) %>%
#   # filter(!(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO"))) %>%
#   select(adm0_name, adm1_name, Region, period, vaccinetype, source, U_mOPV2, ES) %>%
#   left_join(matrix(nrow=))
# 
# 

tt_conv_data <- data_province %>% 
  filter(is.na(U_mOPV2) == F) %>%
  # filter(!(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO"))) %>%
  # filter(period < 2021.750) %>%
  select(adm0_name, adm1_name, Region, period, vaccinetype, source, U_mOPV2, ES)

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

# tt_conv_data %>% filter(adm0_name %in% c("INDONESIA", "MALAYSIA", "PHILLIPINES")) %>% View() # Check for NA's

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

# Plot points and lines for U_d by period
ggplot() +
  geom_vline(xintercept = 2023.750, size = 1) +
  geom_vline(xintercept = 2021.25, color = "red", alpha = 0.25, size = 1) +
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

# Calculate cdf before defined dates
daily_U_conv %>% 
  # filter(period <= 2023.583) %>%
  # filter(period <= 2023.750) %>%
  # filter(period > 2023.750) %>%
  # filter(period >= 2020, period < 2021) %>%
  # filter(period >= 2021, period < 2022) %>%
  # filter(period >= 2022, period < 2023) %>%
  # filter(period >= 2023, period < 2024) %>%
  # filter(period >= 2023, period <= 2023.50) %>%
  group_by(source) %>%
  summarize(sum(U_mOPV2))

# Convert to quarterly figure
daily_U_conv$year <- floor(daily_U_conv$period)
daily_U_conv$quarter <- Func_period_to_quarter(daily_U_conv$period)
daily_U_conv$year.quarter <- daily_U_conv$year + daily_U_conv$quarter
viruses_count_period$year <- floor(viruses_count_period$period)
viruses_count_period$quarter <- Func_period_to_quarter(viruses_count_period$period)
viruses_count_period$year.quarter <- viruses_count_period$year + viruses_count_period$quarter
viruses_count_quarter <- viruses_count_period %>%
  group_by(source, year.quarter) %>%
  summarize(emergences = sum(emergences))

# Plot points and lines for U_d by quarter
ggplot() +
  geom_vline(xintercept = 2023.750, size = 1) +
  geom_vline(xintercept = 2021.25, color = "red", alpha = 0.25, size = 1) +
  geom_line(data = daily_U_conv, aes(x = year.quarter, y = U_mOPV2, color = source), size = 1) +
  # geom_point(data = daily_U_conv, aes(x = period, y = U_mOPV2, color = source), size = 1) +
  geom_point(data = viruses_count_quarter %>% filter(emergences > 0), aes(x = year.quarter, y = emergences, color = source, shape = source), size = 2) +
  theme_bw() +
  ylab("Expected cVDPV2 Emergences\n(Assuming Seeding at mOPV2 Rate)") +
  xlab("Quarter") +
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

# Plot cumsum points and lines for U_d by period
temp <- left_join(daily_U_conv, viruses_count_period, by = c("period", "source"))
temp[is.na(temp$emergences), "emergences"] <- 0
temp <- temp %>% ungroup() %>% group_by(source) %>% arrange(period) %>% mutate(emergences_cumsum = cumsum(emergences))
temp[temp$period > 2023.75, c("emergences", "emergences_cumsum")] <- NA
temp <- temp %>% 
  select(c("period", "source", "U_mOPV2_cumsum", "emergences_cumsum")) %>%
  pivot_longer(cols = !c("period", "source"))
temp <- temp %>% filter(!(source %in%  c("Sabin2") & name %in% c("U_mOPV2_cumsum")))

ggplot() +
  geom_vline(xintercept = 2023.750, alpha = 0.25, size = 1) +
  geom_vline(xintercept = 2021.25, color = "red", alpha = 0.25, size = 1) +
  geom_line(data = temp, aes(x = period, y = value, color = source, linetype = name), size = 1) +
  theme_bw() +
  ylab("Cumulative cVDPV2 Emergences") +
  scale_x_continuous(limits = c(2016, 2027), breaks = seq(2016, 2027, 2), name = "") +
  scale_linetype_discrete(name = "Emergences",
                          labels = c("Observed", "Expected based\non mOPV2 Rate")) +
  # theme(legend.position = c(0.85, 0.3)) +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(5, "in")) +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")
ggsave("figures/Cumulative Emergences.png", device = "png", units = "in", width = 7, height = 5)

#### Repeat, but only for DRC ####
polis_pops %>% filter(adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>% group_by(source)%>% filter(start_date < today()) %>% summarize(sum = sum(target_pop))

viruses_count_period_DRC <- viruses %>% 
  filter(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  filter(seeding_date > "2016-03-01") %>%
  group_by(period, source) %>% 
  summarize(emergences = sum(index_isolate==TRUE))

viruses %>% 
  filter(admin0name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  filter(source == "Sabin2", index_isolate == "TRUE", seeding_date > "2016-03-01") %>% 
  filter(region_who_code %in% c("AFRO", "EMRO", "EURO")) %>% View() # 17 Sabin2-emergences
alpha_DRC <- 6.6e-06 # based on 17 DRC emergences from Sabin2 use

data_province_DRC <- data_province %>% 
  filter(adm0_name == "DEMOCRATIC REPUBLIC OF THE CONGO") %>%
  mutate(U_mOPV2 = Func_u(p=population_total_sum,
                          q = immunity_weighted,
                          alpha = alpha_DRC))

tt_conv_data <- data_province_DRC %>% 
  filter(is.na(U_mOPV2) == F) %>%
  # filter(!(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO"))) %>%
  # filter(period < 2021.750) %>%
  select(adm0_name, adm1_name, Region, period, vaccinetype, source, U_mOPV2, ES)

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

# Plot points and lines for U_d
ggplot() +
  geom_vline(xintercept = 2023.750, size = 1) +
  # geom_vline(xintercept = 2021.667, color = "red", alpha = 0.25, size = 1) +
  geom_line(data = daily_U_conv, aes(x = period, y = U_mOPV2, color = source), size = 1) +
  geom_point(data = viruses_count_period_DRC %>% filter(emergences > 0), aes(x = period, y = emergences, color = source, shape = source), size = 2) +
  theme_bw() +
  ylab("Expected cVDPV2 Emergences\n(Assuming Seeding at mOPV2 Rate)") +
  xlab("Month") +
  scale_y_continuous(limits = c(0, 2.2)) +
  scale_shape_manual(values = c(19,1)) +
  theme(legend.position = c(0.85, 0.85)) +
  ggtitle("DRC Only") +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")

# Calculate cdf before defined dates
daily_U_conv %>% 
  filter(period <= 2023.750) %>%
  # filter(period > 2023.750) %>%
  group_by(source) %>%
  summarize(sum(U_mOPV2))

#### Repeat, but only for Nigeria ####
polis_pops %>% filter(adm0_name == "NIGERIA") %>% group_by(source)%>% filter(start_date < today()) %>% summarize(sum = sum(target_pop))
# Big question about what population to use for DRC for nOPV2. 30m or 23m?

viruses_count_period_DRC <- viruses %>% 
  filter(admin0name == "NIGERIA") %>%
  filter(seeding_date > "2016-03-01") %>%
  group_by(period, source) %>% 
  summarize(emergences = sum(index_isolate==TRUE))

viruses %>% 
  filter(admin0name == "NIGERIA") %>%
  filter(source == "Sabin2", index_isolate == "TRUE", seeding_date > "2016-03-01") %>% 
  filter(region_who_code %in% c("AFRO", "EMRO", "EURO")) %>% View() # 10 Sabin2-emergences
alpha_Nigeria <- 1.58e-06 # based on 10 NIGERIA emergences from Sabin2 use

data_province_NIE <- data_province %>% 
  filter(adm0_name == "NIGERIA") %>%
  mutate(U_mOPV2 = Func_u(p=population_total_sum,
                          q = immunity_weighted,
                          alpha = alpha_Nigeria))

tt_conv_data <- data_province_NIE %>% 
  filter(is.na(U_mOPV2) == F) %>%
  # filter(!(adm0_name %in% c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO"))) %>%
  # filter(period < 2021.750) %>%
  select(adm0_name, adm1_name, Region, period, vaccinetype, source, U_mOPV2, ES)

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

# Plot points and lines for U_d
ggplot() +
  geom_vline(xintercept = 2023.750, size = 1) +
  # geom_vline(xintercept = 2021.25, color = "red", alpha = 0.25, size = 1) +
  geom_line(data = daily_U_conv, aes(x = period, y = U_mOPV2, color = source), size = 1) +
  geom_point(data = viruses_count_period_DRC %>% filter(emergences > 0), aes(x = period, y = emergences, color = source, shape = source), size = 2) +
  theme_bw() +
  ylab("Expected cVDPV2 Emergences\n(Assuming Seeding at mOPV2 Rate)") +
  xlab("Month") +
  scale_y_continuous(limits = c(0, 3)) +
  scale_shape_manual(values = c(19,1)) +
  theme(legend.position = c(0.85, 0.85)) +
  ggtitle("NIGERIA Only") +
  scale_color_discrete(labels = c("nOPV2", "Sabin2"), name = "Vaccine Type")

# Calculate cdf before defined dates (eg, March 1)
daily_U_conv %>% 
  filter(period <= 2023.750) %>%
  # filter(period > 2023.750) %>%
  group_by(source) %>%
  summarize(sum(U_mOPV2))

#### Estimate uncertainty ####
expand_grid(count = c(43, 66)) %>% 
  +   mutate(lower = qpois(0.025,count), upper = qpois(0.975,count))

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

#### Spread of nOPV2 Emergences ####
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

viruses_supplemented %>% filter(active == TRUE) %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(max = max(AFP_cumsum),
            source = unique(source)) %>%
  filter(max < 10)

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
  ggplot(aes(x = AFP_cumsum>=10)) +
  geom_bar(aes(fill = source), stat="count", width = 1, position = "dodge") +
  # geom_density() +
  xlab("10 or more AFP cases in Emergence Group") +
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

  