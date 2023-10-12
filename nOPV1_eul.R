# Analyses to support nOPV1 EUL
rm(list=ls())

#### Load Libraries ####
# library(ggridges)
library(tidyverse)
library(lubridate)
library(PolisAPI)
# library(ggrepel)
# library(GGally)

#### Import virus data ####
token = readLines("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/data_local/token.txt")[1]
viruses_raw_1 = get_polis_virus(token = token, min_date = "2003-01-01", virus_id = 1)
viruses_raw_2 = get_polis_virus(token = token, min_date = "2003-01-01", virus_id = 4)
viruses_raw <- rbind(viruses_raw_1, viruses_raw_2)

# Clean data
viruses = viruses_raw %>%
  select(id, epid, virus_date, surveillance_type_name, region_who_code,
         admin0name, country_iso3code, admin1name, admin1guid, admin2name, admin2guid,
         virus_type_name, vdpv_classification_name, vdpv_emergence_group_name, 
         vdpv_nt_changes_from_sabin, surveillance_type_name, po_ns_seq_date) %>%
  mutate(virus_date = ymd(as.Date(virus_date)),
         id = as.numeric(id),
         dot_name = paste(admin0name, admin1name, sep = ":"),
         dot_year_month = tolower(paste(year(virus_date), month(virus_date), sep = ":"))) %>% 
  filter(virus_type_name %in% c("cVDPV1", "cVDPV2"),
         !is.na(vdpv_emergence_group_name)) %>%
  group_by(vdpv_emergence_group_name) %>%
  mutate(index_date = min(virus_date, na.rm=T),
         most_recent_date = max(virus_date, na.rm=T),
         index_isolate = (virus_date == index_date),
         most_recent_isolate = (virus_date == most_recent_date),
         active = most_recent_date >= today()-months(6),
         time_since_index = virus_date - index_date)

# Check emergence group names
novel_emergences <- c("RDC-SKV-1", "RDC-TAN-2", "RDC-KOR-1", "CAF-KEM-1", "NIE-KBS-1", "RDC-HKA-2","CAF-BNG-3", "BOT-FRA-1", "EGY-NOR-1")
novel_emergences %in% c(viruses$vdpv_emergence_group_name %>% unique())

# If there are ties for the most recent isolate, select one
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
viruses$source = "Sabin"
sort(unique(viruses$vdpv_emergence_group_name))
viruses[viruses$vdpv_emergence_group_name %in% novel_emergences, "source"] <- "nOPV2"

# Calculate period and quarter
viruses$week <- year(viruses$virus_date) + (week(viruses$virus_date)-1)/52
viruses$period <- year(viruses$virus_date) + (month(viruses$virus_date)-1)/12
viruses$quarter <- year(viruses$virus_date) + (quarter(viruses$virus_date)-1)/4

# Dataset for labeling ends
viruses_ends <- viruses %>%
  filter(time_since_index <= 365) %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(desc(virus_date)) %>%
  slice(1) %>%
  ungroup()

# Calculate cumsum of AFP cases
viruses$AFP_count <- 0
viruses[viruses$surveillance_type_name == "AFP", "AFP_count"] <- 1
viruses <- viruses %>%
  group_by(vdpv_emergence_group_name) %>%
  arrange(virus_date) %>%
  mutate(AFP_cumsum = cumsum(AFP_count)) %>%
  ungroup()

# Label date 10th case
viruses$date_case_10 <- NA
for (i in viruses$vdpv_emergence_group_name){
  if (!(is.na(i))){
    if (max(viruses[viruses$vdpv_emergence_group_name %in% i, "AFP_cumsum"]) >= 10){
      date_case_10 <- viruses[viruses$vdpv_emergence_group_name %in% i &
                                       viruses$surveillance_type_name %in% c("AFP") &
                                       viruses$AFP_cumsum == 10, "virus_date"]
      viruses[viruses$vdpv_emergence_group_name %in% i, "date_case_10"] <- date_case_10
    }
  }
}

# Identify 'sputtering' groups
viruses %>% filter(active == TRUE) %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(max = max(AFP_cumsum),
            source = unique(source)) %>%
  filter(max < 10)

# Create virus_type_group
group_cvdpv2_pre <- viruses[viruses$virus_type_name %in% c("cVDPV2") &
                              viruses$index_isolate == TRUE &
                              viruses$seeding_date < "2016-05-01",]$vdpv_emergence_group_name
group_cvdpv2_post <- viruses[viruses$virus_type_name %in% c("cVDPV2") &
                              viruses$index_isolate == TRUE &
                              viruses$seeding_date >= "2016-05-01",]$vdpv_emergence_group_name
viruses$virus_type_group <- NA
viruses[viruses$virus_type_name %in% c("cVDPV1"), "virus_type_group"] <- "cVDPV1"
viruses[viruses$vdpv_emergence_group_name %in% group_cvdpv2_pre, "virus_type_group"] <- "cVDPV2 Pre-Switch"
viruses[viruses$vdpv_emergence_group_name %in% group_cvdpv2_post, "virus_type_group"] <- "cVDPV2 Post-Switch"

#### Figures ####
# Summary of emergences
viruses %>% 
  filter(index_isolate == TRUE) %>% 
  ungroup() %>% 
  group_by(admin0name, virus_type_group, source) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = reorder(admin0name, -count), y = count, fill = source)) +
  facet_grid(virus_type_group~.) +
  geom_col() + coord_flip() + 
  xlab("") + ylab("Number of Emergences since 2003")
ggsave(filename = "figures/number-of-emergences-by-country.png", device = "png", units = "in", width = 7, height = 9)

# Duration of emergence groups
viruses_summary <-
  viruses %>%
  group_by(vdpv_emergence_group_name, virus_type_name) %>%
  summarize(count = n(),
            cases = sum(surveillance_type_name == "AFP"),
            index_date = unique(index_date),
            most_recent_date = unique(most_recent_date),
            date_case_10 = unique(date_case_10),
            duration = most_recent_date - index_date,
            time_to_to = date_case_10 - index_date,
            seeding_date = min(seeding_date, na.rm=T),
            virus_type_group = unique(virus_type_group))

viruses_summary %>%
  ggplot(aes(x = duration/365, fill = virus_type_group)) +
  # facet_grid(virus_type_name~.) +
  geom_histogram(aes(y = ..count..), position = "dodge", bins = 5) +
  ylab("Count") + xlab("Duration of Emergence Group (years)")
ggsave("figures/histogram duration count.png", device = "png", units = "in", width = 5, height = 5)

viruses_summary %>%
  ggplot(aes(x = duration/365, fill = virus_type_group)) +
  # facet_grid(virus_type_name~.) +
  geom_density(aes(y = ..density..), position = "dodge", alpha = 0.5) +
  ylab("Density") + xlab("Duration of Emergence Group (years)")
ggsave("figures/histogram duration density.png", device = "png", units = "in", width = 5, height = 5)

viruses_summary %>%
  ggplot() +
  geom_boxplot(aes(x = virus_type_group, fill = virus_type_group, y = duration/365)) +
  ylab("Duration of Emergence Group (years)") +
  xlab("") +
  theme(legend.position = "none")
ggsave("figures/boxplot duration.png", device = "png", units = "in", width = 5, height = 5)


# Cases of emergence groups
viruses_summary %>%
  ggplot(aes(x = cases, fill = virus_type_group)) +
  # facet_grid(virus_type_name~.) +
  geom_histogram(position = "dodge", bins = 5) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  ylab("Count") + xlab("Cases in Emergence Group")
ggsave("figures/histogram case count.png", device = "png", units = "in", width = 5, height = 5)

viruses_summary %>%
  ggplot(aes(x = cases, fill = virus_type_group)) +
  # facet_grid(virus_type_name~.) +
  geom_density(position = "dodge", alpha = 0.5) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  ylab("Density") + xlab("Cases in Emergence Group")
ggsave("figures/histogram case density.png", device = "png", units = "in", width = 5, height = 5)

viruses_summary %>%
  ggplot(aes(y= cases, x = virus_type_group, fill = virus_type_group)) +
  geom_boxplot() +
  scale_y_log10() +
  xlab("") + ylab("Cases in Emergence Group") +
  theme(legend.position = "none")
ggsave("figures/boxplot case count.png", device = "png", units = "in", width = 5, height = 5)

# Days to 10 cases
viruses_summary %>%
  ggplot(aes(x = (date_case_10 - index_date)/365, fill = virus_type_group)) +
  # facet_grid(virus_type_name~.) +
  geom_histogram(position = "dodge", bins = 5) +
  ylab("Count") + xlab("Years until 10 cases in Emergence Group")
ggsave("figures/histogram days to 10 cases count.png", device = "png", units = "in", width = 5, height = 5)

viruses_summary %>%
  ggplot(aes(x = (date_case_10 - index_date)/365, fill = virus_type_group)) +
  # facet_grid(virus_type_name~.) +
  geom_density(position = "dodge", alpha = 0.5) +
  ylab("Density") + xlab("Years until 10 cases in Emergence Group")
ggsave("figures/histogram days to 10 cases density.png", device = "png", units = "in", width = 5, height = 5)

viruses_summary %>%
  ggplot(aes(y = (date_case_10 - index_date)/365, x = virus_type_group, fill = virus_type_group)) +
  geom_boxplot() +
  xlab("") + ylab("Years until 10 cases in Emergence Group")+
  theme(legend.position = "none")
ggsave("figures/boxplot days to 10 cases count.png", device = "png", units = "in", width = 5, height = 5)

# Growth trajectories
ggplot(data = viruses %>% ungroup(),
       aes(x = time_since_index, y = AFP_cumsum, 
           group = vdpv_emergence_group_name, color = virus_type_group)) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(0, 365), name = "Days Since Index Isolate\n(first year of circulation)") +
  scale_y_log10(name = "Cumulative AFP Cases Reported") +
  scale_color_manual(values = c("red", "#99999950", "black"), name = "Emergence Category") +
  # facet_grid(virus_type_group~.) +
  theme_bw() 
ggsave("figures/growth trajectories.png", device = "png", units = "in", width = 7, height = 7)

temp <- viruses %>%
  filter(vdpv_emergence_group_name %in% viruses_summary[viruses_summary$duration > 100,]$vdpv_emergence_group_name) %>% #lasted at least 100 days
  filter(time_since_index <= 100) %>%
  group_by(vdpv_emergence_group_name) %>%
  summarize(cases_at_100 = max(AFP_cumsum),
            virus_type_name = unique(virus_type_name))

wilcox.test(temp$cases_at_100 ~ temp$virus_type_name,
            exact = FALSE)

ggplot(temp) +
  geom_density(aes(x = cases_at_100, fill = virus_type_name), position = "dodge", alpha = 0.5) +
  scale_x_continuous(name = "AFP Cases in Emergence Group at Day 100") +
  ylab("Density")
ggsave("figures/cases at 100 days.png", device = "png", units = "in", width = 5, height = 5)

ggplot(temp) +
  geom_boxplot(aes(y = cases_at_100, x = virus_type_name, fill = virus_type_name)) +
  xlab("AFP Cases in Emergence Group at Day 100") +
  ylab("") +
  theme(legend.position = "none")
ggsave("figures/boxplot cases at 100 days.png", device = "png", units = "in", width = 5, height = 5)

