# Figures for SC Meeting 

# begin with update_cases.R and aggregate_cases.R

# Epi curve, all cases:
detections_months = detections %>% mutate(period = year(date) + (month(date)-1)/12) %>%
  group_by(source,final_class,period, vdpv_emergence_group) %>% tally()

detections_months$vdpv_emergence_group_simple <- ifelse(detections_months$vdpv_emergence_group == "NIE-JIS-1", "NIE-JIS-1", 
                                                        ifelse(detections_months$vdpv_emergence_group == "PAK-GB-1", "PAK-GB-1",
                                                               ifelse(detections_months$vdpv_emergence_group == "NIE-ZAS-1", "NIE-ZAS-1",
                                                                             ifelse(detections_months$vdpv_emergence_group == "CHA-NDJ-1", "CHA-NDJ-1",
                                                                                    ifelse(detections_months$vdpv_emergence_group == "YEM-TAI-1", "YEM-TAI-1",
                                                                                           ifelse(detections_months$vdpv_emergence_group == "RDC-MAN-3", "RDC-MAN-3",
                                                                      "Other Emergence Groups"))))))
detections_months$vdpv_emergence_group_simple <- factor(detections_months$vdpv_emergence_group_simple, 
                                                        levels = c("NIE-JIS-1", "CHA-NDJ-1", "PAK-GB-1", "Other Emergence Groups", "NIE-ZAS-1", "RDC-MAN-3", "YEM-TAI-1"))
#Detections by month, country, and emergence group
detections_months_2 = detections %>% mutate(period = year(date) + (month(date)-1)/12) %>%
  group_by(source,final_class,period, vdpv_emergence_group, iso3_code) %>% tally()

detections_months_2 <- detections_months_2 %>% mutate(vdpv_emergence_group_simple2 = "Other Emergence Groups")
detections_months_2[detections_months_2$iso3_code == "NGA" & detections_months_2$vdpv_emergence_group == "NIE-ZAS-1", "vdpv_emergence_group_simple2"] <- "NIE-ZAS-1_Nigeria"
detections_months_2[detections_months_2$iso3_code %in% c("TGO", "NER", "GIN", "GHA", "CIV") & detections_months_2$vdpv_emergence_group == "NIE-ZAS-1", "vdpv_emergence_group_simple2"] <-
  "NIE-ZAS-1_West"
detections_months_2[detections_months_2$iso3_code %in% c("CMR", "TCD", "CAF") & detections_months_2$vdpv_emergence_group == "NIE-ZAS-1", "vdpv_emergence_group_simple2"] <- "NIE-ZAS-1_East"
detections_months_2[detections_months_2$iso3_code == "MOZ", "vdpv_emergence_group_simple2"] <- "Mozambique"
detections_months_2[detections_months_2$iso3_code == "COD" & detections_months_2$vdpv_emergence_group %in% c("RDC-MAN-2","RDC-MAN-3","RDC-MAN-4", "RDC-MAN-5", "RDC-MAN-3, RDC-MAN-4"), "vdpv_emergence_group_simple2"] <-
  "RDC-MAN"
detections_months_2[detections_months_2$iso3_code == "YEM", "vdpv_emergence_group_simple2"] <- "Yemen"

detections_months_2$vdpv_emergence_group_simple2 <- factor(detections_months_2$vdpv_emergence_group_simple2, 
                                                        levels = c("NIE-ZAS-1_Nigeria", "NIE-ZAS-1_West", "NIE-ZAS-1_East", "RDC-MAN", "Mozambique", "Yemen", "Other Emergence Groups"))
#Detections by month, country, and emergence group
detections_months_3 = detections %>% mutate(period = year(date) + (month(date)-1)/12) %>%
  group_by(source,final_class,period, vdpv_emergence_group, iso3_code) %>% tally()

detections_months_3 <- detections_months_3 %>% mutate(vdpv_emergence_group_simple3 = "Other")
detections_months_3[detections_months_3$iso3_code == "NGA", "vdpv_emergence_group_simple3"] <- "Nigeria"
detections_months_3[detections_months_3$iso3_code %in% c("AFG", "PAK"), "vdpv_emergence_group_simple3"] <-
  "Afg/Pak"
detections_months_3[detections_months_3$iso3_code == "COD", "vdpv_emergence_group_simple3"] <-
  "DRC"
detections_months_3[detections_months_3$iso3_code == "YEM", "vdpv_emergence_group_simple3"] <- "Yemen"

detections_months_3$vdpv_emergence_group_simple3 <- factor(detections_months_3$vdpv_emergence_group_simple3, 
                                                           levels = c("Nigeria", "Afg/Pak", "DRC","Yemen", "Other"))

#Affected countries by month
countries_affected = detections %>% filter(final_class == "cVDPV2") %>% mutate(period = year(date) + (month(date)-1)/12) %>% 
  count(iso3_code, period) %>% count(period)

# number of affected countries over time
plot(countries_affected)
date_from = "2016-01-01"
ggplot() + 
  geom_col(aes(x=period, y=n), data = countries_affected %>%
             filter(period>= period_from), orientation = 'x') +
  geom_hline(yintercept=0) +
  theme_bw() +
  xlab(NULL) 

detections_months_3 <- detections_months_3 %>% mutate(vdpv_emergence_group_simple3 = "Other")
detections_months_3[detections_months_3$iso3_code == "NGA", "vdpv_emergence_group_simple3"] <- "Nigeria"
detections_months_3[detections_months_3$iso3_code %in% c("AFG", "PAK"), "vdpv_emergence_group_simple3"] <-
  "Afg/Pak"
detections_months_3[detections_months_3$iso3_code == "COD", "vdpv_emergence_group_simple3"] <-
  "DRC"
detections_months_3[detections_months_3$iso3_code == "YEM", "vdpv_emergence_group_simple3"] <- "Yemen"

detections_months_3$vdpv_emergence_group_simple3 <- factor(detections_months_3$vdpv_emergence_group_simple3, 
                                                           levels = c("Nigeria", "Afg/Pak", "DRC","Yemen", "Other"))

#epi curve, recent cases:
date_from = "2016-01-01"
period_from = year(date_from) + (month(date_from)-1)/12
ggplot() + 
  geom_col(aes(x=period,y=n),data = detections_months %>% 
             filter(period >= period_from,source=='AFP',
                    final_class == "cVDPV2")
           ,orientation='x') +
  # geom_col(aes(x=period,y=-n,fill = final_class),data = detections_months %>% filter(period >= period_from,source=='ES'),orientation='x') + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlab(NULL) + 
  scale_y_continuous(NULL,labels = function(x){if_else(x==0,'0',if_else(x>0,str_c(x, ' AFP'),str_c(abs(x), ' ES')))}) + 
  scale_fill_discrete(NULL)+
  theme(legend.position = "none")
ggsave('results/epicurve_recent_cVDPV2.png',width=4,height=1.5,units='in',dpi=600)

#epi curve, recent cases, by emergence group:
date_from = "2020-01-01"
period_from = year(date_from) + (month(date_from)-1)/12
ggplot() + 
  geom_col(aes(x=period,y=n,fill = vdpv_emergence_group_simple3),data = detections_months_3 %>% 
             filter(period >= period_from,source=='AFP',
                    final_class == "cVDPV2")
           ,orientation='x') + 
  # geom_text(aes(x=period, y=n, label = n), data = detections_months %>% 
  #             filter(period >= period_from,source=='AFP',
  #                    final_class == "cVDPV2")) +
  # geom_col(aes(x=period,y=-n,fill = final_class),data = detections_months %>% filter(period >= period_from,source=='ES'),orientation='x') + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlab(NULL) + 
  scale_y_continuous(NULL,labels = function(x){if_else(x==0,'0',if_else(x>0,str_c(x, ' AFP'),str_c(abs(x), ' ES')))}) + 
  scale_fill_brewer(type = 'qual', palette = 6) +
  # facet_wrap(vars(vdpv_emergence_group_simple3),ncol=1) +
  theme(legend.position = "none") 
ggsave('results/epicurve_recent_byemergencegroup.png',width=4,height=8,units='in',dpi=600)

# Table of recent cases
temp <- detections_months %>% filter(final_class == "cVDPV2", period > 2019, source == "AFP") %>%
  group_by(vdpv_emergence_group) %>% mutate(sum(n))
View(temp)

temp <- detections_months %>% filter(final_class == "cVDPV2", period > 2019, source == "AFP") %>%
  group_by(vdpv_emergence_group_simple) %>% mutate(sum(n))
View(temp)


# Focus on JIS-1
date_from = "2018-05-01"   ##Double-check this
period_from = year(date_from) + (month(date_from)-1)/12
detections_months$JIS1_response <- ifelse(detections_months$period <= period_from, "Pre-Response", "Post-Response")

date_from = "2017-01-01"
period_from = year(date_from) + (month(date_from)-1)/12
ggplot() + 
  geom_col(aes(x=period,y=n, fill = JIS1_response),data = detections_months %>% 
             filter(period >= period_from,source=='AFP',
                    final_class == "cVDPV2",
                    vdpv_emergence_group == "NIE-JIS-1")
           ,orientation='x') + 
  # geom_col(aes(x=period,y=-n,fill = final_class),data = detections_months %>% filter(period >= period_from,source=='ES'),orientation='x') + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlab(NULL) + 
  scale_y_continuous(limits = c(0,85),labels = function(x){if_else(x==0,'0',if_else(x>0,str_c(x, ' AFP'),str_c(abs(x), ' ES')))}) + 
  scale_fill_discrete(NULL)+
  theme(legend.position = "none")
ggsave('results/epicurve_recent_JIS1.png',width=4,height=1.5,units='in',dpi=600)


# Focus on ZAS-1
date_from = "2021-06-28"   ##Double-check this
period_from = year(date_from) + (month(date_from)-1)/12
detections_months$ZAS1_response <- ifelse(detections_months$period <= period_from, "Pre-Response", "Post-Response")

date_from = "2019-01-01"
period_from = year(date_from) + (month(date_from)-1)/12
ggplot() + 
  geom_col(aes(x=period,y=n, fill = ZAS1_response),data = detections_months %>% 
             filter(period >= period_from,source=='AFP',
                    final_class == "cVDPV2",
                    vdpv_emergence_group == "NIE-ZAS-1")
           ,orientation='x') + 
  # geom_col(aes(x=period,y=-n,fill = final_class),data = detections_months %>% filter(period >= period_from,source=='ES'),orientation='x') + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlab(NULL) + 
  scale_y_continuous(limits = c(0,85),labels = function(x){if_else(x==0,'0',if_else(x>0,str_c(x, ' AFP'),str_c(abs(x), ' ES')))}) + 
  scale_fill_discrete(NULL)+
  theme(legend.position = "none")
ggsave('results/epicurve_recent_JIS1.png',width=4,height=1.5,units='in',dpi=600) 