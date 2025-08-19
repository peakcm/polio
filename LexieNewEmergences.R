# Updates requested by Lexie
# LOAD VIRUSES first (from VDPV2n_analyses.R)


viruses %>% filter(virus_date > "2022-01-01") %>%
  arrange(vdpv_emergence_group_name,-index_isolate,seeding_date) %>%
  select(vdpv_emergence_group_name, 
         quarter, index_isolate, # get quarter of index case#virus_date, index_date, 
         #seeding_date, 
         source, 
         surveillance_type_name, 
         admin0name) %>%
  mutate(
    year = floor(quarter),  # Extract the year
    quarter = case_when(
      round(quarter, 2) %% 1 == 0 ~ "Q1",   # Detects .0 and .00 as Q1
      round(quarter, 2) %% 1 == 0.25 ~ "Q2",
      round(quarter, 2) %% 1 == 0.50 ~ "Q3",
      round(quarter, 2) %% 1 == 0.75 ~ "Q4"
    ),
    period_emerge = paste(quarter, year)  # Combine Q# and year
  )  %>% select(-year,-quarter) -> dp


# index info - emergence post 2022 Jan
dp %>% filter(index_isolate == "TRUE") %>% 
  select(-surveillance_type_name, -index_isolate) %>%
  rename_with(~ paste0(., "_index")) -> index_info

# number of countries spread to
dp %>% filter(surveillance_type_name %in% c("AFP","Environmental")) %>%
  group_by(vdpv_emergence_group_name) %>%
  summarise(countries_spread_count_total = n_distinct(admin0name), 
            countries_spread_list_total = paste(unique(admin0name), collapse = ", "), 
            .groups = "drop")  -> spread_info_total

# number of AFP vs. Environmental cases
dp %>% filter(surveillance_type_name %in% c("AFP","Environmental")) %>%
  group_by(vdpv_emergence_group_name, surveillance_type_name) %>% 
  summarise(case_counts = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = surveillance_type_name,
    values_from = case_counts,
    names_glue = "{surveillance_type_name}_{.value}"
  ) -> case_info

# number of countries spread to by AFP vs Environmental
dp %>% filter(surveillance_type_name %in% c("AFP","Environmental")) %>%
  group_by(vdpv_emergence_group_name, surveillance_type_name) %>%
  summarise(countries_spread_count = n_distinct(admin0name), 
            countries_spread_list = paste(unique(admin0name), collapse = ", "), 
            .groups = "drop") %>%
  filter(surveillance_type_name %in% c("AFP","Environmental")) %>% 
  pivot_wider(
    names_from = surveillance_type_name,  
    values_from = c(countries_spread_count, countries_spread_list), 
    names_glue = "{surveillance_type_name}_{.value}"
  ) -> spread_info

index_info %>% left_join(spread_info_total,
                         by = c("vdpv_emergence_group_name_index" = "vdpv_emergence_group_name")) %>% 
  left_join(spread_info,
            by = c("vdpv_emergence_group_name_index" = "vdpv_emergence_group_name")) %>%
  left_join(case_info,
            by = c("vdpv_emergence_group_name_index" = "vdpv_emergence_group_name")) %>%
  select("vdpv_emergence_group_name_index",
         "source_index",
         "admin0name_index", "period_emerge_index",
         "countries_spread_count_total", "countries_spread_list_total", 
         "AFP_case_counts", "Environmental_case_counts",
         "AFP_countries_spread_count",  "AFP_countries_spread_list",
         "Environmental_countries_spread_count","Environmental_countries_spread_list") -> output_share 

dp_afp_wide <- dp %>% 
  filter(vdpv_emergence_group_name %in% output_share$vdpv_emergence_group_name_index) %>%
  filter(surveillance_type_name == "AFP") %>%
  group_by(vdpv_emergence_group_name, source, admin0name) %>%
  summarise(case_counts = n(), .groups = "drop")# %>%
#pivot_wider(names_from = admin0name, values_from = case_counts, values_fill = 0)

dp_env_wide <- dp %>% 
  filter(vdpv_emergence_group_name %in% output_share$vdpv_emergence_group_name_index) %>%
  filter(surveillance_type_name == "Environmental") %>%
  group_by(vdpv_emergence_group_name, source, admin0name) %>%
  summarise(case_counts = n(), .groups = "drop") #%>%
#pivot_wider(names_from = admin0name, values_from = case_counts, values_fill = 0)

write.csv(output_share,"cvdpv2spread.csv") # this is the file I share with Lexie
write.csv(dp_afp_wide,"dp_afp_wide.csv")
write.csv(dp_env_wide,"dp_env_wide.csv")
