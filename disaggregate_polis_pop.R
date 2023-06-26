rm(list=ls())
library(PolisAPI)
library(httr)


get_polis_url2 <- function(my_url, verbose = T, clean = T, cast=TRUE, ntries = 3, count = FALSE, api_token = get_token()){
  all_results = NULL
  initial_query = my_url
  while(!is.null(my_url)){
    
    result = 404L
    try_count = 1
    
    while(http_error(result) & try_count <= ntries){
      result = GET(my_url, add_headers("authorization-token" = api_token))
      try_count = try_count + 1
    }
    stop_for_status(result, task = 'pull data from POLIS')
    result_content = content(result, type='text', encoding = 'UTF-8') 
    assertthat::assert_that(jsonlite::validate(result_content), msg = 'POLIS data not in JSON format')
    result_content = jsonlite::fromJSON(result_content)
    
    
    all_results = bind_rows(all_results,mutate_all(result_content$value,as.character))
    my_url =result_content$odata.nextLink
    if(verbose) cat('.')
  }
  if(!is.null(result_content$odata.count)){
    if(nrow(all_results) != as.numeric(result_content$odata.count)){
      warning(paste0('Expected ',result_content$odata.count, ' results, returned ',nrow(all_results))) 
    }
  }
  if(clean) all_results = all_results %>% PolisAPI:::clean_polis()
  if(cast) all_results = cast_polis(all_results, PolisAPI:::polis_table_name(initial_query))
  attr(all_results,'query') = initial_query
  return(all_results)
}

set_token("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/data_local/token.txt")

sia_url <- "https://extranet.who.int/polis/api/v2/SubActivity?$filter=(DateFrom ge DateTime'2014-01-01') and (VaccineType eq 'mOPV2' or VaccineType eq 'tOPV' or VaccineType eq 'nOPV2' or VaccineType eq 'OPV2*')"
sia_polis <- get_polis_url2(URLencode(sia_url))

sia_tsir <- read_csv('C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub//polio-immunity-mapping/results/sia_district_rows.csv') %>%
  janitor::clean_names()
data_tsir <- read_csv('C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub//polio-immunity-mapping/results/tsir_data.csv')
pop_tsir <- data_tsir %>% filter(period == 2023) %>% distinct(guid, population)
sia_polis_target_pop <- sia_polis %>% 
  group_by(childactivitycode = sia_sub_activity_code) %>% 
  summarise(target_pop_total = sum(unique(calculated_target_population)), .groups = 'keep')

sia <- left_join(sia_tsir, pop_tsir) %>% left_join(sia_polis_target_pop)
sia <- sia %>% group_by(childactivitycode) %>%
  mutate(target_pop_frac = fraction*population/sum(population),
         target_pop = target_pop_frac* target_pop_total)

#no tsir population data for philippines, indonesia, malaysia. Perhaps ignore?
sia %>% ungroup %>% filter(start_date > ymd('2014-01-01'), str_detect(vaccinetype, 'tOPV|OPV2')) %>%
  filter(is.na(target_pop)) %>% group_by(childactivitycode, adm0_name, start_date) %>% tally() %>%
  arrange(desc(start_date))


sia %>% ungroup %>% filter(start_date > ymd('2014-01-01'), str_detect(vaccinetype, 'tOPV|OPV2')) %>%
  write_rds('C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/sia_polis_target_pop.rds', compress = 'gz')

sia %>% ungroup %>% 
  filter(vaccinetype == 'nOPV2', status == 'Done') %>% 
  summarise(total = sum(target_pop, na.rm=T)/1e6)
