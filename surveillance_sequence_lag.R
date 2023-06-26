library(tidyverse)
library(lubridate)
library(purrr)
#library(brms)
#options(mc.cores = 4)
setwd("~/GitHub/polio-immunity-mapping")
source("R/pipeline/load_args.R")
args <- load_args()

source('R/utils/polis_download_api_v2.R')
token = readLines(config::get("filename_token"),n=1L)

serotype_use <- ifelse(is.null(args$serotype_use),
                       config::get("serotype_use"),
                       args$serotype_use)
default_serotype_use <- config::get("default_serotype_use")
filename_fit <- config::get("filename_sequence_lag")

if(serotype_use %in% c(2,4)){
  #pull VDPV2 from POLIS:
  virus = get_polis_virus(token=token,virus_id = c(4),updated_since = ymd('2016-06-01'))
  
  virus= virus %>% clean_names() %>% 
    mutate(across(c("vdpv_reported_to_hq_date", "virus_date"), ~as_date(ymd_hms(.)))) %>%
    # mutate(across(matches('date$'),~as_date(ymd_hms(.)))) %>%
    mutate(seq_lag = as.numeric(vdpv_reported_to_hq_date-virus_date))
} else {
  filename_fit <- glue::glue(tools::file_path_sans_ext(filename_fit), "_type{serotype_use}.csv")
  #Pull WILD1 from POLIS... use reported to HQ date to determine lag (mostly missing!)
  virus_raw = get_polis_afp(token=token,updated_since = ymd('2001-01-01'),filter = " and WILD1")
  
  virus = virus_raw %>% clean_names() %>% 
    mutate(across(matches('date'),~as_date(ymd_hms(.x)))) %>%
    mutate(seq_lag = as.numeric(date_notificationto_hq-paralysis_onset_date))
}
#estimate distribution with empirical bayes----

# get initial value from pooled data:
max_lag = 365
xx = virus %>% filter(!is.na(seq_lag),seq_lag > 0, seq_lag < max_lag) %>% pull(seq_lag)
est = c(log(mean(xx)), log(1))
est = optim(est, function(x){ -sum(dgamma(xx, rate = exp(x[2] - x[1]),shape = exp(x[2]),log=TRUE)) })$par
est

#curve(dgamma(x,rate = exp(est[2]-est[1]), shape = exp(est[2])),from=0,to=150)

# fit separately to each country:
fit_eb = virus %>% filter(!is.na(seq_lag), seq_lag > 0, seq_lag <max_lag)%>% 
  group_by(admin0name) %>% mutate(n = n()) %>% ungroup %>% 
  filter(n >= 2) %>% 
  split(.$admin0name) %>%
  map_df(~{
    par = optim(est, function(x){ -sum(dgamma(.x$seq_lag, rate = exp(x[2] - x[1]),shape = exp(x[2]),log=TRUE)) })$par
    tibble(lmu = par[1], lshape = par[2])
  }, .id = 'admin0name')

# point estimate of hyper parameters:
lmu_mean = mean(fit_eb$lmu)
lmu_sd = sd(fit_eb$lmu)

lshape_mean = mean(fit_eb$lshape)
lshape_sd = sd(fit_eb$lshape)

#iteratively get posterior mode and update hyper parameter point estimates:
for(i in 1:10){
  fit_eb = virus %>% filter(!is.na(seq_lag), seq_lag >0, seq_lag < max_lag) %>%
    group_by(admin0name) %>% mutate(n = n()) %>% ungroup %>% 
    split(.$admin0name) %>%
    map_df(~{
      par = optim(est, function(x){ 
        -sum(dgamma(.x$seq_lag, rate = exp(x[2] - x[1]), shape = exp(x[2]),log=TRUE)) + 
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
  geom_histogram(aes(y=after_stat(density),fill=admin0name,x=seq_lag),data=virus %>% filter(!is.na(seq_lag)))+
  geom_line(aes(x=day,y=dseq),colour='black',linetype=2,data=df_fit_eb) + 
  facet_wrap(vars(admin0name),scales = 'free_y')+
  theme(legend.position = 'none') + 
  scale_x_continuous('Days after collection',limits=c(0,max_lag),oob=scales::squish)
# ggsave(glue::glue('results/seq_lag_fit_type{serotype_use}.png'),width=12,height=10,units='in')

ggplot() + 
  geom_step(aes(colour=admin0name,x=seq_lag),stat='ecdf',data=virus %>% filter(!is.na(seq_lag), seq_lag < 365, seq_lag > 0))+
  geom_line(aes(x=day,y=pseq),colour='black',linetype=2,data=df_fit_eb) + 
  facet_wrap(vars(admin0name))+
  theme(legend.position = 'none') + 
  scale_y_continuous('Proportion of Sequences Available',breaks = c(0,0.5,1))+
  scale_x_continuous('Days after collection', limits=c(0,max_lag),oob=scales::squish)
# ggsave(glue::glue('results/seq_lag_fit_ecdf_type{serotype_use}.png'),width=10,height=8,units='in')

# Plot subset
countries <- c("NIGERIA", "DEMOCRATIC REPUBLIC OF THE CONGO", "CENTRAL AFRICAN REPUBLIC", "YEMEN", "ETHIOPIA")
ggplot() + 
  geom_histogram(aes(y=after_stat(density),fill=admin0name,x=seq_lag),data=virus %>% filter(!is.na(seq_lag),
                                                                                            admin0name %in% countries))+
  geom_line(aes(x=day,y=dseq),colour='black',linetype=2,data=df_fit_eb %>% filter(admin0name %in% countries)) + 
  facet_wrap(vars(admin0name),scales = 'free_y')+
  theme(legend.position = 'none') + 
  scale_x_continuous('Days after collection',limits=c(0,max_lag),oob=scales::squish)

# Plot all countries
ggplot() + 
  geom_line(aes(x=day,y=dseq, group=admin0name),colour='black', linetype=1,data=df_fit_eb) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = 'none') +
  ylab("") +
  scale_x_continuous('Days after collection',limits=c(0,200),oob=scales::squish)

# Full bayesian treatment (not really neccesary since only using point estimates):
# fit <- brm(brmsformula(seq_lag ~ (1|admin0name), shape ~ (1|admin0name)),
#            cores = 4, 
#            data = virus %>% filter(!is.na(seq_lag),seq_lag >=0, seq_lag <365), 
#            family = Gamma(link='log'))
# 
# #summarize results:
# df_mu = coef(fit)$admin0name[,,'Intercept'] %>% as_tibble(rownames = 'admin0name') %>%
#   add_row(admin0name = 'other', 'Estimate' = fixef(fit)['Intercept','Estimate'])
# df_shape = coef(fit)$admin0name[,,"shape_Intercept"] %>% as_tibble(rownames = 'admin0name') %>%
#   add_row(admin0name = 'other', Estimate = fixef(fit)['shape_Intercept','Estimate'])
# 
# #compare to empirical bayes. Looks close enough:
# df_compare = left_join(fit_eb, 
#           df_mu %>% transmute(admin0name,lmu_brm = Estimate)) %>%
#   left_join(df_shape %>% transmute(admin0name, lshape_brm = Estimate)) 
# df_compare %>% filter(admin0name=='other')
# df_compare %>% arrange(desc(abs(lmu_brm-lmu)))
# df_compare %>% arrange(desc(abs(lshape_brm-lshape)))
# 
# 
# df_fit = expand_grid(admin0name =  c('other',unique(df_mu$admin0name)), day=0:200)
# 
# df_fit = df_fit %>% left_join(df_mu %>% select(admin0name,mu=Estimate)) %>%
#   left_join(df_shape %>% select(admin0name,shape=Estimate) ) %>% 
#   mutate(dseq = dgamma(day,rate = exp(shape)/exp(mu), shape=exp(shape)),
#          pseq = pgamma(day,rate = exp(shape)/exp(mu), shape=exp(shape)))
# 
# #plot
# ggplot() + 
#   geom_histogram(aes(y=after_stat(density),fill=admin0name,x=seq_lag),data=virus %>% filter(!is.na(seq_lag)))+
#   geom_line(aes(x=day,y=dseq),colour='black',linetype=2,data=df_fit) + 
#   facet_wrap(vars(admin0name),scales = 'free_y')+
#   theme(legend.position = 'none') + 
#   scale_x_continuous('Days after collection',limits=c(0,200),oob=scales::squish)
# ggsave('results/seq_lag_brms_fit.png',width=12,height=10,units='in')
# 
# ggplot() + 
#   geom_step(aes(colour=admin0name,x=seq_lag),stat='ecdf',data=virus %>% filter(!is.na(seq_lag), seq_lag < 365, seq_lag > 0))+
#   geom_line(aes(x=day,y=pseq),colour='black',linetype=2,data=df_fit) + 
#   facet_wrap(vars(admin0name))+
#   theme(legend.position = 'none') + 
#   scale_y_continuous('Proportion of Sequences Available',breaks = c(0,0.5,1))+
#   scale_x_continuous('Days after collection', limits=c(0,200),oob=scales::squish)
# ggsave('results/seq_lag_brms_fit_ecdf.png',width=10,height=8,units='in')

write_csv(df_fit_eb, file = filename_fit)
