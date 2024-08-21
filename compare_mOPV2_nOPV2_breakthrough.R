library(tidyverse)
library(scales)
library(sf)
library(zoo)
library(janitor)
library(lubridate)
library(units)
library(ggrepel)
library(Matrix)
library(RColorBrewer)   # improved plotting
library(cowplot)        # subplots
library(PolisAPI)
library(countrycode)
library(stanley)
breakthrough_days = 28

setwd("~/GitHub/polio-immunity-mapping/")
source('R/utils/polis_shapes_query_sf.R')
shape2 = st_read_geodb(config::get("filename_shapes"), 2, current_only = TRUE)
key2 = shape2 %>% st_drop_geometry %>% clean_names %>% select(guid,adm0_name,adm1_name,adm2_name,startdate,enddate)

as_of="2024-08-01"

sias = read_csv(config::get("filename_sia")) %>% clean_names()
nopv2_countries = sias %>% filter(vaccinetype=='nOPV2',start_date < as.Date(as_of)-months(1)) %>% .$adm0_name %>% unique
mopv2_countries = sias %>% filter(vaccinetype=='mOPV2',start_date < as.Date(as_of)-months(1)) %>% .$adm0_name %>% unique
topv_countries = sias %>% filter(vaccinetype=='tOPV',start_date < as.Date(as_of)-months(1),start_date > as.Date("2016-05-01")) %>% .$adm0_name %>% unique

nopv2_countries = nopv2_countries[!nopv2_countries%in%topv_countries]
mopv2_countries = mopv2_countries[!mopv2_countries%in%topv_countries]

sias%>%filter(adm0_name%in%intersect(nopv2_countries,mopv2_countries))%>%group_by(adm0_name)%>%
  summarise(min_nopv2=min(start_date[vaccinetype=="nOPV2"]),max_mopv2=max(start_date[vaccinetype=="mOPV2"]))

afp <- read_csv(config::get('filename_linelist_afp'))
es <- read_csv(config::get('filename_linelist_es'))
afp<-afp[afp$updated_date<=as.Date(as_of),]
es<-es[es$updated_date<=as.Date(as_of),]

viruses <- bind_rows(
  es %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = collection_date) %>% mutate(source = "ES"),
  afp %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = donset) %>% mutate(source = "AFP")
) %>%
  filter(!is.na(date)) %>%
  distinct(epid, date, .keep_all = TRUE)

viruses <- filter(viruses, adm0_name %in% nopv2_countries, date > ymd('2021-03-01'))

# Problem with Tajikistan SIAs 
st_intersects(st_point(c(68.8,  38.5)),shape2%>%filter(iso_3_code=="TJK"))
d<-merge(shape2%>%filter(iso_3_code=="TJK"),sias[sias$adm0_name=="TAJIKISTAN"&sias$vaccinetype=="nOPV2",])
ggplot()+geom_sf(data=shape2%>%filter(iso_3_code=="TJK"))+geom_sf(data=d,fill="red")+geom_sf(data=shape2%>%filter(adm2_name=="RUDAKI"),fill="blue",alpha=0.5)+facet_wrap(~start_date)

viruses$guid[viruses$guid=="{057000C8-26B9-4E9B-AF4A-86469C828BAA}"]<-"{802F778A-F8EB-45D1-B4BA-91FA9E31E5AB}"

opv2_tally = viruses %>%  
  left_join(sias %>% filter(vaccinetype %in% c('nOPV2')),by='guid') %>%
  group_by(epid) %>% summarise(opv2_before = sum(start_date < date-breakthrough_days & start_date >date %m-%months(12),na.rm=T),
                               opv2_after = sum(start_date > date-breakthrough_days,na.rm=T))

viruses = left_join(viruses,opv2_tally)
viruses = viruses %>% mutate(adm0_name_plot = recode(adm0_name,'DEMOCRATIC REPUBLIC OF THE CONGO'='DRC','CENTRAL AFRICAN REPUBLIC'='CAR'))

bwidth=30
max.SIAs = 4
wk_seq = seq(min(viruses$date),today(),by=bwidth)
afp_wk = viruses %>% filter(source=='AFP') %>% transmute(adm0_name,adm0_name_plot,adm1_name,
                                                         date=wk_seq[findInterval(date,wk_seq)],
                                                         breakthrough=factor(pmin(opv2_before,max.SIAs),levels=c(0:max.SIAs),labels=paste0(0:max.SIAs,ifelse(0:max.SIAs==max.SIAs,"+","")),ordered = TRUE))
es_wk = viruses %>% filter(source == 'ES') %>% transmute(adm0_name,adm0_name_plot,adm1_name,
                                                         date=wk_seq[findInterval(date,wk_seq)],
                                                         breakthrough=factor(pmin(opv2_before,max.SIAs),levels=c(0:max.SIAs),labels=paste0(0:max.SIAs,ifelse(0:max.SIAs==max.SIAs,"+","")),ordered = TRUE))


maxn = afp_wk %>% group_by(adm0_name,date) %>% tally() %>% summarise(maxn=max(n))
minn = es_wk %>% group_by(adm0_name,date) %>% tally() %>% summarise(maxn=max(n))
sia_adm0 = sias %>% filter(adm0_name %in% nopv2_countries,vaccinetype=='nOPV2') %>% distinct(adm0_name,start_date) %>%
  left_join(maxn) %>% left_join(minn) %>% replace_na(list(maxn=1,minn=0)) %>%
  mutate(adm0_name_plot =  recode(adm0_name,'DEMOCRATIC REPUBLIC OF THE CONGO'='DRC','CENTRAL AFRICAN REPUBLIC'='CAR'))

key<-summarise(group_by(bind_rows(afp_wk,es_wk),adm0_name),n=max(breakthrough),latest=max(date))

sel<-key$adm0_name[key$n==0]
p1<-ggplot()+
  geom_segment(aes(x=start_date,xend=start_date,yend=pmax(2,maxn),y=pmax(2,maxn)*1.1),size=0.5,colour='red',alpha=0.8,
               data=sia_adm0%>%filter(adm0_name%in%sel),arrow=arrow(length=unit(0.3,'cm'),type = 'closed',angle = 20))+
  geom_histogram(data=afp_wk %>% filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(x=date,fill=breakthrough),
                 binwidth = bwidth,boundary=mdy('01012021'))+
  geom_histogram(data=es_wk %>%filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(y=-pmin(50,stat(count)),x=date,fill=breakthrough),
                 binwidth = bwidth,alpha=0.7,boundary=mdy('01012021'))+
  theme_bw()+
  scale_x_date(NULL,date_labels = '%y',date_breaks = '1 year', limits = c(min(viruses$date),today()))+
  scale_fill_brewer('Preceding\nnOPV2 SIAs',palette = 'RdYlBu',direction = -1, drop=F)+
  scale_y_continuous('', labels = function(x){abs(x)})+
  expand_limits(y=c(-2,2))+
  facet_wrap(~adm0_name_plot,scales = 'free_y',nrow=1) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = today(),linetype=2,colour='grey50')

sel<-key$adm0_name[key$n==1]
p2<-ggplot()+
  geom_segment(aes(x=start_date,xend=start_date,yend=pmax(2,maxn),y=pmax(2,maxn)*1.1),size=0.5,colour='red',alpha=0.8,
               data=sia_adm0%>%filter(adm0_name%in%sel),arrow=arrow(length=unit(0.3,'cm'),type = 'closed',angle = 20))+
  geom_histogram(data=afp_wk %>% filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(x=date,fill=breakthrough),
                 binwidth = bwidth,boundary=mdy('01012021'))+
  geom_histogram(data=es_wk %>%filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(y=-pmin(50,stat(count)),x=date,fill=breakthrough),
                 binwidth = bwidth,alpha=0.7,boundary=mdy('01012021'))+
  theme_bw()+guides(fill="none")+
  scale_x_date(NULL,date_labels = '%y',date_breaks = '1 year', limits = c(min(viruses$date),today()))+
  scale_fill_brewer('Preceding\nnOPV2 SIAs',palette = 'RdYlBu',direction = -1, drop=F)+
  scale_y_continuous('', labels = function(x){abs(x)})+
  expand_limits(y=c(-2,2))+
  facet_wrap(~adm0_name_plot,scales = 'free_y',nrow=1) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = today(),linetype=2,colour='grey50')

sel<-key$adm0_name[key$n==2]
p3<-ggplot()+
  geom_segment(aes(x=start_date,xend=start_date,yend=pmax(2,maxn),y=pmax(2,maxn)*1.1),size=0.5,colour='red',alpha=0.8,
               data=sia_adm0%>%filter(adm0_name%in%sel),arrow=arrow(length=unit(0.3,'cm'),type = 'closed',angle = 20))+
  geom_histogram(data=afp_wk %>% filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(x=date,fill=breakthrough),
                 binwidth = bwidth,boundary=mdy('01012021'))+
  geom_histogram(data=es_wk %>%filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(y=-pmin(50,stat(count)),x=date,fill=breakthrough),
                 binwidth = bwidth,alpha=0.7,boundary=mdy('01012021'))+
  theme_bw()+guides(fill="none")+
  scale_x_date(NULL,date_labels = '%y',date_breaks = '1 year', limits = c(min(viruses$date),today()))+
  scale_fill_brewer('Preceding\nnOPV2 SIAs',palette = 'RdYlBu',direction = -1, drop=F)+
  scale_y_continuous('', labels = function(x){abs(x)})+
  expand_limits(y=c(-2,2))+
  facet_wrap(~adm0_name_plot,scales = 'free_y',nrow=1) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = today(),linetype=2,colour='grey50')

sel<-key$adm0_name[key$n==3]
p4<-ggplot()+
  geom_segment(aes(x=start_date,xend=start_date,yend=pmax(2,maxn),y=pmax(2,maxn)*1.1),size=0.5,colour='red',alpha=0.8,
               data=sia_adm0%>%filter(adm0_name%in%sel),arrow=arrow(length=unit(0.3,'cm'),type = 'closed',angle = 20))+
  geom_histogram(data=afp_wk %>% filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(x=date,fill=breakthrough),
                 binwidth = bwidth,boundary=mdy('01012021'))+
  geom_histogram(data=es_wk %>%filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(y=-pmin(50,stat(count)),x=date,fill=breakthrough),
                 binwidth = bwidth,alpha=0.7,boundary=mdy('01012021'))+
  theme_bw()+guides(fill="none")+
  scale_x_date(NULL,date_labels = '%y',date_breaks = '1 year', limits = c(min(viruses$date),today()))+
  scale_fill_brewer('Preceding\nnOPV2 SIAs',palette = 'RdYlBu',direction = -1, drop=F)+
  scale_y_continuous('', labels = function(x){abs(x)})+
  expand_limits(y=c(-2,2))+
  facet_wrap(~adm0_name_plot,scales = 'free_y',nrow=1) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = today(),linetype=2,colour='grey50')

sel<-key$adm0_name[key$n=="4+"]
p5<-ggplot()+
  geom_segment(aes(x=start_date,xend=start_date,yend=pmax(2,maxn),y=pmax(2,maxn)*1.1),size=0.5,colour='red',alpha=0.8,
               data=sia_adm0%>%filter(adm0_name%in%sel),arrow=arrow(length=unit(0.3,'cm'),type = 'closed',angle = 20))+
  geom_histogram(data=afp_wk %>% filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(x=date,fill=breakthrough),
                 binwidth = bwidth,boundary=mdy('01012021'))+
  geom_histogram(data=es_wk %>%filter(adm0_name%in%sel)%>%
                   mutate(adm0_name = if_else(adm0_name=='DEMOCRATIC REPUBLIC OF THE CONGO','DRC',adm0_name)),
                 aes(y=-pmin(50,stat(count)),x=date,fill=breakthrough),
                 binwidth = bwidth,alpha=0.7,boundary=mdy('01012021'))+
  theme_bw()+guides(fill="none")+
  scale_x_date(NULL,date_labels = '%y',date_breaks = '1 year', limits = c(min(viruses$date),today()))+
  scale_fill_brewer('Preceding\nnOPV2 SIAs',palette = 'RdYlBu',direction = -1, drop=F)+
  scale_y_continuous('', labels = function(x){abs(x)})+
  expand_limits(y=c(-2,2))+
  facet_wrap(~adm0_name_plot,scales = 'free_y',nrow=1) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept = today(),linetype=2,colour='grey50')

table(key$n)
key%>%filter(latest>Sys.Date()-90)

ggsave('results/nopv2_breakthrough_cases_epicurve1.png',p1,width=17.5,height=1.8,units='in',dpi=300)
ggsave('results/nopv2_breakthrough_cases_epicurve2.png',p2,width=18,height=1.8,units='in',dpi=300)
ggsave('results/nopv2_breakthrough_cases_epicurve3.png',p3,width=18/10*6,height=1.8,units='in',dpi=300)
ggsave('results/nopv2_breakthrough_cases_epicurve4.png',p4,width=5.5,height=1.8,units='in',dpi=300)
ggsave('results/nopv2_breakthrough_cases_epicurve5.png',p5,width=2.1,height=1.8,units='in',dpi=300)

# Compare nOPV2, mOPV2
viruses <- bind_rows(
  es %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = collection_date) %>% mutate(source = "ES"),
  afp %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = donset) %>% mutate(source = "AFP")
) %>%
  filter(!is.na(date)) %>%
  distinct(epid, date, .keep_all = TRUE)
viruses$guid[viruses$guid=="{057000C8-26B9-4E9B-AF4A-86469C828BAA}"]<-"{802F778A-F8EB-45D1-B4BA-91FA9E31E5AB}"

viruses <- filter(viruses, adm0_name %in% nopv2_countries, date > ymd('2021-03-01'))
opv2_tally = viruses %>%  
  left_join(sias %>% filter(vaccinetype %in% c('nOPV2')),by='guid') %>%
  group_by(epid) %>% summarise(nearest_after = min(as.numeric(start_date-(date-breakthrough_days))[start_date > (date-breakthrough_days)]), 
                               opv2_before = sum(start_date < date-breakthrough_days & start_date >date %m-%months(12),na.rm=T),
                               opv2_after = sum(start_date > date-breakthrough_days,na.rm=T))

viruses = left_join(viruses,opv2_tally)
viruses = viruses %>% mutate(adm0_name_plot = recode(adm0_name,'DEMOCRATIC REPUBLIC OF THE CONGO'='DRC','CENTRAL AFRICAN REPUBLIC'='CAR'))

key_nOPV2<-summarise(group_by(viruses,adm0_name),n=max(opv2_before),after=min(nearest_after[opv2_before==n]))

viruses <- bind_rows(
  es %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = collection_date) %>% mutate(source = "ES"),
  afp %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = donset) %>% mutate(source = "AFP")
) %>%
  filter(!is.na(date)) %>%
  distinct(epid, date, .keep_all = TRUE)

viruses <- filter(viruses, adm0_name %in% mopv2_countries, date < ymd('2021-03-01'), date > ymd('2016-05-01'))
opv2_tally = viruses %>%  
  left_join(sias %>% filter(vaccinetype %in% c('mOPV2')),by='guid') %>%
  group_by(epid) %>% summarise(nearest_after = min(as.numeric(start_date-(date-breakthrough_days))[start_date > (date-breakthrough_days)]), 
                               opv2_before = sum(start_date < date-breakthrough_days & start_date >date %m-%months(12),na.rm=T),
                               opv2_after = sum(start_date > date-breakthrough_days,na.rm=T))

viruses = left_join(viruses,opv2_tally)
viruses = viruses %>% mutate(adm0_name_plot = recode(adm0_name,'DEMOCRATIC REPUBLIC OF THE CONGO'='DRC','CENTRAL AFRICAN REPUBLIC'='CAR'))
key_mOPV2<-summarise(group_by(viruses,adm0_name),n=max(opv2_before),after=min(nearest_after[opv2_before==n]))

k<-rbind(cbind(key_mOPV2,vac="mOPV2"),cbind(key_nOPV2,vac="nOPV2"))
k%>%filter(after>180)
k$n_num<-ifelse((k$n+1)>4,5,(k$n+1))
k$n_char<-ifelse((k$n+1)>4,"5+",(k$n+1))
summarise(group_by(k,n_num,n_char,vac),countries=n())%>%group_by(vac)%>%mutate(total=sum(countries))%>%
  ggplot()+geom_col(aes(x=n_char,y=countries/total,fill=vac),position=position_dodge(width=0.5,preserve="single"))+
  geom_text(aes(x=n_num+ifelse(vac=="mOPV2",-0.2,0.2),y=0.02+countries/total,label=countries))+theme_bw()+theme(legend.position=c(0.8,0.8))+
  scale_y_continuous(labels=scales::percent)+labs(x="SIAs before interruption of cVDPV2",y="Proportion of total countries (%)",fill="Vaccine")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"PRGn")[c(1,3)])
ggsave('results/compare_nopv2_mopv2_breakthrough.png',width=4,height=4,units='in',dpi=300)

library(ggrepel)
library(ggbeeswarm)
key_mOPV2$n_num<-ifelse((key_mOPV2$n+1)>4,5,(key_mOPV2$n+1))
key_mOPV2$n_char<-ifelse((key_mOPV2$n+1)>4,"5+",(key_mOPV2$n+1))
key_nOPV2$n_num<-ifelse((key_nOPV2$n+1)>4,5,(key_nOPV2$n+1))
key_nOPV2$n_char<-ifelse((key_nOPV2$n+1)>4,"5+",(key_nOPV2$n+1))

d<-merge(key_mOPV2,key_nOPV2,by="adm0_name",suffixes=c(".m",".n"),all=T)%>%
  mutate(adm0_name_plot =  recode(str_to_title(adm0_name),'Democratic Republic Of The Congo'='DRC','Central African Republic'='CAR',"Côte D’ivoire"="Côte d’Ivoire"),n.m=as.numeric(n.m),n.n=as.numeric(n.n))
d$class<-ifelse(is.na(d$n.m),"nOPV2 only",ifelse(is.na(d$n.n),"mOPV2 only","nOPV2 & mOPV2"))
ggplot()+geom_beeswarm(data=d%>%filter(class=="nOPV2 & mOPV2"),aes(x=n_char.m,y=n_char.n,color=adm0_name_plot))+
  geom_text_repel(data=d%>%filter(class=="nOPV2 & mOPV2"),aes(x=n_num.m,y=n_num.n,label=adm0_name_plot,color=adm0_name_plot))+
  labs(x="mOPV2 SIAs before interruption of cVDPV2",y="nOPV2 SIAs before interruption of cVDPV2")+guides(color="none")+theme_bw()
ggsave('results/compare_nopv2_mopv2_breakthrough_scatter.png',width=4,height=4,units='in',dpi=300)


d$class<-factor(d$class,names(table(d$class)),paste0(names(table(d$class))," (n=",as.numeric(table(d$class)),")"))

e<-merge(rbind(cbind(key_mOPV2,vac="mOPV2"),cbind(key_nOPV2,vac="nOPV2")),d[,c("class","adm0_name")])
e<-summarise(group_by(e,n,vac,class),countries=n())%>%group_by(vac,class)%>%mutate(total=sum(countries),cum=cumsum(countries))


ggplot()+geom_col(data=e,aes(x=ifelse(as.numeric(n+1)>4,"5+",as.numeric(n+1)),y=countries/total,fill=vac),position=position_dodge(width=0.5,preserve="single"))+
  geom_text(data=e,aes(x=ifelse(as.numeric(n+1)>4,5,as.numeric(n+1)),y=countries/total+0.02,label=countries),size=3)+theme_bw()+theme(legend.position=c(0.85,0.8))+
  scale_y_continuous(labels=scales::percent)+labs(x="SIAs before interruption of cVDPV2",y="Proportion of total countries (%)",fill="Vaccine")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"PRGn")[c(1,3)])+facet_wrap(~class)
ggsave('results/compare_nopv2_mopv2_breakthrough_hist.png',width=5,height=3,units='in',dpi=300)

e<-merge(rbind(cbind(key_mOPV2,vac="mOPV2"),cbind(key_nOPV2,vac="nOPV2")),d[,c("class","adm0_name")])
load("~/polio-immunity-mapping/data_local/adm0data.Rdata")
map<-st_as_sf(adm0data[adm0data$ENDDATE>as.Date("2030-01-01"),])
map<-merge(map,summarise(group_by(e,adm0_name),n_both=max(n_num)),by.x="ADM0_NAME",by.y="adm0_name",all.x=T)
map<-merge(map,summarise(group_by(e%>%filter(vac=="mOPV2"),adm0_name),n_mopv2=max(n_num)),by.x="ADM0_NAME",by.y="adm0_name",all.x=T)
map<-merge(map,summarise(group_by(e%>%filter(vac=="nOPV2"),adm0_name),n_nopv2=max(n_num)),by.x="ADM0_NAME",by.y="adm0_name",all.x=T)
map<-pivot_longer(map,cols=c(n_both,n_mopv2,n_nopv2),values_to="n")
b<-st_bbox(map[!is.na(map$n),])
map$name<-factor(map$name,levels=c("n_mopv2","n_nopv2","n_both"),labels=c("mOPV2","nOPV2","nOPV2 or mOPV2"),ordered=T)
# map$n[map$n==5]<-"5+"
map$n_cat<-cut(map$n,c(1,3,5,6),right=F,labels=c("1-2","3-4","5+"))
p<-ggplot()+geom_sf(data=map%>%filter(name%in%"nOPV2 or mOPV2"),aes(fill=as.factor(n_cat)))+coord_sf(xlim=b[c(1,3)],ylim=b[c(2,4)])+facet_wrap(~name)+labs(fill="SIAs before\ninterruption\nof cVDPV2")
ggsave("map_sias_before_int.png",p,width=10,height=10)

# GUID level
viruses <- bind_rows(
  es %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = collection_date) %>% mutate(source = "ES"),
  afp %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = donset) %>% mutate(source = "AFP")
) %>%
  filter(!is.na(date)) %>%
  distinct(epid, date, .keep_all = TRUE)

viruses <- filter(viruses, adm0_name %in% nopv2_countries, date > ymd('2021-03-01'))
opv2_tally = viruses %>%  
  left_join(sias %>% filter(vaccinetype %in% c('nOPV2')),by='guid') %>%
  group_by(epid) %>% summarise(nearest_after = min(as.numeric(start_date-(date-breakthrough_days))[start_date > (date-breakthrough_days)]), 
                               opv2_before = sum(start_date < date-breakthrough_days & start_date >date %m-%months(12),na.rm=T),
                               opv2_after = sum(start_date > date-breakthrough_days,na.rm=T))

viruses = left_join(viruses,opv2_tally)
viruses = viruses %>% mutate(adm0_name_plot = recode(adm0_name,'DEMOCRATIC REPUBLIC OF THE CONGO'='DRC','CENTRAL AFRICAN REPUBLIC'='CAR'))

key_nOPV2<-summarise(group_by(viruses,adm0_name,guid),n=max(opv2_before),after=min(nearest_after[opv2_before==n]))

viruses <- bind_rows(
  es %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = collection_date) %>% mutate(source = "ES"),
  afp %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = donset) %>% mutate(source = "AFP")
) %>%
  filter(!is.na(date)) %>%
  distinct(epid, date, .keep_all = TRUE)

viruses <- filter(viruses, adm0_name %in% mopv2_countries, date < ymd('2021-03-01'), date > ymd('2016-05-01'))
opv2_tally = viruses %>%  
  left_join(sias %>% filter(vaccinetype %in% c('mOPV2')),by='guid') %>%
  group_by(epid) %>% summarise(nearest_after = min(as.numeric(start_date-(date-breakthrough_days))[start_date > (date-breakthrough_days)]), 
                               opv2_before = sum(start_date < date-breakthrough_days & start_date >date %m-%months(12),na.rm=T),
                               opv2_after = sum(start_date > date-breakthrough_days,na.rm=T))

viruses = left_join(viruses,opv2_tally)
viruses = viruses %>% mutate(adm0_name_plot = recode(adm0_name,'DEMOCRATIC REPUBLIC OF THE CONGO'='DRC','CENTRAL AFRICAN REPUBLIC'='CAR'))
key_mOPV2<-summarise(group_by(viruses,adm0_name,guid),n=max(opv2_before),after=min(nearest_after[opv2_before==n]))

k<-rbind(cbind(key_mOPV2,vac="mOPV2"),cbind(key_nOPV2,vac="nOPV2"))
k$n_num<-ifelse((k$n+1)>4,5,(k$n+1))
k$n_char<-ifelse((k$n+1)>4,"5+",(k$n+1))
summarise(group_by(k,n_num,n_char,vac),countries=n())%>%group_by(vac)%>%mutate(total=sum(countries))%>%
  ggplot()+geom_col(aes(x=n_char,y=countries/total,fill=vac),position=position_dodge(width=0.5,preserve="single"))+
  geom_text(aes(x=n_num+ifelse(vac=="mOPV2",-0.2,0.2),y=0.02+countries/total,label=countries))+theme_bw()+theme(legend.position=c(0.8,0.8))+
  scale_y_continuous(labels=scales::percent)+labs(x="SIAs before interruption of cVDPV2",y="Proportion of total countries (%)",fill="Vaccine")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"PRGn")[c(1,3)])

key_mOPV2$n_num<-ifelse((key_mOPV2$n+1)>4,5,(key_mOPV2$n+1))
key_mOPV2$n_char<-ifelse((key_mOPV2$n+1)>4,"5+",(key_mOPV2$n+1))
key_nOPV2$n_num<-ifelse((key_nOPV2$n+1)>4,5,(key_nOPV2$n+1))
key_nOPV2$n_char<-ifelse((key_nOPV2$n+1)>4,"5+",(key_nOPV2$n+1))

d<-merge(key_mOPV2,key_nOPV2,by=c("guid","adm0_name"),suffixes=c(".m",".n"),all=T)%>%
  mutate(adm0_name_plot =  recode(str_to_title(adm0_name),'Democratic Republic Of The Congo'='DRC','Central African Republic'='CAR',"Côte D’ivoire"="Côte d’Ivoire"),n.m=as.numeric(n.m),n.n=as.numeric(n.n))
d$class<-ifelse(is.na(d$n.m),"nOPV2 only",ifelse(is.na(d$n.n),"mOPV2 only","nOPV2 & mOPV2"))
ggplot()+geom_jitter(data=d%>%filter(class=="nOPV2 & mOPV2"),aes(x=n_char.m,y=n_char.n,color=adm0_name_plot))+
  # geom_text_repel(data=d%>%filter(class=="nOPV2 & mOPV2"),aes(x=n_num.m,y=n_num.n,label=adm0_name_plot,color=adm0_name_plot))+
  labs(x="mOPV2 SIAs before interruption of cVDPV2",y="nOPV2 SIAs before interruption of cVDPV2")+guides(color="none")+theme_bw()


d$class<-factor(d$class,names(table(d$class)),paste0(names(table(d$class))," (n=",as.numeric(table(d$class)),")"))

e<-merge(rbind(cbind(key_mOPV2,vac="mOPV2"),cbind(key_nOPV2,vac="nOPV2")),d[,c("class","guid","adm0_name")])
e<-summarise(group_by(e,n_num,n_char,vac,class),countries=n())%>%group_by(vac,class)%>%mutate(total=sum(countries),cum=cumsum(countries))


ggplot()+geom_col(data=e,aes(x=n_char,y=countries/total,fill=vac),position=position_dodge(width=0.5,preserve="single"))+
  geom_text(data=e,aes(x=n_num,y=countries/total+0.02,label=countries),size=3)+theme_bw()+theme(legend.position=c(0.85,0.8))+
  scale_y_continuous(labels=scales::percent)+labs(x="SIAs before interruption of cVDPV2",y="Proportion of total ADM2 (%)",fill="Vaccine")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"PRGn")[c(1,3)])+facet_wrap(~class)

e<-merge(rbind(cbind(key_mOPV2,vac="mOPV2"),cbind(key_nOPV2,vac="nOPV2")),d[,c("class","guid","adm0_name")])
load("~/polio-immunity-mapping/data_local/adm2data.Rdata")
map<-st_as_sf(adm2data[adm2data$ENDDATE>as.Date("2030-01-01")&adm2data$ADM0_NAME%in%e$adm0_name,])
map<-merge(map,summarise(group_by(e,guid,adm0_name),n_both=max(n)),by.x="GUID",by.y="guid",all.x=T)
map$n_cat<-cut(map$n_both+1,c(1,3,5,9),right=F,labels=c("1-2","3-4","5-8"))
b<-st_bbox(map[!is.na(map$n_both),])
map0<-st_as_sf(adm0data[adm0data$ENDDATE>as.Date("2030-01-01"),])
p<-ggplot()+geom_sf(data=map0,fill="grey")+geom_sf(data=map,aes(fill=as.factor(n_cat)),linetype=0)+geom_sf(data=map0,fill="transparent")+coord_sf(xlim=b[c(1,3)],ylim=b[c(2,4)])+labs(fill="SIAs before\ninterruption\nof cVDPV2")
ggsave("map_sias_before_int_adm2.png",p,width=10,height=10)


# ADM1 level
viruses <- bind_rows(
  es %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = collection_date) %>% mutate(source = "ES"),
  afp %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = donset) %>% mutate(source = "AFP")
) %>%
  filter(!is.na(date)) %>%
  distinct(epid, date, .keep_all = TRUE)

viruses <- filter(viruses, adm0_name %in% "NIGERIA", date > ymd('2021-03-01'))
opv2_tally = viruses %>%  
  left_join(sias %>% filter(vaccinetype %in% c('nOPV2')),by='guid') %>%
  group_by(epid) %>% summarise(nearest_after = min(as.numeric(start_date-(date-breakthrough_days))[start_date > (date-breakthrough_days)]), 
                               opv2_before = sum(start_date < date-breakthrough_days & start_date >date %m-%months(12),na.rm=T),
                               opv2_after = sum(start_date > date-breakthrough_days,na.rm=T))

viruses = left_join(viruses,opv2_tally)
viruses = viruses %>% mutate(adm0_name_plot = recode(adm0_name,'DEMOCRATIC REPUBLIC OF THE CONGO'='DRC','CENTRAL AFRICAN REPUBLIC'='CAR'))

key_nOPV2<-summarise(group_by(viruses,adm0_name,adm1_name),n=max(opv2_before),after=min(nearest_after[opv2_before==n]))

viruses <- bind_rows(
  es %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = collection_date) %>% mutate(source = "ES"),
  afp %>% filter(final_class == "cVDPV2") %>% select(epid, guid, matches("adm[0-2]_name"), date = donset) %>% mutate(source = "AFP")
) %>%
  filter(!is.na(date)) %>%
  distinct(epid, date, .keep_all = TRUE)

viruses <- filter(viruses, adm0_name %in% "NIGERIA", date < ymd('2021-01-01'), date > ymd('2016-05-01'))
opv2_tally = viruses %>%  
  left_join(sias %>% filter(vaccinetype %in% c('mOPV2')),by='guid') %>%
  group_by(epid) %>% summarise(nearest_after = min(as.numeric(start_date-(date-breakthrough_days))[start_date > (date-breakthrough_days)]), 
                               opv2_before = sum(start_date < date-breakthrough_days & start_date >date %m-%months(12),na.rm=T),
                               opv2_after = sum(start_date > date-breakthrough_days,na.rm=T))

viruses = left_join(viruses,opv2_tally)
viruses = viruses %>% mutate(adm0_name_plot = recode(adm0_name,'DEMOCRATIC REPUBLIC OF THE CONGO'='DRC','CENTRAL AFRICAN REPUBLIC'='CAR'))

key_mOPV2<-summarise(group_by(viruses,adm0_name,adm1_name),n=max(opv2_before),after=min(nearest_after[opv2_before==n]))

k<-rbind(cbind(key_mOPV2,vac="mOPV2"),cbind(key_nOPV2,vac="nOPV2"))
k$n_num<-ifelse((k$n+1)>4,5,(k$n+1))
k$n_char<-ifelse((k$n+1)>4,"5+",(k$n+1))
summarise(group_by(k,n_num,n_char,vac),countries=n())%>%group_by(vac)%>%mutate(total=sum(countries))%>%
  ggplot()+geom_col(aes(x=n_char,y=countries/total,fill=vac),position=position_dodge(width=0.5,preserve="single"))+
  geom_text(aes(x=n_num+ifelse(vac=="mOPV2",-0.2,0.2),y=0.02+countries/total,label=countries))+theme_bw()+theme(legend.position=c(0.8,0.8))+
  scale_y_continuous(labels=scales::percent)+labs(x="SIAs before interruption of cVDPV2",y="Proportion of total countries (%)",fill="Vaccine")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"PRGn")[c(1,3)])

key_mOPV2$n_num<-ifelse((key_mOPV2$n+1)>4,5,(key_mOPV2$n+1))
key_mOPV2$n_char<-ifelse((key_mOPV2$n+1)>4,"5+",(key_mOPV2$n+1))
key_nOPV2$n_num<-ifelse((key_nOPV2$n+1)>4,5,(key_nOPV2$n+1))
key_nOPV2$n_char<-ifelse((key_nOPV2$n+1)>4,"5+",(key_nOPV2$n+1))

d<-merge(key_mOPV2,key_nOPV2,by=c("adm1_name","adm0_name"),suffixes=c(".m",".n"),all=T)%>%
  mutate(adm0_name_plot =  recode(str_to_title(adm0_name),'Democratic Republic Of The Congo'='DRC','Central African Republic'='CAR',"Côte D’ivoire"="Côte d’Ivoire"),n.m=as.numeric(n.m),n.n=as.numeric(n.n))
d$class<-ifelse(is.na(d$n.m),"nOPV2 only",ifelse(is.na(d$n.n),"mOPV2 only","nOPV2 & mOPV2"))
d$adm1<-tolower(d$adm1_name)
d<-merge(d,read.csv("zone_key.csv"))
ggplot()+geom_beeswarm(data=d%>%filter(class=="nOPV2 & mOPV2"),aes(x=n_char.m,y=n_char.n))+geom_abline()+
  # geom_text_repel(data=d%>%filter(class=="nOPV2 & mOPV2"),aes(x=n_num.m,y=n_num.n,label=adm0_name_plot,color=adm0_name_plot))+
  labs(x="mOPV2 SIAs before interruption of cVDPV2",y="nOPV2 SIAs before interruption of cVDPV2")+theme_bw()

ggplot()+geom_beeswarm(data=d%>%filter(class=="nOPV2 & mOPV2"),aes(x=n.m,y=n.n,color=Zone))+geom_abline()+
  # geom_text_repel(data=d%>%filter(class=="nOPV2 & mOPV2"),aes(x=n_num.m,y=n_num.n,label=adm0_name_plot,color=adm0_name_plot))+
  labs(x="mOPV2 SIAs before interruption of cVDPV2",y="nOPV2 SIAs before interruption of cVDPV2")+theme_bw()



d%>%filter(class=="nOPV2 & mOPV2",n_num.n>3,n_num.m<3)
d%>%filter(class=="nOPV2 & mOPV2",n_num.n>3,n_num.m>3)
viruses%>%filter(adm1_name=="LAGOS")
sias %>% filter(adm1_name %in% c('LAGOS'),vaccinetype %in% c('mOPV2'))%>%group_by(adm1_name,start_date)%>%summarise()
viruses%>%filter(adm1_name=="SOKOTO")
sias %>% filter(adm1_name %in% c('SOKOTO'),vaccinetype %in% c('mOPV2'))%>%group_by(adm1_name,start_date)%>%summarise()
viruses%>%filter(adm1_name=="BORNO")
sias %>% filter(adm1_name %in% c('BORNO'),vaccinetype %in% c('mOPV2'))%>%group_by(adm1_name,start_date)%>%summarise()

d$class<-factor(d$class,names(table(d$class)),paste0(names(table(d$class))," (n=",as.numeric(table(d$class)),")"))

e<-merge(rbind(cbind(key_mOPV2,vac="mOPV2"),cbind(key_nOPV2,vac="nOPV2")),d[,c("class","guid","adm0_name")])
e<-summarise(group_by(e,n_num,n_char,vac,class),countries=n())%>%group_by(vac,class)%>%mutate(total=sum(countries),cum=cumsum(countries))


ggplot()+geom_col(data=e,aes(x=n_char,y=countries/total,fill=vac),position=position_dodge(width=0.5,preserve="single"))+
  geom_text(data=e,aes(x=n_num,y=countries/total+0.02,label=countries),size=3)+theme_bw()+theme(legend.position=c(0.85,0.8))+
  scale_y_continuous(labels=scales::percent)+labs(x="SIAs before interruption of cVDPV2",y="Proportion of total countries (%)",fill="Vaccine")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"PRGn")[c(1,3)])+facet_wrap(~class)

