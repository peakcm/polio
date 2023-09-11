full <- "C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub"
path_to_polis_data<-"/polio-immunity-mapping/results"
afp<-read.csv(paste0(full, path_to_polis_data,"/linelist_afp.csv"))
es<-read.csv(paste0(full, path_to_polis_data,"/linelist_es.csv"))
sia<-read.csv(paste0(full, path_to_polis_data,"/sia_district_rows.csv"))
afp$type<-"AFP"
es$type<-"ES"
library(tidyverse)
virus<-bind_rows(afp,es)
virus<-virus[!is.na(virus$vdpv_emergence_group)&virus$vdpv_emergence_group!="",]
virus$date<-as.Date(ifelse(is.na(virus$case_date),virus$collection_date,virus$case_date))
a<-virus[!grepl(",",virus$vdpv_emergence_group),]
b<-virus[grepl(",",virus$vdpv_emergence_group),]
b<-lapply(1:nrow(b),function(e){merge(b[e,which(names(b)!="vdpv_emergence_group")],cbind.data.frame(vdpv_emergence_group=strsplit(b$vdpv_emergence_group[e],", ")[[1]]))})
b<-do.call("rbind",b)
virus<-rbind(a,b)
#data as of 9/11 missing "cVDPV2" in polio_virus_types for BOT-FRA-1
virus[virus$vdpv_emergence_group == "BOT-FRA-1", "polio_virus_types"] <- "cVDPV2"
virus<-filter(virus,!iso3_code%in%c("AFG","PAK"),grepl("cVDPV2",polio_virus_types),vdpv_emergence_group!="cVDPV2")

virus<-virus[order(virus$date),]
virus<-mutate(group_by(as_tibble(virus),vdpv_emergence_group),start=min(date),end=max(date),n=n())

key<-summarise(group_by(as_tibble(virus),vdpv_emergence_group),iso3=unique(iso3_code[which(date==start)]))
key$iso3[key$vdpv_emergence_group=="SOM-AWL-1"]<-"SOM"

virus$grp<-paste(virus$vdpv_emergence_group,virus$iso3_code)
sub<-virus[virus$grp%in%(summarise(group_by(as_tibble(virus),grp),date=min(date))%>%filter(date>as.Date("2010-01-01")))$grp,]
summ<-summarise(group_by(as_tibble(sub),grp,adm0_name),start=min(date),end=max(date))
sub<-mutate(group_by(as_tibble(sub),grp),n=n())
sia_sub<-unique(sia[which(sia$ADM0_NAME%in%sub$adm0_name&sia$start_date>as.Date("2010-01-01")&sia$status!="Planned"&sia$vaccinetype%in%c("mOPV2",'nOPV2',"tOPV")),c("start_date","ADM0_NAME","vaccinetype")])
sia_sub<-merge(sia_sub,summ,by.x="ADM0_NAME",by.y="adm0_name",all=T)
sia_sub<-sia_sub[which(as.Date(sia_sub$start_date)>=sia_sub$start&as.Date(sia_sub$start_date)<=(sia_sub$end+6*30)),]

sub<-merge(sub,key)
sub$country_of_origin<-ifelse(sub$iso3==sub$iso3_code,"Inside country of origin","Outside country of origin")
sia_sub<-merge(sia_sub,unique(sub[,c("country_of_origin","grp")]),all.x=T)
df<-summarise(group_by(as_tibble(sub),grp,country_of_origin),start=min(date),end=max(date))
df$grp<-factor(df$grp,levels=df$grp[order(df$start)],ordered=T)
sub$grp<-factor(sub$grp,levels=df$grp[order(df$start)],ordered=T)
sia_sub$grp<-factor(sia_sub$grp,levels=df$grp[order(df$start)],ordered=T)

rect<-merge(cbind.data.frame(xmin=Sys.Date()-6*30,xmax=Sys.Date(),ymin=-Inf,ymax=Inf),cbind.data.frame(country_of_origin=unique(df$country_of_origin)),all=T)
p<-ggplot()+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="gray")+
  geom_segment(data=df,aes(x=start,xend=end,y=grp,yend=grp))+
  geom_point(data=sub[sub$n==1,],aes(x=date,y=grp))+
  geom_text(data=sia_sub,aes(x=as.Date(start_date),y=grp,label="|",color=vaccinetype))+
  theme(legend.position=c(0.1,0.9))+
  labs(x="Date",y="Country-lineage outbreak",color="SIA")+facet_wrap(~country_of_origin,scales="free_y")

library(ggpubr)
p1<-ggplot()+geom_histogram(data=df%>%filter(end<(Sys.Date()-6*30)),aes(x=as.numeric(end-start)/30),binwidth=1)+labs(y="Outbreaks",x="Duration (months)")+scale_y_continuous(limits=c(0,30))
d<-merge(sia_sub%>%filter(end<(Sys.Date()-6*30))%>%group_by(grp)%>%summarise(n=n()),cbind.data.frame(grp=(df%>%filter(end<(Sys.Date()-6*30)))$grp),all=T)
p2<-ggplot()+geom_histogram(data=d,aes(x=ifelse(is.na(n),0,n)),binwidth=1)+labs(y="Outbreaks",x="SIAs")+scale_y_continuous(limits=c(0,30))
q<-ggarrange(p1,p2,ncol=1)
ggarrange(p,q,ncol=2,widths=c(3,1))
ggsave("neo_epoch.png",width=10,height=8)

d<-merge(sia_sub%>%filter(end<(Sys.Date()-12*30))%>%group_by(grp)%>%summarise(n=n()),df%>%filter(end<(Sys.Date()-12*30)),all=T)
ggplot(data=d,aes(x=as.numeric((end+6*30)-start)/30,y=ifelse(is.na(n),0,n)))+geom_point(alpha=0.5)+geom_smooth(method="lm")+labs(x="Duration (months)",y="SIAs")

# Statements for pull up point analysis
df %>% View()
df %>% filter(end > Sys.Date()-365) %>% View()

# Number of lineage-country outbreaks active during each year
df <- df %>%
  mutate(
    active_2016 = if_else(start <= "2017-01-01" & end >= "2016-01-01",TRUE, FALSE),
    active_2017 = if_else(start <= "2018-01-01" & end >= "2017-01-01",TRUE, FALSE),
    active_2018 = if_else(start <= "2019-01-01" & end >= "2018-01-01",TRUE, FALSE),
    active_2019 = if_else(start <= "2020-01-01" & end >= "2019-01-01",TRUE, FALSE),
    active_2020 = if_else(start <= "2021-01-01" & end >= "2020-01-01",TRUE, FALSE),
    active_2021 = if_else(start <= "2022-01-01" & end >= "2021-01-01",TRUE, FALSE),
    active_2022 = if_else(start <= "2023-01-01" & end >= "2022-01-01",TRUE, FALSE),
    active_2023 = if_else(start <= "2024-01-01" & end >= "2023-01-01",TRUE, FALSE))

# Duration of lineage-country outbreaks
df$duration <- df$end - df$start
df %>% filter(duration == 0) %>% View()
median(df$duration)

df %>% 
  ungroup() %>%
  # group_by(duration >0) %>%
  # group_by(end > Sys.Date()-365) %>%
  group_by(duration >0, end > Sys.Date()-365) %>%
  summarize(
    count = n(),
    q1 = quantile(duration, 0.25),
    mid = quantile(duration, 0.5),
    q3 = quantile(duration, 0.75))

# Number of SIAs in response to outbreaks
  # Add a field for end date for that country-lineage outbreak
sia_sub$end <- NA
for (i in unique(sia_sub$grp)){
  sia_sub[sia_sub$grp == i, "end"] <- df[df$grp == i, "end"]
}

sia_sub %>% 
  filter(start_date <= end + 365) %>%
  filter(start_date < "2023-09-01") %>%
  group_by(vaccinetype) %>%
  group_by(year(start_date)) %>%
  summarize(
    count = n())

