library(PolisAPI)
library(lubridate)
library(tidyverse)
# Vector of emergence groups linked to nOPV2
n<-c("RDC-SKV-1", "RDC-TAN-2", "RDC-KOR-1",
      "CAF-KEM-1", "NIE-KBS-1", "RDC-HKA-2",
      "CAF-BNG-3", "BOT-FRA-1", "EGY-NOR-1",
      "CAE-EXT-1", "ZIM-HRE-1", "NIE-KTS-1",
      "MOZ-MAN-1") # holding out RSS-WEQ-1 until confirmed


set_token("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/data_local/token.txt")
viruses_raw = get_polis_virus(min_date = '2014-01-01',virus_id = 4)
viruses_raw_save<-viruses_raw

viruses_raw$year<-year(viruses_raw$virus_date)
viruses_raw$half_year<-floor_date(as.Date(viruses_raw$virus_date),"6 months")
viruses_raw<-viruses_raw[viruses_raw$vdpv_classification_name=="Circulating"&!is.na(viruses_raw$vdpv_emergence_group_name),]
viruses_raw$origin<-as.Date(viruses_raw$virus_date)-30.5*as.numeric(viruses_raw$nt_changes)
viruses_raw$vacc<-ifelse(viruses_raw$vdpv_emergence_group_name%in%n,"nOPV2",'Sabin2')
a<-mutate(group_by(viruses_raw[order(viruses_raw$virus_date),],vdpv_emergence_group_name,vacc),rank=1:n())%>%filter(rank<4)%>%summarise(origin=mean(origin,na.rm=T))
a$year<-year(a$origin)
a$half_year<-floor_date(as.Date(a$origin),"6 months")

full <- "C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub"
path_to_polis_data<-"/polio-immunity-mapping/results"
sia<-read.csv(paste0(full, path_to_polis_data,"/sia_district_rows_raw.csv"))
sia<-sia[sia$VaccineType%in%c("mOPV2","nOPV2","tOPV")&sia$Status%in%c("Done"),]
summary(sia$CalculatedDosages)
summary(sia[is.na(sia$CalculatedDosages),])
sia$dose<-sia$CalculatedDosages
sia$dose[is.na(sia$CalculatedDosages)]<-(sia$CalculatedTargetPopulation*sia$WastageFactor)[is.na(sia$CalculatedDosages)]

# afp<-read.csv(paste0(path_to_polis_data,"/linelist_afp.csv"))
# es<-read.csv(paste0(path_to_polis_data,"/linelist_es.csv"))
#
# a<-es[!grepl(",",es$vdpv_emergence_group),]
# b<-es[grepl(",",es$vdpv_emergence_group),]
# b<-lapply(1:nrow(b),function(e){merge(b[e,which(!names(b)%in%c("vdpv_emergence_group","vdpv_nt_changes"))],
#                                       cbind.data.frame(vdpv_emergence_group=strsplit(b$vdpv_emergence_group[e],", ")[[1]],
#                                                        vdpv_nt_changes=strsplit(b$vdpv_nt_changes[e],", ")[[1]]))})
# b<-do.call("rbind",b)
# es<-rbind(a,b)
#
# a<-afp[!grepl(",",afp$vdpv_emergence_group),]
# b<-afp[grepl(",",afp$vdpv_emergence_group),]
# b<-lapply(1:nrow(b),function(e){
#   pvt=strsplit(b$polio_virus_types[e],", ")[[1]]
#   pvt<-pvt[!grepl("VACCINE",pvt)]
#   nt=strsplit(b$vdpv_nt_changes[e],", ")[[1]]
#   nt<-nt[grepl("VDPV",nt)]
# merge(b[e,which(!names(b)%in%c("vdpv_emergence_group","vdpv_nt_changes","polio_virus_types"))],
#                                       cbind.data.frame(polio_virus_types=pvt,
#                                                        vdpv_emergence_group=strsplit(b$vdpv_emergence_group[e],", ")[[1]],
#                                                        vdpv_nt_changes=nt))})
# b<-do.call("rbind",b)
# afp<-rbind(a,b)
#
# afp$date<-as.Date(afp$case_date)
# es$date<-as.Date(es$collection_date)
# afp$type<-"AFP"
# es$type<-"ES"
# viruses<-bind_rows(afp%>%select(final_class,vdpv_emergence_group,date,iso3_code,vdpv_nt_changes,y,x,type,adm1_name,adm2_name),es%>%select(final_class,vdpv_emergence_group,date,iso3_code,vdpv_nt_changes,y,x,type,adm1_name,adm2_name))
# viruses$vdpv_emergence_group_name<-viruses$vdpv_emergence_group
# viruses$year<-year(viruses$date)
# viruses$half_year<-floor_date(viruses$date,"6 months")
# viruses<-viruses[viruses$final_class=="cVDPV2"&!is.na(viruses$vdpv_emergence_group)&viruses$date>=as.Date("2016-05-01"),]
# y<-summarise(group_by(summarise(group_by(as_tibble(viruses),vdpv_emergence_group_name),year=min(as.numeric(year))),year),new_emerg=n())
# z<-summarise(group_by(as_tibble(viruses),year),countries=length(unique(iso3_code)))

x<-unique(sia[,c("VaccineType","ParentActivityCode","ParentPlannedDateFrom","dose")])
x$date<-as.Date(x$ParentPlannedDateFrom,"%d/%m/%Y")
x$date[is.na(x$date)]<-as.Date(paste0("01/",x$ParentPlannedDateFrom[is.na(x$date)]),"%d/%b/%Y")
x<-x[x$date>=as.Date("2016-05-01"),]
x$year<-format(x$date,"%Y")
x$half_year<-floor_date(x$date,"6 months")
x<-summarise(group_by(x[x$half_year>as.Date("2016-01-01"),],half_year,VaccineType),n=sum(dose))

summarise(group_by(as_tibble(viruses_raw),vdpv_emergence_group_name,vacc),half_year=min(half_year))%>%filter(half_year==as.Date("2023-01-01"))

# exc<-c("NIE-KTS-1","MOZ-NPL-2", "RSS-WEQ-1")
exc <- c("")
y<-summarise(group_by(summarise(group_by(as_tibble(viruses_raw%>%filter(!vdpv_emergence_group_name%in%exc)),vdpv_emergence_group_name,vacc),half_year=min(half_year)),half_year,vacc),new_emerg=n())
y<-merge(cbind.data.frame(vacc="nOPV2",half_year=unique(viruses_raw$half_year[viruses_raw$half_year>=as.Date("2021-01-01")])),y,all=T)
y$y<-y$new_emerg
y$y[is.na(y$y)]<-0
y<-y[y$half_year>as.Date("2016-01-01"),]
p<-ggplot()+geom_col(data=x,aes(x=half_year+365/4,y=n/1e6,fill=VaccineType))+
  geom_point(data=y,aes(x=half_year+365/4,y=12*y,color=vacc))+
  geom_line(data=y,aes(x=half_year+365/4,y=12*y,color=vacc))+
  # geom_label(data=y,aes(x=half_year+365/4,y=12*(y+ifelse(half_year=="2023-01-01",ifelse(vacc=="nOPV2",-1,0.9),0.9)),label=y,color=vacc),size=2)+
  geom_label(data=y,aes(x=half_year+365/4,y=12*(y+0.9),label=y,color=vacc),size=2)+
  labs(x="",y="Estimated doses (millions)",linetype="",fill="Doses",color="Linked to",caption="*Data as of 6 Feb 2024.")+
  theme(legend.position=c(0.1,0.68))+scale_x_date(date_breaks="year",date_labels="%b %y")+
  scale_fill_manual(values=c("dark blue","orange","gray"))+scale_color_manual(values=c("blue","red"))+scale_linetype_manual(values=c(1,2))+
  scale_y_continuous(sec.axis=sec_axis(trans=~./12,name="New emergences (date of detection)"))
p
ggsave("figures/nopv2_new_emergence_plot_alt.png",p,width=6,height=4)
