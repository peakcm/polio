library(lubridate)
library(tidyverse)
# Vector of emergence groups linked to nOPV2
# n<-c("RDC-TAN-2","RDC-KOR-1","RDC-SKV-1","CAF-KEM-1","NIE-KBS-1","CAF-BNG-3","RDC-HKA-2","CAF-cVDPV2")
n<-c("RDC-TAN-2","RDC-KOR-1","RDC-SKV-1","CAF-KEM-1","NIE-KBS-1","CAF-BNG-3","RDC-HKA-2", "EGY-NOR-1")
# Countries
iso<-c("NGA","COD","CAF", "EGY")
ct<-c("NIGERIA","DEMOCRATIC REPUBLIC OF THE CONGO","CENTRAL AFRICAN REPUBLIC", "EGYPT")

path_to_polis_data<-"~/polio-immunity-mapping/results/"

upd<-read.csv(paste0("upd.csv"),sep=";")
names(upd)<-gsub(".","_",tolower(names(upd)),fixed=T)
upd$vdpv_emergence_group_name<-upd$emergence_group_s_
upd$admin1name<-upd$place_admin_1
upd$nt_changes<-as.character(upd$nt_changes)
upd$virus_date<-as.character(as.Date(upd$virus_date,"%d/%m/%Y"))
virus<-read.csv(paste0(path_to_polis_data,"linelist_virus.csv"))%>%filter(virus_type_name%in%c("aVDPV2","VDPV2","cVDPV2"),virus_date>=as.Date("2020-01-01"))
virus<-bind_rows(virus,upd%>%select(vdpv_emergence_group_name,nt_changes,virus_date,admin1name))
imm<-bind_rows(read.csv(paste0(path_to_polis_data,"immunity_trace.csv")),read.csv(paste0(path_to_polis_data,"immunity_trace_0.5.csv")))%>%filter(ISO_3_CODE%in%iso)
imm<-pivot_longer(data=imm,cols=starts_with("p2_"),names_prefix="p2_",names_to=c("y","m"),names_sep="_")
imm$date<-as.Date(paste0(imm$y,"-",imm$m,"-01"))
imm<-imm%>%filter(date>=as.Date("2020-01-01"),date<Sys.Date(),summ_level==1)

sia<-read.csv(paste0(path_to_polis_data,"sia_district_rows.csv"))%>%filter(vaccinetype%in%c("mOPV2","nOPV2","tOPV"),status%in%c("Done"),ADM0_NAME%in%ct,start_date>=as.Date("2020-01-01"))

origin<-virus%>%filter(vdpv_emergence_group_name%in%n)
fn<-function(nt,date,days=NA){
  sub<-cbind.data.frame(nt=as.numeric(nt),date=as.Date(date))
  sub<-sub[!is.na(sub$nt)&!is.na(sub$date),]
  if(!is.na(days)){
    sub$time<-as.numeric(sub$date-min(sub$date))
    sub<-sub[sub$time<days,]}
  return(quantile(unlist(foreach(i=1:nrow(sub))%do%(sub$date[i]-rgamma(10000, shape=sub$nt[i], rate = 1/35))),probs=c(0.025,0.5,0.975),na.rm=T))
}
library(foreach)
emerg<-foreach(E=unique(origin$vdpv_emergence_group_name),.combine="rbind.data.frame")%do%fn(nt=origin$nt_changes[origin$vdpv_emergence_group_name==E],date=origin$virus_date[origin$vdpv_emergence_group_name==E])
for(j in 1:3){emerg[,j]<-as.Date(emerg[,j],origin = "1970-01-01")}
names(emerg)<-paste0("date_emergence_",c("lwr","med","upr"))
emerg<-cbind.data.frame(vdpv_emergence_group_name=unique(origin$vdpv_emergence_group_name),emerg)
emerg<-merge(emerg,summarise(group_by(origin,vdpv_emergence_group_name),start=min(virus_date),admin1name=admin1name[which(virus_date==start)]))

p<-NULL
for(i in 1:nrow(emerg)){
  a<-emerg[i,]
  a$x<-a$date_emergence_med
  a$type<-"Emergence"
  b<-origin%>%filter(admin1name%in%emerg$admin1name[i],vdpv_emergence_group_name==emerg$vdpv_emergence_group_name[i])
  b$x<-as.Date(b$virus_date)
  b$type<-"Detection"

p[[i]]<-ggplot()+geom_line(data=imm%>%filter(ADM1_NAME==emerg$admin1name[i]),aes(x=date,y=value,linetype=coverage))+
  geom_segment(data=sia%>%filter(ADM1_NAME==emerg$admin1name[i]),aes(x=as.Date(start_date),xend=as.Date(start_date),y=100,yend=105))+
  geom_errorbarh(data=a,aes(xmin=date_emergence_lwr,xmax=date_emergence_upr,y=95,color=type))+
  geom_point(data=bind_rows(a,b),aes(x=x,y=95,color=type))+
  labs(x="",y="Immunity (6-36 months)",color="Type",linetype="Coverage",shape="Emergence",title=paste0(emerg$admin1name[i],", ",emerg$vdpv_emergence_group_name[i]))
}
library(ggpubr)
ggarrange(plotlist=p,common.legend=T)
ggsave("nopv2_emergence_timeseries.png",width=10,height=8)


load("data_local/adm1data.Rdata")
library(sf)
adm1data<-st_as_sf(adm1data[adm1data$WHO_REGION%in%c("AFRO"),])
adm1data<-adm1data[adm1data$ENDDATE>as.Date("2030-01-01"),]
adm1data<-adm1data[st_is_valid(adm1data),]
ref<-st_is_within_distance(adm1data[adm1data$ADM1_NAME%in%emerg$admin1name,],adm1data,100*1000)
names(ref)<-adm1data$ADM1_NAME[adm1data$ADM1_NAME%in%emerg$admin1name]

p<-NULL
for(i in 1:nrow(emerg)){
  a<-emerg[i,]
  a$x<-a$date_emergence_med
  a$type<-"Emergence"
  b<-origin%>%filter(vdpv_emergence_group_name==emerg$vdpv_emergence_group_name[i])
  b$x<-as.Date(b$virus_date)
  b$type<-"Detection"
  sel<-c(adm1data$ADM1_NAME[ref[[a$admin1name]]],b$admin1name)
     p[[i]]<-ggplot()+
    geom_text(data=sia%>%filter(ADM1_NAME%in%sel,vaccinetype=="nOPV2"),aes(x=as.Date(start_date),y=ADM1_NAME,label="|"))+
    geom_errorbarh(data=a,aes(xmin=date_emergence_lwr,xmax=date_emergence_upr,y=admin1name,color=type))+
    geom_point(data=bind_rows(a,b),aes(x=x,y=admin1name,color=type))+scale_x_date(limits=as.Date(c("2020-01-01","2023-08-01")))+
    labs(x="",y="",color="Type",linetype="Coverage",shape="Emergence",title=emerg$vdpv_emergence_group_name[i])
}

ggarrange(plotlist=p,common.legend=T)
ggsave("nopv2_emergence_carpet.png",width=10,height=8)

