# 03a ventilator mode, outcome-definition and selection data (MIMIC-IV)
suppressMessages({library(data.table);library(DBI);library(mice);library(logistf);library(WeightIt)})
src<-readLines("../R/01_cohort.R"); eval(parse(text=src[1:22]))   # MIMIC up to DP extraction (m = ventilated on arrival, w = DP measurements)
V<-copy(m); V[,hasdp:=as.integer(stay_id%in%w[h<=6]$stay_id)]
dth<-as.POSIXct(V$deathtime,tz="UTC"); dd<-as.POSIXct(ifelse(is.na(V$dod),NA,paste(as.character(V$dod),"23:59:00")),tz="UTC"); dth[is.na(dth)]<-dd[is.na(dth)]; V[,death_h:=pmax(H(dth,intime),0)]
V<-V[is.na(death_h)|death_h>24]; V[,`:=`(pmv=as.integer(ve>24),male=as.integer(gender=="M"),ptype=fifelse(cabg&valve,"comb",fifelse(cabg,"cabg","valve")),d28=as.integer(!is.na(death_h)&death_h<=672))]
ad<-rcsv(paste0(MM,"hosp/admissions.csv.gz"),"hadm_id,admission_type"); ad[,hadm_id:=as.integer(hadm_id)]; V[,emerg:=as.integer(ad$admission_type[match(hadm_id,ad$hadm_id)]%in%c("EW EMER.","DIRECT EMER.","URGENT"))]
O<-list(); O$sel<-V[,.(n=.N,age=median(age),male=mean(male),cabg=mean(ptype=="cabg"),valve=mean(ptype=="valve"),comb=mean(ptype=="comb"),emerg=mean(emerg),pmv=mean(pmv),mvh=median(ve),d28=mean(d28)),keyby=hasdp]
sp<-glm(hasdp~age+male+factor(ptype)+emerg,binomial,V); V[,psel:=fitted(sp)]
## mode at time of DP measurement (MIMIC)
md<-Q("select stay_id,charttime,value from w.ce2 where itemid=223849"); md<-merge(md[stay_id%in%w$stay_id],m[,.(stay_id,intime)],by="stay_id"); md[,hm:=H(charttime,intime)]
ctrl<-c("CMV/ASSIST/AutoFlow","CMV/ASSIST","CMV/AutoFlow","PRVC/AC","PCV+Assist","VOL/AC","CMV","APV (cmv)","PCV+","PRVC/SIMV")
x<-merge(w[h<=6,.(stay_id,h,dp,peep)],md[,.(stay_id,hm,mode=value)],by="stay_id",allow.cartesian=TRUE)[abs(h-hm)<=1][order(abs(h-hm))][,.SD[1],by=.(stay_id,h)]
O$modes<-x[,.N,by=mode][order(-N)][1:12]; O$modecov<-c(withmode=nrow(x),all=nrow(w[h<=6]))
x[,ctl:=mode%in%ctrl]; cm<-x[ctl==TRUE,.(dpc=median(dp)),by=stay_id]; spn<-unique(x[ctl==FALSE]$stay_id)
saveRDS(list(V=V,cm=cm,spn=spn),"r1_mimic.rds")
## harmonised (eICU-like) outcome in MIMIC: ventilator settings charted 25-30 h
c3<-Q("select stay_id,itemid,charttime from w.ce2 where itemid in (220339,224684,224685,224696,223849)"); c3<-merge(c3[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); c3[,h:=H(charttime,intime)]
md2<-Q("select stay_id,charttime,value from w.ce2 where itemid=223849"); c3<-merge(c3,md2,by=c("stay_id","charttime"),all.x=TRUE); h25<-unique(c3[h>25&h<=30&(itemid!=223849|!value%in%c("Standby"))]$stay_id)
saveRDS(h25,"r1_h25.rds"); saveRDS(O,"r1_O.rds"); cat("mimic part ok\n")
