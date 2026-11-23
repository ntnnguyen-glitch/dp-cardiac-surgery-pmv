# manuscript-level summaries: Table 1, absolute risks, clustered SE, counts
suppressMessages({library(data.table);library(DBI);library(sandwich)})
D<-as.data.table(readRDS("D.rds")); A<-as.data.table(readRDS("Dd.rds")); B<-readRDS("R.rds"); set.seed(2026)
D[,x:=as.integer(dp>=15)]; D[,`:=`(d28=as.integer(!is.na(death_h)&death_h<=672),y1=as.integer(!is.na(death_h)&death_h<=8760))]
M2<-"age+male+wt+emerg+cabg+comb+chf+copd+ckd+dm+af+mi+pvd+htn+liver+obese+ph_"
O<-list()
## Table 1 (observed values, not imputed)
md<-function(v) sprintf("%.1f (%.1f–%.1f)",median(v,na.rm=T),quantile(v,.25,na.rm=T),quantile(v,.75,na.rm=T)); np<-function(v) sprintf("%d (%.1f)",sum(v==1,na.rm=T),100*mean(v==1,na.rm=T))
pv<-function(p) if(is.na(p)) "" else if(p<0.001) "<0.001" else sprintf("%.3f",p)
cont<-c(age="Age, years",wt="Weight, kg",dp="Driving pressure, cmH₂O",pplat="Plateau pressure, cmH₂O",peep="PEEP, cmH₂O",vt="Tidal volume, mL",vtk="Tidal volume, mL/kg PBW",pf="PaO₂/FiO₂, mmHg",paco2="PaCO₂, mmHg",sev="Severity score",lact="Lactate, mmol/L",creat="Creatinine, mg/dL",hb="Haemoglobin, g/dL",wbc="White cell count, ×10⁹/L",plt="Platelets, ×10⁹/L",glu="Glucose, mg/dL",hr="Heart rate, beats/min",map="Mean arterial pressure, mmHg")
bin<-c(male="Male sex",emerg="Urgent or emergency admission",cabg_="Isolated CABG",valve_="Isolated valve surgery",comb_="CABG plus valve surgery",chf="Heart failure",copd="COPD",ckd="Chronic kidney disease",dm="Diabetes mellitus",af="Atrial fibrillation",mi="Myocardial infarction",pvd="Peripheral vascular disease",htn="Hypertension",liver="Liver disease",obese="Obesity",ph_="Pulmonary hypertension",vaso="Vasopressor, 0–6 h",ino="Inotrope, 0–6 h",pmv="Ventilation >24 h",reint="Reintubation within 7 days",d28="28-day mortality",hd="In-hospital mortality")
D[,`:=`(cabg_=as.integer(ptype=="cabg"),valve_=as.integer(ptype=="valve"),comb_=as.integer(ptype=="comb"))]
T1<-list(); for(k in c("MIMIC-IV","eICU-CRD")){d<-D[db==k]; r<-list(data.table(Characteristic="n",Overall=format(nrow(d),big.mark=","),`DP <15 cmH₂O`=format(sum(d$x==0),big.mark=","),`DP ≥15 cmH₂O`=format(sum(d$x==1),big.mark=","),`P value`=""))
 for(v in names(cont)){z<-d[[v]]; if(all(is.na(z))||mean(is.na(z))>0.5&&!v%in%c("vtk")) next; if(k=="eICU-CRD"&&v=="lact") next; lab<-cont[v]; if(v=="sev") lab<-if(k=="MIMIC-IV") "SOFA score, day 1 (non-respiratory, non-neurological)" else "APACHE IV score"
  r[[length(r)+1]]<-data.table(Characteristic=lab,Overall=md(z),`DP <15 cmH₂O`=md(z[d$x==0]),`DP ≥15 cmH₂O`=md(z[d$x==1]),`P value`=pv(tryCatch(wilcox.test(z~d$x)$p.value,error=function(e) NA)))}
 for(v in names(bin)){z<-d[[v]]; if(all(is.na(z))) next; r[[length(r)+1]]<-data.table(Characteristic=bin[v],Overall=np(z),`DP <15 cmH₂O`=np(z[d$x==0]),`DP ≥15 cmH₂O`=np(z[d$x==1]),`P value`=pv(tryCatch(suppressWarnings(chisq.test(table(z,d$x))$p.value),error=function(e) NA)))}
 T1[[k]]<-rbindlist(r); O$miss[[k]]<-round(100*sapply(d[,.(pf,lact,creat,hb,wbc,plt,glu,hr,map,sev,vt,vtk,wt)],function(z) mean(is.na(z))),1)}
O$T1<-T1
## counts
O$cnt<-D[,.(n=.N,hi=sum(x),pmv0=sum(pmv[x==0]),n0=sum(x==0),pmv1=sum(pmv[x==1]),n1=sum(x==1),re0=sum(reint[x==0],na.rm=T),re1=sum(reint[x==1],na.rm=T),d280=sum(d28[x==0]),d281=sum(d28[x==1]),hd0=sum(hd[x==0]),hd1=sum(hd[x==1]),y10=sum(y1[x==0]),y11=sum(y1[x==1]),age=median(age),male=mean(male),dpm=median(dp),dpq1=quantile(dp,.25),dpq3=quantile(dp,.75),mvh=as.numeric(median(mvh,na.rm=T)),mvq1=as.numeric(quantile(mvh,.25,na.rm=T)),mvq3=as.numeric(quantile(mvh,.75,na.rm=T)),ndp=as.numeric(median(ndp))),by=db]
## standardised risks and hospital-clustered estimate are taken from the multiply imputed analysis (02_analysis.R)
for(k in c("MIMIC-IV","eICU-CRD")) O$ar[[k]]<-B[[k]]$ar; O$cl<-c(B[["eICU-CRD"]]$cl,B[["eICU-CRD"]]$nhosp)
## DP distribution & correlation with Vt, PEEP
O$cor<-D[,.(r_vt=cor(dp,vt,use="complete.obs"),r_pf=cor(dp,pf,use="complete.obs"),prop_hi=mean(x)),by=db]
saveRDS(O,"M.rds"); print(O$cnt); print(O$ar); print(round(O$cl,4)); print(O$miss); print(O$cor); print(T1[[1]][1:6])
