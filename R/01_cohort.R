# 01 cohort: cardiac surgery, ventilated on ICU arrival, driving pressure 0-6 h
suppressMessages({library(DBI);library(data.table)})
con<-dbConnect(duckdb::duckdb(shared_home=FALSE),":memory:"); invisible(dbExecute(con,sprintf("ATTACH '%s' AS w (READ_ONLY)",Sys.getenv("DUCKDB","w.duckdb")))); invisible(dbExecute(con,"PRAGMA threads=8"))
Q<-function(s) as.data.table(dbGetQuery(con,s)); H<-function(a,b) as.numeric(difftime(a,b,units="hours")); fl<-list()
MM<-paste0(Sys.getenv("MIMIC_DIR"),"/"); EE<-paste0(Sys.getenv("EICU_DIR"),"/")
rcsv<-function(p,sel="*",wh="") Q(sprintf("select %s from read_csv('%s',all_varchar=true,ignore_errors=true) %s",sel,p,wh))
# ================= MIMIC-IV =================
m<-Q("select * from w.st"); fl$m0<-nrow(m)
pr<-rcsv(paste0(MM,"hosp/procedures_icd.csv.gz"),"hadm_id,icd_code,icd_version,chartdate"); pr[,hadm_id:=as.integer(hadm_id)]; pr[,chartdate:=as.Date(chartdate)]
pr[,cabg:=(icd_version=="9"&grepl("^361",icd_code))|(icd_version=="10"&grepl("^021[0-3]0",icd_code))]
pr[,valve:=(icd_version=="9"&grepl("^35(1|2)",icd_code))|(icd_version=="10"&grepl("^02(R|Q)[FGHJ]0",icd_code))]
pr[,excl:=(icd_version=="9"&grepl("^(3965|375)",icd_code))|(icd_version=="10"&grepl("^(5A15|02YA|02HA0)",icd_code))]
cs<-pr[cabg|valve,.(cabg=any(cabg),valve=any(valve),sdate=min(chartdate)),by=hadm_id]; m<-merge(m,cs,by="hadm_id"); fl$m_cs<-nrow(m)
m<-m[as.numeric(as.Date(intime)-sdate)%in%c(0,1)]; fl$m_postop<-nrow(m)
m<-m[!hadm_id%in%pr[excl==TRUE]$hadm_id]; fl$m_ecmo<-nrow(m)
pe<-Q("select stay_id,starttime,endtime from w.pe where itemid=225792"); pe<-merge(pe[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); pe[,`:=`(s=H(starttime,intime),e=H(endtime,intime))]; setorder(pe,stay_id,s)
a0<-pe[s<=6&e>s][,.SD[1],by=stay_id]; m<-m[stay_id%in%a0$stay_id]; fl$m_vent<-nrow(m)
m[,ve:=a0$e[match(stay_id,a0$stay_id)]]
c2<-Q("select stay_id,itemid,charttime,valuenum v from w.ce2 where itemid in (220339,224696,224685,224684)"); c2<-merge(c2[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); c2[,h:=H(charttime,intime)]
w<-dcast(c2[h>=0&h<=24],stay_id+charttime+h~itemid,value.var="v",fun.aggregate=function(x)x[1]); setnames(w,c("220339","224696","224685","224684"),c("peep","pplat","vt","vts"),skip_absent=TRUE)
w<-w[!is.na(peep)&!is.na(pplat)&peep>=0&peep<=25&pplat>=5&pplat<=50&pplat>peep]; w[,dp:=pplat-peep]; w[,vtx:=fifelse(!is.na(vt)&vt>=100&vt<=1200,vt,fifelse(!is.na(vts)&vts>=100&vts<=1200,vts,NA_real_))]
w6<-w[h<=6,.(dp=median(dp),peep=median(peep),pplat=median(pplat),vt=median(vtx,na.rm=TRUE),ndp=.N),by=stay_id]; w24<-w[,.(dp24=median(dp)),by=stay_id]; wf<-w[order(h)][,.(dpfirst=dp[1]),by=stay_id]
m<-merge(m,w6,by="stay_id"); m[,`:=`(dp24=w24$dp24[match(stay_id,w24$stay_id)],dpfirst=wf$dpfirst[match(stay_id,wf$stay_id)])]; fl$m_dp<-nrow(m)
dth<-as.POSIXct(m$deathtime,tz="UTC"); dd<-as.POSIXct(ifelse(is.na(m$dod),NA,paste(as.character(m$dod),"23:59:00")),tz="UTC"); dth[is.na(dth)]<-dd[is.na(dth)]; m[,death_h:=pmax(H(dth,intime),0)]
m<-m[is.na(death_h)|death_h>24]; fl$m_d24<-nrow(m)
ad<-rcsv(paste0(MM,"hosp/admissions.csv.gz"),"hadm_id,admission_type"); ad[,hadm_id:=as.integer(hadm_id)]; m[,emerg:=as.integer(ad$admission_type[match(hadm_id,ad$hadm_id)]%in%c("EW EMER.","DIRECT EMER.","URGENT"))]
om<-rcsv(paste0(MM,"hosp/omr.csv.gz"),"subject_id,result_name,result_value","where result_name like 'Height%'"); om[,ht:=suppressWarnings(as.numeric(result_value))*2.54]; om<-om[ht>=130&ht<=210,.(ht=median(ht)),by=subject_id]; om[,subject_id:=as.integer(subject_id)]; m[,ht:=om$ht[match(subject_id,om$subject_id)]]
dx<-Q("select hadm_id,icd_code from w.dx"); fx<-function(p) as.integer(m$hadm_id%in%dx[grepl(p,icd_code)]$hadm_id)
m[,`:=`(chf=fx("^(428|I50)"),copd=fx("^(49[0-2]|496|J4[0-4])"),ckd=fx("^(585|N18)"),dm=fx("^(250|E1[0-4])"),af=fx("^(42731|42732|I48)"),mi=fx("^(410|I21|I22)"),pvd=fx("^(4439|440|I70|I739)"),htn=fx("^(40[1-5]|I1[0-5])"),liver=fx("^(571|K70|K7[2-6])"),obese=fx("^(2780|E66)"),pna=fx("^(48[0-6]|J1[2-8])"),ards=fx("^(51882|J80)"),ph_=fx("^(416|I27)"),male=as.integer(gender=="M"))]
L<-Q("select stay_id,specimen_id,itemid,charttime,v,value from w.lab"); L<-L[stay_id%in%m$stay_id]; L<-merge(L,m[,.(stay_id,intime)],by="stay_id"); L[,h:=H(charttime,intime)]
gv<-function(it,lo,hi,fun="near",a=-6,z=6){x<-L[itemid%in%it&v>=lo&v<=hi&h>=a&h<=z]; x<-if(fun=="near") x[order(abs(h))][,.(v=v[1]),by=stay_id] else x[,.(v=if(fun=="max") max(v) else min(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(lact=gv(50813,.1,30),bun=gv(51006,1,250),creat=gv(50912,.1,25),hb=gv(51222,3,20),wbc=gv(51301,.1,200),plt=gv(51265,1,1500),glu=gv(50931,20,1500),na=gv(50983,100,180),bili=gv(50885,.1,60,"max",-6,24),pltmin=gv(51265,1,1500,"min",0,24),crmax=gv(50912,.1,25,"max",0,24))]
sp<-L[itemid==52033,.(specimen_id,st=value)]; g<-dcast(L[itemid%in%c(50818,50820,50821,50816)&h>=0&h<=6],stay_id+specimen_id+h~itemid,value.var="v",fun.aggregate=function(x)x[1]); setnames(g,c("50818","50820","50821","50816"),c("paco2","ph","pao2","fl"),skip_absent=TRUE)
g[,st:=sp$st[match(specimen_id,sp$specimen_id)]]; g<-g[(is.na(st)|st=="ART.")&!is.na(pao2)&pao2>=20&pao2<=700]
ce<-Q("select stay_id,itemid,charttime,v from w.ce"); ce<-merge(ce[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); ce[,h:=H(charttime,intime)]
fi<-ce[itemid==223835&v>=21&v<=100]; pc2<-merge(g[,.(stay_id,specimen_id,ha=h)],fi[,.(stay_id,hf=h,f=v)],by="stay_id",allow.cartesian=TRUE)[abs(ha-hf)<=2][order(abs(ha-hf))][,.SD[1],by=specimen_id]
g[,fc:=pc2$f[match(specimen_id,pc2$specimen_id)]]; g[,f:=fifelse(!is.na(fl)&fl>=21&fl<=100,fl,fc)]; pp<-g[!is.na(f)][order(h)][,.SD[1],by=stay_id]
m<-merge(m,pp[,.(stay_id,pao2,paco2,ph,fio2=f)],by="stay_id",all.x=TRUE); m[,pf:=pao2/(fio2/100)]
cm<-function(its,lo,hi,fun=mean,a=0,z=6){x<-ce[itemid%in%its&v>=lo&v<=hi&h>=a&h<=z,.(v=fun(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(hr=cm(220045,20,250),map=cm(c(220052,220181,225312),20,200),mapmin=cm(c(220052,220181,225312),20,200,min,0,24),rr=cm(c(220210,224690),3,70))]
w8<-ce[itemid%in%c(226512,224639)&v>=30&v<=250][order(h)][,.SD[1],by=stay_id]; m[,wt:=w8$v[match(stay_id,w8$stay_id)]]
ie<-Q(sprintf("select stay_id,itemid,starttime from read_csv('%sicu/inputevents.csv.gz',all_varchar=true,ignore_errors=true) where itemid in ('221906','221289','221662','222315','221749','221653','221986') and stay_id in (%s)",MM,paste(m$stay_id,collapse=",")))
ie[,stay_id:=as.integer(stay_id)]; ie<-merge(ie,m[,.(stay_id,intime)],by="stay_id"); ie[,h:=H(as.POSIXct(starttime,tz="UTC"),intime)]; ie<-ie[h>=-6&h<=6]
m[,`:=`(vaso=as.integer(stay_id%in%ie[itemid%in%c("221906","222315","221749","221662")]$stay_id),ino=as.integer(stay_id%in%ie[itemid%in%c("221289","221653","221986")]$stay_id))]
sc<-function(x,br) fifelse(is.na(x),0L,as.integer(findInterval(x,br)))
m[,sev:=as.numeric(fifelse(is.na(pltmin),0L,4L-sc(pltmin,c(20,50,100,150)))+sc(bili,c(1.2,2,6,12))+fifelse(vaso==1|ino==1,3L,fifelse(!is.na(mapmin)&mapmin<70,1L,0L))+sc(crmax,c(1.2,2,3.5,5)))]
ri<-merge(pe,m[,.(stay_id,ve)],by="stay_id")[s>ve+1&s<=ve+168]; r1<-ri[,.(t=min(s)),by=stay_id]; m[,`:=`(reint=as.integer(stay_id%in%ri$stay_id),reint_h=r1$t[match(stay_id,r1$stay_id)])]
tv<-merge(pe,m[,.(stay_id)],by="stay_id")[s<672,.(vh=sum(pmin(e,672)-pmax(s,0))),by=stay_id]; m[,venth28:=tv$vh[match(stay_id,tv$stay_id)]]
m[,`:=`(mvh=ve,los=H(outtime,intime),fu_h=pmin(fifelse(is.na(death_h),365*24,death_h),365*24),hd=hospital_expire_flag,ptype=fifelse(cabg&valve,"comb",fifelse(cabg,"cabg","valve")),db="MIMIC-IV",id=stay_id)]
# ================= eICU-CRD =================
e<-Q("select cast(patientunitstayid as bigint) s,uniquepid,gender,age,try_cast(hospitaldischargeoffset as int) hdo,hospitaldischargestatus hds,try_cast(unitdischargeoffset as int) udo,try_cast(unitvisitnumber as int) uvn,hospitaldischargeyear yr,apacheadmissiondx adx,hospitalid hosp,try_cast(admissionweight as double) wt,try_cast(admissionheight as double) ht from w.ept")
e<-e[uvn==1][order(yr,s)][,.SD[1],by=uniquepid]; e[,age:=as.numeric(ifelse(age=="> 89","90",age))]; e<-e[!is.na(age)&age>=18&hds%in%c("Alive","Expired")&!is.na(hdo)]
e<-e[grepl("^CABG|valve replacement|valve repair|Valve replacement|Tricuspid valve surgery|double valve",adx)]; fl$e_cs<-nrow(e)
r2<-Q("select s,o,l,v from w.erc2 where l in ('PEEP','Plateau Pressure','Exhaled TV (machine)','Tidal Volume Observed (VT)','Tidal Volume (set)','FiO2','FIO2 (%)','Vent Rate') and o between -60 and 1800"); r2<-r2[s%in%e$s]; r2[,v:=suppressWarnings(as.numeric(v))]; r2<-r2[!is.na(v)]
ven<-unique(r2[l%in%c("PEEP","Vent Rate","Tidal Volume (set)")&o>=0&o<=60]$s); e<-e[s%in%ven]; fl$e_vent<-nrow(e)
r2[l%in%c("Exhaled TV (machine)","Tidal Volume Observed (VT)","Tidal Volume (set)"),l:="vt"]
w<-dcast(r2[l%in%c("PEEP","Plateau Pressure","vt")&o>=0&o<=1440],s+o~l,value.var="v",fun.aggregate=function(x)x[1]); setnames(w,c("PEEP","Plateau Pressure"),c("peep","pplat"))
w<-w[!is.na(peep)&!is.na(pplat)&peep>=0&peep<=25&pplat>=5&pplat<=50&pplat>peep]; w[,dp:=pplat-peep]; w[,vtx:=fifelse(!is.na(vt)&vt>=100&vt<=1200,vt,NA_real_)]
w6<-w[o<=360,.(dp=median(dp),peep=median(peep),pplat=median(pplat),vt=median(vtx,na.rm=TRUE),ndp=.N),by=s]; w24<-w[,.(dp24=median(dp)),by=s]; wf<-w[order(o)][,.(dpfirst=dp[1]),by=s]
e<-merge(e,w6,by="s"); e[,`:=`(dp24=w24$dp24[match(s,w24$s)],dpfirst=wf$dpfirst[match(s,wf$s)])]; fl$e_dp<-nrow(e)
e[,death_h:=fifelse(hds=="Expired",pmax(hdo/60,0),NA_real_)]; e<-e[is.na(death_h)|death_h>24]; fl$e_d24<-nrow(e)
v24<-unique(r2[l%in%c("PEEP","Vent Rate","vt","Plateau Pressure")&o>1440&o<=1800]$s); e[,pmv:=as.integer(s%in%v24)]
pv<-rcsv(paste0(EE,"apachePredVar.csv.gz"),"patientunitstayid,electivesurgery"); pv[,s:=as.numeric(patientunitstayid)]; e[,emerg:=as.integer(pv$electivesurgery[match(s,pv$s)]%in%c("0"))]
edx<-Q("select s,d from w.edx"); eph<-Q("select s,d from w.eph"); fe<-function(p) as.integer(e$s%in%c(edx[grepl(p,d)]$s,eph[grepl(p,d)]$s)); fd<-function(p) as.integer(e$s%in%edx[grepl(p,d)]$s)
e[,`:=`(chf=fe("congestive heart failure|chf"),copd=fe("copd|chronic obstructive|emphysema|chronic bronchitis"),ckd=fe("chronic kidney|chronic renal|renal insufficiency"),dm=fe("diabetes"),af=fe("atrial fibrillation"),mi=fe("myocardial infarction"),pvd=fe("peripheral vascular|peripheral arter|claudication"),htn=fe("hypertension"),liver=fe("cirrhosis|chronic liver"),obese=fe("obesity"),pna=fd("pneumonia"),ards=fd("ards|acute respiratory distress"),ph_=fe("pulmonary hypertension"),male=as.integer(gender=="Male"))]
EL<-Q("select s,o,n,v from w.elab"); EL<-EL[s%in%e$s&!is.na(v)]
ge<-function(nm,lo,hi,fun="near",a=-720,z=720){x<-EL[n==nm&v>=lo&v<=hi&o>=a&o<=z]; x<-if(fun=="near") x[order(abs(o))][,.(v=v[1]),by=s] else x[,.(v=if(fun=="max") max(v) else min(v)),by=s]; x$v[match(e$s,x$s)]}
e[,`:=`(lact=ge("lactate",.1,30),bun=ge("BUN",1,250),creat=ge("creatinine",.1,25),hb=ge("Hgb",3,20),wbc=ge("WBC x 1000",.1,200),plt=ge("platelets x 1000",1,1500),glu=ge("glucose",20,1500),na=ge("sodium",100,180))]
po<-EL[n=="paO2"&v>=20&v<=700&o>=0&o<=360][order(o)][,.SD[1],by=s]; fi<-rbind(r2[l%in%c("FiO2","FIO2 (%)")&v>=21&v<=100,.(s,of=o,f=v)],EL[n=="FiO2"&v>=21&v<=100,.(s,of=o,f=v)])
pp<-merge(po[,.(s,o,pao2=v)],fi,by="s",allow.cartesian=TRUE)[of<=o+15&o-of<=240][order(s,o-of)][,.SD[1],by=s]; e[,`:=`(pao2=pp$pao2[match(s,pp$s)],fio2=pp$f[match(s,pp$s)])]; e[,pf:=pao2/(fio2/100)]
pc<-EL[n=="paCO2"&v>=10&v<=150&o>=0&o<=360][order(o)][,.SD[1],by=s]; e[,paco2:=pc$v[match(s,pc$s)]]
aps<-Q("select cast(patientunitstayid as bigint) s,try_cast(heartrate as double) hr,try_cast(meanbp as double) map,try_cast(respiratoryrate as double) rr from w.eaps"); a<-aps[match(e$s,aps$s)]
e[,`:=`(hr=fifelse(a$hr>0,a$hr,NA_real_),map=fifelse(a$map>0,a$map,NA_real_),rr=fifelse(a$rr>0,a$rr,NA_real_))]
ap<-Q("select * from w.eapr"); e[,sev:=as.numeric(ap$apache[match(s,ap$s)])]; e[sev<0,sev:=NA]
inf<-Q(sprintf("select cast(patientunitstayid as bigint) s,try_cast(infusionoffset as int) o,lower(drugname) d from read_csv('%sinfusionDrug.csv.gz',all_varchar=true,ignore_errors=true) where cast(patientunitstayid as bigint) in (%s)",EE,paste(e$s,collapse=","))); inf<-inf[o>=-360&o<=360]
e[,`:=`(vaso=as.integer(s%in%inf[grepl("norepi|levophed|vasopressin|phenyleph|neo-syn|dopamine",d)]$s),ino=as.integer(s%in%inf[grepl("epineph|dobut|milrin|primacor",d)&!grepl("norepi",d)]$s))]
m[,hosp:=NA_character_]; e[,`:=`(los=udo/60,fu_h=fifelse(is.na(death_h),pmax(hdo/60,0.01),death_h),hd=as.integer(hds=="Expired"),ptype=fifelse(grepl("^CABG",adx)&grepl("valve",adx),"comb",fifelse(grepl("^CABG",adx),"cabg","valve")),db="eICU-CRD",id=s,mvh=NA_real_,reint=NA_integer_,reint_h=NA_real_,venth28=NA_real_)]
m[,pmv:=as.integer(mvh>24)]
cols<-c("id","db","hosp","age","male","wt","ht","ptype","emerg","dp","dp24","dpfirst","ndp","peep","pplat","vt","pao2","paco2","fio2","pf","lact","bun","creat","hb","wbc","plt","glu","na","hr","map","rr","sev","vaso","ino","chf","copd","ckd","dm","af","mi","pvd","htn","liver","obese","pna","ards","ph_","pmv","mvh","reint","reint_h","venth28","death_h","fu_h","hd","los")
D<-rbind(m[,..cols],e[,..cols]); D[,pbw:=fifelse(male==1,50,45.5)+0.91*(ht-152.4)]; D[!is.na(ht)&(ht<130|ht>210),pbw:=NA]; D[,vtk:=vt/pbw]
saveRDS(D,"D.rds"); saveRDS(fl,"fl.rds"); print(unlist(fl))
print(D[,.(n=.N,dp=median(dp),hi=mean(dp>=15),pmv=mean(pmv),d28=mean(!is.na(death_h)&death_h<=672),hd=mean(hd),reint=mean(reint)),by=db])
print(D[,lapply(.SD,function(x) round(mean(is.na(x)),2)),by=db,.SDcols=c("pf","lact","creat","hb","hr","map","sev","vt","vtk","emerg","wt")])
