# 02 statistics: driving pressure 0-6 h and prolonged mechanical ventilation after cardiac surgery
# Missing data: covariates with >20% missing values in a database are not used; remaining missing values are
# handled by multiple imputation by chained equations (20 datasets); estimates are pooled with Rubin's rules.
suppressMessages({library(data.table);library(survival);library(splines);library(WeightIt);library(MatchIt);library(mice);library(pROC);library(EValue);library(cmprsk);library(sandwich);library(logistf)})
RAW<-as.data.table(readRDS("D.rds")); set.seed(2026); MIMP<-20
RAW[,`:=`(x=as.integer(dp>=15),d28=as.integer(!is.na(death_h)&death_h<=672),y1=as.integer(!is.na(death_h)&death_h<=8760),lpf=log(pmin(pmax(pf,30),700)),cabg=as.integer(ptype=="cabg"),comb=as.integer(ptype=="comb"))]
RAW[,dpc:=factor(fifelse(dp<10,"lt10",fifelse(dp<15,"10_14","ge15")),levels=c("lt10","10_14","ge15"))]
RAW[db=="MIMIC-IV",vfd:=fifelse(d28==1,0,pmax(28-fifelse(is.na(venth28),0,venth28)/24,0))]
M2<-"age+male+wt+emerg+cabg+comb+chf+copd+ckd+dm+af+mi+pvd+htn+liver+obese+ph_"
cand<-c("sev","vaso","ino","creat","hb","wbc","plt","glu","lact","hr","map","peep","lpf")
miss<-RAW[,lapply(.SD,function(z) mean(is.na(z))),by=db,.SDcols=c("wt",cand)]
M3v<-function(dbx){m<-miss[db==dbx]; c(cand[sapply(cand,function(v) m[[v]]<=0.20)])}
M3f<-function(dbx) paste(M2,"+",paste(sub("^sev$","sevz",M3v(dbx)),collapse="+"))
pr<-function(est,se,df=Inf){q<-qnorm(.975); c(exp(est),exp(est-q*se),exp(est+q*se),2*pnorm(-abs(est/se)))}
rub<-function(b,v){m<-length(b); qb<-mean(b); W<-mean(v); Bv<-var(b); Tt<-W+(1+1/m)*Bv; c(qb,sqrt(Tt))}
pl<-function(fits,tr){b<-sapply(fits,function(f) coef(f)[tr]); v<-sapply(fits,function(f) vcov(f)[tr,tr]); r<-rub(b,v); pr(r[1],r[2])}
plr<-function(lst){b<-sapply(lst,`[`,1); v<-sapply(lst,`[`,2)^2; r<-rub(b,v); pr(r[1],r[2])}
gl<-function(f,d,w=NULL) suppressWarnings(glm(as.formula(f),if(is.null(w)) binomial else quasibinomial,d,weights=w))
rob<-function(fit,tr="x",cl=NULL){V<-if(is.null(cl)) vcovHC(fit,type="HC0") else vcovCL(fit,cluster=cl); c(coef(fit)[tr],sqrt(V[tr,tr]))}
## imputation per database
IMP<-list(); for(dbx in c("MIMIC-IV","eICU-CRD")){d<-copy(RAW[db==dbx]); vv<-unique(c("pmv","x","dp",strsplit(gsub(" ","",M2),"\\+")[[1]],M3v(dbx)))
 im<-mice(d[,..vv],m=MIMP,method="pmm",printFlag=FALSE,seed=2026); IMP[[dbx]]<-lapply(1:MIMP,function(i){z<-copy(d); cc<-as.data.table(complete(im,i)); for(v in vv) set(z,j=v,value=cc[[v]]); z[,sevz:=(sev-mean(sev))/sd(sev)]; z})}
saveRDS(IMP,"Imp.rds")
## early-stop gate: primary analysis in MIMIC-IV
G<-pl(lapply(IMP[["MIMIC-IV"]],function(d) gl(paste("pmv~x+",M3f("MIMIC-IV")),d)),"x"); cat("GATE MIMIC",round(G,4),"\n"); if(G[4]>=0.05){cat("STOP: primary NS\n"); quit(status=3)}
an<-function(dbx){L<-IMP[[dbx]]; M3<-M3f(dbx); mi<-dbx=="MIMIC-IV"; o<-list(); d0<-RAW[db==dbx]; o$n<-nrow(d0); o$t1<-d0[,.(n=.N,ev=sum(pmv)),keyby=x]; o$M3<-M3
 F<-function(f) lapply(L,function(d) gl(f,d))
 o$crude<-pl(F("pmv~x"),"x"); o$m2<-pl(F(paste("pmv~x+",M2)),"x"); o$m3<-pl(F(paste("pmv~x+",M3)),"x"); o$cont<-pl(F(paste("pmv~dp+",M3)),"dp")
 fc<-F(paste("pmv~dpc+",M3)); o$cat<-rbind(pl(fc,"dpc10_14"),pl(fc,"dpcge15")); o$trend<-pl(F(paste("pmv~as.numeric(dpc)+",M3)),"as.numeric(dpc)")
 kn<-quantile(d0$dp,c(.1,.5,.9)); bk<-range(d0$dp); fr<-lapply(L,function(d) gl(sprintf("pmv~ns(dp,knots=c(%s),Boundary.knots=c(%s))+%s",paste(kn,collapse=","),paste(bk,collapse=","),M3),d)); f0<-F(paste("pmv~",M3)); fl<-F(paste("pmv~dp+",M3))
 o$rcs<-c(D1(as.mira(fr),as.mira(f0))$result[4],D1(as.mira(fr),as.mira(fl))$result[4])
 f1<-F(paste("pmv~x+",M3)); o$inc_lr<-D1(as.mira(f1),as.mira(f0))$result[4]
 au<-t(mapply(function(a,b,cc,d){r0<-roc(d$pmv,fitted(a),quiet=T); r1<-roc(d$pmv,fitted(b),quiet=T); r2<-roc(d$pmv,fitted(cc),quiet=T); c(auc(r0),auc(r1),roc.test(r0,r1)$p.value,auc(r2),roc.test(r0,r2)$p.value)},f0,f1,fl,L))
 o$inc<-c(o$inc_lr,apply(au,2,median)[1:3],apply(au,2,median)[4:5])
 ## propensity score methods within each imputation, pooled
 ps1<-lapply(L,function(d){fps<-as.formula(paste("x~",M3)); d<-copy(d); w1<-weightit(fps,d,method="glm",estimand="ATE"); d$ipw<-pmin(w1$weights,quantile(w1$weights,.99)); d$oww<-weightit(fps,d,method="glm",estimand="ATO")$weights
  a<-rob(gl("pmv~x",d,d$ipw)); b<-rob(gl("pmv~x",d,d$oww)); cdr<-rob(gl(paste("pmv~x+",M3),d,d$ipw))
  mt<-matchit(fps,d,method="nearest",caliper=.2,std.caliper=TRUE); md<-match.data(mt); e<-rob(gl("pmv~x",md),cl=~subclass)
  list(a=a,b=b,c=cdr,e=e,np=sum(md$x),smd=c(max(abs(summary(mt)$sum.all[-1,3])),max(abs(summary(mt)$sum.matched[-1,3]))))})
 o$iptw<-plr(lapply(ps1,`[[`,"a")); o$ow<-plr(lapply(ps1,`[[`,"b")); o$dr<-plr(lapply(ps1,`[[`,"c")); o$psm<-c(plr(lapply(ps1,`[[`,"e")),round(median(sapply(ps1,`[[`,"np"))))
 o$smd<-rowMeans(sapply(ps1,`[[`,"smd")); o$ev<-evalues.OR(o$m3[1],o$m3[2],o$m3[3],rare=mean(d0$pmv)<.15)[2,1]
 ## standardised risks (g-computation; 100 bootstrap samples per imputed dataset)
 g<-function(dd){m<-gl(paste("pmv~x+",M3),dd); p0<-mean(predict(m,transform(dd,x=0),type="response")); p1<-mean(predict(m,transform(dd,x=1),type="response")); c(p0,p1,p1-p0)}
 est<-rowMeans(sapply(L,g)); bs<-do.call(rbind,lapply(L,function(dd) t(replicate(100,g(dd[sample(.N,replace=TRUE)]))))); o$ar<-rbind(est=est,lo=apply(bs,2,quantile,.025),hi=apply(bs,2,quantile,.975))
 if(!mi){o$cl<-plr(lapply(L,function(d) rob(gl(paste("pmv~x+",M3),d),cl=~hosp))); o$nhosp<-uniqueN(d0$hosp)}
 ## secondary outcomes (propensity score logit for rare outcomes)
 L<-lapply(L,function(d){d<-copy(d); d[,lps:=predict(glm(as.formula(paste("x~",M3)),binomial,d))]; d}); IMP[[dbx]]<<-L
 o$rate<-c(pmv=mean(d0$pmv),d28=mean(d0$d28),hd=mean(d0$hd),los=median(d0$los/24))
 o$s_d28<-pl(lapply(L,function(d) gl("d28~x+lps",d)),"x"); o$s_hd<-pl(lapply(L,function(d) gl("hd~x+lps",d)),"x")
 o$s_los<-pl(lapply(L,function(d) lm(as.formula(paste("log(los)~x+",M3)),d[d28==0&los>0])),"x")
 if(mi){o$rate<-c(o$rate,reint=mean(d0$reint),y1=mean(d0$y1),vfd=median(d0$vfd),mvh=median(d0$mvh))
  o$s_re<-pl(lapply(L,function(d) gl("reint~x+lps",d)),"x"); o$s_reM<-pl(lapply(L,function(d) gl(paste("reint~x+",M2),d)),"x")
  o$s_y1<-pl(lapply(L,function(d){d$t2<-pmin(d$fu_h,8760)/24; coxph(Surv(t2,y1)~x+lps,d)}),"x"); o$s_y1km<-1-pchisq(survdiff(Surv(pmin(fu_h,8760)/24,y1)~x,d0)$chisq,1)
  s<-pl(lapply(L,function(d) lm(as.formula(paste("vfd~x+",M3)),d)),"x"); o$s_vfd<-c(log(s[1:3]),s[4])
  o$s_mvh<-pl(lapply(L,function(d) lm(as.formula(paste("log(mvh)~x+",M3)),d)),"x")
  st<-function(d){d<-copy(d); d[,st:=fifelse(!is.na(death_h)&death_h<=pmin(mvh,672),2L,fifelse(mvh<=672,1L,0L))]; d[,ti:=pmin(fifelse(st==2,death_h,mvh),672)]; d}
  L<-lapply(L,st); IMP[[dbx]]<<-L
  o$s_ext<-plr(lapply(L,function(d){X<-model.matrix(as.formula(paste("~x+",M3)),d)[,-1]; cr<-crr(d$ti,d$st,X); s<-summary(cr)$coef["x",]; c(s[1],s[3])})); d1<-L[[1]]; o$s_extg<-cuminc(d1$ti,d1$st,d1$x)$Tests[1,2]}
 ## sensitivity analyses
 sen<-function(sub=NULL,f=M3,tr="x",y="pmv",mut=NULL){LL<-lapply(L,function(d){d<-copy(d); if(!is.null(mut)) mut(d); if(!is.null(sub)) d<-d[eval(sub)]; d}); dd<-LL[[1]]
  if(nrow(dd)<200||sum(dd[[y]])<25||sum(dd[[tr]]==1)<20) return(NULL); c(pl(lapply(LL,function(d) gl(paste(y,"~",tr,"+",f),d)),tr),nrow(dd))}
 S<-list(); S[["dp24"]]<-sen(mut=function(d) d[,x:=as.integer(dp24>=15)]); S[["dpfirst"]]<-sen(mut=function(d) d[,x:=as.integer(dpfirst>=15)])
 S[["t14"]]<-sen(mut=function(d) d[,x:=as.integer(dp>=14)]); S[["t16"]]<-sen(mut=function(d) d[,x:=as.integer(dp>=16)])
 S[["noards"]]<-sen(quote(ards==0&pna==0)); S[["nocopd"]]<-sen(quote(copd==0)); S[["pf300"]]<-sen(quote(exp(lpf)>=300))
 S[["nopf"]]<-sen(f=sub("\\+lpf","",M3)); S[["vtk"]]<-sen(quote(!is.na(vtk)&vtk>2&vtk<15),paste(M3,"+vtk"))
 S[["cabg"]]<-sen(quote(cabg==1),sub("\\+cabg\\+comb","",M3)); S[["elective"]]<-sen(quote(emerg==0),sub("\\+emerg","",M3)); S[["ndp2"]]<-sen(quote(ndp>=2))
 S[["novaso"]]<-sen(quote(vaso==0&ino==0),sub("\\+vaso\\+ino","",M3))
 S[["median_imp"]]<-{d<-copy(d0); for(v in c("wt",M3v(dbx))) if(v%in%names(d)) set(d,j=v,value=fifelse(is.na(d[[v]]),median(d[[v]],na.rm=T),d[[v]])); d[,sevz:=(sev-mean(sev))/sd(sev)]; c(pr(coef(gl(paste("pmv~x+",M3),d))["x"],sqrt(vcov(gl(paste("pmv~x+",M3),d))["x","x"])),nrow(d))}
 S[["cc"]]<-{d<-na.omit(d0[,c("pmv","x",strsplit(gsub(" ","",M2),"\\+")[[1]],M3v(dbx)),with=FALSE]); d[,sevz:=(sev-mean(sev))/sd(sev)]; f<-gl(paste("pmv~x+",M3),d); c(pr(coef(f)["x"],sqrt(vcov(f)["x","x"])),nrow(d))}
 S[["ps_adj"]]<-c(pl(lapply(L,function(d) gl("pmv~x+lps",d)),"x"),nrow(d0))
 Pm<-"age+male+wt+cabg+comb+sevz+peep+lpf"; S[["pars"]]<-c(pl(lapply(L,function(d) gl(paste("pmv~x+",Pm),d)),"x"),nrow(d0))
 S[["firth"]]<-c(plr(lapply(L,function(d){fz<-logistf(as.formula(paste("pmv~x+",M3)),data=as.data.frame(d),plconf=2); c(coef(fz)[2],sqrt(diag(vcov(fz)))[2])})),nrow(d0))
 if(mi){S[["p48"]]<-sen(y="p48",mut=function(d) d[,p48:=as.integer(mvh>48)]); S[["comp"]]<-sen(y="pr",mut=function(d) d[,pr:=as.integer(pmv==1|reint==1)])}
 o$S<-S
 sgv<-list(age65=quote(age>=65),male=quote(male==1),cabg=quote(cabg==1),emerg=quote(emerg==1),chf=quote(chf==1),copd=quote(copd==1),ckd=quote(ckd==1),obese=quote(obese==1),pf200=quote(exp(lpf)<200),ino=quote(ino==1))
 o$sg<-lapply(sgv,function(z){LL<-lapply(L,function(d){d<-copy(d); d[,zz:=as.integer(eval(z))]; d}); pint<-tryCatch(pl(lapply(LL,function(d) gl(paste("pmv~x*zz+",M3),d)),"x:zz")[4],error=function(e) NA)
  g<-function(lv) tryCatch({if(sum(LL[[1]][zz==lv]$pmv)<15) matrix(NA,1,4) else matrix(pl(lapply(LL,function(d) gl(paste("pmv~x+",M3),d[zz==lv])),"x"),1)},error=function(e) matrix(NA,1,4))
  list(yes=g(1),no=g(0),nyes=sum(LL[[1]]$zz==1),nno=sum(LL[[1]]$zz==0),pint=pint)})
 o$desc<-d0[,.(n=.N,age=mean(age),male=mean(male),dp=median(dp),peep=median(peep,na.rm=T),pplat=median(pplat),vt=median(vt,na.rm=T),pf=median(pf,na.rm=T),sev=mean(sev,na.rm=T),vaso=mean(vaso),ino=mean(ino),cabg=mean(cabg),valve=mean(ptype=="valve"),comb=mean(comb),emerg=mean(emerg),chf=mean(chf),copd=mean(copd),ckd=mean(ckd),dm=mean(dm),af=mean(af),lact=median(lact,na.rm=T),creat=median(creat,na.rm=T),pmv=mean(pmv),d28=mean(d28),hd=mean(hd)),keyby=x]
 o}
R<-list(); for(kk in c("MIMIC-IV","eICU-CRD")){cat(kk,"\n"); R[[kk]]<-an(kk)}
R$miss<-miss; R$M3<-list(m=M3f("MIMIC-IV"),e=M3f("eICU-CRD"))
saveRDS(R,"R.rds"); saveRDS(IMP,"Imp.rds"); saveRDS(rbind(IMP[[1]][[1]],IMP[[2]][[1]],fill=TRUE),"Dd.rds")
for(kk in c("MIMIC-IV","eICU-CRD")){o<-R[[kk]]; cat("\n==",kk,o$n,o$M3,"\n"); print(round(rbind(crude=o$crude,m2=o$m2,m3=o$m3,cont=o$cont,o$cat,trend=o$trend,iptw=o$iptw,ow=o$ow,dr=o$dr,psm=o$psm[1:4]),4)); print(round(c(o$rcs,o$inc,o$smd,o$ev),4)); print(round(o$ar,4))
 for(s in grep("^s_",names(o),value=T)) cat(s,round(o[[s]],4),"\n"); print(round(do.call(rbind,o$S),4)); print(sapply(o$sg,function(z) round(c(z$yes[1],z$yes[4],z$no[1],z$no[4],z$pint),3)))}
