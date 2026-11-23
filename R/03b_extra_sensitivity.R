# 02b additional sensitivity analyses (ventilator mode, harmonised outcome definition, selection weighting) in MIMIC-IV
suppressMessages({library(data.table);library(sandwich)})
R<-readRDS("R.rds"); IMP<-readRDS("Imp.rds"); L1<-readRDS("r1_mimic.rds"); h25<-readRDS("r1_h25.rds"); M3<-R[["MIMIC-IV"]]$M3
gl<-function(f,d,w=NULL) suppressWarnings(glm(as.formula(f),if(is.null(w)) binomial else quasibinomial,d,weights=w))
rub<-function(b,v){m<-length(b); Tt<-mean(v)+(1+1/m)*var(b); c(mean(b),sqrt(Tt))}; pr<-function(e,s) c(exp(e),exp(e-1.96*s),exp(e+1.96*s),2*pnorm(-abs(e/s)))
run<-function(fun){z<-lapply(IMP[["MIMIC-IV"]],fun); r<-rub(sapply(z,`[`,1),sapply(z,`[`,2)^2); c(pr(r[1],r[2]),z[[1]][3])}
cm<-copy(L1$cm); setnames(cm,"dpc","dpmode"); w<-L1$V[,.(id=stay_id,psel)]
S<-R[["MIMIC-IV"]]$S
S[["mode_ctl"]]<-run(function(d){d<-merge(d,cm,by.x="id",by.y="stay_id"); d[,x:=as.integer(dpmode>=15)]; f<-gl(paste("pmv~x+",M3),d); c(coef(f)["x"],sqrt(vcov(f)["x","x"]),nrow(d))})
S[["mode_nospont"]]<-run(function(d){d<-d[!id%in%L1$spn]; f<-gl(paste("pmv~x+",M3),d); c(coef(f)["x"],sqrt(vcov(f)["x","x"]),nrow(d))})
S[["harm"]]<-run(function(d){d<-copy(d); d[,p25:=as.integer(id%in%h25)]; f<-gl(paste("p25~x+",M3),d); c(coef(f)["x"],sqrt(vcov(f)["x","x"]),nrow(d))})
S[["selw"]]<-run(function(d){d<-merge(d,w,by="id"); f<-gl(paste("pmv~x+",M3),d,1/d$psel); V<-vcovHC(f,type="HC0"); c(coef(f)["x"],sqrt(V["x","x"]),nrow(d))})
R[["MIMIC-IV"]]$S<-S
d0<-IMP[["MIMIC-IV"]][[1]]; d0[,p25:=as.integer(id%in%h25)]; R$agree<-d0[,.(n=.N,both=sum(p25==1&pmv==1),p25only=sum(p25==1&pmv==0),pmvonly=sum(p25==0&pmv==1))]
R$sel<-L1$V[,.(n=.N,age=median(age),male=mean(male),cabg=mean(ptype=="cabg"),valve=mean(ptype=="valve"),comb=mean(ptype=="comb"),emerg=mean(emerg),pmv=mean(pmv),mvh=median(ve),d28=mean(d28)),keyby=hasdp]
md<-readRDS("r1_O.rds"); R$modecov<-md$modecov; R$modes<-md$modes
two<-copy(d0[ndp>=2]); R$two<-list(n=nrow(two),ev=sum(two$pmv),pmv1=mean(d0[ndp==1]$pmv),pmv2=mean(two$pmv),cor=cor(two$dpfirst,two$dp),crude=pr(coef(gl("pmv~x",two))["x"],sqrt(vcov(gl("pmv~x",two))["x","x"])))
saveRDS(R,"R.rds"); for(s in c("mode_ctl","mode_nospont","harm","selw")) cat(s,round(S[[s]],4),"\n"); print(R$agree); print(round(R[["eICU-CRD"]]$cl,4)); cat("hosp",R[["eICU-CRD"]]$nhosp,"\n"); str(R$two)
