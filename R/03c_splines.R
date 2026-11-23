# 02c dose-response: restricted cubic splines pooled across imputations
suppressMessages({library(data.table);library(mice)}); source("../R/rcs.R")
R<-readRDS("R.rds"); IMP<-readRDS("Imp.rds")
for(k in c("MIMIC-IV","eICU-CRD")){L<-IMP[[k]]; M3<-R[[k]]$M3; kn<-quantile(L[[1]]$dp,c(.05,.35,.65,.95),type=1); kn<-unique(kn); if(length(kn)<4) kn<-quantile(L[[1]]$dp,c(.05,.35,.65,.95)) 
 L<-lapply(L,function(d){d<-copy(d); b<-rcsb(d$dp,kn); d[,`:=`(s1=b[,1],s2=b[,2])]; d})
 F<-function(f) as.mira(lapply(L,function(d) glm(as.formula(f),binomial,d)))
 fr<-F(paste("pmv~dp+s1+s2+",M3)); f0<-F(paste("pmv~",M3)); fl<-F(paste("pmv~dp+",M3))
 R[[k]]$rcs<-c(D1(fr,f0)$result[1,4],D1(fr,fl)$result[1,4]); R[[k]]$knots<-kn; cat(k,"knots",kn,"P overall",R[[k]]$rcs[1],"P nonlinear",R[[k]]$rcs[2],"\n")}
saveRDS(R,"R.rds")
