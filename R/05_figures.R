# Figures 1-4 (TIFF 300 dpi, Times New Roman, white background)
suppressMessages({library(data.table);library(ggplot2);library(splines);library(cmprsk)})
B<-readRDS("R.rds"); A<-as.data.table(readRDS("Dd.rds")); fl<-readRDS("fl.rds"); O<-readRDS("M.rds"); FT<-"Times New Roman"; dir.create("fig",showWarnings=FALSE)
sv<-function(p,f,w,h){ggsave(file.path("fig",paste0(f,".tiff")),p,width=w,height=h,dpi=300,compression="none",bg="white"); ggsave(file.path("fig",paste0(f,".png")),p,width=w,height=h,dpi=150,bg="white")}
fm<-function(n) format(n,big.mark=",")
## ---- Figure 1: flow chart (arrow heads stop short of boxes)
bx<-function(x,y,w,h,lab,id) data.table(id=id,x=x,y=y,w=w,h=h,lab=lab)
mk<-function(x0,db,main,exc){B1<-list(); E1<-list(); n<-length(main); ys<-seq(10,10-2*(n-1),by=-2)
 for(i in 1:n) B1[[i]]<-bx(x0,ys[i],3.2,1.05,main[i],paste0(db,"m",i))
 for(i in seq_along(exc)) if(!is.na(exc[i])) E1[[length(E1)+1]]<-bx(x0+3.55,ys[i]-1,2.9,0.9,exc[i],paste0(db,"e",i))
 list(b=rbindlist(B1),e=rbindlist(E1),ys=ys)}
M<-mk(0,"M",c(sprintf("MIMIC-IV v3.1\nFirst ICU stays of adults\nn = %s",fm(fl$m0)),sprintf("Open CABG and/or valve surgery\nn = %s",fm(fl$m_cs)),sprintf("ICU admission on day of surgery\nor next day\nn = %s",fm(fl$m_postop)),sprintf("Invasively ventilated on ICU arrival\nn = %s",fm(fl$m_vent)),sprintf("Plateau pressure and PEEP\nrecorded within 0–6 h\nn = %s",fm(fl$m_dp)),sprintf("Analysed cohort\nn = %s",fm(fl$m_d24))),
 c(NA,sprintf("Not admitted directly\nafter surgery: n = %s",fm(fl$m_cs-fl$m_postop)),sprintf("ECMO, heart transplant or\nVAD: n = %s; not ventilated\non arrival: n = %s",fm(fl$m_postop-fl$m_ecmo),fm(fl$m_ecmo-fl$m_vent)),sprintf("No plateau pressure or\nPEEP within 0–6 h: n = %s",fm(fl$m_vent-fl$m_dp)),sprintf("Death within 24 h:\nn = %s",fm(fl$m_dp-fl$m_d24))))
E<-mk(8,"E",c(sprintf("eICU-CRD v2.0\nFirst ICU stays of adults with\nCABG and/or valve surgery\nn = %s",fm(fl$e_cs)),sprintf("Invasively ventilated on ICU arrival\nn = %s",fm(fl$e_vent)),sprintf("Plateau pressure and PEEP\nrecorded within 0–6 h\nn = %s",fm(fl$e_dp)),sprintf("Analysed cohort\nn = %s",fm(fl$e_d24))),
 c(sprintf("Not ventilated on arrival:\nn = %s",fm(fl$e_cs-fl$e_vent)),sprintf("No plateau pressure or\nPEEP within 0–6 h: n = %s",fm(fl$e_vent-fl$e_dp)),sprintf("Death within 24 h:\nn = %s",fm(fl$e_dp-fl$e_d24))))
c1<-O$cnt; gb<-rbind(bx(-0.95,-2.2,1.5,0.9,sprintf("DP <15 cmH₂O\nn = %s",fm(c1[db=="MIMIC-IV"]$n0)),"Mg0"),bx(0.95,-2.2,1.5,0.9,sprintf("DP ≥15 cmH₂O\nn = %s",fm(c1[db=="MIMIC-IV"]$n1)),"Mg1"),
 bx(7.05,1.8,1.5,0.9,sprintf("DP <15 cmH₂O\nn = %s",fm(c1[db=="eICU-CRD"]$n0)),"Eg0"),bx(8.95,1.8,1.5,0.9,sprintf("DP ≥15 cmH₂O\nn = %s",fm(c1[db=="eICU-CRD"]$n1)),"Eg1"))
BX<-rbind(M$b,M$e,E$b,E$e,gb); g<-0.12
ar<-list(); vs<-function(b) for(i in 1:(nrow(b)-1)) ar[[length(ar)+1]]<<-data.table(x=b$x[i],y=b$y[i]-b$h[i]/2,xe=b$x[i+1],ye=b$y[i+1]+b$h[i+1]/2+g)
vs(M$b); vs(E$b)
hz<-function(b,e){for(i in 1:nrow(e)){yy<-e$y[i]; ar[[length(ar)+1]]<<-data.table(x=b$x[1],y=yy,xe=e$x[i]-e$w[i]/2-g,ye=yy)}}
hz(M$b,M$e); hz(E$b,E$e)
sp<-function(b,gx){lb<-b[.N]; for(i in 1:nrow(gx)){ar[[length(ar)+1]]<<-data.table(x=lb$x,y=lb$y-lb$h/2,xe=gx$x[i],ye=gx$y[i]+gx$h[i]/2+g)}}
sp(M$b,gb[1:2]); sp(E$b,gb[3:4]); AR<-rbindlist(ar)
p1<-ggplot()+geom_rect(data=BX,aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),fill="white",colour="black",linewidth=.4)+geom_text(data=BX,aes(x,y,label=lab),family=FT,size=2.6,lineheight=.95)+
 geom_segment(data=AR,aes(x=x,y=y,xend=xe,yend=ye),arrow=arrow(length=unit(1.6,"mm"),type="closed"),linewidth=.35)+coord_cartesian(xlim=c(-1.9,13.2),ylim=c(-2.8,10.7),expand=FALSE)+theme_void(base_family=FT)+theme(plot.background=element_rect(fill="white",colour=NA))
sv(p1,"Figure1",8.2,7.2)
## ---- Figure 2: RCS
IMP<-readRDS("Imp.rds")
source("../R/rcs.R")
mkr<-function(k){L<-IMP[[k]]; M3<-B[[k]]$M3; kn<-B[[k]]$knots; d0<-L[[1]]; gg<-seq(quantile(d0$dp,.02),quantile(d0$dp,.98),length.out=120)
 C<-cbind(gg-10,sweep(rcsb(gg,kn),2,rcsb(10,kn)))
 z<-lapply(L,function(d){d<-copy(d); bb<-rcsb(d$dp,kn); d[,`:=`(s1=bb[,1],s2=bb[,2])]; f<-glm(as.formula(paste("pmv~dp+s1+s2+",M3)),binomial,d); j<-c("dp","s1","s2"); list(e=C%*%coef(f)[j],v=diag(C%*%vcov(f)[j,j]%*%t(C)))})
 E<-sapply(z,`[[`,"e"); V<-sapply(z,`[[`,"v"); m<-length(L); qb<-rowMeans(E); se<-sqrt(rowMeans(V)+(1+1/m)*apply(E,1,var))
 o<-B[[k]]; data.table(dp=gg,or=exp(qb),lo=exp(qb-1.96*se),hi=exp(qb+1.96*se),db=sprintf("%s\nP overall %s; P nonlinear = %s",k,ifelse(o$rcs[1]<0.001,"< 0.001",sprintf("= %.3f",o$rcs[1])),sprintf("%.3f",o$rcs[2])))}
fd<-rbindlist(lapply(c("MIMIC-IV","eICU-CRD"),mkr)); fd[,db:=factor(db,levels=unique(db))]
p2<-ggplot(fd,aes(dp,or))+geom_ribbon(aes(ymin=lo,ymax=hi),fill="grey80")+geom_line(linewidth=.7)+geom_hline(yintercept=1,linetype=2,linewidth=.4)+geom_vline(xintercept=15,linetype=3,linewidth=.4)+facet_wrap(~db,scales="free_y")+scale_y_log10()+
 labs(x=expression("Driving pressure, median 0–6 h after ICU admission (cm"*H[2]*"O)"),y="Adjusted odds ratio for ventilation >24 h")+theme_classic(base_size=10,base_family=FT)+theme(strip.background=element_blank(),strip.text=element_text(face="bold"))
sv(p2,"Figure2",7.5,3.6)
## ---- Figure 3: cumulative incidence of extubation (MIMIC-IV)
d<-IMP[["MIMIC-IV"]][[1]]; cu<-cuminc(d$ti,d$st,d$x); cf<-rbindlist(lapply(c("0 1","1 1"),function(nm){z<-cu[[nm]]; data.table(t=z$time/24,y=z$est,g=ifelse(nm=="1 1","DP ≥15 cmH₂O","DP <15 cmH₂O"))}))
p3<-ggplot(cf[t<=7],aes(t,y,linetype=g))+geom_step(linewidth=.6)+scale_linetype_manual(values=c("solid","dashed"),name=NULL)+scale_y_continuous(labels=function(z) paste0(100*z,"%"),limits=c(0,1))+scale_x_continuous(breaks=0:7)+
 labs(x="Days after ICU admission",y="Cumulative incidence of extubation")+annotate("text",x=3.0,y=.35,label=sprintf("Gray's test P < 0.001\nAdjusted SHR %.2f (95%% CI %.2f–%.2f)",B[["MIMIC-IV"]]$s_ext[1],B[["MIMIC-IV"]]$s_ext[2],B[["MIMIC-IV"]]$s_ext[3]),family=FT,size=3.2,hjust=0)+
 theme_classic(base_size=10,base_family=FT)+theme(legend.position=c(.8,.6))
sv(p3,"Figure3",6,4.2)
## ---- Figure 4: subgroup forest (both databases)
en<-c(age65="Age ≥65 years",male="Male sex",cabg="Isolated CABG",emerg="Urgent or emergency admission",chf="Heart failure",copd="COPD",ckd="Chronic kidney disease",obese="Obesity",pf200="PaO₂/FiO₂ <200 mmHg",ino="Inotrope use")
fp<-rbindlist(lapply(c("MIMIC-IV","eICU-CRD"),function(k){sg<-B[[k]]$sg; rbindlist(lapply(names(sg),function(s){z<-sg[[s]]; rbind(data.table(db=k,var=en[s],lv="Yes",n=z$nyes,or=z$yes[1,1],lo=z$yes[1,2],hi=z$yes[1,3],pint=z$pint),data.table(db=k,var=en[s],lv="No",n=z$nno,or=z$no[1,1],lo=z$no[1,2],hi=z$no[1,3],pint=z$pint))}))}))
fp<-fp[db=="MIMIC-IV"]; fp[,lab:=paste0(var,": ",tolower(lv))]; fp[,lab:=factor(lab,levels=rev(lab))]; fp[,txt:=sprintf("%.2f (%.2f–%.2f)",or,lo,hi)]; fp[,pt:=ifelse(lv=="Yes",sprintf("%.2f",pint),"")]
p4<-ggplot(fp,aes(or,lab))+geom_vline(xintercept=1,linetype=2,linewidth=.4)+geom_errorbar(aes(xmin=lo,xmax=hi),width=.25,orientation="y",linewidth=.4)+geom_point(shape=15,size=1.8)+
 geom_text(aes(x=6.2,label=txt),family=FT,size=2.8,hjust=0)+geom_text(aes(x=15,label=pt),family=FT,size=2.8,hjust=0)+annotate("text",x=c(6.2,15),y=nrow(fp)+.9,label=c("OR (95% CI)","P for\ninteraction"),family=FT,size=2.9,hjust=0,fontface="bold")+
 scale_x_log10(breaks=c(.5,1,2,4))+coord_cartesian(clip="off",xlim=c(.4,5),ylim=c(1,nrow(fp)+1))+labs(x="Adjusted odds ratio for ventilation >24 h\n(DP ≥15 vs <15 cmH₂O)",y=NULL)+theme_classic(base_size=9.5,base_family=FT)+theme(plot.margin=margin(5,190,5,5))
sv(p4,"Figure4",7.6,6)
cat("figs ok\n")
