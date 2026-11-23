# restricted cubic spline basis (Harrell), knots at the 5th, 35th, 65th and 95th percentiles
rcsb<-function(x,kn){k<-length(kn); p<-function(z) pmax(z,0)^3; sapply(1:(k-2),function(j) (p(x-kn[j])-p(x-kn[k-1])*(kn[k]-kn[j])/(kn[k]-kn[k-1])+p(x-kn[k])*(kn[k-1]-kn[j])/(kn[k]-kn[k-1]))/(kn[k]-kn[1])^2)}
