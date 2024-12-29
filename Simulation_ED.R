rm(list=ls())
library(tweedie)

set.seed(0)

theta<-c(0.132680378342590,1.47311344871203,1.10950195745342,0.457467403417009,1.82688627901974)
u<-theta[1]
p<-theta[2]
lam<-theta[3]
b<-theta[4]
a<-theta[5]

ff=100
t<-seq(from=0,to=3000,by=ff)

J<-3
n<-30
m<-length(t)

S<-c(65,85,100)
s<-rep(NA,J)
S0<-40
SH<-100
for(k in 1:J){
  s[k]=(1/(S0+273.15)-1/(S[k]+273.15))/(1/(S0+273.15)-1/(SH+273.15))
}


tau<-c(0,0,1000,2000)
l1<-1000/ff
l2<-1000/ff
l3<-1000/ff
lll<-c(l1,l2,l3)
ll<-c(0,l1,l1+l2,l1+l2+l3)
ss<-c(0,0,s[1],s[2],s[3])


y<-array(NA,dim=c(J,n,max(lll)+1))
P1<-array(NA,dim=c(J,max(lll)))
P2<-array(NA,dim=c(J,max(lll)))
rou<-array(NA,dim=c(J,max(lll)))
y[1,,1]<-rep(0,n)

for(k in 1:J){
    P1[k,1:lll[k]]<-(t[(2+ll[k]):(ll[k+1]+1)]-tau[k+1])*exp(a*ss[k+2]/b)+exp(a*ss[k+1]/b)*(tau[k+1]-tau[k])+exp(a*ss[k]/b)*tau[k]
    P2[k,1:lll[k]]<-(t[(1+ll[k]):(ll[k+1])]-tau[k+1])*exp(a*ss[k+2]/b)+exp(a*ss[k+1]/b)*(tau[k+1]-tau[k])+exp(a*ss[k]/b)*tau[k]
    rou[k,1:lll[k]]<-P1[k,1:lll[k]]^b-P2[k,1:lll[k]]^b
}

for(k in 1:J){
  j<-1
  if(k>1){
    y[k,,1]<-y[k-1,,max(lll)+1] 
  }
  for(j in 1:max(lll)){
    dy<-rtweedie(n, mu=u*rou[k,j], phi=1/(lam*rou[k,j]^(p-1)), power=p)
    while(any(dy==0)){
      dy<-rtweedie(n, mu=u*rou[k,j], phi=1/(lam*rou[k,j]^(p-1)), power=p)
    }
    y[k,,j+1]<-y[k,,j]+dy
  }
}

plot(t[1:(l1+1)],y[1,1,1:(l1+1)],col="red",type="b",xlim=c(0,3000),ylim=c(0,30))
for(i in 1:n){
  points(t[1:(l1+1)],y[1,i,1:(l1+1)],col="red",type="b")
  points(t[(l1+1):(l1+l2+1)],y[2,i,1:(l2+1)],col="blue",type="b")
  points(t[(l1+l2+1):(l1+l2+l3+1)],y[3,i,1:(l3+1)],col="green",type="b")
}

write.csv(y[1,,],"y1.csv")
write.csv(y[2,,],"y2.csv")
write.csv(y[3,,],"y3.csv")
