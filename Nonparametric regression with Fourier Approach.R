library(readxl)
library(gtools)
library(MASS) 
library(pracma)

data <- read_excel("Dataset.xlsx")
data

#Response variable
y=as.matrix(data[,5]) 
#Predictor variables 
x=as.matrix(data[,c(2,3,4)])
n=nrow(data)                #number of observations

#Scatterplots
par(mfrow = c(1, 3))
plot(data$Permintaan, data$Harga,pch=19, col="red", xlab="Permintaan",ylab="Harga")
plot(data$Persediaan, data$Harga,pch=19, col="red", xlab="Persediaan",ylab="Harga")
plot(data$Produksi, data$Harga,pch=19, col="red", xlab="Produksi",ylab="Harga")

#Linear Parametric Regression for Permintaan
Permintaan=as.matrix(data[,2])
reg<-lm(y~Permintaan, data)
summary(reg)

#Nonparametric Regression with Fourier Series Approach for Variable BOPO
fourierseries<-function(y, x, K){ 
  n=length(y) 
  a<-(1*(K+1))+1 
  Z<-matrix(0, n, a) 
  result<-matrix(0, K, 4) 
  for (k in 1:K) 
  { 
    for(i in 1:n) 
    { 
      for(j in 1:k) 
      { 
        Z[i,1]=1/2 
        Z[i,2]<-x[i] 
        Z[i,2+j]=cos(j*x[i])
      } 
    } 
    I<-diag(1,n,n) 
    M<- Z%*%ginv(t(Z)%*%Z)%*%t(Z) 
    Ytop<-M%*%y
    W<-(1/n)*I 
    MSE<-t(y-Ytop) %*% W %*% (y-Ytop) 
    Q<-(n^-1*sum(diag(I-M)))^2 
    GCV<-MSE/Q
    SSE<-sum((y-Ytop)^2)
    SST<-sum((y-mean(y))^2)
    Rsq=(1-(SSE/(SST)))*100
    result[k,1]<-k 
    result[k,2]<-GCV
    result[k,3]<-Rsq
    result[k,4]<-MSE
  } 
  print(result) 
  
  GCV2<-min(result[,2]) 
  s<-1
  repeat{ 
    if(result[s,2]==GCV2) 
    { 
      kOpt<-result[s,1] 
      GCVOpt<-GCV2 
      break 
    } 
    else s<-s+1 
  } 
  cat("optimal number of K is \t",kOpt,"\n")
  cat("the smallest GCV value is \t",GCV2,"\n")
  print(GCV2)
} 
A=fourierseries(y, Permintaan, K=10)

##MAIN RESULT USING ALL THE PREDICTOR VARIABLES
#Combinations Number of K
p=ncol(x)                   #number of predictor variables
K=3                       #number of K
K.com=permutations(K, p, c(1:K), repeats.allowed = TRUE)

#GCV Value and Hypothesis testing for all the combinations of K
Best=data.frame(matrix(0,1,15)) 
colnames(Best)=c("Var Permintaan","Var Persediaan","Var Produksi", "K.opt Var Permintaan","K.opt Var Persediaan", "K.opt Var Produksi", "GCV", "MSE", "R.sq (%)", "d1", "d2", "Fvalue", "Ftable", "Pvalue", "Decisions")
X.first=c(1:(K+1))
Best[,1:p]=colnames(x) 
h=0
for (k in 1:nrow(K.com)){
  X=matrix(0,nrow(x)) 
  for (l in 1:p)
  {
    X1=matrix(0,nrow(x)) 
    for (f in 1:K.com[k,l])
    {
      cosine=cos(X.first[f]*x[,l])
      X1=cbind(X1,cosine)
    }
    X1=X1[,-1] 
    X1=cbind(1/2,x[,l],X1)
    X=cbind(X,X1)
  }
  X=X[,-1]
  GCV=0; R.sq=0; MSE=0; d1=0; d2=0; Fvalue=0; Ftable=0; Pvalue=0; Decision=0
  
  I=diag(n)
  V=X%*%ginv((t(X)%*%X))%*%t(X) 
  MSE=(n^-1)*t(y)%*%t(I-V)%*%(I-V)%*%y 
  GCV.all=MSE/(((n^-1)*sum(diag(I-V)))^2) 
  B=ginv((t(X)%*%X))%*%t(X)%*%y
  y.pred=V%*%y
  R2=(1-(sum((y-y.pred)^2)/sum((y-mean(y))^2)))*100
  mse=sum((y-y.pred)^2)/(n)
  MSE=rbind(MSE,mse) 
  R.sq=rbind(R.sq,R2) 
  GCV=rbind(GCV,GCV.all)
  
  d1.all=sum(diag(V))
  d2.all=n-sum(diag(V))
  F1=t(B)%*%t(X)%*%y/d1.all
  F2=t(y-X%*%B)%*%(y-X%*%B)/d2.all
  Fvalue.all=F1/F2
  alfa=0.05
  Ftable.all=qf(alfa,d1.all,d2.all, lower.tail = FALSE)
  Pvalue.all=pf(Fvalue.all,d1.all,d2.all, lower.tail= FALSE)
  Decision.all=ifelse(Fvalue.all >= Ftable.all, "H0 Rejected", "H0 Accepted")
  d1=rbind(d1,d1.all) 
  d2=rbind(d2,d2.all) 
  Fvalue=rbind(Fvalue,Fvalue.all)
  Ftable=rbind(Ftable,Ftable.all) 
  Pvalue=rbind(Pvalue,Pvalue.all)
  Decision=rbind(Decision,Decision.all)
  
  list=cbind(GCV[-1],MSE[-1],R.sq[-1],d1[-1],d2[-1],Fvalue[-1],Ftable[-1],Pvalue[-1],Decision[-1]) 
  best_comb=which(list[,1]==min(list[,1])) 
  hh=list[best_comb,]
  h=rbind(h,hh)
}
Optimal=h[-1,] 
rr=matrix(Optimal,nrow=nrow(Optimal)) 
list2=cbind(K.com,rr)
colnames(list2)=c("K.opt Var Permintaan","K.opt Var Persediaan", "K.opt Var Produksi", "GCV", "MSE", "R.sq (%)", "d1", "d2", "Fvalue", "Ftable", "Pvalue", "Decisions")

#The Optimum Combinations of K Based on the Smallest GCV Value
K.optimum=which(list2[,4]==min(list2[,4])) 
ss=as.vector(list2[K.optimum,])
Best[,4:15]=ss 
Best_combination=Best

#Parameters Estimation for the Optimal Combinations of K
x.opt=matrix(0,ncol=1,nrow=nrow(data)) 
for (xb in 1:p){
  tb=data[Best_combination[,xb]] 
  x.opt=cbind(x.opt,tb)
}
x.opt=as.matrix(x.opt[-1])

K.opt=matrix(0,ncol=1,nrow=1) 
for (zb in 4:6){
  Kb=Best_combination[,zb] 
  K.opt=cbind(K.opt,Kb)
}
K.opt=as.matrix(K.opt[-1]) 
rownames(K.opt)=c("K.opt Var Permintaan","K.opt Var Persediaan","K.opt Var Produksi")

X.opt=matrix(0,nrow(x.opt)) 
for (ll in 1:p)
{
  X1.opt=matrix(0,nrow(x.opt)) 
  for (ff in 1:K.opt[ll])
  {
    cosine=cos(X.first[ff]*x.opt[,ll]) 
    X1.opt=cbind(X1.opt,cosine)
  }
  X1.opt=X1.opt[,-1] 
  X1.opt=cbind(1/2,x.opt[,ll],X1.opt) 
  X.opt=cbind(X.opt,X1.opt)
}
X.opt=X.opt[,-1] 

V.opt=X.opt%*%ginv(t(X.opt)%*%X.opt)%*%t(X.opt) 
MSE.opt=n^-1*t(y)%*%t(I-V.opt)%*%(I-V.opt)%*%y
GCV.opt=MSE.opt/((n^-1*sum(diag(I-V.opt)))^2)
y.hat=V.opt%*%y
R2.opt=(1-(sum((y-y.hat)^2)/sum((y-mean(y))^2)))*100
mse.opt=sum((y-y.hat)^2)/(n)
B.opt=ginv(t(X.opt)%*%X.opt)%*%t(X.opt)%*%y
colnames(B.opt)=c("Estimated Parameters")
B.opt

#Hypothesis Testing for the Optimal Combinations of K
d1.opt=sum(diag(V.opt))
d2.opt=n-sum(diag(V.opt))
F1.opt=t(B.opt)%*%t(X.opt)%*%y/d1.opt
F2.opt=t(y-X.opt%*%B.opt)%*%(y-X.opt%*%B.opt)/d2.opt
Fvalue.opt=F1.opt/F2.opt
Ftable.opt=qf(alfa,d1.opt,d2.opt, lower.tail = FALSE)
Pvalue.opt=pf(Fvalue.opt,d1.opt,d2.opt, lower.tail= FALSE)
Hypothesis_table <- data.frame(GCV = as.numeric(GCV.opt),
                               df1 = d1.opt,
                               df2 = d2.opt,
                               Ftest = Fvalue.opt,
                               Ftable = Ftable.opt,
                               PValue = Pvalue.opt,
                               R2 = R2.opt)
colnames(Hypothesis_table) <- c("GCV","df1","df2","Ftest",
                                "Ftable","PValue","R2")
rownames(Hypothesis_table) <- NULL
Hypothesis_table
cat("\n Optimal Combinatios of K =",ss[1:p],"\n")

#Decision
if (Fvalue.opt >= Ftable.opt) "The null hypothesis is rejected (at least one of the parameters is not zero)" else "The null hypothesis fails to be rejected (all the parameters are zero)"