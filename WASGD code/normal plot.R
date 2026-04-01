library(mvtnorm) 
#############Squared Loss######################
N=450
#you can try N=200 if N=450 is beyond the limitation
n=100000
d=5
alpha=0.505
#We choose step size constant to be 1, as same as the simulation settings in Zhu et al.(2021)
Sigma=diag(1,5)
xt=c(1,-2,0,0,4)
eta=3
w=matrix(0,d,N)
ws=matrix(0,d,N)
ww=matrix(0,d,N)

#functions
grad_comp = function(x_old, a_new, b_new, type = "Linear"){
  if (type == "Linear"){
    grad = t(a_new)*c(a_new%*%x_old - b_new)
  }else{
    grad = -b_new*t(a_new)/c(1+exp(b_new*a_new%*%x_old))
  }
  return (grad)
}


generate_data = function(n,d, Sigma, sigma=1, x_true, type = "Linear"){
  
  linear=matrix(0,d+1,n)
  for(i in 1:n){
    a = rmvnorm(1, sigma = Sigma)
    b = c(a%*%x_true) + rnorm(1, sd = sigma) 
    linear[,i]=c(b,a)
  }
  
  
  return(linear)
}

iterate=function(M,alpha){
  n=dim(M)[2]
  d=dim(M)[1]-1
  X=matrix(0,d,n)
  X[,1]=rep(0,d)
  for(i in 2:n){
    X[,i]=X[,i-1]-i^(-alpha)*grad_comp(X[,i-1],M[-1,i],M[1,i]
                                       ,type = "Linear")
  }
  return(X)
}

#generate data
datasquare=list()

for(j in 1:N){
  set.seed(j)
  datasquare[[j]]=generate_data(n,d,Sigma,1,xt,type="Linear")
}


ite=list()

for(j in 1:N){
  ite[[j]]=iterate(datasquare[[j]],alpha)
}


#Compute weighted averaged SGD

#Polynomial Decay
  for(j in 1:N){
    x=rep(0,d)
    for(i in 2:n){
      w[,j]=(1-(eta+1)/(eta+i))*w[,j]+(eta+1)/(eta+i)*ite[[j]][,i]
    }
  }


#Suffix
for(j in 1:N){
  ws[,j]=apply(ite[[j]][,(n/2+1):n],1,mean)
}



#Optimal Weights
  for(j in 1:N){
    ww[,j]=(1-2^(alpha-1))*ite[[j]][,1]+2^(alpha-1)*ite[[j]][,2]
    for(i in 3:n){
      ww[,j]=(i-1)/i*ww[,j]+
        (1-i^alpha)*ite[[j]][,i-1]/i+i^(alpha-1)*ite[[j]][,i]
    }
  }




compo=2
plot(density((w[compo,]-xt[compo])*
               sqrt(n*(2*eta+1))/(eta+1)),xlab="",ylab="",col="blue",main="",
     xlim=c(-6,6))

lines(density((ws[compo,]-xt[compo])*
                sqrt(n*0.5)),col="green")
lines(density((ww[compo,]-xt[compo])*
                sqrt(n)),col="red")

curve(dnorm,add=TRUE,lwd=1.5,lty=4)
lines(density((w[compo,]-xt[compo])*
                sqrt(n)),col="blue",lty=2)
lines(density((ws[compo,]-xt[compo])*
                sqrt(n)),col="green",lty=2)

legend("topleft", c("standard normal distribution",
                    "Polynomial-Decay (correctly standardized)",
                    "0.5-suffix (correctly standardized)",
                    "Optimal Weights",
                    "Polynomial-Decay (unstandardized)",
                    "0.5-suffix (unstandardized)"),
       lty=c(1,1,1,2,2),col=c("black","blue","green","red","blue","green"),
       lwd=c(1.5,1,1,1,1),cex=0.43)




#############Logit Loss######################
grad_comp = function(x_old, a_new, b_new, type = "Linear"){
  if (type == "Linear"){
    grad = t(a_new)*c(a_new%*%x_old - b_new)
  }else{
    grad = -b_new*t(a_new)/c(1+exp(b_new*a_new%*%x_old))
  }
  return (grad)
}




generate_data = function(n,d, Sigma, sigma=1, x_true, type = "Logistic"){
  
  
  
  logis=matrix(0,d+1,n)
  for(i in 1:n){
    a = rmvnorm(1, sigma = Sigma)
    p = 1/c(1+exp(-a%*%x_true)) 
    b = 2*rbinom(1, 1, p)-1
    logis[,i]=c(b,a)
  }
  
  
  return(logis)
}




iterate=function(M,alpha){
  n=dim(M)[2]
  d=dim(M)[1]-1
  X=matrix(0,d,n)
  X[,1]=rep(0,d)
  for(i in 2:n){
    X[,i]=X[,i-1]-i^(-alpha)*grad_comp(X[,i-1],M[-1,i],M[1,i]
                                       ,type = "Logistic")
  }
  return(X)
}

#Hessian Matrix
heis_comp = function(x_old, a_new, b_new, type = "Linear"){
  if (type == "Linear"){
    heis = t(a_new)%*%a_new 
  }else{
    heis =  t(a_new)%*%a_new/c(1+exp(a_new%*%x_old))/c(1+exp(-a_new%*%x_old))
  }
  return (heis)
}
H=matrix(0,d,d)
for(i in 1:1000000){
  a = rmvnorm(1, sigma = Sigma)
  p = 1/c(1+exp(-a%*%xt)) 
  b = rbinom(1, 1, p)
  H=H+heis_comp(xt,a,b,type="Logistic")
}
H=H/1000000
H=solve(H)


#generate data
datalg=list()

for(j in 1:N){
  set.seed(j)
  datalg[[j]]=generate_data(n,d,Sigma,1,xt,type="Logistic")
}


ite=list()

for(j in 1:N){
  ite[[j]]=iterate(datalg[[j]],alpha)
}







#Compute weighted averaged SGD

#Polynomial Decay
for(j in 1:N){
  x=rep(0,d)
  for(i in 2:n){
    w[,j]=(1-(eta+1)/(eta+i))*w[,j]+(eta+1)/(eta+i)*ite[[j]][,i]
  }
}


#Suffix
for(j in 1:N){
  ws[,j]=apply(ite[[j]][,(n/2+1):n],1,mean)
}



#Optimal Weights
for(j in 1:N){
  ww[,j]=(1-2^(alpha-1))*ite[[j]][,1]+2^(alpha-1)*ite[[j]][,2]
  for(i in 3:n){
    ww[,j]=(i-1)/i*ww[,j]+
      (1-i^alpha)*ite[[j]][,i-1]/i+i^(alpha-1)*ite[[j]][,i]
  }
}





compo=4
plot(density((w[compo,]-xt[compo])*
               sqrt(n*(2*eta+1))/(eta+1)/sqrt(H[compo,compo])), 
     main="",xlab="",ylab="",col="blue",
     xlim=c(-6,6),ylim=c(0,0.45))

lines(density((ws[compo,]-xt[compo])*
                sqrt(n*0.5)/sqrt(H[compo,compo])),col="red")
lines(density((ww[compo,]-xt[compo])*
                sqrt(n)/sqrt(H[compo,compo])),col="green")


curve(dnorm,add=TRUE,lwd=1.5)
lines(density((w[compo,]-xt[compo])*
                sqrt(n)/sqrt(H[compo,compo])),col="blue",lty=2)
lines(density((ws[compo,]-xt[compo])*
                sqrt(n)/sqrt(H[compo,compo])),col="red",lty=2)

legend("topleft", c("standard normal distribution",
                    "Polynomial-Decay (correctly standardized)",
                    "0.5-suffix (correctly standardized)",
                    "Optimal Weights (correctly standardized)",
                    "Polynomial-Decay (unstandardized)",
                    "0.5-suffix (unstandardized)"),
       lty=c(1,1,1,1,2,2),col=c("black","blue","red","green","blue","red"),
       lwd=c(1.5,1,1,1,1,1),cex=0.43)


