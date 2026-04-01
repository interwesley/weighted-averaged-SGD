##########Squared Loss#######################
library(mvtnorm) 
N=400
xt=c(1,-2,0,0,4)
d=5
n=10000
Sigma=diag(1,5)
alpha=0.8
burn=50
index=c(1950,3950,5950,7950,9950)



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



#ASGD
mseasgd=matrix(0,n,N)
solmse=matrix(0,d,N)
for(j in 1:N){
  mseasgd[1,j]=mean((solmse[,j]-xt)^2)
  for(i in 2:(n-burn)){
    solmse[,j]=(1-1/i)*solmse[,j]+1/i*ite[[j]][,i+burn]
    mseasgd[i,j]=mean((solmse[,j]-xt)^2)
  }
}
ASGD=apply(mseasgd[1:(n-burn),],1,mean)
Asd=apply(mseasgd[1:(n-burn),],1,sd)






#Optimal Weights
mseow=matrix(0,n,N)
solmse=matrix(0,d,N)


  
  for(j in 1:N){
    mseow[1,j]=mean((solmse[,j]-xt)^2)
    solmse[,j]=2^(alpha-1)*ite[[j]][,2]
    mseow[2,j]=mean((solmse[,j]-xt)^2)
    for(i in 3:(n-burn)){
      solmse[,j]=(i-1)/i*solmse[,j]+
        (1-i^alpha)*ite[[j]][,i+burn-1]/i+i^(alpha-1)*ite[[j]][,i+burn]
      mseow[i,j]=mean((solmse[,j]-xt)^2)
    }
  }
OW=apply(mseow[1:(n-burn),],1,mean)
OWsd=apply(mseow[1:(n-burn),],1,sd)





#Suffix
msesuf=matrix(0,n,N)
solmse=matrix(0,d,N)

  
  for(j in 1:N){
    msesuf[1,j]=mean((solmse[,j]-xt)^2)
    msesuf[2,j]=mean((ite[[j]][,2]-xt)^2)
    msesuf[3,j]=mean((ite[[j]][,3]-xt)^2)
    solmse[,j]=apply(ite[[j]][,3:4],1,mean)
    msesuf[4,j]=mean((solmse[,j]-xt)^2)
    for(i in 5:(n-burn)){
      if(i%%2==0){
        solmse[,j]=solmse[,j]*(i-2)/i+ite[[j]][,i+burn]*2/i
      }else{
        solmse[,j]=solmse[,j]-2/(i-1)*ite[[j]][,(i/2+0.5+burn/2)]+
          2/(i-1)*ite[[j]][,i+burn]
      }
      msesuf[i,j]=mean((solmse[,j]-xt)^2)
    }
  }


SUFFIX=apply(msesuf[1:(n-burn),],1,mean)
Ssd=apply(msesuf[1:(n-burn),],1,sd)



#Polynoial Decay

eta=3
msepoly=matrix(0,n,N)
solmse=matrix(0,d,N)


for(j in 1:N){
  msepoly[1,j]=mean((solmse[,j]-xt)^2)
  for(i in 2:(n-burn)){
    solmse[,j]=(1-(eta+1)/(eta+i))*solmse[,j]+(eta+1)/(eta+i)*ite[[j]][,i+burn]
    msepoly[i,j]=mean((solmse[,j]-xt)^2)
  }
}
PD=apply(msepoly[1:(n-burn),],1,mean)
Psd=apply(msepoly[1:(n-burn),],1,sd)





#Last
msel=matrix(0,n,N)

for(j in 1:N){
  msel[,j]=apply((ite[[j]][,1:n]-xt)^2,2,mean)
}



LAST=apply(msel[(burn+1):n,],1,mean)
Lsd=apply(msel[(burn+1):n,],1,sd)



OW[index]
ASGD[index]
SUFFIX[index]
PD[index]
LAST[index]

OWsd[index]
Asd[index]
Ssd[index]
Psd[index]
Lsd[index]





#######Mean Estimation###########

MEdata=function(n,N,alpha){
X=matrix(0,n,N)
  for(j in 1:N){
    set.seed(1000+j)
    X[1,j]=rnorm(1)
    for(i in 2:n){
      X[i,j]=(1-i^(-alpha))*X[i-1,j]+i^(-alpha)*rnorm(1)
      
    }
  }
return(X)
}

n=1600
N=400
alpha=0.505
X=MEdata(n,N,alpha)

index=c(400,800,1200,1600)

#Optimal Weight

mseow=matrix(0,n,N)
solmse=rep(0,N)

system.time({
  
  for(j in 1:N){
    mseow[1,j]=X[1,j]^2
    solmse[j]=(1-2^(alpha-1))*X[1,j]+2^(alpha-1)*X[2,j]
    mseow[2,j]=solmse[j]^2
    for(i in 3:n){
      solmse[j]=(i-1)/i*solmse[j]+(1-i^alpha)*X[i-1,j]/i+i^(alpha-1)*X[i,j]
      mseow[i,j]=solmse[j]^2
    }
  }
})

OW=apply(mseow,1,mean)
OWsd=apply(mseow,1,sd)









#ASGD

mseasgd=matrix(0,n,N)
solmse=rep(0,N)

system.time({
  
  for(j in 1:N){
    mseasgd[1,j]=X[1,j]^2
    solmse[j]=X[1,j]
    for(i in 2:n){
      solmse[j]=(1-1/i)*solmse[j]+1/i*X[i,j]
      mseasgd[i,j]=solmse[j]^2
    }
  }
})

ASGD=apply(mseasgd,1,mean)
Asd=apply(mseasgd,1,sd)










#Poly Decay
eta=3
msepd=matrix(0,n,N)
solmse=rep(0,N)

system.time({
  
  for(j in 1:N){
    msepd[1,j]=X[1,j]^2
    solmse[j]=X[1,j]
    for(i in 2:n){
      solmse[j]=(1-(eta+1)/(eta+i))*solmse[j]+(eta+1)/(eta+i)*X[i,j]
      msepd[i,j]=solmse[j]^2
    }
  }
})

PD=apply(msepd,1,mean)
Psd=apply(msepd,1,sd)



#suffix

msesuf=matrix(0,n,N)
solmse=rep(0,N)

system.time({
  
  for(j in 1:N){
    msesuf[1,j]=X[1,j]^2
    msesuf[2,j]=X[2,j]^2
    msesuf[3,j]=X[3,j]^2
    solmse[j]=(X[3,j]+X[4,j])/2
    msesuf[4,j]=solmse[j]^2
    for(i in 5:n){
      if(i%%2==0){
        solmse[j]=solmse[j]*(i-2)/i+X[i,j]*2/i
      }else{
        solmse[j]=solmse[j]-2/(i-1)*X[(i/2+0.5),j]+
          2/(i-1)*X[i,j]
      }
      msesuf[i,j]=solmse[j]^2
    }
  }
})

SUFFIX=apply(msesuf,1,mean)
Ssd=apply(msesuf,1,sd)


#Last

LAST=apply(X[1:n,]^2,1,mean)
Lsd=apply(X[1:n,]^2,1,sd)


OW[index]
ASGD[index]
SUFFIX[index]
PD[index]
LAST[index]

OWsd[index]
Asd[index]
Ssd[index]
Psd[index]
Lsd[index]




#change the value of alpha___________

alpha=0.8
X=MEdata(n,N,alpha)



#Optimal Weight

mseow=matrix(0,n,N)
solmse=rep(0,N)

system.time({
  
  for(j in 1:N){
    mseow[1,j]=X[1,j]^2
    solmse[j]=(1-2^(alpha-1))*X[1,j]+2^(alpha-1)*X[2,j]
    mseow[2,j]=solmse[j]^2
    for(i in 3:n){
      solmse[j]=(i-1)/i*solmse[j]+(1-i^alpha)*X[i-1,j]/i+i^(alpha-1)*X[i,j]
      mseow[i,j]=solmse[j]^2
    }
  }
})

OW=apply(mseow,1,mean)
OWsd=apply(mseow,1,sd)








#ASGD

mseasgd=matrix(0,n,N)
solmse=rep(0,N)

system.time({
  
  for(j in 1:N){
    mseasgd[1,j]=X[1,j]^2
    solmse[j]=X[1,j]
    for(i in 2:n){
      solmse[j]=(1-1/i)*solmse[j]+1/i*X[i,j]
      mseasgd[i,j]=solmse[j]^2
    }
  }
})

ASGD=apply(mseasgd,1,mean)
Asd=apply(mseasgd,1,sd)









#Poly Decay
eta=3
msepd=matrix(0,n,N)
solmse=rep(0,N)

system.time({
  
  for(j in 1:N){
    msepd[1,j]=X[1,j]^2
    solmse[j]=X[1,j]
    for(i in 2:n){
      solmse[j]=(1-(eta+1)/(eta+i))*solmse[j]+(eta+1)/(eta+i)*X[i,j]
      msepd[i,j]=solmse[j]^2
    }
  }
})

PD=apply(msepd,1,mean)
Psd=apply(msepd,1,sd)









#suffix

msesuf=matrix(0,n,N)
solmse=rep(0,N)

system.time({
  
  for(j in 1:N){
    msesuf[1,j]=X[1,j]^2
    msesuf[2,j]=X[2,j]^2
    msesuf[3,j]=X[3,j]^2
    solmse[j]=(X[3,j]+X[4,j])/2
    msesuf[4,j]=solmse[j]^2
    for(i in 5:n){
      if(i%%2==0){
        solmse[j]=solmse[j]*(i-2)/i+X[i,j]*2/i
      }else{
        solmse[j]=solmse[j]-2/(i-1)*X[(i/2+0.5),j]+
          2/(i-1)*X[i,j]
      }
      msesuf[i,j]=solmse[j]^2
    }
  }
})

SUFFIX=apply(msesuf,1,mean)
Ssd=apply(msesuf,1,sd)




#Last

LAST=apply(X[1:n,]^2,1,mean)
Lsd=apply(X[1:n,]^2,1,sd)



OW[index]
ASGD[index]
SUFFIX[index]
PD[index]
LAST[index]

OWsd[index]
Asd[index]
Ssd[index]
Psd[index]
Lsd[index]

