#########Plot Adaptive Weights for Expectile Regression###############

library(expectreg)
SGD_expectile = function(p, eta, a,polypara,n,root){
  
  xs_SGD = rep(0, n) 
  xs_ASGD = rep(0, n) 
  xs_optimalSGD = rep(0, n) 
  xs_polySGD = rep(0, n) 
  xs_rootSGD = rep(0, n) 
  xs_suffixSGD = rep(0, n) 
  x = x_bar = x_tilde=x_poly= x_root = x_suffic = 0
  
  #Since we focus on small finite-sample performance, we will not use burn in
  x = rnorm(1, 0, 0.5)
  x_root = x
  v=x
  m = S_0 = S_1 = 0
  for (rep in 1:n){
    y = rnorm(1)
    g = 2*(x - y)*abs(p - (y<x))
    
    
    x_old = x
    x_rootold=x_root
    
    g_old = 2*(x_rootold - y)*abs(p - (y<x_rootold))
    g_root = 2*(x_root - y)*abs(p - (y<x_root))
    
    v = g_root+ (rep-1)/rep*(v-g_old)
    
    x = x - eta*n**(-a)*g
    x_bar = (x + (rep-1)*x_bar)/rep
    x_tilde = ((rep-1)*x_tilde + (1 - rep**a)*x_old + rep**a*x)/rep
    x_poly=(1-(polypara+1)/(polypara+rep))*x_poly+(polypara+1)/(polypara+rep)*x
    
    
    if (rep >2^m){
      m=m+1
      S_0=S_1
      S_1=x
    }else{
      S_1=S_1+x
    }
    x_suffix=(S_0+S_1)/(rep- floor(2^(m-2)))
    
    x_root=x_root-root*v
    
    xs_SGD[rep] = x
    xs_ASGD[rep] =  x_bar
    xs_optimalSGD[rep] = x_tilde
    xs_polySGD[rep] = x_poly
    xs_rootSGD[rep]= x_root
    xs_suffixSGD[rep]=x_suffix
  }
  
  ans= list()
  ans$SGD = xs_SGD
  ans$ASGD = xs_ASGD
  ans$OSGD = xs_optimalSGD
  ans$PSGD = xs_polySGD
  ans$RSGD = xs_rootSGD
  ans$SSGD = xs_suffixSGD
  ans
}


experiment_weight = function(nrep, p,eta, a,n,root){
  x_true = enorm(p) 
  Sigma = matrix(0, nrow = n, ncol = n)
  for (rep in 1:nrep){
    one_trial = SGD_expectile(p, eta, a,3,n,root)
    one_trial$SGD 
    Sigma = Sigma + (one_trial$SGD - x_true)%*%t(one_trial$SGD - x_true)
  } 
  Sigma = Sigma/nrep
  Sigma_inv = solve(Sigma)
  colSums(Sigma_inv)/sum(Sigma_inv)
} 


polyweight=function(n,para){
  w=1
  for(i in 2:n){
    w=c((1-(para+1)/(para+i))*w,(para+1)/(para+i))
  }
  w
}




set.seed(123)

n = 50
nrep = 50000
root=0.0002
optimal_weight = experiment_weight(nrep, 0.8, 1, 0.505, n,root)
mean_optimal_weight =  c((1 + (1:(n-1))**0.505 -(2:n)**0.505)/n, n**(-0.495))

plot(1:n, optimal_weight,type = 'l', col = 'red',lty=4,lwd=2,
     xla="n",ylab="Weights")
lines(1:n, mean_optimal_weight, col = 'black')
lines(1:n, rep(1/n, n),col="pink")
lines(1:n,c(rep(0,n/2),rep(2/n,n/2)),col="green")
lines(1:n,polyweight(n,3),col="blue")


legend("top", legend = c('oracle','adaptive', 'ASGD','suffix','polynomial-decay'), 
       lwd = c(2,1,1,1,1),
       lty=c(4,1,1,1,1),col = c('red','black','pink','green','blue') )


