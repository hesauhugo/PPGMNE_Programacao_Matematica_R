## Projeto: Modelo Poisson com efeito aleatório --------------------------------
## Prof. Wagner Hugo Bonat LEG/UFPR --------------------------------------------
## Data: 18/06/2021 ------------------------------------------------------------
## Hesau,Ingrid

##---------------Código Melhorado-----------------------------------------------
## Fixando valores do parâmetros
install.packages("Rcpp")
install.packages("profvis")
install.packages("RcppArmadillo")
install.packages("std")

require(Rcpp)
require(profvis)
require(RcppArmadillo)
require(Matrix)
require(expm)
require(std)
require(rbenchmark)

sourceCpp("RcppPoisson.cpp")

Tsigma2 <- 0.5
Tbeta0 <- log(10)

#Tamanho da amostra
set.seed(123)
Tn = 100

## Simular o efeito aleatório b_i
Tbi <- rnorm(n = Tn, mean = 0, sd = sqrt(Tsigma2))

## lambda
Tlambdai <- exp(Tbeta0 + Tbi)

## Simulando os dados Yi's

Tyi <- rpois(n = Tn, lambda = Tlambdai)

#Integrandos
sum_cpp<-cppFunction(
  "double sum_cpp(NumericVector x) {
  const int n = x.size();
  double y  ; 
  for (int i=1; i < n; ++i) {
    y = y+x[i];
  }
  return y;
}")

Tintegrando<-function(Tbi, Tbeta0, Tsigma2, Ty){
  a= Exp_dpois(Ty,exp(Tbeta0+Tbi))
  b= Ext_dnorm(Tbi, mean = 0, sd = sqrt(Tsigma2))
  out <- exp(log(a) + log(b))
  return(out)
}

Tintegrando2<-function(Tbi, Tbeta0, Tsigma2, Ty){
  dpois(Ty, exp(Tbeta0 + Tbi))*(1/((sqrt(2 * pi * Tsigma2))*exp((Tbi^2)/(2*Tsigma2))))
}

Tintegrando3<-function(Tbi, Tbeta0, Tsigma2, Ty){
  Ext_Integrando_NumericVector(Tbi,Tbeta0,Tsigma2,Ty) * dnorm(Tbi, mean = 0, sd = sqrt(Tsigma2))
}

TintegrandoNaMao<-function(Tbi, Tbeta0, Tsigma2, Ty){
  lambda <- exp(Tbeta0+Tbi)
  x <- log(lambda^Ty)
  w <- Tbi^2
  v <- (2*Tsigma2)
  e <- -lambda - (w/v)
  f <- log(sqrt(2*pi*Tsigma2))
  z <- log(factorial(Ty))
  exp(x + e - z -f)
}

TintegrandoSapply<-function(Tbi, Tbeta0, Tsigma2, Ty){
  dpois(Ty, exp(Tbeta0 + Tbi))*(1/((sqrt(2 * pi * Tsigma2))*exp((Tbi^2)/(2*Tsigma2))))
}

integrando <- function(bi, beta0, sigma2, y) {
  lambda <- exp(beta0 + bi)
  out <- dpois(y, lambda = lambda)*dnorm(bi, mean = 0, sd = sqrt(sigma2))
  return(out)
}

Ext_Integrando_Completo(Tbi,Tbeta0,Tsigma2,Tyi[1])
Ext_Integrando(Tbi,Tbeta0,Tsigma2,Tyi[1])

integrando(Tbi,Tbeta0,Tsigma2,Tyi[1])
Tintegrando(Tbi,Tbeta0,Tsigma2,Tyi[1])
Tintegrando2(Tbi,Tbeta0,Tsigma2,Tyi[1])
Tintegrando3(Tbi,Tbeta0,Tsigma2,Tyi[1])
TintegrandoNaMao(Tbi,Tbeta0,Tsigma2,Tyi[1])

my_ll <- function(par, y) {
  integral <- c()
  for(i in 1:length(y)) {
    integral[i] <- integrate(f = integrando, lower = -Inf, upper = Inf, 
                             beta0 = par[1], sigma2 = exp(par[2]), y = y[i])$value
  }
  #print(c(round(par, 2)))
  ll <- sum(log(integral))
  return(-ll)
}

Tmy_ll <- function(par, y){
  
-sum(sapply(y,function(y) log(integrate(f = Tintegrando, lower = -Inf, upper = Inf, Tbeta0 = par[1], Tsigma2 = exp(par[2]), Ty = y)$value)))

}

Tmy_ll2 <- function(par, y){
  integral<- c()
  for(i in 1:length(y)){
    integral[i] <- integrate(f =  Tintegrando2, lower = -Inf, upper = Inf, Tbeta0 = par[1], Tsigma2 = exp(par[2]), Ty = y[i])$value
  }
  -sum(log(integral))
}


Tmy_ll3 <- function(par, y){
  integral<- c()
  for(i in 1:length(y)){
    integral[i] <- integrate(f =  Ext_Integrando_NumericVector, lower = -Inf, upper = Inf, Tbeta0 = par[1], Tsigma2 = exp(par[2]), Ty = y[i])$value
  }
  -sum_cpp(log(integral))
}

Tmy_ll4 <- function(par, y){
  integral<- c()
  for(i in 1:length(y)){
    integral[i] <- integrate(f =  TintegrandoNaMao, lower = 0, upper = Inf, Tbeta0 = par[1], Tsigma2 = exp(par[2]), Ty = y[i])$value
  }
  -sum_cpp(log(integral))
}

Tmy_ll5 <- function(par, y){
  
  -sum(sapply(-Inf:Inf,function(y) log(sum(Tintegrando2(par[1],exp(par[2])))), Ty = y))
  
}

Tmy_ll_sapply <- function(par, y){
  -sum(sapply(y,function(y) log(integrate(f = Tintegrando, lower = -Inf, upper = Inf, Tbeta0 = par[1], Tsigma2 = exp(par[2]), Ty = y)$value)))
}

temp <- optim(par = c(log(8), log(0.3) ), fn = my_ll, y = Tyi,method = "Nelder-Mead")
temp
Ttemp <- optim(par = c(log(8), log(0.3)), fn = Tmy_ll, y = Tyi, method = "Nelder-Mead")
Ttemp
Ttemp2 <- optim(par = c(log(8), log(0.3)), fn = Tmy_ll2, y = Tyi, method = "Nelder-Mead")
Ttemp2
Ttemp3 <- optim(par = c(log(8), log(0.3)), fn = Tmy_ll3, y = Tyi, method = "Nelder-Mead")
Ttemp3
Ttemp4 <- optim(par = c(log(8), log(0.3)), fn = Tmy_ll4, y = Tyi, method = "Nelder-Mead")
Ttemp4
Ttemp5 <- optim(par = c(log(8), log(0.3)), fn = Tmy_ll5, y = Tyi, method = "Nelder-Mead")
Ttemp5
Ttemp6 <- optim(par = c(log(8), log(0.3)), fn = Tmy_ll_sapply, y = Tyi, method = "Nelder-Mead")
Ttemp6


benchmark("Tmy_ll" = optim(par = c(log(8), log(0.3)), fn = Tmy_ll, y = Tyi, method = "Nelder-Mead"),
          
          "Tmy_ll2" =optim(par = c(log(8), log(0.3)), fn = Tmy_ll2, y = Tyi, method = "Nelder-Mead"),
          
          "Tmy_ll_sapply" =optim(par = c(log(8), log(0.3)), fn = Tmy_ll_sapply, y = Tyi, method = "Nelder-Mead"),
          
          "my_ll" =optim(par = c(log(8), log(0.3)), fn = Tmy_ll_sapply, y = Tyi, method = "Nelder-Mead"),
          
          replications = 5)




