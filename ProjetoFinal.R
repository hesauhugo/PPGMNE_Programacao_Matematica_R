## Projeto: Modelo Poisson com efeito aleatório --------------------------------
## Prof. Wagner Hugo Bonat LEG/UFPR --------------------------------------------
## Data: 18/06/2021 ------------------------------------------------------------
## Hesau,Ingrid

##---------------Código Melhorado-----------------------------------------------
## Fixando valores do parâmetros
install.packages("profvis")
require(profvis)

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

## Integrando
Tintegrando<-function(Tbi, Tbeta0, Tsigma2, Ty){dpois(Ty, exp(Tbeta0 + Tbi))*(1/((sqrt(2 * pi * Tsigma2))*exp((Tbi^2)/(2*Tsigma2))))}
## verossimilhança
Tmy_ll <- function(par, y){-sum(sapply(y,function(y) log(integrate(f = Tintegrando, lower = -Inf, upper = Inf, Tbeta0 = par[1], Tsigma2 = exp(par[2]), Ty = y)$value)))}
## resultado
Ttemp <- optim(par = c(log(8), log(0.3)), fn = Tmy_ll, y = Tyi, method = "Nelder-Mead")
Ttemp
##---------------Código Antigo-----------------------------------------------

## Fixando valores do parâmetros
sigma2 <- 0.5
beta0 <- log(10)

# Tamanho da amostra
set.seed(123)

n = 100

## Simular o efeito aleatório b_i
bi <- rnorm(n = n, mean = 0, sd = sqrt(sigma2))

## Preditor linear
lambdai <- exp(beta0 + bi)

## Simulando os dados Yi's
yi <- rpois(n = n, lambda = lambdai)
hist(yi)

## FIM DA PARTE DE SIMULAÇÃO ---------------------------------------------------

integrando <- function(bi, beta0, sigma2, y) {
  lambda <- exp(beta0 + bi)
  out <- dpois(y, lambda = lambda)*dnorm(bi, mean = 0, sd = sqrt(sigma2))
  return(out)
}

## Gráfico do integrando
b <- seq(-3, 3, l = 100)

plot(integrando(bi = b, beta0 = log(10), sigma2 = 0.5, y = yi[1]) ~ b, type = "l")
for(i in 2:100) {
  lines(b, integrando(bi = b, beta0 = log(10), sigma2 = 0.5, y = yi[i]))
}

## Resolvendo a integral
integral <- c()
for(i in 1:100) {
  integral[i] <- integrate(f = integrando, lower = -Inf, upper = Inf, 
                           beta0 = beta0, sigma2 = sigma2, y = yi[i])$value
}

ll <- sum(log(integral))
ll

## Criando uma função: log-verossimilhança

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

temp <- optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi,method = "Nelder-Mead")
temp
##----------------------BenchMark------------------------------------------------

 
 require(bench)
 require(rbenchmark)

benchmark("Professor" = optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi,method = "Nelder-Mead"),
           
           "Meu" = optim(par = c(log(8), log(0.3)), fn = Tmy_ll, y = Tyi, method = "Nelder-Mead"),
           replications = 10)


benchmark("Professor" = optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi,method = "BFGS"),
          
          "Meu" = optim(par = c(log(8), log(0.3)), fn = Tmy_ll, y = Tyi, method = "BFGS"),
          replications = 10)

benchmark("Professor" = optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi,method = "SANN"),
          
          "Meu" = optim(par = c(log(8), log(0.3)), fn = Tmy_ll, y = Tyi, method = "SANN"),
          replications = 10)

benchmark("Professor" = optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi,method = "CG"),
          
          "Meu" = optim(par = c(log(8), log(0.3)), fn = Tmy_ll, y = Tyi, method = "CG"),
          replications = 10)
