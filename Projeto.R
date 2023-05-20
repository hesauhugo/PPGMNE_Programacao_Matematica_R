## Projeto: Modelo Poisson com efeito aleatório --------------------------------
## Prof. Wagner Hugo Bonat LEG/UFPR --------------------------------------------
## Data: 18/06/2021 ------------------------------------------------------------

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

## FIM DA PARTE DE SIMULAÇÃO ---------------------------------------------------
integrando <- function(bi, beta0, sigma2, y) {dpois(y, lambda = exp(beta0 + bi)) * dnorm(bi, mean = 0, sd = sqrt(sigma2)}

## Criando uma função: log-verossimilhança

my_ll <- function(par, y) {
  integral <- c()
    for (i in 1:length(y)) {
      
    integral[i] <- integrate(f = integrando, lower = -Inf, upper = Inf, 
                             beta0 = par[1], sigma2 = exp(par[2]), y = y[i])$value
  }
  
  ll <- sum(log(integral))
    return(-ll)
}


my_ll(par = c(beta0, log(sigma2) ), y = yi)
 
temp <- optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi, method = "Nelder-Mead")
temp
system.time(temp <- optim(par = c(log(8), log(0.3)), fn = my_ll, y = yi, method = "Nelder-Mead"))
# temp2 <- optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi, 
#               method = "BFGS")
# temp2
# temp3 <- optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi, 
#                method = "SANN")
# 
# temp3 <- optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi, 
#                method = "CG")
# temp3
# 
system.time(temp <- optim(par = c(log(8), log(0.3)), fn = my_ll, y = yi,
               method = "Nelder-Mead"))



