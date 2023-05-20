## Projeto: Modelo Poisson com efeito aleatorio --------------------------------
## Michely Castro dos Santos Sena ------------------------------------------------------------
require (pracma)

#Constantes
n = 100 #N amostras
## Fixando valores do parâmetros
## bi ~ N(media, sigma^2)
media <- 0
sigma2 <- 0.5

# Tamanho da amostra
set.seed(123) #pseudo aleatório

## Simular o efeito aleatório b_i
bi <- rnorm(n = n, mean = media, sd = sqrt(sigma2))

#onde lambdaI =  e^{beta0 + bi} ## Preditor linear
beta0 <- log(10)
lambdaI <- exp(beta0 + bi)

## Simulando os dados Yi's
#Yi|bi ~ P(lambdaI)
yi <- rpois(n = n, lambda = lambdaI)
hist(yi)

## FIM DA PARTE DE SIMULACAO ---------------------------------------------------

## DECLARANDO FUNÇÕES


calcInt <- function(bi, beta0, sigma2, y){
  lambda <- exp(beta0 + bi)
  
  #c <- (lambda^y*exp(-lambda))/factorial(y)   
  #Não foi possível abrir a distrição de
  #Poisson, pois retorna números inválidos
  #print(c)
  #x <- dpois(y, lambda = lambda)
  #print(x)
  
  d <- sqrt(2 * 3.14 * sigma2)
  e <- (bi^2)/(2*sigma2)
  f <- exp(e)
  #resultado <- (c * (1/(d*f)))
  resultado <- dpois(y, lambda = lambda)*(1/(d*f))
  return(resultado)
}

### Criando uma função: log-verossimulhanança

my_ll <- function(par, y) {
  integralMi <- c()
  for(i in 1:length(y)) {
    integralMi[i] <- integrate(f = calcInt, lower = -Inf, upper = Inf,
                               beta0 = beta0, sigma2 = sigma2, y = yi[i])$value
  }
  ll <- sum(log(integralMi))
  return(-ll)
}



## Gráfico do integrando
b <- seq(-3, 3, l = n)
plot(calcInt(bi = b, beta0 = log(10), sigma2 = sigma2, y = yi[1]) ~ b, type = "l")
for(i in 2:n) {
  lines(b, calcInt(bi = b, beta0 = log(10), sigma2 = sigma2, y = yi[i]))
}


## Resolvendo a integral

integralMi <- c()
for(i in 1:n) {
  integralMi[i] <- integrate(f = calcInt, lower = -Inf, upper = Inf,
                             beta0 = beta0, sigma2 = sigma2, y = yi[i])$value
}

ll <- sum(log(integralMi))


system.time(temp <- optim(par = c(log(8), log(0.3) ), fn = my_ll, y = yi,method = "Nelder-Mead"))
