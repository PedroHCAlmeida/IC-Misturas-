library(tidyverse)
source("scriptsEM/regEM.R")
source("scriptsEM/auxfuncs.R")

# Função auxiliar
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

precision = function(x){
  p = sum(diag(x))/sum(x)
  if(p<0.5){
    p = sum(diag(x[,c(2,1)]))/sum(x)
  }
  p
}

# Simulando dados

set.seed(123)
n = 3000
X = cbind(rep(1, n), rnorm(n, 5, 3), rnorm(n, 10, 5))
x = X[,-1]

R = cbind(rep(1, n), runif(n, -2, 1), runif(n, -1, 1))
r = R[,-1]

beta01 <- c(5, 6, -1)
alpha01 <- c(0.7, 1, 2)
sigma01 <- 100
mu01 <- X%*%beta01

beta02 <- c(-10, 1, 4)
sigma02 <- 30
mu02 <- X%*%beta02

alpha <- matrix(alpha01, byrow = T, nrow = 3)
P <- matrizP(alpha, R)
table(apply(P, 1, which.max))

obs = do.call(rbind, 
              lapply(1:n,
                     function(i){
                       arg1 <- list(mu = mu01[i,], sigma2 = sigma01, shape = 0)
                       arg2 <- list(mu = mu02[i,], sigma2 = sigma02, shape = 0)
                       
                       obs <- rmix(1, pii = P[i,], family = 'Normal', 
                                   arg = list(arg1, arg2), cluster = T)
                       
                       return(obs)
                     }))

rm(beta01, alpha01, sigma01, mu01, beta02, sigma02, mu02, alpha)
y = unlist(obs[,1])
grupo = unlist(obs[,2]) 

table(grupo)
ng = 2

# Analisando Grupos

data.frame(grupo = as.factor(grupo), y = y) %>%
  ggplot(aes(x = y, fill = grupo)) +
  geom_histogram(bins = 15, alpha=0.6, position = 'identity') +
  theme_minimal()

# Regressão usual

ml = lm(y ~ ., data.frame(cbind(x, y = y)))
summary(ml)

# Misturas

regMoEN = regEM(y, x, r = r, g = ng, tol = 1E-6,
                grupoReal = grupo,
                family = "MoENormal", initGrupo = "KMeans", min_iter = 1)

regMoEN$tabela
precision(regMoEN$tabela)
regMoEN$l
regMoEN$metricas
regMoEN$Parâmetros

regMixN = regEM(y, x, g = ng, tol = 1E-6, grupoReal = grupo, family = "MixNormal")
regMixN$tabela
regMixN$metricas
precision(regMixN$tabela)
regMixN$Parâmetros

# Análise resíduos

dadosMetodos = data.frame(do.call(cbind.fill, regMoEN$residuos)) %>%
  pivot_longer(everything()) %>% 
  filter(!is.na(value)) %>%
  mutate(Metodo = "MoE") %>%
  bind_rows(
    data.frame(do.call(cbind.fill, regMixN$residuos)) %>%
      pivot_longer(everything()) %>% 
      filter(!is.na(value)) %>%
      mutate(Metodo = "Mix")) %>%
  rename(Grupo = name) %>%
  mutate(Grupo = str_replace(Grupo, "X", ""))

dadosMetodos %>%
  ggplot(aes(value, ..density.., fill = Grupo)) +
  facet_grid(vars(Grupo), vars(Metodo)) +
  geom_histogram() +
  theme_minimal()

dadosMetodos %>%
  ggplot(aes(Grupo, value, fill = Grupo)) +
  facet_wrap(~Metodo, nrow = 2) +
  geom_boxplot() +
  theme_minimal()

#################################################

library(mixtools)

regmix = regmixEM(y, x)
reghhme = hmeEM(y, x)

################################################
