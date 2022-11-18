library(tidyverse)
library(AppliedPredictiveModeling)
source("scriptsEM/regEM.R")

 # Função auxiliar
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}


# AED

data(abalone)
dim(abalone)
head(abalone)

abalone$Age = abalone$Rings +1.5
abalone = subset(abalone, select = -Rings)

abalone = abalone %>%
  mutate(age_box_cox = box_cox(Age, 0.3))

abalone %>%
  ggplot(aes(box_cox(abalone$Age, 0.3), fill = Type)) +
  geom_histogram(bins = 15, alpha=0.6, position = 'identity') +
  theme_minimal()

abalone %>%
  filter(Type != "F") %>%
  ggplot(aes(Age, fill = Type)) +
  geom_histogram(bins = 15, alpha=0.6, position = 'identity') +
  theme_minimal()

abalone %>%
  filter(Type != "M") %>%
  ggplot(aes(Age, fill = Type)) +
  geom_histogram(bins = 15, alpha=0.6, position = 'identity') +
  theme_minimal()

abalone %>%
  filter(Type != "I") %>%
  ggplot(aes(Age, fill = Type)) +
  geom_histogram(bins = 15, alpha=0.6, position = 'identity') +
  theme_minimal()

abalone %>%
  pivot_longer(-Type) %>%
  ggplot(aes(value, fill = Type)) +
  facet_wrap(~name, scales = 'free') +
  geom_histogram(bins = 10, alpha=0.6, position = 'identity') +
  theme_minimal()

abalone %>%
  pivot_longer(-Type) %>%
  ggplot(aes(value)) +
  facet_wrap(~name, scales = 'free') +
  geom_histogram(bins = 10) +
  theme_minimal()

abalone %>%
  ggplot(aes(LongestShell, Age, colour = Type)) +
  geom_point() +
  theme_minimal()

help(abalone)

# Regressão 

ml = lm(Age ~ ., abalone)
summary(ml)

# Misturas 

estimativasMoEAbalone = regEM(y = abalone$age_box_cox,
      x = abalone %>% select(-Age, -age_box_cox, -Type, -LongestShell) %>% select_if(is.numeric) %>% as.matrix(),
      r = abalone %>% select(-Age, -age_box_cox, -Type) %>% select_if(is.numeric) %>% as.matrix(),
      g = 3, family = "MoENormal", tol = 1E-6, grupoReal = abalone$Type, initGrupo = "KMeans")

estimativasMoEAbalone
estimativasMoEAbalone$l
estimativasMoEAbalone$tabela
estimativasMoEAbalone$metricas
estimativasMoEAbalone$metricasTotais

estimativasMixAbalone = regEM(y = abalone$age_box_cox,
                              x = abalone %>% select(-Age, -age_box_cox, -Type, -LongestShell) %>% select_if(is.numeric) %>% as.matrix(),
                              g = 3, family = "MixNormal", tol = 1E-6, grupoReal = abalone$Type, initGrupo = "KMeans")

estimativasMixAbalone
estimativasMixAbalone$l
estimativasMixAbalone$tabela
estimativasMixAbalone$metricas
estimativasMixAbalone$metricasTotais

# Testando outro R

estimativasMoEAbaloneTeste = regEM(y = abalone$Age,
                                   x = abalone %>% select(-Age, -Type, -LongestShell) %>% select_if(is.numeric) %>% as.matrix(),
                                   r = abalone %>% select(LongestShell) %>% select_if(is.numeric) %>% as.matrix(),
                                   g = 3, family = "MoENormal", tol = 1E-6, grupoReal = abalone$Type, initGrupo = "KMeans")

estimativasMoEAbaloneTeste
estimativasMoEAbaloneTeste$l
estimativasMoEAbaloneTeste$tabela
estimativasMoEAbaloneTeste$metricas
estimativasMoEAbaloneTeste$metricasTotais

# Análise Resíduos

dadosAbaloneMetodos = data.frame(do.call(cbind.fill, estimativasMixAbalone$residuos)) %>%
  pivot_longer(everything()) %>% 
  filter(!is.na(value)) %>%
  mutate(Metodo = "Mix") %>%
  bind_rows(
    data.frame(do.call(cbind.fill, estimativasMoEAbalone$residuos)) %>%
      pivot_longer(everything()) %>% 
      filter(!is.na(value)) %>%
      mutate(Metodo = "MoE")) %>%
  rename(Grupo = name) %>%
  mutate(Grupo = str_replace(Grupo, "X", ""))

dadosAbaloneMetodos %>%
  ggplot(aes(value, ..density.., fill = Grupo)) +
  facet_grid(vars(Grupo), vars(Metodo)) +
  geom_histogram() +
  theme_minimal()

dadosAbaloneMetodos %>%
  ggplot(aes(Grupo, value, fill = Grupo)) +
  facet_wrap(~Metodo, nrow = 2) +
  geom_boxplot() +
  theme_minimal()
