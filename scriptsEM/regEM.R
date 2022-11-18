library(mixsmsn)
library(codetools)

regEM = function(y, x, ..., tol = 1E-10, family = "MoENormal", 
                 grupoReal = NULL, max_iter = 1000, min_iter = 5){
  
  args = list(...)
  env = new.env(parent = environment())

  list2env(args, envir = env)
  
  setwd("~/ufjf/IC")
  source("scriptsEM/classes.R", local = env)
  source("scriptsEM/funcoes.R", local = env)
  source("scriptsEM/MetodosNormal.R", local = env)
  source("scriptsEM/MetodosMixNormal.R", local = env)
  source("scriptsEM/MetodosMoENormal.R", local = env)
  source("scriptsEM/resultadosEM.R", local = env)
  source("scriptsEM/auxFuncs.R", local = env)
  
  evalq({
    n = length(y)
    
    X = cbind(rep(1, n), x)
    p = ncol(X)
    try({
      R = cbind(rep(1, n), r)
      #list2env(list(R = R), envir = environment())
    })
      
    y = eval(parse(text = family))(y)
    X = eval(parse(text = family))(X)
      
    paramsAtual = chuteInicial(y, X, n, p)
    #paramsAtual = chuteInicial(Y = y, X = X, n = n, p = p, R = R, g = g)
    medias = estimaMedia(X, paramsAtual$params)
    crit = 1
    it = 0
    
    while(((crit > tol) & (it < max_iter)) | (it < min_iter)){
      # Calculando Verossimilhança
      vero0 = vero(y, X, paramsAtual)
      print(vero0)
      # Etapa E
      U = etapaE(y, X, paramsAtual, medias)
      #U = etapaE.MoENormal(Y = y, params = paramsAtual, medias = medias)
      
      # Etapa M
      paramsNovo = etapaM(y, X, U, paramsAtual)
      #paramsNovo = etapaM.MoENormal(Y = y, X = X, U = U, g = g, R = R, params = paramsAtual)
    
      # Estimando valores esperados
      medias = estimaMedia(X, paramsAtual$params)
        
      paramsAtual = paramsNovo
        
      # Calculando critério
      crit = abs(vero(y, X, paramsAtual)/vero0 - 1)
      print(crit)
      it = it+1
    }
      
    gruposEM = apply(U$Z, 1, which.max)
    
    if(!is.null(grupoReal)){
      
      ordemReg = order(table(gruposEM))
      ordemReal = order(table(grupoReal))
      
      tabela = table(grupoReal, replace(gruposEM, ordemReg, ordemReal),
                     dnn = list('Real', 'Modelo'))
      
      #tabela = table(grupoReal, gruposEM, dnn = list('Real', 'Modelo'))
      }
    else tabela = NA
      
    rownames(paramsNovo$params) = 1:nrow(paramsNovo$params)
      
    metricas = lapply(1:g, function(j) calculaMetricas(y[gruposEM == j], 
                                                       medias[gruposEM == j, j]))
    
    resultados = list(
      `Iterações` = it,
      l = vero(y, X, paramsNovo),
      `Parâmetros` = t(paramsNovo$params),
      U = U,
      tabela = tabela,
      medias = medias,
      metricas = metricas,
      params = paramsNovo,
      metricasTotais = colSums(do.call(rbind,
                                       lapply(metricas,
                                              function(x) sapply(2:4, function(i) x$n*x[[i]]/n)))),
      residuos = lapply(1:g, function(j) medias[gruposEM == j, j] - y[gruposEM == j])
    )
    class(resultados) = c("resultadosEM", family)
  },
  envir = env)
  
  #rm(list = ls(envir = env), envir = env)
  
  return(env$resultados)
}

