MoENormal = function(x){
  class(x) = c("MoENormal", "MixNormal", class(x))
  return(x)
}

dMix.MoENormal = function(y, X, beta, sigma, P){
  total = lapply(1:nrow(beta),
                 function(j) P[,j]*dnorm(y, X%*%beta[j,], sigma[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MoENormal", dMix.MoENormal)

vero.MoENormal = function(Y, X, params){
  sum(log(dMix(Y, X, 
               params$params[, startsWith(colnames(params$params), "beta")], 
               params$params[,"sigma"], params$P)))
}
.S3method("vero", "MoENormal", vero.MoENormal)

estimaTeta.MoENormal = function(Y, X, Z, R, alpha, P){
  
  beta = solve(t(X)%*%diag(Z)%*%X)%*%(t(X)%*%(Z*y))
  sigma = sqrt(sum(Z*(Y - (X%*%beta))^2)/sum(Z))
  alphaNovo = alpha + 4*solve(t(R)%*%R)%*%(t(R)%*%(Z - P))
  
  list(params = c(beta = beta, sigma = sigma, alpha = alphaNovo))
}
.S3method("estimaTeta", "MoENormal", estimaTeta.MoENormal)

estimaMedia.MoENormal = function(X, params){
  apply(params[, startsWith(colnames(params), "beta")], 1,
        function(betaI) estimaMedia.Normal(X = X, beta = betaI))
}
.S3method("estimaMedia", "MoENormal", estimaMedia.MoENormal)

chuteInicial.MoENormal = function(Y, X, n, p){
  
  k = ncol(R)
  dados = cbind(Y, X[, -1])
  
  grupos = tryCatch({
    if(initGrupo == "Aleatorio") sample(1:g, n, replace = T)
    else if(initGrupo == "KMeans") kmeans(dados, centers = g)$cluster
    else stop()
  },
  error = function(err) kmeans(dados, centers = g)$cluster)
  
  dadosGrupos = lapply(list("X" = X, "Y" = y), 
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]), 
                       grupos = grupos)
  
  params = do.call(rbind, mapply(estimaTeta.Normal, 
                                 dadosGrupos$Y, 
                                 dadosGrupos$X, 
                                 SIMPLIFY = F,
                                 MoreArgs = list(n = n, p = p)))
  
  P = matrix(rep(c(prop.table(table(grupos))), n), byrow = T, ncol = g)
  alpha = matrix(c(rep(0, (g-1)*k), rep(NA, k)), nrow = g, ncol = k, byrow = T)
  #alpha = matrix(c(nnet::multinom(grupos ~ R[,-1])$wts[-1], rep(NA, p)), nrow = g, ncol = k, byrow = T)
  colnames(alpha) = paste0("alpha", 1:k)
  params = cbind(params, alpha = alpha)
  
  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoENormal", chuteInicial.MoENormal)

etapaE.MoENormal = function(y, X, params, medias, ...){
  
  Z = sapply(seq_len(g), 
             function(j) params$P[,j]*dnorm(y, medias[,j], params$params[j,"sigma"]))
  
  Z = t(sapply(1:n, function(i) Z[i,]/sum(Z[i,])))
  return(list(Z = Z))
}
.S3method("etapaE", "MoENormal", etapaE.MoENormal)

etapaM.MoENormal = function(Y, X, U, params){
  
  paramsNovo = do.call(
    rbind, 
    sapply(1:g, 
          function(j) estimaTeta.MoENormal(Y = Y, X = X, Z = U$Z[,j], R = R,
                                           alpha = params$params[j,startsWith(colnames(params$params), "alpha")],
                                           P = params$P[,j])))

  P = matrizP(t(paramsNovo)[startsWith(colnames(paramsNovo), "alpha"), 1:(g-1)], R)

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MoENormal", etapaM.MoENormal)

predictReg.MoENormal = function(reg, x, r, type){
  
  source("scriptsEM/classes.R")
  source("scriptsEM/funcoes.R")
  source("scriptsEM/MetodosNormal.R")
  source("scriptsEM/MetodosMoENormal.R")
  source("scriptsEM/auxFuncs.R")
  
  n = nrow(x)
  
  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(x)), r)
  
  alpha = reg$Parâmetros %>%
    t() %>%
    as_tibble() %>%
    dplyr::select(starts_with("alpha")) %>%
    t()
  
  P = matrizP(alpha[,-ncol(alpha)], R)
  medias = estimaMedia.MoENormal(X, reg$params$params)
  if(type == 1) y = apply(P*medias, 1, sum)
  else{
    grupos = apply(P, 1, which.max)
    y = sapply(1:n, function(i) medias[i, grupos[i]])
  }
  return(y)
}
.S3method("predictReg", "MoENormal", predictReg.MoENormal)







