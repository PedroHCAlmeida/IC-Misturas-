MixNormal = function(x){
  class(x) = c("MixNormal", "Normal", class(x))
  return(x)
}

dMix.MixNormal = function(y, X, beta, sigma, P){
  total = lapply(1:nrow(beta),
                 function(j) P[j]*dnorm(y, X%*%beta[j,], sigma[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MixNormal", dMix.MixNormal)

vero.MixNormal = function(Y, X, params){
  sum(log(dMix.MixNormal(Y, X, 
                         params$params[, startsWith(colnames(params$params), "beta")], 
                         params$params[,"sigma"], params$P)))
}
.S3method("vero", "MixNormal", vero.MixNormal)

estimaTeta.MixNormal = function(Y, X, Z){
  
  beta = solve(t(X)%*%diag(Z)%*%X)%*%(t(X)%*%(Z*y))
  sigma = sqrt(sum(Z*(Y - (X%*%beta))^2)/sum(Z))
  
  c(beta = beta, sigma = sigma)
}
.S3method("estimaTeta", "MixNormal", estimaTeta.MixNormal)

estimaMedia.MixNormal = function(X, params){
  apply(params[, startsWith(colnames(params), "beta")], 1,
        function(betaI) estimaMedia.Normal(X = X, beta = betaI))}
.S3method("estimaMedia", "MixNormal", estimaMedia.MixNormal)

chuteInicial.MixNormal = function(Y, X, n, p){
 
  dados = cbind(Y, X[, -1])

  grupos = tryCatch({
    if(initGrupo == "Aleatorio") sample(1:g, n, replace = T)
    else if(initGrupo == "Kmeans") kmeans(dados, centers = g)$cluster
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
  
  P = prop.table(table(grupos))
  
  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixNormal", chuteInicial.MixNormal)

etapaE.MixNormal = function(y, X, params, medias, ...){
  
  Z = sapply(seq_len(g), 
              function(j) params$P[j]*dnorm(y, medias[,j], params$params[j,"sigma"]))
  
  Z = t(sapply(1:n, function(i) Z[i,]/sum(Z[i,])))
  
  return(list(Z = Z))
}
.S3method("etapaE", "MixNormal", etapaE.MixNormal)

etapaM.MixNormal = function(Y, X, U, params){
  
  paramsNovo = do.call(
    rbind, 
    lapply(1:g, 
           function(j) estimaTeta.MixNormal(Y = Y, X = X, Z = U$Z[,j])))
  
  
  P = colMeans(U$Z)
  
  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MixNormal", etapaM.MixNormal)

