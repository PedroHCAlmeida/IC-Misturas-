library(mixsmsn)

dMix.SN = function(Y, X, beta, sigma){
  total = mapply(function(X, y, sigma, beta) dsn(y, X%*%beta, sigma), 
                 beta = beta, sigma = sigma, MoreArgs = list(y = Y, X = X))
  return(rowSums(total))
}
.S3method("dMix", "SN", dMix.SN)

vero.SN = function(Y, X, params){
  sum(log(dMix.SN(Y, X, params$beta, params$sigma)))
}
.S3method("vero", "SN", vero.SN)

estimaTeta.SN = function(Y, X, n, p){
  
  Xl = t(X)
  beta = solve(Xl%*%X)%*%(Xl%*%Y)
  
  sigma = sqrt(sum((Y - (X%*%beta))^2)/(n - p))
  
  c(beta = beta, sigma = sigma)
  
  parte1Beta <- solve(t(X)%*%diag(U$Z)%*%X)
  parte2Beta <- t(X)%*%(U$Z*y - params$S2[,j]*params$delta)
  
  beta <- parte1Beta%*%parte2Beta
}
.S3method("estimaTeta", "SN", estimaTeta.SN)


