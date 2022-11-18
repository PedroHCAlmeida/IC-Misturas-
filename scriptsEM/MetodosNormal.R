estimaTeta.Normal = function(Y, X, n, p){
  
  Xl = t(X)
  beta = solve(Xl%*%X)%*%(Xl%*%Y)
  
  sigma = sqrt(sum((Y - (X%*%beta))^2)/(n - p))
  
  c(beta = beta, sigma = sigma)
}
.S3method("estimaTeta", "Normal", estimaTeta.Normal)

estimaMedia.Normal = function(X, beta){
  X%*%beta
}
.S3method("estimaMedia", "Normal", estimaMedia.Normal)

chuteInicial.Normal = function(Y, X){
  params = estimaTeta.Normal(Y, X)
  return(params)
  }
.S3method("chuteInicial", "Normal", chuteInicial.Normal)



# chuteInicial.MoENormal = function(Y, X, R, g){
#             
#   n = length(y)
#   k = ncol(R)
#   dados = cbind(y, X[, -1])
#             
#             
#   #grupos = kmeans(dados, centers = g)$cluster
#   grupos = sample(1:g, n, replace = T)
#             
#   dadosGrupos = lapply(list("X" = X, "Y" = Y), 
#                        function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]), 
#                        grupos = grupos
#   )
#             
#   beta = mapply(estimaBeta.Normal, dadosGrupos$Y, dadosGrupos$X, SIMPLIFY = F)
#   p = matrix(rep(as.vector(prop.table(sapply(dadosGrupos$X, nrow))), n), 
#              byrow = T, ncol = g)
#   sigma = mapply(estimaSigma.Normal, dadosGrupos$Y, dadosGrupos$X, beta, SIMPLIFY = F)
#   alpha = matrix(0, nrow = g-1, ncol = k)
#   
#   return(list(beta = beta,
#               sigma = sigma,
#               p = p, 
#               grupos = grupos,
#               alpha = alpha))
#           }