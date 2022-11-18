library(mixsmsn)

# Função

MoESkewNormal <- function(dados, X, R, nGrupos = 2, maxIter = 100, 
                          tol = 1e-6, ordem, grupo){
  # Pacotes exigidos e funções auxiliares
  
  library(moments)
  library(sn)
  library(mixtools)
  
  matrizP <- function(alpha, R){
    P1 <- exp(R%*%t(alpha))/(1 + rowSums(exp(R%*%t(alpha))))
    P <- cbind(P1, 1 - rowSums(P1))
    return(P)
  }
  
  dMixSN <- function(y, X, beta, s, a, p, D){
    total <- matrix(NA, nrow = n, ncol = g)
    for(j in 1:g){
      media <- X%*%beta[j,] + b*D[j]
      total[,j] <- p[,j]*dsn(y, media, s[j], a[j])
    }
    return(rowSums(total))
  }
  
  logVero <- function(beta, s, a, p, D){
    vetor <- log(dMixSN(y, X, beta, s, a, p, D))
    return(sum(vetor))
  }
  
  erroPadraoMoESN2 <- function(y, X, R, beta, s, a, p){
    g <- nrow(beta)
    q <- ncol(beta)
    k <- ncol(R)
    n <- length(y)
    
    S <- matrix(NA, nrow = n, ncol = g*(q + 2) + (g-1)*k)
    
    dMixSN <- function(i, beta, s, a, p, D){
      vetor <- numeric(g)
      for(j in 1:g){
        media <- t(X[i,])%*%beta[j,] + b*D[j]
        vetor[j] <- p[j]*dsn(y[i], media, s[j], a[j])
      }
      return(sum(vetor))
    }
    
    for(i in 1:n){
      scoreBeta <- matrix(NA, nrow = g, ncol = q)
      scoreAlpha <- matrix(NA, nrow = g-1, ncol = k)
      scoreS <- scoreA <- numeric(g)
      
      yi <- y[i]
      x <- X[i,]
      pi <- p[i,]
      
      for(j in 1:g){
        delta <- a/sqrt(1+a^2)
        D <- s*delta
        z <- as.numeric((yi - t(x)%*%beta[j,] - b*D[j])/s[j])
        
        f <- dMixSN(i, beta, s, a, pi, D)
        
        scoreBeta[j,] <- (pi[j]/f)*(2*x/s[j]^2)*dnorm(z)*(z*pnorm(a[j]*z) - a[j]*dnorm(a[j]*z))
        
        parte1S <- pnorm(a[j]*z)*(z*((yi - t(x)%*%beta[j,])/s[j]) - 1)
        parte2S <- dnorm(a[j]*z)*(a[j]*((yi - t(x)%*%beta[j,])/s[j]))
        scoreS[j] <- (pi[j]/f)*(2/s[j]^2)*dnorm(z)*(parte1S - parte2S)
        
        parte1A <- z*b*pnorm(a[j]*z)*(((1 + a[j]^2)^0.5 - (a[j]^2)*(1 + a[j]^2)^(-0.5))/(1 + a[j]^2))
        parte2A <- (dnorm(a[j]*z)/s[j])*(yi - t(x)%*%beta[j,] - b*s[j]*a[j]*((2*(1 + a[j]^2)^0.5 - (a[j]^2)*(1 + a[j]^2)^(-0.5))/(1 + a[j]^2)))
        scoreA[j] <- (pi[j]/f)*(2/s[j])*dnorm(z)*(parte1A + parte2A)
        
        if(j != g){
          parte1Alpha <- pi[j]*(1 - pi[j])*dsn(yi, (t(x)%*%beta[j,] + b*D[j]), s[j], a[j])*R[i,]
          parte2Alpha <- pi[g]*(1 - pi[g])*dsn(yi, (t(x)%*%beta[g,] + b*D[g]), s[g], a[g])*R[i,]
          scoreAlpha[j,] <- (1/f)*(parte1Alpha - parte2Alpha)
        }
      }
      
      S[i, 1:(g*q)] <- as.vector(t(scoreBeta))
      S[i, (g*q + 1):(g*q + g)] <- scoreS
      S[i, (g*q + g + 1):(g*q + 2*g)] <- scoreA
      S[i, (g*q + 2*g + 1):(g*(q + 2) + (g-1)*k)] <- as.vector(t(scoreAlpha))
    }
    
    erros <- round(sqrt(diag(solve(t(S)%*%S, tol = 1e-161))), 6)
    nomes <- character(g*(q + 3) - 1)
    nomes[1:(g*q)] <- paste0('beta', rep(1:g, rep(q, g)), rep(0:(q-1), g))
    nomes[(g*q + 1):(g*q + g)] <- paste0('sigma', 1:g)
    nomes[(g*q + g + 1):(g*q + 2*g)] <- paste0('ass', 1:g)
    nomes[(g*q + 2*g + 1):(g*(q + 2) + (g-1)*k)] <- paste0('alpha', rep(1:(g-1), rep(k, g-1)), rep(0:(k-1), g-1))
    names(erros) <- nomes
    return(erros)
  }
  
  # Chute inicial
  
  y <- dados
  dados <- as.data.frame(cbind(y, X[,-1]))
  colnames(dados)[1] <- 'Y'
  
  p <- ncol(X)
  k <- ncol(R)
  g <- nGrupos
  n <- length(y)
  
  km <- kmeans(dados, centers = g, iter.max = 100, nstart = 50)
  aux <- as.numeric(names(sort(table(km$cluster))))
  
  betaNova <- matrix(NA, nrow = g, ncol = p)
  alphaNova <- matrix(0, nrow = g-1, ncol = k)
  pNova <- matrix(NA, nrow = n, ncol = g)
  dpNova <- assNova <- numeric(g)
  
  for(i in 1:g){
    dados0 <- dados[km$cluster == aux[ordem[i]],]
    modelo <- lm(Y ~ ., data = dados0)
    betaNova[i,] <- as.numeric(modelo$coefficients)
    dpNova[i] <- sqrt(sum((dados0$Y - modelo$fitted.values)^2)/(nrow(dados0) - p))
    assNova[i] <- skewness(modelo$residuals)
    pNova[,i] <- rep(nrow(dados0)/n, n) 
  }
  
  # Algoritmo
  
  crit <- 1
  c <- 1
  
  while(crit > tol & c < maxIter){
    # Parâmetros atuais
    
    pAtual <- pNova
    alphaAtual <- alphaNova
    betaAtual <- betaNova
    dpAtual <- dpNova
    assAtual <- assNova
    
    deltaAtual <- assAtual/sqrt(1 + assAtual^2)
    
    DeltaAtual <- dpAtual*deltaAtual
    
    GamaAtual <- (dpAtual^2)*(1 - deltaAtual^2)
    
    # Etapa E
    
    b <- -sqrt(2/pi)
    Z <- S2 <- S3 <- Mu <- Tau <- matrix(0, nrow = n, ncol = g)
    M2T <- GamaAtual/(GamaAtual + DeltaAtual^2)
    MT <- sqrt(M2T)
    
    for(j in 1:g){
      media <- as.vector(X%*%betaAtual[j,] + b*DeltaAtual[j])
      Z[,j] <- pAtual[,j]*dsn(y, media, dpAtual[j], 
                              assAtual[j])/dMixSN(y, X, betaAtual, dpAtual, 
                                                  assAtual, pAtual, DeltaAtual)
      
      Mu[,j] <- (DeltaAtual[j]/(GamaAtual[j] + DeltaAtual[j]^2))*(y - media)
      
      aux <- (pnorm(Mu[,j]/MT[j]) == 0)
      vetorAux <- numeric(n)
      vetorAux[aux] <- (dnorm(Mu[,j]/MT[j])[aux])/(.Machine$double.xmin)
      vetorAux[!aux] <- (dnorm(Mu[,j]/MT[j])[!aux])/(pnorm(Mu[,j]/MT[j])[!aux])
      
      Tau[,j] <- vetorAux
      
      S2[,j] <- Z[,j]*(Mu[,j] + b + MT[j]*Tau[,j])
      
      S3[,j] <- Z[,j]*((Mu[,j] + b)^2 + M2T[j] + (Mu[,j] + 2*b)*MT[j]*Tau[,j])
      
    }
    
    # Etapa M
    
    DeltaNova <- GamaNova <- dpNova <- assNova <- numeric(g)
    betaNova <- matrix(0, nrow = g, ncol = p)
    alphaNova <- matrix(0, nrow = g-1, ncol = k)
    
    for(j in 1:g){
      ## Proporção
      
      if(j != g){
        alphaNova[j,] <- alphaAtual[j,] + 4*solve(t(R)%*%R)%*%(t(R)%*%(Z[,j] - pAtual[,j]))
      }
      
      ## Coeficientes
      
      parte1Beta <- solve(t(X)%*%diag(Z[,j])%*%X)
      parte2Beta <- t(X)%*%(Z[,j]*y - S2[,j]*DeltaAtual[j])
      
      betaNova[j,] <- parte1Beta%*%parte2Beta
      
      ## Delta
      
      parte1Delta <- sum(S2[,j]*(y - X%*%betaAtual[j,]))
      parte2Delta <- sum(S3[,j])
      
      DeltaNova[j] <- parte1Delta/parte2Delta
      
      ## Gama
      
      parte11Gama <- Z[,j]*(y - X%*%betaAtual[j,])^2
      parte12Gama <- 2*DeltaAtual[j]*S2[,j]*(y - X%*%betaAtual[j,])
      parte13Gama <- S3[,j]*(DeltaAtual[j]^2)
      
      parte1Gama <-  sum(parte11Gama - parte12Gama + parte13Gama)
      parte2Gama <- sum(Z[,j])
      
      GamaNova[j] <- parte1Gama/parte2Gama
      
      ## Desvio padrão e Assimetria
      
      dpNova[j] <- sqrt(GamaNova[j] + (DeltaNova[j]^2))
      assNova[j] <- DeltaNova[j]/sqrt(GamaNova[j])
    }
    
    pNova <- matrizP(alphaNova, R)
    
    # Critério de parada
    
    crit <- abs(logVero(betaNova, dpNova, assNova, pNova, DeltaNova)/logVero(betaAtual, dpAtual, assAtual, pAtual, DeltaAtual) - 1)
    c <- c + 1
  }
  
  # Resultados
  
  erros <- erroPadraoMoESN2(y, X, R, betaNova, dpNova, assNova, pNova)
  
  agrup <- as.numeric(apply(pNova, 1, which.max))
  tabela <- table(grupo, agrup, dnn = list('Real', 'Modelo'))
  
  resultados <- list(`Iterações` = c, `Coeficientes` = betaNova, 
                     `Desvio-padrão` = dpNova, Assimetria = assNova, 
                     `Coeficientes da proporção` = alphaNova, 
                     `Erros-padrão` = erros, 
                     `Classificação` = tabela)
  return(resultados)
}