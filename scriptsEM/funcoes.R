dMix = function(Y, X, params, ...){
  UseMethod("dMix")
}

vero = function(Y, X, params, ...){
  UseMethod("vero")
}

estimaTeta = function(Y, X, ...){
  UseMethod("estimaTeta")
}

estimaMedia = function(X, beta, ...){
  UseMethod("estimaMedia")
}

estimaSigma = function(Y, X, beta, ...){
  UseMethod("estimaSigma")
}

chuteInicial = function(Y, X, ...){
  UseMethod("chuteInicial")
}

etapaE = function(Y, params, medias, ...){
  UseMethod("etapaE")
}

etapaM = function(Y, X, U, params, ...){
  UseMethod("etapaM")
}

calculaMetricas = function(y, medias){
  UseMethod("calculaMetricas")
}

predictReg = function(reg, ...){
  UseMethod("predictReg")
}
