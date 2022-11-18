setClassUnion("Y",   c("numeric",   "matrix", "array"))

setClass("Normal",
         representation(Y = "Y"))

setClass("MixNormal",
         representation(Y = "Y"),
         contains="Normal")

setClass("MoENormal",
         representation(Y = "Y"),
         contains = c("MixNormal"))

setClass("MoEKernelNormal",
         representation(Y = "Y"),
         contains = c("MixNormal"))

setClass("SN",
         representation(Y = "Y"))

setClass("MixSN",
         representation(Y = "Y"),
         contains = c("MixNormal", "Normal"))

setClass("MoESN",
         representation(Y = "Y"),
         contains = c("MixSN", "SN"))

setClass("resultadosEM",
         representation(resultados = "list"))

