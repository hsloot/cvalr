#' @include allClass.R
NULL

setGeneric("getDimension",
  function(object) {
    standardGeneric("getDimension")
  })
setGeneric("setDimension<-",
  function(object, value) {
    standardGeneric("setDimension<-")
  })

setGeneric("getQMatrix",
  function(object) {
    standardGeneric("getQMatrix")
  })
setGeneric("setQMatrix<-",
  function(object, value) {
    standardGeneric("setQMatrix<-")
  })

setGeneric("getExIntensities",
  function(object) {
    standardGeneric("getExIntensities")
  })
setGeneric("setExIntensities<-",
  function(object, value) {
    standardGeneric("setExIntensities<-")
  })

setGeneric("getBernsteinFunction",
  function(object) {
    standardGeneric("getBernsteinFunction")
  })
setGeneric("setBernsteinFunction<-",
 function(object, value) {
   standardGeneric("setBernsteinFunction<-")
 })

setGeneric("getLambda",
  function(object) {
    standardGeneric("getLambda")
  })
setGeneric("setLambda<-",
  function(object, value) {
    standardGeneric("setLambda<-")
  })

setGeneric("getNu",
  function(object) {
    standardGeneric("getNu")
  })
setGeneric("setNu<-",
  function(object, value) {
    standardGeneric("setNu<-")
  })

setGeneric("getRho",
  function(object) {
    standardGeneric("getRho")
  })
setGeneric("setRho<-",
  function(object, value) {
    standardGeneric("setRho<-")
  })

setGeneric("getTau",
  function(object) {
    standardGeneric("getTau")
  })
setGeneric("setTau<-",
  function(object, value) {
    standardGeneric("setTau<-")
  })

setGeneric("getAlpha",
  function(object) {
    standardGeneric("getAlpha")
  })
setGeneric("setAlpha<-",
  function(object, value) {
    standardGeneric("setAlpha<-")
  })

setGeneric("invRho",
  function(object, value) {
    standardGeneric("invRho")
  })
setGeneric("invTau",
  function(object, value) {
    standardGeneric("invTau")
  })
setGeneric("invAlpha",
  function(object, value) {
    standardGeneric("invAlpha")
  })

setGeneric("constructBernsteinFunction",
  function(object, ...) {
    standardGeneric("constructBernsteinFunction")
  })
