#' Invisible Fence model selection (Linear Model)
#'
#' Invisible Fence model selection (Linear Model)
#'
#' @param full formula of full model
#' @param data data
#' @param B number of bootstrap sample, parametric for lm
#' @param cpus number of parallel computers
#' @param lftype subtractive measure type, e.g., absolute value of coefficients, p-value, t-value, etc.
#' @details This method (Jiang et. al, 2011) is motivated by computational expensive in complex and high dimensional problem.
#' The idea of the method--there is the best model in each dimension (in model space).  The boostrapping determines the coverage
#' probability of the selected model in each dimensions. The parsimonious model is the selected model with the highest coverage probabily
#' (except the one for the full model, always probability of 1.)
#' @return
#' \item{full}{list the full model}
#' \item{B}{list the number of bootstrap samples that have been used}
#' \item{freq}{list the coverage probabilities of the selected model for each dimension}
#' \item{size}{list the number of variables in the parsimonious model}
#' \item{term}{list variables included in the full model}
#' \item{model}{list the variables selected in-the-order in the parsimonious model}
#'  @note The current Invisible Fence focuses on variable selection. The current routine is applicable to the case in which
#'  the subtractive measure is the absolute value of the coefficients, p-value, t-value.
#'  However, the method can be extended to other subtractive measures.  See Jiang et. al (2011) for more details.
#' @author Jiming Jiang  Jianyang Zhao  J. Sunil Rao  Thuan Nguyen
#' @references
#' \itemize{
#'  \item{Jiang J., Rao J.S., Gu Z., Nguyen T. (2008),  Fence Methods for Mixed Model Selection. The Annals of Statistics, 36(4): 1669-1692}
#'  \item{Jiming Jiang, Thuan Nguyen and J. Sunil Rao (2011), Invisible fence methods and the identification of differentially expressed gene sets. Statistics and Its Interface, Volume 4, 403-415.}
#' }
#'
#' @examples
#' library(fence)
#' library(MASS)
#' library(snow)
#' r =1234; set.seed(r)
#' p=10; n=300; rho = 0.6
#' R = diag(p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'      R[i,j] = rho^(abs(i-j))
#'   }
#' }
#' R = 1*R
#' x=mvrnorm(n, rep(0, p), R)
#' colnames(x)=paste('x',1:p, sep='')
#' X = cbind(rep(1,n),x)
#' tbetas = c(1,1,1,0,1,1,0,1,0,0,0)  # non-zero beta 1,2,4,5,7
#' epsilon = rnorm(n)
#' y = as.matrix(X)%*%tbetas + epsilon
#' colnames(y) = 'y'
#' data = data.frame(cbind(X,y))
#' full = y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10
#' # Takes greater than 5 seconds (~`17 seconds) to run
#' # obj1 = IF.lm(full = full, data = data, B = 100, lftype = "abscoef")
#' # sort((names(obj1$model$coef)[-1]))  
#' # obj2 = IF.lm(full = full, data = data, B = 100, lftype = "pvalue")
#' # sort(setdiff(names(data[c(-1,-12)]), names(obj2$model$coef)))
#' 
#' @export


IF.lm = function(
  full, data, B = 100, cpus = 2, lftype = c("abscoef", "pvalue")) {

  full = cleanformula(full) # remove random effect
  #my first comment
  # model fit function
  mf = lm
  # lack of fit function
  lftype = match.arg(lftype)
  lf = switch(lftype,
              abscoef = function(res) abs(coef(res))[-1],
              pvalue = function(res) summary(res)$coefficients[-1,4])

  # bootstrap sample
  bs = bootstrap.lm(B, full, data)

  sfInit(parallel = TRUE, cpus = cpus)
  sfExportAll()
  invisiblefence(mf, full, data, lf, bs)
}

RIF.lm = function(
  full, data, B = 100, cpus = 2) {
  ans = IF.lm(full, data, B, cpus)

  B = ans$B
  times = ans$freq * B
  size = 1:length(ans$freq)
  pvalues = sapply(size, function(i) pbirthday(B, choose(length(ans$freq), i), times[i]))

  ans$freq = NULL
  ans$pvalue = pvalues
  ans$size = which.min(pvalues)
  ans$model = ans$terms[1:ans$size]
  class(ans) = "RIF"
  ans
}

ad_hoc_relative = function(B, times, nogs, size) {
  hugepbirthday <- chooseZ <- NULL
  rm(hugepbirthday, chooseZ)
  ans = hugepbirthday(B, chooseZ(nogs, size), times)
  names(ans) = NULL
  ans
}


bootstrap.lm = function(B, full, data) {
  X = model.matrix(full, data)
  model = lm(full, data)
  sigma = summary(model)$sigma
  beta = coef(model)
  bootsmp = matrix(as.vector(X %*% beta) + rnorm(nrow(X) * B, 0, sigma), nrow = nrow(X))

  ans = replicate(B, data, FALSE)
  for (i in 1:B) {
    ans[[i]][,deparse(full[[2]])] = bootsmp[,i]
  }
  ans
}

# TODO
cleanformula = function(full) {
  full=as.formula(full)
  resp = as.character(full)[2]
  tms = attributes(terms(full))$term.labels
  fr = grepl("\\|", tms)
  fixs = tms[!fr]
  cleanfull=paste(resp, '~', paste(fixs, sep='', collapse='+'), sep='')
  return(as.formula(cleanfull))
}
