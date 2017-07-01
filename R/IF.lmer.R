#' Invisible Fence model selection (Linear Mixed Model)
#'
#' Invisible Fence model selection (Linear Mixed Model)
#'
#' @param full formula of full model
#' @param data data
#' @param B number of bootstrap sample, parametric for lmer
#' @param REML Restricted maximum likelihood estimation
#' @param method choose either marginal (e.g., GEE) or conditional model
#' @param lftype subtractive measure type, e.g., absolute value of coefficients, p-value, t-value, etc.
#' @param cpus Number of parallel computers
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
#' require(fence)
#' library(snow)
#' library(MASS)
#' data("X.lmer")
#' data = data.frame(X.lmer)
#' # non-zero beta I.col.2, I.col.3a, I.col.3b, V5, V7, V8, V9
#' beta = matrix(c(0, 1, 1, 1, 1, 0, 0.1, 0.05, 0.25, 0), ncol = 1)
#' set.seed(1234)
#' alpha = rep(rnorm(100), each = 3)
#' mu = alpha + as.matrix(data[,-1]) %*% beta
#' data$id = as.factor(data$id)
#' data$y = mu + rnorm(300)
#' raw = "y ~ (1|id)+I.col.2+I.col.3a+I.col.3b"
#' for (i in 5:10) {
#'     raw = paste0(raw, "+V", i)
#' }
#' full = as.formula(raw)
#' # The following output takes more than 5 seconds (~70 seconds) to run. 
#' 
#' # obj1.lmer = IF.lmer(full = full, data = data, B = 100, method="conditional",lftype = "abscoef")
#' # sort(obj1.lmer$model) 
#' 
#' # obj2.lmer = IF.lmer(full = full, data = data, B = 100, method="conditional",lftype = "tvalue")
#' # sort(obj2.lmer$model)
#' 
#' # Similarly, the following scenarios can be run
#' 
#' # obj2.lmer = IF.lmer(full = full, data = data, B = 100, method="conditional",lftype = "tvalue")
#' # sort(obj2.lmer$model)
#' # obj1.lm = IF.lmer(full = full, data = data, B = 100, method="marginal", lftype = "abscoef")
#' # sort(names(obj1.lm$model$coefficients[-1]))
#' # obj2.lm = IF.lmer(full = full, data = data, B = 100, method="marginal", lftype = "tvalue")
#' # sort(names(obj2.lm$model$coefficients[-1]))
#' 
#' @export
#' @import lme4

IF.lmer = function(
  full, data, B = 100, REML = TRUE, method = c("marginal", "conditional"),
  cpus = parallel::detectCores(), lftype = c("abscoef", "tvalue")) {
  lme4 <- NULL
  rm(lme4)

  method = match.arg(method)
  if (method == "marginal") {
    if(lftype =="abscoef") {
      return(IF.lm(full, data, B, lftype ="abscoef", cpus))
    }else{
      return(IF.lm(full, data, B, lftype ="pvalue", cpus))
    }
  }

  # model fit function
  mf = function(formula, data) lmer(formula, data, REML = REML)
  # lack of fit function
  lftype = match.arg(lftype)
  lf = switch(lftype,
              abscoef = function(res) abs(fixef(res))[-1],
              tvalue = function(res) abs(summary(res)$coefficients[-1,3]))

  # bootstrap sample
  bs = bootstrap.lmer(B, full, data, REML)

  sfInit(parallel = TRUE, cpus = cpus)
  sfExportAll()
  #sfLibrary(glmmADMB)
  sfLibrary(lme4)
  invisiblefence_old(mf, full, data, lf, bs)
}
