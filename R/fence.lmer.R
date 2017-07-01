#' Fence model selection (Linear Mixed Model)
#'
#' Fence model selection (Linear Mixed Model)
#'
#' @param full formula of full model
#' @param data data
#' @param B number of bootstrap samples, parametric bootstrap is used
#' @param grid grid for c
#' @param cpus Number of parallel computers 
#' @param fence a procedure of the fence method to be used.  
#' It's suggested to choose nonadaptive procedure if c is known; otherwise nonadaptive must be chosen
#' @param cn cn value for nonadaptive
#' @param REML Restricted Maximum Likelihood approach
#' @param bandwidth bandwidth for kernel smooth function
#' @details In Jiang et. al (2008), the adaptive c value is chosen from the highest peak in the p* vs. c plot.  
#' In Jiang et. al (2009), 95\% CI is taken into account while choosing such an adaptive choice of c.
#' In Thuan Nguyen et. al (2014), the adaptive c value is chosen from the first peak. This approach works better in the 
#' moderate sample size or weak signal situations.  Empirically, the first peak becomes highest peak when sample size 
#' increases or signals become stronger 
#' @return 
#' \item{models}{list all model candidates in the model space}
#' \item{B}{list the number of bootstrap samples that have been used}
#' \item{lack_of_fit_matrix}{list a matrix of Qs for all model candidates (in columns). Each row is for each bootstrap sample}
#' \item{Qd_matrix}{list a matrix of QM - QM.tilde for all model candidates. Each row is for each bootrap sample}
#' \item{bandwidth}{list the value of bandwidth}
#' \item{model_mat}{list a matrix of selected models at each c values in grid (in columns). Each row is for each bootstrap sample}
#' \item{freq_mat}{list a matrix of coverage probabilities (frequency/smooth_frequency) of each selected models for a given c value (index)}
#' \item{c}{list the adaptive choice of c value from which the parsimonious model is selected}
#' \item{sel_model}{list the selected (parsimonious) model given the adaptive c value}
#'  @note The current Fence package focuses on variable selection. 
#'  However, Fence methods can be used to select other parameters of interest, e.g., tunning parameter, variance-covariance structure, etc.
#' @author Jiming Jiang  Jianyang Zhao  J. Sunil Rao  Thuan Nguyen
#' @references 
#' \itemize{
#'  \item{Jiang J., Rao J.S., Gu Z., Nguyen T. (2008),  Fence Methods for Mixed Model Selection. The Annals of Statistics, 36(4): 1669-1692}
#'  \item{Jiang J., Nguyen T., Rao J.S. (2009), A Simplified Adaptive Fence Procedure. Statistics and Probability Letters, 79, 625-629}
#'  \item{Thuan Nguyen, Jie Peng, Jiming Jiang (2014), Fence Methods for Backcross Experiments.  Statistical Computation and Simulation, 84(3), 644-662}
#' }
#' 
#' @examples
#' require(fence)
#' library(snow)
#'
#' #### Example 1 #####
#' data(iris)
#' full = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + (1|Species)
#' # Takes greater than 5 seconds to run
#' # test_af = fence.lmer(full, iris)
#' # test_af$c
#' # test_naf = fence.lmer(full, iris, fence = "nonadaptive", cn = 12)
#' # plot(test_af)
#' # test_af$sel_model
#' # test_naf$sel_model
#' @export
#' @import stats MASS compiler fields ggplot2 lme4 sae snowfall snow
#' @importFrom grDevices dev.off

fence.lmer = function(
  full, data, B = 100, grid = 101, fence = c("adaptive", "nonadaptive"),
  cn = NA, REML = TRUE, bandwidth = NA, cpus = parallel::detectCores()) {
  lme4 <- NULL
  rm(lme4)

  fence = match.arg(fence)
  if (fence == "adaptive" & !is.na(cn) |
      fence == "nonadaptive" & is.na(cn)) {
    stop("Adaptive agreement doesn't match!")
  }

  # find all candidate submodels
  ms = findsubmodel.lmer(full)
  # model fit function
  mf = function(formula, data) lmer(formula, data, REML = REML)
  # lack of fit function
  lf = function(res) -logLik(res)
  # pick up function
  pf = size.lmer

  if (fence == "nonadaptive") {
    return(nonadaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      cn = cn))
  }
  if (fence == "adaptive")    {
    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    sfLibrary(lme4)
    return(   adaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      bs = bootstrap.lmer(B, full, data, REML), grid = grid, bandwidth = bandwidth))
  }
}

findsubmodel.lmer = function(full) {
  resp = as.character(full)[2]
  tms = attributes(terms(full))$term.labels
  fr = grepl("\\|", tms)
  fixs = tms[!fr]
  rans = tms[fr]
  rans = paste("(", rans, ")", sep = "")
  res = paste(resp, "~", rans, sep = "")
  for (fix in fixs) {
    res = as.vector(sapply(res, function(x) paste(x, c("", paste("+", fix, sep = "")), sep = "")))
  }
  lapply(res, as.formula)
}

bootstrap.lmer = function (B, f, data, REML) {
  full = lmer(formula = f, data = data, REML = REML)
  ans = replicate(B, data, FALSE)
  X = full@pp$X
  beta = fixef(full)
  Z = t(as.matrix(full@pp$Zt))
  # random intercept only
  # can't handle random slope or interactive term for now
  # TODO: sigma = full@pp$Lambdat * tau
  #       sigma = sigma %*% t(sigma)
  tau = attr(VarCorr(full), "sc")
  sigmas = diag(as.matrix(full@pp$Lambdat) * tau)
  n = nrow(X)

  # generate
  fe = X %*% beta
  a = matrix(rnorm(B * length(sigmas), 0, sigmas), nrow = length(sigmas))
  re = Z %*% a
  bootsmp = as.matrix(as.vector(fe) + re + rnorm(n * B, 0, tau))
  for (i in 1:B) {
    ans[[i]][,deparse(f[[2]])] = bootsmp[,i]
  }
  ans
}

size.lmer = function(res) {
  length(fixef(res))
}
