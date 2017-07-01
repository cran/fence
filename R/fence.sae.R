#' Fence model selection (Small Area Estmation)
#'
#' Fence model selection (Small Area Estmation)
#'
#' @param full formular of full model
#' @param data data
#' @param B number of bootstrap sample, parametric for lmer
#' @param grid grid for c
#' @param fence fence method to be used, e.g., adaptive, or nonadaptive.
#' It's suggested to choose nonadaptive procedure if c is known; otherwise nonadaptive must be chosen
#' @param cn cn for nonadaptive
#' @param REML Restricted Maximum Likelihood approach
#' @param bandwidth bandwidth for kernel smooth function
#' @param method Select method to use
#' @param D vector containing the D sampling variances of direct estimators for each domain. The values must be sorted as the variables in formula. Only used in FH model
#' @param cpus Number of parallel computers 
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
#' @note 
#' \itemize{
#' \item{The current Fence package focuses on variable selection. 
#'  However, Fence methods can be used to select other parameters of interest, e.g., tunning parameter, variance-covariance structure, etc.}
#' \item{The number of bootstrap samples is suggested to be increased, e.g., B=1000 when the sample size is small, or signals are weak}
#' }
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
#' ### example 1 ####
#' data("kidney")
#' data = kidney[-which.max(kidney$x),]     # Delete a suspicious data point #
#' data$x2 = data$x^2
#' data$x3 = data$x^3
#' data$x4 = data$x^4
#' data$D = data$sqrt.D.^2
#' plot(data$y ~ data$x)
#' full = y~x+x2+x3+x4
#' # Takes more than 5 seconds to run
#' # testfh = fence.sae(full, data, B=100, fence="adaptive", method="F-H", D = D)
#' # testfh$sel_model
#' # testfh$c
#' @export


fence.sae = function(
  full, data, B = 100, grid = 101, fence = c("adaptive", "nonadaptive"),
  cn = NA, method = c("F-H", "NER"), D = NA, REML = FALSE, bandwidth = NA, 
  cpus = parallel::detectCores()) {
  sae <- NULL
  rm(sae)

  fence = match.arg(fence)
  if (fence == "adaptive" & !is.na(cn) |
      fence == "nonadaptive" & is.na(cn)) {
    stop("Adaptive agreement doesn't match!")
  }

  method = match.arg(method)
  if (method == "NER") {
    return(fence.lmer(full, data, B, grid, fence, cn, REML, bandwidth, cpus))
  }
  # if (method == "F-H") {
  #   return(fence.fh(full, data, B, grid, fence, cn, D, bandwidth, cpus))
  # }
  # find all candidate submodels
  ms = findsubmodel.fh(full)
  # model fit function
  mf = function(m, b) eblupFH(formula = m, vardir = D, data = b, method = "FH")
  # lack of fit function
  lf = function(res) -res$fit$goodness[1]
  # pick up function
  pf = function(res) nrow(res$fit$estcoef) - 1

  if (fence == "nonadaptive") {
    return(nonadaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      cn = cn))
  }

  if (fence == "adaptive")    {
    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    sfLibrary(sae)
    return(   adaptivefence.fh(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      bs = bootstrap.fh(B, full, data, D), grid = grid, bandwidth = bandwidth, method=method))
  }
}

findsubmodel.fh = function(full) {
  resp = as.character(full)[2]
  tms = attributes(terms(full))$term.labels
  res = paste(resp, "~", sep = "")
  for (tm in tms) {
    res = as.vector(sapply(res, function(x) paste(x, c("", paste("+", tm, sep = "")), sep = "")))
  }
  res = gsub("~ +", "~", res)
  res = res[res != "y~"]
  lapply(res, as.formula)
}

bootstrap.fh = function(B, full, data, D) {
  X = model.matrix(full, data)
  model = eblupFH(formula = full, vardir = D, data = data, method = "FH")
  beta = model$fit$estcoef[,1]
  tau = model$fit$refvar

  ans = replicate(B, data, FALSE)
  bootsmp = as.vector(X %*% beta) + replicate(B, rnorm(nrow(X), 0, sqrt(tau + data[,deparse(substitute(D))])))
  for (i in 1:B) {
    ans[[i]][,deparse(full[[2]])] = bootsmp[,i]
  }
  ans
}