#' Nonadaptive Fence model selection
#'
#' Nonadaptive Fence model selection
#'
#' @param mf function for fitting the model
#' @param f formula of full model
#' @param ms list of formula of candidates models
#' @param d data
#' @param lf measure lack of fit (to minimize)
#' @param pf model selection criteria, e.g., model dimension
#' @param cn given a specific c value
#' @return 
#' \item{models}{list all model candidates in the model space}
#' \item{lack_of_fit}{list a vector of Qs for all model candidates}
#' \item{formula}{list the model of the selected parsimonious model}
#' \item{sel_model}{list the selected (parsimonious) model given the adaptive c value}
#' @author Jiming Jiang  Jianyang Zhao  J. Sunil Rao  Thuan Nguyen
#' @references 
#' \itemize{
#'  \item{Jiang J., Rao J.S., Gu Z., Nguyen T. (2008),  Fence Methods for Mixed Model Selection. The Annals of Statistics, 36(4): 1669-1692}
#'  \item{Jiang J., Nguyen T., Rao J.S. (2009), A Simplified Adaptive Fence Procedure. Statistics and Probability Letters, 79, 625-629}
#'  \item{Thuan Nguyen, Jie Peng, Jiming Jiang (2014), Fence Methods for Backcross Experiments.  Statistical Computation and Simulation, 84(3), 644-662}
#' }
#' @examples
#' \dontrun{
#' require(fence)
#'
#' #### Example 1 #####
#' data(iris)
#' full = Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + (1|Species)
#' test_naf = fence.lmer(full, iris, fence = "nonadaptive", cn = 12)
#' test_naf$sel_model
#' }
#' @export

nonadaptivefence = function(
  mf, f, ms, d, lf, pf, cn) {
  ans = list()

  mf = cmpfun(mf)

  if (missing(ms)) {
    stop("No candidate models specified!")
  }
  ans$models = ms

  eval_models = lapply(ms, function(m) {
    mf(m, d)
  })

  lack_of_fit = sapply(eval_models, lf)
  ans$lack_of_fit = lack_of_fit

  Q = lack_of_fit - min(lack_of_fit)
  pick = sapply(eval_models, pf)

  infence = Q < cn
  pick[!infence] = Inf
  inpick = pick == min(pick)
  Q[!inpick] = Inf
  index = which.min(Q)

  ans$formula = ms[[index]]
  ans$sel_model = mf(ans$formula, d)

  class(ans) = "NAF"
  return(ans)
}

summary.NAF = function(res) {
  print(res$sel_model)
}
