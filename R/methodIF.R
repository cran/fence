#' Invisible Fence model selection
#'
#' Invisible Fence model selection
#' 
#' @param mf Call function, for example: default calls: function(m, b) eblupFH(formula = m, vardir = D, data = b, method = "FH")
#' @param f Full model
#' @param d Dimension number
#' @param lf Measures lack of fit using function(res) -res$fit$goodness[1]
#' @param bs Bootstrap
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
#' \dontrun{
#' data("X.lmer")
#' data = data.frame(X.lmer)
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
#' obj1.lmer = IF.lmer(full = full, data = data, B = 100, method="conditional",lftype = "abscoef")
#' obj1.lmer$model$coefficients
#' obj2.lmer = IF.lmer(full = full, data = data, B = 100, method="conditional",lftype = "tvalue")
#' obj2.lmer$model$coefficients
#' obj1.lm = IF.lmer(full = full, data = data, B = 100, method="marginal", lftype = "abscoef")
#' obj1.lm$model$coefficients
#' obj2.lm = IF.lmer(full = full, data = data, B = 100, method="marginal", lftype = "tvalue")
#' obj2.lm$model$coefficients
#' }
#' @export


invisiblefence = function(
  # model and lack of fit related
  mf, f, d, lf,
  # bootstrap sample
  bs) {

  ans = list(full = f)
  mf = cmpfun(mf)

  if (missing(bs)) {
    stop("No bootstrap sample specified!")
  }

  boot_evaluations = sfClusterApplyLB(bs, function(b) {
    try(mf(f, b), silent = TRUE)
  })
  sfStop()

  bm = sapply(boot_evaluations, function(x) class(x)[1])
  bb = sum(bm != "try-error")
  if (bb != length(bs)) {
    warning(paste0("Some bootstrap sample are not avaiable, new bootstrap size is ", sum(bb)))
  }
  B = bb
  ans$B = bb

  boot_evaluations = boot_evaluations[bm != "try-error"]

  s_mat = sapply(boot_evaluations, lf)
  orders_mat = apply(s_mat, 2, function(x) order(x, decreasing = TRUE))

  freq = sapply(1:nrow(orders_mat), function(size) {
    max(table(apply(orders_mat, 2, function(x) do.call(paste, as.list(sort(x[1:size]))))))
  })
  freq = freq / B

  size = peakglobal(freq[1:(nrow(orders_mat) - 1)]) 
  orlof = lf(mf(f, d))
  terms = names(sort(orlof, decreasing = TRUE))
  model = terms[1:size]
  ans$modeltest = names(sort(table(apply(orders_mat, 2, function(x) do.call(paste, as.list(sort(x[1:size]))))))[size])

  ans$freq = freq
  ans$size = size
  ans$terms = terms
  ans$model = mf(as.formula(paste0(as.character((f)[2]), "~", gsub(", ", "+", toString(model)))), d)
  class(ans) = "IF"
  return(ans)
}

invisiblefence_old = function(
  # model and lack of fit related
  mf, f, d, lf,
  # bootstrap sample
  bs) {

  ans = list(full = f)
  mf = cmpfun(mf)

  if (missing(bs)) {
    stop("No bootstrap sample specified!")
  }

  boot_evaluations = sfClusterApplyLB(bs, function(b) {
    try(mf(f, b), silent = TRUE)
  })
  sfStop()

  bm = sapply(boot_evaluations, function(x) class(x)[1])
  bb = sum(bm != "try-error")
  if (bb != length(bs)) {
    warning(paste0("Some bootstrap sample are not avaiable, new bootstrap size is ", sum(bb)))
  }
  B = bb
  ans$B = bb

  boot_evaluations = boot_evaluations[bm != "try-error"]

  s_mat = sapply(boot_evaluations, lf)
  orders_mat = apply(s_mat, 2, function(x) order(x, decreasing = TRUE))

  freq = sapply(1:nrow(orders_mat), function(size) {
    max(table(apply(orders_mat, 2, function(x) do.call(paste, as.list(sort(x[1:size]))))))
  })
  freq = freq / B

  size = peakglobal(freq[1:(nrow(orders_mat) - 1)])
  orlof = lf(mf(f, d))
  terms = names(sort(orlof, decreasing = TRUE))
  model = terms[1:size]
  ans$modeltest = names(sort(table(apply(orders_mat, 2, function(x) do.call(paste, as.list(sort(x[1:size]))))))[size])

  ans$freq = freq
  ans$size = size
  ans$terms = terms
  ans$model = model
  class(ans) = "IF"
  return(ans)
}
