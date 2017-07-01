#' Adaptive Fence model selection (Small Area Estmation)
#' 
#' Adaptive Fence model selection (Small Area Estmation)
#'
#' @param grid grid for c
#' @param bandwidth bandwidth for kernel smooth function
#' @param mf Call function, for example: default calls: function(m, b) eblupFH(formula = m, vardir = D, data = b, method = "FH")
#' @param f Full Model
#' @param ms find candidate model, findsubmodel.fh(full)
#' @param d Dimension number
#' @param lf Measures lack of fit using function(res) -res$fit$goodness[1]
#' @param pf Dimensions of model
#' @param bs Bootstrap
#' @param method Method to be used. Fay-Herriot method is the default.
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
#' \dontrun{
#' require(fence)
#' ### example 1 ####
#' data("kidney")
#' data = kidney[-which.max(kidney$x),]     # Delete a suspicious data point #
#' data$x2 = data$x^2
#' data$x3 = data$x^3
#' data$x4 = data$x^4
#' data$D = data$sqrt.D.^2
#' plot(data$y ~ data$x)
#' full = y~x+x2+x3+x4
#' testfh = fence.sae(full, data, B=1000, fence="adaptive", method="F-H", D = D)
#' testfh$sel_model
#' testfh$c
#' }
#' @export

adaptivefence.fh = function(
  # model and lack of fit related
  mf, f, ms, d, lf, pf,
  # bootstrap sample
  bs,
  # fence related
  grid = 101, bandwidth, method) {

  ans = list(full = f, models = ms, pickfunc = pf)
  mf = cmpfun(mf)

  if (missing(ms)) {
    stop("No candidate models specified!")
  }

  if (missing(bs)) {
    stop("No bootstrap sample specified!")
  }

  eval_models = sfClusterApplyLB(ms, function(m) {
    lapply(bs, function(b) {
      try(mf(m, b), silent = TRUE)
    })
  })
  sfStop()

  em = sapply(eval_models, function(eval_model) sapply(eval_model, class))
  eb = rowSums(em == "try-error") == 0
  if (sum(eb) != length(bs)) {
    warning(paste0("Some bootstrap sample are not avaiable, new bootstrap size is ", sum(eb)))
  }
  B = sum(eb)
  ans$B = sum(eb)
  for (i in 1:length(ms)) {
    eval_models[[i]] = eval_models[[i]][eb]
  }
  
if(method=='F-H'){
  em2 = sapply(eval_models, function(eval_model) sapply(eval_model, function(mod){mod$fit$convergence==TRUE}))
  eb2= rowSums(em2) == length(ms)
  if (sum(eb2) != length(bs)) {
    warning(paste0("Some bootstrap sample resulted in unconverged eblupFH, new bootstrap size is ", sum(eb2)))
  }
  B = sum(eb2)
  ans$B = sum(eb2)
  for (i in 1:length(ms)) {
    eval_models[[i]] = eval_models[[i]][eb2]
  }
}

  mi = 0
  bi = 0
  lack_of_fit_matrix = replicate(length(ms), {
    mi <<- mi + 1
    bi <<- 0
    replicate(B, {
      bi <<- bi + 1
      lf(eval_models[[mi]][[bi]])
    })
  })
  ans$lack_of_fit_matrix = lack_of_fit_matrix

  mi = 0
  bi = 0
  pick_matrix = replicate(length(ms), {
    mi <<- mi + 1
    bi <<- 0
    replicate(B, {
      bi <<- bi + 1
      pf(eval_models[[mi]][[bi]])
    })
  })
  ans$pick_matrix = pick_matrix

  rm(mi, bi)

  Q_m = sweep(lack_of_fit_matrix, 1, apply(lack_of_fit_matrix, 1, min), '-')
  ans$Qd_matrix = Q_m
  lof_lower = 0
  lof_upper = max(Q_m)
  cs = seq(lof_lower, lof_upper, length.out = grid)

  if (is.na(bandwidth)) {
    bandwidth = (cs[2] - cs[1]) * 3
  }
  ans$bandwidth = bandwidth * 1

  model_mat = matrix(NA, nrow = B, ncol = grid)
  for (i in 1:length(cs)) {
    infence_matrix = Q_m <= cs[i]
    for (bi in 1:B) {
      b_infence = infence_matrix[bi,]
      b_lack = lack_of_fit_matrix[bi,]
      b_pick = pick_matrix[bi,]
      b_pick[!b_infence] = Inf
      b_pick = which(b_pick == min(b_pick))
      model_mat[bi, i] = b_pick[which.min(b_lack[b_pick])]
    }
  }
  ans$model_mat = model_mat

  # if two models have same frequency, this frequency must
  # be lower than 0.5, so maybe we don't have to worry about
  # this case too much?

  freq_mat = apply(model_mat, 2, function(l) {
    tab = sort(table(l), decreasing = TRUE)
    c(as.numeric(names(tab)[1]), tab[1])
  })
  freq_mat[2,] = freq_mat[2,] / B
  freq_mat = rbind(freq_mat, ksmooth(cs, freq_mat[2,], kernel = "normal", bandwidth = bandwidth, x.points = cs)$y)

  colnames(freq_mat) = cs
  rownames(freq_mat) = c("index", "frequency", "smooth_frequency")
  ans$freq_mat = freq_mat

  cindex = peakw(cs, freq_mat[3,], 2)
  ans$c = cs[cindex]

  if (!is.na(cindex)) {
    ans$formula = ms[[freq_mat[1,cindex]]]
    ans$sel_model = mf(ans$formula, d)
  } else {
    ans$formula = NA
    ans$sel_model = NA
  }
  class(ans) = "AF"
  return(ans)
}

