#' Adaptive Fence model selection
#'
#' Adaptive Fence model selection
#'
#' @param mf function for fitting the model
#' @param f formula of full model
#' @param ms list of formula of candidates models
#' @param d data
#' @param lf measure lack of fit (to minimize)
#' @param pf model selection criteria, e.g., model dimension
#' @param bs bootstrap samples
#' @param grid grid for c
#' @param bandwidth bandwidth for kernel smooth function
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
#' test_af = fence.lmer(full, iris)
#' plot(test_af)
#' test_af$sel_model
#' 
#' #### Example 2 #####
#' r =1234; set.seed(r)  
#' p=8; n=150; rho = 0.6
#' id = rep(1:50,each=3)
#' R = diag(p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'      R[i,j] = rho^(abs(i-j))
#'   }
#' } 
#' R = 1*R
#' x=mvrnorm(n, rep(0, p), R)  # all x's are time-varying dependence #
#' colnames(x)=paste('x',1:p, sep='')
#' tbetas = c(0,0.5,1,0,0.5,1,0,0.5)  # non-zero beta 2,3,5,6,8
#' epsilon = rnorm(150)
#' y = x%*%tbetas + epsilon
#' colnames(y) = 'y'
#' data = data.frame(cbind(x,y,id))
#' full = y ~ x1+x2+x3+x4+x5+x6+x7+x8+(1|id)
#' #X = paste('x',1:p, sep='', collapse='+')
#' #full = as.formula(paste('y~',X,'+(1|id)', sep=""))  #same as previous one
#' fence_obj = fence.lmer(full,data)   # it takes 3-5 min #
#' plot(fence_obj)
#' fence_obj$sel_model
#' }
#' @export

adaptivefence = function(
  # model and lack of fit related
  mf, f, ms, d, lf, pf,
  # bootstrap sample
  bs,
  # fence related
  grid = 101, bandwidth) {

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
  }
  else {
    ans$formula = NA
    ans$sel_model = NA
  }
  class(ans) = "AF"
  return(ans)
}

#' Plot Adaptive Fence model selection
#' 
#' @param x Object to be plotted
#' @param ... Additional arguments. CNot currently used.
#'
#' @export
plot.AF = function(x = res, ...) {
  res <- m <- sp <- NULL
  rm(res,m,sp)
  tmp = data.frame(c = as.numeric(colnames(x$freq_mat)),
                   p = x$freq_mat[2, ],
                   sp = x$freq_mat[3, ],
                   m = as.factor(x$freq_mat[1, ]))
  p = ggplot(tmp) +
    geom_point(aes(x = c, y = p, colour = m)) +
    geom_line(aes(x = c, y = sp), linetype="dashed") +
    geom_line(aes(x = c, y = sp, colour = m)) +
    ylim(0, 1)
  p
}

#' Summary Adaptive Fence model selection
#' 
#' @param object Object to be summarized
#' @param ... addition arguments. Not currently used
#'
#' @export
summary.AF = function(object = res, ...) {
  res <- m <- sp <- NULL
  rm(res,m,sp)
  print(object$sel_model)
}
