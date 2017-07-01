#' Adaptive Fence model selection (Restricted Fence)
#'
#' Adaptive Fence model selection (Restricted Fence)
#' 
#' @param full formula of full model
#' @param data data
#' @param groups A list of formulas of (full) model in each bins (groups) of variables
#' @param B number of bootstrap sample, parametric for lmer
#' @param grid grid for c
#' @param bandwidth bandwidth for kernel smooth function
#' @param method either marginal (GEE) or conditional approach is selected
#' @param id Subject or cluster id variable
#' @param plot Plot object
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
#' @note bandwidth = (cs[2] - cs[1]) * 3. So it's chosen as 3 times grid between two c values.
#' @references 
#' \itemize{
#'  \item{Jiang J., Rao J.S., Gu Z., Nguyen T. (2008),  Fence Methods for Mixed Model Selection. The Annals of Statistics, 36(4): 1669-1692}
#'  \item{Jiang J., Nguyen T., Rao J.S. (2009), A Simplified Adaptive Fence Procedure. Statistics and Probability Letters, 79, 625-629}
#'  \item{Thuan Nguyen, Jiming Jiang (2012), Restricted fence method for covariate selection in longitudinal data analysis. Biostatistics, 13(2), 303-314}
#'  \item{Thuan Nguyen, Jie Peng, Jiming Jiang (2014), Fence Methods for Backcross Experiments.  Statistical Computation and Simulation, 84(3), 644-662}
#' }
#' 
#' @examples
#' \dontrun{
#' r =1234; set.seed(r)
#' n = 100; p=15; rho = 0.6
#' beta = c(1,1,1,0,1,1,0,1,0,0,1,0,0,0,0)  # non-zero beta 1,2,3,V6,V7,V9,V12
#' id = rep(1:n,each=3)
#' V.1 = rep(1,n*3)
#' I.1 = rep(c(1,-1),each=150)
#' I.2a = rep(c(0,1,-1),n)
#' I.2b = rep(c(0,-1,1),n)
#' x = matrix(rnorm(n*3*11), nrow=n*3, ncol=11)
#' x = cbind(id,V.1,I.1,I.2a,I.2b,x)
#' R = diag(3)
#' for(i in 1:3){
#'  for(j in 1:3){
#'    R[i,j] = rho^(abs(i-j))
#'  }
#' } 
#' e=as.vector(t(mvrnorm(n, rep(0, 3), R)))  
#' y = as.vector(x[,-1]%*%beta) + e
#' data = data.frame(x,y)
#' raw = "y ~ V.1 + I.1 + I.2a +I.2b"
#' for (i in 6:16) { raw = paste0(raw, "+V", i)}; full = as.formula(raw)
#' bin1="y ~ V.1 + I.1 + I.2a +I.2b"
#' for (i in 6:8) { bin1 = paste0(bin1, "+V", i)}; bin1 = as.formula(bin1)
#' bin2="y ~ V9"
#' for (i in 10:16){ bin2 = paste0(bin2, "+V", i)}; bin2 = as.formula(bin2)
#' # May take longer than 30 min since there are two stages in this RF procedure
#' obj1.RF = RF(full = full, data = data, groups = list(bin1,bin2), method="conditional")
#' obj1.RF$sel_model
#' obj2.RF = RF(full = full, data = data, groups = list(bin1,bin2), B=100, method="marginal")
#' obj2.RF$sel_model
#' }
#' @export

RF = function(full, data, groups, B = 100, grid = 101, bandwidth = NA, plot = FALSE, method = c("marginal", "conditional"), id = "id", cpus = parallel::detectCores()) {
  #if (is.na(bandwidth)) {
  # stop("Must assign a bandwidth.")
  #}

  method = match.arg(method)

  # 1st step
  cands = character(0)
  
  resp = as.character(full)[2]
  y = as.matrix(data[,resp])

  for (group in groups) {
    bs = bootstrap.RF(B, full, group, data)
    ms = findsubmodel.RF(full, group, data) 
    mf = function(model, bootsmp) {
      if (is.data.frame(bootsmp)) {return(model$model)}
      list(Q = t(bootsmp) %*% model$meat %*% bootsmp, 
           size = length(model$model))
    }
    lf = function(x) {x$Q}
    pf = function(x) {x$size}

    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    groupaf = adaptivefence(mf = mf, f = group, ms = ms, d = data, lf = lf, pf = pf, bs = bs, grid = grid, bandwidth = bandwidth)
    if (plot) {
      print(plot(groupaf))
      dev.off()
    }    
    groupc = groupaf$c

    Qs = sapply(ms, function(m) t(y) %*% m$meat %*% y)
    Qs = Qs - min(Qs)
    infence = Qs < groupc
    sizes = sapply(ms, function(x) length(x$model))
    sizes[!infence] = Inf
    Qs[sizes != min(sizes)] = Inf
    #cands = c(cands, ms[[which.min(Qs)]]$model)
    cands = c(cands, ms[[which.min(Qs)]]$model[-length(ms[[which.min(Qs)]]$model)])
  }
  cands = unique(cands)

  # 2nd step
  if (method == "marginal") {
    resp = as.character(full)[2]
    full2nd = as.formula(paste0(resp, "~", gsub(" ", "+", do.call(paste, as.list(cands)))))
    bs = bootstrap.RF2(B, full2nd, data)
    ms = findsubmodel.RF2(resp, cands)
    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    res = adaptivefence(mf = lm, f = full2nd, ms = ms, d = data, lf = function(x) -logLik(x), 
      pf = function(x) length(attributes(terms(x))$term.labels), bs = bs, bandwidth = NA)
  }

  if (method == "conditional") {
    resp = as.character(full)[2]
    full2nd = as.formula(paste0(resp, "~", gsub(" ", "+", do.call(paste, as.list(cands))), "+(1|", id, ")"))
    res = fence.lmer(full2nd, data, B = B, grid = grid)
  }

  class(res) = "RF"
  res
}

bootstrap.RF = function(B, full, group, data) {
  resp = as.character(full)[2]
  fulltms = attributes(terms(full))$term.labels
  grouptms = attributes(terms(group))$term.labels
  groupxtms = fulltms[!(fulltms %in%  grouptms)]

  y = data[,resp]

  grouptms = as.list(grouptms)
  grouptms$sep = "+"
  groupxtms = as.list(groupxtms)
  groupxtms$sep = "+"

  X1 = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, grouptms))), data)
  X2 = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, groupxtms))), data)

  PX2O = diag(nrow(X2)) - X2 %*% ginv(t(X2) %*% X2) %*% t(X2)
  PX2OMX1 = PX2O - PX2O %*% X1 %*% ginv(t(X1) %*% PX2O %*% X1) %*% t(X1) %*% PX2O

  df = nrow(data) - length(fulltms)
  beta1 = ginv(t(X1) %*% PX2O %*% X1) %*% t(X1) %*% PX2O %*% matrix(y, ncol = 1)

  sigma = matrix(y, nrow = 1) %*% PX2OMX1 %*% matrix(y, ncol = 1) / df
  base = X1 %*% beta1
  lapply(1:B, function(i) matrix(base + rnorm(nrow(data), 0, sqrt(sigma)), ncol = 1))
}

findsubmodel.RF = function(full, group, data) {
  resp = as.character(full)[2]
  fulltms = attributes(terms(full))$term.labels
  grouptms = attributes(terms(group))$term.labels
  groupxtms = fulltms[!(fulltms %in%  grouptms)]

  y = data[,resp]

  grouptms = as.list(grouptms)
  grouptms$sep = "+"
  groupxtms = as.list(groupxtms)
  groupxtms$sep = "+"

  X1 = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, grouptms))), data)
  X2 = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, groupxtms))), data)
  
  PX2O = diag(nrow(X2)) - X2 %*% ginv(t(X2) %*% X2) %*% t(X2)

  res = ""
  for (i in 1:length(grouptms)) {
    res = append(lapply(res, function(x) c(x, grouptms[i])), lapply(res, function(x) x))
  }
  res = lapply(res, function(x) x[-1])[-length(res)]
  for (i in 1:length(res)) {
    X1 = as.matrix(data[,unlist(res[[i]][-length(res[[i]])])])
    if (ncol(X1) == 0) {
      next
    }
    res[[i]] = list(model = res[[i]], meat = PX2O - PX2O %*% X1 %*% ginv(t(X1) %*% PX2O %*% X1) %*% t(X1) %*% PX2O)
  }
  #res = res[!(sapply(res,is.null))]
  res = res[!(sapply(res, function(l){is.null(l$meat)}))]
  res
}

bootstrap.RF2 = function(B, full2nd, data) {
  resp = as.character(full2nd)[2]
  m = lm(full2nd, data)
  base = fitted(m)
  sigma = summary(m)$sigma
  newresp = lapply(1:B, function(i) base + rnorm(length(base), 0, sigma))
  ans = list()
  for (i in 1:B) {
    data[,resp] = newresp[[i]]
    ans[[i]] = data
  }
  ans
}

findsubmodel.RF2 = function(resp, tms) {
  res = paste(resp, "~1", sep = "")
  for (fix in tms) {
    res = as.vector(sapply(res, function(x) paste(x, c("", paste("+", fix, sep = "")), sep = "")))
  }
  lapply(res, as.formula)
}

