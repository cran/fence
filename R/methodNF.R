#' Fence model selection (Nonparametric Model)
#'
#' Fence model selection (Noparametric Model)

#' @param full formula of full model
#' @param data data
#' @param spline variable needed for spline terms
#' @param ps order of power
#' @param qs number of knots
#' @param B number of bootstrap sample, parametric for lmer
#' @param grid grid for c
#' @param bandwidth bandwidth for kernel smooth function
#' @param lambda A grid of lambda values
#' @return 
#' \item{models}{list all model candidates with p polynomial degrees and q knots in the model space}
#' \item{Qd_matrix}{list a matrix of QM - QM.tilde for all model candidates. Each row is for each bootrap sample}
#' \item{bandwidth}{list the value of bandwidth}
#' \item{model_mat}{list a matrix of selected models at each c values in grid (in columns). Each row is for each bootstrap sample}
#' \item{freq_mat}{list a matrix of coverage probabilities (frequency/smooth_frequency) of each selected models for a given c value (index)}
#' \item{c}{list the adaptive choice of c value from which the parsimonious model is selected}
#' \item{lambda}{penalty (or smoothing) parameter estimate given selected p and q}
#' \item{sel_model}{list the selected (parsimonious) model given the adaptive c value}
#' \item{beta.est.u}{A list of coefficient estimates given a lambda value}
#' \item{f.x.hat}{A vector of fitted values obtained from a given lambda value and beta.est.u}
#'  @note The current Fence method in Nonparametric model focuses on one spline variable. 
#'  This method can be extended to a general case with more than one spline variables, and includes non-spline variables.
#' @author Jiming Jiang Jianyang Zhao J. Sunil Rao Bao-Qui Tran Thuan Nguyen
#' @references 
#' \itemize{
#'  \item{Jiang J., Rao J.S., Gu Z., Nguyen T. (2008),  Fence Methods for Mixed Model Selection. The Annals of Statistics, 36(4): 1669-1692}
#'  \item{Jiang J., Nguyen T., Rao J.S. (2009), A Simplified Adaptive Fence Procedure. Statistics and Probability Letters, 79, 625-629}
#'  \item{Jiang J., Nguyen T., Rao J.S. (2010), Fence Method for Nonparametric Small Area Estimation. Survey Methodology, 36, 1, 3-11}
#' }
#' 
#' @examples
#' \dontrun{
#' require(fence)
#' n = 100
#' set.seed(1234)
#' x=runif(n,0,3)
#' y = 1-x+x^2- 2*(x-1)^2*(x>1) + 2*(x-2)^2*(x>2) + rnorm(n,sd=.2)
#' lambda=exp((c(1:60)-30)/3)
#' data=data.frame(cbind(x,y))   
#' test_NF = fence.NF(full=y~x, data=data, spline='x', ps=c(1:3), qs=c(2,5), B=1000, lambda=lambda)
#' plot(test_NF)
#' summary <- summary(test_NF) 
#' model_sel <- summary[[1]]
#' model_sel
#' lambda_sel <- summary[[2]]
#' lambda_sel
#' }
#' @export

fence.NF = function(full, data, spline, ps = 1:3, qs = NA, B = 100, grid = 101, bandwidth = NA, lambda) {
  n = nrow(data)
  if (all(is.na(qs))) {
    lower = ceiling(n/5)
    upper = floor(n/4) 
    if (n < 50) {
      lower = 2
    }
    if (n > 500) {
      lower = 100
      upper = 125
    }
    qs = lower:upper
  }
  
  y = matrix(data[,as.character(full)[2]], ncol = 1)
  baseX = model.matrix(full, data)
  nX= ncol(baseX)-1
  #### STOP IF MORE THAN 1 COVARIATES ###
  if (nX >= 2) {
    stop(paste0("Model input has ", ncol(baseX)-1, ' predictors. Current fence package can only considers 1 predictor at a time.'))
  }
  ####
  
  splineX=which(colnames(baseX)[-1] %in% spline)
  
  expand.grid.X = function(nX,ps, qs){
    psX=replicate(nX, c(0,ps), simplify=FALSE)
    qsX=replicate(nX, c(0,qs), simplify=FALSE)
    qpX = c(psX, qsX)
    model.cand <- expand.grid(qpX)
    colnames(model.cand) <- paste(c(rep('p',nX),rep('q',nX)), rep(1:nX, times=2), sep='')
    for (i in c(1:nX)){
      #remove candidate models where p=0 and q>0
      model.cand <- model.cand[!(model.cand[,i]==0 & model.cand[,nX+i]>0),]
      #remove candidate models where q>0 for X's not in spline list
      if(!i %in% splineX) model.cand <- model.cand[model.cand[,nX+i]==0,]
    } 
    rownames(model.cand) <- NULL
    return(model.cand)
  }
  
  model.cand <- expand.grid.X(nX, ps, qs)
  
  #Get X^p and knots for all p,q,x
  addtionXlist =apply(as.matrix(baseX[,-1],nrow=n),2, genaddX, ps=ps, qs=qs)
  
  # Get matrix of X's and Z's for each model with known p and q
  make.eachX = function(addtionXlist, xindex, p, q){
    psm = addtionXlist[[xindex]]$pmatrix #matrix where each column is X^p for different p
    qsm = addtionXlist[[xindex]]$qmatrix #list contaning knots column for different number of knots    
    if(p==0 & q==0) Xall=as.matrix(baseX[,1], ncol=1); Zall=list(NULL)
    if(p==1) Xall = baseX[,c(1,xindex+1)]
    if(p>1) Xall=cbind(baseX[,c(1,xindex+1)], psm[,2:p])
    if(q==0) Zall = list(NULL)
    if(q>0) {
      whichq=which(qs==q)
      Zall = qsm[[p]][[whichq]]
    }
    return(list('X'=Xall,'Z'= Zall))
  }
  
  #X and Z list for all candidate models
  Xall=Zall=list()
  for (i in 1:nrow(model.cand)){
    m=model.cand[i,]
    #get X's columns and Z's column of each x covariate
    eachmod=lapply(1:nX, FUN=function(xi) {make.eachX(addtionXlist, xi, p=as.numeric(m[xi]), q=as.numeric(m[nX+xi]))})
    # extract X's columns for each x and combine into a matrix
    Xmat=matrix(unlist(lapply(eachmod, function(eachXls){eachXls$X})), nrow=n, byrow = FALSE)
    if(nX>1) Xmat=Xmat[,-(which(colSums(Xmat==1)==n)[-1])] #remove repeated intercept-column of 1's
    # extract Z's columns for each x and combine into a matrix
    Zmat=unlist(lapply(eachmod, function(eachXls){eachXls$Z}))
    if(length(Zmat)==0) {Zmat=list()} else{Zmat = matrix(Zmat, nrow=n, byrow = FALSE)} 
    
    Xall[[i]]=Xmat
    Zall[[i]]=Zmat
  }
  Wall=mapply(cbind,Xall, Zall)
  
  ##Calculate lack-of-fit  
  Q.comb=sapply(Wall, function(W){
    W=matrix(unlist(W),nrow=n)
    temp.b=solve(t(W)%*%W)%*%t(W)%*%y
    y.hat = W%*%temp.b
    list(t(y - y.hat)%*%(y - y.hat), temp.b)
  })	
  
  Q.hat=unlist(Q.comb[1,])
  beta.hat=Q.comb[2,]
  Qdiff = Q.hat - min(Q.hat)
  
  # Chose model with minimum lack-of-fit  
  index.min = which.min(Q.hat)
  model.min = model.cand[index.min,]
  
  # Bootstrap based on chosen model (not full model) 
  bs = genNFbs(B, y, Wall[[index.min]], beta.hat[[index.min]]) 
  
  # Recalculate lack-of-fit for all models & bootstraps  
  Q.star=sapply(Wall, function(W){
    W=matrix(unlist(W),nrow=n)
    temp.b=solve(t(W)%*%W)%*%t(W)%*%bs
    y.hat = W%*%temp.b
    colSums((bs-y.hat)^2)
  })	
  
  Qdiff.star = sweep(Q.star, 1, apply(Q.star, 1, min), '-') #B x number of modelcandidate
  
  # Find peak Cn  
  cs = seq(0, max(Qdiff.star), length.out = grid)
  if (is.na(bandwidth)) {
    bandwidth = (cs[2] - cs[1]) * 3
  }
  
  # Optimal model < Cn
  model.search=function(model.cand, infence){
    m.infence=model.cand[infence,]
    index.infence=which(infence)
    sump=apply(m.infence, 1, function(r){sum(r[1:nX])})
    sumq=apply(m.infence, 1, function(r){sum(r[(nX+1):(nX+nX)])})
    return(min(index.infence[order(sumq, sump)]))
  }
  model_mat = matrix(NA, nrow = B, ncol = grid)
  for (i in 1:length(cs)) {
    infence_matrix = Qdiff.star <= cs[i]
    model_mat[, i] = apply(infence_matrix, 1, function(x){model.search(model.cand,x)})
  }
  
  freq_mat = apply(model_mat, 2, function(l) {
    tab = sort(table(l), decreasing = TRUE)
    c(as.numeric(names(tab)[1]), tab[1])
  })
  freq_mat[2,] = freq_mat[2,] / B
  freq_mat = rbind(freq_mat, ksmooth(cs, freq_mat[2,], kernel = "normal", bandwidth = bandwidth, x.points = cs)$y)
  colnames(freq_mat) = cs
  rownames(freq_mat) = c("index", "frequency", "smooth_frequency")
  
  cindex = peakw(cs, freq_mat[3,], 2)
  Cn.peak = cs[cindex]
  
  # Find p and q that give Qdiff < Cn.peak
  index.infence=which(Qdiff < Cn.peak)
  correct.m=model.cand[index.infence,]
  sump.m=apply(correct.m, 1, function(r){sum(r[1:nX])})
  sumq.m=apply(correct.m, 1, function(r){sum(r[(nX+1):(nX+nX)])})
  index.m=min(index.infence[order(sumq.m, sump.m)])
             
  # Get y fit
  get.y.fit.out = get.y.fit(index.ind=index.m,beta.hat,X=Xall,Z=Zall,y,Q.hat,index=index.min, Cn.ind=Cn.peak, lambda, model.cand)
  lambda.fit = get.y.fit.out[[1]]
  beta.est.u = get.y.fit.out[[2]]
  f.x.hat = get.y.fit.out[[3]]
  
  ans = list(full = full, models = model.cand, 
             model_mat =model_mat, freq_mat =freq_mat,
             Qd_matrix = Qdiff,  bandwidth =bandwidth,beta.est.u=beta.est.u,  
             f.x.hat =f.x.hat, sel_model=model.cand[index.m,], c=Cn.peak, lambda = lambda.fit)
  class(ans) = "NF"
  return(ans)
}

############# HELPER FUNCTIONS #######################
genaddX = function(svalue, ps, qs) {
  psm = outer(svalue, ps, '^')
  qsm = lapply(ps, function(p) {
    res = lapply(qs, function(q) {
      knots = cover.design(R=as.matrix(svalue, ncol = 1), nd=q)
      drop(outer(svalue, knots$design, function(x, y) ifelse(x < y, 0, x - y)) ^ p)
    })
    names(res) = paste0("q", qs)
    res
  })
  names(qsm) = paste0("p", ps)
  list(pmatrix = psm, qmatrix = qsm)
}

genNFbs = function(B, y, W0, beta0) {
  sigma = sqrt(sum((W0 %*% beta0 - y)^2) / (nrow(W0) - ncol(W0)))
  ybase = W0 %*% beta0
  matrix(as.vector(ybase) + rnorm(B * nrow(W0), 0, sigma), nrow = nrow(W0))
}


Q.cal.fit = function(lambda.ind,temp.x,temp.z, y,len.q){
  para.est = para.fit(lambda.ind,temp.x,temp.z,y,len.q)
  x.beta.fit = para.est[[1]]
  u.est.fit = para.est[[2]]
  temp.e = y - temp.x%*%x.beta.fit - temp.z%*%u.est.fit
  t(temp.e)%*%temp.e
}

lambda.cal = function(k.index){
  exp((k.index-30)/3)
}

para.fit = function(lambda.est,temp.x,temp.z,y,len.q){
  n <- NULL
  rm(n)
  v.est = diag(n) + lambda.est^(-1)*(temp.z%*%t(temp.z))
  v.est.inv = solve(v.est)
  x.v.inv = t(temp.x)%*%v.est.inv
  beta.est.fit = solve(x.v.inv%*%temp.x)%*%x.v.inv%*%y
  x.beta.fit = temp.x%*%beta.est.fit
  u.est.fit = lambda.est^(-1)*solve(diag(len.q)+ lambda.est^(-1)*(t(temp.z)%*%temp.z))%*%t(temp.z)%*%(y - x.beta.fit)
  list(beta.est.fit,u.est.fit)
}

get.y.fit = function(index.ind,beta.hat,X,Z,y,Q.hat,index,Cn.ind,lambda,model.cand){
  model.m=model.cand[index.ind,]
  beta.m = beta.hat[[index.ind]]
  nX=ncol(model.cand)/2
  len.q = sum(model.m[(nX+1):(nX+nX)])
  
  if(len.q==0){
    f.x.hat = X[[index.ind]]%*%beta.m
    beta.est.u = list(beta.m, NULL)  
    lambda.fit = 'no knots'
  }else{
    temp.x = X[[index.ind]]
    temp.z = Z[[index.ind]]
    
    Q.hat.est = rep(NA,length(lambda))
    for(k in 1:length(lambda)){
      Q.hat.est[k] = Q.cal.fit(lambda[k],temp.x,temp.z,y,len.q)
      
    }
    Q.diff.est = Q.hat.est - Q.hat[index]
    k.index = which(Q.diff.est <= Cn.ind)
    k.index = k.index[length(k.index)]
    
    lambda.u = lambda.cal(k.index+1)
    lambda.l = lambda.cal(k.index)
    lambda.diff = lambda.u - lambda.l
    
    error = lambda.diff - 10^(-6)
    
    while(error > 0){
      mid.lambda = lambda.diff/2  + lambda.l
      Q.est = Q.cal.fit(mid.lambda,temp.x,temp.z,y,len.q)
      Q.tilte = Q.hat[index]
      Q.diff.fit = Q.est - Q.tilte
      if(Q.diff.fit <= Cn.ind){
        lambda.l = mid.lambda
        lambda.diff = lambda.u - lambda.l
      }else{
        lambda.u = mid.lambda
        lambda.diff = lambda.u - lambda.l				
      }
      error = lambda.diff - 10^(-6)
    }
    
    lambda.fit = lambda.l
    beta.est.u = para.fit(lambda.fit,temp.x,temp.z,y,len.q)
    f.x.hat = temp.x%*%beta.est.u[[1]]+temp.z%*%beta.est.u[[2]]
  }
  list(lambda.fit,beta.est.u,f.x.hat)
}

#' Plot Nonparametric Fence model selection
#'
#' @param x Object to be plotted
#' @param ... Additional arguments. CNot currently used.
#' @export
#' @importFrom utils globalVariables
plot.NF = function(x = res, ...) {
  res <- sp <- m <- NULL
  rm(res,sp,m)
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

#' Summary Nonparametric Fence model selection
#' 
#' @param object Object to be summarized
#' @param ... addition arguments. Not currently used
#' @export
summary.NF = function(object = res, ...) {
  res <- NULL
  rm(res)
  print(list(object$sel_model, object$lambda))
}
