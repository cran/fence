% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methodAF.fh.R
\name{adaptivefence.fh}
\alias{adaptivefence.fh}
\title{Adaptive Fence model selection (Small Area Estmation)}
\usage{
adaptivefence.fh(mf, f, ms, d, lf, pf, bs, grid = 101, bandwidth, method)
}
\arguments{
\item{mf}{Call function, for example: default calls: function(m, b) eblupFH(formula = m, vardir = D, data = b, method = "FH")}

\item{f}{Full Model}

\item{ms}{find candidate model, findsubmodel.fh(full)}

\item{d}{Dimension number}

\item{lf}{Measures lack of fit using function(res) -res$fit$goodness[1]}

\item{pf}{Dimensions of model}

\item{bs}{Bootstrap}

\item{grid}{grid for c}

\item{bandwidth}{bandwidth for kernel smooth function}

\item{method}{Method to be used. Fay-Herriot method is the default.}
}
\value{
\item{models}{list all model candidates in the model space}
\item{B}{list the number of bootstrap samples that have been used}
\item{lack_of_fit_matrix}{list a matrix of Qs for all model candidates (in columns). Each row is for each bootstrap sample}
\item{Qd_matrix}{list a matrix of QM - QM.tilde for all model candidates. Each row is for each bootrap sample}
\item{bandwidth}{list the value of bandwidth}
\item{model_mat}{list a matrix of selected models at each c values in grid (in columns). Each row is for each bootstrap sample}
\item{freq_mat}{list a matrix of coverage probabilities (frequency/smooth_frequency) of each selected models for a given c value (index)}
\item{c}{list the adaptive choice of c value from which the parsimonious model is selected}
\item{sel_model}{list the selected (parsimonious) model given the adaptive c value}
}
\description{
Adaptive Fence model selection (Small Area Estmation)
}
\details{
In Jiang et. al (2008), the adaptive c value is chosen from the highest peak in the p* vs. c plot.  
In Jiang et. al (2009), 95\% CI is taken into account while choosing such an adaptive choice of c.
In Thuan Nguyen et. al (2014), the adaptive c value is chosen from the first peak. This approach works better in the 
moderate sample size or weak signal situations.  Empirically, the first peak becomes highest peak when sample size 
increases or signals become stronger
}
\note{
\itemize{
\item{The current Fence package focuses on variable selection. 
 However, Fence methods can be used to select other parameters of interest, e.g., tunning parameter, variance-covariance structure, etc.}
\item{The number of bootstrap samples is suggested to be increased, e.g., B=1000 when the sample size is small, or signals are weak}
}
}
\examples{
\dontrun{
require(fence)
### example 1 ####
data("kidney")
data = kidney[-which.max(kidney$x),]     # Delete a suspicious data point #
data$x2 = data$x^2
data$x3 = data$x^3
data$x4 = data$x^4
data$D = data$sqrt.D.^2
plot(data$y ~ data$x)
full = y~x+x2+x3+x4
testfh = fence.sae(full, data, B=1000, fence="adaptive", method="F-H", D = D)
testfh$sel_model
testfh$c
}
}
\references{
\itemize{
 \item{Jiang J., Rao J.S., Gu Z., Nguyen T. (2008),  Fence Methods for Mixed Model Selection. The Annals of Statistics, 36(4): 1669-1692}
 \item{Jiang J., Nguyen T., Rao J.S. (2009), A Simplified Adaptive Fence Procedure. Statistics and Probability Letters, 79, 625-629}
 \item{Thuan Nguyen, Jie Peng, Jiming Jiang (2014), Fence Methods for Backcross Experiments.  Statistical Computation and Simulation, 84(3), 644-662}
}
}
\author{
Jiming Jiang  Jianyang Zhao  J. Sunil Rao  Thuan Nguyen
}
