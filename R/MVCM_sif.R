#' MVCM_sif is to implement Zhu's (2010) method of smoothing individual function without preselected bandwidth in MVCM
#'
#' Input:
#'     arclength   - col vector of the arclength from one end to the other end
#'     ResYdesign  - a n x L0 x m matrix of difference between fiber bundle diffusion properties and fitted fiber bundle diffusion properties
#'     kstr        - Kernel function
#' Output:
#'     ResEtas     - a n x L0 x m matrix of of difference between ResYdesign and fitted eta
#'     efitEtas    - a n x L0 x m matrix of estimated eta
#'     eSigEta     - a m x m x L0 x L0 matrix of covariance matrix of etas
#' ###################################################################################################################
#' Please run
#'     [efitBetas,InvSigmats,efitYdesign] = MVCM_lpks_wb1(NoSetup,arclength,Xdesign,Ydesign,mh,kstr)
#' before you use MVCM_sif
#' ###################################################################################################################
#' March 28, 2010 @ AA
#'
#' @references Hongtu Zhu. Runze Li. Linglong Kong. "Multivariate varying coefficient model for functional responses." Ann. Statist. 40 (5) 2634 - 2666, October 2012. https://doi.org/10.1214/12-AOS1045
#' @import pracma
#' @export

MVCM_sif <- function(arclength, ResYdesign, kstr='exp(-.5*t^2)'){

  n <- dim(ResYdesign)[1]
  L0 <- dim(ResYdesign)[2]
  m <- dim(ResYdesign)[3]
  xrange <- base::max(arclength) - base::min(arclength)
  nh <- base::max(30, floor(L0 / 2))
  hmin <- 1.0 * xrange / L0
  hmax <- xrange / 8
  vh <- pracma::logspace(log10(hmin), log10(hmax), nh)
  Tmat0 <- arclength %*% matrix(1, 1, L0) - matrix(1, L0, 1) %*% t(arclength)
  GCVs <- matrix(0, nh, m)

  for (nhii in 1:nh){
    h <- vh[nhii]
    Tmat <- Tmat0 / h
    t <- Tmat
    Kmat <- eval(parse(text = kstr)) / h
    S0 <- matrix(1, L0, 1) %*% matrix(apply(Kmat, 2, sum), ncol=dim(Kmat)[2])
    S1 <- matrix(1, L0, 1) %*% matrix(apply(Kmat * Tmat0, 2, sum), ncol=dim(Kmat * Tmat0)[2])
    S2 <- matrix(1, L0, 1) %*% matrix(apply(Kmat * Tmat0^2, 2, sum), ncol=dim(Kmat * Tmat0^2)[2])
    Smat <- t(Kmat * (S2 - S1 * Tmat0) / (S2 * S0 - S1^2))
    for (mii in 1:m){
      GCVs[nhii, mii] <- sum(sum((t(ResYdesign[, , mii]) - Smat %*% t(ResYdesign[, , mii]))^2)) / (1 - sum(diag(Smat)) / L0)^2 / n
    }
  }

  flag <- apply(GCVs, 2, which.min)
  mh <- vh[flag]

  efitEtas <- array(0, c(n, L0, m))
  ResEtas <- array(0, c(n, L0, m))
  for (mii in 1:m){
    h <- mh[mii]
    Tmat <- Tmat0 / h
    t <- Tmat
    Kmat <- eval(parse(text = kstr)) / h
    S0 <- matrix(1, L0, 1) %*% matrix(apply(Kmat, 2, sum), ncol=dim(Kmat)[2])
    S1 <- matrix(1, L0, 1) %*% matrix(apply(Kmat * Tmat0, 2, sum), ncol=dim(Kmat * Tmat0)[2])
    S2 <- matrix(1, L0, 1) %*% matrix(apply(Kmat * Tmat0^2, 2, sum), ncol=dim(Kmat * Tmat0^2)[2])
    Smat <- Kmat * (S2 - S1 * Tmat0) / (S2 * S0 - S1^2)
    efitEtas[, , mii] <- ResYdesign[, , mii] %*% Smat
    ResEtas[, , mii] <- ResYdesign[, , mii] - efitEtas[, , mii]
  }

  eSigEta <- array(0, c(m, m, L0, L0))
  for (L0ii in 1:L0){
    for (L0jj in 1:L0){
      eSigEta[, , L0ii, L0jj] <- t(efitEtas[, L0ii, ]) %*% efitEtas[, L0jj, ] / n
    }
  }

  return(list(ResEtas=ResEtas, efitEtas=efitEtas, eSigEta=eSigEta))
}
