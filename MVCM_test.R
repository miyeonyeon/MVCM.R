#############################################################################################
## This is a test code to implement Zhu's (2010) method of FMVCM. Note here we use FA and MD.
#############################################################################################

root <- "/MVCM.R/"
devtools::load_all(root, quiet = TRUE)
devtools::load_all()

library(R.matlab, warn.conflicts = FALSE) # For loading mat files

tractdata <- readMat("Data/tractdata.mat")$tractdata # L0 x 3 matrix
designdata <- readMat("Data/designdata.mat")$designdata # n x p matrix
m <- dim(readMat("Data/diffusiondata.mat")$diffusionFiles)[1]
diffusionFiles <- vector("list", m) # m x 1 list; L0 x n matrix
diffusionFiles[[1]] <- readMat("Data/diffusiondata.mat")$diffusionFiles[[1]][[1]]
diffusionFiles[[2]] <- readMat("Data/diffusiondata.mat")$diffusionFiles[[2]][[1]]


## ---- Analysis Flu Study
result0 <- MVCM_read(tractdata, designdata, diffusionFiles, nofeatures=m)
NoSetup <- result0$NoSetup
arclength <- result0$arclength
Xdesign <- result0$Xdesign
Ydesign <- result0$Ydesign

(n <- NoSetup[1])
(L0 <- NoSetup[2])
(p <- NoSetup[3])
# m <- NoSetup[4]


## plot tract
Mnames <- list('FA', 'MD')

for (mii in 1:m) {
  png(file=sprintf("Output/%s_tract.png", Mnames[[mii]]), height=300*2, width=300*2, res=100)
  par(mar=c(5,5,4,2)+0.1)
  plot(NULL, xlim=c(0, length(arclength)-1), ylim=range(Ydesign), xlab='arc-length', ylab=Mnames[[mii]], type='n')

  for (nii in 1:n) {
    lines(arclength, Ydesign[nii, , mii], col='red', lwd=2)
  }

  dev.off()
}


## ---- 2.1 estimate betas using local polynomial kernel smoothing
result1 <- MVCM_lpks_wob(NoSetup, arclength, Xdesign, Ydesign)
mh <- result1$mh

result2 <- MVCM_lpks_wb1(NoSetup, arclength, Xdesign, Ydesign, mh)
efitBetas <- result2$efitBetas


## plot betas
for (mii in 1:m) {
  png(fil=sprintf("Output/%s_beta.png", Mnames[[mii]]), height=300*2, width=300*2, res=100)
  par(mar=c(5,5,4,2)+0.1)
  plot(arclength, efitBetas[1, , mii], type='l', col='blue', lwd=2, xlab='arc-length', ylab=Mnames[[mii]], xlim=c(0, length(arclength)-1), ylim=range(efitBetas))
    lines(arclength, efitBetas[2, , mii], col='red', lwd=2)
    lines(arclength, efitBetas[3, , mii], col='green', lwd=2)
    # lines(arclength, efitBetas[4, , mii], col='black', lwd=2)
    # legend('topleft', legend=c('Intercept', 'Gender', 'Flu', 'Age'), col=c('blue', 'red', 'green', 'black'), lwd=2, bty='n', cex=1, text.font=2)
    legend('topleft', legend=c('Intercept', 'Gender', 'Age'), col=c('blue', 'red', 'green'), lwd=2, bty='n', cex=1, text.font=2)

  dev.off()
}


## ---- 2.2 smoothing individual function
ResYdesign <- Ydesign-result2$efitYdesign
result3 <- MVCM_sif(arclength, ResYdesign)
