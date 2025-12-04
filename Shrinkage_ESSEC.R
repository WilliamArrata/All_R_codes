
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   ################

require("pacman")
pacman::p_load("tseries","readxl","dplyr", "tidyr", "data.table", "ggplot2")

#####################   DATA DOWNLOAD AND COMPUTATION OF EXPECTED RETURNS AND COVARIANCES   ################

#I load the data
returns <- as.matrix(read_excel("stock_prices.xlsx") %>%  select_if(is.numeric) %>%  mutate_all(~ ( (.) - shift(.))/(.)) %>% 
                       na.omit() %>% rename_with(~gsub(" Equity","", (.)) ))      #daily historical returns
mean <- 252*matrix(colMeans(returns))                                                #annualized expected returns
sig <- 252*cov(returns)                                                              #annualized covariances


################################# SIGN CONSTRAINED EFFICIENT FRONTIER   #####################################

EF = function (returns, nports, shorts, wmax){
  max_ret<-max(mean)
  #max_ret<-(1+as.numeric(shorts)*0.5)*max(mean)      #la cible de renta maximale
  target <- seq(-max_ret, max_ret, len= nports)       #on définit les cibles de renta via nports et maxret
  reslow <- rep(-as.numeric(shorts), length(mean))    #vecteur de poids minimum
  reshigh <- rep(wmax,length(mean))                   #vecteur de poids max
  output <- list()
  for (i in seq_along(target)){
    sol <- NULL
    try(sol <- portfolio.optim(returns, pm = target[i]/252, reshigh = reshigh, reslow = reslow, shorts = shorts), silent = T)
    if(!is.null(sol)){
      output[[i]] <- c(i, sqrt(252)*sol$ps, 252*sol$pm, sol$pw)
      names(output[[i]]) <- c("i", "vol", "return", paste0("w", 1:length(mean)))}
  }
  output<-as.data.frame(do.call(rbind,output))
  rownames(output)<-output$i
  return(output)
}

nports <- 300   #nb of ptf, thus we have 300 target expected returns

#Efficient frontier when short selling is forbidden
shorts <- F
wmax <- 1

ptfs_no_s <- EF(returns = returns, nports = nports, shorts = shorts, wmax = wmax)
low_no_s <- which.min(ptfs_no_s$vol)
high_no_s <- which.max(ptfs_no_s$return)
effi_no_s <- ptfs_no_s[low_no_s:high_no_s,]


#######################################   RESAMPLING HISTORICAL RETURNS   #####################################

require(epca)

shrink_cov <- shrinkage(sig, gamma = 2, shrink="soft", epsilon = 1e-11)

#Simulating n_samp samples of length n_tirages for the 6 assets
require(MASS)
set.seed(33)
n_tirages <- nrow(returns)                                       #length of each sample

#daily simulated returns for the six stocks
estim <- mvrnorm(n_tirages, mean/252, shrink_cov/252, tol = 1e-06, empirical = F)

###########################   EFFICIENT FRONTIER WITH SHRUNK ESTIMATE   ##############################

#I run the optimization for the 1000 simulated sets of returns
ptfs_no_s_shrunk <- EF(returns = estim, nports = nports, shorts = shorts, wmax = wmax)
low_no_s_shrunk <- which.min(ptfs_no_s_shrunk$vol)
high_no_s_shrunk <- which.max(ptfs_no_s_shrunk$return)
effi_no_s_shrunk <- ptfs_no_s_shrunk[low_no_s_shrunk:high_no_s_shrunk,]


#graph of the Markowitz frontier and the average of the resampled efficient frontiers
col<-c("darkblue","indianred")
par(mar=c(7, 6, 4, 4),xpd=T)
plot(100*ptfs_no_s$vol, 100*ptfs_no_s$return,col="darkblue", lwd=2,xlab="standard deviation (%)",
     ylab="expected return (%)",las=1, type="l",pch=20,ylim=100*range(c(ptfs_no_s$return, ptfs_no_s$return), na.rm = T),
     xlim=100*range(c(ptfs_no_s$vol, ptfs_no_s_shrunk$vol), na.rm = T))
lines(100*ptfs_no_s_shrunk$vol, 100*ptfs_no_s_shrunk$return,col="indianred", lwd=2)
legend("bottom", horiz = T,inset = c(0,-0.4),text.col=col,pch=rep(NA,3),lty=rep(1,3),col=col, bty="n",
       legend= c("Markowitz efficient frontier","Efficient frontier with shrunk covariance matrix"))

#Graph
weights <- ptfs_no_s_shrunk[,-c(1:3)]
at_1 = seq(0, 1, 0.25)
at_2 = seq(1, ncol(weights), length.out = 7)
colvector <- rainbow(6)

cex <- 0.8
par(mar = c(8,4,4,4) + 0.1, xpd = T, cex.axis = cex)
for (i in 1:nrow(weights)){
  plot(1:ncol(weights),weights[1+nrow(weights)-i,], xlab="",ylab="", ylim=0:1,
       xlim=c(0,ncol(weights)),las=1, col=colvector[i],pch=20, axes=F)
  polygon(c(1:ncol(weights),ncol(weights):1), c(rep(0,ncol(weights)),rev(weights[1+nrow(weights)-i,])),
          col=colvector[i])
  par(new=T)}
axis(1, at=at_2, labels=round(100*ptfs_no_s_shrunk$return,1)[at_2], cex.axis =cex)
axis(2, at=at_1, labels=at_1*100,cex.axis=cex)
mapply(mtext, c("expected return (%)", "weights (%)"), side=c(1,2), line = rep(2.5,2))
legend("bottom",ncol=3,inset = c(0,-0.5),legend=rev(colnames(returns)),text.col=colvector,col=colvector,
       pch=c(15), bty="n")
box()