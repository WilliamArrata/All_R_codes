
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


#Simulating n_samp samples of length n_tirages for the 6 assets
require(MASS)
set.seed(33)
n_samp <- 1000                                                       #number of samples
n_tirages <- 250                                                     #length of each sample
#daily simulated returns for the six stocks in each simu
estim <- replicate(n_samp , mvrnorm(n_tirages, mean, 252*sig, tol = 1e-06, empirical = F)/252)
resampm <- colMeans(estim, dims = 1)                                #the average return for each asset in each simu

#Distribution of daily returns for 3 random simulations for a given asset
alea <- sort(sample(n_samp, 3))                                  
dens_ex <- 10000*data.frame(estim[, 6, alea])/n_samp                   #I pick 3 simulations randomly for asset 6
dens <- apply(dens_ex, 2, density)                                   #I represent their return distrib with density

#Graph
par(mar=c(7,5,4,3),xpd=T)
plot(NA, xlim = range(sapply(dens, "[", "x")), ylim = range(sapply(dens, "[", "y")), xlab = "in %",
     ylab = "frequency of observation (in %)", main = "A few resampled distributions of Hermès daily stock returns")
mapply(lines, dens, col = 1:length(dens))
legend("bottom", horiz = T,inset = c(0, -0.3), text.col = 1:length(dens), pch = c(NA, NA), lty = rep(1, 3),
       col = 1:length(dens), bty = "n",legend = paste("simulation", alea))

#Distribution of mean returns across all simu for a given asset class
dens_moy <- density( 100*resampm[6, ])
dens_moy$y <- 100*dens_moy$y/n_samp

#Graph
plot(NA, xlim = range(dens_moy$x), ylim = range(dens_moy$y), xlab = "in %", ylab = "frequency of observation (in %)",
     col = "darkblue", main = "Distribution of expected returns of Hermès daily stock returns")
lines(dens_moy)

#######################################   RESAMPLED EFFICIENT FRONTIER   #####################################

#I write a new function which has resampled series of returns as input
EF2 = function (nports, shorts, wmax){
  return(EF( returns = estim[,,i], nports, shorts, wmax))}

#I run the optimization for the 1000 simulated sets of returns
resampw <- list()
for (i in 1:n_samp){
  output <- EF2(nports = nports, shorts = shorts, wmax = wmax)
  if (nrow(output) == 0){                                #weights are all equal to 0 when no solution
    output <- matrix(c(0, rep(NA, 2), rep(0, ncol(returns))), nports, 3 + ncol(returns), byrow = T,
                   dimnames = list((1:nports), c("i", "vol", "return", paste0("w", 1:ncol(returns)))))}
  else {
    output <- rbind(matrix( c(0, rep(NA, 2), rep(0, ncol(returns))), nports - nrow(output), ncol(output), byrow = T,
                         dimnames = list((1:nports)[-output$i], colnames(output))), output)}
  output <- output[order(as.numeric(rownames(output))),]
  resampw[[i]] <- output[ ,grep("w", colnames(output))]
}

#I average weights across each target return for all optimizations
aveweight <- as.matrix(Reduce("+", resampw)/n_samp)
#I rescale for the number of simulations where the target return not reached and weights are 0 by default
aveweight <- aveweight/replicate(ncol(aveweight), rowSums(aveweight))

#I apply average weights to initial parameters to get robust efficient frontier
resamp <- data.frame(vol = diag(sqrt(aveweight%*%sig%*%t(aveweight))), return = aveweight%*%mean)

#graph of the Markowitz frontier and the average of the resampled efficient frontiers
col<-c("darkblue","indianred")
par(mar=c(7, 6, 4, 4),xpd=T)
plot(100*ptfs_no_s$vol,100*ptfs_no_s$return,col="darkblue", lwd=2,xlab="standard deviation (%)",
     ylab="expected return (%)",las=1, type="l",pch=20,ylim=100*range(c(resamp$return,ptfs_no_s$return), na.rm = T),
     xlim=100*range(c(resamp$vol,ptfs_no_s$vol), na.rm = T))
lines(100*resamp$vol,100*resamp$return,col="indianred", lwd=2)
legend("bottom", horiz = T,inset = c(0,-0.4),text.col=col,pch=rep(NA,3),lty=rep(1,3),col=col, bty="n",
       legend= c("Markowitz efficient frontier","resampled efficient frontier"))

#weights by asset across the frontier (transition map)
cum_ave_w <- apply(aveweight, 1, cumsum)

#Graph
at_1 = seq(0, 1, 0.25)
at_2 = seq(1, ncol(cum_ave_w), length.out = 7)
colvector <- rainbow(6)

cex <- 0.8
par(mar = c(8,4,4,4) + 0.1, xpd = T, cex.axis = cex)
for (i in 1:nrow(cum_ave_w)){
  plot(1:ncol(cum_ave_w),cum_ave_w[1+nrow(cum_ave_w)-i,], xlab="",ylab="", ylim=0:1,
       xlim=c(0,ncol(cum_ave_w)),las=1, col=colvector[i],pch=20, axes=F)
  polygon(c(1:ncol(cum_ave_w),ncol(cum_ave_w):1), c(rep(0,ncol(cum_ave_w)),rev(cum_ave_w[1+nrow(cum_ave_w)-i,])),
          col=colvector[i])
  par(new=T)}
axis(1, at=at_2, labels=round(100*resamp$return,1)[at_2], cex.axis =cex)
axis(2, at=at_1, labels=at_1*100,cex.axis=cex)
mapply(mtext, c("expected return (%)", "weights (%)"), side=c(1,2), line = rep(2.5,2))
legend("bottom",ncol=3,inset = c(0,-0.5),legend=rev(colnames(returns)),text.col=colvector,col=colvector,
       pch=c(15), bty="n")
box()



#######################################   BLACK LITTERMAN MODEL   #####################################

# relative weights of assets in benchmark portfolio
w_m <- read_excel("stock_prices_2.xlsx",2) %>% mutate_if(is.character, as.numeric) %>% na.omit() %>% 
  select_if(is.numeric) %>% mutate(across()/rowSums(across()))

#expected return on the market portfolio
return_m <- read_excel("stock_prices_2.xlsx",3) %>% rename(date="...1") %>% 
  mutate_if(is.character, as.numeric) %>% mutate(date = as.Date(date, origin = "1899-12-30")) %>%
  filter(date >= as.Date("2019-01-01")) %>%                                               #load stock prices
  mutate_if(is.numeric, ~ ( (.) - shift(.))/(.)) %>% na.omit() %>% select_if(is.numeric)  #daily historical returns
r_m <- 252*colMeans(return_m)
var_m <- 252*var(return_m)

r_f <- 0.01      #the riskfree rate

#Participation Matrix: Vue relative sur MC vs OR; vue absolue sur Bouygues
P <- matrix(c(-1,0,1,0,0,0,0,1,0,0,0,0), ncol = ncol(return), byrow=T)

#Vector of subjective expected returns means: MC FP va surperformer OR FP de 7%; Bouygues va avoir une perf de -5%
Q <- c(0.20,-0.20)

Omega <- matrix(c(0.005,0,0,0.02), nrow = length(Q))      #Matrix of views uncertainties
tau <- 0.025    #shrinkage coefficient

lambda_m <- (r_m - r_f)/var_m                  #market risk aversion coefficient: expected return on mkt ptf in excess of rf divided by market ptf variance

pi <- c(lambda_m)*sig%*%t(w_m) + r_f    #Vector of prior mean expected returns

#Vector of posterior mean expected returns:
post_m <- pi + tau*sig%*%t(P)%*%solve(tau*P%*%sig%*%t(P) + Omega)%*%(Q - P%*%pi)

#a series of returns in the posterior distribution is obtained from the mean of expected returns
require(MASS)
set.seed(33)
n_tirages <- nrow(return)
return_post <- mvrnorm(n_tirages, post_m, Sigma = 252*sig, tol = 1e-06, empirical = FALSE)/252

#I plot the efficient frontier
ptf_BL <- EF(returns = return_post, nports = nports, shorts = shorts, wmax = wmax)
low_BL <- which.min(ptf_BL$vol)                   
high_BL <- which.max(ptf_BL$return)
effi_BL <- ptf_BL[low_BL:high_BL,]

#Graph efficient frontier
col<-c("darkblue","indianred")
par(mar=c(7,6,4,4),xpd=T)
plot(100*ptfs_no_s$vol,100*ptfs_no_s$return,col="darkblue", lwd=2,xlab="standard deviation (%)",
     ylab="expected return (%)",las=1, type="l",pch=20,ylim=100*range(c(ptf_BL$return,ptfs_no_s$return)),
     xlim=100*range(c(ptf_BL$vol,ptfs_no_s$vol)))
lines(100*ptf_BL$vol,100*ptf_BL$return,col="indianred", lwd=2)
legend("bottom", horiz=T,inset = c(0,-0.4),text.col=col,pch=rep(NA,2),lty=rep(1,2),col=col, bty="n",
       legend= c("Markowitz efficient frontier","BL efficient frontier"))

#I plot weightings for all portfolios on the frontier
cum_w_BL<-apply(ptf_BL[,grep("w",colnames(ptf_BL))],1,cumsum)
colvector<-rainbow(ncol(return_post))
at_BL_1=seq(1,ncol(cum_w_BL), length.out=7)
at_BL_2<-seq(0,1,0.25)

cex<-0.8
par(mar=c(8,4,4,4) + 0.1,xpd=T, cex.axis=cex)
for (i in 1:nrow(cum_w_BL)){
  plot(1:ncol(cum_w_BL),cum_w_BL[1+nrow(cum_w_BL)-i,],xlab="",ylab="", ylim=range(cum_w_BL),xlim=c(0,ncol(cum_w_BL)),
       las=1,col=colvector[i],pch=20,axes=F)
  polygon(c(1:ncol(cum_w_BL), ncol(cum_w_BL):1),c(rep(0,ncol(cum_w_BL)),rev(cum_w_BL[nrow(cum_w_BL)+1-i,])),
          col=colvector[i])
  par(new=T)}
axis(1, at=at_BL_1, labels=round(100*ptf_BL$ret,1)[at_BL_1], cex.axis = cex)
axis(2, at=at_BL_2,labels=100*at_BL_2,cex.axis =cex)
mapply(mtext, c("expected return (%)", "weights (%)"), side=c(1,2), line = rep(2.5,2))
legend("bottom",ncol=3,inset = c(0,-0.45),legend=rev(colnames(returns)),text.col=colvector,col=colvector,
       lty=1, bty="n")
box()