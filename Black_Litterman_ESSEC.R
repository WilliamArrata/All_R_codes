
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