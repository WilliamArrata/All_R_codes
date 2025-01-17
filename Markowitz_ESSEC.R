
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   ################

require("pacman")
pacman::p_load("tseries","readxl","dplyr", "tidyr", "data.table", "ggplot2", "viridis", "stringr")

#####################   DATA DOWNLOAD AND COMPUTATION OF EXPECTED RETURNS AND COVARIANCES   ################

#I load the data
setwd("Z://1_Service/1.2_Agents_Service/William/Enseignement/Mon cours 2022/mean variance")
returns <- read_excel("stock_prices.xlsx") %>%  select_if(is.numeric) %>%  mutate_all(~ ( (.) - shift(.))/(.)) %>% 
  na.omit() %>% rename_with(~word(., 1)) %>% as.matrix
annu <- 252
mean <- annu*matrix(colMeans(returns))                             #annualized expected returns
sig <- annu*cov(returns)                                           #annualized covariances

###########################   OPTIMAL PORTFOLIOS FOR SOME TARGET EXPECTED RETURNS   #########################

#1. short selling allowed

#do not put annualized returns into the optimizer as will compute cov matrix based on them
#unless cov matrix is also provided. however, same weights whether annualized or not

#Global Minimum Variance Portfolio
gmvp <- portfolio.optim(returns, shorts = T)         #only risk minimized when no return specified
w_gmvp <- gmvp$pw
gmvp_assets <- data.frame(ptf = "GMVP", stock = colnames(returns), weight = w_gmvp)

print(c(w_gmvp, annu*gmvp$pm, sqrt(annu)*gmvp$ps))     #weights, annualized expected return and stddev

ggplot(data = gmvp_assets, aes(x = stock, y = weight)) + geom_bar(stat = "identity", aes(fill = stock)) + 
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) +
  scale_fill_viridis(discrete = T, option = "E") + scale_y_continuous(labels = scales::percent) + 
  theme(legend.position = "none",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#Portfolio with an 20% target expected return
target_1 <- portfolio.optim(returns, pm = 0.20/annu, shorts = T)
w_1_s <- target_1$pw
assets_1 <- data.frame(ptf = "20% return", stock = colnames(returns), weight = w_1_s)

ggplot(data = assets_1, aes(x = stock, y = weight)) + geom_bar(stat = "identity", aes(fill = stock)) + 
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) + 
  scale_fill_viridis(discrete = T, option = "E") + scale_y_continuous(labels = scales::percent) + 
  theme(legend.position = "none",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#Portfolio with a 10% target expected return
target_2 <- portfolio.optim(returns, pm = 0.10/annu, shorts = T)
w_2_s <- target_2$pw
assets_2 <- data.frame(ptf = "10% return", stock = colnames(returns), weight = w_2_s)

ggplot(data = assets_2, aes(x = stock, y = weight)) + geom_bar(stat = "identity", aes(fill = stock)) + 
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) + 
  scale_fill_viridis(discrete = T, option = "E") + scale_y_continuous(labels = scales::percent) + 
  theme(legend.position = "none",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#Comparaison between weights in the three portfolios
#NB: the higher the target expected return, the higher the weight on the asset class (if expected return positive)

comp_short <- bind_rows(gmvp_assets, assets_1, assets_2 )

ggplot(comp_short, aes(fill = ptf, y = weight, x = stock)) + geom_bar(position = "dodge", stat = "identity") + 
  scale_fill_viridis(discrete = T, option = "E") + scale_y_continuous(labels = scales::percent) + 
  theme(legend.position = "bottom",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#with facet_wrap
library(hrbrthemes)
ggplot(comp_short, aes(fill = ptf, y = weight, x = stock)) + geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) +
  facet_wrap(~ptf) +  scale_y_continuous(labels = scales::percent) + theme_ipsum() + scale_fill_viridis(discrete = T, option = "E") +
  theme(legend.position = "bottom",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#base plot
c_y <- 1.3
col <- c("darkred", "darkblue", "darkgrey")
cex <- 0.8
par(mar = c(8,4,4,4) + 0.1, cex.axis = cex, xpd = T)
xx <- barplot(rbind(w_gmvp, w_1_s, w_2_s), ylab = "%", col = col, las = 1, beside = T, 
              ylim = c_y*range(w_gmvp, w_1_s, w_2_s), names.arg = colnames(returns))
text(x = xx, y = c(rbind(w_gmvp, w_1_s, w_2_s)), paste(round(100*c(rbind(w_gmvp, w_1_s, w_2_s))), "%", sep=""), 
     pos = 3, font = 3)
legend("bottom", horiz = T, inset = c(0,-0.35), legend = c("GMVP","10% target return","20% target return"), 
       text.col = col, pch=c(15), col = col, bty="n")
box()

#2. short selling forbidden

#The Global Minimum Variance Portfolio
gmvp_no <- portfolio.optim(returns) #GMVP. By default, short selling is banned
w_gmvp_n <- gmvp_no$pw
gmvp_assets_n <- data.frame(ptf = "GMVP", stock = colnames(returns), weight = w_gmvp_n)

ggplot(data = gmvp_assets_n, aes(x = stock, y = weight)) + geom_bar(stat = "identity", aes(fill = stock)) + 
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) +
  scale_fill_viridis(discrete = T, option = "E") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#A portfolio with an expected return of 20%
target_1_no <- portfolio.optim(returns, pm = 0.20/annu)
w_1_n <- target_1_no$pw
assets_1_n <- data.frame(ptf = "10% return", stock = colnames(returns), weight = w_1_n)

ggplot(data = assets_1_n, aes(x = stock, y = weight)) + geom_bar(stat = "identity", aes(fill = stock)) + 
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) +
  scale_fill_viridis(discrete = T, option = "E") + scale_y_continuous(labels = scales::percent) + 
  theme(legend.position = "none",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#A portfolio with an expected return of 10%
target_2_no <- portfolio.optim(returns, pm = 0.10/annu)
w_2_n <- target_2_no$pw
assets_2_n <- data.frame(ptf = "20% return", stock = colnames(returns), weight = w_2_n)

ggplot(data = assets_2_n, aes(x = stock, y = weight)) + geom_bar(stat = "identity", aes(fill = stock)) + 
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) +
  scale_fill_viridis(discrete = T, option = "E") + scale_y_continuous(labels = scales::percent) + 
  theme(legend.position = "none",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#Comparison of weights between the three portfolios
comp_no_short <- bind_rows(gmvp_assets_n, assets_1_n, assets_2_n )

ggplot(comp_no_short, aes(fill = ptf, y = weight, x = stock)) + geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) +
  facet_wrap(~ptf) +  scale_y_continuous(labels = scales::percent) + theme_ipsum() + scale_fill_viridis(discrete = T, option = "E") +
  theme(legend.position = "bottom",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#3. comparison short selling- no short selling

#Comparison between GMVPs
comp_gmvp <- bind_rows(gmvp_assets_n %>% mutate_at("ptf", ~"no short selling"), gmvp_assets %>% mutate_at("ptf", ~"short selling"))

ggplot(comp_gmvp, aes(fill = ptf, y = weight, x = stock)) + geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) +
  facet_wrap(~ptf) +  scale_y_continuous(labels = scales::percent) + theme_ipsum() + scale_fill_viridis(discrete = T, option = "E") +
  theme(legend.position = "bottom",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#Comparison of portfolios with 20% target expected returns : negative yielding assets now have zero weight
comp_2 <- bind_rows(assets_2_n %>% mutate_at("ptf", ~"no short selling"), assets_2 %>% mutate_at("ptf", ~"short selling"))

ggplot(comp_2, aes(fill = ptf, y = weight, x = stock)) + geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = paste( round(weight * 100), "%"), vjust = ifelse(weight >= 0, -0.7, 1.2))) +
  facet_wrap(~ptf) +  scale_y_continuous(labels = scales::percent) + theme_ipsum() + 
  scale_fill_viridis(discrete = T, option = "E") + theme(legend.position = "bottom",  plot.margin = margin(.8,.5,.8,.5, "cm"))

#######################################   FINDING THE EFFICIENT FRONTIER   #####################################

EF = function (returns, nports, shorts, wmax){
  max_ret <- max(mean)
  #max_ret<-(1+as.numeric(shorts)*0.5)*max(mean)     #la cible de renta maximale
  target <- seq(-max_ret, max_ret, len= nports)       #on définit les cibles de renta via nports et maxret
  reslow <- rep(-as.numeric(shorts), length(mean))    #vecteur de poids minimum
  reshigh <- rep(wmax,length(mean))                   #vecteur de poids max
  output <- list()
  for (i in seq_along(target)){
    sol <- NULL
    try(sol <- portfolio.optim(returns, pm = target[i]/annu, reshigh = reshigh, reslow = reslow, shorts=shorts), 
        silent = T)
    if(!is.null(sol)){
      output[[i]] <- c(i, sqrt(annu)*sol$ps, annu*sol$pm, sol$pw)
      names(output[[i]]) <- c("i", "vol", "return", paste0("w", 1:length(mean)))}
  }
  output <- as.data.frame(do.call(rbind, output))
  rownames(output) <- output$i
  return(output)
}

nports <- 300   #nb of ptf, thus we have 300 target expected returns

#1. Efficient frontier when short selling is allowed
shorts <- T
wmax <- 1       #max weight on each asset class

#minimum variance frontier and efficient frontier
ptfs <- EF(returns = returns, nports = nports, shorts = shorts, wmax = wmax)  #ptfs expected returns, variances and weights
low <- which.min(ptfs$vol)                                                    #variance minimale
high <- which.max(ptfs$return)                                                #expected return max
effi <- ptfs[low:high,]                                                       #coordonnées de chaque pf sur l'EF

#Graph minimum variance frontier and efficient frontier
par(mar = c(7,5,4,3), xpd = T)
plot(ptfs$vol[ c(low, high) ], ptfs$return[ c(low, high) ], las=1, xlab = "standard deviation", ylab = "expected return",
     ylim = c_y*range(ptfs$return), xlim = c(0.9, 1.1)*range(ptfs$vol), col = "lightblue", pch = 19)
lines(effi$vol, effi$return, col = "lightblue", lwd = 2)
lines(ptfs$vol, ptfs$return, col = "indianred", pch = 20)
legend("bottom", horiz = T, inset = c(0, -0.35), legend = c("Minimum variance frontier","Efficient frontier"),
       text.col = c("indianred","lightblue"), col = c("indianred","lightblue"), bty = "n", lty = 1)

#Graph minimum variance frontier and efficient frontier with ggplot2
ggplot(ptfs, aes(vol, return)) +  geom_point(aes(color = "Minimum Variance Frontier")) +
  geom_line(data = effi, aes(color = "Efficient frontier"), size = 1) +
  geom_point(data = ptfs[c(low, high),], aes(vol, return), size = 2) +
  annotate("text", x = ptfs$vol[ c(low, high)], y = ptfs$return[ c(low, high)], hjust = -0.1, vjust = 0.2, label = 
             as.expression(sapply(split(ptfs[ c(low, high),], ptfs[ c(low, high),]$vol), function(x)
               bquote(mu == .(round(x[[3]]*100, 1)) ~ "% ; "~ sigma == .(round(x[[2]]*100, 1)) ~ "%")))) +
  scale_y_continuous(labels = scales::percent) +  scale_x_continuous(labels = scales::percent) +
  labs(x = "standard deviation", y = "expected return") + 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#another try with geom_segment
ggplot() +
  geom_segment(data = ptfs, aes(x = vol, xend = dplyr::lead(vol), y = return, yend = dplyr::lead(return)), 
               size = 2.5, color="lightblue") +
  geom_segment(data = effi, aes(x = vol, xend = dplyr::lead(vol), y = return, yend = dplyr::lead(return)), 
               size = 1, color="indianred") +
  geom_point(data = ptfs[c(low,high),], aes(vol, return), color="black", size = 3) +
  annotate("text", x = ptfs$vol[c(low,high)], y = ptfs$return[c(low,high)], hjust = -0.1, vjust = 0.2, label = 
             as.expression(sapply(split(ptfs[c(low,high),], ptfs[c(low,high),]$vol), function(x)
               bquote(mu == .(round(x[[3]]*100, 1)) ~ "% ; "~ sigma == .(round(x[[2]]*100, 1)) ~ "%")))) +
  labs(x="standard deviation", y="expected return") + theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm")) +
  scale_y_continuous(labels = scales::percent) +  scale_x_continuous(labels = scales::percent)


#2. Efficient frontier when short selling is forbidden

shorts <- F

ptfs_no_s <- EF(returns = returns, nports = nports, shorts = shorts, wmax = wmax)     #some returns not attainable with sign constrained optim
low_no_s <- which.min(ptfs_no_s$vol)
high_no_s <- which.max(ptfs_no_s$return)
effi_no_s <- ptfs_no_s[low_no_s:high_no_s,]

#Graph minimum variance frontier and efficient frontier
par(mar=c(7,5,4,3),xpd=T)
plot(ptfs_no_s$vol[c(low_no_s,high_no_s)],ptfs_no_s$ret[c(low_no_s,high_no_s)],las=1, xlab="standard deviation",
     ylab="expected return",ylim=1.3*range(ptfs_no_s$return), xlim=c(0.8,1.3)*range(ptfs_no_s$vol),col="black", pch=19)
lines(effi_no_s$vol,effi_no_s$return,col="darkblue", lwd=2.0)
lines(ptfs_no_s$vol,ptfs_no_s$return,col="indianred",pch=20)
legend("bottom",horiz=T,inset = c(0,-0.35),legend=c("Minimum variance frontier","Efficient frontier"),
       text.col=c("indianred","darkblue"),col=c("indianred","darkblue"),lty=1, bty="n")

#Graph minimum variance frontier and efficient frontier with ggplot2
ggplot(ptfs_no_s, aes(vol,return)) +  geom_point(aes(color = "Minimum Variance Frontier")) +
  geom_line(data = effi_no_s, aes(color = "Efficient frontier"), size = 1) +
  geom_point(data = ptfs_no_s[c(low_no_s,high_no_s),], aes(vol, return), size = 2) +
  annotate("text", x = ptfs_no_s$vol[c(low_no_s,high_no_s)], y = ptfs_no_s$return[c(low_no_s,high_no_s)], hjust = -0.1, vjust = 0.2, label = 
             as.expression(sapply(split(ptfs_no_s[c(low_no_s,high_no_s), ], ptfs_no_s$vol[c(low_no_s,high_no_s)]), 
          function(x) bquote(mu == .(round(x[[3]]*100, 1)) ~ "% ; "~ sigma == .(round(x[[2]]*100, 1)) ~ "%")))) +
  labs(x = "standard deviation", y = "expected return") + scale_y_continuous(labels = scales::percent) +  
  scale_x_continuous(labels = scales::percent) + 
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm")) 

#Weights for each target return
cum_w <- apply( ptfs_no_s[ ,grep("w", colnames(ptfs_no_s))], 1, cumsum)

#Transition map of weights for all target returns
colvector <- rainbow(6)
at = seq(1, ncol(cum_w), length.out = 7)

cex <- 0.8
par( mar = c(8, 4, 4, 4) + 0.1, cex.axis = cex, xpd = T)
for (i in 1:nrow(cum_w)){
  plot(1:ncol(cum_w),cum_w[1+nrow(cum_w)-i,], xlab="",ylab="", xlim=c(0,ncol(cum_w)),ylim=range(cum_w),
       las=1,col=colvector[i],pch=20,axes=F)
  polygon(c(1:ncol(cum_w), ncol(cum_w):1), c(rep(0,ncol(cum_w)),rev(cum_w[nrow(cum_w)+1-i,])),col=colvector[i])
  par(new=T)}
axis(1, at = at, labels = round(100*ptfs_no_s$ret,1)[at], cex.axis = cex)
axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), cex.axis = cex)
mapply(mtext, c("expected return (%)", "weights (%)"), side=c(1,2), line = rep(2.5,2))
legend("bottom", ncol = 3, inset = c(0,-0.5), legend = rev(colnames(returns)), text.col = colvector,
       col = colvector, pch=c(15), bty="n")
box()


#Comparaison of frontiers short selling allowed short selling forbidden
col_no <- c("lightblue","blue","indianred","red")
par(mar = c(8, 5, 4, 3), xpd = T)
plot( c( ptfs$vol[low], ptfs_no_s$vol[low_no_s]), c(ptfs$return[low], ptfs_no_s$return[low_no_s]), las = 1,
     xlab = "standard deviation", ylab = "expected return",
     ylim = 1.2*range( c(ptfs_no_s$return, ptfs$return)), xlim = c(0.8, 1.1)*range( c(ptfs$vol, ptfs_no_s$vol)),
     col = col_no[ c(2, 4)], pch = 19)
lines(effi$vol, effi$return, col = col_no[2], lwd = 2.0)
lines(ptfs$vol, ptfs$return, col = col_no[1], pch = 20)
lines(effi_no_s$vol, effi_no_s$return, col = col_no[4], lwd = 2.0)
lines(ptfs_no_s$vol, ptfs_no_s$return, col = col_no[3], pch = 20)
legend("bottom", ncol = 2, inset = c(0,-0.5), text.col = col_no, col = col_no, lty = 1, bty = "n",
       legend=c("MV frontier with short", "EF with short", "MV frontier w/o short", "EF w/o short"))


#for every level of expected return, the ptf with the short sales constraint incurs a higher risk
compar<-merge(x = ptfs, y = ptfs_no_s, by = "i", all=T)
compar<-compar[apply(!is.na(compar),1,all),]
compar$diff_vol<-compar$vol.y-compar$vol.x

#3. Efficient frontier when individual asset weights are capped at 25%

wmax<-0.25
ptfs_25<-EF(returns=returns,nports=nports,shorts=shorts,wmax=wmax)
low_25<-which.min(ptfs_25$vol)
high_25<-which.max(ptfs_25$return) 
effi_25<-ptfs_25[low_25:high_25,]

#Graph efficient frontiers short selling allowed, no short selling and weights capped at 25%
col_f<-c("darkgrey","darkblue","indianred")
par(mar=c(7,5,4,3),xpd=T)
plot(effi_25$vol,effi_25$return, las=1, pch=20,col=col_f[1],
     xlab="", ylab="expected return", ylim=c(0.8,1.1)*range(c(effi_25$return,effi_no_s$return,effi$return)),
     xlim=c(0.9,1.1)*range(c(effi_25$vol,effi_no_s$vol,effi$vol)))
lines(effi$vol,effi$return, col=col_f[3], pch=20)
lines(effi_no_s$vol,effi_no_s$return,col=col_f[2], lwd=2.0)
title(sub="standard deviation",adj=1,line=2)
legend("bottom",horiz=T,inset = c(0,-0.35), text.col=col_f,col=col_f,lty=1, bty="n",
       legend=c("Cap on weights at 25%","Short selling forbidden","short selling allowed"))
box()

#Weights for all target returns
cum_w_25<-apply(ptfs_25[,grep("w",colnames(ptfs_25))],1,cumsum)

#Transition map of weights for all target returns
at_25=seq(1,ncol(cum_w_25), length.out=7)

cex<-0.8
par(mar=c(8,4,4,4) + 0.1, cex.axis=cex, xpd=T)
for (i in 1:nrow(cum_w_25)){
  plot(1:ncol(cum_w_25),cum_w_25[1+nrow(cum_w_25)-i,], xlab="",ylab="", xlim=c(0,ncol(cum_w_25)),
       ylim=range(cum_w_25), las=1,col=colvector[i],pch=20,axes=F)
  polygon(c(1:ncol(cum_w_25), ncol(cum_w_25):1), c(rep(0,ncol(cum_w_25)),rev(cum_w_25[nrow(cum_w_25)+1-i,])),
          col=colvector[i])
  par(new=T)}
axis(1, at=at_25, labels=round(100*ptfs_25$ret,1)[at_25], cex.axis=cex)
axis(2, at=seq(0,1,0.25),labels=seq(0,100,25),cex.axis=cex)
mapply(mtext, c("expected return (%)", "weights (%)"), side=c(1,2), line = rep(2.5,2))
legend("bottom",ncol=3,inset = c(0,-0.5),legend=rev(colnames(returns)),text.col=colvector,col=colvector,
       pch=c(15), bty="n")
box()