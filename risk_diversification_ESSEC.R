require("pacman")
pacman::p_load("data.table","dplyr","ggplot2")

#######################   VARYING THE CORRELATION COEFFICIENT BETWEEN TWO ASSETS    ###############################

m <- c(3,7)*1e-2                                                         #Expected returns on assets A and B
V <- c(1.44, 6.25)*1e-2                                                  #Variance on assets A and B
rho <- c(-1, 1, 0, 0.25)                                                 #Various values for correlation btw A and B
cov <- apply(cbind(V[1],replicate(2,sqrt(prod(V))*rho),V[2]),1,list)     #all covariances matrix

w <- seq(0, 1, length.out = 300)
w <- cbind(w, rev(w))                                                    #a range of weights for both assets

mean <- 100*w%*%m                                                        #portfolio expected returns
sd <- pair <- list()                               
for (i in seq_along(rho)){
  sd[[i]] <- 100*sqrt(diag(w%*%matrix(cov[[i]][[1]], nrow = 2)%*%t(w)))  #portfolio standard deviations
  pair[[i]] <- data.frame(stdev = sd[[i]], mean = mean)}                

#Graph of the 300 portfolios in the mean standard deviation space
xlim <- 1.1*range(sd)
ylim <- c(0.8,1.1)*range(mean)
colv <- c(rainbow(length(rho)-1),"black")
lty <- rep(1:2,c(length(rho)-1,1))

par(mar=c(7, 4, 4, 3),xpd=T)
plot(pair[[1]][c(1,length(mean)),], las=1, xlab="", ylab="", xlim = xlim, ylim = ylim, pch = 20)
mapply(lines, pair, col = colv, lty = lty)
mapply(mtext, c(expression(sigma(r[P]) ("%"), E(r[P]) ("%"))), side=c(1,2), line = rep(2.5,2))
text(x = 1.1*pair[[1]]$stdev[c(1,length(mean))], y = pair[[1]]$mean[c(1,length(mean))], label = c("B","A"), font =2)
legend("bottom", horiz=T, inset = c(0,-0.4), text.col = colv, pch = rep(NA,4), lty = c(rep(1,3),3), col = colv, 
       bty="n", legend = parse(text = sprintf("rho[A][B] == %s",rho)))


#######################################   COMBINING ASSET BY PAIRS   ##################################

mu <- list(c(0.03,0.05), c(0.05,0.065), c(0.065,0.075),             #expected returns on pairs of assets
           c(0.075,0.11), c(0.11,0.15), c(0.03,0.15))
sig <- list(c(0.05,0.08,0.03), c(0.08,0.11,-0.03), c(0.11,0.12,0.06), #covariances between pairs of assets
          c(0.12,0.16,0.07), c(0.16,0.21,0.12), c(0.05,0.21,-0.08))

#assets coordinates
assets <- 100*data.frame(stdev = unique(unlist(lapply(sig, function(x) sqrt(head(x,2))))), 
                         mean = unique(unlist(mu)),
                         correl = unlist(lapply(sig, function(x) last(x)/sqrt(prod(head(x,2))))))

w_a <- seq(0, 1, length.out = 300)
w_a <- cbind(w_a, rev(w_a))       #a range of weights for each asset in the portfolio, varying by steps of 1/300

mu_tot <- sd_tot <- pair_tot <- list() #Expected returns & standard deviations on combinations of 2 consecutive assets
for (i in seq_along(mu)){
  mu_tot[[i]] <- 100*w_a%*%mu[[i]]
  sd_tot[[i]] <- 100*sqrt(diag(w_a%*%matrix(sig[[i]][c(1,3,3,2)], nrow=2)%*%t(w_a)))
  pair_tot[[i]] <- data.frame( stdev = sd_tot[[i]], mean = mu_tot[[i]])}

#the resulting curves
curve <- do.call(rbind, head(pair_tot, -1)) %>% arrange(mean)
curve_2 <- tail(pair_tot, 1)[[1]]

#Graph of curves formed by consecutive pairs of assets
ggplot(curve, aes(stdev, mean)) +
  geom_point(color="blue", size = 0.5) + geom_point(data = assets, aes(stdev, mean), size=3) +
  geom_text(data = assets, aes(x=stdev, y=mean, label = LETTERS[seq_along(mu)]), hjust = -0.6, vjust = 1.2, fontface = 2) +
  xlim(c(0.9,1.1)*range(curve$stdev)) + ylim(c(0.9,1.1)*range(curve$mean)) +
  labs(x = "standard deviation (%)", y = "expected return (%)") + theme(plot.margin = margin(1.2,.5,1.2,.5, "cm"))

#Adding in the graph the curve between most distant assets
asset_pair <- paste0(LETTERS[seq_along(mu)], LETTERS[c(tail(seq_along(mu), -1), 1)])
lab <- paste0(expression(rho~ plain("_") ),"~",  asset_pair, expression(~ plain("=") ), "~", round(assets$correl/100,2))

ggplot() + geom_point(data = assets, aes(x = stdev, y = mean, fill = as.factor(correl)), size = 3) +
  annotate(geom = "text", x = assets$stdev,  y = assets$mean, label = LETTERS[seq_along(mu)], 
           fontface = 2, hjust = -0.6, vjust = 1.2) +
  geom_segment(data = curve, aes(x = stdev, xend = dplyr::lead(stdev), y = mean, yend = dplyr::lead(mean)),
               size = 1, color="indianred") +
  geom_segment(data = curve_2, aes(x = stdev, xend = dplyr::lead(stdev), y = mean, yend = dplyr::lead(mean)), size = 1) +
  scale_fill_manual(values = rep("black", 6), labels = parse(text = sprintf(lab))) +
  labs(x = "standard deviation (%)", y = "expected return (%)", fill=NULL) +
  xlim(c(0.9,1.1)*range(curve_2$stdev)) + ylim(c(0.9,1.1)*range(curve_2$mean)) +
  theme(legend.position = "bottom", plot.margin = margin(1.2,.5,1.2,.5, "cm"))