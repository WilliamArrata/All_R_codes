require("pacman")
pacman::p_load("stringr", "Hmisc", "stats", "readxl", "data.table", "zoo", "dplyr", "tidyr", "janitor", "ggplot2")

##########################################   DOWNLOAD DATA    ##########################################

#1. Options prices
options <- read_excel("inputs/ERA_options_31_mai_2023.xlsx")  %>% row_to_names(row_number = 1) %>% 
  clean_names() %>% select(contains(c("strike", "last"))) %>% mutate_if(is.character, ~replace_na(.,"matu")) %>% 
  rename_with(~c(outer(c("call_", "put_"), c("strike", "price"), paste0)))

#2. Futures contracts prices and maturities
charac <- options %>% mutate(mat = row_number()) %>% filter(if_any(everything(), ~ grepl('matu',.))) %>% 
  mutate(option_matu = word(call_strike, 1, 3), fut_price = as.numeric( word(call_strike, -1))) %>% 
  mutate(terms = as.numeric(gsub('[^0-9.-]','', word(option_matu, 2)))/365, fut_contract = word(call_strike,-2)) %>%
  select(-colnames(options)) %>% mutate_at("option_matu", ~as.Date(gsub("\\).*","",word(.,-1)), format = "%m/%d/%y"))

#graph option prices for the most remote maturity
last_mat <- options %>% mutate_if(is.character, as.numeric) %>% slice((last(charac$mat)+1):nrow(options))

cex <- 0.8
col <- c("lightblue","indianred")
par(mar = c(6, 4, 4, 4) + 0.1, xpd = T, cex.axis = cex)
plot(last_mat$call_strike, last_mat$call_price, xlim = range(c(last_mat$call_strike, last_mat$put_strike)),
     ylim = range( c(last_mat$call_price, last_mat$put_price) ), col=col[1], type="l", pch=20, xlab=" ",
     main = paste(last(charac$option_matu),"Euribor 3 mth option prices at all strikes, 05/31/2023",sep=" "),
     ylab = "option premium (EUR)")
lines(last_mat$put_strike, last_mat$put_price, col=col[2])
title(xlab="strike price (EUR)",adj=1)
legend("bottom", horiz=T, bty="n",inset=c(-0.05,-0.35),legend=c("calls","puts"),lty=1,text.col=col,col=col)

#3. Riskfree rates at options' maturities (discount prices)
rates <- read_excel("inputs/EUR_rates.xlsx") %>% mutate_if(is.character, as.numeric)

#get by linear extrapolation a risk free rate at each option maturity
rates_n <- approxExtrap(rates$term, rates$Yield, xout=charac$terms, method = "linear", n = 50, rule = 2, f = 0, 
                        ties = "ordered", na.rm = FALSE)$y/100

###############################  CALIBRATION OF PARAMETERS  ##########################################

call <- function(x, KC){                          #call price in the B&S model
  d1_C <- (x[1] + x[2]^2 - log(KC))/x[2]
  d2_C <- d1_C - x[2]
  call <- exp(-r*T)*(exp(x[1] + (x[2]^2/2))*pnorm(d1_C) - KC*pnorm(d2_C))}

esp <- function(x){exp(x[1]+(x[2]^2/2))}          #expected value for a lognormal distribution

#European call & put prices, expected spot price as a function of a and b for a sum of 2 lognormals in B&S model

call_2_log <- function(x, KC){ x[5]*call(x[c(1,3)], KC) + (1-x[5])*call(x[c(2,4)], KC)}
put_2_log <- function(x, KP){ call_2_log(x,KP) + exp(-r*T)*(KP - FWD)}                  #put call parity
esp_2_log <- function(x){ x[5]*esp(x[c(1,3)]) + (1-x[5])*esp(x[c(2,4)])}

#Function to minimize over 7 parameters

MSE_2_log <- function(x){
  C_INF <- pmax(esp_2_log(x) - KC,call_2_log(x,KC))
  C_SUP <- exp(r*T)*call_2_log(x,KC)
  P_INF <- pmax(KP - esp_2_log(x), put_2_log(x,KP))
  P_SUP <- exp(r*T)*put_2_log(x,KP)
  A <- as.numeric(KC<=esp_2_log(x))
  B <- as.numeric(KP>=esp_2_log(x))
  w_call <- A*x[6] + (1-A)*x[7]
  w_put <- B*x[6] + (1-B)*x[7]
  CALL <- w_call*C_INF + (1 - w_call)*C_SUP
  PUT <- w_put*P_INF + (1 - w_put)*P_SUP
  RES_C <- sum((C - CALL)^2, na.rm=T)
  RES_P <- sum((P - PUT)^2, na.rm=T)
  RES_F <- (FWD - esp_2_log(x))^2
  MSE_2_log <- RES_C + RES_P + RES_F
  return(MSE_2_log)
}

#weights on itm and otm options fixed for the moment at 0.5 each thus 1st optim on first 5 parameters
PR <- seq(0.1, 0.49, 0.01)

objective <- function(x){ MSE_2_log(c(x[1:4], PR[i], rep(0.5, 2))) }

#Calibration of the 7 parameters using market data
mat <- c(charac$mat,nrow(options))                         #adding one last term to mat for the loop
params <- CV <- PX <- range_px <- nb_opt <- list()
x_axis <- c(0.9, 1.05)

for (m in 1:length(charac$terms)){
  
  #Elements of the option price function which are not random variables
  T <- charac$terms[m]                                             #maturity m
  r <- rates_n[m]                                                  #discount rate for maturity m
  prices <- options %>% select(-put_strike) %>% slice(mat[m]:mat[m+1]) %>% 
    mutate_if(is.character, as.numeric) %>% na.omit %>% mutate_all(funs(./100))
  C <- prices$call_price                                           #prices of calls for maturity m
  P <- prices$put_price                                            #prices of puts for maturity m
  KC <- KP <- prices$call_strike                                   #strikes of options for maturity m
  FWD <- charac$fut_price[m]/100                                   #future price for maturity m
  range_px[[m]] <- x_axis*range(KC, na.rm = T)                    #the augmented range of strike for matu m
  PX[[m]] <- Reduce(seq, 1e4*range_px[[m]])*1e-4                   #values of x to comput PDF and CDF
  nb_opt[[m]] <- nrow(prices)                                      #number of options for matu m
  
  #1st optimization over 6 parameters to get initialization values for second optim
  PARA <- matrix(nrow = length(PR), ncol = 8, dimnames =
                   list(c(), c(paste0("m", seq(2)), paste0("s", seq(2)), "p", paste0("w", seq(2)), "SCE")))
  start <- rep(c(log(FWD),0.2), each = 2)
  lower <- rep(c(-10,1e-6), each = 2)
  upper <- rep(c(10,0.9), each = 2)
  
  for (i in 1:length(PR)){
    sol <- nlminb(start = start, objective = objective, lower = lower, upper = upper, 
                  control = list(iter.max=500))
    PARA[i, 1:4] <- sol$par
    PARA[i, 8] <- sol$objective
  }
  
  PARA[, 5] <- PR
  PARA[, 6:7] <- 0.5
  
  param <- PARA[which.min(PARA[,8]), -8]
  param[param==0] <- 1e-6
  
  #2nd optimization over 8 parameters
  L <- U <- rep(0, length(param))
  L[sign(param) == -1] <- 2*param[sign(param) == -1]
  L[sign(param) == 1] <- 1e-2*param[sign(param) == 1]
  U[sign(param) == -1] <- 1e-2*param[sign(param) == -1]
  U[sign(param) == 1] <- 2*param[sign(param) == 1]
  CI <- c(L, -U)
  UI <- rbind(diag(length(L)), -diag(length(L)))
  
  solu <- constrOptim(param, MSE_2_log, NULL, ui = UI, ci = CI, mu = 1e-05, control = list(iter.max = 2000), 
                      method = "Nelder-Mead")
  CV[[m]] <- solu$convergence
  
  #conversion of (a,b) into (mu, sigma)
  params[[m]] <- c(log(FWD) + (solu$par[1:2] - log(FWD))/T, 
                   solu$par[3:4]/sqrt(T),
                   solu$par[5])
}

#European call & put prices, expected spot price as a function of a and b for a sum of 3 lognormals in B&S model

call_3_log <- function(x, KC){
  x[7]*call(x[c(1,4)], KC) + x[8]*call(x[c(2,5)], KC) + (1-sum(x[7:8]))*call(x[c(3,6)], KC)}
put_3_log <- function(x,KP){ call_3_log(x,KP) + exp(-r*T)*(KP-FWD)}
esp_3_log <- function(x){ x[7]*esp(x[c(1,4)]) + x[8]*esp(x[c(2,5)]) + (1-sum(x[7:8]))*esp(x[c(3,6)])}

#function to minimize over 10 parameters

MSE_3_log <- function(x){
  C_INF <- pmax(esp_3_log(x) - KC,call_3_log(x,KC))
  C_SUP <- exp(r*T)*call_3_log(x,KC)
  P_INF <- pmax(KP - esp_3_log(x),put_3_log(x,KP))
  P_SUP <- exp(r*T)*put_3_log(x,KP)
  A <- as.numeric(KC<=esp_3_log(x))
  B <- as.numeric(KP>=esp_3_log(x))
  w_call <- A*x[9] + (1-A)*x[10]
  w_put <- B*x[9] + (1-B)*x[10]
  CALL <- w_call*C_INF + (1-w_call)*C_SUP
  PUT <- w_put*P_INF + (1-w_put)*P_SUP
  RES_C <- sum((C-CALL)^2, na.rm=T)
  RES_P <- sum((P-PUT)^2, na.rm=T)
  RES_F <- (FWD-esp_3_log(x))^2
  MSE_3_log <- RES_C + RES_P + RES_F
  return(MSE_3_log)
}

#weights on itm and otm options fixed for the moment at 0.5 each thus 1st optim on 8 parameters
PR <- seq(0.1,1,0.01)                  #range of weights on the first two densities
PR <- expand.grid(c(rep(list(PR), 2)))
PR <- PR[rowSums(PR)<0.9,]             #sum of the weights on the first two densities capped at 90%

objective<-function(x){ MSE_3_log(c(x[1:6],PR[i,1],PR[i,2],rep(0.5,2)))}

mat<-c(charac$mat,nrow(options))

params <- CV <- PX <- range_px <- nb_opt <- list()

for (m in 1:length(charac$terms)){
  
  #Elements of the option price function which are not random variables
  T <- charac$terms[m]                                             #maturity m
  r <- rates_n[m]                                                  #discount rate for maturity m
  prices <- options %>% select(-put_strike) %>% slice(mat[m]:mat[m+1]) %>% 
    mutate_if(is.character, as.numeric) %>% na.omit %>% mutate_all(funs(./100))
  C <- prices$call_price                                           #prices of calls for maturity m
  P <- prices$put_price                                            #prices of puts for maturity m
  KC <- KP <- prices$call_strike                                   #strikes of options for maturity m
  FWD <- charac$fut_price[m]/100                                   #future price for maturity m
  range_px[[m]] <- c(0.9, 1.05)*range(KC, na.rm = T)               #the range of strike for matu m
  PX[[m]] <- Reduce(seq, 1e4*range_px[[m]])*1e-4                   #values of x to comput PDF and CDF
  nb_opt[[m]] <- nrow(prices)                                      #number of options for matu m
  
  #Thus 1st optimization over first 8 parameters to get initialization values for second optim
  PARA<-matrix(nrow=nrow(PR),ncol=12,dimnames=
                 list(c(),c(paste0("m",seq(3)),paste0("s",seq(3)),paste0("p",seq(2)),paste0("w",seq(2)),"p1+p2","SCE")))
  lower<-rep(c(-10,1e-6),each=3)
  upper<-rep(c(10,0.8),each=3)
  start<-rep(c(log(FWD),0.2),each=3)
  
  for (i in 1:nrow(PR)){
    sol<-nlminb(start=start,objective=objective,lower=lower, upper = upper, control=list(iter.max=500))
    PARA[i,1:6]<-sol$par
    PARA[i,12]<-sol$objective
  }
  PARA[,7]<-PR[,1]
  PARA[,8]<-PR[,2]
  PARA[,9:10]<-0.5
  PARA[,11]<-rowSums(PR)
  param<-PARA[which.min(PARA[,12]),-12]
  param[param==0]<-1e-6
  
  #2nd optimization over 10 parameters
  L<-U<-rep(0,length(param))
  L[sign(param)==-1]<-2*param[sign(param)==-1]
  L[sign(param)==1]<-1e-2*param[sign(param)==1]
  U[sign(param)==-1]<-1e-2*param[sign(param)==-1]
  U[sign(param)==1]<-2*param[sign(param)==1]
  CI<-c(L,-U)
  UI<-rbind(diag(length(L)),-diag(length(L)))
  
  solu <- constrOptim(param, MSE_3_log, NULL, ui = UI, ci = CI, mu = 1e-05, control = list(iter.max = 2000),
                      method = "Nelder-Mead")
  
  CV[[m]] <- solu$convergence
  
  #conversion of (a,b) into (mu, sigma)
  params[[m]]<-c(log(FWD) + (solu$par[1:3] - log(FWD))/T,
                 solu$par[4:6]/sqrt(T),
                 solu$par[7:8])
}

###############################  GRAPH OF RISK NEUTRAL DENSITIES       ########################################

sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }

PDF <- function(x, y){
  ifelse(unique(lengths(params))==5,
         return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1-x[5]), y) ),
         return(sub(x[c(1, 4, 7)], y) + sub(x[c(2, 5, 8)], y) + sub( c(x[c(3, 6)], 1-sum(x[7:8])), y))) }

DNR <- mapply(PDF, params, PX)
mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR, PX)   #check that integral of PDF*dPX is worth 1

#Graph of risk neutral densities for Euribor futures prices
co <- rainbow(nrow(charac))
xlim <- range(PX)
ylim <- range(DNR)
series <- mapply(cbind, PX, DNR)

nb_log <- 2
nb_log[unique(lengths(params))!=5] <- 3

cex <- 0.8
par(mar = c(8, 4, 4, 4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "", ylab = "density", xlim = xlim, ylim = ylim, las = 1,
     main = paste("RNDs from a mixture of",nb_log,"lognormals"))
mapply(lines, series, col = co)
title(sub = "3 mth Euribor future price (EUR)", adj = 1, line = 2)
legend("bottom", inset = c(-0.05,-0.4), legend = charac$option_matu, ncol = 6, col = co, lty = 1, bty = "n")

#Graph of risk neutral densities of prices with ggplot2
series <- lapply(mapply(cbind, series, charac$option_matu), data.frame)
series <- do.call(rbind, series) %>% rename_with(~c("price", "density", "maturity")) %>%
  mutate_at("maturity", ~as.Date(.))

ggplot() + geom_line(data = series, aes(x = price, y = density, color = maturity)) +
  labs(title = "OAT futures prices (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "futures prices (% of par)") + 
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#Graph in base R of risk neutral densities for Euribor rates
xlim_r <- 1 - rev(xlim)                          #rates derived from prices
yields <- sapply(PX, function(x) 1 - rev(x))     #we use rev to display rates in increasing order
DNR_rev <- sapply(DNR, rev)                      #to be consistent with the reordering of yields
series_y <- mapply(cbind, yields, DNR_rev)

par(mar = c(8,4,4,4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "", ylab = "density", xlim = xlim_r, ylim = ylim, las = 1,
     main = paste("RNDs from a mixture of",nb_log,"lognormals"))
mapply(lines, series_y, col = co)
title(sub = "3 mth Euribor future rate", adj = 1, line = 2)
legend("bottom", inset = c(-0.05,-0.45), legend = charac$option_matu, ncol = 6, col = co, lty = 1, bty = "n")


#Ggplot2 graph of RNDs with RND maturity on x axis
x0 <- sapply(DNR_rev, function(x) ceiling(max(x)))  #the max of probability density value per RND
x0 <- c(0, cumsum(x0))                              #new xaxis : cumulative max probability densities
y0 <- c(0, charac$terms)                            #RND's terms
scale <- exp(diff(log(x0[-1])))/exp(diff(log(y0[-1])))   #ratio of consecutive growth rates

#a transformation of x0 which makes them proportional to options' terms
z0 <- y0*1.3*max(scale)*x0[which.min(scale)]/y0[which.min(scale)]

print(exp(diff(log(z0[-1])))/exp(diff(log(y0[-1]))))      #check that DNR max values now proportional to terms
print(cumsum(diff(z0))/cumsum(diff(y0)))

#The value of each RND following the first are shifted by a constant to allow for a representation proportional to terms
path <- mapply(function(x, y, z) cbind(density = x + y, yield = z),  DNR_rev,  cumsum(diff(z0)), yields)
path <- mapply(rbind, path, 0)
path <- mapply(cbind, path, charac$terms)
path <- lapply(path, data.frame)
path <- lapply(path, setNames, nm =c("density", "yield", "maturity"))
path <- do.call(rbind, path)

yield_min <- max(sapply(yields, function(x) min(x)))
yield_max <- min(sapply(yields, function(x) max(x)))

ggplot() +
  geom_path(data = path, aes(x = density, y = yield, colour = maturity)) +
  labs(x = "options' maturity (years)", y = 'Euribor 3 month values', title = "3-month Euribor RNDs") +
  scale_x_continuous(labels = function(x) round(x/(max(z0)/max(y0)), 2) , 
                     breaks = scales::pretty_breaks(n = 6), limits = c(0, 1.1*round(max(z0)))) +
  scale_y_continuous(labels = scales::percent, limits = c(yield_min/10, 0.5*yield_max) )  +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))

yield_min <- max(sapply(yields, function(x) min(x)))
yield_max <- min(sapply(yields, function(x) max(x)))

ggplot() +
  geom_path(data = do.call(rbind, path_2), aes(x = density, y = yield, colour = maturity)) +
  labs(x = "options' maturity (years)", y = 'Euribor 3 month values', title = "3-month Euribor RNDs") +
  scale_x_continuous(labels = function(x) round(x/(max(z0)/max(y0)), 2) ,
                     breaks = scales::pretty_breaks(n = 6), limits = c(0, 1.1*round(max(z0)))) +
  scale_y_continuous(labels = scales::percent, limits = c(yield_min/10, 0.5*yield_max) )  +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))

#Cumulative Density Function for any maturity for a sum of 2 or 3 lognormals
sub_2 <- function(x, y){ x[3]*plnorm(y, meanlog = x[1], sdlog = x[2]) }

CDF <- function(x, y){
  ifelse(unique(lengths(params))==5,
         return(sub_2(x[c(1, 3, 5)], y) + sub_2(c(x[c(2, 4)], 1-x[5]), y) ),
         return(sub_2(x[c(1, 4, 7)], y) + sub_2(x[c(2, 5, 8)], y) + sub_2( c(x[c(3, 6)], 1-sum(x[7:8])), y)) ) }

#Graph of cumulative density functions for contract prices and rates
NCDF <- mapply(CDF, params, PX)
NCDF_rev <- sapply(NCDF, rev)
series_CDF <- mapply(cbind, PX, NCDF)
series_CDF_rev <- mapply(cbind, yields, NCDF_rev)

#prices
par(mar = c(8,6,4,4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "", ylab = "cumulative probability", las = 1, xlim = xlim, ylim = 0:1, 
     main = paste("RNDs from a mixture of",nb_log,"lognormals"))
mapply(lines, series_CDF, col = co)
title(sub = "3 mth Euribor futures prices (%)", adj = 1, line = 2)
legend("bottom", inset = c(-0.05,-0.4), legend = format(as.yearmon(charac$option_matu), "%b %y"),
       ncol = 5, col = co, lty = 1, bty = "n")

#rates
par(mar = c(8,6,4,4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "", ylab = "cumulative probability", las = 1, xlim = xlim_r, ylim = 0:1, 
     main = paste("RNDs from a mixture of",nb_log,"lognormals"))
mapply(lines, series_CDF_rev, col = co)
title(sub = "3 mth Euribor rate (%)", adj = 1, line = 2)
legend("bottom", inset = c(-0.05,-0.5), legend = format(as.yearmon(charac$option_matu), "%b %y"),
       ncol = 5, col = co, lty = 1, bty = "n")

#mean, standard deviation, skewness and kurtosis for the distribution at each options' maturity
E_y <- 1 - mapply(function(x, y) sum(rollmean(x*y, 2)*diff(x)), PX, DNR)

moments <- function(x){
  return(mapply(function(x, y, z, t) sum(rollmean( ( (1 - t - y)^x )*z, 2)*diff(t)), x, E_y, DNR, PX))}

SD_y <- sqrt(moments(2))
SK_y <- moments(3)/SD_y^3
KU_y <- moments(4)/SD_y^4

charac <- charac %>% select(-c(mat, fut_contract)) %>% mutate(fut_rate = 100 - fut_price) %>%
  bind_cols(t(100*(1-sapply(range_px, rev)/rev(x_axis))), nb_opt = unlist(nb_opt), 100*E_y, 100*SD_y, SK_y, KU_y) %>%
  rename_at(c(5,6,8:11), ~c("min_strike (%)", "max_strike (%)", "mean (%)", "stddev (%)", "skewness", "kurtosis"))

#a few quantiles
nb_q <- 100
thres <- c(1, 5, 25, 50, 75, 95, 99)/nb_q
quantiles <- list()
for (i in 1:length(params)){
  quantiles[[i]] <- list()
  for (j in 1:length(thres)){
    quantiles[[i]][[j]] <- 1 - mean(PX[[i]][c(min(which(NCDF[[i]] > thres[j] - 1e-5)),
                                            max(which(NCDF[[i]] < thres[j] + 1e-5)))]) }
  quantiles[[i]] <- unlist(quantiles[[i]])
}

mean_r <- data.frame(term = c(0, charac$terms), mean = c(eur_spot, E_y))

#graph of quantiles through time with shaded areas
quantiles_2 <- bind_cols(c(0, charac$terms), rbind(eur_spot, do.call(rbind, quantiles))) %>%
  rename_with(~c("term", paste0("q", nb_q*thres)))

ggplot(quantiles_2, aes(x = term)) +
  geom_ribbon(aes(ymin = q1, ymax = q99, fill = "min-max")) +
  geom_ribbon(aes(ymin = q1, ymax = q95, fill = "1st decile - 9th décile")) +
  geom_ribbon(aes(ymin = q1, ymax = q75, fill = "1st quartile - 3rd quartile")) +
  geom_ribbon(aes(ymin = q1, ymax = q25), fill = "mistyrose1") +
  geom_ribbon(aes(ymin = q1, ymax = q5), fill = "lightblue") +
  geom_line(aes(y = q50, color = "median"), size = 0.6) +
  geom_line(aes(y = mean_r$mean, color = "mean"), size = 0.6) +
  scale_x_continuous( breaks = scales::pretty_breaks(n=5)) + theme_light() +
  labs(x = "term (years)", y = "Euribor rate (%)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "bottom", legend.title=element_blank(),
        legend.box = "vertical", plot.margin = margin(.5, .5, 1.2, .5, "cm")) +
  scale_fill_manual(values = c("mistyrose1", "plum3", "lightblue")) +
  scale_color_manual(values = c("darkred", "darkgreen"))

#graph of quantiles through time unshaded #1
ggplot(left_join(quantiles_2, mean_r), aes(x = term)) +
  geom_line(aes(y = q1)) +
  geom_line(aes(y = q5)) +
  geom_line(aes(y = q25)) +
  geom_line(aes(y = q50, color = "median")) +
  geom_line(aes(y = mean, color = "mean")) +
  geom_line(aes(y = q75)) +
  geom_line(aes(y = q95)) +
  geom_line(aes(y = q99)) +
  labs(x = "term", y = "Euribor rate (%)", color = c("median" = "deepskyblue",  "mean" = "coral1")) +
  scale_color_manual(values = c("median" = "deepskyblue", "mean" = "coral1")) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "bottom", legend.title=element_blank(), plot.margin = margin(1.2,.5,1.2,.5, "cm"))

#graph of quantiles through time unshaded #2
quantiles_3 <- bind_cols(rep(c(0, charac$terms), c(unique(lengths(quantiles)), lengths(quantiles)) ),
                         c(rep(eur_spot, unique(lengths(quantiles))), unlist(quantiles)),
                         rep(rev(paste0("q", nb_q*thres)), 1 + length(quantiles))) %>% 
  rename_all(~c("term", "value", "quantile"))


ggplot(quantiles_3, aes(term, value, color = quantile)) +  geom_line() +
  geom_line(data = mean_r, aes(term, mean), color = "coral1")+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "bottom", legend.title=element_blank(), 
        plot.margin = margin(1.2,.5,1.2,.5, "cm")) +
  labs(x = "term", y = "Euribor rate (%)", color = c("q50" = "deepskyblue")) +
  scale_color_manual(values = c("q50" = "deepskyblue"))
  
#graph of quantile of order q for maturity d
q <- 90
d <- 6
cutoff <- mean(PX[[d]][c(min(which(NCDF[[d]] > q/100 - 1e-5)),
                         max(which(NCDF[[d]] < q/100 + 1e-5)))])
dnr_q <- data.frame(x = PX[[d]], y = DNR[[d]]) %>% mutate(area = x > cutoff)

ggplot(data = dnr_q, aes(x = x)) + geom_ribbon(aes(ymin = 0, ymax = y, fill = area)) +
  geom_line(aes(y = y)) +
  annotate(geom = 'text', x = cutoff, y = -1, label = paste0("q",q), hjust = 0.5) +
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n =6)) +
  labs(x = 'Euribor 3 mth future price (%)', y = 'probability density', 
       title = paste0("Euribor 3 month RND and quantile of order ", q, "%"),
       subtitle = paste0("Probability density ", charac$option_matu[d])) +
  theme_light() +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))
