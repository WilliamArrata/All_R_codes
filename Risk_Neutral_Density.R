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
last_mat <- options %>% mutate_if(is.character, as.numeric) %>% slice( (last(charac$mat) + 1) : nrow(options) ) %>% 
  ggplot() + geom_line(aes(x = call_strike, y = call_price, color = "calls" ) )+
  geom_line(aes(x = call_strike, y = put_price, color = "puts") ) +
  labs(title = paste0(last(charac$option_matu)," Euribor 3-month futures options prices at all strikes"), subtitle = "31th May 2023") +
  labs(y = "option premium (EUR)", x = "option strike (EUR)") + scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#3. Riskfree rates at options' maturities (discount prices)
rates <- read_excel("inputs/EUR_rates.xlsx") %>% mutate_if(is.character, as.numeric)

#get by linear extrapolation a risk free rate at each option maturity
rates_n <- approx(rates$term, rates$Yield, xout=charac$terms, method = "linear", n = 50, rule = 2, f = 0, 
                        ties = "ordered", na.rm = FALSE)$y/100

###############################  CALIBRATION OF PARAMETERS  ##########################################

#European call & put prices, expected spot price as a function of transformed parameters a and b
#for a sum of 2 or 3 lognormals in B&S model

call <- function(x, KC){                          #call price in the B&S model
  d1_C <- (x[1] + x[2]^2 - log(KC))/x[2]
  d2_C <- d1_C - x[2]
  call <- exp(-r*T)*(exp(x[1] + (x[2]^2/2))*pnorm(d1_C) - KC*pnorm(d2_C))}

esp <- function(x){ exp(x[1] + (x[2]^2/2))}          #expected value for a lognormal distribution

#call price when the underlying asset follows a mixture of lognormal laws
call_log <- function(x, KC){
ifelse(length(x) == 7, return(x[5]*call(x[c(1, 3)], KC) + (1 - x[5])*call(x[c(2, 4)], KC) ),
       return(x[7]*call(x[c(1, 4)], KC) + x[8]*call(x[c(2, 5)], KC) + (1 - sum(x[7:8]))*call(x[c(3, 6)], KC))) }

#expected spot price when the underlying asset follows a mixture of lognormal laws
esp_log <- function(x){
  ifelse(length(x) == 7, x[5]*esp(x[c(1, 3)]) + (1 - x[5])*esp(x[c(2, 4)]),
         x[7]*esp(x[c(1, 4)]) + x[8]*esp(x[c(2, 5)]) + (1 - sum(x[7:8]))*esp(x[c(3, 6)]) )}

put_log <- function(x, KP){ call_log(x, KP) + exp(-r*T)*(KP - FWD)}   #put call parity

#The function to minimize over 7 or 10 parameters

MSE_log <- function(x){
  C_INF <- pmax(esp_log(x) - KC, call_log(x, KC))
  C_SUP <- exp(r*T)*call_log(x, KC)
  P_INF <- pmax(KP - esp_log(x), put_log(x, KP))
  P_SUP <- exp(r*T)*put_log(x, KP)
  A <- as.numeric(KC <= esp_log(x))
  B <- as.numeric(KP >= esp_log(x))
  w_call <- A*first(tail(x, 2)) + (1 - A)*last(x)
  w_put <- B*first(tail(x, 2)) + (1 - B)*last(x)
  CALL <- w_call*C_INF + (1 - w_call)*C_SUP
  PUT <- w_put*P_INF + (1 - w_put)*P_SUP
  RES_C <- sum((C - CALL)^2, na.rm = T)
  RES_P <- sum((P - PUT)^2, na.rm = T)
  RES_F <- (FWD - esp_log(x))^2
  MSE_log <- RES_C + RES_P + RES_F
  return(MSE_log)
}

#The probability density function

sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }
PDF <- function(x, y){
  ifelse(length(x) == 5, return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1 - x[5]), y) ),
         return(sub(x[c(1, 4, 7)], y) + sub(x[c(2, 5, 8)], y) + sub( c(x[c(3, 6)], 1-sum(x[7:8])), y))) }

#Calibration of parameters

#weights on itm and otm options fixed for the moment at 0.5 each thus 1st optim on some param only
PR <- seq(0.1, 0.49, 0.01)

#objective function to be minimized
objective <- function(x){
  ifelse( length(PR) !=2, MSE_log( c(x[1:4], PR[i], rep(0.5, 2))),
          MSE_log( c(x[1:6], PR[i, 1], PR[i, 2], rep(0.5, 2)))) }

#the density function
sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }
PDF <- function(x, y){
  ifelse(unique(lengths(params))==5,
         return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1-x[5]), y) ),
         return(sub(x[c(1, 4, 7)], y) + sub(x[c(2, 5, 8)], y) + sub( c(x[c(3, 6)], 1-sum(x[7:8])), y))) }

#Calibration of the 7 parameters using market data
mat <- c( charac$mat, nrow(options))                         #adding one last term to mat for the loop
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

  #1st optimization over 6 parameters to get initialization values for second optim
  PARA <- matrix(nrow = length(PR), ncol = 8, dimnames =
                   list(c(), c(paste0("m", seq(2)), paste0("s", seq(2)), "pr", paste0("w", seq(2)), "SCE")))
  start <- rep(c(log(FWD), 0.2), each = 2)
  lower <- rep(c(-10, 1e-6), each = 2)
  upper <- rep(c(10, 0.9), each = 2)

   for (i in 1:length(PR)){
    sol <- nlminb(start = start, objective = objective, lower = lower, upper = upper, 
                  control = list(iter.max = 500))
    PARA[i, grep(paste(c("m", "s"), collapse = "|"), colnames(PARA))] <- sol$par
    PARA[i, ncol(PARA)] <- sol$objective
  }
  
  PARA[, grep("pr", colnames(PARA))] <- PR
  PARA[, grep("w", colnames(PARA))] <- 0.5

  param <- PARA[which.min(PARA[, ncol(PARA)]), -ncol(PARA)]
  param[param==0] <- 1e-6

  #2nd optimization over 8 parameters
  L <- U <- rep(0, length(param))
  L[sign(param) == -1] <- 2*param[sign(param) == -1]
  L[sign(param) == 1] <- 1e-2*param[sign(param) == 1]
  U[sign(param) == -1] <- 1e-2*param[sign(param) == -1]
  U[sign(param) == 1] <- 2*param[sign(param) == 1]
  CI <- c(L, -U)
  UI <- rbind(diag(length(L)), -diag(length(L)))
  
  solu <- constrOptim(param, MSE_log, NULL, ui = UI, ci = CI, mu = 1e-05, control = list(iter.max = 2000), 
                      method = "Nelder-Mead")
  CV[[m]] <- solu$convergence
  
  #conversion of (a,b) into (mu, sigma)
  params[[m]] <- c(log(FWD) + (solu$par[1:2] - log(FWD))/T, solu$par[3:4]/sqrt(T), solu$par[5])

  nb_opt[[m]] <- nrow(prices)                                      #number of options for matu m
  range_px[[m]] <- range(KC, na.rm = T)                            #the range of strike for matu m
  PX[[m]] <- Reduce(seq, 1e4*range_px[[m]])*1e-4                   #values of x to comput PDF and CDF
}

DNR <- mapply(PDF, params, PX)
integ <- mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR, PX)   #check that integral of PDF*dPX is worth 1
integ_fail <- integ < 0.99
x_axis <- c(0.98, 1.02)

for (m in which(integ_fail)){
  while (sum(rollmean(PDF(params[[m]], PX[[m]]), 2)*diff(PX[[m]]), na.rm = T) < 0.99){
    range_px[[m]] <- x_axis*range_px[[m]]
    PX[[m]] <- Reduce(seq, 1e4*range_px[[m]])*1e-4 }          #augmented values of x to comput PDF and CDF
}

DNR <- mapply(PDF, params, PX)
mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR, PX)

###############################  GRAPH OF RISK NEUTRAL DENSITIES       ########################################

#Graph of risk neutral densities for Euribor futures prices
co <- rainbow(nrow(charac))
xlim <- range(PX)
ylim <- range(DNR)
series <- mapply(cbind, PX, DNR)

nb_log <- 2    #nb of lognormals by default
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
df_p <- do.call(rbind, series) %>% rename_with(~c("price", "density", "maturity")) %>%
  mutate_at("maturity", as.Date)

ggplot() + geom_line(data = df_p, aes(x = price, y = density, group = maturity, color = maturity)) +
  labs(title = "Euribor 3-month prices (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "futures prices (% of par)") + 
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "none", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))


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

series_y <- lapply(mapply(cbind, yields, DNR_rev, charac$option_matu), data.frame)
df_y <- do.call(rbind, series_y) %>% rename_with(~c("rate", "density", "maturity")) %>%
  mutate_at("maturity", as.Date)

ggplot() + geom_point(data = df_y, aes(x = rate, y = density, group = maturity, color = maturity), size = 0.5) +
  labs(title = "Euribor 3-month rate (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "futures rate") + 
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "none", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

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
  scale_y_continuous(labels = scales::percent, limits = c(yield_min/10, yield_max) )  +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))

#Cumulative Density Function for any maturity for a sum of 2 or 3 lognormals
sub_2 <- function(x, y){ x[3]*plnorm(y, meanlog = x[1], sdlog = x[2]) }

CDF <- function(x, y){
  ifelse(unique(lengths(params))==5,
         return(sub_2(x[c(1, 3, 5)], y) + sub_2(c(x[c(2, 4)], 1-x[5]), y) ),
         return(sub_2(x[c(1, 4, 7)], y) + sub_2(x[c(2, 5, 8)], y) + sub_2( c(x[c(3, 6)], 1-sum(x[7:8])), y)) ) }

#Graph of cumulative density functions for contract pricers and rates
NCDF <- mapply(CDF, params, PX)
NCDF_rev <- sapply(NCDF, rev)
series_CDF <- mapply(cbind, PX, NCDF)
series_CDF_rev <- mapply(cbind, yields, NCDF_rev)


par(mar = c(8,6,4,4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "", ylab = "cumulative probability", las = 1, xlim = xlim, ylim = 0:1, 
     main = paste("RNDs from a mixture of",nb_log,"lognormals"))
mapply(lines, series_CDF, col = co)
title(sub = "3 mth Euribor rate (%)", adj = 1, line = 2)
legend("bottom", inset = c(-0.05,-0.4), legend = format(as.yearmon(charac$option_matu), "%b %y"),
       ncol = 5, col = co, lty = 1, bty = "n")


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

nb_q <- 1000
thres <- seq(nb_q)/nb_q
quantiles <- list()
for (i in 1:length(params)){
  quantiles[[i]] <- list()
  for (j in 1:length(thres)){
    quantiles[[i]][[j]] <- 1 - mean(PX[[i]][c(min(which(NCDF[[i]] > thres[j] - 1e-5)),
                                            max(which(NCDF[[i]] < thres[j] + 1e-5)))]) }
  quantiles[[i]] <- unlist(quantiles[[i]])
}

eur_spot <- 0.0347

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


#graph of quantiles through time unshaded
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

# quantiles_2 <- bind_cols(rep(c(0, charac$terms), each = unique(lengths(quantiles))),
#                          c(rep(eur_spot, unique(lengths(quantiles))), unlist(quantiles)),
#                          rep(paste0("q", nb_q*thres), 1 + length(quantiles))) %>% rename_with( ~c("term", "rate", "quantile"))


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
  
#graph of quantile of order q for d^th maturity
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