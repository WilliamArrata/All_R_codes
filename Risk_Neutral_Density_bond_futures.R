require("pacman")
pacman::p_load("stringr","Hmisc","stats","readxl","data.table","dplyr","tidyr","zoo", "janitor", "ggplot2")

##########################################   DOWNLOAD DATA    ##########################################

#1. Options prices
options <- read_excel("inputs/OATA_options_31_mai_2023.xlsx",1) %>% row_to_names(row_number = 1) %>% 
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
  labs(title = paste0(last(charac$option_matu)," OAT futures options prices at all strikes"), subtitle = "31th May 2023") +
  labs(y = "option premium (EUR)", x = "option strike (EUR)") + scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#3. Riskfree rates at options' maturities
rates <- read_excel("inputs/EUR_rates.xlsx") %>% mutate_if(is.character, as.numeric)

#get by extrapolation a risk free rate for each option maturity
rates_n <- approx(rates$term, rates$Yield, xout=charac$terms, method = "linear", n = 50, rule = 2, f = 0, 
                        ties = "ordered", na.rm = F)$y/100


###############################  CALIBRATION OF PARAMETERS  ##########################################

#European call & put prices, expected spot price as a function of transformed parameters a and b
#for a sum of 2 or 3 lognormals in B&S model

nb_log <- 3    #choose number of lognormal laws in the mixture: 2 or 3

call <- function(x, KC){                          #call price in the B&S model
  d1_C <- (x[1] + x[2]^2 - log(KC))/x[2]
  d2_C <- d1_C - x[2]
  call <- exp(-r*T)*(exp(x[1] + (x[2]^2/2))*pnorm(d1_C) - KC*pnorm(d2_C))}

esp <- function(x){ exp(x[1] + (x[2]^2/2))}          #expected value for a lognormal distribution

#call price when the underlying asset follows a mixture of lognormal laws
call_mix <- function(x, KC){
  ifelse(length(x) == 7, return(x[5]*call(x[c(1, 3)], KC) + (1 - x[5])*call(x[c(2, 4)], KC) ),
         return(x[7]*call(x[c(1, 4)], KC) + x[8]*call(x[c(2, 5)], KC) + (1 - sum(x[7:8]))*call(x[c(3, 6)], KC))) }

#expected spot price when the underlying asset follows a mixture of lognormal laws
esp_mix <- function(x){
  ifelse(length(x) == 7, x[5]*esp(x[c(1, 3)]) + (1 - x[5])*esp(x[c(2, 4)]),
         x[7]*esp(x[c(1, 4)]) + x[8]*esp(x[c(2, 5)]) + (1 - sum(x[7:8]))*esp(x[c(3, 6)]) )}

put_mix <- function(x, KP){ call_mix(x, KP) + exp(-r*T)*(KP - FWD)}   #put call parity

#The function to minimize over 7 or 10 parameters

MSE_mix <- function(x){
  C_INF <- pmax(esp_mix(x) - KC, call_mix(x, KC))
  C_SUP <- exp(r*T)*call_mix(x, KC)
  P_INF <- pmax(KP - esp_mix(x), put_mix(x, KP))
  P_SUP <- exp(r*T)*put_mix(x, KP)
  A <- as.numeric(KC <= esp_mix(x))
  B <- as.numeric(KP >= esp_mix(x))
  w_call <- A*first(tail(x, 2)) + (1 - A)*last(x)
  w_put <- B*first(tail(x, 2)) + (1 - B)*last(x)
  CALL <- w_call*C_INF + (1 - w_call)*C_SUP
  PUT <- w_put*P_INF + (1 - w_put)*P_SUP
  RES_C <- sum((C - CALL)^2, na.rm = T)
  RES_P <- sum((P - PUT)^2, na.rm = T)
  RES_F <- (FWD - esp_mix(x))^2
  MSE_mix <- RES_C + RES_P + RES_F
  return(MSE_mix)
}

#The probability density function

sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }
PDF <- function(x, y){
  ifelse(length(x) == 5, return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1 - x[5]), y) ),
         return(sub(x[c(1, 4, 7)], y) + sub(x[c(2, 5, 8)], y) + sub( c(x[c(3, 6)], 1-sum(x[7:8])), y))) }

#Calibration of parameters

#weights on itm and otm options fixed for the moment at 0.5 each thus 1st optim on some param only
PR <- seq(0.1, 0.49, 0.01)

if(nb_log != 2){
  PR <- seq(0.1, 1, 0.01)                           #range of weights on the first 2 densities
  PR <- expand.grid(c(rep(list(PR), 2)))
  PR <- PR[rowSums(PR) < 0.9, ] }                   #sum of weights on the first 2 densities capped at 90%

#objective function to be minimized
objective <- function(x){
  ifelse( length(PR) !=2, MSE_mix( c(x[1:4], PR[i], rep(0.5, 2))),
          MSE_mix( c(x[1:6], PR[i, 1], PR[i, 2], rep(0.5, 2)))) }

#the density function
sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }
PDF <- function(x, y){
  ifelse(unique(lengths(params)) == 5,
         return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1 - x[5]), y) ),
         return(sub(x[c(1, 4, 7)], y) + sub(x[c(2, 5, 8)], y) + sub( c(x[c(3, 6)], 1 - sum(x[7:8])), y))) }

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
  m1 <- m2 <- m3 <- s1 <- s2 <- s3 <- SCE <- NA
  PARA <- as.matrix(data.frame(m1, m2, s1, s2, pr = PR, w1 = 0.5, w2 = 0.5, SCE))
  if(nb_log != 2){PARA <- as.matrix(data.frame(m1, m2, m3, s1, s2, s3, pr1 = PR[, 1], pr2 = PR[, 2], 
                                               w1 = 0.5, w2 = 0.5, p1_p2 = rowSums(PR), SCE))}
  start <- rep(c(log(FWD), 0.2), each = nb_log)
  lower <- rep(c(-10, 1e-6), each = nb_log)
  upper <- rep(c(10, 0.9), each = nb_log)
  for (i in 1:length(PR)){
    sol <- nlminb(start = start, objective = objective, lower = lower, upper = upper, 
                  control = list(iter.max = 500))
    PARA[i, grep(paste(c("m", "s"), collapse = "|"), colnames(PARA))] <- sol$par
    PARA[i, ncol(PARA)] <- sol$objective
  }
  
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
  
  solu <- constrOptim(param, MSE_mix, NULL, ui = UI, ci = CI, mu = 1e-05, control = list(iter.max = 2000), 
                      method = "Nelder-Mead")
  CV[[m]] <- solu$convergence
  
  #conversion of (a,b) into (mu, sigma)
  params[[m]] <- c(log(FWD) + (solu$par[1:2] - log(FWD))/T, solu$par[3:4]/sqrt(T), solu$par[5])
  if(nb_log != 2){params[[m]] <- c(log(FWD) + (solu$par[1:3] - log(FWD))/T, solu$par[4:6]/sqrt(T),
                                   solu$par[7:8])}
  nb_opt[[m]] <- nrow(prices)                                      #number of options for matu m
  range_px[[m]] <- range(KC, na.rm = T)                            #the range of strike for matu m
  PX[[m]] <- Reduce(seq, 1e4*range_px[[m]])*1e-4                   #values of x to comput PDF and CDF
}

#check that integral of PDF*dPX is worth 1, and if necessary augment the range of PX
DNR <- mapply(PDF, params, PX)
integ <- mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR, PX)
integ_fail <- integ < 0.99
x_axis <- c(0.98, 1.02)

for (m in which(integ_fail)){
  while (sum(rollmean(PDF(params[[m]], PX[[m]]), 2)*diff(PX[[m]]), na.rm = T) < 0.99){
    range_px[[m]] <- x_axis*range_px[[m]]
    PX[[m]] <- Reduce(seq, 1e4*range_px[[m]])*1e-4 }          #augmented values of x to comput PDF and CDF
}

DNR <- mapply(PDF, params, PX)
round(mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR, PX))

###############################  GRAPH OF RISK NEUTRAL DENSITIES       ########################################

#removing maturities where fit did not work
bad_fit <- which(round(mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR, PX)) > 1)
if(length(bad_fit) > 0){
  PX <- PX[-bad_fit]
  DNR <- DNR[-bad_fit]
  params <- params[-bad_fit]
  charac <- charac[-bad_fit, ]}

#Graph of risk neutral densities for Euribor futures prices
co <- rainbow(nrow(charac))
xlim <- range(PX)
ylim <- range(DNR)
series <- mapply(cbind, PX, DNR)

nb_log[unique(lengths(params))!=5] <- 3

cex <- 0.8
par(mar = c(8,4,4,4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "", ylab = "density", main = paste("RNDs from a mixture of", nb_log, "lognormals"), 
     xlim = xlim, ylim = ylim, las = 1)
mapply(lines, series, col = co)
title(sub = "OAT future price (% of par)", adj = 1, line = 2)
legend("bottom", inset = c(-0.05, -0.4), legend = format(as.yearmon(charac$option_matu), "%b %y"),
       horiz = T, col = co, lty = 1, bty = "n")

#Graph of risk neutral densities of prices with ggplot2
series <- lapply(mapply(cbind, series, as.character(charac$option_matu)), data.frame)
series <- do.call(rbind, series) %>% rename_with(~c("price", "density", "maturity")) %>%
  mutate_at(c("price", "density"), as.numeric)

ggplot() + geom_line(data = series, aes(x= price, y = density,  colour =  maturity)) +
  labs(title = "OAT futures prices (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "futures prices (% of par)") + 
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#######################  CALCULATION OF ACCRUED COUPON OF CTDs AT OPTION MATURITY ###########################

#Loading futures contracts characteristics and merging with options characteristics

bond_fut <- read_excel("inputs/OATA_fut_characteristics.xlsx", 1) %>%  
  rename_with(~c("fut_contract", "conv_factor", "ctd_cp", "ctd_matu")) %>% 
  mutate_at("fut_contract", ~word(., 1)) %>%  filter(fut_contract%in%charac$fut_contract) %>% 
  mutate_at("ctd_matu", as.Date, format = "%d/%m/%Y") %>% left_join(charac) %>% 
  mutate(prev_cp_dt = as.Date(paste0(format(option_matu, "%Y"),"-",format(ctd_matu, "%m-%d")))) %>% 
  mutate(prev_cp_dt = ifelse(as.numeric(option_matu - prev_cp_dt) < 0, 
                             as.Date(paste0(as.numeric(format(option_matu, "%Y")) - 1, "-", format(ctd_matu, "%m-%d"))),
                             prev_cp_dt)) %>% mutate_at("prev_cp_dt", as.Date) %>% 
  mutate(cc = as.numeric((option_matu - prev_cp_dt )/365), acc = ctd_cp*cc, 
         PX_liv = fut_price*conv_factor + acc, years = as.numeric(ctd_matu - option_matu)/365 )

####################  CONVERSION OF FUTURES PRICES INTO CTD PRICES THEN YIELDS AT MATU #######################

#conversion des prix futures en prix de CtD à matu de l'option
P <- mapply(function(x, y, z) 100*x*y + z, PX, bond_fut$conv_factor, bond_fut$acc)

N <- 100 + bond_fut$ctd_cp                                      #le flux payé à maturité par chaque CtD

years_c <- trunc(bond_fut$years)                                #le nb d'années pleines de paiement cp/ppal
full_y_c <- sapply(years_c, seq, from = 0)                      #toutes les années pleines intermédiaires
if(length(unique(years_c))==1){full_y_c <- as.list(as.data.frame(full_y_c))}
term <- mapply("+", full_y_c, bond_fut$years - years_c )         #le terme de tous les flux par CtD
term <- apply(term, 2, list)

cf <- split(rep(bond_fut$ctd_cp, years_c), 
            rep(seq_along(years_c), years_c))                    #les coupons (sauf le final) par CtD

#le YTM par obligation à partir de son prix, pour tous les prix possibles de chaque distribution
require('tvm')

tri <- list()
for (i in 1:nrow(bond_fut)){
  tri[[i]] <- list()
  for (j in 1:length(P[[i]])){
    tri[[i]][[j]] <- xirr(cf = c(-P[[i]][[j]], cf[[i]], N[i]), tau = c(0, term[[i]][[1]]), comp_freq = 1, 
                          interval = c(0, 10))}
  tri[[i]]<-unlist(tri[[i]])}

#############################  STATISTICS OF THE DISTRIBUTION #############################

#mean of futures prices at each options's maturity
E <- mapply(function(x, y, z) sum(rollmean(x*y, 2)*diff(z)), PX, DNR, PX)

#mean of CtD ytm from mean of futures prices at each option's maturity
E_y <- mapply(function(x, y, z, t, u, v) 
  xirr(cf = c(-(100*x*y+z), t, u), tau = c(0, v[[1]]), comp_freq = 1, interval = c(0, 10)),
  E, bond_fut$conv_factor, bond_fut$acc, cf, N, term)

#consistency check : ytm implicit to future contracts prices
y_fut <- mapply(function(x, y, z, t, u, v) 
  xirr(cf = c(-(x*y+z), t, u), tau = c(0, v[[1]]), comp_freq = 1, interval = c(0, 10)),
  bond_fut$fut_price, bond_fut$conv_factor, bond_fut$acc, cf, N, term)

#graph of implied ytm densities in base R
co <- rainbow(nrow(charac_1))
xlim_r <- range(tri, na.rm = T)
series_rev <- mapply(cbind, tri, DNR)

par(mar=c(7,4,4,4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "OAT future yield (%)", ylab = "density", xlim = xlim_r, ylim = ylim, las = 1, 
     main = "RNDs from a mixture of 2 lognormals")
mapply(lines, series_rev, col = co)
legend("bottom", inset = c(-0.05,-0.2), legend = charac$option_matu, horiz = T, col=co, lty = 1, bty = "n")

#graph of implied ytm densities using ggplot
series_rev <- lapply(mapply(cbind, series_rev, as.character(charac$option_matu)), data.frame)
series_rev <- do.call(rbind, series_rev) %>% rename_with(~c("ytm", "density", "maturity")) %>%
  mutate_at(c("ytm", "density"), as.numeric)

ggplot() + geom_line(data = series_rev, aes(x= ytm, y = density,  colour =  maturity)) +
  labs(title = "OAT future ytm (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "Cheapest-to-Deliver ytm") + 
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#Cumulative Density Function
sub_2 <- function(x, y){ x[3]*plnorm(y, meanlog = x[1], sdlog = x[2]) }

CDF <- function(x, y){
  ifelse(unique(lengths(params))==5,
         return(sub_2(x[c(1, 3, 5)], y) + sub_2(c(x[c(2, 4)], 1-x[5]), y) ),
         return(sub_2(x[c(1, 4, 7)], y) + sub_2(x[c(2, 5, 8)], y) + sub_2( c(x[c(3, 6)], 1-sum(x[7:8])), y)) ) }

#Graph of cumulative density functions for rates
NCDF <- mapply(CDF, params, PX)
series_CDF <- mapply(cbind, PX, NCDF)
series_CDF_tri <- mapply(cbind, tri, NCDF)

par(mar = c(7, 6, 4, 4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "", ylab = "cumulative probability", xlim = xlim, ylim = 0:1, las = 1,
     main = paste("RNDs from a mixture of", nb_log, "lognormals"))
mapply(lines, series_CDF, col = co)
title(sub = "OAT future price (% of par)", adj = 1, line = 2)
legend("bottom", inset = c(-0.05,-0.35), legend = charac$option_matu, horiz = T,col = co, lty = 1, bty = "n")

par(mar = c(7, 6, 4, 4) + 0.1, xpd = T, cex.axis = cex)
plot(NA, pch = 20, xlab = "", ylab = "cumulative probability", xlim = xlim_r, ylim = 0:1, las = 1,
     main = paste("RNDs from a mixture of", nb_log, "lognormals"))
mapply(lines, series_CDF_tri, col = co)
title(sub = "OAT future price (% of par)", adj = 1, line = 2)
legend("bottom", inset = c(-0.05,-0.35), legend = charac$option_matu, horiz = T,col = co, lty = 1, bty = "n")


#Standard deviation, skewness and kurtosis for the distribution at each options' maturity
moments <- function(x){
  return(mapply(function(x, y, z, t) sum(rollmean( z*(t - y)^x , 2)*diff(t)), x, E, DNR, PX))}

SD <- sqrt(moments(2))
SK <- moments(3)/SD^3
KU <- moments(4)/SD^4

#all statistics at a glance
desc_stats <- bond_fut %>% bind_cols(t(sapply(range_px, function(x) x/x_axis))) %>%
  rename_at(14:15, ~c("min_strike", "max_strike")) %>% 
  mutate_at(vars(contains("strike")), ~ .*conv_factor + acc) %>% 
  mutate_at(vars(contains("strike")), ~ 100*xirr(cf = c(-., rep(ctd_cp, trunc(years)), 100 + ctd_cp), 
                                                 tau = c(0, (0:years) + years - trunc(years)),
                                                 comp_freq = 1, interval = c(0, 10))) %>% 
  select(c(fut_contract, option_matu, fut_rate, min_strike, max_strike)) %>% 
  bind_cols(100*E_y, SD, SK, KU, nb_opt = unlist(nb_opt)) %>%
  rename_at(6:9, ~c("mean (%)", "stddev (%)", "skewness", "kurtosis")) 

#a few quantiles
nb_q <- 10
thres <- seq(nb_q)/nb_q
quantiles <- list()
for (i in 1:nrow(charac)){
  quantiles[[i]] <- list()
  for (j in 1:length(thres)){
    quantiles[[i]][[j]] <- mean(tri[[i]][c(min(which(NCDF[[i]] > thres[j] - 1e-5)),
                                           max(which(NCDF[[i]] < thres[j] + 1e-5)))])
  }
  quantiles[[i]] <- unlist(quantiles[[i]])
}

rate_spot <- 0.02848

#graph of quantiles through time with shaded areas
quantiles_2 <- bind_cols(c(0, charac$terms), rbind(rate_spot, do.call(rbind, quantiles))) %>%
  rename_with(~c("term", paste0("q", nb_q*thres)))

ggplot(quantiles_2, aes(x = term)) +
  geom_ribbon(aes(ymin = q1, ymax = q9, fill = "80%-90%"), alpha=0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q8, fill = "70%-80%"), alpha=0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q7, fill = "60%-70%"), alpha=0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q6, fill = "50%-60%"), alpha=0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q5, fill = "40%-50%"), alpha=0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q4, fill = "30%-40%"), alpha=0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "20%-30%"), alpha=0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q2, fill = "10%-20%"), alpha=0.2) +
  labs(x = "term (years)", y = "10Y OAT rate (%)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "bottom", legend.title=element_blank(), plot.margin = margin(1.2,.5,1.2,.5, "cm"))

#graph of quantiles through time unshaded
quantiles_3 <- bind_cols(rep(c(0, charac$terms), c(unique(lengths(quantiles)), lengths(quantiles)) ),
                         c(rep(rate_spot, unique(lengths(quantiles))), unlist(quantiles)),
                         rep(rev(paste0("q", nb_q*thres)), 1+ length(quantiles))) %>% 
  rename_all(~c("term", "value", "quantile"))

ggplot(quantiles_3, aes(term, value, color = quantile)) +  geom_line() +
  labs(x = "term (years)", y = "10Y OAT rate (%)") + 
  theme(legend.position= "bottom", legend.title=element_blank(), plot.margin = margin(1.2,.5,1.2,.5, "cm")) +
  scale_y_continuous(labels = scales::percent) 


#graph of quantile of order q for a given maturity
q <- 90
d <- 2
cutoff <- mean(PX[[d]][c(min(which(NCDF[[d]] > q/100 - 1e-5)),
                         max(which(NCDF[[d]] < q/100 + 1e-5)))])
dnr_q <- data.frame(x = PX[[d]], y = DNR[[d]]) %>% mutate(area = x > cutoff)

ggplot(data = dnr_q, aes(x = x)) + geom_ribbon(aes(ymin = 0, ymax = y, fill = area)) +
  geom_line(aes(y = y)) +
  annotate(geom = 'text', x = cutoff, y = -1, label = paste0("q",q), hjust = 0.5) +
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n =6)) +
  labs(x = '10-Y OAT contract future price (%)', y = 'probability density', 
       title = paste0("10-Y OAT contract RND and quantile of order ", q, "%"),
       subtitle = paste0("Probability density ", charac$option_matu[d])) +
  theme_light() +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))