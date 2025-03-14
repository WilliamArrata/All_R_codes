require("pacman")
pacman::p_load("stringr", "Hmisc", "stats", "readxl", "data.table", "zoo", "dplyr", "tidyr",
               "janitor", "ggplot2", "lubridate")

##########################################   DOWNLOAD DATA    ##########################################

#1. Options prices
options <- read_excel("inputs/OATA_options_31_mai_2023.xlsx")  %>% row_to_names(row_number = 1) %>% 
  clean_names() %>% select(contains(c("strike", "last"))) %>% mutate_if(is.character, ~replace_na(.,"matu")) %>% 
  rename_with(~c(outer(c("call_", "put_"), c("strike", "price"), paste0))) %>% 
  mutate(option_matu = ifelse(put_price == "matu", mdy(gsub("\\).*", "", word(call_strike, 3))), NA )) %>% 
  fill("option_matu", .direction = "down") %>% mutate_at("option_matu", as.Date)

#2. Futures contracts prices and maturities
charac <- options %>% filter(if_any(everything(), ~ grepl('matu',.))) %>% 
  mutate(fut_price = as.numeric( word(call_strike, -1)), terms = as.numeric(gsub('[^0-9.-]','', word(call_strike, 2)))/365) %>% 
  mutate(fut_contract = word(call_strike, - 2)) %>% select(c("fut_price", "terms", "fut_contract", "option_matu")) 

options <- options %>% group_by(option_matu) %>% slice(-1) %>% ungroup() %>% mutate_at(-ncol(.), as.numeric) 

#option prices for all maturities on one plot
ggplot() +  geom_line(data = options, aes(x = call_strike, y = call_price, group = option_matu, color = option_matu)) +
  geom_line(data = options, aes(x = call_strike, y = put_price, group = option_matu, color = option_matu)) + 
  labs(x = "option strike (EUR)", y = "option premium (EUR)") +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm")) +
  guides(color = guide_legend(title = "maturity", title.position = "top", title.hjust = 0.5))

#options prices, one plot per maturity
ggplot() +  geom_line(data = options, aes(x = call_strike, y = call_price, group = option_matu, color = "calls")) + 
  geom_line(data = options, aes(x = call_strike, y = put_price, group = option_matu, color = "puts")) +
  geom_bar(stat = "identity") + labs(x = "option strike (EUR)", y = "option premium (EUR)") + facet_wrap(~option_matu) + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#3. Riskfree rates at options' maturities (discount prices)
rates <- read_excel("inputs/EUR_rates.xlsx") %>% mutate_if(is.character, as.numeric)

#get by linear extrapolation a risk free rate at each option maturity
rates_n <- approx(rates$term, rates$Yield, xout = charac$terms, method = "linear", n = 50, rule = 2, f = 0, 
                  ties = "ordered", na.rm = F)$y/100

#####################################     INTRODUCTION OF THE FUNCTIONS TO OPTIMIZE ########################################


nb_log <- 2    #choose number of lognormal laws in the mixture: 2 or 3

#European call & put prices, expected spot price as a function of transformed parameters a and b
#for a sum of 2 or 3 lognormals in B&S model

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
  C_INF <- pmax(esp_mix(x) - KC, call_mix(x, KC))  #upper and lower bounds for American option prices
  C_SUP <- exp(r*T)*call_mix(x, KC)
  P_INF <- pmax(KP - esp_mix(x), put_mix(x, KP))
  P_SUP <- exp(r*T)*put_mix(x, KP)
  A <- as.numeric(KC <= esp_mix(x))   #indicator function worth 1 if call option itm, 0 otherwise
  B <- as.numeric(KP >= esp_mix(x))   # indicator function worth 1 if put option itm, 0 otherwise
  w_call <- A*first(tail(x, 2)) + (1 - A)*last(x)
  w_put <- B*first(tail(x, 2)) + (1 - B)*last(x)
  CALL <- w_call*C_INF + (1 - w_call)*C_SUP      #American call price estimator
  PUT <- w_put*P_INF + (1 - w_put)*P_SUP         #American put pric estimator
  MSE_mix <- sum((C - CALL)^2, na.rm = T) + sum((P - PUT)^2, na.rm = T) + (FWD - esp_mix(x))^2
  return(MSE_mix)
}

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

#The probability density function
sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }
PDF <- function(x, y){
  ifelse(unique(lengths(params)) == 5,
         return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1 - x[5]), y) ),
         return(sub(x[c(1, 4, 7)], y) + sub(x[c(2, 5, 8)], y) + sub( c(x[c(3, 6)], 1 - sum(x[7:8])), y))) }

#The cumulative Density Function
sub_2 <- function(x, y){ x[3]*plnorm(y, meanlog = x[1], sdlog = x[2]) }
CDF <- function(x, y){
  ifelse(unique(lengths(params)) == 5,
         return(sub_2(x[c(1, 3, 5)], y) + sub_2(c(x[c(2, 4)], 1 - x[5]), y) ),
         return(sub_2(x[c(1, 4, 7)], y) + sub_2(x[c(2, 5, 8)], y) + sub_2( c(x[c(3, 6)], 1 - sum( x[7:8])), y)) ) }


###############################  CALIBRATION OF PARAMETERS  ##########################################

#Calibration of the 7 parameters using market data
params <- CV <- PX <- range_px <- nb_opt <- list()
x_axis <- c(0.98, 1.02)

for (m in 1:length(charac$option_matu)){
  
  #Elements of the option price function which are not random variables
  T <- charac$terms[m]                                                   #maturity m
  r <- rates_n[m]                                                        #discount rate for maturity m
  prices <- options %>% filter(option_matu == charac$option_matu[m])     #isolating data for maturity m
  C <- prices$call_price                                                 #prices of calls for maturity m
  P <- prices$put_price                                                  #prices of puts for maturity m
  KC <- KP <- prices$call_strike                                         #strikes of options for maturity m
  FWD <- charac$fut_price[m]                                             #future price for maturity m
  
  #1st optimization excluding weights on itm and otm options, to get initialization values for second optim
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
    PARA[i, grep( paste( c("m", "s"), collapse = "|"), colnames(PARA))] <- sol$par
    PARA[i, "SCE"] <- sol$objective
  }
  PARA <- PARA[ !is.na(PARA[, "m1"]), ]
  
  if(length(PARA) > 0 ){ 
    param <- PARA[which.min(PARA[, "SCE"]), -ncol(PARA)]
    param[param == 0] <- 1e-6
    
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
  }
  else(params[[m]] <- PX[[m]] <- 0)
  nb_opt[[m]] <- nrow(prices)                                       #number of options for matu m
  range_px[[m]] <- range(KC)                                        #the range of strike for matu m
  PX[[m]] <- Reduce(seq, 1e2*range_px[[m]])*1e-2                   #values of x to comput PDF and CDF
  
}


#check that integral of PDF*dPX is worth 1, and if necessary augment the range of PX
DNR <- mapply(PDF, params, PX)
integ <- mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR, PX)
integ_fail <- integ < 0.99
print(lapply(PX, range))
PX_2 <- PX
range_px_2 <- range_px
for (z in 1:length(params)){
  counter <- 0
  while (counter < 50 & sum(rollmean(PDF(params[[z]], PX_2[[z]]), 2)*diff(PX_2[[z]]), na.rm = T) < 0.99){
    counter <- 1 + counter
    range_px_2[[z]] <- x_axis*range_px_2[[z]]
    PX_2[[z]] <- Reduce(seq, 1e2*range_px_2[[z]])*1e-2
    print(counter)
  }   
}
rg_1 <- lapply(PX, range)
rg_2 <- lapply(PX_2, range)
extension <- mapply("/", rg_2, rg_1)
keep <- which(extension[2, ] < 5)
PX_2 <- PX_2[keep]
params_2 <- params[keep]
charac_2 <- charac[keep, ]
DNR_2 <- mapply(PDF, params_2, PX_2)
bad_fit <- round(mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR_2, PX_2))
bad_fit <- which(bad_fit > 1 | is.na(bad_fit)) 

if(length(bad_fit) > 0){
  PX_2 <- PX_2[-bad_fit]
  DNR_2 <- DNR_2[-bad_fit]
  params_2 <- params_2[-bad_fit]
  charac_2 <- charac_2[-bad_fit, ]}

NCDF <- mapply(CDF, params_2, PX_2)

#check that now all sum to 1
print(round(mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR_2, PX_2)))


#######################  CALCULATION OF ACCRUED COUPON OF CTDs AT OPTION MATURITY ###########################

#Loading futures contracts characteristics and merging with options characteristics

bond_fut <- read_excel("inputs/OATA_fut_characteristics.xlsx", 1) %>%  
  rename_with(~c("fut_contract", "conv_factor", "ctd_cp", "ctd_matu")) %>% 
  mutate_at("fut_contract", ~word(., 1)) %>%  filter(fut_contract%in%charac_2$fut_contract) %>% 
  mutate_at("ctd_matu", as.Date, format = "%d/%m/%Y") %>% left_join(charac_2) %>% 
  mutate(prev_cp_dt = as.Date(paste0(format(option_matu, "%Y"),"-",format(ctd_matu, "%m-%d")))) %>% 
  mutate(prev_cp_dt = ifelse(as.numeric(option_matu - prev_cp_dt) < 0, 
                             as.Date(paste0(as.numeric(format(option_matu, "%Y")) - 1, "-", format(ctd_matu, "%m-%d"))),
                             prev_cp_dt)) %>% mutate_at("prev_cp_dt", as.Date) %>% 
  mutate(cc = as.numeric((option_matu - prev_cp_dt )/365), acc = ctd_cp*cc, 
         PX_liv = fut_price*conv_factor + acc, years = as.numeric(ctd_matu - option_matu)/365 )

####################  CONVERSION OF FUTURES PRICES INTO CTD PRICES THEN YIELDS AT MATU #######################

#conversion des prix futures en prix de CtD à matu de l'option
P <- mapply(function(x, y, z) x*y + z, PX_2, bond_fut$conv_factor, bond_fut$acc)

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
    tri[[i]][[j]] <- xirr(cf = c(-P[[i]][[j]], cf[[i]], N[i]), tau = c(0, unlist(term[[i]]) ), comp_freq = 1, 
                          interval = c(0, 10))}
  tri[[i]] <- unlist(tri[[i]])}

#############################  STATISTICS OF THE DISTRIBUTION #############################

#mean of futures prices at each options's maturity
E <- mapply(function(x, y, z) sum(rollmean(x*y, 2)*diff(z)), PX_2, DNR_2, PX_2)

#mean of CtD ytm from mean of futures prices at each option's maturity
E_y <- mapply(function(x, y, z, t, u, v) 
  xirr(cf = c(-(x*y + z), t, u), tau = c(0, v[[1]]), comp_freq = 1, interval = c(0, 10)),
  E, bond_fut$conv_factor, bond_fut$acc, cf, N, term)

#consistency check : ytm implicit to future contracts prices
y_fut <- mapply(function(x, y, z, t, u, v) 
  xirr(cf = c(-(x*y + z), t, u), tau = c(0, v[[1]]), comp_freq = 1, interval = c(0, 10)),
  bond_fut$fut_price, bond_fut$conv_factor, bond_fut$acc, cf, N, term)

#Standard deviation, skewness and kurtosis for the distribution at each options' maturity
moments <- function(x){
  return(mapply(function(x, y, z, t) sum(rollmean( z*(t - y)^x , 2)*diff(t)), x, E, DNR, PX))}

SD <- sqrt(moments(2))
SK <- moments(3)/SD^3
KU <- moments(4)/SD^4

#all statistics at a glance
desc_stats <- bond_fut %>% bind_cols(t(sapply(range_px_2, function(x) x/x_axis))) %>%
  rename_at( (ncol(.) - 1): ncol(.), ~c("min_strike", "max_strike")) %>% 
  mutate_at(vars(contains("strike")), ~ .*conv_factor + acc) %>% 
  mutate_at(vars(contains("strike")), ~ 100*xirr(cf = c(-., rep(ctd_cp, trunc(years)), 100 + ctd_cp), 
                                                 tau = c(0, (0:years) + years - trunc(years)),
                                                 comp_freq = 1, interval = c(0, 10))) %>% 
  select(c(fut_contract, option_matu, fut_price, min_strike, max_strike)) %>% 
  bind_cols(mean = 100*E_y, stddev = SD, skewness = SK, kurtosis = KU, nb_opt = unlist(nb_opt)) 


######################  SEVERAL GRAPHS AROUND THE DISTRIBUTIONS:PDFs, CDFs, QUANTILES... ################

#Graph of risk neutral densities of prices
df_p <- mapply(cbind, price = PX_2, density = DNR_2, maturity = charac_2$option_matu)
df_p <- do.call(rbind, df_p) %>% data.frame %>% mutate_at("maturity", as.Date)

ggplot() + geom_line(data = df_p, aes(x = price, y = density, group = maturity, color = maturity)) +
  labs(title = "OAT future prices (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "futures prices (% of par)") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm")) +
  guides(color = guide_legend(title = "maturity", title.position = "top", title.hjust = 0.5))

#Graph of risk neutral densities of irr
df_y <- df_p %>% mutate(price = unlist(tri))

ggplot() + geom_line(data = df_y, aes(x = price, y = density, group = maturity, color = maturity)) +
  labs(title = "OAT future ytm (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "future rates") + 
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) + 
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm"))  +
  guides(color = guide_legend(title = "maturity", title.position = "top", title.hjust = 0.5))


#Graph of cumulative density functions for contract prices
ncdf_p <- mapply(cbind, price = PX_2, cdf = NCDF, maturity = charac_2$option_matu)
ncdf_p <- do.call(rbind, ncdf_p) %>% data.frame %>% mutate_at("maturity", as.Date)

ggplot() + geom_point(data = ncdf_p, aes(x = price, y = cdf, group = maturity, color = maturity), size = 0.5) +
  labs(title = "Euribor 3-month future contract (% of par)", subtitle = "Cumulative Density Function") +
  labs(y = "cumulative density", x = "futures contract price (% of par)") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm")) +
  guides(color = guide_legend(title = "maturity", title.position = "top", title.hjust = 0.5))


#a few quantiles
nb_q <- 1000
thres <- seq(nb_q)/nb_q
quantiles <- list()
for (i in 1:nrow(charac)){
  quantiles[[i]] <- list()
  for (j in 1:length(thres)){
    quantiles[[i]][[j]] <- mean(tri[[i]][c(min(which(NCDF[[i]] > thres[j] - 1e-5)),
                                           max(which(NCDF[[i]] < thres[j] + 1e-5)))])
  }
  quantiles[[i]] <- t(data.frame(c(charac_2$terms[i], unlist(quantiles[[i]])))) %>% data.frame %>% 
    rename_with(~c("term", paste0("q", nb_q*thres)))
}

quantiles_1 <- lapply(quantiles, pivot_longer, cols =! "term", names_to = "quantile", values_to = 'value')
quantiles_1 <- lapply(quantiles_1, fill, "value", .direction ="downup")
quantiles_1 <- do.call(rbind, quantiles_1)

#Graph of all quantiles
ggplot() +  geom_line(data = quantiles_1, aes(term, value, color = quantile)) + scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "none", plot.margin = margin(1.2,.5,1.2,.5, "cm")) +
  labs(x = "term", y = "Euribor rate (%)", color = c("q500" = "deepskyblue"))

#graph of specific quantiles with shaded areas
quantiles_2 <- do.call(rbind, quantiles)

ggplot(quantiles_2, aes(x = term)) +
  geom_ribbon(aes(ymin = q1, ymax = q999, fill = "min-max")) +
  geom_ribbon(aes(ymin = q1, ymax = q950, fill = "1st decile - 9th décile")) +
  geom_ribbon(aes(ymin = q1, ymax = q750, fill = "1st quartile - 3rd quartile")) +
  geom_ribbon(aes(ymin = q1, ymax = q250), fill = "mistyrose1") +
  geom_ribbon(aes(ymin = q1, ymax = q50), fill = "lightblue") +
  geom_line(aes(y = q50, color = "median"), size = 0.6) +
  scale_x_continuous( breaks = scales::pretty_breaks(n = 5)) + theme_light() +
  labs(x = "term (years)", y = "Euribor rate (%)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "bottom", legend.title=element_blank(), 
        legend.box = "vertical", plot.margin = margin(.5, .5, 1.2, .5, "cm")) +
  scale_fill_manual(values = c("mistyrose1", "plum3", "lightblue")) +
  scale_color_manual(values = c("darkred", "darkgreen"))


#graph of quantiles through time unshaded
ggplot(quantiles_2, aes(x = term)) +
  geom_line(aes(y = q1)) +
  geom_line(aes(y = q5)) +
  geom_line(aes(y = q25)) +
  geom_line(aes(y = q50, color = "median")) +
  geom_line(aes(y = q75)) +
  geom_line(aes(y = q95)) +
  geom_line(aes(y = q99)) +
  labs(x = "term", y = "Euribor rate (%)", color = c("median" = "deepskyblue",  "mean" = "coral1")) +
  scale_color_manual(values = c("median" = "deepskyblue", "mean" = "coral1")) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "bottom", legend.title=element_blank(), plot.margin = margin(1.2,.5,1.2,.5, "cm"))


ggplot(quantiles_2, aes(x = term)) +
  geom_ribbon(aes(ymin = q1, ymax = q900, fill = "80%-90%"), alpha = 0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q800, fill = "70%-80%"), alpha = 0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q700, fill = "60%-70%"), alpha = 0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q600, fill = "50%-60%"), alpha = 0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q500, fill = "40%-50%"), alpha = 0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q400, fill = "30%-40%"), alpha = 0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q300, fill = "20%-30%"), alpha = 0.2) +
  geom_ribbon(aes(ymin = q1, ymax = q200, fill = "10%-20%"), alpha = 0.2) +
  labs(x = "term (years)", y = "10Y OAT rate (%)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "bottom", legend.title=element_blank(), plot.margin = margin(1.2,.5,1.2,.5, "cm"))

#graph of quantile of order q for a given maturity
q <- 90
d <- 2
cutoff <- mean(PX_2[[d]][c(min(which(NCDF[[d]] > q/100 - 1e-5)),
                         max(which(NCDF[[d]] < q/100 + 1e-5)))])
dnr_q <- data.frame(x = PX_2[[d]], y = DNR_2[[d]]) %>% mutate(area = x > cutoff)

ggplot(data = dnr_q, aes(x = x)) + geom_ribbon(aes(ymin = 0, ymax = y, fill = area)) +
  geom_line(aes(y = y)) +
  annotate(geom = 'text', x = cutoff, y = -1, label = paste0("q",q), hjust = 0.5) +
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n =6)) +
  labs(x = '10-Y OAT contract future price (%)', y = 'probability density', 
       title = paste0("10-Y OAT contract RND and quantile of order ", q, "%"),
       subtitle = paste0("Probability density ", charac$option_matu[d])) +
  theme_light() +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))