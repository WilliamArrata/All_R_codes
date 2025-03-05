require("pacman")
pacman::p_load("stringr", "Hmisc", "stats", "readxl", "data.table", "zoo", "dplyr", "tidyr", 
               "janitor", "ggplot2", "lubridate")

##########################################   DOWNLOAD DATA    ##########################################

#1. Options prices
options <- read_excel("inputs/ERA_options_31_mai_2023.xlsx")  %>% row_to_names(row_number = 1) %>% 
  clean_names() %>% select(contains(c("strike", "last"))) %>% mutate_if(is.character, ~replace_na(.,"matu")) %>% 
  rename_with(~c(outer(c("call_", "put_"), c("strike", "price"), paste0))) %>% 
  mutate(option_matu = ifelse(put_price == "matu", mdy(gsub("\\).*", "", word(call_strike, 3))), NA )) %>% 
  fill("option_matu", .direction = "down") %>% mutate_at("option_matu", as.Date)

#2. Futures contracts prices and maturities
charac <- options %>% filter(if_any(everything(), ~ grepl('matu',.))) %>% 
  mutate(fut_price = as.numeric( word(call_strike, -1)), terms = as.numeric(gsub('[^0-9.-]','', word(call_strike, 2)))/365) %>% 
  mutate(fut_contract = word(call_strike, - 2)) %>% select(c("fut_price", "terms", "fut_contract", "option_matu")) 

#graph option premia by maturity
ggplot() + geom_line(data = options, aes(x = call_strike, y = call_price, group = option_matu, color = "calls")) +  #group/fill
  geom_line(data = options, aes(x = call_strike, y = put_price, group = option_matu, color = "puts")) +
  geom_bar(stat = "identity") + labs(x = "option strike (EUR)", y = "option premium (EUR)") + facet_wrap(~option_matu) + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

#3. Riskfree rates at options' maturities (discount rates)
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

params <- CV <- PX <- range_px <- nb_opt <- list()
x_axis <- c(0.98, 1.02)

for (m in 1:length(charac$option_matu)){
    
    #Elements of the option price function which are not random variables
    T <- charac$terms[m]                                                   #maturity m
    r <- rates_n[m]                                                        #discount rate for maturity m
    prices <- options %>% filter(option_matu == charac$option_matu[m]) %>% 
      mutate_if(is.character, as.numeric) %>% na.omit 
    C <- prices$call_price                                                 #prices of calls for maturity m
    P <- prices$put_price                                                  #prices of puts for maturity m
    KC <- KP <- prices$call_strike                                         #strikes of options for maturity m
    FWD <- charac$fut_price[m]                                             #future price for maturity m
  
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
for (i in 1:length(options)){
  print(round(mapply(function(x,y) sum(rollmean(x, 2)*diff(y)), DNR_2, PX_2)))
}


######################  SEVERAL GRAPHS AROUND THE DISTRIBUTIONS:PDFs, CDFs, QUANTILES... ################

#Graph of risk neutral densities of prices
df_p <- mapply(cbind, price = PX_2, density = DNR_2, maturity = charac_2$option_matu)
df_p <- do.call(rbind, df_p) %>% data.frame %>% mutate_at("maturity", as.Date)

ggplot() + geom_line(data = df_p, aes(x = price, y = density, group = maturity, color = maturity)) +
  labs(title = "Euribor 3-month prices (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "futures prices (% of par)") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm")) +
  guides(color = guide_legend(title = "maturity", title.position = "top", title.hjust = 0.5))

#Graph of risk neutral densities of yields
df_y <- df_p %>% mutate_at("price", ~1 - ./100)

ggplot() + geom_line(data = df_y, aes(x = price, y = density, group = maturity, color = maturity)) +
  labs(title = "Euribor 3-month prices (%)", subtitle = "Probability density functions") +
  labs(y = "probability density", x = "future rates") + 
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) + 
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm"))  +
  guides(color = guide_legend(title = "maturity", title.position = "top", title.hjust = 0.5))

#Graph of RNDs with RND maturity on x axis
yields <- sapply(PX_2, function(x) 1 - x/100)
x0 <- sapply(DNR_2, function(x) ceiling(max(x)))  #the max of probability density value per RND
x0 <- c(0, cumsum(x0))                              #new xaxis : cumulative max probability densities
y0 <- c(0, charac_2$terms)                            #RND's terms
scale <- exp(diff(log(x0[-1])))/exp(diff(log(y0[-1])))   #ratio of consecutive growth rates

#a transformation of x0 which makes them proportional to options' terms
z0 <- y0*1.3*max(scale)*x0[which.min(scale)]/y0[which.min(scale)]

print(exp(diff(log(z0[-1])))/exp(diff(log(y0[-1]))))      #check that DNR max values now proportional to terms
print(cumsum(diff(z0))/cumsum(diff(y0)))

#The value of each RND following the first are shifted by a constant to allow for a representation proportional to terms
path <- mapply(function(x, y, z) cbind(density = x + y, yield = z),  DNR_2,  cumsum(diff(z0)), yields)
path <- mapply(rbind, path, 0)
path <- mapply(cbind, path, charac_2$terms)
path <- lapply(lapply(path, data.frame), setNames, nm =c("density", "yield", "maturity"))
path <- do.call(rbind, path)

yield_min <- max(sapply(yields, function(x) min(x)))
yield_max <- min(sapply(yields, function(x) max(x)))

ggplot() +  geom_path(data = path, aes(x = density, y = yield, colour = maturity)) +
  labs(x = "options' maturity (years)", y = 'Euribor 3 month values', title = "3-month Euribor RNDs") +
  scale_x_continuous(labels = function(x) round(x/(max(z0)/max(y0)), 2) , 
                     breaks = scales::pretty_breaks(n = 6), limits = c(0, 1.1*round(max(z0)))) +
  scale_y_continuous(labels = scales::percent, limits = c(yield_min/10, yield_max) )  +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))


#Graph of cumulative density functions for contract pricers and rates
ncdf_p <- mapply(cbind, price = PX_2, cdf = NCDF, maturity = charac_2$option_matu)
ncdf_p <- do.call(rbind, ncdf_p) %>% data.frame %>% mutate_at("maturity", as.Date)

ggplot() + geom_point(data = ncdf_p, aes(x = price, y = cdf, group = maturity, color = maturity), size = 0.5) +
  labs(title = "Euribor 3-month future contract (% of par)", subtitle = "Cumulative Density Function") +
  labs(y = "cumulative density", x = "futures contract price (% of par)") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm")) +
  guides(color = guide_legend(title = "maturity", title.position = "top", title.hjust = 0.5))

#mean, standard deviation, skewness and kurtosis for the distribution at each options' maturity
E_y <- 100 - mapply(function(x, y) sum(rollmean(x*y, 2)*diff(x)), PX_2, DNR_2)

moments <- function(x){
  return(mapply(function(x, y, z, t) sum(rollmean( ( (100 - t - y)^x )*z, 2)*diff(t)), x, E_y, DNR_2, PX_2))}

SD_y <- sqrt(moments(2))
SK_y <- moments(3)/SD_y^3
KU_y <- moments(4)/SD_y^4

charac_2 <- charac_2 %>% select(-c(fut_contract)) %>% mutate(fut_rate = 100 - fut_price) %>%
  bind_cols(t((100 - sapply(range_px, rev)/rev(x_axis))), nb_opt = unlist(nb_opt), E_y, SD_y, SK_y, KU_y) %>%
  rename_at(c(5,6,8:11), ~c("min_strike (%)", "max_strike (%)", "mean (%)", "stddev (%)", "skewness", "kurtosis"))

#quantiles of order 0.1%
nb_q <- 1000
thres <- seq(nb_q)/nb_q
quantiles <- list()
for (i in 1:length(params_2)){
  quantiles[[i]] <- list()
  for (j in 1:length(thres)){
    quantiles[[i]][[j]] <- 100 - mean(PX_2[[i]][c(min(which(NCDF[[i]] > thres[j] - 1e-5)),
                                            max(which(NCDF[[i]] < thres[j] + 1e-5)))]) }
  quantiles[[i]] <- t(data.frame(c(charac_2$terms[i], unlist(quantiles[[i]])))) %>% data.frame %>% 
    rename_with(~c("term", paste0("q", nb_q*thres)))
}

#Graph of all quantiles
quantiles_1 <- lapply(quantiles, pivot_longer, cols =! "term", names_to = "quantile", values_to = 'value')
quantiles_1 <- do.call(rbind, quantiles_1)

ggplot(quantiles_1, aes(term, value, color = quantile)) +  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position= "bottom", legend.title=element_blank(), 
        plot.margin = margin(1.2,.5,1.2,.5, "cm")) +
  labs(x = "term", y = "Euribor rate (%)", color = c("q500" = "deepskyblue")) +
  scale_color_manual(values = c("q500" = "deepskyblue"))

#Graph of selected quantiles
quantiles_2 <- do.call(rbind, quantiles)

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

#Graph of specific quantiles with shaded areas
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


#graph of quantile of order q for d^th maturity
q <- 90
d <- 6
cutoff <- mean(PX[[d]][c(min(which(NCDF[[d]] > q/100 - 1e-5)),
                         max(which(NCDF[[d]] < q/100 + 1e-5)))])
dnr_q <- data.frame(x = PX[[d]], y = DNR[[d]]) %>% mutate(area = x > cutoff)

ggplot(data = dnr_q, aes(x = x)) + geom_ribbon(aes(ymin = 0, ymax = y, fill = area)) +
  geom_line(aes(y = y)) +
  annotate(geom = 'text', x = cutoff, y = -1, label = paste0("q", q), hjust = 0.5) +
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 6)) +
  labs(x = 'Euribor 3 mth future price (%)', y = 'probability density', 
       title = paste0("Euribor 3 month RND and quantile of order ", q, "%"),
       subtitle = paste0("Probability density ", charac$option_matu[d])) +
  theme_light() +
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))