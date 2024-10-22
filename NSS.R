
################################   WILLIAM ARRATA - NELSON SIEGEL SVENSSON MODEL  ################################

require("pacman")
pacman::p_load("nloptr","readxl","dplyr","tidyr", "ggplot2")

#########################            MATRIX OF PARAMETERS           ############################

#The matrix of the 6 parameters
coeff <- data.frame(matrix(c(0.02, -0.02, "", 0.06, "", 13), nrow = 4, ncol = 6, byrow = T, dimnames =
                             list( c(), c(paste0("beta_", 1:4), paste0("lambda_", 1:2)))))
coeff_2 <- coeff_3 <- coeff

#varying the values of beta_3
coeff_2 <- coeff_2 %>% mutate(beta_3 = c(-8, -2, 2, 10)/100, lambda_1 = 1) %>% mutate_all(as.numeric)
coeff_2 <- split(coeff_2, coeff_2$beta_3)

#varying the values of lambda_1
coeff_3 <- coeff_3 %>% mutate(lambda_1 = c(0.2,0.5,1,2), beta_3 = 0.1) %>% mutate_all(as.numeric)
coeff_3 <- split(coeff_3, coeff_3$lambda_1)

####################################         NS MODEL        ######################################

#the function describing the hump of the curve
hump <- function(x, y) { x[,1]*((1-exp(-y/x[,2]))/(y/x[,2])-exp(-y/x[,2])) }

#the NS function, with one hump
NS <- function(x, y){
  curve <- x[,1] + x[,2]*(1-exp(-y/x[,4]))/(y/x[,4]) + hump(x[,3:4], y)
  return(curve)
}

col <- rainbow(4)
matu <- (1:240)/12

#The NS function with the 4 coefficients from coeff_2 with 4 possible values for beta_3
NS_1 <- lapply(coeff_2, function(x) 100*NS(x[c(1:3,5)], matu))

graph_ns <- data.frame(matu = rep(matu, length(NS_1)), taux = unlist(NS_1),
                       lambda = rep(names(coeff_2), each = lengths(NS_1)))

#Influence of beta_3 on the shape of the curve:
p <- ggplot(graph_ns, aes(matu, taux, color = lambda)) + geom_line(size = 1) +  
  labs(x = "term (years)", y = "ytm (%)") + 
  theme(legend.position = "bottom", plot.margin = margin(.8,.5,.8,.5, "cm")) +
  guides(color = guide_legend(title = expression(beta[3]), title.hjust = 1.5))

#The NS function with the 4 coefficients from coeff_3 with 4 possible values for lambda_1
NS_2 <- lapply(coeff_3, function(x) 100*NS(x[c(1:3,5)], matu))

#Influence of lambda_1 on the shape of the curve:
graph_ns_2 <- data.frame(matu = rep(matu, length(NS_2)), taux = unlist(NS_2),
                      lambda = rep(names(coeff_3), each = lengths(NS_2)))

lab <- as.expression(bquote(lambda == .(as.numeric(names(lapply(coeff_2, function(x) grep("lambda_1", x)))))))
lab <- parse(text = sprintf("lambda==%s", names(lapply(coeff_2, function(x) grep("lambda_1", x)))))

ggplot(graph_ns_2, aes(matu, taux, color = lambda)) + geom_line(size = 1) + labs(x = "term (years)", y = "ytm (%)") + 
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
  annotate(geom = 'text', x = 4.5, y = seq(2.25, 3.75, 0.5), label = lab, hjust = 0.025, 
           color = unique(sapply(ggplot_build(p)$data,'[[','colour')), vjust = -1, size = 4.5)

####################################          NSS MODEL         ###################################

#The NSS function, with two humps
NSS <- function(x, y){
  curve <- x[,1] + x[,2]*(1-exp(-y/x[,5]))/(y/x[,5]) + hump(x[,c(3,5)], y)+ hump(x[,c(4,6)], y)
  return(curve)
}

NSS_1 <- 100*mapply(NSS, coeff_2, list(matu))
NSS_1 <- split(t(NSS_1), 1: length(coeff_2))

plot(NA, type="l", xlim = range(matu), ylim = range(NSS_1), xlab = "term (years)", ylab = "ytm (%)")
lapply(NS_1, function(x, t) lines(matu, x), t=t)
text(x = matu[30], y = lapply(NSS_1, function(x) x[30]), pos = 1, col = col, cex = 1.5, 
     label = parse(text = sprintf("beta[3]==%s", names(lapply(coeff_2, function(x) grep("beta_3", x))))))

NSS_2 <- lapply(coeff_3, function(x) 100*NSS(x, matu))

plot(NA, type="l", xlim = range(matu), ylim = range(NSS_2), xlab="term (years)", ylab="ytm (%)")
lapply(NS_2, function(x, t) lines(matu, x), t=t)
text(x = matu[60], y = lapply(NSS_2, function(x) x[60]), pos = 1, col = col, cex = 1.5, 
     label = parse(text = sprintf("lambda[1]==%s", names(lapply(coeff_3, function(x) grep("lambda_1", x))))))


##########################   CORRELATION BETWEEN 2nd AND 3nd FACTOR LOADINGS   ###########################


loading_1 <- function(x, y){ (1-exp(-y/x))/(y/x) }
loading_2 <- function(x, y) { (loading_1(x ,y) - exp(-y/x)) }

lambda <- seq(0.1, 15, 0.01)       #we create a series of possible values for lambda_1

matu_1 <- head(matu, -60)             #a subset of maturities
matu_2 <- head(matu, -100)          #an even smaller subset of maturities
all_matu <- apply(mapply(paste0, apply(round(sapply(list(matu, matu_1, matu_2), range)), 2, as.list), "Y"),
              2, paste, collapse="-")

#values of 2nd and third factor loadings for all tested values of lambda_1 and corresp correl
slope_0 <- lapply(lambda, function(x) loading_1(x, matu))
hump_0 <- lapply(lambda, function(x) loading_2(x, matu))
correl_0 <- mapply(cor, slope_0, hump_0)
correl_0 <- data.frame(lambda = lambda, correlation = correl_0, matu = all_matu[1])

#values of 2nd and third factor loadings for all tested values of lambda_1 for a smaller set of maturities
#and corresponding correlations:
slope_1 <- lapply(lambda, function(x) loading_1(x, matu_1))
hump_1 <- lapply(lambda, function(x) loading_2(x, matu_1))
correl_1 <- mapply(cor, slope_1, hump_1)
correl_1 <- data.frame(lambda = lambda, correlation = correl_1, matu = all_matu[2])

#values of 2nd and third factor loadings for all tested values of lambda_1 for an even smaller set of maturities
#and corresponding correlationq
slope_2 <- lapply(lambda, function(x) loading_1(x, matu_2))
hump_2 <- lapply(lambda, function(x) loading_2(x, matu_2))
correl_2 <- mapply(cor, slope_2, hump_2)
correl_2 <- data.frame(lambda = lambda, correlation = correl_2, matu = all_matu[3])

#graph
correl_all <- rbind(correl, correl_1, correl_2)
ggplot() +
  geom_point(data = correl_all, aes(lambda, correlation, color = matu), size = 0.5) +
  labs(y = "correlation between slope and curvature factors", x = expression(lambda[1])) +
  theme(legend.position = "bottom", plot.margin = margin(.5,.5,.5,.5, "cm")) + 
  guides(color = guide_legend(title = "maturity spectrum", title.position = "top", title.hjust = 0.5))


########################   NON LINEAR ESTIMATION PROCEDURE WITH RANDOM INITIAL VALUES    ###########################


#I load the data
data <- read_excel("French_bonds_06_10_2023.xlsx",1) %>% filter(!Series%in%c("OATe","OATi")) %>%
  select(c("Maturity","Ask Yield to Maturity")) %>% rename_all(~c("term","ytm")) %>% 
  mutate_at("term", ~(as.Date(., format= "%d/%m/%Y") - as.Date("2023-10-06"))/365) %>% 
  mutate(term = as.numeric(term), ytm = ytm/100) %>% filter(term <= 10)


#A function to compute the sum of squares between theoretical and market rates for different combinations of param
GSS1 <- function(x, m, r){
  NSS1 <- SQ1 <- array(dim = c(nrow(data), 1, dim(h)[1]))
  SSQ1 <- matrix(nrow = dim(SQ1)[2], ncol = dim(SQ1)[3])
  
  for (i in 1:dim(NSS1)[1]){                  #theoretical rate by maturity
    NSS1[i,,] <- x[1,,] + x[2,,]*(1 - exp(-m[i]/x[5,,]))/(m[i]/x[5,,])+
      x[3,,]*((1 - exp(-m[i]/x[5,,]))/(m[i]/x[5,,]) - exp(-m[i]/x[5,,]))+
      x[4,,]*((1 - exp(-m[i]/x[6,,]))/(m[i]/x[6,,]) - exp(-m[i]/x[6,,]))
    for (j in 1:dim(SQ1)[2]){
      SQ1[i,j,]<-(r[i]-NSS1[i,j,])^2       #squared diff between theoretical and market rate by maturity
    }
    SSQ1 <- colSums(SQ1, dims = 1, na.rm=T)   #sum of squared diff across maturities
  }
  return(SSQ1)
}

#The six parameters are drawn randomly 10000 times from predefined intervals
n_samp <- 10000
h <- do.call(cbind, lapply(list((0:600)/10000, (-350:350)/10000, (-160:160)/1000, 
                                (-100:100)/1000, (330:630)/100, (30:240)/100),
                           sample, n_samp, replace = T))
h[h==0] <- 1e-6

#rearrange such that lower bound always inferior to upper bound
lower <- 0.9*h
upper <- 1.1*h
replace <- upper<lower
reorder_1 <- lower[replace]
reorder_2 <- upper[replace]
lower[replace] <- reorder_2
upper[replace] <- reorder_1

#Objective function to be optimized.
sub_3 <- function(x, m){ x[1]*((1-exp(-m/x[2]))/(m/x[2])-exp(-m/x[2])) }

GSS <- function(x, m, r){
  NSS <- x[1] + x[2]*(1-exp(-m/x[5]))/(m/x[5]) + sub_3(x[c(3,5)], m)+ sub_3(x[c(4,6)], m)
  SSQ <- sum((r-NSS)^2,na.rm=T)
  return(SSQ)
}

objective <- function(x){
  return(GSS(x, m = data$term, r = as.matrix(data$ytm)))
}

#Calibration of the 6 parameters
CI <- cbind(lower, -upper)
UI <- rbind(diag(6), -diag(6))

#changer l'expression de la valeur initiale
para <- CV <- error <- list()
for (i in 1:nrow(h)){
  sol <- constrOptim(h[i, ], objective, NULL, ui = UI, ci = CI[i, ], mu = 1e-05, 
                     control = list(iter.max = 2000), method = "Nelder-Mead")
  para[[i]] <- sol$par
  CV[[i]] <- sol$convergence
  error[[i]] <- objective(x =  para[[i]])
}

e1 <- which.min(error)

para_nls <- data.frame(matrix(unlist(para[e1]), nrow = 1, ncol = 6))  %>% rename_with(~colnames(coeff))
curve_nls <- data.frame(matu = data$term, rates = NSS(para_nls, data$term))

#graph
ggplot() +
  geom_point(data = curve_nls, aes(matu, rates, color = "fitted rates"), size = 1) +
  geom_point(data = data, aes(term, ytm, color = "market rates"), size = 1) +
  geom_line(data = curve_nls, aes(matu, rates, color = "fitted rates"), size = 0.5) +
  geom_line(data = data, aes(term, ytm, color = "market rates"), size = 0.5) +
  labs(y = "ytm (%)", x = "term (years)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))


#########################     NON LINEAR ESTIMATION PROCEDURE WITH GRID SEARCH    ########################


# Determination of the initial values to start non linear optimization: first grid search

#The six parameters are described with a grid
param <- mapply(seq, list(0, -3.5, -16, -10, 3.3, 0.3), list(6, 3.5, 16, 10, 6.3, 2.3))
param <- mapply("/", param, list(100, 100, 100, 100, 1, 1))
h <- expand.grid(param)

comb <- array(dim = c(6, 1, dim(h)[1]))                #All possible param combinations
comb[1:dim(comb)[1],,] <- t(h)                         #Last four parameters to be determined by grid search

#Retrieving squared differences for all param combinations and sorting
SSQRA <- GSS1(x = comb, m = data$term, r = data$ytm)

#finding out the combination with the lowest SSQ
cmatrix <- comb[,,apply(SSQRA, 1, which.min)] %>% data.frame %>%  mutate(across(, ~ifelse(.==0, 1e-6, .)))
cmatrix <- split(unlist(cmatrix), 1:6)

param_2 <- mapply(c, param, cmatrix)
param_2 <- data.frame(unlist(cmatrix), do.call(rbind, lapply(param_2, function(x) head(sort(diff(sort(x))),2)))) %>% 
  rename_all(~c("central", "lower", "upper"))

range_2 <- param_2 %>% rowwise() %>% mutate(max = max(across(2:3)))
range_2 <- param_2$central + cbind(- range_2$max,  range_2$max)


# Determination of the initial values to start non linear optimization: second grid search

#Creation of a narrower set of values around best values obtained in the first grid search
spread <- apply(range_2, 1, diff)/(lengths(param) - 1)
new_set <- mapply("*", apply(range_2/spread, 1, function(x) Reduce(seq, x)), as.list(spread))
hh <- expand.grid(new_set)

comb2 <- comb                                     #First 2 parameters do not change
comb2[1:dim(comb2)[1], , ] <- t(hh)               #Next 4 param are given by the combination matrix

SSQRB <- GSS1(x = comb2, m = data$term, r = data$ytm)
print(min(SSQRB)/min(SSQRA)<1)                    #Check that new parameters helped minimize SSQ

cmatrix_2 <- comb2[, , apply(SSQRB, 1, which.min)] %>% data.frame %>%  mutate(across(, ~ifelse(.==0, 1e-6, .)))
cmatrix_3 <- split(unlist(cmatrix_2), 1:6)

param_3 <- mapply(c, new_set, cmatrix_3)
param_3 <- data.frame(unlist(cmatrix_3), do.call(rbind,lapply(param_3, function(x) head(sort(diff(unique(sort(x)))), 2)))) %>% 
  rename_all(~c("central", "lower", "upper"))

range_3 <- param_3 %>% rowwise() %>% mutate(max = max(across(2:3)))
range_3 <- param_3$central + data.frame(lower = -range_3$max, start = 0, upper = range_3$max)


#Calibration of the 6 parameters
CI <- c(range_3$lower, -range_3$upper)
UI <- rbind(diag(6), -diag(6))

para_nlls <- CV <- error <- list()
for (i in 1:ncol(data[,-1])){
  sol <- constrOptim(cmatrix_2[,i], objective, NULL, ui = UI, ci = CI, mu = 1e-05, 
                     control = list(iter.max = 2000), method = "Nelder-Mead")
  para_nlls[[i]] <- sol$par
  CV[[i]] <- sol$convergence
  error[[i]] <- objective(x =  para[[i]])
}

e1 <- unlist(error)

para_nlls <- data.frame(matrix(do.call(rbind, para_nlls), nrow = 1, ncol = 6)) %>% rename_with(~colnames(coeff))
curve_nlls <- data.frame(matu = data$term, rates = NSS(para_nlls, data$term))

#graph
ggplot() +
  geom_point(data = curve_nlls, aes(matu, rates, color = "fitted rates"), size = 1) +
  geom_point(data = data, aes(term, ytm, color = "market rates"), size = 1) +
  geom_line(data = curve_nlls, aes(matu, rates, color = "fitted rates"), size = 0.5) +
  geom_line(data = data, aes(term, ytm, color = "market rates"), size = 0.5) +
  labs(y = "ytm (%)", x = "term (years)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))


##########################          LINEARIZING THE MODEL : GRID SEARCH OLS        ###################

#combinations of values of lambda_1 and lambda_2
lambda_ols <- seq(1, 15, 0.01)
lambda_ols <- expand.grid(lambda_ols, lambda_ols) %>% rename_with(~c("lambda_1", "lambda_2"))
  
#values of the different factor loadings for the different values of lambdas for each maturity
load_2 <- lapply(lambda_ols$lambda_1, function(x) loading_1(x, data$term))
load_3 <- lapply(lambda_ols$lambda_1, function(x) loading_2(x, data$term))
load_4 <- lapply(lambda_ols$lambda_2, function(x) loading_2(x, data$term))

#Grid search based OLS
RSS <- list()     #the RSS for each combination of lambdas in the linear model OLS estimate
for (i in 1:nrow(lambda_ols)){
    RSS[[i]] <- summary(lm(data$ytm ~ load_2[[i]] + load_3[[i]] + load_4[[i]]))$sigma
  }

lambda_opt <- lambda_ols[which.min(RSS), ]   #Index positioning of optimal lambda values

#linearized model with the best goodness of fit
reg <- data.frame(y = data$ytm, x_2 = load_2[[which.min(RSS)]],
                  x_3 = load_3[[which.min(RSS)]], x_4 = load_4[[which.min(RSS)]])

model <- lm(y ~ x_2 + x_3 + x_4, data = reg)

#associated estimated beta parameters with lambda parameters
para <- data.frame(matrix(c(summary(model)$coefficients[,1], unlist(lambda_opt)), nrow =1, ncol = 6)) %>% 
  rename_with(~colnames(coeff))
curve <- data.frame(matu = data$term, rates = NSS(para, data$term))

#graph
ggplot() +
  geom_point(data = curve, aes(matu, rates, color = "fitted rates"), size = 1) +
  geom_point(data = data, aes(term, ytm, color = "market rates"), size = 1) + 
  geom_line(data = curve, aes(matu, rates, color = "fitted rates"), size = 0.5) +
  geom_line(data = data, aes(term, ytm, color = "market rates"), size = 0.5) + 
  labs(y = "ytm (%)", x = "term (years)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))


################            LINEARIZING THE MODEL: FIXING THE LAMBDAS             ####################


#which values of lambda_1 generate a low correlation between the two factor loadings and will be tested
lambda_target <- lambda[ abs(correl$correlation) < 0.8]

lam_1 <- tail(lambda_target,  trunc(length(lambda_target)/2))  #tested values for lambda_1
lam_2 <- head(lambda_target,  trunc(length(lambda_target)/2))  #tested values for lambda_2

lambda_fix <- expand.grid(lam_1, lam_2) %>% rename_with(~c("lambda_1", "lambda_2"))

#values of the different factor loadings for the different values of lambdas
load_2 <- lapply(lambda_fix$lambda_1, function(x) loading_1(x, data$term))
load_3 <- lapply(lambda_fix$lambda_1, function(x) loading_2(x, data$term))
load_4 <- lapply(lambda_fix$lambda_2, function(x) loading_2(x, data$term))

#OLS
RSS_fixed <- list()
for (i in 1:nrow(lambda_fix)){
  RSS_fixed[[i]] <- summary(lm(data$ytm ~load_2[[i]] + load_3[[i]] + load_4[[i]]))$sigma
}

lambda_fix_opt <- lambda_fix[which.min(RSS_fixed), ]

reg_fix <- data.frame(y = data$ytm, x_2 = load_2[[which.min(RSS_fixed)]], x_3 = load_3[[which.min(RSS_fixed)]],
                        x_4 = load_4[[which.min(RSS_fixed)]])

model_fix <- lm ( y ~ ., data = reg_fix)

para_fix <- data.frame(matrix(c(summary(model_fix)$coefficients[ ,1], unlist(lambda_fix_opt)), nrow = 1, ncol = 6))
colnames(para_fix) <- colnames(coeff)
curve_fix <- data.frame(matu = data$term, rates = NSS(para_fix, data$term))

#graph
ggplot() +
  geom_point(data = curve_fix, aes(matu, rates, color = "fitted rates"), size = 1) +
  geom_point(data = data, aes(term, ytm, color = "market rates"), size = 1) +
  geom_line(data = curve_fix, aes(matu, rates, color = "fitted rates"), size = 0.5) +
  geom_line(data = data, aes(term, ytm, color = "market rates"), size = 0.5) +
  labs(y = "ytm (%)", x = "term (years)") + scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))