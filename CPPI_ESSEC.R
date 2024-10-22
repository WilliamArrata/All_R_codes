require("pacman")
pacman::p_load("readxl", "data.table", "dplyr", "tidyr", "ggplot2", "janitor", "stringr")

##############################  WILLIAM ARRATA - CPPI - ESSEC WINTER 2023  ####################################

risky <- read_excel("SXXE.xlsx") %>% rename_with(~c("Date", "price")) %>%   #Loading risky asset price history
  mutate_all(as.numeric) %>% mutate_at("Date", as.Date, origin='1899-12-30') %>% drop_na()

#I define the parameters
F <- 90                                                                #Floor at inception
m <- 4                                                                 #Multiplier
chg <- which(risky$Date == "2019-12-31")                                #Date when riskfree rates changes
F1 <- 1.2                                                              #Additional Floors if cliquet
F2 <- 1.1

####################################  BACKTESTING WITH HISTORICAL DATA  #######################################

#DAILY REBALANCING, MAXIMUM EXPOSURE, NO CLIQUET 

#Initialization
backt <- data.frame(Date = risky$Date, r = rep( c(2.75,2.25)/100, c(chg, nrow(risky) - chg)), Vpf = 100, 
                    Px = risky$price, brebal = NA, rf = rep(c(2.75, 2.25)/100, c(chg, nrow(risky) - chg))) %>%
  mutate(Floor = F*cumprod(ifelse(row_number()!=1, 1 + shift(r)*as.numeric(Date - shift(Date))/365, 1)) ) %>% 
  mutate(Cushion = Vpf - Floor, Vrisky = m*Cushion, Vrf = Vpf - Vrisky, nsh = Vrisky/Px)

#Iteration
for (i in 2:nrow(backt)){
  #values of risky asset and cushion before rebalancing
  backt$brebal[i] <- backt$Px[i]*backt$nsh[i-1]
  backt$Cushion[i] <- (backt$Vrf*(1 + backt$rf*diff(backt$Date)/365))[i-1] + (backt$brebal - backt$Floor)[i]
  #rebalancing actions
  backt$Vrisky[i] <- ifelse(backt$Cushion[i] > backt$Floor[i]/(m-1), backt$brebal[i],
                            ifelse(backt$Cushion[i] > 0, m*backt$Cushion[i], backt$Vrisky[i]))     
  backt$nsh[i] <- ifelse(backt$Cushion[i] > backt$Floor[i]/(m-1), backt$nsh[i-1],
                       ifelse(backt$Cushion[i] > 0, (backt$Vrisky/backt$Px)[i], 0))
  backt$Vrf[i] <- (backt$Vrf*(1 + backt$rf*diff(backt$Date)/365) )[i-1] + 
    ifelse(backt$Cushion[i] > backt$Floor[i]/(m-1), 0,
           ifelse(backt$Cushion[i] > 0, (backt$brebal - backt$Vrisky)[i], backt$Px[i]*backt$nsh[i-1]))
  backt$Vpf[i] <- ifelse(backt$Cushion[i] > 0, (backt$Vrisky + backt$Vrf)[i],
                         (backt$Vrf*(1 + backt$rf*diff(backt$Date)/365))[i-1] + backt$Px[i]*backt$nsh[i-1])
}

ggplot() + geom_line(data = backt, aes(x = Date, y = Vpf)) +  labs(x = "year", y = 'risky value') #graph


################       SIMULATIONS OF CPPI VALUES WITH THE GIVEN PARAMETERS      #################

#Simulation of 1000 paths for the price of the risky asset over 500 days
R <- 0.0225
t <- 0:500
sig <- sqrt(252)*sd(diff(risky$price)/risky$price[-last(risky$price)])     #historical volatility used for simu
nsim <- 1000
set.seed(123)

#Brownian Motion simulation
dW <- split(rnorm(n = nsim*(length(t) - 1), sd = sig), 1:nsim)                    #mean=0 thus no need to specify
W <- lapply(dW, function(x) c(0, cumsum(x)))                                      #integration and inception at 0
plot(NA, type = "l", xlim = range(t), ylim = range(W), xlab="time", ylab="brownian motion") #graph
lapply(W, lines)

#asset price simulation
S <- lapply(W, function(x) 100*exp( (R - 0.5*sig^2)*t + sig*x ) )
plot(NA, type = "l", xlim = range(t), ylim = range(S), xlab="time", ylab="asset value") #graph
lapply(S, lines)

#CPPI values

simu_tot <- list()

for (j in 1:length(S)){

#Initialization
    simu <- data.frame(Date = t, Px = S[[j]], Vpf = 100, brebal = NA,
                      Floor = F*c(1, cumprod(1 + as.numeric( tail(R, -1)*diff(t)/365)))) %>% 
    mutate(Cushion = Vpf - Floor, Vrisky = m*Cushion, Vrf = Vpf - Vrisky, nsh = Vrisky/Px)
  #Iteration
  for (i in 2:nrow(simu)){
    #values of risky asset and cushion before rebalancing
    simu$brebal[i] <- simu$Px[i]*simu$n[i-1]
    simu$Cushion[i] <- (simu$Vrf*(1 + R*diff(simu$Date)/365))[i-1] + (simu$brebal - simu$Floor)[i]
    #rebalancing actions
    simu$Vrisky[i] <- ifelse(simu$Cushion[i] > simu$Floor[i]/(m-1), simu$brebal[i],
                              ifelse(simu$Cushion[i] > 0, m*simu$Cushion[i], simu$Vrisky[i]))     
    simu$nsh[i] <- ifelse(simu$Cushion[i] > simu$Floor[i]/(m-1), simu$nsh[i-1],
                           ifelse(simu$Cushion[i] > 0, (simu$Vrisky/simu$Px)[i], 0))
    simu$Vrf[i] <- (simu$Vrf*(1 + R*diff(simu$Date)/365) )[i-1] + 
      ifelse(simu$Cushion[i] > simu$Floor[i]/(m-1), 0,
             ifelse(simu$Cushion[i] > 0, (simu$brebal - simu$Vrisky)[i], simu$Px[i]*simu$nsh[i-1]))
    simu$Vpf[i] <- ifelse(simu$Cushion[i] > 0, (simu$Vrisky + simu$Vrf)[i],
                           (simu$Vrf*(1 + R*diff(simu$Date)/365))[i-1] + simu$Px[i]*simu$nsh[i-1])
  }
  simu_tot[[j]] <- simu$Vpf
}

plot(NA, type="l", xlim = range(t), ylim = round(range(simu_tot)), xlab="time", ylab="risky value",)  #graph
lapply(simu_tot, lines)


###########################     FINDING OUT THE OPTIMAL MULTIPLIER VALUE      ############################

mul <- seq(2, 8, 0.1)     #trying different muyltipliers values

simu_tot <- list()

for (j in 1:length(S)){
  
  simu_tot[[j]] <- list()
  
  for (z in 1:length(mul)){
  
  #Initialization
  simu <- data.frame(Date = t, Px = S[[j]], Vpf = 100, brebal = NA,
                     Floor = F*c(1, cumprod(1 + as.numeric( tail(R, -1)*diff(t)/365)))) %>% 
    mutate(Cushion = Vpf - Floor, Vrisky = mul[z]*Cushion, Vrf = Vpf - Vrisky, nsh = Vrisky/Px)
  #Iteration
  for (i in 2:nrow(simu)){
    #values of risky asset and cushion before rebalancing
    simu$brebal[i] <- simu$Px[i]*simu$n[i-1]
    simu$Cushion[i] <- (simu$Vrf*(1 + R*diff(simu$Date)/365))[i-1] + (simu$brebal - simu$Floor)[i]
    #rebalancing actions
    simu$Vrisky[i] <- ifelse(simu$Cushion[i] > simu$Floor[i]/(m-1), simu$brebal[i],
                             ifelse(simu$Cushion[i] > 0, mul[z]*simu$Cushion[i], simu$Vrisky[i]))     
    simu$nsh[i] <- ifelse(simu$Cushion[i] > simu$Floor[i]/(mul[z]-1), simu$nsh[i-1],
                          ifelse(simu$Cushion[i] > 0, (simu$Vrisky/simu$Px)[i], 0))
    simu$Vrf[i] <- (simu$Vrf*(1 + R*diff(simu$Date)/365) )[i-1] + 
      ifelse(simu$Cushion[i] > simu$Floor[i]/(mul[z]-1), 0,
             ifelse(simu$Cushion[i] > 0, (simu$brebal - simu$Vrisky)[i], simu$Px[i]*simu$nsh[i-1]))
    simu$Vpf[i] <- ifelse(simu$Cushion[i] > 0, (simu$Vrisky + simu$Vrf)[i],
                          (simu$Vrf*(1 + R*diff(simu$Date)/365))[i-1] + simu$Px[i]*simu$nsh[i-1])
  }
  simu_tot[[j]][[z]] <- simu$Vpf
  }
}

#we draw one risky asset price trajectory and examine CPPI value for the different multiplier levels.
sc <- sample(length(S), 1)
print(sc)
sim <- mapply(c, mul, simu_tot[[sc]]) %>% row_to_names(1) %>% data.frame %>% rename_with(~paste0(., "mul")) %>% 
  bind_cols(t) %>% pivot_longer(cols = ends_with('mul'), names_to = 'multiplier', values_to = 'CPPI_value') %>% 
  mutate_at("multiplier", ~str_squish(gsub('[a-zA-Z]',' ', .))) %>% rename_at(1, ~"Date")

ggplot() + geom_line(data = sim, aes(x = Date, y = CPPI_value, color = multiplier)) + 
  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm"))