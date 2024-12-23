
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   #####################

require("pacman")
pacman::p_load("readxl", "ggplot2", "tibble", "dplyr", "data.table", "stringr")

#############################################  FEASIBLE SET #######################################################

#load data and calculate annualized historical returns and covariances
ret <- read_excel("stock_prices.xlsx") %>%  select_if(is.numeric) %>% mutate_all(~ ( (.) - shift(.))/(.)) %>% na.omit()
moy <- 252*colMeans(ret)                                                            # Annualized expected returns
cov_mat <- 252*cov(ret)                                                             # Annualized covariances
assets <- data.frame(mean = moy, stddev = sqrt(diag(cov_mat))) %>% arrange(mean)    # Coordinates of the six assets

#construction of portfolios with multiple combinations of the six assets
w <- list(seq(from = 0, to = 1, by = 0.1))                                # Create sets of weights for the 6 assets
w <- expand.grid(rep(w, length(moy)))
w <- as.matrix( w[rowSums(w) == 1,] )                                     # Weights sum to 1
means <- w%*%moy                                                          # expected returns by portfolio
stdev <- sqrt(diag(w%*%cov_mat%*%t(w)))                                   # stddev by portfolio
portfolio <- data.frame(Mean = means, Risk = stdev)                       # coordinates by portfolio

#Graph in the mean standard deviation space of all created portfolios
ggplot() + geom_point(data = portfolio, aes(x = Risk, y = Mean), color = "indianred", size = 0.5) +
  geom_point(data = assets, aes(x = stddev, y = mean)) + labs(x = 'Standard Deviations', y = 'Expected Returns') +
  annotate("text", x = coord$stddev, y = coord$mean, hjust = c(1.1, 0.4, rep(-0.1, 3), 1), 
           vjust = c(0.5, 0.2)[c(1, 2, 2, 2, 1, 1)], size = 4, 
           label = as.expression( sapply( split(assets, assets$mean), function(x)
             bquote(mu == .(round(x[[1]]*100, 1)) ~ "% ; "~ sigma == .(round(x[[2]]*100, 1)) ~ "%")))) +
  scale_y_continuous(labels = scales::percent) +  scale_x_continuous(labels = scales::percent)