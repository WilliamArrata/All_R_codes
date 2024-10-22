
#####################   WILLIAM ARRATA - ESSEC PORTFOLIO MANAGEMENT COURSE WINTER 2023   #####################

require("pacman")
pacman::p_load("readxl", "ggplot2", "tibble", "dplyr", "data.table", "stringr")

#############################################  FEASIBLE SET #######################################################

#load data and calculate annualized historical returns and covariances
ret <- read_excel("stock_prices.xlsx") %>%  select_if(is.numeric) %>% mutate_all(~ ( (.) - shift(.))/(.)) %>% na.omit()
moy <- 252*colMeans(ret)                                             # Annualized expected returns
cov_mat <- 252*cov(ret)                                              # Annualized covariances

w <- list(seq(from = 0, to = 1, by = 0.1))                           # Create sets of weights for the 6 assets
w <- expand.grid(rep(w, length(moy)))
w <- as.matrix(w[rowSums(w)==1,])                                    # Weights sum to 1

means <- rowSums(w%*%moy)                                            # expected returns for each portfolio
stdev <- diag(sqrt(w%*%cov_mat%*%t(w)))                              # stddev for each portfolio
portfolio <- data.frame(Mean = means, Risk = stdev)                  # portfolio's coordinates

#Graph in the mean standard deviation space of all created portfolios
coord <- data.frame(mean = moy, stddev = sqrt(diag(cov_mat))) %>% arrange(mean)

ggplot(portfolio, aes(x = Risk, y = Mean)) + geom_point(color="indianred", size = 0.5) +
  geom_point(data = coord, aes(stddev, mean)) +
  annotate("text", x = coord$stddev, y = coord$mean, hjust = c(1.1, 0.4, rep(-0.1, 3), 1), size = 4, 
           vjust = c(0.5, 0.2)[c(1,2,2,2,1,1)], 
           label = as.expression(sapply(split(coord, coord$mean), function(x)
             bquote(mu == .(round(x[[1]]*100, 1)) ~ "% ; "~ sigma == .(round(x[[2]]*100, 1)) ~ "%")))) +
  scale_y_continuous(labels = scales::percent) +  scale_x_continuous(labels = scales::percent) +
  labs(x = 'Standard Deviations', y = 'Expected Returns')