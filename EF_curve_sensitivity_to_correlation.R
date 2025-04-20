##########################   WILLIAM ARRATA 2023 - ESSEC PORTFOLIO MANAGEMENT COURSE   ################################

require("pacman")
pacman::p_load("stringr", "stats", "data.table", "zoo", "dplyr", "tidyr", "janitor", "tibble",
               "lubridate", "ggplot2", "roll")


##################  EFFICIENT FRONTIER - SENSITIVITY OF PORTFOLIOS DIVERSIFICATION TO CORRELATION  ########################

#Sensitivity of the weights of portfolios on the efficient frontier to correlation coefficients between risky assets
n <- 7
set.seed(123)
mu <- abs(rnorm(7))*1e-1   #the vector of expected returns

#1. Impact of correlation coefficients on GMVP weights

#case of zero correlation
sig <- diag(n)
diag(sig) <- runif(n)/100    #the covariances matrix
gmvp <- c(solve(sig)%*%matrix(1, nrow = n))/
  (matrix(1, ncol = n)%*%solve(sig)%*%matrix(1, nrow = n))

#case of perfect positive correlation
sig_2 <- diag(n)
diag(sig_2) <- runif(n)/100
cov <- sqrt(expand.grid(diag(sig_2), diag(sig_2))) %>% data.frame %>% rename_with(~c("sig_i", "sig_j")) %>% 
  mutate(cov = sig_i*sig_j) %>% select(cov) %>% c
cov <- cov$cov
dim(cov) <- c(n, n) 

#The covariances matrix is not invertible!
gmvp_2 <- c(solve(cov)%*%matrix(1, nrow = n))/(matrix(1, ncol = n)%*%solve(cov)%*%matrix(1, nrow = n))

#case of perfect negative correlation
mat_2 <- -cov
diag(mat_2) <- diag(cov)   #the covariances matrix
gmvp_3 <- c(solve(mat_2)%*%matrix(1, nrow = n))/
  (matrix(1, ncol = n)%*%solve(mat_2)%*%matrix(1, nrow = n))

compar_gmvp <- data.frame(asset = LETTERS[1:7], corr0 = gmvp, corrmoins1 = gmvp_3) %>% 
  pivot_longer(cols =! asset, names_to = "correlation", values_to = "weight") %>% 
  mutate_at("correlation", ~as.factor(-as.numeric(gsub("[^0-9-]", "", .))))

#graph of compared weights depending on correlation coefficients
ggplot() +  geom_bar(data = compar_gmvp, aes( x = asset, y = weight, fill = correlation ), 
                     position = "dodge", stat="identity") +
  scale_y_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
  guides(fill= guide_legend(title.position = "left", title = "weight for different correlation coefficients")) +
  theme(legend.position = "bottom",  plot.margin = margin(.8,.5,.8,.5, "cm"))


#2. Impact of correlation coefficients on a portfolio on the efficient frontier other than the GMVP
target_mu <- 0.12  #an example of target expected return

#case of zero correlation
#the four scalars obtained from FOCs
A <- colSums(solve(sig))%*%mu
B <- t(mu)%*%solve(sig)%*%mu
C <- sum(solve(sig))
D <- B*C - A^2

g <- (matrix(B*rowSums(solve(sig)), nrow = n) - matrix( A, nrow = n)*solve(sig)%*%mu)/matrix(D, nrow = n)
h <- (matrix(C, nrow = n)*solve(sig)%*%mu - matrix(A*rowSums(solve(sig)), nrow = n)) /matrix(D, nrow = n)

#optimal weights for the target expected return
w_ef <- g + h*target_mu


#case of perfect negative correlation
#the four scalars obtained from FOCs
A_3 <- colSums(solve(mat_2))%*%mu
B_3 <- t(mu)%*%solve(mat_2)%*%mu
C_3 <- sum(solve(mat_2))
D_3 <- B_3*C_3 - A_3^2

g_3 <- (matrix(B_3*rowSums(solve(mat_2)), nrow = n) - matrix( A_3, nrow = n)*solve(mat_2)%*%mu)/matrix(D_3, nrow = n)
h_3 <- (matrix(C_3, nrow = n)*solve(mat_2)%*%mu - matrix(A_3*rowSums(solve(mat_2)), nrow = n)) /matrix(D_3, nrow = n)

#optimal weights for the target expected return
w_ef_3 <- g_3 + h_3*target_mu


compar_ef <- data.frame(asset = LETTERS[1:7], corr0 = w_ef, corrmoins1 = w_ef_3) %>% 
  pivot_longer(cols =! asset, names_to = "correlation", values_to = "weight") %>% 
  mutate_at("correlation", ~as.factor(-as.numeric(gsub("[^0-9-]", "", .))))

#graph of compared weights depending on the correlation coefficient
ggplot() +  geom_bar(data = compar_ef, aes( x = asset, y = weight, fill = correlation ), 
                     position = "dodge", stat="identity") +
  scale_y_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
  guides(fill= guide_legend(title.position = "left", title = "weight for different correlation coefficients")) +
  theme(legend.position = "bottom",  plot.margin = margin(.8,.5,.8,.5, "cm"))


###########################  EFFICIENT FRONTIER - PORTFOLIO WEIGHTS ACROSS THE CURVE  ##########################


#the covariances matrix
sig <- matrix( c(0.16, 0, 0.24, 0, 1.44, 0, 0.24, 0, 0.64)/100, nrow = 3)

#the inverse of the covarianes matrix
sig_inv <- solve(sig)

#the expected return vector
mu <- ( c(2, 4, 3) )/100

#the four scalars obtained from FOCs
A <- colSums(sig_inv)%*%mu
B <- t(mu)%*%sig_inv%*%mu
C <- sum(sig_inv)
D <- B*C - A^2

g <- (matrix(B*rowSums(sig_inv), nrow = 3) - matrix( A, nrow = 3)*sig_inv%*%mu)/matrix(D, nrow = 3)
h <- (matrix(C, nrow = 3)*sig_inv%*%mu - matrix(A*rowSums(sig_inv), nrow = 3)) /matrix(D, nrow = 3)

#a example
target_mu <- 0.04

#optimal weights
w <- g + (h)*target_mu

#generalization
targets <- seq(0, 0.06, 1e-3)

all_w <- matrix(mapply("+", g, mapply("*", list(h), as.list(targets))), nrow =3)

#stddev for each portfolio
std <- sqrt(diag(t(all_w)%*%sig%*%all_w))

coord <- bind_cols(targets, std) %>% rename_with(~c("mu", "sigma"))

ggplot() + geom_segment(data = coord, aes( x = std, xend = dplyr::lead(std), y = mu, yend = dplyr::lead(mu)),
                        size = 1.5, color="plum") +
  scale_y_continuous(labels = scales::percent) +  scale_x_continuous(labels = scales::percent) +
  labs(x = "standard deviation", y = "expected return") + theme(plot.margin = margin(.8,.5,.8,.5, "cm")) 


all_w_2 <- bind_cols(LETTERS[1:3], all_w) %>% data.frame %>% rename_with(~c("asset", as.character(targets))) %>% 
  pivot_longer(cols =! asset, names_to = "exp_return_target", values_to = "weight") %>% arrange(exp_return_target) 

ggplot() +  geom_bar(data = all_w_2, aes( x = asset, y = weight, fill = exp_return_target ), 
                     position = "dodge", stat="identity") +
  scale_y_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 7)) +
  guides(fill= guide_legend(title.position = "left", title = "weight for different correlation coefficients")) +
  theme(legend.position = "none",  plot.margin = margin(.8,.5,.8,.5, "cm"))