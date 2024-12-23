
##################  EFFICIENT FRONTIER - THE CLOSED FORM SOLUTION WITH SHORT SELLING  ###################

require("pacman")
pacman::p_load("stringr","stats","readxl","data.table","zoo","dplyr","tidyr", "janitor", "tibble",
               "lubridate", "zoo", "textclean", "ggplot2", "writexl")

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

std <- sqrt(diag(t(all_w)%*%sig%*%all_w))   #stddev for each portfolio

coord <- bind_cols(targets, std) %>% rename_with(~c("mu", "sigma"))

ggplot() + geom_segment(data = coord, aes( x = std, xend = dplyr::lead(std), y = mu, yend = dplyr::lead(mu)),
                       size = 1.5, color="plum") +
  scale_y_continuous(labels = scales::percent) +  scale_x_continuous(labels = scales::percent) +
  labs(x = "standard deviation", y = "expected return") + theme(plot.margin = margin(.8,.5,.8,.5, "cm")) 
