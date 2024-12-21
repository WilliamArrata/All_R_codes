
##################   WILLIAM ARRATA 2024 -  EXPOSURES TO COMMODITY FUTURES CONTRACTS   ################

require("pacman")
pacman::p_load("writexl", "varhandle", "readxl", "stringr", "data.table", "nloptr", "textreadr", 
               "dplyr", "ggplot2", "tidyr")

prices <- read_excel("commodities_futures.xlsx", 1) %>% slice(-1) %>% mutate_all(as.numeric) %>% 
  rename_at(1, ~"Date") %>% mutate_at("Date", ~as.Date (., origin = '1899-12-30')) %>% rename_at( ncol(.), ~"r_f") %>% 
  mutate_at("r_f", ~./36500)

matu <- read_excel("commodities_futures.xlsx", 2) %>% rename_with(~c("contract", "matu", "start")) %>% 
  mutate_at(-1, ~as.Date(., format = "%d/%m/%Y"))
  

################################           UNFUNDED POSITION              ##############################

#here, the number of contracts is not adjusted according to the change in the future price value
#thus the collateral is not equal to the price of a future contract anymore through time, position is either underfunded or overfunded

unfunded <- prices %>% select(c("Date", "CLG5", "r_f")) %>% rename(fut_v = CLG5) %>% mutate(p_l = fut_v - shift(fut_v)) %>% 
  replace(is.na(.), 0) %>% mutate(capi_collat = first(fut_v)*cumprod(1 + r_f) , 
                                  capi_p_l = cumsum(p_l*cumprod(1 + r_f)/ (1 + first(r_f) ) ),  Collat = capi_collat + capi_p_l)

ggplot(unfunded) + geom_line(aes( x = Date, y = fut_v, color = "January 2025 WTI future") ) +
  geom_line(aes( x = Date, y = Collat, color = "collateralized position") ) + labs(y = "value (USD)", x = "Date") + 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))

################################           FULLY FUNDED POSITION          ##############################
         
#here, the number of contracts is dynamically adjusted according to the change in the future price value
#thus the collateral is always equal to the price of a future contract anymore through time

fully_funded <- prices %>% select(c("Date", "CLG5", "r_f")) %>% rename(fut_v = CLG5) %>%
  mutate(adj_fut =  1 + r_f*fut_v/shift(fut_v) ) %>% replace(is.na(.), 1) %>% mutate_at("adj_fut", ~cumprod(.) ) %>% 
  mutate(fut_v_new = fut_v *adj_fut) %>% mutate(p_l = fut_v_new - shift(fut_v_new)) %>% replace(is.na(.), 0) %>% 
  mutate(capi_collat = first(fut_v_new)*cumprod(1 + r_f) , capi_p_l = cumsum(p_l*cumprod(1 + r_f)/ (1 + first(r_f) ) ),  
         Collat = capi_collat + capi_p_l)

ggplot(fully_funded) + geom_line(aes( x = Date, y = fut_v, color = "January 2025 WTI future") ) +
  geom_line(aes( x = Date, y = Collat, color ="collateralized position") ) + labs(y = "value (USD)", x = "Date") + 
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))


####################           CONSTANT TERM FUTURES FULLY FUNDED             ###########################

#two contracts are combined and their relative weights adjusted to maintain the term of the portfolio constant

constant_term <- prices %>% select(c("Date", "CLG4", "r_f")) %>% rename(fut_v = CLG4) %>%  
  mutate(new_fut = 1 + r_f*fut_v/shift(fut_v) ) %>% replace(is.na(.), 1) %>% mutate(cum_fut = cumprod(new_fut) ) %>% 
  mutate(fut_v_new = fut_v *cum_fut) %>% mutate(p_l = fut_v_new - shift(fut_v_new)) %>% replace(is.na(.), 0) %>% 
  mutate(capi_collat = first(fut_v_new)*cumprod(1 + r_f) , capi_p_l = cumsum(p_l*cumprod(1 + r_f)/ (1 + first(r_f) ) ),  
         Collat = capi_collat + capi_p_l) %>% select(c("Date", "Collat", "fut_v")) %>% rename(fut_jan_24 = fut_v, collat_jan_24 = Collat) %>% 
  left_join(ff_ffunded) %>% rename(fut_jan_25 = fut_v, collat_jan_25 = Collat) %>% select(c(Date, fut_jan_24, collat_jan_24, fut_jan_25, collat_jan_25)) %>% 
  mutate(term_jan_25 = pmax(0, as.numeric((as.Date("2025-01-21") - Date)/365)), term_jan_24 = pmax(0, as.numeric((as.Date("2024-01-22") - Date)/365)) )

#A 2-year term is achieved using the WTI Jan 24 contract and the WTI Jan 25 contract

constant_term <- constant_term %>% mutate(q_jan_25 = (2 - term_jan_24)/(term_jan_25 - term_jan_24), q_jan_24 = (term_jan_25 - 2)/(term_jan_25 - term_jan_24) ) %>% 
  mutate( term = q_jan_25*term_jan_25  + q_jan_24*term_jan_24) %>% mutate(collat_2y = q_jan_25*collat_jan_25 + q_jan_24*collat_jan_24) %>% 
  filter(q_jan_25 >= 0 & q_jan_24 >= 0 ) %>% left_join(ff_ffunded)

ggplot(constant_term) + geom_line(aes( x = Date, y = collat_2y, color = "2-year position") ) +
  geom_line(aes( x = Date, y = Collat, color ="January 2025 future") ) + labs(y = "value (USD)", x = "Date") + 
  theme(legend.position = "bottom", legend.title=element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm"))
