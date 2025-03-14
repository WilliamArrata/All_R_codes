require("pacman")
pacman::p_load("writexl", "varhandle", "readxl", "stringr", "data.table", "dplyr", "janitor", "ggplot2", 
                 "tibble", "readxlsb", "purrr", "zoo", "tidyr", "PerformanceAnalytics", "stats", "scales",
               "tseries")

#BKXE filtering according to strategy' s criteria: sequential filtering
port <- read_excel("exemple_investment_project.xlsx", 1) %>% mutate(Weight = ifelse(is.na(Ticker), Name, NA)) %>% 
  fill(Weight, .direction = "down") %>% rename(sector = Weight) %>% filter(!is.na(Ticker)) %>% select(-Shares) %>% 
  mutate_at(-c(1:3), as.numeric) %>% select(-Price) %>%
  rename_at(-c(1:3), ~c("PEG", "EPS_g_1y", "PER", "FCF_2020", "FCF_2021", "FCF_2022", "FCF_g", "FCF_y-3", 
                        "OPM_y", "OPM_y-1", "PEG_new", "sales_g_5y", "Mkt_cap")) %>%
  group_by(sector) %>% mutate(OPM_sector = mean(OPM_y, na.rm = T)) %>% ungroup() %>% 
  mutate_at("sector", ~str_squish(gsub("\\(.*", "", .))) %>% mutate_at("Mkt_cap", ~.*1e-9) %>% 
  filter(EPS_g_1y >= 7 & EPS_g_1y <= 20, PEG < 1, sales_g_5y > 0, OPM_y > OPM_y-1, OPM_y > OPM_sector, FCF_g > 0) %>% 
  select_if( is.character)

#this for portfolio construction: returns at weekly frequency
returns_p <- read_excel("extract_prices_ptf.xlsx", 1) %>% row_to_names(2) %>% data.frame %>% mutate_all(as.numeric) %>% 
  mutate_at(1, ~as.Date(., origin = "1899-12-30")) %>% rename_at(1, ~"Date")

returns_p_1 <- returns_p %>% filter(Date < "2023-03-01") %>%  mutate_at(-1, ~ (. - shift(.))/(.)) %>% drop_na() %>% 
  select(-1) %>% as.matrix

optim_no_short <- portfolio.optim(returns_p_1, pm = 0.025/252, reshigh = rep(0.2, ncol(returns_p_1)))
weights <- optim_no_short$pw

port <- port %>% bind_cols(weights_p = weights)


#this for perf attribution - useful only if ptf securities are not all in the benchmark
returns_p_2 <- returns_p %>% filter(Date >= "2023-03-01") %>%  filter(Date %in% as.Date(c("2023-03-03", "2024-11-15") )) %>% 
  mutate_at(-1, ~ (. - shift(.))/(.)) %>% drop_na() %>% select(-1) %>% unlist()

port <- port %>% bind_cols(return = returns_p_2)

#portfolio return, at the security level
r_p <- sum(port$weights_p*port$return)

#strategy's benchmark
bck <- read_excel("exemple_investment_project.xlsx", 2) %>% filter( !is.na("free float market cap")) %>% 
  rename_at( c(5, 7:9), ~ c("weight_bck", "sector_name_1", "global_sector", "sector")) %>% 
  mutate_at("weight_bck", ~./sum(., na.rm = T)) %>% select( c(1, 2, 8, 9, 5) )

corres <- bck %>% select(c(global_sector, sector)) %>% unique

port <- left_join(port, corres) %>% 
  mutate(global_sector = ifelse(row_number() == 4, "Banks",  ifelse(row_number() %in% c(8, 11), "Health Care Equipment & Services", global_sector))) %>% 
  select(-sector) %>% relocate(global_sector, .before =weights_p)

returns_b <- read_excel("extract_prices_bck.xlsx", 1) %>% row_to_names(2) %>% data.frame %>% mutate_all(as.numeric) %>% 
  mutate_at(1, ~as.Date(., origin = "1899-12-30")) %>% rename_at(1, ~"Date") %>% filter(Date %in% as.Date(c("2023-03-03", "2024-11-15") )) %>% 
  mutate_at(-1, ~ (. - shift(.))/(.)) %>% drop_na() %>% select(-1) %>% unlist()

bck <- bck %>% bind_cols(return = returns_b) %>% select(-sector)

#benchmark return, at the security level
r_b <- sum(bck$weight_bck*bck$return, na.rm = T)

er <- r_p - r_b


#BF Model
miss_sec_bck <- port[which(!port$Ticker%in%bck$Ticker), ] %>% bind_cols(weight_bck = 0) %>% relocate(weight_bck, .before = return) %>% 
  relocate(weights_p, .after = return) %>% select(-Ticker)

bf_model <- left_join( bck[, -1], port[, -c(1, 3)]) %>% bind_rows(miss_sec_bck ) %>% data.frame %>% replace(is.na(.), 0) %>% 
  group_by(global_sector) %>% mutate(w_sector_bck = sum(weight_bck), w_sector_ptf = sum(weights_p), r_sector_bck = sum(weight_bck*return), r_sector_ptf = sum(weights_p*return) )
  
bf <- bf_model %>% select(-c(Name, return, weight_bck, weights_p)) %>% unique %>% arrange(global_sector) %>% 
  mutate(ae = (w_sector_ptf - w_sector_bck)*r_sector_bck, sse = w_sector_bck*(r_sector_ptf - r_sector_bck),
         ie = (w_sector_ptf - w_sector_bck)*(r_sector_ptf - r_sector_bck) ) 

sum(bf$w_sector_ptf*bf$r_sector_ptf) - sum(bf$w_sector_bck*bf$r_sector_bck)

sum(bf$ae) + sum(bf$sse) + sum(bf$ie)