library(quadprog)
library(wpp2019)
library(tidyverse)
library(gridExtra)
source("R/QOSplit.R")

# main_funtion

# age survivor spline
lx_fun <- function(data, age){
  splinefun(data$x, data$lx, method = "monoH.FC")(age)
}

# function to calculate all the measures
CS_results <- function(LA, country, period, a){   # LA is the data base
  
  # example
  # country = "Brazil"; period = "1970-1975"; a = 70
  
  # filter country
  country = LA %>% filter(name == country, Period==period)%>% 
                    mutate(dx = c(-diff(lx),lx[101]),
                           Lx = ifelse(x==0, Lx_n, # assume uniform deaths within interval
                                       ifelse(x==100, ex, 
                                              lx-dx/2)),
                           Lx = as.numeric(Lx), # in case
                           ex = rev(cumsum(rev(Lx)))/lx)
  
  # calculate: 
    # assume that birth is in exact age (not in the middle). That allow us to work with lx, instead of Lx.
    # consider that i=a is age a-1.
  CS = sum(unlist(sapply(1:a, function(i) country$asfr[i]*country$lx[a+2-i]/100000)), na.rm = T)
  Fa = sum(country$asfr[1:a], na.rm = T)
  CS_prob = CS/Fa
  k = sum(country$asfr[1:a]*0:(a-1), na.rm = T)/Fa # must be kappa but...
  CS_aprox = Fa * country$lx[a+1-k]/100000
  error_aprox = (CS-CS_aprox)/CS*100 # relative error approx
  MTSL = 0
  for(i in 1:a){
    for(j in 1:(a+1-i)){
      MTSL = sum(MTSL,
                 country$asfr[i] * country$dx[j]/100000 * sum(country$Lx[j:a])/country$lx[j],
                 na.rm=T)
    }
  }
  MTSA = 0
  for(i in 1:a){
    for(j in 1:(a+1-i)){
      MTSA = sum(MTSA,
                 country$asfr[i] * sum(country$Lx[j:a])/country$lx[j],
                 na.rm=T)
    }
  }
  ITL = MTSL/MTSA * 100
  m_MAD = 0
  for(i in 1:a){
    for(j in 1:(a+1-i)){
      m_MAD = sum(m_MAD,
                  country$asfr[i] * country$dx[j]/100000 * (j-1),
                  na.rm=T)
    }
  }
  MAL = k + m_MAD/Fa
  k_CS = sum(unlist(sapply(1:a, function(i) country$asfr[i]*country$lx[a+2-i]/100000*(i-1))), na.rm = T)/CS
  abs_change = -(a-k_CS)
  
  # Results
  return(data.frame(CS,Fa,CS_prob,k,CS_aprox,error_aprox,MTSL,ITL,MAL))
}

# Get wpp locations
data(UNlocations)

# countries of LA
paises_LA_seleccion <- c(904, # promedio regional
                         192,214,332,388,630,188,222,320,340,484,558,591,32,68,76,152,170,218,600,604,858,862)

# get fertility data
data(tfr)
data(percentASFR)
# tfr
tfr_paises_LA <- tfr %>% select(-last.observed) %>% 
                  gather(Period,tfr,-country_code,-name) %>% 
                  filter(country_code %in% paises_LA_seleccion)
# asfr join tfr
asfr_paises_LA <- percentASFR %>% gather(Period, asfr, -name, -country_code, -age) %>% 
                  mutate(x = as.numeric(substr(age,1,2)), asfr=asfr/100) %>% 
                  left_join(tfr_paises_LA, by=c("country_code","name","Period")) %>% 
                  mutate(asfr = asfr/5 * tfr) %>% 
                  select(-age) %>% 
                  filter(country_code %in% paises_LA_seleccion)

# lt must downloaded manually (maybe too heavy for CRAN)
lt_paises_LA <- readxl::read_xlsx("Data/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", 
                                  sheet = 1, range = "A17:T77381") %>% 
                  rename(country_code=5,
                        x = "Age (x)", 
                        mx = "Central death rate m(x,n)",
                        qx = "Probability of dying q(x,n)",
                        dx = "Number of deaths d(x,n)",
                        lx = "Number of survivors l(x)",
                        Lx_n = "Number of person-years lived L(x,n)",
                        ex = "Expectation of life e(x)")  %>% 
                  select(country_code,Period,x, mx, qx, dx, lx, Lx_n, ex) %>% 
                  filter(country_code %in% paises_LA_seleccion)

# smooth both components
# initial df
asfr_paises_LA_s <- data.frame(country_code=0, Period = " ", x = 0, asfr = 0, stringsAsFactors = F)
lx_paises_LA_s   <- data.frame(country_code=0, Period = " ", x = 0, lx = 0,   stringsAsFactors = F)
# fill it
countries = unique(asfr_paises_LA$country_code)
for(c in 1:length(countries)){ # loop countries
  for(p in unique(asfr_paises_LA$Period)){ # loop periods
    # p = "1950-1955"; c = 1
    cp = asfr_paises_LA %>% filter(country_code==countries[c], Period==p) %>% select(asfr)
    if(any(is.na(cp$asfr))){next}
    out <- tryCatch(QOSplit(c(0,cp$asfr), L=seq(10,45,5), AgeInt=rep(5,8)),
                    error = data.frame(Age=10:49, ASFR=rep(NA,40)))
    asfr_paises_LA_s <- rbind(asfr_paises_LA_s, 
                              data.frame(country_code = countries[c], Period = p, 
                                         x = out$Age, asfr = out$ASFR))
    cp <- lt_paises_LA %>% filter(country_code==countries[c], Period==p) %>% select(x,lx)
    out <- lx_fun(cp, age=0:100)
    lx_paises_LA_s <- rbind(lx_paises_LA_s, data.frame(country_code=countries[c], 
                                                       Period = p, x = 0:100, lx = out))
  }
}

# join everything
# assuming
im = 100/204
# db with smooth components
LA <- lx_paises_LA_s %>% slice(-1) %>% mutate(country_code=as.integer(country_code)) %>% 
                          left_join(asfr_paises_LA_s %>% slice(-1), by=c("country_code","Period","x")) %>%
                          left_join(tfr_paises_LA, by=c("country_code","Period")) %>% 
                          left_join(lt_paises_LA %>% 
                                      select(country_code, Period, x, Lx_n, ex),
                                    by = c("country_code","Period","x")) %>% 
                          mutate(name = as.character(name),
                                 tfr = tfr*im,
                                 asfr = asfr*im) %>% 
                          left_join(
                            lx_paises_LA_s %>% slice(-1) %>% group_by(country_code, Period) %>% 
                              mutate(dx=c(-diff(lx),min(lx)),
                                     qx=dx/lx) )%>% ungroup() %>% 
                          mutate(Year = as.integer(substr(Period, 1, 4)) + 2.5)

# Bolivia name
LA$name[LA$name == "Bolivia (Plurinational State of)"] = "Bolivia"

# get results for every country/period/age 
CS_outs <- as.data.frame(expand.grid(name = unique(LA$name),
                                     Age = c(30,40,50,70),
                                     Period = unique(LA$Period))) %>% 
                                     mutate(Year = as.integer(substr(Period, 1, 4)) + 2.5)
for(i in 1:nrow(CS_outs)){
        CS_outs[i,5:13] = CS_results(LA, country = CS_outs[i,1], 
                                      period = CS_outs[i,3], 
                                      a = CS_outs[i,2])}

# CS abs and relative change

# biiiig loop, could be optimized with map/purr
age <- 50 # choose age 50, but could be 30,40,50,70
deltas <- c(.0001, .0001, 0.005, 0.01)
CS_change_app <- as.data.frame(expand.grid(
                                name = unique(LA$name),
                                Age = age,
                                Period = unique(LA$Period),
                                delta = deltas))
CS_fun = function(Age, mx, lx){
  sum(unlist(sapply(1:Age, function(j)
    mx[j] * lx[Age+2-j]/100000)),na.rm = T)
}
for(i in 1:nrow(CS_change_app)){
        
        # i = 1
  
        # loop
        country = CS_change_app[i,1]
        period = CS_change_app[i,3]
        age = CS_change_app[i,2]
        delta = CS_change_app[i,4]
        
        # country loop
        country_i = LA %>% 
                        filter(name == country, Period == period) %>%
                        mutate(Lx = ifelse(x==0, Lx_n, # assume uniform deaths within interval
                                           ifelse(x==100, ex,
                                                  lx-dx/2)),
                               Lx = as.numeric(Lx),
                               Mx = dx/Lx)
        CS = CS_outs %>% filter(name == country, Period == period, Age == age) %>% select(CS)
        kappa = CS_outs %>% filter(name == country, Period == period, Age == age) %>% select(k)

        # observed changes
        CS_changed_abs = sum(unlist(sapply(1:age, function(j)
                              country_i$asfr[j] * country_i$lx[age+2-j] * exp(-delta*(age+2-j)) /100000)),
                              na.rm = T)
        change_obs_abs = (CS_changed_abs - CS)/CS
    
        CS_changed_rel = sum(unlist(sapply(1:age, function(j)
                              country_i$asfr[j] * (country_i$lx[age+2-j]/100000)^(1+delta))),
                              na.rm = T)
      
        change_obs_rel = (CS_changed_rel - CS)/CS

        # analytical
        change_app_abs = -(age - kappa) * delta
        change_app_rel = 0


        
        change_app_rel = sum(sapply(1:age, function(m){
          -country_i$Mx[m] * CS_fun(age - m + 2, country_i$asfr, country_i$lx)
                        }))/CS * delta

        # precision
        precision_abs = (change_app_abs-change_obs_abs)/change_obs_abs * 100
        precision_rel = (change_app_rel-change_obs_rel)/change_obs_rel * 100

        # collect
        CS_change_app[i,5:10] = data.frame( change_obs_abs,
                                            change_obs_rel,
                                            change_app_abs,
                                            change_app_rel,
                                             precision_abs,
                                             precision_rel)
}

CS_change_app <- CS_change_app %>% rename(change_obs_abs = 5,
                                          change_obs_rel = 6,
                                          change_app_abs = 7,
                                          change_app_rel = 8,
                                          precision_abs = 9,
                                          precision_rel = 10)

# end
save(LA,
     CS_outs,
     im,
     CS_change_app,
     file = "R/CS_outs.RData")
