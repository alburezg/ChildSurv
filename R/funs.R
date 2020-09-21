library(quadprog)
library(wpp2019)
library(tidyverse)
library(gridExtra)
source("R/QOSplit.R")

# main functions ----------------------------------------------------------

# age survivor split spline
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
                           Lx = as.numeric(Lx), # just in case
                           ex = rev(cumsum(rev(Lx)))/lx)
  
  # calculate: 
    # assume that birth is in exact age (not in the middle). That allow us to work with lx, instead of Lx.
    # consider that i=a is age a-1.
  CS = sum(unlist(sapply(1:a, function(i) country$asfr[i]*country$lx[a+2-i]/100000)), na.rm = T)
  Fa = sum(country$asfr[1:a], na.rm = T)
  CS_prob = CS/Fa
  k = sum(country$asfr[1:a]*0:(a-1), na.rm = T)/Fa # should be noted as kappa but...
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


# get data and smooth -----------------------------------------------------

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

# lt data must be downloaded manually (maybe too heavy for CRAN)
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
# assuming at birth
im = 100/204
# db with smooth components
LA <- lx_paises_LA_s %>% slice(-1) %>% 
                          mutate(country_code=as.integer(country_code)) %>% 
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

# change Bolivia name for graphs
LA$name[LA$name == "Bolivia (Plurinational State of)"] = "Bolivia"

# get results for every country/period/age 
CS_outs <- as.data.frame(expand.grid(name = unique(LA$name),
                                     Age = c(30,40,50,70),
                                     Period = unique(LA$Period))) %>% 
                                     mutate(Year = as.integer(substr(Period, 1, 4)) + 2.5)
for(i in 1:nrow(CS_outs)){
        CS_outs[i,5:13] = CS_results(LA, country = CS_outs[i,1], 
                                      period = CS_outs[i,3], 
                                      a = CS_outs[i,2])
        print(round(i/nrow(CS_outs)*100,1))
        }


# Ch Ch Ch Change ---------------------------------------------------------
# CS abs and relative change
# biiiig loop, could be optimized with map/purr
age <- 50 # choose age 50, but could be 30,40,50,70
# proposed changes
deltas <- c(.0001, .001, 0.01, .1)
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
        change_app_rel = sum(sapply(1:age, function(m){
          -country_i$Mx[m] * CS_fun(age - m + 2, country_i$asfr, country_i$lx)
                        }))/CS * delta

        # precision
        precision_abs = (change_app_abs-change_obs_abs)/change_obs_abs * 100
        precision_rel = (change_app_rel-change_obs_rel)/change_obs_rel * 100

        # collect
        CS_change_app[i,5:13] = data.frame( change_obs_abs,
                                            change_obs_rel,
                                            change_app_abs,
                                            change_app_rel,
                                             precision_abs,
                                             precision_rel,
                                            CS,
                                            CS_changed_abs,
                                            CS_changed_rel)
        print(round(i/nrow(CS_change_app)*100,1))
}

CS_change_app <- CS_change_app %>% rename(change_obs_abs = 5,
                                          change_obs_rel = 6,
                                          change_app_abs = 7,
                                          change_app_rel = 8,
                                          precision_abs = 9,
                                          precision_rel = 10,
                                          CS = 11,
                                          CS_changed_abs = 12,
                                          CS_changed_rel = 13)



# mentions in text --------------------------------------------------------
LA %>% filter(x==0, name=="Latin America and the Caribbean") %>% select(qx, ex, Period)
CS_outs %>% filter(name=="Latin America and the Caribbean", Age==50) %>% 
            select(Age, Period, CS, Fa, CS_prob, k)


# get graphs --------------------------------------------------------------
# TFR & q0
library(ggpubr)

a <- ggplot(LA %>% filter(x==0, name!="Latin America and the Caribbean"), 
       aes(x=qx, y=tfr))+
      geom_point(aes(color=Year)) +
      geom_smooth(se = F, color = "grey", method = "lm",
                  formula = y ~ log(x), linetype = "dashed") +
      geom_line(data = LA %>% filter(x==0, name=="Latin America and the Caribbean"),
                aes(x=qx, y=tfr, color=Year), size=1)+
      theme_classic() + labs(x="q(0)", y="GRR") + 
      theme(legend.position = "none") 

# CS
b <- ggplot(CS_outs %>% filter(name!="Latin America and the Caribbean"), 
            aes(x=Year, y=CS, color=factor(Age)))+
      geom_point(alpha = .3) + 
      geom_line(data = CS_outs %>% filter(name=="Latin America and the Caribbean"), 
                aes(x=Year, y=CS, color=factor(Age))) +
      theme_classic() +
      labs(y = "CSa") + theme(legend.position = "none")
c <- ggplot(CS_outs %>% filter(name!="Latin America and the Caribbean"), 
            aes(x=Year, y=CS_prob, color=factor(Age)))+
      geom_point(alpha = .3) + 
      geom_line(data = CS_outs %>% filter(name=="Latin America and the Caribbean"), 
                aes(x=Year, y=CS_prob, color=factor(Age))) +
      theme_classic() +
      labs(y = "CSa/Fa") + 
      scale_colour_discrete(name="Age")

pdf("R/GRR_q0_CS.pdf")
ggarrange(a,                                                 
          ggarrange(b, c, ncol = 2, labels = c("B", "C")),           
          nrow = 2, labels = "A") 
dev.off()

# CS approx
CS_outs %>% filter(name=="Latin America and the Caribbean", Age==50) %>% 
  select(Age, Period, CS, CS_aprox, error_aprox)

age = 30
CS_app_30 <- ggplot(CS_outs %>% filter(Age == age, name!="Latin America and the Caribbean"), 
                 aes(x=CS, y=error_aprox, color=Year))+
  geom_point(alpha = .6) + 
  geom_point(data = CS_outs %>% filter(name=="Guatemala",Age == age),
             aes(x=CS, y=error_aprox), color="red") +
  geom_line(data = CS_outs %>% filter(Age == age, name=="Latin America and the Caribbean"), 
            aes(x=CS, y=error_aprox), size=2) +
  geom_hline(yintercept = 0, color="grey", linetype = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "%", x="CS30")+
  ylim(-2, 2)
age = 40
CS_app_40 <- ggplot(CS_outs %>% filter(Age == age, name!="Latin America and the Caribbean"), 
                    aes(x=CS, y=error_aprox, color=Year))+
  geom_point(alpha = .6) + 
  geom_point(data = CS_outs %>% filter(name=="Guatemala",Age == age),
             aes(x=CS, y=error_aprox), color="red") +
  geom_line(data = CS_outs %>% filter(Age == age, name=="Latin America and the Caribbean"), 
            aes(x=CS, y=error_aprox), size=2) +
  geom_hline(yintercept = 0, color="grey", linetype = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "%", x="CS40")+
  ylim(-2, 2)
age = 50
CS_app_50 <- ggplot(CS_outs %>% filter(Age == age, name!="Latin America and the Caribbean"), 
                    aes(x=CS, y=error_aprox, color=Year))+
  geom_point(alpha = .6) + 
  geom_point(data = CS_outs %>% filter(name=="Guatemala",Age == age),
             aes(x=CS, y=error_aprox), color="red") +
  geom_line(data = CS_outs %>% filter(Age == age, name=="Latin America and the Caribbean"), 
            aes(x=CS, y=error_aprox), size=2) +
  geom_hline(yintercept = 0, color="grey", linetype = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "%", x="CS50")+
  ylim(-2, 2)
Guate_lt <- LA %>% filter(name=="Guatemala") %>% ggplot() +
  geom_line(aes(x=x, y=lx/100000, color=Period)) +
  theme_classic() +
  scale_color_grey()+
  theme(legend.position = "none") +
  labs(y = "lx",x="Age")
Guate_mx <- LA %>% filter(name=="Guatemala", between(x, 15,49)) %>% ggplot() +
  geom_line(aes(x=x, y=asfr/tfr*100, color=Period)) +
  theme_classic() +
  scale_color_grey()+
  theme(legend.position = "none") +
  labs(y = "%",x="Age")


pdf("R/CS_app.pdf")
ggarrange(ggarrange(CS_app_30, CS_app_40, CS_app_50, ncol = 3, labels = c("A", "B", "C")),                                                 
          ggarrange(Guate_lt,Guate_mx,ncol=2,labels = c("D","E")),           
          nrow = 2) 
dev.off()

# changes approximation
CS_change_app %>% filter(name=="Latin America and the Caribbean") %>% 
                  select(Period, CS_changed_abs, CS,
                         change_obs_abs, change_app_abs, precision_abs, delta)
2.2365080-2.2413889


ch_abs <-  CS_change_app %>%
  ggplot() +
  geom_point(aes(x=change_obs_abs, y=change_app_abs, color=delta)) +
  geom_abline(slope=1, intercept = 0) +
  scale_x_continuous("Empirical") +
  scale_y_continuous("Approximated") +
  scale_color_continuous("Delta") +
  theme_classic()+
  theme(legend.position = "none") 

CS_change_app %>% filter(name=="Latin America and the Caribbean") %>% 
                  select(Period, CS_changed_rel, CS,
                         change_obs_rel, change_app_rel, precision_rel, delta)

ch_rel <-  CS_change_app %>%
  ggplot() +
  geom_point(aes(x=change_obs_rel, y=change_app_rel, color=delta)) +
  geom_abline(slope=1, intercept = 0) +
  scale_x_continuous("Empirical") +
  scale_y_continuous("") +
  scale_color_continuous("Delta") +
  theme_classic()

pdf("R/change_app.pdf")
ggarrange(ch_abs, ch_rel, labels = c("A","B"), ncol = 2) 
dev.off()

# S-shape
S_shape <- data.frame(country = " ",
                        age = 0,
                        CS_CS = 0)

for(country in c("Guatemala", "Argentina", "Haiti")){
  S_shape_i = data.frame(country = country, 
                         age = 15:Age,
                         CS_CS = sapply(15:Age, function(age)
                                       CS_results(country = country, 
                                               period = "1950-1955", 
                                               a = age, LA = LA)$CS)/
                                       CS_results(country = country, 
                                               period = "1950-1955", 
                                               a = Age, LA = LA)$CS)
  S_shape = rbind(S_shape, S_shape_i)
}

S_shape = S_shape[-1,]

g_shape <- ggplot(S_shape) + 
            geom_line(aes(age,CS_CS,color=country), size=2)+
            geom_hline(yintercept = 1, linetype="dashed")+
            theme_classic() +
            theme(legend.position = "right") +
            labs(y = "CSx/CS50", x="Age")+
            scale_color_discrete(name = "Country")
pdf("R/S_Shape.pdf")
g_shape
dev.off()


# ITL & MAL   
  ITL <- ggplot(CS_outs %>% filter(name!="Latin America and the Caribbean"), 
                  aes(x=Fa, y=ITL, color=factor(Age)), alpha = .3) + # alpha not working
      geom_point() +
      geom_line(data = CS_outs %>% filter(name=="Latin America and the Caribbean"), 
                aes(x=Fa, y=ITL, color=factor(Age)), size =1)+
      theme_classic() +
      theme(legend.position = "right") +
      labs(y = "%")+
      scale_color_discrete(name = "Age")
  
age = 30
MAL_30 <- ggplot(CS_outs %>% filter(Age==age, name!="Latin America and the Caribbean"))+
    geom_point(aes(x=Fa, y=MAL, color=Year), alpha = .3) +
    geom_line(data = CS_outs %>% filter(Age==age, name=="Latin America and the Caribbean"),
              aes(x=Fa, y=MAL, color=Year),size= 1) +
    theme_classic()+
    theme(legend.position = "none") +
    labs(y = "Mean Age at Lost")
age = 50
MAL_50 <- ggplot(CS_outs %>% filter(Age==age, name!="Latin America and the Caribbean"))+
  geom_point(aes(x=Fa, y=MAL, color=Year), alpha = .3) +
  geom_line(data = CS_outs %>% filter(Age==age, name=="Latin America and the Caribbean"),
            aes(x=Fa, y=MAL, color=Year),size= 1) +
  theme_classic()+
  theme(legend.position = "none") +
  labs(y = "Mean Age at Loss")


pdf("R/itl_mal.pdf")
ggarrange(ITL, MAL_50, labels = c("A","B"),ncol=2)
dev.off()


CS_outs %>% filter(name=="Latin America and the Caribbean", Age==30) %>% 
              select(MAL,Period,Age)
  
  
# end
save(LA,
     CS_outs,
     im,
     CS_change_app,
     file = "R/CS_outs.RData")
