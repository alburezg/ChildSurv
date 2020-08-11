# age survivor
lx_fun <- function(data, age){
  splinefun(data$x, data$lx, method = "monoH.FC")(age)
}
# p-y
Lx_fun <- function(data, agei, agen){
  integrate(lx_fun, agei, agen, data=data)
  }

# split asfr
install.packages("quadprog")
library(quadprog)
source("split/QOSplit.R")
int = QOSplit (Fx= c(0.00008, 0.02655, 0.10538, 0.09983, 0.05105, 0.01904, 0.00359,
               0.00033, 0.00000), L=seq(10,50,5), AgeInt=rep (5,9))

plot(L,Fx,t="s")
lines(int,col=2)

# data
library(tidyverse)
library(quadprog)
source("split/QOSplit.R")
library(wpp2019)

# Get wpp data
data(UNlocations)
data(tfr)
data(percentASFR)
paises_LA <- UNlocations %>% filter(area_name == "Latin America and the Caribbean") %>% 
              select(country_code, name) %>% distinct()
paises_LA_seleccion <- c(192,214,332,388,630,188,222,320,340,484,558,591,32,68,76,152,170,218,600,604,858,862)
lt_paises_LA <- readxl::read_xlsx("WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", sheet = 1, 
                range = "A17:T77381") %>% rename(country_code=5,
                                      x = "Age (x)", 
                                      mx = "Central death rate m(x,n)",
                                      qx = "Probability of dying q(x,n)",
                                      dx = "Number of deaths d(x,n)",
                                      lx = "Number of survivors l(x)",
                                      Lx_n = "Number of person-years lived L(x,n)",
                                      ex = "Expectation of life e(x)")  %>% 
  select(country_code,Period,x, mx, qx, dx, lx, Lx_n, ex) %>% 
  filter(country_code %in% paises_LA_seleccion)
tfr_paises_LA <- tfr %>% select(-last.observed) %>% gather(Period,tfr,-country_code,-name) %>% 
  filter(country_code %in% paises_LA_seleccion)
asfr_paises_LA <- percentASFR %>% gather(Period, asfr, -name, -country_code, -age) %>% 
  mutate(x = as.numeric(substr(age,1,2)), asfr=asfr/100) %>% 
  left_join(tfr_paises_LA, by=c("country_code","name","Period")) %>% 
  mutate(asfr = asfr/5 * tfr) %>% 
  select(-age) %>% 
  filter(country_code %in% paises_LA_seleccion)

# smooth
asfr_paises_LA_s <- data.frame(country_code=0, Period = " ", x = 0, asfr = 0, stringsAsFactors = F)
lx_paises_LA_s <- data.frame(country_code=0, Period = " ", x = 0, lx = 0, stringsAsFactors = F)
countries = unique(asfr_paises_LA$country_code)
for(c in 1:length(countries)){
  for(p in unique(asfr_paises_LA$Period)){
    # c = unique(asfr_paises_LA$country_code)[1]; c=32
    # p = "1950-1955"
    cp = asfr_paises_LA %>% filter(country_code==countries[c], Period==p) %>% select(asfr)
    if(any(is.na(cp$asfr))){next}
    out <- tryCatch(QOSplit(c(0,cp$asfr), L=seq(10,45,5), AgeInt=rep(5,8)),
                    error = data.frame(Age=10:49, ASFR=rep(NA,40)))
    asfr_paises_LA_s <- rbind(asfr_paises_LA_s, data.frame(country_code = countries[c], Period = p, x = out$Age, asfr = out$ASFR))
    cp = lt_paises_LA %>% filter(country_code==countries[c], Period==p) %>% select(x,lx)
    out <- lx_fun(cp, age=0:100)
    lx_paises_LA_s <- rbind(lx_paises_LA_s, data.frame(country_code=countries[c], Period = p, x = 0:100, lx = out))
    # plot(seq(10,45,5),c(0,cp$asfr),t="s");lines(out,col=2)
  }
  print(c)
}

LA <- lx_paises_LA_s %>% slice(-1) %>% mutate(country_code=as.integer(country_code)) %>% 
          left_join(asfr_paises_LA_s %>% slice(-1), by=c("country_code","Period","x")) %>%
          left_join(tfr_paises_LA, by=c("country_code","Period")) %>% 
          left_join(lt_paises_LA %>% 
                      select(country_code, Period, x, Lx_n, ex),
                    by = c("country_code","Period","x")) %>% 
          mutate(q = (1-lx/100000),
                 name = as.character(name),
                 tfr = tfr*1/2.04,
                 asfr = asfr*1/2.04)
LA$name[LA$name == "Bolivia (Plurinational State of)"] = "Bolivia"

# graph tfr paises
plot_tfr_q015 <- ggplot(LA %>% filter(x==15), aes(x=tfr, y=q, color=Period))+
                  geom_point() + geom_smooth(se = F, color = "grey") + 
                  theme_classic() + labs(x="TFR", y="q(0,15)")

# CS results
CS_results <- function(country, period, a){
  
  # filter country
  # a = 40; country = "Guatemala"; period = "1950-1955"
  country = LA %>% filter(name == country, Period==period)%>% 
    mutate(dx = c(-diff(lx),lx[101]),
           Lx = ifelse(x==0, Lx_n,
                       ifelse(x==100, ex, 
                              lx-dx/2)),
           Lx = as.numeric(Lx),
           ex = rev(cumsum(rev(Lx)))/lx)
  
  # calculate
  CS = sum(unlist(sapply(1:a, function(i) country$asfr[i]*country$lx[a-(i-1)]/100000)), na.rm = T)
  Fa = sum(country$asfr[1:a], na.rm = T)
  CS_prob = CS/Fa
  k = sum(country$asfr[1:a]*0:(a-1), na.rm = T)/Fa
  CS_aprox = Fa * country$lx[a-k+1]/100000
  error_aprox = (CS-CS_aprox)/CS*100
  
  MTSL = 0
  for(i in 1:a){
    for(j in 1:(a-i)){
      MTSL = sum(MTSL,
                 country$asfr[i] * country$dx[j]/100000 * sum(country$Lx[j:a])/country$lx[j],
                 na.rm=T)
    }
  }
  MTSA = 0
  for(i in 1:a){
    for(j in 1:(a-i)){
      MTSA = sum(MTSA,
                 country$asfr[i] * sum(country$Lx[j:a])/country$lx[j],
                 na.rm=T)
    }
  }
  ITL = MTSL/MTSA * 100
  m_MAD = 0
  for(i in 1:a){
    for(j in 1:(a-i)){
      m_MAD = sum(m_MAD,
                  country$asfr[i] * country$dx[j]/100000 * (j-1),
                  na.rm=T)
    }
  }
  MAL = k + m_MAD/Fa
  k_CS = sum(unlist(sapply(1:a, function(i) country$asfr[i]*country$lx[a-(i-1)]/100000*i)), na.rm = T)/CS
  abs_change = -(a-k_CS)
  
  # Results
  return(data.frame(CS,Fa,CS_prob,k,CS_aprox,error_aprox,MTSL,ITL,MAL))
}

# get results
CS_outs <- as.data.frame(expand.grid(name = unique(LA$name),
                                     Age = c(30,40,50),
                                     Period = unique(LA$Period)))
for(i in 1:nrow(CS_outs)){
  CS_outs[i,4:12] = CS_results(country = CS_outs[i,1], 
                               period = CS_outs[i,3], 
                               a = CS_outs[i,2])}
CS_outs <- CS_outs %>% mutate(Year = as.integer(substr(Period, 1, 4)) + 2.5)

# plot CS

plot_CS <- ggplot(CS_outs, aes(x=Year, y=CS_prob, color=factor(Age)))+
              geom_point(alpha = .3) + 
              geom_smooth(se = F) +
              theme_classic() +
              labs(y = "CS")
      
# plot Aprox

plot_CS_aprox <- ggplot(CS_outs, aes(x=CS, y=CS_aprox, color=factor(Period)))+
                    geom_point(alpha = .3) + 
                    geom_smooth(se = F) +
                    theme_classic()

# plot ITL

plot_ITL <- ggplot(CS_outs %>% filter(Age==50))+
                    geom_smooth(aes(x=Fa, y=ITL), col=1, se = F)+
                    geom_point(aes(x=Fa, y=ITL, color=factor(Period)), alpha = .3) +
                    theme_classic()

# plot MAL

plot_MAL <- ggplot(CS_outs %>% filter(Age==50))+
                    geom_smooth(aes(x=Fa, y=MAL), col=1, se = F)+
                    geom_point(aes(x=Fa, y=MAL, color=factor(Period)), alpha = .3) +
                    theme_classic()





# library(tidyverse)
# library(quadprog)
# source("split/QOSplit.R")
# library(wpp2019)
# 
# # Get wpp data
# data(UNlocations)
# data(tfr)
# data(percentASFR)
# paises_LA <- UNlocations %>% filter(area_name == "Latin America and the Caribbean") %>% select(country_code, name) %>% distinct()
# lt_paises_LA <- readxl::read_xlsx("WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", sheet = 1, range = "A17:T77381") %>% rename(country_code=5,
#                   x = "Age (x)", 
#                  mx = "Central death rate m(x,n)",
#                  qx = "Probability of dying q(x,n)",
#                  dx = "Number of deaths d(x,n)",
#                  lx = "Number of survivors l(x)",
#                  ex = "Expectation of life e(x)")  %>% 
#       select(country_code,Period,x, mx, qx, dx, lx, ex)
# tfr_paises_LA <- tfr %>% select(-last.observed) %>% gather(Period,tfr,-country_code,-name)
# asfr_paises_LA <- percentASFR %>% gather(Period, asfr, -name, -country_code, -age) %>% 
#   mutate(x = as.numeric(substr(age,1,2)), asfr=asfr/100) %>% 
#   left_join(tfr_paises_LA, by=c("country_code","name","Period")) %>% 
#   mutate(asfr = asfr * tfr) %>% 
#   select(-age)
# 
# # smooth
# # smooth
# asfr_paises_LA_s <- data.frame(country_code=, Period = " ", x = 0, asfr = 0, stringsAsFactors = F)
# lx_paises_LA_s <- data.frame(country_code=0, Period = " ", x = 0, lx = 0, stringsAsFactors = F)
# countries = unique(asfr_paises_LA$country_code)
# for(c in 1:length(countries)){
#   for(p in unique(asfr_paises_LA$Period)){
#     # c = unique(asfr_paises_LA$country_code)[1]
#     # p = "1950-1955"
#     cp = asfr_paises_LA %>% filter(country_code==countries[c], Period==p) %>% select(asfr)
#     if(any(is.na(cp$asfr))){next}
#     out <- QOSplit(c(0,cp$asfr), L=seq(10,45,5), AgeInt=rep(5,8))
#     asfr_paises_LA_s <- rbind(asfr_paises_LA_s, data.frame(country_code = countries[c], Period = p, x = out$Age, asfr = out$ASFR))
#     cp = lt_paises_LA %>% filter(country_code==countries[c], Period==p) %>% select(x,lx)
#     out <- lx_fun(cp, age=0:100)
#     lx_paises_LA_s <- rbind(lx_paises_LA_s, data.frame(country_code=countries[c], Period = p, x = 0:100, lx = out))
#     # plot(seq(10,45,5),c(0,cp$asfr),t="s");lines(out,col=2)
#   }
#   print(c)
# }
# 
# # final data
# LA <- paises_LA %>% 
#       left_join(lx_paises_LA_s,by=c("country_code")) %>% 
#       left_join(asfr_paises_LA_s, by=c("country_code","Period","x")) %>%
#       left_join(tfr_paises_LA,by=c("country_code","Period"))
# 
# # grafs
# tfr_q0 = ggplot(LA %>% filter(x==1), aes(x=tfr, y=(1-lx/100000))) + geom_smooth(se = T, color=1) + geom_point( aes(color=Period)) + theme_classic()




LA %>% filter(name=="Guatemala") %>% ggplot() +
  geom_line(aes(x=x, y=asfr/tfr*100, color=Period))
LA %>% filter(name=="Guatemala") %>% ggplot() +
  geom_line(aes(x=x, y=lx, color=Period))



LA %>% filter(name=="Argentina") %>% ggplot() +
  geom_line(aes(x=x, y=asfr/tfr*100, color=Period))
CS_outs %>% filter(name=="Argentina")


