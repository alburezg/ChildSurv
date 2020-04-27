
source("QOSplit.R")

# main_funtion
# age survivor spline

lx_fun <- function(data, age){
  splinefun(data$x, data$lx, method = "monoH.FC")(age)
}

# function to calculate all the measures

CS_results <- function(country, period, a){
  
  # filter country
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

# Get wpp locations
data(UNlocations)

# countries of LA
paises_LA_seleccion <- c(192,214,332,388,630,188,222,320,340,484,558,591,32,68,76,152,170,218,600,604,858,862)

# get fertility data
data(tfr)
data(percentASFR)
tfr_paises_LA <- tfr %>% select(-last.observed) %>% gather(Period,tfr,-country_code,-name) %>% filter(country_code %in% paises_LA_seleccion)

asfr_paises_LA <- percentASFR %>% gather(Period, asfr, -name, -country_code, -age) %>% 
  mutate(x = as.numeric(substr(age,1,2)), asfr=asfr/100) %>% 
  left_join(tfr_paises_LA, by=c("country_code","name","Period")) %>% 
  mutate(asfr = asfr/5 * tfr) %>% 
  select(-age) %>% 
  filter(country_code %in% paises_LA_seleccion)

# lt must downloaded manually
lt_paises_LA <- readxl::read_xlsx("WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", sheet = 1, range = "A17:T77381") %>% 
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

# smooth both
asfr_paises_LA_s <- data.frame(country_code=0, Period = " ", x = 0, asfr = 0, stringsAsFactors = F)
lx_paises_LA_s <- data.frame(country_code=0, Period = " ", x = 0, lx = 0, stringsAsFactors = F)
countries = unique(asfr_paises_LA$country_code)
for(c in 1:length(countries)){
  for(p in unique(asfr_paises_LA$Period)){
    cp = asfr_paises_LA %>% filter(country_code==countries[c], Period==p) %>% select(asfr)
    if(any(is.na(cp$asfr))){next}
    out <- tryCatch(QOSplit(c(0,cp$asfr), L=seq(10,45,5), AgeInt=rep(5,8)),
                    error = data.frame(Age=10:49, ASFR=rep(NA,40)))
    asfr_paises_LA_s <- rbind(asfr_paises_LA_s, data.frame(country_code = countries[c], Period = p, x = 
                                                             out$Age, asfr = out$ASFR))
    cp <- lt_paises_LA %>% filter(country_code==countries[c], Period==p) %>% select(x,lx)
    out <- lx_fun(cp, age=0:100)
    lx_paises_LA_s <- rbind(lx_paises_LA_s, data.frame(country_code=countries[c], Period = p, x = 0:100, lx = out))
  }
}

# join everything
# assume
im = 100/204

LA <- lx_paises_LA_s %>% slice(-1) %>% mutate(country_code=as.integer(country_code)) %>% 
  left_join(asfr_paises_LA_s %>% slice(-1), by=c("country_code","Period","x")) %>%
  left_join(tfr_paises_LA, by=c("country_code","Period")) %>% 
  left_join(lt_paises_LA %>% 
              select(country_code, Period, x, Lx_n, ex),
            by = c("country_code","Period","x")) %>% 
  mutate(q = (1-lx/100000),
         name = as.character(name),
         tfr = tfr*im,
         asfr = asfr*im)
LA$name[LA$name == "Bolivia (Plurinational State of)"] = "Bolivia"


# get_results
CS_outs <- as.data.frame(expand.grid(name = unique(LA$name),
                                     Age = c(30,40,50),
                                     Period = unique(LA$Period)))
for(i in 1:nrow(CS_outs)){
  CS_outs[i,4:12] = CS_results(country = CS_outs[i,1], 
                               period = CS_outs[i,3], 
                               a = CS_outs[i,2])}
CS_outs <- CS_outs %>% mutate(Year = as.integer(substr(Period, 1, 4)) + 2.5)

