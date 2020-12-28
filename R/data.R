get_data <- function(source = c("WPP", "HMD_HFD"),
                     country = "Sweden",
                     years = 1950:2020){
  
  
  lx_fun <- function(data, age){
    splinefun(data$x, data$lx, method = "monoH.FC")(age)
  }
  
  ################# WPP
  library(wpp2019)
  # get fertility data
  data(tfr)
  data(tfrprojMed)
  data(percentASFR)
  
  # tfr
  tfr_country <- tfr %>% bind_cols(tfrprojMed %>% select(-1,-2)) %>% 
                    select(-last.observed) %>% 
                    gather(Period, tfr, -country_code, -name) %>% 
                    filter(name %in% country)
  # asfr join tfr
  asfr_country <- percentASFR %>% gather(Period, asfr, -name, -country_code, -age) %>% 
                    mutate(x = as.numeric(substr(age,1,2)), asfr=asfr/100) %>% 
                    inner_join(tfr_country, by=c("country_code","name","Period")) %>% 
                    mutate(asfr = asfr/5 * tfr) %>% 
                    select(-age)
  
  # lt data must be downloaded manually (maybe too heavy for CRAN)
  lt_paises_LA <- readxl::read_xlsx("C:/Proyectos/ChildSurv/Data/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", 
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
  
  
} 