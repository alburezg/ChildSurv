library(tidyverse)

get_kins <- function(U_list, F_list, w=100, pop_mother){
  # x ego´s age at time t + n
  # U_list survive matrix
  # F_list fertility matrix
  
  asfrt = LA %>% filter(name=="Argentina", Period=="1950-1955") %>%
              mutate(asfr = ifelse(is.na(asfr),0,asfr)) %>% pull(asfr)
  Pt = LA %>% filter(name=="Argentina", Period=="1950-1955") %>% 
              mutate(px = 1-qx) %>% pull(px)
  ages = length(Pt)
  Ut = Ft = matrix(0, nrow=ages, ncol=ages)
  Ut[row(Ut)-1 == col(Ut)] <- Pt[-101]
  Ut[ages,ages]=Pt[101]
  Ft[1,] = asfrt
  
  kins = data.frame(x=0, a=a, b=b, c=c, n=n)
  
  a = b = c = n = rep(0,ages)
  # vectors of kins by age
  d = runif(ages)
  d[c(1:15,50:ages)] = 0
  d = d/sum(d)
  g = d
  m = d
  
  e = rep(0,ages)
  
  age_ego = 70
  for(x in 2:age_ego){
    
    # daughters
    e[x] = 1 
    a = Ut %*% a + Ft %*% e
    
    # granddaughters
    b = Ut %*% b + Ft %*% a
    
    # mother
    d = Ut %*% d
    
    # younger sisters
    n = Ut %*% n + Ft %*% d
    
    # grandmother
    g = Ut %*% g
    
    # older sisters
    m = Ut %*% m
    
    e = rep(0,ages)
    }
    
  plot(a,ylim=c(0,.1))
  lines(b)  
  plot(n)
  
    # get results
    kins = rbind(kins, data.frame(x=x, a=a, b=b, c=c, d=d, g=g, m=m, n=n))
  }
  
  
  






F_m = matrix(c(.3,.3,.3,
               0,0,0,
               0,0,0),nrow = 3,byrow = T)


ex = c(0,1,0)

F_m%*%t(ex)
