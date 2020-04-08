
library(HMDHFDplus)
user = "ivanwilliams1985@gmail.com"
pwHMD ="volveroman"
pwHFD = 52962

Austria_F <- readHFDweb("AUT","asfrRR", user, pwHFD)
Austria_M <- readHMDweb("AUT","fltper_1x1", user, pwHMD)
alpha = 15
a = 40
range = a - alpha
year = 1951
Austria_F2017 <- subset(Austria_F, Year==year & Age %in% alpha:a)
Austria_M2017 <- subset(Austria_M, Year==year & Age %in% 0:range)
plot(Austria_F2017$Age, Austria_F2017$ASFR, t="l")
par(new = T)
plot(Austria_M2017$Age, rev(Austria_M2017$lx)/100000, ylim = c(.9,1), 
     col=2, axes=F, xlab=NA, ylab=NA, t="l")



# theorical cases

# big hump case

#input
par(mfrow=c(1,1),mar = c(5,5,2,5))
f_young = data.frame(x=c(15,20,25,30,35,40), f=c(0,.3,.3,.05,.01,.01))
f_old = data.frame(  x=c(15,20,25,30,35,40), f=c(0,.04,.06,.25,.1,.01))
f_young_p = predict(loess(f~x, data = f_young), 15:40); sum(f_young_p)
f_old_p = predict(loess(f~x, data = f_old), 15:40); sum(f_old_p)
f_young_mean = sum(f_young_p*15:40/sum(f_young_p)); f_young_mean
f_old_mean = sum(f_old_p*15:40/sum(f_old_p)); f_old_mean
library(latex2exp)
m_si_hump = data.frame(x= c(0,1,5,10,15,20,25), l=c(1,.98,.97,.97,.96,.9,.85))
m_si_hump_sp = splinefun(m_si_hump, method = "monoH.FC")
#graph
plot(15:40, f_young_p, 
     t="l",xlab="x (mother age)", ylab= TeX(sprintf("$m_x$"))) 
lines(15:40, f_old_p, lty=2) 
abline(v=f_young_mean, lty=2, col="grey")
abline(v=f_old_mean, lty=2, col="grey")
par(new = T)
plot(0:25, rev(m_si_hump_sp(0:25)), axes=F, xlab=NA, ylab=NA, t="l", col=2)
axis(side = 4)
mtext(side = 4, line = 3, TeX(sprintf("$l_{40-x}$")))


H_young = -sum(f_young_p*rev(m_si_hump_sp(0:25))*log(rev(m_si_hump_sp(0:25))))/
  sum(f_young_p*rev(m_si_hump_sp(0:25)))

H_old = -sum(f_old_p*rev(m_si_hump_sp(0:25))*log(rev(m_si_hump_sp(0:25))))/
  sum(f_old_p*rev(m_si_hump_sp(0:25)))
