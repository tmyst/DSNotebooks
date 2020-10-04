library(tidyverse)
library(rms)
library(survival)
library(Hmisc)
library(final)
library(finalfit)

n <- 200
p <- 20
set.seed(6)

xx <- matrix(rnorm(n * p), nrow = n, ncol = p)
y <- runif(n)

units(y) <- "Year"

colnames(xx) <- str_c("x", seq(1,20,1))
e <- c(rep(0, n/2), rep(1, n/2))
d <- data.frame(cbind(xx, y, e))

f <- cph(Surv(y, e) ~ xx, x = TRUE, y = TRUE, time.inc = 0.5, surv = TRUE)
g <- coxph(Surv(y, e) ~ xx)

# ブートストラプサンプルを使ってクロスバリデーションを行う
# overfitによるバイアスを加味した「予測vs観測」の推定を得るため。
# 各インターバルに対し部分的予測を用いるか、
# またはノンパラメトリックなスムージングを行う
# 予測値は定点における予測生存率、観測はKM法による確率
# 何もcmethodに入れなければ以下のようにpolspline を利用してスムージングした
# ものとなる。KMを指定した場合は、
cal <- calibrate(f, u=.5, B=200, cmethod = "hare")
cal2 <- calibrate(f, u=.5, B=300)

attributes(cal)
attributes(cal2)

plot(cal, ylim=c(.4, 1), subtitles=T)
plot(cal2, ylim=c(.4, 1), subtitles=T)

calkm <- calibrate(f, u=.5, m=50, cmethod= "KM", B=200)
# 青い点はKMの修正バージョン、黒はKMそのまま
plot(calkm)


# Proportional Hazard model -----------------------------------------------

colnames(xx) <- str_c("x", 1:ncol(xx))
vars <- colnames(xx)
surv <- "Surv(y,e)"
fit <- finalfit(dependent = surv, explanatory = vars, data = data.frame(y,e,xx))
fit %>% htmlTable::htmlTable() # to save the result, use kableExtra::save_kable(file = "xxx.html")

# plot HR
data.frame(y,e,xx) %>% hr_plot(explanatory = vars, dependent = surv)


# estimation and validation -----------------------------------------------

library(rms)
library(survminer)
library(ggsci)
library(unikn)

sobj <- with(veteran,Surv(time,status))

line_c <- pal_simpsons()(2)[1]
ci_c <- pal_simpsons()(2)[2]

km <- survfit(sobj ~ 1, data = veteran)
plot(km, ylab = "survival", xlab = "time(/days)", lty = c(1,2,2), lw = c(2,1,1), col= c(line_c, ci_c, ci_c))

# not so much related
km2 <- survfit(sobj ~ trt, data = veteran)
ggsurvplot(km2, conf.int = T)

# important
km3 <- survfit(sobj ~ veteran$celltype, data = veteran)
ggsurvplot(km3, conf.int = T)

# important
km4 <- survfit(sobj ~ veteran$karno, data = veteran)
ggsurvplot(km4, conf.int = T)

# not so much related
km5 <- survfit(sobj ~ veteran$age, data = veteran)
ggsurvplot(km5, conf.int = T)


# cox ph ------------------------------------------------------------------
# use cut of age
df <- veteran %>% mutate(age_c = cut(age, breaks = seq(0, 100, 5)) %>% droplevels())
model <- coxph(sobj ~ strata(trt) + karno + age_c + celltype, data = df) %>% step(direction = "both")
summary(model)
zph  <- cox.zph(model)
modelf <- finalfit(.data = df, explanatory = c("celltype", "karno", "age_c", "strata(trt)"), dependent = "Surv(time, status)")
hr_plot(.data = df, explanatory = c("celltype", "karno", "age_c", "strata(trt)"), dependent = "Surv(time, status)")
modelf %>% htmlTable::htmlTable()

# use raw age
model2 <- coxph(sobj ~ strata(trt) + karno + age + celltype, data = df) %>% step(direction = "both")
summary(model2)
zph2  <- cox.zph(model2)
modelf2 <- finalfit(.data = df, explanatory = c("celltype", "karno", "age", "strata(trt)"), dependent = "Surv(time, status)")
hr_plot(.data = df, explanatory = c("celltype", "karno", "age", "strata(trt)"), dependent = "Surv(time, status)")
modelf2 %>% htmlTable::htmlTable()


# from Harrell's book -----------------------------------------------------

set.seed(109090)
train_ind <- sample(x = 1:nrow(veteran), size = floor(nrow(veteran)*0.6), replace=F) %>% sort
is_train <- 1:nrow(veteran) %in% train_ind
sobj_train <- with(veteran[is_train,], Surv(time, status))

cox <- cph(sobj~celltype + karno, 
           data=veteran, 
           x=TRUE, 
           y=TRUE, 
           surv=TRUE,
           time.inc=5*365)

cox_train <- cph(sobj_train~celltype + karno, 
           data=veteran[train_ind,], 
           x=TRUE, 
           y=TRUE, 
           surv=TRUE,
           time.inc=5*365)

# External data
test_dat <- data.frame(trt=replicate(500,NA), 
                       celltype=replicate(500,NA), 
                       time=replicate(500,NA), 
                       status=replicate(500,NA), 
                       karno=replicate(500,NA), 
                       diagtime=replicate(500,NA), 
                       age=replicate(500,NA), 
                       prior=replicate(500,NA))
for(i in seq(8)){
  test_dat[,i]=sample(veteran[,i],500,replace=T)
}

#
test_veteran <- veteran[!is_train,]

#
sobj_test_dat     <- with(test_dat,Surv(time,status))
sobj_test_veteran <- with(test_veteran, Surv(time, status))

# Create estimation
estimates_test_dat     <- survest(cox, newdata = test_dat, times=5*365)$surv
estimates_test_veteran <- survest(cox, newdata = veteran[!is_train,], times = 5*777)$surv
estimates_train_veteran <- survest(cox, newdata = veteran[is_train,], times = 5*365)$surv

# Concordance
rcorr.cens(x = estimates_test_dat, S = sobj_test_dat)
rcorr.cens(x = estimates_test_veteran, S = sobj_test_veteran)
rcorr.cens(x = estimates_train_veteran, S = sobj_train)

#
latex(validate(f, B=200), digits=3, file= '', caption= '', table.env=TRUE, label= 'tab:cox-val-random') 

# P488 make dataset
n <- 2000

set.seed (3)
age <- 50 + 12 * rnorm(n)
label(age) <- 'Age'

factor(c(T, F), 0:1, c("huhuhu", "ipoip"))
sex  <-  factor(1 + (runif(n) <= .4), 1:2, c( 'Male', 'Female'))
cens <-  15 * runif(n)
h    <-  0.02 * exp(0.04 * (age - 50) + 0.8 * (sex == 'Female'))
ft   <-  -log(runif(n)) / h
g1 <- cbind(age, sex, h) %>% data.frame() %>% ggplot() + 
  geom_point(aes(x = age, y = h, color = factor(sex))) +
  theme_light(base_size = 14) + scale_color_aaas() + 
  theme(legend.title = element_blank()) 

g2 <- cbind(age, sex, ft) %>% data.frame() %>% ggplot() + 
  geom_point(aes(x = age, y = ft, color = factor(sex))) +
  theme_light(base_size = 14) + scale_color_aaas() + 
  theme(legend.title = element_blank()) 

e  <-  ifelse(ft <= cens, 1, 0)
print (table (e))

ft <- pmin(ft, cens)
units(ft) <- 'Year'

Srv <- Surv(ft, e)

age.dec <- cut2(age, g=10, levels.mean=TRUE)
label(age.dec) <- 'Age'

dd <- datadist(age, sex, age.dec); 
options(datadist= 'dd') 
f.np <- cph(Srv ~ strat(age.dec) + strat(sex), surv=TRUE)
km <- survfit(Srv~age.dec, data = data.frame(cbind(ft, e, age.dec)))
ggsurvplot(km, conf.int = T, palette = pal_simpsons()(10))

# surv=TRUE speeds up computations, and confidence limits when 
# there are no covariables are still accurate.
p <- Predict(f.np, age.dec, sex, time=3, loglog=TRUE)

# Treat age.dec as a numeric variable (means within deciles) 
p$age.dec <- as.numeric(as.character(p$age.dec))
ggplot(p, ylim=c(-5, -.5)) + theme_bw()


f.noia <- cph(Srv ~ rcs(age,4) + strat(sex), x=TRUE, y=TRUE)

# Get accurate C.L. for any age by specifying x=TRUE y=TRUE
# Note : for evaluating shape of regression , we would not
# ordinarily bother to get 3-year survival probabilities -
# would just use X * beta
# We do so here to use same scale as nonparametric estimates
library(texPreview)
w <- latex(f.noia , inline =TRUE , digits =3)
latex(anova(f.noia), table.env=FALSE, file='')

ploglog <-  Predict(f.noia, age, sex, time=3, loglog=TRUE)
p <-  Predict(f.noia, age, sex, time=3, loglog=FALSE)

ggplot(ploglog, ylim=c(-5, -.5))
ggplot(p)

x1 <- seq(0, 100, 0.5)
yhat <- x1
lower <- x1-5
upper <- x1+5


## 
f.ia <-  cph(Srv ~ rcs(age,4) * strat(sex), x=TRUE, y=TRUE, surv=TRUE)
w <- latex(f.ia , inline = TRUE, digits = 3)
latex(anova(f.ia), table.env=FALSE, file= '')
ploglog <-  Predict(f.ia, age, sex, time = 3, loglog = T)
p <-  Predict(f.ia, age, sex, time = 3, loglog = F)
ggplot(p)
ggplot(ploglog, ylim = c(-5, -0.5))



# kihon -------------------------------------------------------------------
x1 <- seq(0, 100, 0.1)
srv_a <- (1/10)*exp(-x1/10)
plot(srv_a ~ x1)

ht <- -log(srv_a)
plot(ht~x1)

lam <- (diff(srv_a))/(srv_a[2:length(srv_a)])
plot(lam~x1[2:length(x1)], ylim = c(-1, 1))


# Spline ------------------------------------------------------------------
rcspline.eval(seq(0, 1, 0.01), knots = seq(0.05, 0.95, length=5), inclx = T) %>% matplot(.[,1], ., lty=1, type="l")
x <- Hmisc::rcspline.eval(seq(0,1,.01), knots = seq(.05,.95,length=5), inclx=T)
xm <- x
graphics::matplot(x[,1], x, lty=1, type="l")

xm[xm > .0106] <- NA
graphics::matplot(x[,1], xm, type="l", ylim=c(0,.01), xlab=expression(X), ylab= '', lty=1) 
graphics::matplot(x[,1], x,  type="l", xlab=expression(X), ylab= '', lty=1)

x <- seq(0, 1, length =300)
for(nk in 3:6){
  set.seed(nk)
  knots <- seq(0.05, 0.95, length=nk)
  xx    <- rcspline.eval(x, knots=knots, inclx=T) 
  for(i in 1:(nk - 1)) xx[,i] <- (xx[,i] - min(xx[,i])) / (max(xx[,i]) - min(xx[,i]))
  for(i in 1:20){ 
    beta <- 2*runif(nk-1) - 1
    xbeta <- xx %*% beta + 2*runif(1) - 1 
    xbete <- (xbeta - min(xbeta)) / (max(xbeta) - min(xbeta))
    
    if(i == 1){
      plot(x, xbeta, type="l", lty=1, xlab=expression(X), ylab='', bty="l", ylim=c(0,1)) 
      title(sub=paste(nk,"knots"), adj=0, cex=.75)
      for(j in 1:nk) arrows(knots[j], 0.04, knots[j], -0.03, angle=20, length=0.07, lwd=1.5)
    }
    else lines(x, xbeta, col=i)
  }
}





