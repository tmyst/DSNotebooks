#rethinking
devtools::install_github("rmcelreath/rethinking")

#ggplot2 loading error in library(retihinking) was fixed by installing "scales" again
library(rethinking)
library(tidyverse)

data(WaffleDivorce) 
d <- WaffleDivorce

theme_set(theme_bw(base_size = 16))
# points
d %>% ggplot() + 
  geom_point(aes(x = MedianAgeMarriage, y = Divorce))+
  geom_errorbar(aes(x = MedianAgeMarriage, 
                    ymin = Divorce-Divorce.SE, 
                    ymax = Divorce+Divorce.SE), 
                width = .1,linetype = 1) +
  xlab("Median age marriage") +
  ylab("Divorce rate") 

d %>% as_tibble()

dlist <- list(
  D_obs = standardize( d$Divorce ),
  D_sd =  d$Divorce.SE / sd( d$Divorce ), 
  M = standardize( d$Marriage ),
  A = standardize( d$MedianAgeMarriage ), 
  N = nrow(d)
)

m15.1 <- ulam( alist(
  D_obs ~ dnorm( D_true , D_sd ), 
  vector[N]:D_true ~ dnorm( mu , sigma ), 
  mu <- a + bA*A + bM*M,
  a ~ dnorm(0,0.2),
  bA ~ dnorm(0,0.5),
  bM ~ dnorm(0,0.5),
  sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 )

precis( m15.1 , depth=2 )
library(vcd)

#install.packages("Rcpp", repos="https://rcppcore.github.io/drat")
