library(tidyverse)
library(JM)
library(survival)
library(survminer)
library(unikn)
library(scico)
library(ggsci)
library(wesanderson)
library(colorspace)
library(ggpubr)
library(ggrepel)
library(patchwork)
#library(ochRe)
library(naniar)
library(eechidna)
library(oceanmap)
library(marmap)
library(pals)
library(yarrr)
library(viridis)
library(RColorBrewer)

theme_set(new = theme_light(base_size = 14))

aids %>% ggplot() + geom_line(aes(x = obstime , y = CD4, group = patient, color = gender)) + 
  facet_grid(.~drug) + 
  scale_color_nejm()

surv_data <- with(data = aids.id, Surv(time = Time, event = death))

# random intercept model for CD4 ------------------------------------------

lme_fit <- lme(CD4~obstime, random = ~1 | patient, data = aids)
summary(lme_fit)

#Vi =σb21ni1⊤ni +σ2Ini,
margCov1 <- getVarCov(lme_fit, individuals = 10, type = "marginal")
margCov1
cor1 <- cov2cor(margCov1[[1]])
cor1

#ただし、どの時間帯のペアをとっても相関が不変というのは不自然

# random intercept & slope model ------------------------------------------
lme_slope <- lme(CD4~obstime, random = ~ obstime | patient, data = aids)
summary(lme_slope)
margCov2 <- getVarCov(lme_slope, individuals = 10, type = "marginal")
margCov2
cor2 <- cov2cor((margCov2[[1]]))
cor2

