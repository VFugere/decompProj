### Code for experimental project on drivers of decomposition in Kibale streams ###
### Vincent Fugère, 2019

rm(list=ls())
par(family='sans')
# options(scipen=999)

library(tidyverse)
library(RColorBrewer)
library(glmmTMB)
library(magrittr)
library(plotrix)

cols <- brewer.pal(3, 'Dark2')

#load and format data

data <- read_csv('~/Google Drive/Recherche/PhD/manuscripts/caterpillar/decompProj/leafbags.csv') %>% filter(!is.na(rep.nb))

dmg <-  read_csv('~/Google Drive/Recherche/PhD/manuscripts/caterpillar/decompProj/caterpillar.csv') %>%
  select(leaf.nb,damage.area)

data <- inner_join(data, dmg, by = 'leaf.nb') #loosing one replicate - the leaf fragment with a missing picture
rm(dmg)

data$prop <- data$prct.mass.rem/100
data$prop.decomp <- 1-data$prop

#transformation based on: https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
data$rv <- 1-(data$prop*(nrow(data)-1)+0.5) / nrow(data)
#1- because we want prop decomposed instead of prop left

#data$dmg <- scale(data$damage.area)
data$dmg <- rescale(data$damage.area, c(0,1))
data$lu <- as.factor(data$land.use)
data$lu <- relevel(data$lu, ref = 'forest')
#data$lu <- scale(as.numeric(data$lu))
data$ha <- as.factor(data$leaf.deployment)
data$ha <- relevel(data$ha, ref = 'home')
#data$ha <- scale(as.numeric(data$ha))
data$wks <- scale(data$weeks)
data$wtr <- as.factor(data$watershed)

## fine mesh bags

fm <- filter(data, mesh.type == 'fine')

m0 <- glmmTMB(rv ~ 1 + wks + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))
m1 <- glmmTMB(rv ~ 1 + wks + wks:lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))
m2 <- glmmTMB(rv ~ 1 + wks + wks:ha + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))
m3 <- glmmTMB(rv ~ 1 + wks + wks:dmg + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))
m4 <- glmmTMB(rv ~ 1 + wks + wks:dmg + wks:ha + wks:lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))
m5 <- glmmTMB(rv ~ 1 + wks + wks:ha*lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))
m6 <- glmmTMB(rv ~ 1 + wks + wks:ha*dmg + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))
m7 <- glmmTMB(rv ~ 1 + wks + wks:dmg*lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))

models.fm <- list(m0,m1,m2,m3,m4,m5,m6,m7)
map(models.fm, summary)

f0 <- formula(rv ~ 1 + wks + watershed + dmg + ha + lu)
f1 <- update(f0, . ~ wks * .)
f2 <- update(f1, . ~ . + (wks|site) + (wks|leaf) + (wks|leaf.origin))
m4 <- glmmTMB(f2, fm, family=beta_family(link = "logit"))

m4 <- glmmTMB(rv ~ 0 + weeks + weeks:wtr + weeks:dmg + weeks:ha + weeks:lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), fm, family=beta_family(link = "logit"))

summary(m4)

#Fig. 2a) caterpillar plot for FM bags

#getting coefs
coefs <- rownames_to_column(as.data.frame(summary(m4)$coefficients[1])) %>%
  rename(coef = rowname, value = cond.Estimate, se = cond.Std..Error) %>% select(coef:se) %>%
  filter(coef != '(Intercept)', coef != 'weeks')
coefs %<>% mutate('lwr' = value-1.96*se, 'upr' = value+1.96*se)

xmin <- min(coefs$lwr)
xmax <- max(coefs$upr)

pdf('~/Desktop/Fig2.pdf',width=7,height = 3.8,pointsize = 8)

layout(rbind(c(1,1,2,3,4),c(5,5,6,7,8)))

par(mar=c(4,7,1,1),cex=1)

plot(0,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(xmin,xmax),ylim=c(0.5,3.5))
axis(2,cex.axis=1,lwd=0,lwd.ticks=0,at=1:3,labels = rev(c('farm vs. forest','away vs. home','leaf damage')),las=1)
abline(v=0,lty=2,lwd=1)
title(xlab='parameter estimate')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
for(w in 1:3){
  pt <- coefs[w,'value']
  lwr <- coefs[w,'lwr']
  upr <- coefs[w,'upr']
  if(0 < upr & 0 > lwr){
    alph <- 0.3
  } else {
    alph <- 1
  }
  segments(x0=lwr,x1=upr,y0=w,y1=w,col=alpha(rev(cols)[w],alph),lwd=3)
  points(x=pt,y=w,pch=16,col=alpha(rev(cols)[w],alph),cex=2)
}
legend('topright',bty='n',fill=0,border=0,legend=expression(italic('fine mesh')))
#legend('topright',bty='n',legend = '(a)', cex =1.2)

par(mar=c(4,4,1,1))

plot(prop.decomp~weeks,fm,type='n',yaxt='n',xaxt='n',xlim=c(1,4),ylim=range(fm$prop.decomp),ann=F,bty='l')
title(xlab='days in stream')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:4,labels = seq(7,28,by=7))
title(ylab='proportion decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5))
means <- aggregate(prop.decomp~weeks*land.use, fm, FUN = 'mean')
se <- aggregate(prop.decomp~weeks*land.use, fm, FUN = 'sd')
means$se <- 1.96*se$prop.decomp/sqrt(nrow(se))
means$land.use <- as.factor(means$land.use)
coltmp <- c(cols[1],1)
for(i in 1:2){
  data.tmp <- subset(means, land.use == levels(means$land.use)[i])
  plot.data <- as.data.frame(approx(x = data.tmp$weeks, y = data.tmp$prop.decomp, n = 1000))
  plot.data$se <- as.data.frame(approx(x = data.tmp$weeks, y = data.tmp$se, n = 1000))[,2]
  polygon(x=c(plot.data$x,rev(plot.data$x)), y=c(plot.data$y+plot.data$se,rev(plot.data$y-plot.data$se)), border=NA, col=alpha(coltmp[i],0.2))
  lines(y ~ x, plot.data, lwd=2, col = alpha(coltmp[i],0.8))
}
rm(data.tmp,i)
legend('topleft',bty='n',legend = c('forest','farm'), lty=1, col=rev(coltmp),lwd=2,seg.len = 0.8)
#legend('topright',bty='n',legend = '(b)', cex =1.2)

plot(prop.decomp~weeks,fm,type='n',yaxt='n',xaxt='n',xlim=c(1,4),ylim=range(fm$prop.decomp),ann=F,bty='l')
title(xlab='days in stream')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:4,labels = seq(7,28,by=7))
title(ylab='proportion decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5))
means <- aggregate(prop.decomp~weeks*leaf.deployment, fm, FUN = 'mean')
se <- aggregate(prop.decomp~weeks*leaf.deployment, fm, FUN = 'sd')
means$se <- 1.96*se$prop.decomp/sqrt(nrow(se))
means$leaf.deployment <- as.factor(means$leaf.deployment)
coltmp <- c(cols[2],1)
for(i in 1:2){
  data.tmp <- subset(means, leaf.deployment == levels(means$leaf.deployment)[i])
  plot.data <- as.data.frame(approx(x = data.tmp$weeks, y = data.tmp$prop.decomp, n = 1000))
  plot.data$se <- as.data.frame(approx(x = data.tmp$weeks, y = data.tmp$se, n = 1000))[,2]
  polygon(x=c(plot.data$x,rev(plot.data$x)), y=c(plot.data$y+plot.data$se,rev(plot.data$y-plot.data$se)), border=NA, col=alpha(coltmp[i],0.2))
  lines(y ~ x, plot.data, lwd=2, col = alpha(coltmp[i],0.8))
}
rm(data.tmp,i)
legend('topleft',bty='n',legend = c('home','away'), lty=1, col=rev(coltmp),lwd=2,seg.len = 0.8)
#legend('topright',bty='n',legend = '(c)', cex =1.2)

tp4 <- subset(fm, weeks == 4)
tp4$damage.area <- tp4$damage.area*100
plot(prop.decomp~damage.area,tp4,type='n',yaxt='n',xaxt='n',xlim=range(tp4$damage.area),ylim=range(tp4$prop.decomp),ann=F,bty='l')
title(xlab='leaf damage (%)')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='proportion decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5))
#points(prop.decomp~damage.area,subset(tp4, weeks == 3),col=alpha(1,0.5),pch=16)
points(prop.decomp~damage.area,tp4,col=alpha(cols[3],0.5),pch=16)
#legend('topright',bty='n',legend = c('week 3','week 4'), pch=15, col=c(1,cols[3]))
legend('topright',bty='n',legend=expression(italic('day 28')))
#legend('topright',bty='n',legend = '(d)', cex =1.2)

#### coarse mesh bags ####

cm <- filter(data, mesh.type == 'coarse')

m0 <- glmmTMB(rv ~ 1 + wks + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), cm, family=beta_family(link = "logit"))
m1 <- glmmTMB(rv ~ 1 + wks + wks:lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), cm, family=beta_family(link = "logit"))
m2 <- glmmTMB(rv ~ 1 + wks + wks:ha + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), cm, family=beta_family(link = "logit"))
m3 <- glmmTMB(rv ~ 1 + wks + wks:dmg + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), cm, family=beta_family(link = "logit"))
m4 <- glmmTMB(rv ~ 1 + wks + wks:dmg + wks:ha + wks:lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), cm, family=beta_family(link = "logit"))
m5 <- glmmTMB(rv ~ 1 + wks + wks:ha*lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), cm, family=beta_family(link = "logit"))
m6 <- glmmTMB(rv ~ 1 + wks + wks:ha*dmg + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), cm, family=beta_family(link = "logit"))
m7 <- glmmTMB(rv ~ 1 + wks + wks:dmg*lu + (wks-1|site) + (wks-1|leaf) + (wks-1|leaf.origin), cm, family=beta_family(link = "logit"))

models.cm <- list(m0,m1,m2,m3,m4,m5,m6,m7)
map(models.cm, summary)
summary(m4)

#Fig. 2a) caterpillar plot for cm bags

#getting coefs
coefs <- rownames_to_column(as.data.frame(summary(m4)$coefficients[1])) %>%
  rename(coef = rowname, value = cond.Estimate, se = cond.Std..Error) %>% select(coef:se) %>%
  filter(coef != '(Intercept)', coef != 'weeks')
coefs %<>% mutate('lwr' = value-1.96*se, 'upr' = value+1.96*se)

xmin <- min(coefs$lwr)
xmax <- max(coefs$upr)

par(mar=c(4,7,1,1))
plot(0,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(xmin,xmax),ylim=c(0.5,3.5))
axis(2,cex.axis=1,lwd=0,lwd.ticks=0,at=1:3,labels = rev(c('farm vs. forest','away vs. home','leaf damage')),las=1)
abline(v=0,lty=2,lwd=1)
title(xlab='parameter estimate')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
for(w in 1:3){
  pt <- coefs[w,'value']
  lwr <- coefs[w,'lwr']
  upr <- coefs[w,'upr']
  if(0 < upr & 0 > lwr){
    alph <- 0.3
  } else {
    alph <- 1
  }
  segments(x0=lwr,x1=upr,y0=w,y1=w,col=alpha(rev(cols)[w],alph),lwd=3)
  points(x=pt,y=w,pch=16,col=alpha(rev(cols)[w],alph),cex=2)
}
legend('topright',bty='n',fill=0,border=0,legend=expression(italic('coarse mesh')))
#legend('topright',bty='n',legend = '(b)', cex =1.2)

par(mar=c(4,4,1,1))

plot(prop.decomp~weeks,cm,type='n',yaxt='n',xaxt='n',xlim=c(1,4),ylim=range(cm$prop.decomp),ann=F,bty='l')
title(xlab='days in stream')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:4,labels = seq(7,28,by=7))
title(ylab='proportion decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5))
means <- aggregate(prop.decomp~weeks*land.use, cm, FUN = 'mean')
se <- aggregate(prop.decomp~weeks*land.use, cm, FUN = 'sd')
means$se <- 1.96*se$prop.decomp/sqrt(nrow(se))
means$land.use <- as.factor(means$land.use)
coltmp <- c(cols[1],1)
for(i in 1:2){
  data.tmp <- subset(means, land.use == levels(means$land.use)[i])
  plot.data <- as.data.frame(approx(x = data.tmp$weeks, y = data.tmp$prop.decomp, n = 1000))
  plot.data$se <- as.data.frame(approx(x = data.tmp$weeks, y = data.tmp$se, n = 1000))[,2]
  polygon(x=c(plot.data$x,rev(plot.data$x)), y=c(plot.data$y+plot.data$se,rev(plot.data$y-plot.data$se)), border=NA, col=alpha(coltmp[i],0.2))
  lines(y ~ x, plot.data, lwd=2, col = alpha(coltmp[i],0.8))
}
rm(data.tmp,i)
#legend('topleft',bty='n',legend = c('forest','farm'), lty=1, col=rev(coltmp),lwd=2,seg.len = 0.5)
#legend('topright',bty='n',legend = '(f)', cex =1.2)

plot(prop.decomp~weeks,cm,type='n',yaxt='n',xaxt='n',xlim=c(1,4),ylim=range(cm$prop.decomp),ann=F,bty='l')
title(xlab='days in stream')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:4,labels = seq(7,28,by=7))
title(ylab='proportion decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5))
means <- aggregate(prop.decomp~weeks*leaf.deployment, cm, FUN = 'mean')
se <- aggregate(prop.decomp~weeks*leaf.deployment, cm, FUN = 'sd')
means$se <- 1.96*se$prop.decomp/sqrt(nrow(se))
means$leaf.deployment <- as.factor(means$leaf.deployment)
coltmp <- c(cols[2],1)
for(i in 1:2){
  data.tmp <- subset(means, leaf.deployment == levels(means$leaf.deployment)[i])
  plot.data <- as.data.frame(approx(x = data.tmp$weeks, y = data.tmp$prop.decomp, n = 1000))
  plot.data$se <- as.data.frame(approx(x = data.tmp$weeks, y = data.tmp$se, n = 1000))[,2]
  polygon(x=c(plot.data$x,rev(plot.data$x)), y=c(plot.data$y+plot.data$se,rev(plot.data$y-plot.data$se)), border=NA, col=alpha(coltmp[i],0.2))
  lines(y ~ x, plot.data, lwd=2, col = alpha(coltmp[i],0.8))
}
rm(data.tmp,i)
#legend('topleft',bty='n',legend = c('home','away'), lty=1, col=rev(coltmp),lwd=2,seg.len = 0.5)
#legend('topright',bty='n',legend = '(g)', cex =1.2)

tp4 <- subset(cm, weeks == 4)
tp4$damage.area <- tp4$damage.area*100
plot(prop.decomp~damage.area,tp4,type='n',yaxt='n',xaxt='n',xlim=range(tp4$damage.area),ylim=range(tp4$prop.decomp),ann=F,bty='l')
title(xlab='leaf damage (%)')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='proportion decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,1,length.out = 5),labels=seq(0,1,length.out = 5))
#points(prop.decomp~damage.area,subset(tp4, weeks == 3),col=alpha(1,0.5),pch=16)
points(prop.decomp~damage.area,tp4,col=alpha(cols[3],0.5),pch=16)
#legend('topright',bty='n',legend = c('week 3','week 4'), pch=15, col=c(1,cols[3]))
#legend('topright',bty='n',legend=expression(italic('day 28')))
#legend('topright',bty='n',legend = '(h)', cex =1.2)

dev.off()