### Code for experimental project on drivers of decomposition in Kibale streams ###
### Vincent Fugère, 2019

rm(list=ls())

library(tidyverse)
library(RColorBrewer)
library(glmmTMB)
library(magrittr)
library(plotrix)
library(party)
library(devtools)
library(performance)
source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/vif.R')
source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

cols <- brewer.pal(5, 'Dark2')

#### load and format data ####

data <- read_csv('~/Google Drive/Recherche/PhD/manuscripts/caterpillar/decompProj/leafbags.csv') %>% filter(!is.na(rep.nb))

dmg <-  read_csv('~/Google Drive/Recherche/PhD/manuscripts/caterpillar/decompProj/caterpillar.csv') %>%
  select(leaf.nb,damage.area)

data <- inner_join(data, dmg, by = 'leaf.nb') #loosing one replicate - the leaf fragment with a missing picture
rm(dmg)

#to conduct dmg analysis at whole-leaf or whole-tree level
data <- data %>% group_by(leaf) %>% mutate(avg.dmg = mean(damage.area)) %>% ungroup
data <- data %>% group_by(tree) %>% mutate(tree.avg.dmg = mean(damage.area)) %>% ungroup

data$leaf.lu <- 'forest'
data$leaf.lu[data$leaf.origin %in% c('miranga','bugembe')] <- 'farm'
data$leaf.lu <- as.factor(data$leaf.lu)

data$prop <- data$prct.mass.rem/100
data$prop.decomp <- 1-data$prop

#leaching un-adjusted for k calculation
data$exp.dec <- with(data, (AFDM.after / dry.weight)*100)

#transformation based on: https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
data$rv <- 1-(data$prop*(nrow(data)-1)+0.5) / nrow(data)
#1- because we want prop decomposed instead of prop left

data$dmg <- rescale(data$damage.area, c(0,1))
data$lu <- as.factor(data$land.use)
data$lu <- relevel(data$lu, ref = 'forest')
data$ha <- as.factor(data$leaf.deployment)
data$ha <- relevel(data$ha, ref = 'home')
data$wks <- rescale(data$weeks, c(0,1))
data$wtr <- as.factor(data$watershed)
data$lf.lu <- relevel(data$leaf.lu, ref = 'forest')
data$lf.dmg <- rescale(data$avg.dmg, c(0,1))
data$tree.dmg <- rescale(data$tree.avg.dmg, c(0,1))

#### is leaf mass different among treatments? ####

#data$dm <- data$dry.weight
#
# mean(data$dm)
# sd(data$dm)
# sd(data$dm)/sqrt(nrow(data)) #se
# range(data$dm)
# 
# anova(lm(dm~lu,data))
# anova(lm(dm~ha,data))
# summary(lm(dm~dmg,data))
# plot(dm~dmg,data)
# #relationship driven by high-leverage outlier
# data2 <- data[data$dmg < 0.95,]
# anova(lm(dm~lu,data2))
# anova(lm(dm~ha,data2))
# summary(lm(dm~dmg,data2))

#### setting up models ####

#collinearity

corvif(data[,c('dmg','lu','ha','wks','wtr','lf.lu','lf.dmg','tree.dmg')]) 
#fragment damage and leaf total damage are obviously correlated
#also correlated with tree damage, and leaf land use
#will test separately. No dmg var has an effect; chose fragment damage in final models because cleaner

#model with time

f0 <- formula(rv ~ 1 + wks + wtr + dmg + ha + lu + lf.lu)
f1 <- update(f0, . ~ wks * .)
f2 <- update(f1, . ~ . + (wks|site) + (wks|leaf.origin/tree/leaf))

#model with final time point only

fx <- formula(rv ~ 1 + dmg + wtr + ha + lu + lf.lu) 
fx <- update(fx, . ~ . + (1|site) + (1|leaf.origin/tree/leaf))

#### models ####

data$mass.init <- rescale(data$leaching.adj.dry.weight, c(0,1))

fm <- filter(data, mesh.type == 'fine')
cm <- filter(data, mesh.type == 'coarse')

#calculating k/day
summary(nls(exp.dec~100*(exp(1)^(-k*(weeks*7))),fm,start=list('k'= 0.01)))
summary(nls(exp.dec~100*(exp(1)^(-k*(weeks*7))),cm,start=list('k'= 0.01)))
summary(nls(exp.dec~100*(exp(1)^(-k*(weeks*7))),subset(cm, land.use == 'farm'),start=list('k'= 0.01)))
summary(nls(exp.dec~100*(exp(1)^(-k*(weeks*7))),subset(cm, land.use == 'forest'),start=list('k'= 0.01)))

ts.mod.fm <- glmmTMB(f2, fm, family=beta_family(link = "logit"))
summary(ts.mod.fm)

ts.mod.cm <- glmmTMB(f2, cm, family=beta_family(link = "logit"))
summary(ts.mod.cm)

# #saving model results for Table 1
# rownames_to_column(as.data.frame(summary(ts.mod.fm)$coefficients[1])) %>% write_csv(., '~/Desktop/fmmod.csv')
# rownames_to_column(as.data.frame(summary(ts.mod.cm)$coefficients[1])) %>% write_csv(., '~/Desktop/cmmod.csv')

# #interaction model
# f3 <- formula(rv ~ wks + dmg + ha + lu + (wks | site) + (wks | leaf.origin/tree/leaf) +
#                 wks:dmg + wks:ha + wks:lu + dmg:ha + dmg:lu + ha:lu +
#                 dmg:ha:wks + dmg:lu:wks + ha:lu:wks)
# summary(glmmTMB(f3, fm, family=beta_family(link = "logit"))) #no higher-order interactions
# summary(glmmTMB(f3, cm, family=beta_family(link = "logit")))
# #does not converge, going Bayesian
# library(brms)
# m1 <- brm(f3, cm, family=beta_family(link = "logit"))
# summary(m1)
# plot(m1)
# #no higher-order interactions

fm.end <- filter(fm, weeks == 4)
cm.end <- filter(cm, weeks == 4)

# #regression trees
# fit.fm <- ctree(rv ~ dmg+lu+ha+wks+wtr+leaf.lu+lf.lu, data = fm.end, controls = ctree_control(minsplit = 1, testtype = 'MonteCarlo'))
# plot(fit.fm, inner_panel=node_inner(fit.fm,pval = T),
#      terminal_panel=node_boxplot(fit.fm, width=0.4,fill='white',ylines=3,id=F))
# fit.cm <- ctree(rv ~ dmg+lu+ha+wks+wtr+leaf.lu+lf.lu, data = cm.end, controls = ctree_control(minsplit = 1, testtype = 'MonteCarlo'))
# plot(fit.cm, inner_panel=node_inner(fit.cm,pval = T),
#      terminal_panel=node_boxplot(fit.cm, width=0.4,fill='white',ylines=3,id=F))
# #the only thing that matters is land use for CM bags

t4.mod.fm <- glmmTMB(fx, fm.end, family=beta_family(link = "logit"))
summary(t4.mod.fm)

# #removing site random effect and treating as fixed effect
# fx2 <- formula(rv ~ 1 + dmg + wtr + ha + lu + lf.lu)
# fx2 <- update(fx2, . ~ . + (1|leaf.origin/tree/leaf))
# t4.mod.fm2 <- glmmTMB(fx2, fm.end, family=beta_family(link = "logit"))
# summary(t4.mod.fm2)
# #effect of land use only significant in one watershed
# 
# plot(fitted(t4.mod.fm)~fm.end$rv, pch=16, col=fm.end$site_f)
# plot(fitted(t4.mod.fm2)~fm.end$rv, pch=16, col=fm.end$site_f)
# #fit is identical

t4.mod.cm <- glmmTMB(fx, cm.end, family=beta_family(link = "logit"))
summary(t4.mod.cm)

fm$prop.decomp <- fm$prop.decomp*100
cm$prop.decomp <- cm$prop.decomp*100

# #interaction models with factors of interest (not enough data to add all factors)
# fx.i <- formula(rv ~ dmg + ha + lu + (1 | site) + (1 | leaf.origin/tree/leaf) + 
#                   dmg:ha + dmg:lu + ha:lu)
# summary(glmmTMB(fx.i, fm.end, family=beta_family(link = "logit")))
# summary(glmmTMB(fx.i, cm.end, family=beta_family(link = "logit")))
# #convergence problems again
# m1 <- brm(fx.i, cm.end, family=beta_family(link = "logit"),iter=4000,thin=2)
# summary(m1)
# plot(m1)
# #no significant two-way interactions

# # effect of starting AFDM on decomp
# 
# summary(lm(rv~mass.init*weeks,fm))
# summary(lm(rv~mass.init*weeks,cm))
# summary(lm(rv~mass.init,fm.end))
# summary(lm(rv~mass.init,cm.end))
# 
# pdf('~/Desktop/FigS1.pdf',width=5.5,height = 2.8,pointsize = 8)
# par(mfrow=c(1,2),cex=1)
# plot(rv~mass.init,fm.end,xlab = 'initial mass', ylab = 'proportion decomposed',bty='l')
# title(main='fine mesh',cex=1)
# plot(rv~mass.init,cm.end,xlab = 'initial mass', ylab = 'proportion decomposed',bty='l')
# title(main='coarse mesh',cex=1)
# dev.off()
# par(mfrow=c(1,1))

### wtrshed as random effect as asked by Editor, remove leaf land use

f0 <- formula(rv ~ 1 + wks + dmg + ha + lu)
f1 <- update(f0, . ~ wks * .)
f2 <- update(f1, . ~ . + (1|watershed/site) + (1|leaf.origin/tree/leaf))

ts.mod.fm <- glmmTMB(f2, fm, family=beta_family(link = "logit"))
summary(ts.mod.fm)

ts.mod.cm <- glmmTMB(f2, cm, family=beta_family(link = "logit"))
summary(ts.mod.cm)

#saving model results for Table 1
rownames_to_column(as.data.frame(summary(ts.mod.fm)$coefficients[1])) %>% 
  bind_rows(rownames_to_column(as.data.frame(summary(ts.mod.cm)$coefficients[1]))) %>%
  write_csv(., '~/Desktop/tsmod.csv')

fx<- formula(rv ~ dmg + ha + lu + (1 | wtr/site) + (1 | leaf.origin/tree/leaf))

t4.mod.fm <- glmmTMB(fx, fm.end, family=beta_family(link = "logit"))
summary(t4.mod.fm)

t4.mod.cm <- glmmTMB(fx, cm.end, family=beta_family(link = "logit"))
summary(t4.mod.cm)
plot(resid(t4.mod.cm)~fitted(t4.mod.cm))
plot(fitted(t4.mod.cm) ~ cm.end$rv)

# #### veryfing that I get similar coefficients with a bayesian model with strong-ish priors for random effects
# 
# library(brms)
# get_prior(fx,cm.end,family='beta')
# priorz <- prior_string("cauchy(0, 1)", class = "sd") #helps random effects converging on small sd
# 
# mod1 <- brm(f2, fm, family='beta',cores=2, iter=8000,prior=priorz)
# mod2 <- brm(f2, cm, family='beta',cores=2, iter=8000,prior=priorz)
# mod3 <- brm(fx, fm.end, family='beta',cores=2, iter=8000,prior=priorz)
# mod4 <- brm(fx, cm.end, family='beta',cores=2, iter=8000,prior=priorz)
# 
# plot(mod1)
# summary(mod1);summary(ts.mod.fm)
# plot(fitted(mod1)[,1] ~ fm$rv)
# pp_check(mod1)
# 
# plot(mod2)
# summary(mod2);summary(ts.mod.cm)
# plot(fitted(mod2)[,1] ~ cm$rv)
# pp_check(mod2)
# 
# plot(mod3)
# summary(mod3);summary(t4.mod.fm)
# plot(fitted(mod3)[,1] ~ fm.end$rv)
# pp_check(mod3)
# 
# plot(mod4)
# summary(mod4);summary(t4.mod.cm)
# plot(fitted(mod4)[,1] ~ cm.end$rv)
# pp_check(mod4)

#### Figure 2 ####

#pdf('~/Desktop/Fig2.pdf',width=4.5,height = 3.5,pointsize = 8)
par(mfrow=c(2,3),mar=c(4,4,1,1),cex=1)

#fm

plot(prop.decomp~weeks,fm,type='n',yaxt='n',xaxt='n',xlim=c(1,4),ylim=range(fm$prop.decomp),ann=F,bty='l')
title(xlab='days in stream')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:4,labels = seq(7,28,by=7))
title(ylab='% decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,100,length.out = 5),labels=seq(0,100,length.out = 5))
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
legend('topright',bty='n',legend = '(a)', cex =1.2)

plot(prop.decomp~weeks,fm,type='n',yaxt='n',xaxt='n',xlim=c(1,4),ylim=range(fm$prop.decomp),ann=F,bty='l')
title(xlab='days in stream')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:4,labels = seq(7,28,by=7))
title(ylab='% decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,100,length.out = 5),labels=seq(0,100,length.out = 5))
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

tp4 <- subset(fm, weeks == 4)
tp4$damage.area <- tp4$damage.area*100
plot(prop.decomp~damage.area,tp4,type='n',yaxt='n',xaxt='n',xlim=range(tp4$damage.area),ylim=range(tp4$prop.decomp),ann=F,bty='l')
title(xlab='leaf damage (%)')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='% decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,100,length.out = 5),labels=seq(0,100,length.out = 5))
points(prop.decomp~damage.area,tp4,col=alpha(cols[3],0.5),pch=16)
legend('topright',bty='n',legend=expression(italic('day 28')))

#cm

plot(prop.decomp~weeks,cm,type='n',yaxt='n',xaxt='n',xlim=c(1,4),ylim=range(cm$prop.decomp),ann=F,bty='l')
title(xlab='days in stream')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:4,labels = seq(7,28,by=7))
title(ylab='% decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,100,length.out = 5),labels=seq(0,100,length.out = 5))
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
legend('topright',bty='n',legend = '(b)', cex =1.2)

plot(prop.decomp~weeks,cm,type='n',yaxt='n',xaxt='n',xlim=c(1,4),ylim=range(cm$prop.decomp),ann=F,bty='l')
title(xlab='days in stream')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:4,labels = seq(7,28,by=7))
title(ylab='% decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,100,length.out = 5),labels=seq(0,100,length.out = 5))
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

tp4 <- subset(cm, weeks == 4)
tp4$damage.area <- tp4$damage.area*100
plot(prop.decomp~damage.area,tp4,type='n',yaxt='n',xaxt='n',xlim=range(tp4$damage.area),ylim=range(tp4$prop.decomp),ann=F,bty='l')
title(xlab='leaf damage (%)')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='% decomposed')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(0,100,length.out = 5),labels=seq(0,100,length.out = 5))
points(prop.decomp~damage.area,tp4,col=alpha(cols[3],0.5),pch=16)

#dev.off()

#### Figure 3 ####

pdf('~/Desktop/Fig3.pdf',width=3,height = 2.5,pointsize = 8)
par(mfrow=c(1,1),mar=c(4,7,1,1),cex=1)

plot(0,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(-6,2),ylim=c(0.5,6.5))
#axis(2,cex.axis=1,lwd=0,lwd.ticks=0,at=1:5,labels = rev(c('farm vs. forest','away vs. home','leaf damage','leaf land use','watershed')),las=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=0,at=1:6,labels = rev(c('farm vs. forest','away vs. home','leaf damage','farm vs. forest','away vs. home','leaf damage')),las=1)
abline(v=0,lty=2,lwd=1)
abline(h=3.5,lty=1,lwd=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=3.5,labels = '')
title(xlab='parameter estimate')
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)

#getting coefs: fm
coefs <- rownames_to_column(as.data.frame(summary(t4.mod.fm)$coefficients[1])) %>%
  rename(coef = rowname, value = cond.Estimate, se = cond.Std..Error) %>% select(coef:se) %>%
  filter(coef != '(Intercept)')
coefs %<>% mutate('lwr' = value-1.96*se, 'upr' = value+1.96*se)
#coefs <- coefs[c(2,5,1,3,4),]

for(w in 4:6){
  pt <- coefs[w-3,'value']
  lwr <- coefs[w-3,'lwr']
  upr <- coefs[w-3,'upr']
  if(0 < upr & 0 > lwr){
    alph <- 0.3
  } else {
    alph <- 1
  }
  segments(x0=lwr,x1=upr,y0=w,y1=w,col=alpha(cols[w-3],alph),lwd=3)
  points(x=pt,y=w,pch=16,col=alpha(cols[w-3],alph),cex=2)
}
legend('topleft',bty='n',legend=expression(italic('fine-mesh')))

#getting coefs: cm
coefs <- rownames_to_column(as.data.frame(summary(t4.mod.cm)$coefficients[1])) %>%
  rename(coef = rowname, value = cond.Estimate, se = cond.Std..Error) %>% select(coef:se) %>%
  filter(coef != '(Intercept)') 
coefs %<>% mutate('lwr' = value-1.96*se, 'upr' = value+1.96*se)

for(w in 1:3){
  pt <- coefs[w,'value']
  lwr <- coefs[w,'lwr']
  upr <- coefs[w,'upr']
  if(0 < upr & 0 > lwr){
    alph <- 0.3
  } else {
    alph <- 1
  }
  segments(x0=lwr,x1=upr,y0=w,y1=w,col=alpha(cols[w],alph),lwd=3)
  points(x=pt,y=w,pch=16,col=alpha(cols[w],alph),cex=2)
}
legend('bottomleft',bty='n',legend=expression(italic('coarse-mesh')))

dev.off()