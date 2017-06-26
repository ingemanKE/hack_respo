

##Load data file psData0.5 and psData6.10


merge0.10<-merge(psData0.5[,c(1:2, 28)], psData6.10[,c(1:6, 41)])

psData0.10<-merge0.10[, c(1:2, 4:7, 3, 8)]

names(psData0.10)[8]<- "Tot.Abund6.10"

psData0.10$Tot.Abund0.10<-apply(psData0.10[,c(7,8)], 1, sum)



####Model Prey Abundance 0-5 with Site Means####

plot(psData0.10$Tot.Abund0.5~psData0.10$LF.Abund.ha)
abline(lm(psData0.10$Tot.Abund0.5~psData0.10$LF.Abund.ha))
cor.test(psData0.10$Tot.Abund0.5,psData0.10$LF.Abund.ha)
##NO significant correlation p-value = 0.3927

hist(psData0.10$Tot.Abund0.5)
shapiro.test(psData0.10$Tot.Abund0.5)
##p-value 2.942e-11 not normal

library(lme4)

##Test for best family and link with Site as Random factor
pab.gu.5<-lmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, REML=FALSE)
pab.gl.5<-glmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, family=gaussian("log"), nAGQ=0)
pab.gin.5<-glmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, family=gaussian("inverse"), nAGQ=0)
pab.gmin.5<-glmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, family=Gamma("inverse"), nAGQ=0)
pab.gmi.5<-glmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, family=Gamma("identity"), nAGQ=0)
pab.gml.5<-glmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, family=Gamma("log"), nAGQ=0)
pab.igi.5<-glmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, family=inverse.gaussian("identity"), nAGQ=0)
pab.igl.5<-glmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, family=inverse.gaussian("log"), nAGQ=0)
pab.igin.5<-glmer(Tot.Abund0.5~1+(1|Site), data=psData0.10, family=inverse.gaussian("inverse"), nAGQ=0)

AIC(pab.gu.5, pab.gl.5, pab.gin.5, pab.gmin.5, pab.gmi.5, pab.gml.5, pab.igin.5)

##Based on AIC- Model with Gamma distribution and log link


##Full model
p.abund.mod.s.5<-glmer(Tot.Abund0.5~scale(Year)+scale(LF.Abund.ha)+scale(Rugosity)+(1|Site),
                       data=psData0.10, family=Gamma("log"), nAGQ=0)

##Check model fit
plot(log(psData0.10$Tot.Abund0.5)~predict(p.abund.mod.s.5))
abline(lm(log(psData0.10$Tot.Abund0.5)~predict(p.abund.mod.s.5)))
cor(log(psData0.10$Tot.Abund0.5),predict(p.abund.mod.s.5))

##Check residuals
plot(resid(p.abund.mod.s.5))
plot(resid(p.abund.mod.s.5)~psData0.10$Year) 
plot(resid(p.abund.mod.s.5)~psData0.10$LF.Abund.ha)
plot(resid(p.abund.mod.s.5)~psData0.10$Rugosity) 

##Check VIF
vif.mer(p.abund.mod.s.5)
##No VIf>2 so no correlation

##Check residuals distribution
library(car)
qqPlot(residuals(p.abund.mod.s.5))
##Residuals are normal



####Model Prey Abundance 0-10 with Site Means####

plot(psData0.10$Tot.Abund0.10~psData0.10$LF.Abund.ha)
abline(lm(psData0.10$Tot.Abund0.10~psData0.10$LF.Abund.ha))
cor.test(psData0.10$Tot.Abund0.10,psData0.10$LF.Abund.ha)
##NO significant correlation p-value = 0.4179

hist(psData0.10$Tot.Abund0.10)
shapiro.test(psData0.10$Tot.Abund0.10)
##p-value 5.27e-11 not normal

library(lme4)

##Test for best family and link with Site as Random factor
pab.gu.10<-lmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, REML=FALSE)
pab.gl.10<-glmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, family=gaussian("log"), nAGQ=0)
pab.gin.10<-glmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, family=gaussian("inverse"), nAGQ=0)
pab.gmin.10<-glmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, family=Gamma("inverse"), nAGQ=0)
pab.gmi.10<-glmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, family=Gamma("identity"), nAGQ=0)
pab.gml.10<-glmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, family=Gamma("log"), nAGQ=0)
pab.igi.10<-glmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, family=inverse.gaussian("identity"), nAGQ=0)
pab.igl.10<-glmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, family=inverse.gaussian("log"), nAGQ=0)
pab.igin.10<-glmer(Tot.Abund0.10~1+(1|Site), data=psData0.10, family=inverse.gaussian("inverse"), nAGQ=0)

AIC(pab.gu.10, pab.gl.10, pab.gin.10, pab.gmin.10, pab.gmi.10, pab.gml.10, pab.igin.10)

##Based on AIC- Model with Gamma distribution and log link


##Full model
p.abund.mod.s.10<-glmer(Tot.Abund0.10~scale(Year)+scale(LF.Abund.ha)+scale(Rugosity)+(1|Site),
                        data=psData0.10, family=Gamma("log"), nAGQ=0)

##Check model fit
plot(log(psData0.10$Tot.Abund0.10)~predict(p.abund.mod.s.10))
abline(lm(log(psData0.10$Tot.Abund0.10)~predict(p.abund.mod.s.10)))
cor(log(psData0.10$Tot.Abund0.10),predict(p.abund.mod.s.10))

##Check residuals
plot(resid(p.abund.mod.s.10))
plot(resid(p.abund.mod.s.10)~psData0.10$Year) 
plot(resid(p.abund.mod.s.10)~psData0.10$LF.Abund.ha)
plot(resid(p.abund.mod.s.10)~psData0.10$Rugosity) 

##Check VIF
vif.mer(p.abund.mod.s.10)
##No VIf>2 so no correlation

##Check residuals distribution
library(car)
qqPlot(residuals(p.abund.mod.s.10))
##Residuals are normal
