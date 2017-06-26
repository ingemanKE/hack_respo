
###############################################################################################
# Invasive lionfish did not affect fish community structure across the Belizean Barrier Reef
# Authors: Hackerott, S., A. Valdivia C. E. Cox, N. J. Silbiger6, and J. F. Bruno
# Last edited: 12/2/2016
# This code runs through the all the models in the main text of the maniscript
##############################################################################################

####Main Analysis: Effect of lionfish on prey species at site level####

##Load data file psData6.10
psData6.10<-read.csv('psData6.10.csv')
head(psData6.10)
##Give meaningful rownames of Site and Year
rownames(psData6.10)<-paste(psData6.10$Site, psData6.10$Year-2000, sep=".")

## functions --- calculates variance inflation factors
vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}



####Model Prey Abundance with Site Means####

plot(psData6.10$Tot.Abund~psData6.10$LF.Abund.ha)
abline(lm(psData6.10$Tot.Abund~psData6.10$LF.Abund.ha))
cor.test(psData6.10$Tot.Abund,psData6.10$LF.Abund.ha)
##NO significant correlation p-value = 0.7423

hist(psData6.10$Tot.Abund)
shapiro.test(psData6.10$Tot.Abund)
##p-value 5.583e-06 not normal

library(lme4)

##Test for best family and link with Site as Random factor
pab.gu<-lmer(Tot.Abund~1+(1|Site), data=psData6.10, REML=FALSE)
pab.gl<-glmer(Tot.Abund~1+(1|Site), data=psData6.10, family=gaussian("log"), nAGQ=0)
pab.gin<-glmer(Tot.Abund~1+(1|Site), data=psData6.10, family=gaussian("inverse"), nAGQ=0)
pab.gmin<-glmer(Tot.Abund~1+(1|Site), data=psData6.10, family=Gamma("inverse"), nAGQ=0)
pab.gmi<-glmer(Tot.Abund~1+(1|Site), data=psData6.10, family=Gamma("identity"), nAGQ=0)
pab.gml<-glmer(Tot.Abund~1+(1|Site), data=psData6.10, family=Gamma("log"), nAGQ=0)
pab.igi<-glmer(Tot.Abund~1+(1|Site), data=psData6.10, family=inverse.gaussian("identity"), nAGQ=0)
pab.igl<-glmer(Tot.Abund~1+(1|Site), data=psData6.10, family=inverse.gaussian("log"), nAGQ=0)
pab.igin<-glmer(Tot.Abund~1+(1|Site), data=psData6.10, family=inverse.gaussian("inverse"), nAGQ=0)

AIC(pab.gu, pab.gl, pab.gin, pab.gmin, pab.gmi, pab.gml, pab.igi, pab.igl, pab.igin)

##Based on AIC- Model with Gamma distribution and log link


##Full model
p.abund.mod.s<-glmer(Tot.Abund~scale(Year)+scale(LF.Abund.ha)+scale(Rugosity)+(1|Site),
                     data=psData6.10, family=Gamma("log"), nAGQ=0)

##Check model fit
plot(log(psData6.10$Tot.Abund)~predict(p.abund.mod.s))
abline(lm(log(psData6.10$Tot.Abund)~predict(p.abund.mod.s)))
cor(log(psData6.10$Tot.Abund),predict(p.abund.mod.s))

##Check residuals
plot(resid(p.abund.mod.s))
plot(resid(p.abund.mod.s)~psData6.10$Year) 
plot(resid(p.abund.mod.s)~psData6.10$LF.Abund.ha)
plot(resid(p.abund.mod.s)~psData6.10$Rugosity) 

##Check VIF
vif.mer(p.abund.mod.s)
##No VIf>2 so no correlation

##Check residuals distribution
library(car)
qqPlot(residuals(p.abund.mod.s))
##Residuals are normal



####Model Prey Richness with Site Means####
plot(psData6.10$Rich~psData6.10$LF.Abund.ha)
abline(lm(psData6.10$Rich~psData6.10$LF.Abund.ha))
cor.test(psData6.10$Rich,psData6.10$LF.Abund.ha)
##Significant positive correlation 0.2986292; t = 2.383, df = 58, p-value = 0.02047

hist(psData6.10$Rich)
shapiro.test(psData6.10$Rich)
##p-value 0.0472 normal

##Test for best family and link with Site as Random factor
prich.gu<-lmer(Rich~1+(1|Site), data=psData6.10, REML=FALSE)
prich.gl<-glmer(Rich~1+(1|Site), data=psData6.10, family=gaussian("log"), nAGQ=0)
prich.gin<-glmer(Rich~1+(1|Site), data=psData6.10, family=gaussian("inverse"), nAGQ=0)
prich.gmin<-glmer(Rich~1+(1|Site), data=psData6.10, family=Gamma("inverse"), nAGQ=0)
prich.gmi<-glmer(Rich~1+(1|Site), data=psData6.10, family=Gamma("identity"), nAGQ=0)
prich.gml<-glmer(Rich~1+(1|Site), data=psData6.10, family=Gamma("log"), nAGQ=0)
prich.igi<-glmer(Rich~1+(1|Site), data=psData6.10, family=inverse.gaussian("identity"), nAGQ=0)
prich.igl<-glmer(Rich~1+(1|Site), data=psData6.10, family=inverse.gaussian("log"), nAGQ=0)
prich.igin<-glmer(Rich~1+(1|Site), data=psData6.10, family=inverse.gaussian("inverse"), nAGQ=0)

AIC(prich.gu, prich.gl, prich.gin, prich.gmin, prich.gmi, prich.gml, prich.igl, prich.igin)

##Normal and Based on AIC-Model with Gaussian


##Full model
p.rich.mod.s<-lmer(Rich~scale(Year)+scale(LF.Abund.ha)+scale(Rugosity)+(1|Site), data=psData6.10)

##Check model fit
plot(psData6.10$Rich~predict(p.rich.mod.s))
abline(lm(psData6.10$Rich~predict(p.rich.mod.s)))
cor(psData6.10$Rich,predict(p.rich.mod.s))

##Check residuals
plot(resid(p.rich.mod.s))
plot(resid(p.rich.mod.s)~psData6.10$Year) 
plot(resid(p.rich.mod.s)~psData6.10$LF.Abund.ha)
plot(resid(p.rich.mod.s)~psData6.10$Rugosity) 

##Check VIF
vif.mer(p.rich.mod.s)
##No VIf>2 so no correlation

##Check residuals distribution
library(car)
qqPlot(residuals(p.rich.mod.s))
##Residuals are normal


####Model Prey Pomacentridae with Site Means####

plot(psData6.10$Pomacentridae~psData6.10$LF.Abund.ha)
abline(lm(psData6.10$Pomacentridae~psData6.10$LF.Abund.ha))
cor.test(psData6.10$Pomacentridae,psData6.10$LF.Abund.ha)
##NO significant correlation p-value = 0.6594

hist(psData6.10$Pomacentridae)
shapiro.test(psData6.10$Pomacentridae)
##p-value 5.951e-09 not normal

library(lme4)

##Test for best family and link with Site as Random factor
p.dams.gu<-lmer(Pomacentridae~1+(1|Site), data=psData6.10, REML=FALSE)
p.dams.gl<-glmer(Pomacentridae~1+(1|Site), data=psData6.10, family=gaussian("log"), nAGQ=0)
p.dams.gin<-glmer(Pomacentridae~1+(1|Site), data=psData6.10, family=gaussian("inverse"), nAGQ=0)
p.dams.gmin<-glmer(Pomacentridae~1+(1|Site), data=psData6.10, family=Gamma("inverse"), nAGQ=0)
p.dams.gmi<-glmer(Pomacentridae~1+(1|Site), data=psData6.10, family=Gamma("identity"), nAGQ=0)
p.dams.gml<-glmer(Pomacentridae~1+(1|Site), data=psData6.10, family=Gamma("log"), nAGQ=0)
p.dams.igi<-glmer(Pomacentridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("identity"), nAGQ=0)
p.dams.igl<-glmer(Pomacentridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("log"), nAGQ=0)
p.dams.igin<-glmer(Pomacentridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("inverse"), nAGQ=0)

AIC(p.dams.gu, p.dams.gl, p.dams.gin, p.dams.gmin, p.dams.gmi, p.dams.gml, p.dams.igin)

##Based on AIC-Model with Gamma distribution and log link


##Full model
p.dams.mod.s<-glmer(Pomacentridae~scale(Year)+scale(LF.Abund.ha)+scale(Rugosity)+(1|Site),
                    data=psData6.10, family=Gamma("log"), nAGQ=0)

##Check model fit
plot(log(psData6.10$Pomacentridae)~predict(p.dams.mod.s))
abline(lm(log(psData6.10$Pomacentridae)~predict(p.dams.mod.s)))
cor(log(psData6.10$Pomacentridae),predict(p.dams.mod.s))

##Check residuals
plot(resid(p.dams.mod.s))
plot(resid(p.dams.mod.s)~psData6.10$Year) 
plot(resid(p.dams.mod.s)~psData6.10$LF.Abund.ha)
plot(resid(p.dams.mod.s)~psData6.10$Rugosity) 

##Check VIF
vif.mer(p.dams.mod.s)
##No VIf>2 so no correlation

##Check residuals distribution
library(car)
qqPlot(residuals(p.dams.mod.s))
##Residuals are normal



####Model Prey Labridae with Site Means####

plot(psData6.10$Labridae~psData6.10$LF.Abund.ha)
abline(lm(psData6.10$Labridae~psData6.10$LF.Abund.ha))
cor.test(psData6.10$Labridae,psData6.10$LF.Abund.ha)
##NO significant correlation p-value = 0.8948

hist(psData6.10$Labridae)
shapiro.test(psData6.10$Labridae)
##p-value 1.679e-09 not normal

library(lme4)

##Test for best family and link with Site as Random factor
p.wras.gu<-lmer(Labridae~1+(1|Site), data=psData6.10, REML=FALSE)
p.wras.gl<-glmer(Labridae~1+(1|Site), data=psData6.10, family=gaussian("log"), nAGQ=0)
p.wras.gin<-glmer(Labridae~1+(1|Site), data=psData6.10, family=gaussian("inverse"), nAGQ=0)
p.wras.gmin<-glmer(Labridae~1+(1|Site), data=psData6.10, family=Gamma("inverse"), nAGQ=0)
p.wras.gmi<-glmer(Labridae~1+(1|Site), data=psData6.10, family=Gamma("identity"), nAGQ=0)
p.wras.gml<-glmer(Labridae~1+(1|Site), data=psData6.10, family=Gamma("log"), nAGQ=0)
p.wras.igi<-glmer(Labridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("identity"), nAGQ=0)
p.wras.igl<-glmer(Labridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("log"), nAGQ=0)
p.wras.igin<-glmer(Labridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("inverse"), nAGQ=0)

AIC(p.wras.gu, p.wras.gl, p.wras.gin, p.wras.gmin, p.wras.gmi, p.wras.gml, p.wras.igin)

##Based on AIC-Model with Gamma distribution and log link


##Full model
p.wras.mod.s<-glmer(Labridae~scale(Year)+scale(LF.Abund.ha)+scale(Rugosity)+(1|Site),
                    data=psData6.10, family=Gamma("log"), nAGQ=0)

##Check model fit
plot(log(psData6.10$Labridae)~predict(p.wras.mod.s))
abline(lm(log(psData6.10$Labridae)~predict(p.wras.mod.s)))
cor(log(psData6.10$Labridae),predict(p.wras.mod.s))

##Check residuals
plot(resid(p.wras.mod.s))
plot(resid(p.wras.mod.s)~psData6.10$Year) 
plot(resid(p.wras.mod.s)~psData6.10$LF.Abund.ha)
plot(resid(p.wras.mod.s)~psData6.10$Rugosity) 

##Check VIF
vif.mer(p.wras.mod.s)
##No VIf>2 so no correlation

##Check residuals distribution
library(car)
qqPlot(residuals(p.wras.mod.s))
##Residuals are normal




####Model Prey Scaridae with Site Means####

plot(psData6.10$Scaridae~psData6.10$LF.Abund.ha)
abline(lm(psData6.10$Scaridae~psData6.10$LF.Abund.ha))
cor.test(psData6.10$Scaridae,psData6.10$LF.Abund.ha)
##NO significant correlation p-value = 0.06856

hist(psData6.10$Scaridae)
shapiro.test(psData6.10$Scaridae)
##p-value 0.0007884 not normal

library(lme4)

##Test for best family and link with Site as Random factor
p.parr.gu<-lmer(Scaridae~1+(1|Site), data=psData6.10, REML=FALSE)
p.parr.gl<-glmer(Scaridae~1+(1|Site), data=psData6.10, family=gaussian("log"), nAGQ=0)
p.parr.gin<-glmer(Scaridae~1+(1|Site), data=psData6.10, family=gaussian("inverse"), nAGQ=0)
p.parr.gmin<-glmer(Scaridae~1+(1|Site), data=psData6.10, family=Gamma("inverse"), nAGQ=0)
p.parr.gmi<-glmer(Scaridae~1+(1|Site), data=psData6.10, family=Gamma("identity"), nAGQ=0)
p.parr.gml<-glmer(Scaridae~1+(1|Site), data=psData6.10, family=Gamma("log"), nAGQ=0)
p.parr.igi<-glmer(Scaridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("identity"), nAGQ=0)
p.parr.igl<-glmer(Scaridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("log"), nAGQ=0)
p.parr.igin<-glmer(Scaridae~1+(1|Site), data=psData6.10, family=inverse.gaussian("inverse"), nAGQ=0)

AIC(p.parr.gu, p.parr.gl, p.parr.gin, p.parr.gmin, p.parr.gmi, p.parr.gml, p.parr.igin)

##Based on AIC-Model with Gamma distribution and log link


##Full model
p.parr.mod.s<-glmer(Scaridae~scale(Year)+scale(LF.Abund.ha)+scale(Rugosity)+(1|Site),
                    data=psData6.10, family=Gamma("log"), nAGQ=0)

##Check model fit
plot(log(psData6.10$Scaridae)~predict(p.parr.mod.s))
abline(lm(log(psData6.10$Scaridae)~predict(p.parr.mod.s)))
cor(log(psData6.10$Scaridae),predict(p.parr.mod.s))

##Check residuals
plot(resid(p.parr.mod.s))
plot(resid(p.parr.mod.s)~psData6.10$Year) 
plot(resid(p.parr.mod.s)~psData6.10$LF.Abund.ha)
plot(resid(p.parr.mod.s)~psData6.10$Rugosity) 

##Check VIF
vif.mer(p.parr.mod.s)
##No VIf>2 so no correlation

##Check residuals distribution
library(car)
qqPlot(residuals(p.parr.mod.s))
##Residuals are normal



####Community Composition Analysis- Prey Species with Transect Data####

##Load data file: pData6.10
pData6.10<-read.csv('pData6.10.csv')

library(vegan)

##Permanova: testing the effect of year and lionfish abundance on community composition; stratafied by site

##Remove rows with zero fish
perm.pData<-pData6.10[-c(which(pData6.10$Tot.Abund==0)),]

perm.prey<-adonis(perm.pData[, c(8:40)]~perm.pData$Year+perm.pData$LF.Abund.ha, 
                  strata=perm.pData$Site, permutations=10000)
##p-values- Year: 9.999e-5, Lionfish: 0.5378





