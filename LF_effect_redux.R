###############################################################################################
# RESPONSE: Invasive lionfish did not affect fish community structure across the Belizean Barrier Reef
# Authors: Hackerott, S., A. Valdivia C. E. Cox, N. J. Silbiger6, and J. F. Bruno
# Reanalysis: Ingeman, K. and M. Albins 
#Last edited: 7/1/2017
# NOTE: This code uses the raw count data at the transect level rather than mean density at the Site level 
##############################################################################################

library(lattice)
library(dplyr)
library(lme4)
library(ggplot2)

# convenience functions from b.bolker
plot.lmList <-
  function(object,ord.var=1,...) {
    ## Ord.var indicates which coefficient to 
    ##     sort the data frame of coefficients by.
    require(reshape)
    ## create a data.frame of coefs from list of glm fits.
    cL <- coef(object) 
    ## Adds a column for group ID --
    ## rownames (the variable upon which the fits were conditioned)
    ## -- here genotype. 
    cL$grp  <- rownames(cL) 
    if (is.numeric(ord.var) & length(ord.var)==1) {
      ## identify the ordering coefficient by name
      ord.var <- names(cL)[ord.var+1] 
      if (!ord.var %in% names(cL)) stop("unknown ordering variable")
      ## create a new order of factor levels for group ID
      cL$grp <- reorder(cL$grp,cL[[ord.var]])
    } else  
      ##otherwise just use the original ordering
      cL$grp <- reorder(cL$grp,ord.var) 
    ##"melt" or stack the data frame to 
    ##   allow a dotplot of each coefficient.
    dotplot(grp~value|variable,data=melt(cL),...) 
  }


qqmath.lmList <- function(object,...) {
  require(reshape)
  qqmath(~value|variable,data=melt(coef(object)),
         prepanel = prepanel.qqmathline,
         panel = function(x, ...) {
           panel.qqmathline(x, ...)
           panel.qqmath(x, ...)
         },
         scale=list(y=list(relation="free")))
}

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  (rdf <- nrow(model@frame)-model.df)
  rp <- residuals(model)
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE,log.p=TRUE)
  c(chisq=Pearson.chisq,ratio=prat,p=exp(pval))
}

sprint <- function(m) lme4:::printMer(m,correlation=FALSE)

locscaleplot <- function(model,col="black") {
  f <- fitted(model)
  r <- sqrt(abs(residuals(model)))  ## had better be Pearson resids
  plot(f,r,col=col) 
  L1 <- loess(r~f)
  fvec = seq(min(f),max(f),length.out=150)
  lines(fvec,predict(L1,fvec),col=2)
}

dfun <- function(x) {
  x$AIC <- x$AIC-min(x$AIC)
  names(x)[2] <- "dAIC"
  x
}

badcorr <- function(x) {
  cc <- attr(x,"correlation")
  diag(cc) <- 0
  any(abs((abs(cc)-1))<1e-5)
}
anybadcorr <- function(x) {
  any(sapply(VarCorr(x),badcorr))
}

locscaleplot <- function(model,col="black") {
  f <- fitted(model)
  r <- abs(residuals(model))
  plot(f,r,col=col) 
  L1 <- loess(r~f)
  fvec = seq(min(f),max(f),length.out=150)
  lines(fvec,predict(L1,fvec),col=2)
}

printvc <- function(m,digits=2,ctol=1e-3) {
  v <- VarCorr(m)
  prtfun <- function(x) {
    cc <- attr(x,"correlation")
    diag(cc) <- 0
    corstr <- ifelse(abs(abs(cc)-1)<ctol,"="," ")
    ss <- format(x,digits=digits) ## sprintf(fmt,x)
    ss <- paste(ss,corstr,sep="")
    m <- matrix(ss,nrow=nrow(x))
    m[upper.tri(m)] <- ""
    dimnames(m) <- dimnames(x)
    ## writeLines(apply(m,1,paste,collapse=" "))
    print(m,quote=FALSE)
  }
  for (i in seq_along(v)) {
    cat(names(v)[i],":\n",sep="")
    prtfun(v[[i]])
    cat("\n")
  }
}

plot.ICtab <- function(x,sort=TRUE,within) {
  z <- with(x,data.frame(n=attr(x,"row.names"),dAIC))
  if (sort) z <- transform(z,n=reorder(n,dAIC))
  if (!missing(within)) z$dAIC[z$dAIC>within] <- NA
  dotplot(n~dAIC,data=z)
}

##Load data file pData6.10 (THIS IS COUNT DATA)
pData6.10<-read.csv('pData6.10.csv')

str(pData6.10) # Year and Transect treated as integers instead of factors

# f.Year is factor variable
pData6.10$Year.fac <- factor(pData6.10$Year)

# convert Transect to factor variable
as.factor(pData6.10$Transect) 

# calculate *relative* area of total transects searched at each site each year 
# to offset survey intensities (doesn't need to be actual area)
# using length() for each factor combinations

pData6.10$Area.off <- ave(pData6.10$Transect, pData6.10[,c("Site","Year.fac")], FUN=length)
str(pData6.10)

# some exploratory plots
plot(Tot.Abund~LF.Abund.ha, data=pData6.10)

ggplot(pData6.10,aes(x=LF.Abund.ha,y=Tot.Abund, colour=Year.fac))+
  geom_point(alpha = 1/4, size = 2, position = "jitter") +
  geom_smooth(method=lm, se=FALSE) +
  theme_bw()
# no lionfish and few other fish in 2009
# not a linear relationship; try log abunda

ggplot(pData6.10,aes(x=LF.Abund.ha,y=log(Tot.Abund), colour=Year.fac))+
  geom_point(alpha = 1/4, size = 2, position = "jitter") +
  geom_smooth(method=lm, se=FALSE) +
  theme_bw()

ggplot(pData6.10,aes(x=LF.Abund.ha,y=Tot.Abund)) +
  geom_point(alpha = 1/4, size = 3) +
  geom_smooth(method=lm, se=FALSE) + 
  theme_bw()+
  facet_wrap(~Site, nrow=4)

ggplot(pData6.10,aes(x=LF.Abund.ha,y=log(Tot.Abund))) +
  geom_point(alpha = 1/4, size = 3) +
  geom_smooth(method=lm, se=FALSE) + 
  theme_bw()+
  facet_wrap(~Site, nrow=4)
# no real pattern in overall abundance or (log)abundance at the transect level
# consider dropping CLPA as they are not particularly important in diet and hugely variable
############################################################################


############################################################################
# sum total fish abundance and lionfish density at each SITE for each YEAR # 
############################################################################
grp.dat <- group_by(pData6.10, Site, Year.fac, Area.off)
dat6.10.sum <- summarize(grp.dat, Abund.site = sum(Tot.Abund), Lion.site = sum(LF.Abund.ha))

# visulaize relationship between lionfish and prey abundance 
plot(dat6.10.sum$Abund.site~dat6.10.sum$Lion.site) # large number of observations with no lionfish
# and a couple of very high values 

# visulaize relationship between and site and prey abundance
plot(dat6.10.sum$Abund.site~dat6.10.sum$Site) # unequal variance among sites esp HM

plot(dat6.10.sum$Abund.site~dat6.10.sum$Year.fac) # similar variances; consider outliers esp in 2012-2013; Species = CLPA


# other plots to look at patterns in data

ggplot(dat6.10.sum,aes(x=Lion.site,y=Abund.site))+
  geom_point() +
  geom_smooth(method=lm, se=FALSE) + 
  theme_bw()+
  facet_wrap(~Site, nrow=4) # don't see a clear pattern of lionfish effects at each site
# does suggest random site effects
# also suggests to me that lionfish averages at the site level are not a good predictor of Site-wide abundance 
ggplot(dat6.10.sum,aes(x=Lion.site,y=Abund.site, colour=Year.fac))+
  geom_point() +
  theme_bw()+
  facet_wrap(~Site, nrow=4)

#######################################################################################################
################ Drop CLPA  ###########################################################################
#######################################################################################################

noCLPA<-pData6.10[,-20]

dat.noCLPA<-noCLPA %>%
mutate(Abund.no = rowSums(.[8:39]))

grp.dat.no <- group_by(dat.noCLPA, Site, Year.fac, Area.off)
dat6.10.sum.no <- summarize(grp.dat.no, Abund.site = sum(Abund.no), Lion.site = sum(LF.Abund.ha))

#######################################################################################################
################ Fit models ###########################################################################
#######################################################################################################
# start with lm
m0 <- lm(Abund.site  ~ Year.fac + Lion.site + offset(Area.off), 
            data = dat6.10.sum)
plot(m0)
# obs 22, 23, 39 extreme values
m0.no <- lm(Abund.site  ~ Year.fac + Lion.site + offset(Area.off), 
            data = dat6.10.sum.no)
summary(m0.no)
plot(m0.no)

# fit fixed and random effects
m1 <- glmer(Abund.site  ~ Year.fac + Lion.site + 
                 offset(Area.off) + (1 | Site), 
               data = dat6.10.sum, family = poisson(), nAGQ = 0)
summary(m1) # lionfish effects non-sig and postive

m1.no <- glmer(Abund.site  ~ Year.fac + Lion.site + 
                 offset(Area.off) + (1 | Site), 
               data = dat6.10.sum.no, family = poisson(), nAGQ = 0)
summary(m1.no)

plot(resid(m1)~dat6.10.sum$Lion.site) # variance reverse horn-shaped 
plot(resid(m1)~dat6.10.sum$Year.fac)
plot(resid(m1)~dat6.10.sum$Site)

qqmath(m1) # issues: some observations fall outside the range of distribution
deviance(m1)

# visualize 

locscaleplot(m1)
#######################################################################################################
################ Tests for significance and effect size ################################################
#######################################################################################################


exp(cbind(fixef(m1), 
          confint(m1.sc, parm = "beta_", method = "Wald")))

exp(cbind(fixef(m1.sc), 
          confint(m1.sc, parm = "beta_", method = "boot")))

drop1(m1.sc, test = "Chisq")

