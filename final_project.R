##### 555 Spatial Epidemiology Final Project #####

#### Step 1: Input Data #######

#load data
library(haven) #allows read_dta
hiv.data <- read_dta("C:/Users/steuben/Desktop/555/ZWE_DHS7_2015_HIV_Y2016M12D13.DTA")
hhm.data <- read_dta("C:/Users/steuben/Desktop/555/ZWE_DHS7_2015_HHM_Y2016M12D13.DTA")

#Keep only the needed variables
hhm.var <- c("hhid", "hvidx", "hv001", "hv002", "hv003", "hv004", "hv005","hv021","hv022", "hv023", 
             "hv024", "hv025", "hv105", "hv104")
hhm.names <- c("id", "hivline", "hivclust", "hivnumb", "respline", "area_unit","house_weight", "psu",
               "strata_error", "strata_sample_design", "province", "urban_rural", "age", "sex")
hhm.data <- hhm.data[,hhm.var]
colnames(hhm.data) <- hhm.names

hiv.var <- c("hivclust", "hivnumb", "hivline", "hiv03","hiv05")
hiv.data <- hiv.data[,hiv.var]

#Merge to create one dataset
zwe.data <- merge(hiv.data, hhm.data, by = c("hivclust", "hivnumb", "hivline"), all.x =T, all.y =F)


# Clean data
zwe.data$weight <- ifelse(zwe.data$age <15, 
                          zwe.data$house_weight, zwe.data$hiv05) #See note 70 here https://dhsprogram.com/data/available-datasets.cfm  
zwe.data <- subset(zwe.data, hiv03 == 1 | hiv03 == 0)

# Bring in shapefile
library(rgdal)
zwemap <- readOGR(dsn = "C:/Users/steuben/Desktop/555/ZWE_adm_shp",layer = "ZWE_adm1")

#recode provinces to match the shapefile
library(car)
zwe.data$province <- recode(zwe.data$province, "1 = '3';2='4';
                            3 = '5'; 4 = '6'; 5 = '8'; 6 = '9';
                            7 = '10' ; 8 = '7'; 9 = '2'; 10 = '1'", as.numeric.result = T)

# Step 2 : Modeling and Creating maps

##################
##Naive estimates
##################
library(dplyr)
provinces <- as.numeric(unique(zwe.data$province))
props <- matrix(NA, nrow = length(provinces), ncol = 5)
props <- as.data.frame(props)
colnames(props)<- c("provcode", "p.hat","se.p.hat","y.i","n.i")
props$provcode <- provinces
for (i in props$provcode) {
  props[props$provcode == i, "p.hat"] <- mean(zwe.data$hiv03[zwe.data$province == i])
  props[props$provcode == i, "y.i"] <- sum(zwe.data$hiv03[zwe.data$province == i])
  props[props$provcode == i, "n.i"] <- length(zwe.data$hiv03[zwe.data$province == i])
  naivevar <- props[props$provcode == i, "p.hat"] * (1 -props[props$provcode == i, "p.hat"])/props[props$provcode == i, "n.i"]
  props[props$provcode == i, "se.p.hat"] <- sqrt(naivevar)
}
props <- arrange(props, provcode)

# Map sample sizes
library(SpatialEpi)
zweshapepoly <- SpatialPolygons(zwemap@polygons, proj4string = zwemap@proj4string)
summary(props[,"n.i"])
dev.new()
mapvariable(props[,"n.i"],zweshapepoly, ncut =1000, nlevels = 10, main = "Sample Size by Province")

#Mapping Naive estimates
dev.new()
mapvariable(props[,"p.hat"],zweshapepoly, main = "Naive Estimates")

#####################
#Naive binomial model
#####################

#create .graph file
library(spdep)
zwe.neigh <- poly2nb(zweshapepoly)
library(INLA)
nb2INLA("ZWE_adm_shp/zweNb.graph", zwe.neigh)

#Naive binomial model
props$unstruct <- props$struct <- 1:length(provinces)
formula = y.i ~1+f(struct, model = "besag",
                   constr = T, graph = "ZWE_adm_shp/zweNb.graph") +
                f(unstruct, model = "iid")

mod.smooth.unweighted <- inla(formula, family = "binomial", data = props,
                               Ntrials = n.i, control.predictor = list(compute = TRUE))# Create data set for INLA including design and weights
# Post medians of prevalences
psmoothunwt <- mod.smooth.unweighted$summary.fitted.values[,'0.5quant']
# Post standard deviations of prevalences
psmoothunwtsd <- mod.smooth.unweighted$summary.fitted.values[,'sd']
# Post medians of unstructured random effects
unwtunstruct <- mod.smooth.unweighted$summary.random$unstruct[,"0.5quant"]
# Post medians of spatial random effects
unwtstruct <- mod.smooth.unweighted$summary.random$struct[,"0.5quant"]

#Map unstructured and structured random effects
par(mar = c(1, 1, 1, 1))
mapvariable(unwtunstruct, zweshapepoly, ncut = 1000,
            nlevels = 10, main = "Unstrcutred Random Effects (Naive Binomial)")
par(mar = c(1, 1, 1, 1))
mapvariable(unwtstruct, zweshapepoly, ncut = 1000,
            nlevels = 10, main = "Spatial Random Effects (Naive Binomial)")


#Map Naive Binomial Prevalence
par(mar = c(1, 1, 1, 1))
mapvariable(psmoothunwt, zweshapepoly, main = "Naive Binomial Predicted Prevalence")

# Proportion of variation that is spatial
nareas <- 10
mat.marg <- matrix(NA, nrow = nareas, ncol = 1000)
m <- mod.smooth.unweighted$marginals.random$struct
for (i in 1:nareas) {
  Sre <- m[[i]]
  mat.marg[i, ] <- inla.rmarginal(1000, Sre)
}
var.Sre <- apply(mat.marg, 2, var)
var.eps <- inla.rmarginal(1000, inla.tmarginal(function(x) 1/x,
                                               mod.smooth.unweighted$marginals.hyper$"Precision for unstruct"))
mean(var.Sre)
mean(var.eps)
perc.var.Sre <- mean(var.Sre/(var.Sre + var.eps))

#Compare Naive and Binomial Naive
summary(props[,"p.hat"])
summary(psmoothunwt)
plot(psmoothunwt ~ props[, "p.hat"], pch = 19, xlim = c(0,
               0.15), ylim = c(0, 0.15), col = "blue", cex = 0.5,
               xlab = "Naive estimates", ylab = "Smoothed")
abline(0, 1, col = "green")

dev.new()
plot(psmoothunwtsd ~ props[, "se.p.hat"], pch = 19,
     xlim = c(0.0045, 0.006), ylim = c(0.0045, 0.006),
     col = "blue", cex = 0.5, xlab = "Naive Std Err",
     ylab = "Smoothed Post SD")
abline(0, 1, col = "green")

################################
# Incorporating Design Wieghts
################################

hist(zwe.data$weight, xlab = "Weights", main = "")

library(survey)
zwe.des <- svydesign(ids =~1, weights = ~weight, 
                     strata = ~strata_sample_design, data = zwe.data)
p.i <- svyby(~hiv03, ~province, zwe.des, svymean)$hiv03

dv.i <- svyby(~hiv03, ~province, zwe.des, svymean)$se^2

logit.pi <- log(p.i/(1-p.i))

v.i <- dv.i/(p.i^2 *(1-p.i)^2)

unwtvar <- props[, "se.p.hat"]^2
deff <- dv.i/unwtvar
effss <- props[, "n.i"]/deff
par(mfrow = c(1, 2))
hist(deff, main = "", xlab = "Design Effect")
plot(effss ~ props[, "n.i"], pch = 19, col = "blue",
     cex = 0.5, xlab = "Sample Size", ylab = "Effective Sample Size")
abline(0, 1, col = "green")

#Map weighted basic prevalence
par(mar = c(1,1,1,1), mfrow = c(1,1))
mapvariable(p.i, zweshapepoly, main = "Basic Weighted Prevalence")

par(mar = c(1,1,1,1), mfrow = c(1,1))
mapvariable(dv.i, zweshapepoly, main = "Design Variance")

#######################
#Adjusted Spatial Model
#######################

# Construct data frame for INLA
n.area = 10
data <- matrix(NA, nrow = n.area, ncol = 1)
data <- as.data.frame(data)
colnames(data)[1] <- "unstruct"
data$prov <- as.numeric(unique(zwe.data$province))
data$p.i <- p.i
data$dv.i <- dv.i
data$v.i <- v.i
data$logit.pi <- logit.pi
data$logit.prec <- 1/v.i
data <- data[order(data$prov), ]
data$unstruct <- 1:(n.area)
data$struct <- 1:(n.area)

#Fit model
formula = logit.pi ~ 1 + f(struct, model='besag', constr= T, graph='ZWE_adm_shp/zweNb.graph') + 
                        f(unstruct,model='iid')
mod.smooth <- inla(formula, family = "gaussian", data = data, 
                   control.predictor = list(compute = TRUE),
                   control.family = list( hyper = list(prec = list( initial = log(1), fixed=TRUE))), 
                   scale=logit.prec)

fixed.med <- rep(mod.smooth$summary.fixed[,4],dim(data)[1])
random.iid <- mod.smooth$summary.random$unstruct[,5]
random.smooth <- mod.smooth$summary.random$struct[,5]
linpred <- mod.smooth$summary.fitted.values[,'0.5quant']
pred <- exp(linpred)/(1+exp(linpred))
odds <- exp(linpred)
res <- cbind(data,fixed.med,random.iid,random.smooth,linpred,pred,odds)

# Create Maps
library(SpatialEpi)
par(mar=c(1,1,1,1))
mapvariable(res[,'pred'], zweshapepoly, main = "Prevlance - Adjusted Spatial Model")

par(mar = c(1, 1, 1, 1))
mapvariable(res[, "odds"], zweshapepoly, main = "HIV odds")

#Compare models
summary(p.i)
summary(res[, "pred"])
dev.new()
plot(res[, "pred"] ~ p.i, pch = 19, xlim = c(.07,.12), ylim = c(.07,.12), col = "blue", cex = 0.5, xlab = "Weighted", ylab = "Smoothed estimates")
abline(0, 1, col = "green")

#Post sd of prevalence
expit <- function(x) {
  exp(x)/(exp(x) + 1)
}
n.sim <- 1000
test <- matrix(NA, nrow = n.area, ncol = 2)
test <- as.data.frame(test)
colnames(test) <- c("simulated", "e.marginal")
for (i in 1:n.area) {
  test[i, "simulated"] <- sd(expit(inla.rmarginal(n.sim,
                                                  mod.smooth$marginals.linear.predictor[[i]])))
  expectations <- inla.emarginal(function(x) c(expit(x),
                                               expit(x)^2), mod.smooth$marginals.linear.predictor[[i]])
  test[i, "e.marginal"] <- sqrt(expectations[2] -
                                  expectations[1]^2)
}

plot(test$simulated, test$e.marginal, xlab = "Simulation",
     ylab = "Numerical Integration", pch = 19,
     cex = 0.5, col = "blue")
abline(a = 0, b = 1, col = "red")

plot(test$e.marginal ~ sqrt(dv.i), pch = 19, cex = 0.5,
     col = "blue", ylab = "Wtd Smooth Post S.D.", xlab = "Wtd estimate", xlim = c(.004,.012), ylim = c(.004,.012))
     abline(0, 1, col = "red")
     
     
#####################################
# Age Breakdowns
####################################
setwd("C:/Users/steuben/Desktop/555")
     
# split up by age groups
zwe.data.age1 <- zwe.data[zwe.data$age <= 14, ]
zwe.data.age2 <- subset(zwe.data, age >14 & age <= 24)
zwe.data.age3 <- subset(zwe.data, age >24 & age <= 34)
zwe.data.age4 <- subset(zwe.data, age >34 & age <=49)

pdf("plots.pdf")

###############
#age1
#hist(zwe.data.age1$weight, xlab = "Weights", main = "")

zwe.des <- svydesign(ids =~1, weights = ~weight, 
                     strata = ~strata_sample_design, data = zwe.data.age1)
p.i <- svyby(~hiv03, ~province, zwe.des, svymean)$hiv03

dv.i <- svyby(~hiv03, ~province, zwe.des, svymean)$se^2

logit.pi <- log(p.i/(1-p.i))

v.i <- dv.i/(p.i^2 *(1-p.i)^2)

unwtvar <- props[, "se.p.hat"]^2
deff <- dv.i/unwtvar
effss <- props[, "n.i"]/deff
#par(mfrow = c(1, 2))
#hist(deff, main = "", xlab = "Design Effect")
# plot(effss ~ props[, "n.i"], pch = 19, col = "blue",
#      cex = 0.5, xlab = "Sample Size", ylab = "Effective Sample Size")
# abline(0, 1, col = "green")

#Map weighted basic prevalence
par(mar = c(1,1,1,1), mfrow = c(1,1))
mapvariable(p.i, zweshapepoly, main = "Prevalence 0-14 years", lower = .01, upper = .4)

# Construct data frame for INLA
n.area = 10
data <- matrix(NA, nrow = n.area, ncol = 1)
data <- as.data.frame(data)
colnames(data)[1] <- "unstruct"
data$prov <- as.numeric(unique(zwe.data$province))
data$p.i <- p.i
data$dv.i <- dv.i
data$v.i <- v.i
data$logit.pi <- logit.pi
data$logit.prec <- 1/v.i
data <- data[order(data$prov), ]
data$unstruct <- 1:(n.area)
data$struct <- 1:(n.area)

#Fit model
formula = logit.pi ~ 1 + f(struct, model='besag', constr= T, graph='ZWE_adm_shp/zweNb.graph') + 
  f(unstruct,model='iid')
mod.smooth <- inla(formula, family = "gaussian", data = data, 
                   control.predictor = list(compute = TRUE),
                   control.family = list( hyper = list(prec = list( initial = log(1), fixed=TRUE))), 
                   scale=logit.prec)

fixed.med <- rep(mod.smooth$summary.fixed[,4],dim(data)[1])
random.iid <- mod.smooth$summary.random$unstruct[,5]
random.smooth <- mod.smooth$summary.random$struct[,5]
linpred <- mod.smooth$summary.fitted.values[,'0.5quant']
pred <- exp(linpred)/(1+exp(linpred))
odds <- exp(linpred)
res <- cbind(data,fixed.med,random.iid,random.smooth,linpred,pred,odds)

# Create Maps
par(mar=c(1,1,1,1))
mapvariable(res[,'pred'], zweshapepoly, main = "Prevlance - Adjusted Spatial Model (Age0-14)")


#Compare models
summary(p.i)
summary(res[, "pred"])
dev.new()
plot(res[, "pred"] ~ p.i, pch = 19, xlim = c(.02,.12), ylim = c(.02,.12), col = "blue", cex = 0.5, xlab = "Weighted", ylab = "Smoothed estimates")
abline(0, 1, col = "green")

#################
#age2
#hist(zwe.data.age2$weight, xlab = "Weights", main = "")
zwe.des <- svydesign(ids =~1, weights = ~weight, 
                     strata = ~strata_sample_design, data = zwe.data.age2)
p.i <- svyby(~hiv03, ~province, zwe.des, svymean)$hiv03

dv.i <- svyby(~hiv03, ~province, zwe.des, svymean)$se^2

logit.pi <- log(p.i/(1-p.i))

v.i <- dv.i/(p.i^2 *(1-p.i)^2)

unwtvar <- props[, "se.p.hat"]^2
deff <- dv.i/unwtvar
effss <- props[, "n.i"]/deff
#par(mfrow = c(1, 2))
#hist(deff, main = "", xlab = "Design Effect")
#plot(effss ~ props[, "n.i"], pch = 19, col = "blue",
#     cex = 0.5, xlab = "Sample Size", ylab = "Effective Sample Size")
#abline(0, 1, col = "green")

#Map weighted basic prevalence
par(mar = c(1,1,1,1), mfrow = c(1,1))
mapvariable(p.i, zweshapepoly, main = "Prevalence 15-24 years",lower = .01, upper = .4)

# Construct data frame for INLA
n.area = 10
data <- matrix(NA, nrow = n.area, ncol = 1)
data <- as.data.frame(data)
colnames(data)[1] <- "unstruct"
data$prov <- as.numeric(unique(zwe.data$province))
data$p.i <- p.i
data$dv.i <- dv.i
data$v.i <- v.i
data$logit.pi <- logit.pi
data$logit.prec <- 1/v.i
data <- data[order(data$prov), ]
data$unstruct <- 1:(n.area)
data$struct <- 1:(n.area)

#Fit model
formula = logit.pi ~ 1 + f(struct, model='besag', constr= T, graph='ZWE_adm_shp/zweNb.graph') + 
  f(unstruct,model='iid')
mod.smooth <- inla(formula, family = "gaussian", data = data, 
                   control.predictor = list(compute = TRUE),
                   control.family = list( hyper = list(prec = list( initial = log(1), fixed=TRUE))), 
                   scale=logit.prec)

fixed.med <- rep(mod.smooth$summary.fixed[,4],dim(data)[1])
random.iid <- mod.smooth$summary.random$unstruct[,5]
random.smooth <- mod.smooth$summary.random$struct[,5]
linpred <- mod.smooth$summary.fitted.values[,'0.5quant']
pred <- exp(linpred)/(1+exp(linpred))
odds <- exp(linpred)
res <- cbind(data,fixed.med,random.iid,random.smooth,linpred,pred,odds)

# Create Maps
par(mar=c(1,1,1,1))
mapvariable(res[,'pred'], zweshapepoly, main = "Prevlance - Adjusted Spatial Model (Age 15-24)")

#Compare models
summary(p.i)
summary(res[, "pred"])
dev.new()
plot(res[, "pred"] ~ p.i, pch = 19, xlim = c(.02,.12), ylim = c(.02,.12), col = "blue", cex = 0.5, xlab = "Weighted", ylab = "Smoothed estimates")
abline(0, 1, col = "green")

################
#age3
#hist(zwe.data.age3$weight, xlab = "Weights", main = "")

zwe.des <- svydesign(ids =~1, weights = ~weight, 
                     strata = ~strata_sample_design, data = zwe.data.age3)
p.i <- svyby(~hiv03, ~province, zwe.des, svymean)$hiv03

dv.i <- svyby(~hiv03, ~province, zwe.des, svymean)$se^2

logit.pi <- log(p.i/(1-p.i))

v.i <- dv.i/(p.i^2 *(1-p.i)^2)

unwtvar <- props[, "se.p.hat"]^2
deff <- dv.i/unwtvar
effss <- props[, "n.i"]/deff
#par(mfrow = c(1, 2))
#hist(deff, main = "", xlab = "Design Effect")
# plot(effss ~ props[, "n.i"], pch = 19, col = "blue",
#      cex = 0.5, xlab = "Sample Size", ylab = "Effective Sample Size")
# abline(0, 1, col = "green")

#Map weighted basic prevalence
par(mar = c(1,1,1,1), mfrow = c(1,1))
mapvariable(p.i, zweshapepoly, main = "Prevalence 25-34 years",lower = .01, upper = .4)

# Construct data frame for INLA
n.area = 10
data <- matrix(NA, nrow = n.area, ncol = 1)
data <- as.data.frame(data)
colnames(data)[1] <- "unstruct"
data$prov <- as.numeric(unique(zwe.data$province))
data$p.i <- p.i
data$dv.i <- dv.i
data$v.i <- v.i
data$logit.pi <- logit.pi
data$logit.prec <- 1/v.i
data <- data[order(data$prov), ]
data$unstruct <- 1:(n.area)
data$struct <- 1:(n.area)

#Fit model
formula = logit.pi ~ 1 + f(struct, model='besag', constr= T, graph='ZWE_adm_shp/zweNb.graph') + 
  f(unstruct,model='iid')
mod.smooth <- inla(formula, family = "gaussian", data = data, 
                   control.predictor = list(compute = TRUE),
                   control.family = list( hyper = list(prec = list( initial = log(1), fixed=TRUE))), 
                   scale=logit.prec)

fixed.med <- rep(mod.smooth$summary.fixed[,4],dim(data)[1])
random.iid <- mod.smooth$summary.random$unstruct[,5]
random.smooth <- mod.smooth$summary.random$struct[,5]
linpred <- mod.smooth$summary.fitted.values[,'0.5quant']
pred <- exp(linpred)/(1+exp(linpred))
odds <- exp(linpred)
res <- cbind(data,fixed.med,random.iid,random.smooth,linpred,pred,odds)

# Create Maps
par(mar=c(1,1,1,1))
mapvariable(res[,'pred'], zweshapepoly, main = "Prevlance - Adjusted Spatial Model (Age 25-34)")


#Compare models
summary(p.i)
summary(res[, "pred"])
dev.new()
plot(res[, "pred"] ~ p.i, pch = 19, xlim = c(.02,.12), ylim = c(.02,.12), col = "blue", cex = 0.5, xlab = "Weighted", ylab = "Smoothed estimates")
abline(0, 1, col = "green")

##############
#age4
#hist(zwe.data.age4$weight, xlab = "Weights", main = "")

zwe.des <- svydesign(ids =~1, weights = ~weight, 
                     strata = ~strata_sample_design, data = zwe.data.age4)
p.i <- svyby(~hiv03, ~province, zwe.des, svymean)$hiv03

dv.i <- svyby(~hiv03, ~province, zwe.des, svymean)$se^2

logit.pi <- log(p.i/(1-p.i))

v.i <- dv.i/(p.i^2 *(1-p.i)^2)

unwtvar <- props[, "se.p.hat"]^2
deff <- dv.i/unwtvar
effss <- props[, "n.i"]/deff
#par(mfrow = c(1, 2))
# hist(deff, main = "", xlab = "Design Effect")
# plot(effss ~ props[, "n.i"], pch = 19, col = "blue",
#      cex = 0.5, xlab = "Sample Size", ylab = "Effective Sample Size")
# abline(0, 1, col = "green")

#Map weighted basic prevalence
par(mar = c(1,1,1,1), mfrow = c(1,1))
mapvariable(p.i, zweshapepoly, main = "Prevalence 35-49 years",lower = .01, upper = .4)

# Construct data frame for INLA
n.area = 10
data <- matrix(NA, nrow = n.area, ncol = 1)
data <- as.data.frame(data)
colnames(data)[1] <- "unstruct"
data$prov <- as.numeric(unique(zwe.data$province))
data$p.i <- p.i
data$dv.i <- dv.i
data$v.i <- v.i
data$logit.pi <- logit.pi
data$logit.prec <- 1/v.i
data <- data[order(data$prov), ]
data$unstruct <- 1:(n.area)
data$struct <- 1:(n.area)

#Fit model
formula = logit.pi ~ 1 + f(struct, model='besag', constr= T, graph='ZWE_adm_shp/zweNb.graph') + 
  f(unstruct,model='iid')
mod.smooth <- inla(formula, family = "gaussian", data = data, 
                   control.predictor = list(compute = TRUE),
                   control.family = list( hyper = list(prec = list( initial = log(1), fixed=TRUE))), 
                   scale=logit.prec)

fixed.med <- rep(mod.smooth$summary.fixed[,4],dim(data)[1])
random.iid <- mod.smooth$summary.random$unstruct[,5]
random.smooth <- mod.smooth$summary.random$struct[,5]
linpred <- mod.smooth$summary.fitted.values[,'0.5quant']
pred <- exp(linpred)/(1+exp(linpred))
odds <- exp(linpred)
res <- cbind(data,fixed.med,random.iid,random.smooth,linpred,pred,odds)

# Create Maps

par(mar=c(1,1,1,1))
mapvariable(res[,'pred'], zweshapepoly, main = "Prevlance - Adjusted Spatial Model (Age 35 -49)")

#Compare models
summary(p.i)
summary(res[, "pred"])
dev.new()
plot(res[, "pred"] ~ p.i, pch = 19, xlim = c(.02,.12), ylim = c(.02,.12), col = "blue", cex = 0.5, xlab = "Weighted", ylab = "Smoothed estimates")
abline(0, 1, col = "green")

dev.off()