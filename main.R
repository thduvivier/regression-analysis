library(MASS)
library(glmnet)
library(robustbase)

###############.
# HOSPITAL DATA----
###############.

# Import data and packages
setwd('C:/Users/duviv/Documents/University/KUL/S1/G0S75A - Regression analysis/Project')
hospital.full <- read.table("hospital.txt", header=TRUE)

# Split data
set.seed(0733342)
h.test <- sample(1:dim(hospital.full)[1], 20 )
hospital.test <- hospital.full[h.test, ]
hospital <- hospital.full[-h.test, ]
attach(hospital)

cor(hospital)

##################.
# Model Selection----
##################.

# Backward elimination based on AIC
fit.full <- lm(irisk ~ ., data = hospital)
stepAIC(fit.full, list(upper = ~ lstay + age + cultratio + xrayratio + nbeds + census + nnurse + facil, lower = ~ 1), direction = "backward")
# AIC: -9.98
# irisk ~ lstay + cultratio + xrayratio + facil

# Forward selection based on AIC
fit.null <- lm(irisk ~ 1, data = hospital)
fit.null
stepAIC(fit.null, scope = list(upper = ~ lstay + age + cultratio + xrayratio + nbeds + census + nnurse + facil, lower = ~ 1), direction = "forward")
# AIC = -9.98
# irisk ~ lstay + cultratio + xrayratio + facil

# Stepwise selection based on AIC (started at full model)
stepAIC(fit.full, scope=list(upper = ~ lstay + age + cultratio + xrayratio + nbeds + census + nnurse + facil, lower = ~ 1), direction = "both")
# AIC = -9.98
# irisk ~ lstay + cultratio + xrayratio + facil

# Stepwise selection based on AIC (started at null model)
stepAIC(fit.null, scope=list(upper = ~ lstay + age + cultratio + xrayratio + nbeds + census + nnurse + facil, lower = ~ 1), direction = "both")
# AIC = -9.98
# irisk ~ lstay + cultratio + xrayratio + facil

# Backward elimination based on F-statistic/t-statistic
dropterm(fit.full, test = "F") # if I drop this term, AIC will become ...
# Remove nnurse
fit1 <- update(fit.full, ~ . - nnurse)
dropterm(fit1, test = "F")
# Remove census
fit2 <- update(fit1, ~ . - census)
dropterm(fit2, test = "F")
# Remove nbeds
fit3 <- update(fit2, ~ . - nbeds)
dropterm(fit3, test = "F")
# Remove age
fit4 <- update(fit3, ~ . - age)
dropterm(fit4, test = "F") # if none added, have the lowest AIC so stop here
# irisk ~ lstay + cultratio + xrayratio + facil

# Forward selection based on F-statistic/t-statistic
addterm(fit.null, ~ . + lstay + age + cultratio + xrayratio + nbeds + census + nnurse + facil, test = "F")
# Add cultratio
fit1 <- update(fit.null, ~ . + cultratio)
addterm(fit1, ~ . + lstay + age + xrayratio + nbeds + census + nnurse + facil, test = "F")
# Add lstay
fit2 <- update(fit1, ~ . + lstay)
addterm(fit2, ~. + age + xrayratio + nbeds + census + nnurse + facil, test = "F")
# Add facil
fit3 <- update(fit2, ~ . + facil)
addterm(fit3, ~. + age + xrayratio + nbeds + census + nnurse, test = "F")
# Add xrayratio
fit4 <- update(fit3, ~ . + xrayratio)
addterm(fit4, ~. + age + nbeds + census + nnurse, test = "F")
# Stop here with AIC = -9.98
# irisk ~ lstay + cultratio + xrayratio + facil

# All methods used gave the same model which is 
# irisk ~ lstay + cultratio + xrayratio + facil
# Therefore we will accept it as the best model
# Could make model validation if we had alternative models. Will be checked with models with interaction/quadratic

fit1 <- lm(irisk ~ lstay + cultratio + xrayratio + facil, data = hospital)
summary(fit1)


#########################.
# Gauss Markov Conditions----
#########################.

  # Graphical
fit1.res <- residuals(fit1)
fit1.stdres <- stdres(fit1)
fit1.fittedvalues <- fitted.values(fit1)
par(mfrow = c(1,1))
    # Normality
qqnorm(fit1.stdres)
qqline(fit1.stdres)
    #uncorrelated errors
plot(fit1.res, xlab = "Index", ylab = "Residual")
    #Homoscedasticity
plot(fit1.fittedvalues, fit1.res, xlab = "Fitted value", ylab = "Residual")
lines(lowess(fit1.res ~ fit1.fittedvalues), col = "red")
    #detect influential obs and outliers
plot(fit1.stdres, xlab = "Index", ylab = "Standardized residual", ylim = c(-3,3))
abline(h = -2.5, lty = 2)
abline(h = 2.5, lty = 2)

#################.
#Partial Residual plots to detect curvature----
#################.

fit1.coef <- coefficients(fit1)
fit1.pres.lstay <- fit1.res + fit1.coef[2] * lstay
fit1.pres.cultratio <- fit1.res + fit1.coef[3] * cultratio
fit1.pres.xrayratio <- fit1.res + fit1.coef[4] * xrayratio
fit1.pres.facil <- fit1.res + fit1.coef[5] * facil

par(mfrow = c(2,2))
plot(lstay, fit1.pres.lstay, ylab = "Partial residuals (lstay)")
abline(lm(unname(fit1.pres.lstay) ~ lstay))
lines(lowess(lstay, fit1.pres.lstay), col = "red")

plot(cultratio,fit1.pres.cultratio, ylab = "Partial residuals (cultratio)")
abline(lm(fit1.pres.cultratio ~ cultratio))
lines(lowess(cultratio, fit1.pres.cultratio), col = "red")

plot(xrayratio,fit1.pres.xrayratio, ylab = "Partial residuals (xrayratio)")
abline(lm(fit1.pres.xrayratio ~ xrayratio))
lines(lowess(xrayratio, fit1.pres.xrayratio), col = "red")

plot(facil,fit1.pres.facil, ylab = "Partial residuals (facil)")
abline(lm(fit1.pres.facil ~ facil))
lines(lowess(facil, fit1.pres.facil), col = "red")

#################.
# Transformations----
#################.

# Quadratic term for facil
#####.

fit2 <- lm(irisk ~ lstay + cultratio + facil + I(facil^2) + xrayratio , data = hospital)
summary(fit2)
sum(residuals(fit2)^2)

fit3 <- lm(irisk ~ lstay + cultratio + facil + I(facil^2), data = hospital)
summary(fit3)
sum(residuals(fit3)^2)

#partial residual for quadratic facil

fit3.res <- residuals(fit3)
fit3.stdres <- stdres(fit3)
fit3.fittedvalues <- fitted.values(fit3)
fit3.coef <- coefficients(fit3)
fit3.pres.facil <- fit3.res + fit3.coef[4] * facil + fit3.coef[5] * (facil^2)
fit3.pred <- predict(lm(fit3.pres.facil ~ facil + I(facil^2)))

par(mfrow = c(1,1))
plot(facil, fit3.pres.facil, ylab = "Partial residuals (facil)")
points(facil, fit3.pred, col= "red")
lines(lowess(facil, fit3.pres.facil), col = "red")

#residual vs fitted values for the whole model
plot(fit3.fittedvalues, fit3.res, xlab = "Fitted value", ylab = "Residual")
lines(lowess(fit3.res ~ fit3.fittedvalues), col = "red")

# Logarithmic term for cultratio
.
fit4 <- lm(irisk ~ lstay + cultratio + log(cultratio) + facil + xrayratio, data = hospital)
summary(fit4)
sum(residuals(fit4)^2)

fit5 <- lm(irisk ~ lstay + log(cultratio) + facil, data = hospital)
summary(fit5)
sum(residuals(fit5)^2)

fit5.res <- residuals(fit5)
fit5.stdres <- stdres(fit5)
fit5.fittedvalues <- fitted.values(fit5)
fit5.coef <- coefficients(fit5)

fit5.pres.cultratio <- fit5.res + fit5.coef[3] * log(cultratio)
fit5.pred <- predict(lm(fit5.pres.cultratio ~ log(cultratio)))

#partial residual for logarithmic cultratio
par(mfrow = c(1,1))
plot(cultratio, fit5.pres.facil, ylab = "Partial residuals (facil)")
points(cultratio, fit5.pred, col= "red")
lines(lowess(cultratio, fit5.pres.facil), col = "red")

#residual vs fitted values for the whole model
plot(fit5.fittedvalues, fit5.res, xlab = "Fitted value", ylab = "Residual")
lines(lowess(fit5.res ~ fit5.fittedvalues), col = "red")

# Weighted least squares
.

stdev <- lm(abs(residuals(fit1)) ~ lstay + cultratio + xrayratio + facil)
summary(stdev)
weight.y <- 1/stdev$fitted^2
fit6 <- lm(irisk ~ lstay + cultratio + xrayratio + facil, weights = weight.y)
summary(fit6)
plot(fitted.values(fit1), residuals(fit6)*sqrt(weight.y), ylab = "Weighted residuals")

###################.
# Model validation----
###################.

# PRESS
PRESS1 <- sum((residuals(fit1) / (1 - lm.influence(fit1)$hat))^2)
PRESS2 <- sum((residuals(fit3) / (1 - lm.influence(fit3)$hat))^2)
PRESS3 <- sum((residuals(fit5) / (1 - lm.influence(fit5)$hat))^2)
PRESS4 <- sum((residuals(fit6) / (1 - lm.influence(fit6)$hat))^2)

PRESS <- c(PRESS1, PRESS2, PRESS3, PRESS4)
names(PRESS) <- c("model1", "model2", "model3", "model4")
sort(PRESS)

# MSE
MSE1 <- summary(fit1)$sigma^2
MSE2 <- summary(fit3)$sigma^2
MSE3 <- summary(fit5)$sigma^2
MSE4 <- summary(fit6)$sigma^2
MSE <- c(MSE1, MSE2, MSE3, MSE4)
names(MSE) <- c("model1", "model2", "model3", "model4")
sort(MSE)

# MSEP
MSEP1 <- mean((predict(fit1, newdata = hospital.test) - hospital.test$irisk)^2)
MSEP2 <- mean((predict(fit3, newdata = hospital.test) - hospital.test$irisk)^2)
MSEP3 <- mean((predict(fit5, newdata = hospital.test) - hospital.test$irisk)^2)
MSEP4 <- mean((predict(fit6, newdata = hospital.test) - hospital.test$irisk)^2)
MSEP <- c(MSEP1, MSEP2, MSEP3,MSEP4)
names(MSEP) <- c("model1", "model2", "model3", "model4")
sort(MSEP)

# Results
validation.results <- data.frame(rbind(PRESS, MSE, MSEP), row.names = c("PRESS", "MSE", "MSEP"))
names(validation.results) <- c("model1", "model2", "model3", "model4")
validation.results

##################.
# Ridge Regression----
##################.

# With selected variables
############.

hospital.std <- scale(hospital)
hospital.std <- as.data.frame(hospital.std)
hospital.std <- hospital.std[,c(1,3,4,5,9)]
hospital.std <- hospital.std[,c(2,1,3,4,5)]

y1 <- hospital.std$irisk
f <- as.formula(y1 ~ .*.) #put all the interactions in the model
x1 <- model.matrix(f, hospital.std[,c(2:5)])[, -1]

lambdas <- 10^seq(3, -2, by = -.1) #so glmnet runs the model with many different lambdas
ridge.fit <- glmnet(x = x1, y = y1, alpha = 0, lambda = lambdas)
summary(ridge.fit)

cv_fit <- cv.glmnet(x1, y1, alpha = 0, lambda = lambdas)
opt_lambda <- cv_fit$lambda.min
opt_lambda

best.ridge <- cv_fit$glmnet.fit
summary(best.ridge)

#predict
hospital.test.std <- scale(hospital.test)
hospital.test.std <- as.data.frame(hospital.test.std)
hospital.test.std <- hospital.test.std[,c(1,3,4,5,9)]
hospital.test.std <- hospital.test.std[,c(2,1,3,4,5)]

y2 <- hospital.test.std$irisk
f <- as.formula(y2 ~ .*.) #put all the interactions in the model
x2 <- model.matrix(f, hospital.test.std[,c(2:5)])[, -1]

ridge.predict <- predict(best.ridge, s = opt_lambda, newx = x2, newy = y2)

sst1 <- sum((y2 - mean(y2))^2)
sse1 <- sum((ridge.predict - y2)^2)

# R squared
rsq1 <- 1 - sse1 / sst1
rsq1

# With all variables and interaction
############.

hospital.std.full <- scale(hospital)
hospital.std.full <- as.data.frame(hospital.std.full)

f2 <- as.formula(y1 ~ .*.) #put all the interactions in the model
x2 <- model.matrix(f2, as.data.frame(hospital.std.full[,c(1:2, 4:9)]))[, -1]

lambdas <- 10^seq(3, -2, by = -.1) #so glmnet runs the model with many different lambdas
ridge.fit.full <- glmnet(x = x2, y = y1, alpha = 0, lambda = lambdas)
summary(ridge.fit.full)

cv_fit.full <- cv.glmnet(x2, y1, alpha = 0, lambda = lambdas)
opt_lambda.full <- cv_fit.full$lambda.min
opt_lambda.full

best.ridge.full <- cv_fit.full$glmnet.fit
summary(best.ridge.full)

#predict
hospital.test.std.full <- as.data.frame(scale(hospital.test))

x3 <- model.matrix(f, as.data.frame(hospital.test.std.full[,c(1:2, 4:9)]))[, -1]

ridge.predict.full <- predict(best.ridge.full, s = opt_lambda.full, newx = x3, newy = y2)

sst2 <- sum((y2 - mean(y2))^2)
sse2 <- sum((ridge.predict.full - y2)^2)

rsq2 <- 1 - sse2 / sst2
rsq2



###############.
# FOSSIL DATA----
###############.

setwd('C:/Users/duviv/Documents/University/KUL/S1/G0S75A - Regression analysis/Project')
fossil <- read.table("fossil.txt", header=TRUE)
attach(fossil)

plot(strontium.ratio ~ age)

############.
# Parametric----
############.

# Linear model
fit1 <- lm(strontium.ratio ~ age)
summary(fit1)
abline(fit1, col = "red")

# Quadratic model
fit2 <- lm(strontium.ratio ~ age + I(age^2))
summary(fit2)
fit2.coef <- fit2$coefficients
curve(fit2.coef[1] + fit2.coef[2]*x + fit2.coef[3]*x^2, add = TRUE, col = "green")

# Cubic model
fit3 <- lm(strontium.ratio ~ age + I(age^2) + I(age^3))
summary(fit3)
fit3.coef <- fit3$coefficients
curve(fit3.coef[1] + fit3.coef[2]*x + fit3.coef[3]*x^2 + fit3.coef[4]*x^3, add = TRUE, col = "blue")

############.
# Non-Parametric----
############.

# Local linear regression----
########.

plot(age, strontium.ratio, main = "Local linear regression")
s <- c(2/3, 1/3, 1/10)
colors <- c("red", "green", "blue")
for (i in 1:length(s)) lines(age, predict(loess(strontium.ratio ~ age, span = s[i], degree = 1)), col = colors[i])
legend(95, 0.7073, c("span = 2/3", "span = 1/3", "span = 1/10"), lty = 1, col = colors)

# Local quadratic regression----
########.

plot(age, strontium.ratio, main = "Local quadratic regression")
for (i in 1:length(s)) lines(age, predict(loess(strontium.ratio ~ age, span = s[i], degree = 2)), col = colors[i])
legend(95, 0.7073, c("span = 2/3", "span = 1/3", "span = 1/10"), lty = 1, col = colors)

# Check assumptions----
########.

# Check model assumptions
fit4 <- loess(strontium.ratio ~ age, span = 2/3, degree = 2)
summary(fit4)
par(mfrow=c(2,2))
qqnorm(residuals(fit4))
qqline(residuals(fit4))
scatter.smooth(residuals(fit4), span = 1, degree = 1)
scatter.smooth(fitted(fit4), sqrt(abs(residuals(fit4))), span = 1, degree = 1)
# constant variance is questionable (region with times = 1:14 has small variance, region with times 15:57 has large variance)


###############.
# FUEL DATA----
###############.

# Import data and packages
setwd('C:/Users/thoma/Documents/University/KUL/S1/G0S75A - Regression analysis/Project')
fuel <- read.table("fuel2001.txt", header=TRUE)
attach(fuel)

fit1 <- lm(Fuel ~ Income + Miles + Tax + Dlic)
summary(fit1)
resid1 <- stdres(fit1)

plot(x = Fuel, y = resid1)
abline(h = -2.5, lty = 2)
abline(h = 2.5, lty = 2) # Find two outliers
text(Fuel, resid1, row.names(fuel), cex=0.6, pos=4, col="red")
# We see we need to use Robust LTS

##############.
# Robust LTS----
##############.

# RLTS (50% breakdown value)
RLTS <- ltsReg(Fuel ~ Income + Miles + Tax + Dlic)
summary(RLTS)

# Detection of outliers
plot(RLTS, which = "rindex")
plot(RLTS, which = "rdiag")

# Standardized residuals
RLTS.stdres <- RLTS$residuals/RLTS$scale
plot(RLTS.stdres, ylim = c(-12,12), ylab = "Standardized residuals")
abline(h = c(-2.5,2.5), col = "red")

# Diagnostic plot
p <- length(fuel)
plot(RLTS$RD, RLTS.stdres, ylim = c(-12,12), xlab = "Robust distances", ylab = "Standardized residuals")
abline(v = sqrt(qchisq(0.975, p - 1)), col = "red")
abline(h = c(-2.5,2.5), col = "red")
text(RLTS$RD, RLTS.stdres, row.names(fuel), cex=0.6, pos=4, col="red")

