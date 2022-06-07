# packages
# 
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(nnet)

###############################################################################
# data 
# dat ------------------------------------------------------------
dat <- read.dta("https://stats.idre.ucla.edu/stat/data/ologit.dta")
#The data set has a dependent variable known as apply. It has 3 levels namely “unlikely”, “somewhat likely”, and “very likely”, coded in 1, 2, and 3 respectively. 3 being highest and 1 being lowest. This situation is best for using ordinal regression because of presence of ordered categories. Pared (0/1) refers to at least one parent has a graduate degree; public (0/1) refers to the type of undergraduate institute.
#
# ml -------------------------------------------------------------
ml <- read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")
#The data set contains variables on 200 students. The outcome variable is prog, program type. The predictor variables are social economic status, ses, a three-level categorical variable and writing score, write, a continuous variable.
#
# cases ----------------------------------------------------------
cases <-
  structure(list(Days = c(1L, 2L, 3L, 3L, 4L, 4L, 4L, 6L, 7L, 8L,
                          8L, 8L, 8L, 12L, 14L, 15L, 17L, 17L, 17L, 18L, 19L, 19L, 20L,
                          23L, 23L, 23L, 24L, 24L, 25L, 26L, 27L, 28L, 29L, 34L, 36L, 36L,
                          42L, 42L, 43L, 43L, 44L, 44L, 44L, 44L, 45L, 46L, 48L, 48L, 49L,
                          49L, 53L, 53L, 53L, 54L, 55L, 56L, 56L, 58L, 60L, 63L, 65L, 67L,
                          67L, 68L, 71L, 71L, 72L, 72L, 72L, 73L, 74L, 74L, 74L, 75L, 75L,
                          80L, 81L, 81L, 81L, 81L, 88L, 88L, 90L, 93L, 93L, 94L, 95L, 95L,
                          95L, 96L, 96L, 97L, 98L, 100L, 101L, 102L, 103L, 104L, 105L,
                          106L, 107L, 108L, 109L, 110L, 111L, 112L, 113L, 114L, 115L),
                 Students = c(6L, 8L, 12L, 9L, 3L, 3L, 11L, 5L, 7L, 3L, 8L,
                              4L, 6L, 8L, 3L, 6L, 3L, 2L, 2L, 6L, 3L, 7L, 7L, 2L, 2L, 8L,
                              3L, 6L, 5L, 7L, 6L, 4L, 4L, 3L, 3L, 5L, 3L, 3L, 3L, 5L, 3L,
                              5L, 6L, 3L, 3L, 3L, 3L, 2L, 3L, 1L, 3L, 3L, 5L, 4L, 4L, 3L,
                              5L, 4L, 3L, 5L, 3L, 4L, 2L, 3L, 3L, 1L, 3L, 2L, 5L, 4L, 3L,
                              0L, 3L, 3L, 4L, 0L, 3L, 3L, 4L, 0L, 2L, 2L, 1L, 1L, 2L, 0L,
                              2L, 1L, 1L, 0L, 0L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 0L, 0L,
                              0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)),
            .Names = c("Days", "Students"), 
            class = "data.frame", 
            row.names = c(NA, -109L))
#The data set consists of counts of high school students diagnosed with an infectious disease within a period of days from an initial outbreak.
#
# A ----------------------------------------------------------
A <- structure(list(numeracy = c(6.6, 7.1, 7.3, 7.5, 7.9, 7.9, 8,
                                 8.2, 8.3, 8.3, 8.4, 8.4, 8.6, 8.7, 8.8, 8.8, 9.1, 9.1, 9.1, 9.3,
                                 9.5, 9.8, 10.1, 10.5, 10.6, 10.6, 10.6, 10.7, 10.8, 11, 11.1,
                                 11.2, 11.3, 12, 12.3, 12.4, 12.8, 12.8, 12.9, 13.4, 13.5, 13.6,
                                 13.8, 14.2, 14.3, 14.5, 14.6, 15, 15.1, 15.7), 
                    anxiety = c(13.8, 14.6, 17.4, 14.9, 13.4, 13.5, 
                                13.8, 16.6, 13.5, 15.7, 13.6, 14, 16.1, 10.5, 16.9, 17.4, 13.9, 
                                15.8, 16.4, 14.7, 15, 13.3, 10.9, 12.4, 12.9, 16.6, 16.9, 15.4, 
                                13.1, 17.3, 13.1, 14, 17.7, 10.6, 14.7, 10.1, 11.6, 14.2, 12.1, 
                                13.9, 11.4, 15.1, 13, 11.3, 11.4, 10.4, 14.4, 11, 14, 13.4), 
                    success = c(0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 
                                0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 
                                0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L)),
               .Names = c("numeracy", "anxiety", "success"), 
               row.names = c(NA, -50L), 
               class = "data.frame")

# ------------ data cleaning --------------------------------------
# Load the raw training data and replace missing values with NA
training.data.raw <- read.csv('train.csv',header=T,na.strings=c(""))

# Output the number of missing values for each column
sapply(training.data.raw,function(x) sum(is.na(x)))

# Quick check for how many different values for each feature
sapply(training.data.raw, function(x) length(unique(x)))

# A visual way to check for missing data
library(Amelia)
missmap(training.data.raw, main = "Missing values vs observed")

# Subsetting the data
data <- subset(training.data.raw,select=c(2,3,5,6,7,8,10,12))

# Substitute the missing values with the average value
data$Age[is.na(data$Age)] <- mean(data$Age,na.rm=T)

# Check categorical variables encoding for better understanding of the fitted model
contrasts(data$Sex)
contrasts(data$Embarked)

# Remove rows (Embarked) with NAs
data <- data[!is.na(data$Embarked),]
rownames(data) <- NULL

# ------------ create training / test data -----------------------------
#partition and training data
library(caret)

split <- createDataPartition(y = train$Survived,p = 0.6,list = FALSE)

new_train <- train[split] 
new_test <- train[-split]

##################################################################
# analysis
# ----------- linear regression -----------------------------------
# Multiple Linear Regression Example
fit <- lm(y ~ x1 + x2 + x3, data=mydata)
summary(fit) # show results

# Other useful functions
coefficients(fit) # model coefficients
confint(fit, level=0.95) # CIs for model parameters
fitted(fit) # predicted values
residuals(fit) # residuals
anova(fit) # anova table
vcov(fit) # covariance matrix for model parameters
influence(fit) # regression diagnostics

# diagnostic plots
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(fit)

# compare models
fit1 <- lm(y ~ x1 + x2 + x3 + x4, data=mydata)
fit2 <- lm(y ~ x1 + x2)
anova(fit1, fit2)

# K-fold cross-validation
library(DAAG)
cv.lm(df=mydata, fit, m=3) # 3 fold cross-validation

# ===
# Assessing R2 shrinkage using 10-Fold Cross-Validation

fit <- lm(y~x1+x2+x3,data=mydata)

library(bootstrap)
# define functions
theta.fit <- function(x,y){lsfit(x,y)}
theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}

# matrix of predictors
X <- as.matrix(mydata[c("x1","x2","x3")])
# vector of predicted values
y <- as.matrix(mydata[c("y")])

results <- crossval(X,y,theta.fit,theta.predict,ngroup=10)
cor(y, fit$fitted.values)**2 # raw R2
cor(y,results$cv.fit)**2 # cross-validated R2
# ===

# Stepwise Regression
library(MASS)
fit <- lm(y~x1+x2+x3,data=mydata)
step <- stepAIC(fit, direction="both")
step$anova # display results


# All Subsets Regression ===
library(leaps)
attach(mydata)
leaps<-regsubsets(y~x1+x2+x3+x4,data=mydata,nbest=10)
# view results
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size
library(car)
subsets(leaps, statistic="rsq")
# ===


# Calculate Relative Importance for Each Predictor ===
library(relaimpo)
calc.relimp(fit,type=c("lmg","last","first","pratt"),
   rela=TRUE)

# Bootstrap Measures of Relative Importance (1000 samples)
boot <- boot.relimp(fit, b = 1000, type = c("lmg",
  "last", "first", "pratt"), rank = TRUE,
  diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result
# ===

# ----------- ordinal regression ----------------------------------
# # One of the assumptions underlying ordered logistic (and ordered probit) regression is that the relationship between each pair of outcome groups is the same.There is only one set of coefficients (only one model).

m <- polr(apply ~ pared + public + gpa, data = dat, Hess=TRUE)
summary(m) 
coef(summary(m))
coef(m) # log(confidence interval)
confint(m) # confidence interval
exp(coef(m)) # odds ratio

# ------------------- multinomial ------------------------------
# check data 
with(ml, table(ses, prog))

# check mean -- group by prog and calculate the mean and sd
with(ml, do.call(rbind, tapply(write, prog, function(x) c(M = mean(x), SD = sd(x)))))

# re-order prog and make "academic" first
ml$prog2 <- relevel(ml$prog, ref = "academic") 

#Multinomial logistic regression
test <- multinom(prog2 ~ ses + write, data = ml)

summary(test)
summary(test)$coefficients
z <- summary(test)$coefficients/summary(test)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2  # 2-tailed z test | p value

# odds ration
exp(coef(test))

#prediction?
#calculate predicted probabilities for each of the outcome levels using the fitted function
pp <- fitted(test)

#create small datasets varying one variable while holding the other constant
#holding write at its mean and examining the predicted probabilities for each level of ses.
dses <- data.frame(ses = c("low", "middle", "high"), write = mean(ml$write))
predict(test, newdata = dses, "probs")

#check averaged predicted probabilities for different values of the continuous predictor variable write within each level of ses
dwrite <- data.frame(ses = rep(c("low", "middle", "high"), each = 41), write = rep(c(30:70), 3))

# store the predicted probabilities for each value of ses and write
pp.write <- cbind(dwrite, predict(test, newdata = dwrite, type = "probs", se = TRUE))

# calculate the mean probabilities within each level of ses
by(pp.write[, 3:5], pp.write$ses, colMeans)

# plot
# melt data set to long for ggplot2
lpp <- melt(pp.write, id.vars = c("ses", "write"), value.name = "probability")

# plot predicted probabilities across write values for each level of ses
# facetted by program type
ggplot(lpp, aes(x = write, y = probability, colour = ses)) + 
  geom_line() + 
  facet_grid(variable ~ ., scales = "free")

# ------------- poisson --------------------------------------------------
# The mean and variance are different (actually, the variance is greater). 
# Now we plot the data.
attach(cases)
plot(Days, Students, xlab = "DAYS", ylab = "STUDENTS", pch = 16)
#
# poisson regression
model1 <- glm(Students ~ Days, poisson)
summary(model1)

# the residual deviance is greater than the degrees of freedom, so that we have over-dispersion. 
# This means that there is extra variance not accounted for by the model or by the error structure.
# Now let’s fit a quasi-Poisson model to the same data.
# 
model2 <- glm(Students ~ Days, quasipoisson) 
summary(model2) 
model2$coefficients

# prediction
# 
timeaxis <-seq (0,150,0.1)

Y <- predict(model2, list(Days = timeaxis))
plot(Days, Students, xlab = "DAYS", ylab = "STUDENTS", pch = 16)
lines(timeaxis, exp(Y), lwd = 2, col = "blue")

Z <- predict(model2, list(Days = timeaxis), type = "response")
plot(Days, Students, xlab = "DAYS", ylab = "NUMBER", pch = 16)
lines(timeaxis, Z, lwd = 2, col = "red")

# calculate the impact on the number of cases arising from a one day increase along the time axis. by take the exponential of the coefficients.

coeffs <- exp(coef(model2))

# calculate the 95% confidence interval (upper and lower confidence limits) as follows:

CI <- exp(confint.default(model2))

# the change in number of students presenting with the disease for each additional day:
1 - 0.9826884
# [1] 0.0173116
# The reduction (rate ratio) is approximately 0.02 cases for each additional day.
# 
# ------------ glm logistic ----------------------------------
# 
attach(A)
names(A)

model1 <- glm(success ~ numeracy * anxiety, binomial)
# After the ~, we list the two predictor variables. The * indicates that not only do we want each main effect, but we also want an interaction term between numeracy and anxiety.
# 
summary(model1)

# plotting fits
model_numeracy <- glm(success ~ numeracy, binomial)
summary(model_numeracy)

model_anxiety <- glm(success ~ anxiety, binomial)
summary(model_anxiety)

#find the range of each variable.

range(numeracy)
# [1] 6.6 15.7

range(anxiety)
# [1] 10.1 17.7
# 
# setup x - numveracy 
xnumeracy <- seq(0, 15, 0.01)

# use the predict() function to set up the fitted values. The syntax type = “response” back-transforms from a linear logit model to the original scale of the observed data (i.e. binary).
ynumeracy <- predict(model_numeracy, list(numeracy = xnumeracy),type = "response")

# plotting 
plot(numeracy, success, pch = 16, xlab = "NUMERACY SCORE", ylab = "ADMISSION")
lines(xnumeracy, ynumeracy, col = "red", lwd = 2)

# same with anxiety
xanxiety <- seq(10, 20, 0.1)
yanxiety <- predict(model_anxiety, list(anxiety = xanxiety),type = "response")

plot(anxiety, success, pch = 16, xlab = "ANXIETY SCORE", ylab = "SUCCESS")
lines(xanxiety, yanxiety, col = "blue", lwd = 2)
# Clearly, those who score high on anxiety are unlikely to be admitted, possibly because their admissions test results are affected by their high level of anxiety.
# 
# 
# ------------ logistic regression ---------------------------------
# Model fitting
model <- glm(Survived ~.,family=binomial(link='logit'),data=train)
summary(model)

# Analysis of deviance
anova(model,test="Chisq")

# McFadden R^2
library(pscl)
pR2(model)

# prediction --------
# If prob > 0.5 then 1, else 0. Threshold can be set for better results
fitted.results <- predict(model,newdata=subset(test,select=c(2,3,4,5,6,7,8)),type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

# compare result ----
misClasificError <- mean(fitted.results != test$Survived)
print(paste('Accuracy',1-misClasificError))

# Confusion matrix
library(caret)
confusionMatrix(data=fitted.results, reference=test$Survived)


# ----------- ROC AUC -----------------------------------------------
# library(ROCR)
# ROC and AUC
p <- predict(model, newdata=subset(test,select=c(2,3,4,5,6,7,8)), type="response")
pr <- prediction(p, test$Survived)
# TPR = sensitivity, FPR=specificity
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

# ------------ Titanic data ---------------------------------------
# model training
log_model <- glm(Survived ~ Pclass + Sex + Age +  SibSp + Fare + title, data = new_train[,-c("PassengerId","Name","Ticket")],family = binomial(link="logit"))
summary(log_model)
log_predict <- predict(log_model,newdata = new_test,type = "response")
log_predict <- ifelse(log_predict > 0.5,1,0)

#plot ROC
library(ROCR)
library(Metrics)
pr <- prediction(log_predict,new_test$Survived)
perf <- performance(pr,measure = "tpr",x.measure = "fpr")
plot(perf)
auc(new_test$Survived,log_predict) #0.78

#plot ROC with threshold value
log_predict <- predict(log_model,newdata = new_test,type = "response")
log_predict <- ifelse(log_predict > 0.6,1,0)
pr <- prediction(log_predict,new_test$Survived)
perf <- performance(pr,measure = "tpr",x.measure = "fpr")
plot(perf)
auc(new_test$Survived,log_predict) #0.8008







