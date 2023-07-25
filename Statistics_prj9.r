library(dplyr)
library(reshape2)
# install.packages("dlookr")
library(dlookr)
library(arm)
library(jtools)
library(interactions)
install.packages("broom.mixed")
library(broom.mixed)

HorseKicks <- read.delim("HorseKicks.txt", row.names=NULL, stringsAsFactors=FALSE)
str(HorseKicks)

# transform Year column data type to factor 
HorseKicks <- mutate(HorseKicks, Year = as.factor(HorseKicks$Year)) 

# transform to short df with rowname = 1st column
HorseKicks.s <- data.frame(HorseKicks[,c(2:15)], row.names = HorseKicks$Year)

# count soldiers killed by horse kicks by corp membership and by year
# yearsum 
HorseKicks.s[,c(15:16)] <- NULL
yearsum <- apply(HorseKicks.s, 1, sum) 
yearsum

# corpsum
corpsum <- apply(HorseKicks.s, 2, sum)
corpsum

# average number of death per year
year.average <- apply(HorseKicks.s, 1, mean)
year.average

# plot the data and show your fit
# transform data to long df
HorseKicks2 <- melt(HorseKicks, id.vars = "Year", 
                    variable.name = "corps", value.name = "knumber")
head(HorseKicks2)

# plot the data
hist(HorseKicks2$knumber)

# Shapiro-Wilk normality test to show non-normal distribution
shapiro.test(HorseKicks2$knumber)


# fit the model
hk.pos.model <- glm(knumber ~ Year + corps, data = HorseKicks2, family = poisson(link = "log"))
summary(hk.pos.model)

# extract coefficients and their std.error
coef.hk <- coef(hk.pos.model)
se.coef.hk <- se.coef(hk.pos.model)

# col-combine to show model 
model.hk <- cbind(coef.hk, se.coef.hk)
model.hk

# plot regression coefficients
plot_summs(hk.pos.model, scale = TRUE, exp = TRUE)



