library(ggplot2)
library(magrittr)
library(reshape2)
library(GGally)

# use the data set in weight-height.csv
wh <- read.csv("weight-height.csv")
class(wh); str(wh); head(wh); summary(wh)

# Question1. plot the distribution of heights in the group; is the distribution Gaussian?
ggplot(wh, aes(x = Height))+
  geom_histogram(binwidth = 1, color = "black", fill = "lightblue")+
  geom_freqpoly(binwidth = 1, color = "red")

#access normality
qqnorm(wh$Height)
qqline(wh$Height, col = 3, lwd = 2)


# method 1: round the numbers to integer
h.table <- table(as.integer(wh$Height))
plot(h.table)
plot(h.table, type = "l")


# method 2: histogram
hist(wh$Height)


# method 3: plot density of height
plot(density(wh$Height))


# method 4: using binning function and then plot the data
mn <- floor(min(wh$Height))
mx <- ceiling(max(wh$Height))
bin <- cut(wh$Height, breaks = seq(mn,mx,1))
plot(bin)

# Question 4: plot the distribution of heights of men in black and women in red on the same graph and for the whole group (blue)

wh_male <- wh[which(wh[,1] == "Male"),]
wh_female <- wh[which(wh[,1] == "Female"),]

#using R basic plot
plot(density(wh$Height), col = "blue", ylim = c(0, 0.15))
lines(density(wh_male$Height), col = "black", ylim = c(0, 0.3))
lines(density(wh_female$Height), col = "red", ylim = c(0, 0.3))

#using ggplot
ggplot(data = wh, mapping = aes(x = Height, color = Gender))+
  geom_histogram(binwidth = 1 )+
  geom_freqpoly(binwidth = 1)+
  geom_freqpoly(aes(x = Height), binwidth = 1, color = "blue")
  
  
#====================#
# Question5. repeat 1-4 with weight instead of height
#====================#

# 1. plot the distribution of weights in the group; is the distribution Gaussian?
ggplot(wh, aes(x = Weight))+
  geom_histogram(binwidth = 10, color = "black", fill = "lightblue")+
  geom_freqpoly(binwidth = 10, color = "red")

#access normality
qqnorm(wh$Weight)
qqline(wh$Weight, col = 3, lwd = 2)


# method 1: round the numbers to integer and plot table
w.table <- table(as.integer(wh$Weight))
plot(w.table)
plot(w.table, type = "l")


# method 2: histogram
hist(wh$Weight)

# method 3: plot density of weight
plot(density(wh$Weight))


# method 4: using binning function and then plot the data
mn <- floor(min(wh$Weight))
mx <- ceiling(max(wh$Weight))
bin <- cut(wh$Weight, breaks = seq(mn,mx,10))
plot(bin)


# Question 4: plot the distribution of weights of men in black and women in red on the same graph and for the whole group (blue)

plot(density(wh$Weight), col = "blue", ylim = c(0, 0.023))
lines(density(wh_male$Weight), col = "black", ylim = c(0, 0.023))
lines(density(wh_female$Weight), col = "red", ylim = c(0, 0.023))


ggplot(data = wh, mapping = aes(x = Weight, color = Gender))+
  geom_histogram(binwidth = 10 )+
  geom_freqpoly(binwidth = 10)+
  geom_freqpoly(aes(x = Weight), binwidth = 10, color = "blue")

# Question6. plot height (x-axis) versus weight (y-axis) for the group
plot(x = wh$Height, y = wh$Weight, cex = 0.1, col = "blue")

# Question7. plot height versus weight with man in black and women in red
plot(x = wh$Height, y = wh$Weight, cex = 1, col = "blue")
points(x = wh_male$Height, y = wh_male$Weight, cex = 0.1, col = "black")
points(x = wh_female$Height, y = wh_female$Weight, cex = 0.1, col = "red")

abline(lm(data = wh, Weight ~ Height), lwd = 3, lty = 1, col = "blue")
abline(lm(data = wh_male, Weight ~ Height), lwd = 3, lty = 2, col = "black")
abline(lm(data = wh_female, Weight ~ Height), lwd = 3, lty = 3, col = "red")

# Question8. a) is height a good proxy for weight?
summary(lm(data = wh, Weight ~ Height))

# b) which rate of increase in weight with height is higher, men or women? 
lm(data = wh_male, Weight ~ Height)
lm(data = wh_female, Weight ~ Height)


# c) for which group (men versus women), is height a better proxy for weight?
summary(lm(data = wh_male, Weight ~ Height))
summary(lm(data = wh_female, Weight ~ Height))

# d) better picture is seen using 3-d plot, but not binning data
mn <- floor(min(wh$Weight))
mx <- ceiling(max(wh$Weight))
bin <- cut(wh$Weight, breaks = seq(mn,mx,5))
plot(bin, cex = 1, col = "blue")
abline(lm(data = wh, Weight ~ Height), lwd = 3, lty = 1, col = "blue")

ggplot(data = wh, mapping = aes(x = Height, y = Weight, color = Gender))+
  geom_point(alpha = 0.05, size = 0.5)+
  geom_smooth(method = "lm", se = T )+
  geom_point( aes(x = Height, y = Weight), color = "black", alpha = 0.03, size = 0.3)+
  geom_smooth(aes(x = Height, y = Weight), color = "black", method = "lm", se = T)
  
# Question9 use dataset iris
class(iris); dim(iris); str(iris)

# describe the data through two interesting plot
ggpairs(data = iris)
ggplot(data = iris, mapping = aes(x = Petal.Length, color = Species))+
  geom_histogram(bins = 60)+
  geom_density()

ggplot(data = iris, mapping = aes(x = Petal.Length, y = Petal.Width, color = Species))+
  geom_point()+
  geom_smooth(method = "lm", se = T)

# use the reshape2 to show multiple panels in one plot
iris.long <- melt(iris, id = "Species")
ggplot(data = iris.long, mapping = aes(x = variable, y = value, color = Species))+
  geom_boxplot()

