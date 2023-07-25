rm(list = ls())
graphics.off()

library(dplyr)
wh <- read.csv("weight-height.csv")
swh <- sample_n(wh, 10)

summary(swh)
?dbinom
dbinom(5, 10, prob = .5)

# find out the sample size at the point error less than 5%

set.seed(10000)
n <- 10
swh <- sample_n(wh, 10)
p <- dbinom(dim(swh[swh$Gender == "Male", ])[1], 10, prob = .5)
while(p > 0.05){
  n <- n + 10
  swh <- sample_n(wh, n)
  p <- dbinom(dim(swh[swh$Gender == "Male", ])[1], n, prob = .5)
}

result <- paste("when n = ", n, " the error is less than 5%.") 
result
  
# binom.test
binom.test(5, 18, p = 0.3)
binom.test(6, 18, p = 0.3)
binom.test(4, 18, p = 0.3)


# kolmogorov-smirmov (K-S) test
x <- rnorm(1e5, 1, 2)
ks.test(x, "pnorm")

#==wilcoxon test==#
x <- rnorm(50)
y <- runif(30)
ks.test(x,y)

x2 <- rnorm(50, -1)

plot(ecdf(x))
plot(ecdf(x2), add = T)

t.test(x, x2)
wilcox.test(x, x2)
ks.test (x, x2, alternative = "less")

##== wilcoxon test II ==##

women_weight <- c (38.9, 61.2, 73.3, 21.8, 63.4, 64.6, 48.4, 48.8, 48.5)
men_weight <- c(67.8, 60, 63.4, 76, 89.4, 73.3, 67.3, 61.3, 62.4)

# create a data frame
my_data <- data.frame (group = rep(c("Women", "Man"), each = 9),
                       weight = c(women_weight, men_weight))
boxplot(weight ~ group, data = my_data)

res <- wilcox.test(women_weight, men_weight, exact = F)
res <- wilcox.test(weight ~ group, data = my_data, exact = F)
res

res <- wilcox.test(weight ~ group, data = my_data, exact = F, alternative = "less")

res <- wilcox.test(weight ~ group, data = my_data, exact = F, alternative = "greater")

##== wilcoxon test III ==## repeat observations of same subject
library(MASS)
head(immer)
dim(immer); summary(immer)
wilcox.test(immer$Y1, immer$Y2, paired = T)

#==Mann-Whitney-Wilcoxon test==#
summary (mtcars)
wilcox.test(mpg ~ am, data = mtcars, exact = F)


# == Kruskal-Wallis ==# for multiple group comparison
head(airquality)
summary(airquality)
kruskal.test(Ozone ~ Month, data = airquality)



















