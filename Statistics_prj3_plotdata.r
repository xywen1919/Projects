# rm(list = ls())
# graphics.off()

library(ggplot2)
library(magrittr)
library(reshape2)

x <- c(1,2,3,4,5,6)
y <- x*x

help(plot)
plot(x,y, type ="l")
plot(x,y,type = "b", cex = 2.0, cex.lab = 1.5)

plot(x,y,
     main = "Title",
     sub = "subTitle",
     xlab = "x-axis",
     ylab = "y-axis",
     xlim = c(0,7),
     ylim = c(1,48))


help("par")
par(mar = c(5,4,4,2) +0.1)
plot(x,y,type = "b", pch = 21, col = "red",
     main = "Title", sub = "subtitle",
     xlab = "x-axis",
     ylab = "y-axis",
     yaxt = "n",
     lty = 3)

barplot(y)

# add more data to the plot
z <- x*y
df <- data.frame(x,y,z)

clrs <- c("red", "blue", "green")

barplot(t(as.matrix(df)), beside = TRUE, col = clrs)
legend("topleft", c("x","y","z"), cex = 1.3, bty = "n", fill = clrs)

clrs2 <- c("red", "blue", "green", "cyan", "yellow", "grey")
barplot(as.matrix(df), beside = T, col = clrs2)
legend("topleft", c("1", "2", "3", "4", "5", "6"), cex = 0.9, bty = "n", fill = clrs2)

# save figure
barplot(y)

jpeg("y-bar.jpg")
barplot(y)
dev.off()


#ggplot2

library(ggplot2)
library(magrittr)

a <- c (1,2,3,4,5,6)
b <- a*a
c <- a*b

dfa <- data.frame(a, b, c)
  

ggplot(data = dfa, aes(x = a, y = b))+
  geom_point()+
  geom_line(color = "darkblue")+
  scale_x_continuous(trans = "log2")+
  #scale_y_reverse()
  scale_y_sqrt()

dfa %>% 
  ggplot(aes(x = b, y = c))+
  geom_col(fill = "blue")


dfa %>% 
  ggplot(aes(x = b, y = c))+
  geom_line()+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")


dfa %>% 
  ggplot(aes(x = b, y = c))+
  geom_line()+
  scale_x_continuous(trans = "log", breaks = c(1,10))+
  scale_y_continuous(trans = "log", breaks = c(1, 10, 100))

dfa %>% 
  ggplot(aes(x = b, y = c))+
  geom_line()+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "reverse")

dfa %>% 
  ggplot(aes(x = b, y = c))+
  geom_line()+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "sqrt")

mpg_db <-ggplot2::mpg

ggplot(data = mpg_db, aes(displ, hwy, colour = class))+
  geom_point()

ggplot(data = mpg_db, aes(x = hwy, y = displ, colour = cyl))+
  geom_point()+
  xlab("mpg-hwy")+
  ylab("wt-displ")
  # xlim(c(0,35))+
  # ylim(c(1,6))


ggplot(data = mpg_db, aes(x = hwy, y = displ, colour = factor(cyl)))+
  geom_point()

ggplot(data = mpg_db, aes(x = hwy, y = displ, colour = factor(cyl)), size = 10)+
  geom_point()+
  guides(color = guide_legend(title = "Cylinders"))

car_db <- cars
plot(car_db)

wh_db <- read.csv("~/Biostatassignt/weight-height.csv")
save(car_db, mpg_db, wh_db, file = "DB.Rdata")

#abline
#help("abline")

plot(car_db)
abline(v = 10)
abline(h = 60, col = "red")
abline(a = 20, b =2, col = "green")

# rm(list = ls())
# graphics.off()

plot(car_db)
fit <- lm (dist ~ speed, data = car_db)
abline(fit, lwd = 3, lty = 3, col = "blue")

ft <- round(coefficients(fit), 1)

str <- paste("fit = ", ft[1], "+", ft[2], "* x")
title(main = str, col.main = "blue")


#abline (ggplot)

ggplot(data = car_db, aes(x = speed, y = dist))+
  geom_point()+
  geom_vline(xintercept = 10)+
  geom_hline(yintercept = 60, col = "red")+
  geom_abline(intercept = 20, slope = 2, col = "blue")+
  geom_abline(intercept = ft[1], slope = ft[2], col = "green", cex = 2)+
  geom_smooth(method = "lm", col = "yellow", cex = 1.5)+
  ggtitle(str)+
  theme(plot.title = element_text(color = "green", size = 14, face = "bold.italic"),
        axis.title.x = element_text(color = "blue", size = 14, face ="bold"),
        axis.title.y = element_text(color = "#993333", size = 14, face = "bold"))

#mtcars plot
View(mtcars)
ggplot(data = mtcars, aes(x = mpg, y = hp))+
  geom_point()


ggplot(data = mtcars, aes(x = mpg, y = hp, size = cyl, color = carb))+
  geom_point()


ggplot(data = mtcars, aes(x = mpg, y = hp, size = cyl, color = factor(carb)))+
  geom_point()


ggplot(data = mtcars, aes(x = mpg, y = hp, size = cyl, color = factor(carb)))+
  geom_point()+
  geom_smooth(method = "lm", col = "yellow", cex = 2)


#Quantiles
#QQ plots

#help(rt)
y <- rt(200, df = 5)
qqnorm(y)
qqline(y, col = 2)

#help(rnorm)
y <- rnorm(200, mean = 10, sd = 5)
y2 <- rnorm(100, mean = 0, sd = 2)

qqnorm(y)
qqline(y, col = 3)

qqnorm(y2)
qqline(y2, col = 4)

qqplot(y, y2)

abline(a = -3.6, b = 0.4)

# applicatoin of qq plots
View(precip)
qqnorm(precip, 
       ylab = "Precipitation [in/yr] for 70 US cities")
qqline(precip)

library(ggplot2)
ggplot(data =mpg, aes(x = displ, y = hwy, colour = class))+
  geom_point()


# reshape data

install.packages("reshape2")
library(reshape2)

name <- c("x", "y", "z")
last <- c("A", "B", "C")
age <- c(25, 29, 19)
sex <- c("male", "female", "male")
country <- c("Turkey", "Estonia", "Bosnia")
df_rsp <-data.frame(name, age, sex, country)
write.csv(df_rsp, file = "df_rsp.csv")

df_melt <- melt(df_rsp, id = "name")   #melt
write.csv(df_melt, file = "df_melt.csv")


dcast(df_melt, name ~ variable)        #dcast

df2_rsp <- data.frame(name, last, age, sex, country)
df2_melt <- melt(df2_rsp, id = c("name", "last"), 
                 variable.name = "details")      #melt
dcast(df2_melt, name + last ~ details)          #dcast

#Time series data
help("runif")

x <- seq(0, 4*pi, 0.1)
n <- length(x)

y1 <- 0.5 * runif(n) + sin(x)
y2 <- 0.5 * runif(n) + cos(x) - sin(x)

plot(x, y1, col = "blue", pch = 20)
plot(x, y2, col = "red", pch = 23)

df <- data.frame(x, y1, y2)
df_long <- melt(df, id = "x")  #melt

ggplot(data = df_long, aes(x = x, y = value, color = variable ))+
  geom_point()

ggplot(data = df_long, aes(x = x, y = value, color = variable, size = value ))+
  geom_point()


ggplot(data = df_long, aes(x = x, y = value, color = variable, size = abs(value) ))+
  geom_point()


#Facets

ggplot(data = df_long, aes(x = x, y = value, color = variable ))+
  geom_point()+
  facet_grid(. ~ variable)                 # facet_grid

ggplot(data = df_long, aes(x = x, y = value, color = variable ))+
  geom_point()+
  facet_wrap(~ variable, nrow = 2)         # facet_wrap

ggplot(data = df_long, aes(x = x, y = value, color = variable ))+
  geom_point()+
  facet_wrap ( ~ variable, ncol = 2)          # facet_wrap


# multipanel plots in basic R
mtcars %>% 
par(mfrow = c(2,2)) 

plot(wt, mpg, main = "wt vs mpg")  
plot(wt, disp, main = "wt vs disp")  
hist(wt, main = "wt") 
boxplot(wt, main = "wt")





















