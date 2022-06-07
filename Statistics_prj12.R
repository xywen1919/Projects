require(stats)
require(graphics)
require(jtools)
# require(sandwich)
# require(ggstance)
# require(broom.mixed)

summary(warpbreaks)
opar <- par(mfrow = c(1, 2), oma = c(0, 0, 1.1, 0))
plot(breaks ~ tension, data = warpbreaks, col = "lightgray",
     varwidth = TRUE, subset = wool == "A", main = "Wool A")
plot(breaks ~ tension, data = warpbreaks, col = "lightgray",
     varwidth = TRUE, subset = wool == "B", main = "Wool B")
mtext("warpbreaks data", side = 3, outer = TRUE)
par(opar)

# linear regression model
summary(fm1 <- lm(breaks ~ wool*tension, data = warpbreaks))
anova(fm1)
# plotting linear regression
par(mfrow = c(2,2))
plot(fm1, main = "linear regression")

# Q1 - How is it different from the results of the Poisson distribution?
# glm model
attach(warpbreaks)
hist(breaks)
mean(breaks); var(breaks)

# poisson model
poisson.model <- glm(breaks ~ wool + tension, family = poisson(link = "log"))
summary(poisson.model)
exp(coef(poisson.model)); 1-exp(coef(poisson.model))


# quasipoisson model
quasipoisson.model <- glm(breaks ~ wool + tension, family = quasipoisson(link = "log"))
summary(quasipoisson.model)

# summarize results using jtool "summ"
require(jtools)
summ(poisson.model, confint = TRUE, ci.width = .5, digits = 5)
summ(quasipoisson.model, confint = TRUE, ci.width = .5, digits = 5)

# Q3 - if possible, try to use jtools and plot_summs to plot model variables
# plotting coefficients for poisson model
p.model1 <- glm(breaks ~ tension, family = poisson(link = "log"))
p.model2 <- glm(breaks ~ wool, family = poisson(link = "log"))
p.model3 <- glm(breaks ~ tension * wool, family = poisson(link = "log"))
plot_summs(p.model1, p.model2, p.model3, scale = TRUE)

# plotting coefficients for quasipoisson model
qp.model1 <- glm(breaks ~ tension, family = quasipoisson(link = "log"))
qp.model2 <- glm(breaks ~ wool, family = quasipoisson(link = "log"))
qp.model3 <- glm(breaks ~ tension * wool, family = quasipoisson(link = "log"))
plot_summs(qp.model1, qp.model2, qp.model3, scale = TRUE)

# plotting predicted lines
effect_plot(poisson.model, pred = tension, interval = TRUE, plot.points = TRUE)
effect_plot(poisson.model, pred = wool, interval = TRUE, plot.points = TRUE)

effect_plot(quasipoisson.model, pred = tension, interval = TRUE, plot.points = TRUE)
effect_plot(quasipoisson.model, pred = wool, interval = TRUE, plot.points = TRUE)
