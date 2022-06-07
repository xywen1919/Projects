rm(list = ls())
graphics.off()

#=====#
library(dplyr)
#library(bootstrap)
library(ggplot2)
# library(magrittr)
# library(reshape2)
# library(GGally)
# library(readr)

# 1) Randomly select 30 males and 30 females
rnd <- floor(runif(30, min = 1, max = 5001))
s_wh_male <- wh_male[rnd, ]
s_wh_female <- wh_female[rnd,]
# 2)combine rows - dataframe with equal n of male & female
s_wh <- bind_rows(s_wh_female, s_wh_male)

# 3) Randomly select n rows in dataset regardless of Gender
sample_wh_no_r <- sample_n(wh, 30, replace = FALSE)

#======#
# create subsets & t test

wh <- read.csv("weight-height.csv")
fh <- wh[c(5001:5015), 2]
mh <- wh[c(1:15), 2]
h <- wh[c(1:15, 5001:5015), c(1,2)]


t.test(h[,2] ~ h[,1])
t.test(fh, mh)
t.test(h[c(1:15), 2], mu = 70)    # one sample t test

#plot data
boxplot(Height ~ Gender, data = h, main = "G-H", xlab = "Gender", ylab = "Height in cm")

# create data set & one-way ANOVA test

datafilename ="http://personality-project.org/r/datasets/R.appendix1.data"
data.ex1 <- read.table(datafilename, header = T)
head(data.ex1); str(data.ex1); dim(data.ex1)

boxplot(Alertness ~ Dosage, data = data.ex1)
aov.ex1 <- aov(Alertness ~ Dosage, data = data.ex1)

print(aov.ex1)
summary(aov.ex1)   # check if there is an effect

?TukeyHSD()    # Hukey's Honest Significant Difference, create confidence intervals
TukeyHSD(aov.ex1)

plot(TukeyHSD(aov.ex1))

#====# t test & ANOVA test ==
# t test 
wh <- read.csv("weight-height.csv")
fh <- wh[c(5001:5015), 2]
mh <- wh[c(1:15), 2]
nh <- (fh$Height + mh$Height)/2

t.test(nh, mh)$p.value   # p for neutral and male
t.test(fh, mh)$p.value   # p for female and male


#  ANOVA test
h <- wh[c(1:20, 5001:5020), c(1,2)]
levels(h[,1]) <- c("Female", "Male", "Neutral")

mid <- (as.numeric(h[c(1:20), 2]) + as.numeric(h[c(21:40), 2]))/2
names <- matrix("Neutral", nrow = 20, ncol = 1)
t <- cbind(names, mid)
colnames(t) <- c("Gender", "Height")

h <- rbind(h,t)

h20.aov <- aov(Height ~ Gender, data = h)
summary(h20.aov)
TukeyHSD(h20.aov)



################
#===# Home work ===#

#################

#randomly pick 20 male and 20 female from wh data.
wh <- read.csv("weight-height.csv")
wh_male <- wh[wh[,1] == "Male", ]
wh_female <- wh[wh[,1] == "Female", ]

rnd <- floor(runif(20, min = 1, max = 5001))
s_wh_male <- wh_male[rnd, ]
s_wh_female <- wh_female[rnd,]
# combine rows to small subset
s_wh <- bind_rows(s_wh_female, s_wh_male)

# 1) t test of Height and Weight
ht <-t.test(Height ~ Gender, data = s_wh)
wt <-t.test(Weight ~ Gender, data = s_wh)

ht$p.value
wt$p.value

# 2) add 0.1 to each female height until p increase above 0.05
mh <- s_wh_male$Height
updated_fh <- s_wh_female$Height
updated_ht <- t.test(mh, updated_fh)
updated_hp <- updated_ht$p.value
n = 0

while (updated_hp < 0.05){
 updated_fh <- (updated_fh + 0.1)
 updated_ht <- t.test(mh, updated_fh)
 updated_hp <- updated_ht$p.value
 n = n + 1
}

print(n)                     # 35
print(updated_hp)            # 0.05547221


#add 0.1 to each female height until mean difference falls below 0.05
# mh <- s_wh_male$Height
# updated_fh <- s_wh_female$Height
# df_mean <- abs(mean(mh) - mean(updated_fh))
# n = 0
# 
# while (df_mean > 0.05){
#   updated_fh <- (updated_fh + 0.1)
#   df_mean <- abs(mean(mh) - mean(updated_fh))
#   n = n + 1
# }
# 
# print(n*0.1)                     # 5.6
# print(df_mean)                   # 0.03678393


# 3) chose weight of every female equal to the mean, could you add more before the test showed lack of significance?
#set female weight add 0.1 each time
mw <- s_wh_male$Weight
updated_fw <- s_wh_female$Weight 
updated_wt <- t.test(mw, updated_fw)
updated_wp <- updated_wt$p.value
n = 0

while (updated_wp < 0.05){
  updated_fw <- (updated_fw + 0.1)
  updated_wt <- t.test(mw, updated_fw)
  updated_wp <- updated_wt$p.value
  n = n + 1
}

print(n*0.1)                     # 37.3
print(updated_wp)            # 0.05060767

# set female weight to mean and add 0.1 each time

mw <- s_wh_male$Weight
updated_fw <- rep(mean(s_wh_female$Weight), length(s_wh_female$Weight)) 
updated_wt <- t.test(mw, updated_fw)
updated_wp <- updated_wt$p.value
n = 0

while (updated_wp < 0.05){
  updated_fw <- (updated_fw + 0.1)
  updated_wt <- t.test(mw, updated_fw)
  updated_wp <- updated_wt$p.value
  n = n + 1
}

print(n*0.1)                     # 77.1
print(updated_wp)            # 0.05117562


# 4) add 0.1 to each female height until p increase above 0.05, 
# relation with mean difference, and standard deviation of each group?

mh <- s_wh_male$Height
updated_fh <- s_wh_female$Height
mean_df <- c(abs(mean(mh) - mean(updated_fh)))
msd <- c(sd(mh))
fsd <- c(sd(updated_fh))
wsd <- c(sd(c(mh, updated_fh)))

mvar <- c(var(mh))
fvar <- c(var(updated_fh))

updated_ht <- t.test(mh, updated_fh)
updated_hp <- updated_ht$p.value
n = 0

while (updated_hp < 0.05){
  updated_fh <- (updated_fh + 0.1)
  mean_df <- append(mean_df, abs(mean(mh) - mean(updated_fh)))
  msd <- append(msd, sd(mh))
  fsd <- append(fsd, sd(updated_fh))
  wsd <- append(wsd, sd(c(mh, updated_fh)))
  
  mvar <- append(mvar, var(mh))
  fvar <- append(fvar, var(updated_fh))
  
  updated_ht <- t.test(mh, updated_fh)
  updated_hp <- updated_ht$p.value
  n = n + 1
}

relat_table <- data.frame(mean_df = mean_df, male_sd = msd, female_sd = fsd, whole_sd = wsd, male_var = mvar, female_var = fvar)

print(n)                     # 35
print(updated_hp)            # 0.05547221
View(relat_table)            # n = mean_df/0.1 - (sd1^2 + sd2^2) 

# II) order the male and female by height, pick the top/bottom 20 from male and female

h_male <- wh[wh[,1] == "Male", c(1, 2)]
h_female <- wh[wh[,1] == "Female", c(1, 2)]

mt <- slice_max(as.data.frame(h_male), order_by = Height, n = 20)  
ms <- slice_min(as.data.frame(h_male), order_by = Height, n = 20)
ft <- slice_max(as.data.frame(h_female), order_by = Height, n = 20)
fs <- slice_min(as.data.frame(h_female), order_by = Height, n = 20)

# a) paired t tests
t.test(mt$Height, ms$Height)$p.value
t.test(mt$Height, ft$Height)$p.value
t.test(mt$Height, fs$Height)$p.value

t.test(ms$Height, ft$Height)$p.value
t.test(ms$Height, fs$Height)$p.value

t.test(ft$Height, fs$Height)$p.value


# b) ANOVA test
mt2 <- mutate(mt, Groups = as.factor("Male_tall"))
ms2 <- mutate(ms, Groups = as.factor("Male_short"))
ft2 <- mutate(ft, Groups = as.factor("Female_tall"))
fs2 <- mutate(fs, Groups = as.factor("Female_short"))

t80_h <- rbind(mt2, ms2, ft2, fs2)

t80_h.aov <- aov(Height ~ Groups, data = t80_h)
summary(t80_h.aov)
TukeyHSD(t80_h.aov)


ggplot(data =t80_h, mapping = aes(x = Height, color = Groups))+
  geom_histogram(binwidth = 1)

# 2) how many males would you have to delete from ms to (increase)modulate significance?
ms <- slice_min(as.data.frame(h_male), order_by = Height, n = 20)
fs <- slice_min(as.data.frame(h_female), order_by = Height, n = 20)
pval <- t.test(ms$Height, fs$Height)$p.value
i = 20

while(pval < 0.0001){
  i = i - 1
  ms_updated <- slice_min(as.data.frame(h_male), order_by = Height, n = i)
  pval <- t.test(ms_updated$Height, fs$Height)$p.value
}

print(pval)
print(20 - i)


# 3) how many females would you have to delete from ft to reduce significance?
mt <- slice_max(as.data.frame(h_male), order_by = Height, n = 20)
ft <- slice_max(as.data.frame(h_female), order_by = Height, n = 20)
pval <- t.test(mt$Height, ft$Height)$p.value
i = 20

while(pval < 0.0001){
  i = i - 1
  ft_updated <- slice_max(as.data.frame(h_female), order_by = Height, n = i)
  pval <- t.test(mt$Height, ft_updated$Height)$p.value
}

print(pval)
print(20-i)

#4) 
ms <- slice_min(as.data.frame(h_male), order_by = Height, n = 20)
ft <- slice_max(as.data.frame(h_female), order_by = Height, n = 20)

m_f_ratio <- mean(ms$Height)/mean(ft$Height)
pval <- t.test(ms$Height, ft$Height)$p.value

m_f_ratio
pval


#=====#

m2s <- slice_min(as.data.frame(h_male), order_by = Height, n = 5000)
f2t <- slice_max(as.data.frame(h_female), order_by = Height, n = 5000)

m_f_ratio <- mean(m2s$Height)/mean(f2t$Height)
pval <- t.test(m2s$Height, f2t$Height)$p.value

m_f_ratio
pval


i = 5000
while (m_f_ratio > 1 | pval > 0.05) {
  i = i - 1
  m2s_updated <- slice_min(as.data.frame(m2s), order_by = Height, n = i)
  f2t_updated <- slice_max(as.data.frame(f2t), order_by = Height, n = i)
  m_f_ratio <- mean(m2s_updated$Height)/mean(f2t_updated$Height)
  pval <- t.test(m2s_updated$Height, f2t_updated$Height)$p.value
}

print(m_f_ratio)
print(pval)
print(5000-i)





