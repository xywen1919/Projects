rm(list = ls())   
graphics.off()

#### Part I ##############

###########################

cell_id <- read.csv("cell_id.tab-2.csv", row.names = 1, header = T, stringsAsFactors = F)
str(cell_id)

# 1) How many rows have 293T in their name ?
grep("293T", rownames(cell_id))
cell_id[grep("293T", rownames(cell_id)), ]


# 2) What is the name of marker in column 7 in the table
colnames(cell_id)[6]

# 3) What is the entry in row 2, col 7 ?
cell_id[1, 6]

# 4) For the column labeled "D13S317 " what is most common entry ?
sort(table(cell_id$D13S317), decreasing = T)

#### Part II ##############

###########################
# use Use the weight-height data (weight height.csv)
wh <- read.csv("weight-height.csv")
wh_male <- wh[wh$Gender == "Male", ]
wh_female <- wh[wh$Gender == "Female", ]

library(dplyr)
library(bootstrap)

# 1) Randomly select 30 males and 30 females
rnd <- floor(runif(30, min = 1, max = 5001))
s_wh_male <- wh_male[rnd, ]
s_wh_female <- wh_female[rnd,]

s_wh <- bind_rows(s_wh_female, s_wh_male)

# a) Estimate mean and variance using jackknife 
#for the 30 males, 30 females and all 60 together 
#===== Height =======
# Male
jack_male_height_mn <- jackknife(s_wh_male$Height, theta = mean)
mean_jack_male_height_mn <- round(mean(jack_male_height_mn[[3]]), 2)

jack_male_height_va <- jackknife(s_wh_male$Height, theta = var)
mean_jack_male_height_va <- round(mean(jack_male_height_va[[3]]), 2)

# Female
jack_female_height_mn <- jackknife(s_wh_female$Height, theta = mean)
mean_jack_female_height_mn <- round(mean(jack_female_height_mn[[3]]), 2)

jack_female_height_va <- jackknife(s_wh_female$Height, theta = var)
mean_jack_female_height_va <- round(mean(jack_female_height_va[[3]]), 2)

# All
jack_all_height_mn <- jackknife(s_wh$Height, theta = mean)
mean_jack_all_height_mn <- round(mean(jack_all_height_mn[[3]]), 2)

jack_all_height_va <- jackknife(s_wh$Height, theta = var)
mean_jack_all_height_va <- round(mean(jack_all_height_va[[3]]), 2)


#===== Weight =======
# Male
jack_male_weight_mn <- jackknife(s_wh_male$Weight, theta = mean)
mean_jack_male_weight_mn <- round(mean(jack_male_weight_mn[[3]]), 2)

jack_male_weight_va <- jackknife(s_wh_male$Weight, theta = var)
mean_jack_male_weight_va <- round(mean(jack_male_weight_va[[3]]), 2)

# Female
jack_female_weight_mn <- jackknife(s_wh_female$Weight, theta = mean)
mean_jack_female_weight_mn <- round(mean(jack_female_weight_mn[[3]]), 2)

jack_female_weight_va <- jackknife(s_wh_female$Weight, theta = var)
mean_jack_female_weight_va <- round(mean(jack_female_weight_va[[3]]), 2)

# All
jack_all_weight_mn <- jackknife(s_wh$Weight, theta = mean)
mean_jack_all_weight_mn <- round(mean(jack_all_weight_mn[[3]]), 2)

jack_all_weight_va <- jackknife(s_wh$Weight, theta = var)
mean_jack_all_weight_va <- round(mean(jack_all_weight_va[[3]]), 2)


# b) Estimate mean and variance using bootstrap 
#for the 30 males, 30 females and all 60 together 

s_db <- list(male = s_wh_male, female = s_wh_female, all = s_wh)

#=== Height ===
mean_boot_height_mn <- c(1:3)*0
mean_boot_height_va <- c(1:3)*0

for (i in 1:3){
  boot_height_mn <- bootstrap(s_db[[i]]$Height, 1000, theta = mean)
  mean_boot_height_mn [i] <- round(mean(boot_height_mn$thetastar), 2)
  
  boot_height_va <- bootstrap(s_db[[i]]$Height, 1000, theta = var)
  mean_boot_height_va [i] <- round(mean(boot_height_va$thetastar), 2)
  
  }
mean_boot_height_mn
mean_boot_height_va



#=== Weight ===
mean_boot_weight_mn <- c(1:3)*0
mean_boot_weight_va <- c(1:3)*0

for (i in 1:3){
  boot_weight_mn <- bootstrap(s_db[[i]]$Weight, 1000, theta = mean)
  mean_boot_weight_mn [i] <- round(mean(boot_weight_mn$thetastar), 2)
  
  boot_weight_va <- bootstrap(s_db[[i]]$Weight, 1000, theta = var)
  mean_boot_weight_va [i] <- round(mean(boot_weight_va$thetastar), 2)
  
}
mean_boot_weight_mn
mean_boot_weight_va


# c) Do the estimates change between 1000 and 10000 sampling in bootstrap? 

#=== Height ===
mean_boot_height_mn_4z <- c(1:3)*0
mean_boot_height_va_4z <- c(1:3)*0

for (i in 1:3){
  boot_height_mn_4z <- bootstrap(s_db[[i]]$Height, 10000, theta = mean)
  mean_boot_height_mn_4z [i] <- round(mean(boot_height_mn_4z$thetastar), 2)
  
  boot_height_va_4z <- bootstrap(s_db[[i]]$Height, 10000, theta = var)
  mean_boot_height_va_4z [i] <- round(mean(boot_height_va_4z$thetastar), 2)
  
}
mean_boot_height_mn_4z
mean_boot_height_va_4z



#=== Weight ===
mean_boot_weight_mn_4z <- c(1:3)*0
mean_boot_weight_va_4z <- c(1:3)*0

for (i in 1:3){
  boot_weight_mn_4z <- bootstrap(s_db[[i]]$Weight, 10000, theta = mean)
  mean_boot_weight_mn_4z [i] <- round(mean(boot_weight_mn_4z$thetastar), 2)
  
  boot_weight_va_4z <- bootstrap(s_db[[i]]$Weight, 10000, theta = var)
  mean_boot_weight_va_4z [i] <- round(mean(boot_weight_va_4z$thetastar), 2)
  
}
mean_boot_weight_mn_4z
mean_boot_weight_va_4z


# 2) Using the sample function in R, write a bootstrap routine 
#that outputs the bootstrap mean. Explain the logic by annotating the code.

wh_sample_height <- sample(wh$height, 30, replace = TRUE)  # sampling 30 individuals' height in the weight-height database with replacement. 
                                          #Note that sampling with replacement is used in bootstrap, and without replacement to give a permutation
theta <- function(x){mean(x)}             # set mean as function to be bootstrapped, since mean is a built-in method, this stop could be omitted...
result <- bootstrap(wh_sample_height, 100, theta)# 100 bootstraps of the sample mean
round(mean(result$thetastar), 2)          # the bootstrap values of theta rounded to two decimal places.


# 3) Sample without replacement people's heights and weights from our weight height data. 
#Demonstrate the correlation between height and weight
#=================
# #Randomly Select Rows In R Dataframes
# 
# #method 1
# # r sample dataframe; selecting a random subset in r
# # df is a data frame; pick 5 rows
# 
# df[sample(nrow(df), 5), ]
# 
# #method 2
# # dplyr r sample_n example
# sample_n(df, 10)
#================

library(dplyr)
sample_wh_no_r <- sample_n(wh, 30, replace = FALSE)

#Spearman-rank correlation for the full dataset:
cor(wh[,2], wh[,3], method = c("spearman"))  

#Spearman-rank correlation for the sampled dataset:
cor(sample_wh_no_r[,2], sample_wh_no_r[,3], method = c("spearman"))




