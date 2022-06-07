library(readr)
library(readxl)
library(magrittr)
library(arsenal)
library(dplyr)
library(energy)


# import data
Substracted_Dataset <- read_csv("Substracted Dataset.csv") %>% as.data.frame()
stategundata2 <- read_excel("stategundata2.xlsx")
stategundata <- stategundata2[,c(2,5)] %>% as.data.frame() 
str(Substracted_Dataset); str(stategundata)

plotdata <- round(Substracted_Dataset,1)
colnames(plotdata) <- c("homicideRate", "bradyScore")

paperdata <- round(stategundata2[,c(2,5)],1)
colnames(paperdata) <- c("homicideRate", "bradyScore")
paperdata <- arrange(paperdata, homicideRate)



# compare two data sets
summary(plotdata) 
summary(paperdata)

m <- par(mfrow = c(1,2))
plot(plotdata, pch = 2, col = "blue", 
     xlab = "Homicide rate", ylab = "Brady score", main = "Plotdata")
plot(paperdata, pch = 3, col = "red", 
     xlab = "Homicide rate", ylab = "Brady score", main = "Paperdata")
par(m)

all_equal(plotdata,paperdata)
summary(comparedf(plotdata,paperdata, tol.num.val = 1))[6:7]


# Analyze the data (paper data) using Pearson's correlation
cor(paperdata$homicideRate, paperdata$bradyScore, method = "pearson")
cor.test(paperdata$homicideRate, paperdata$bradyScore, method = "pearson")


# Analyze the data using distance correlation
dcor(paperdata$homicideRate, paperdata$bradyScore)
dcorT.test(paperdata$homicideRate, paperdata$bradyScore)

# #?
# plot(paperdata,pch = 8, xlab = "Homicide rate", ylab = "Brady score")
# abline(b = .0655, a = .45956, col = "blue", lty = 3, lwd = 2)
# abline(b = -0.03, a = -1.052, col = "red", lty = 2, lwd = 2)


