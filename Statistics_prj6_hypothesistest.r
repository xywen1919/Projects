rm(list = ls())
graphics.off()

install.packages("gplots")
install.packages("corrplot")
library(gplots)
library(corrplot)
library(magrittr)
library(dplyr)
install.packages("tidyr")
library(tidyr)


file_path <- "http://www.sthda.com/sthda/RDoc/data/housetasks.txt"
housetasks <- read.delim(file_path, row.names = 1)
dt <- as.table(as.matrix(housetasks))  # convert data to table
dt

balloonplot(t(dt), main = "housetasks", xlab = "", ylab = "", label = F, show.margins = F)
chisq <- chisq.test(housetasks)
chisq
chisq$observed
round(chisq$expected, 2)
round(chisq$residuals, 3)
corrplot(chisq$residuals, is.corr = F)

chisq$statistic
?corrplot

# contribution in percentage
contrib <- 100*chisq$residuals^2/chisq$statistic
corrplot(round(contrib, 3), is.corr = F)


?p.adjust
pvalues <- c(0.002, 0.003, 0.015, 0.113, 0.222, 0.227, 0.454, 0.552, 0.663, 0.751)
p.adjust(pvalues, method = "BY")

#==============================
######### Assignment ##########
#==============================

rm(list = ls())
graphics.off()

###1) sample 20 people from weight-height data
library(dplyr)

wh <- read.csv("weight-height.csv")
wh_s <- sample_n(wh, 20)
wh_s_gender <- table(wh_s$Gender)  #count the number of female and male
nullprobs <- c(.5, .5)  #set up null 

gchq <- chisq.test(wh_s_gender, p = nullprobs)
gchq$p.value

#=======

gchq_pvalue <-c()
gchq_power <- c()

for (i in c(20, 40, 60, 80, 100)){
  wh_s_n <- sample_n(wh, i)
  gchqi <- chisq.test(table(wh_s_n$Gender), p = nullprobs)
  gchq_pvalue <- append(gchq_pvalue, gchqi$p.value)
  gchq_power <- append(gchq_power, power.prop.test(n = i, 
                                                   p1 = prop.table(table(wh_s_n$Gender))[1], 
                                                   p2 = prop.table(table(wh_s_n$Gender))[1])$power)
}
gchq_pvalue
gchq_power

#=== how many people need to be sampled to get a reliable estimate for the fraction of male/female?
wh_malefraction <- prop.table(table(wh$Gender))[2]
cpower <- 0
i <- 100

while (cpower < 0.8) {
  i <- i + 10
  wh_s_n <- sample_n(wh, i)
  cpower <- power.prop.test(n = i, 
                p1 = prop.table(table(wh_s_n$Gender))[2], 
                p2 = wh_malefraction)$power
}
i
cpower


#=== how many people do we need to sample in NYC (population 8 million) for a reliable estimate of male/female in the city?
nyc_male <-array("male", dim = c(8*10^6 * 0.47,1))
nyc_femal <- array("female", dim = c(8*10^6 * 0.53, 1))
nyc_gender <- rbind(nyc_male, nyc_femal)
table(nyc_gender)                     # build nyc gender table 

nyc_male_fraction <- prop.table(table(nyc_gender))[2]    # nyc male proportion

cpower <- 0                   
i <- 1000

while (cpower < 0.8) {
  i <- i + 1000
  nyc_gender_s <- sample(nyc_gender, i)
  cpower <- power.prop.test(n = i, 
                            p1 = prop.table(table(nyc_gender_s))[2], 
                            p2 = nyc_male_fraction)$power
}
i
cpower


###2) expression data from genes from two groups of samples:
library(dplyr)

genexpr <- read.csv("diffexpr-results_no_p.csv", row.names = 1, header = T, stringsAsFactors = F)
dim(genexpr); str(genexpr); head(genexpr); summary(genexpr)

pvalue <- c()
genename <- c()
for (i in c(1:length(genexpr$Gene))){
c <-t.test(genexpr[i, c(2:5)], genexpr[i, c(6:9)])
pvalue <- append(pvalue, c$p.value)
genename <- append(genename, genexpr[i, 1])
}
gene_p <- bind_cols(genename, pvalue)
colnames (gene_p) <- c("gene_name", "p_value")

summary(gene_p)
# head(gene_p)
# str(gene_p)
# summary(complete.cases(gene_p))


library(tidyr)
gene_p_1 <- replace_na(gene_p, list(p_value = 1))
gene_p_1 <- mutate(gene_p_1, FDR_BH = p.adjust(p_value, method = "BH"))
gene_p_1 <- mutate(gene_p_1, bonferroni = p.adjust(p_value, method = "bonferroni"))
summary(gene_p_1)


# select FDR < 10%
gene_p_fdr_10 <- gene_p_1[(gene_p_1$FDR_BH < 0.1),]
summary(gene_p_fdr_10)  

# select FDR < 20%
gene_p_fdr_20 <- gene_p_1[(gene_p_1$FDR_BH < 0.2),]
summary(gene_p_fdr_20)

#select bonferroni < 10%
gene_p_bonferroni_10 <- gene_p_1[(gene_p_1$bonferroni < 0.1),]
summary(gene_p_bonferroni_10)  

#select bonferroni < 20%
gene_p_bonferroni_20 <- gene_p_1[(gene_p_1$bonferroni < 0.2),]
summary(gene_p_bonferroni_20) 

