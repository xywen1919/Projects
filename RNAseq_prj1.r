---
title: "Homework1"
author: "XiaoyanWen"
date: "6/8/2021"
output:
  word_document: default
  html_document: default
  pdf_document: default
---

Homework 1	DUE: 06/18/21

## Part 2 (10 pts)

# 1.	Create a matrix (call it transcriptome) with the values below. The experiments are column names and genes are the row names. (3pts)

#	Control	Nitrogen	Phosphate	Potassium
#GeneA	89	78	77	56
#GeneB	90	99	85	97
#GeneC	78	94	99	87
#GeneD	81	83	80	79
#GeneE	62	51	99	88
```{r}
# create matrix 
transcriptome <- matrix(c(89,90,78,81,62,78,99,94,83,51,77,85,99,80,99,56,97,87,79,88), nrow = 5) 
# add rownames, colnames
rownames(transcriptome) <- paste(rep("Gene", times=5), c("A", "B", "C", "D", "E"), sep = "")
colnames(transcriptome) <- c("Control", "Nitrogen", "Phosphate", "Potassium")
# check result
transcriptome
```

# 2.	Use R code to calculate the average expression for each gene across all experiments. Call the vector (call it expression_average). ( The vector should contain 5 values – one for each gene) (3pts)
```{r}
# compute row means using apply()
expression_avarage <- apply(transcriptome, 1, mean)
# check result
expression_avarage
```
# 3.	Sort the matrix so that the gene with the highest average expression value is on top and save it into a new data frame (call it “sorted_genes”) (4pts)
```{r}
# using order() function to find out sorted index with the highest average on top
# and then sort the transcriptome rows using the order() output.
sorted_genes <- transcriptome[order(expression_avarage, decreasing = TRUE),]
# check result
sorted_genes
```

## Part 3- (20 pts) Fold change calculation:

#We will use the data file “expvalues.txt”. Remember, the first line defines the different experiments and after that every row is a different gene and its expression (RNA transcript abundance) values. The first column is the gene name, followed by the different experiments. First three experiments are replicates of the “control” experiments and the last three are replicates for the “treatments”.

# a.	Load the file expvalues.txt into R. (2pts)
```{r}
# import data
expvalues <- read.csv("C:/Users/wen_x/Downloads/trantomicsData/expvalues.txt", sep="")
# check dimension
dim(expvalues)
```

# b.	The first three columns are “control” and last three columns are “treatment” group. Use a loop or apply functions to calculate the mean of the control samples and the mean of the treatment samples for each gene. (5pts)
```{r}
# create factor for columns
expgroup <- factor(rep(c("control", "treatment"), each=3))
# apply() function to calculate group means for each row and transform the result to data frame
groupmean <- data.frame(t(apply(expvalues, 1, tapply, expgroup, mean)))
# check result wit dimension and print the top 6 rows
dim(groupmean)
head(groupmean)
```

# c.	Calculate the ratio of average treatment/ average control for each gene. (2pts)
```{r}
# calculate ratio and add as a column
groupmean$trtc_ratio <- groupmean$treatment/groupmean$control
# check result by printing top 6 rows
head(groupmean)
```

# d.	Take the log2 of the ratio. You have just calculated log fold change. (2pts)
```{r}
# calculate log2(ratio) and add as a column
groupmean$logfold <- log2(groupmean$trtc_ratio)
# check result
head(groupmean)
```

# e.	How many genes have a log2 fold change > 1 OR < -1 ? (2pts)
```{r}
#  if $logfold > 1 or < -1, assign 1; else assign 0. 
#  sum the number so we know how many genes meet the conditions
sum(ifelse((groupmean$logfold > 1 | groupmean$logfold < -1), 1, 0))
```

# f.	Save the names of the genes that have a log2 fold change > 1 into a file called “Induced_genes.txt” (2pts)
```{r}
# use which() to find out the rows that meet the condition
# extract the subset from expvalues and assign to induced_g
induced_g <- expvalues[which(groupmean$logfold > 1),]
# check subset result
head(induced_g)
dim(induced_g)
# write rownames (gene names) to file
write.csv(rownames(induced_g), file = "Induced_genes.txt", quote = FALSE, row.names = FALSE)
```

# g.	Using the same set of induced genes in the previous question, create a boxplot to show the distribution of values for each induced genes in each experiment. The x-axis should have all six experiments, and the y-axis is the expression level. Save the boxplot as a pdf file called “boxplot.pdf” (5pts)
```{r}
# open a pdf file using pdf() function
pdf(file = "boxplot.pdf")
# draw boxplot
boxplot(induced_g)
# close pdf()
dev.off()
```


