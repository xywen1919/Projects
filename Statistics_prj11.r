library(readr)      # for read.csv
library(tibble)
library(dendextend) # for plot(as.dendrogram(...))
library(ggplot2)
library(pheatmap)  
library(devtools)  # for loading ggbiplot pkg
library(ggbiplot)  # for ggbiplot()
library(Rtsne)
library(MASS)  # for lda()
library(dplyr) # for sampling

# dataset
ppg2008 <- read_csv("ppg2008.csv")
ppg <- column_to_rownames(ppg2008, var = "Name")

# function for z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

##calculate z-scale per each row, create normalized data
ppg_norm <- apply(ppg, 2, cal_z_score)



# hierarchy cluster - cutree k=3
ppg_hclt <- hclust(dist(ppg_norm) , method = "complete")
plot(hclust(dist(ppg_norm) , method = "complete"), hang = -1)

ppg_hclt3 <- cutree(ppg_hclt, k = 3)
ppg_hclt3


# cut_tree 3 x 3 - heat map
pheatmap(ppg_norm, cutree_rows = 3, cutree_cols = 3)


# ===== lda model === to find out contribute par for category

# add category column based on tree cut to the normalized data-set 
ppg_norm2 <- cbind(ppg_norm[order(row.names(ppg_norm)),],
                   cat = ppg_hclt3[order(names(ppg_hclt3))])
ppg_norm2 <- as.data.frame(ppg_norm2)
ppg_norm2$cat <- as.character(ppg_norm2$cat)
summary(ppg_norm2)

# lda model - step - data prepare
# sampling(sample size = n; training = ppg_s; test = ppg_t)
n = 30
set.seed(1000) 
ppg_s <- sample_n(ppg_norm2, n)
ppg_t <- subset(ppg_norm2, !(rownames(ppg_norm2) %in% rownames(ppg_s)))
nrow(ppg_s);nrow(ppg_t)

# lda model - step - train model
ppg.lda.test <- lda(cat ~., data = ppg_s)
ppg.lda.test

# lda model - step - prediction
ppg_new_cat <- predict(ppg.lda.test, ppg_t)$class

# compare result - prediction accurate %
as.numeric(ppg_new_cat);as.numeric(ppg_t$cat)
c <- sum(as.numeric(ppg_new_cat) == as.numeric(ppg_t$cat))
r <-  paste("Pediction model correct rate: ", c , " out of ", (50 - n), " = ", c/(50 - n)*100, "%")

print(r)


# ===== PCA model === to find out contribute par for category
# modeling
ppg.pca <- prcomp(ppg_norm, center = TRUE, scale. = TRUE)
print(ppg.pca); summary(ppg.pca)

# result visualize
labels <- as.character(ppg_hclt3)
ggbiplot(ppg.pca, obs.scale = 1, var.scale = 1, 
         groups = labels, labels = rownames(ppg_norm),
         ellipse = TRUE,circle = TRUE)


# ===== tSNE 

tsne <- Rtsne(ppg_norm,
              dims = 2,
              perplexity = 10,
              pca = FALSE, 
              theta = 0.5)


# display the result of t-sne
colors = rainbow(3)
plot(tsne$Y, t = 'n')
text(tsne$Y, labels = ppg_hclt3, col = colors[ppg_hclt3])
