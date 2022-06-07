#packages
### to install from github
# install.packages("devtools")
# install_github("vqv/ggbiplot")
# install Rtool40 (https://cran.r-project.org/bin/windows/Rtools/)
# install.packages("rgl")

library(readxl)
library(readr) # for read.csv
library(tibble)
library(dplyr)
library(dendextend)  # for plot(as.dendrogram(...))
library(graphics)
library(energy)  # for dcor
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(pheatmap)
library(MASS)     # for lda()
library(MultiRNG) # for generate.point.in.sphere()
library(devtools)  # for loading ggbiplot pkg from gethub_
library(ggbiplot)  # for ggbiplot()
library(Rtsne)     # for Rtsne() plot
library(animation)  # for animation
library(fpc)        # for animation
library(rgl)        # for 3D t-SNE display
library(tidyverse)
library(ggtree)     # for tree plot
library(treeio)     # for tree.data read



############################################
## function
## #########################################
# ======== dist ================
?dist
x <- matrix(rnorm(100), nrow = 5)
dist(x)
dist(x, diag = TRUE)
dist(x, upper = TRUE)
m <- as.matrix(dist(x))
d <- as.dist(m)

## Use correlations between variables "as distance"
str(USJudgeRatings)
dd <- as.dist((1 - cor(USJudgeRatings))/2)
round(1000 * dd) # (prints more nicely)
plot(hclust(dd)) # to see a dendrogram of clustered variables
#
#
# ======== hclust ================
?hclust
### Example 1: Violent crime rates by US state
str(USArrests)
hc <- hclust(dist(USArrests), "ave")
plot(hc)
plot(hc, hang = -1)

## Do the same with centroid clustering and *squared* Euclidean distance,
## cut the tree into ten clusters and reconstruct the upper part of the
## tree from the cluster centers.
hc <- hclust(dist(USArrests)^2, "cen")
memb <- cutree(hc, k = 10)
cent <- NULL
for (k in 1:10) {
  cent <- rbind(cent, colMeans(USArrests[memb == k, , drop = FALSE]))
}
hc1 <- hclust(dist(cent)^2, method = "cen", members = table(memb))
opar <- par(mfrow = c(1, 2))
plot(hc,  labels = FALSE, hang = -1, main = "Original Tree")
plot(hc1, labels = FALSE, hang = -1, main = "Re-start from 10 clusters")
par(opar)


### Example 2: Straight-line distances among 10 US cities
##  Compare the results of algorithms "ward.D" and "ward.D2"
?cmdscale
mds2 <- -cmdscale(UScitiesD)
plot(mds2, type = "n", axes = FALSE, ann = FALSE)
text(mds2, labels = rownames(mds2), xpd = NA)

hcity.D  <- hclust(UScitiesD, "ward.D") # "wrong"
hcity.D2 <- hclust(UScitiesD, "ward.D2")
opar <- par(mfrow = c(1, 2))
plot(hcity.D,  hang = -1)
plot(hcity.D2, hang = -1)
par(opar)
#
#
# ====== phylogenetic tress | ggtree ======
# https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
# library(ggtree)
set.seed(2017-02-16)
tree <- rtree(50)
ggtree(tree)
ggtree(tree, layout="slanted") 
ggtree(tree, layout="circular")
ggtree(tree, layout="fan", open.angle=120)
ggtree(tree, layout="equal_angle")
ggtree(tree, layout="daylight")
ggtree(tree, branch.length='none')
ggtree(tree, branch.length='none', layout='circular')
ggtree(tree, layout="daylight", branch.length = 'none')
# 
# 
# Time-scaled layout
# library(treeio) 
beast_file <- system.file("examples/MCC_FluA_H3.tree", 
                          package="ggtree")
beast_tree <- read.beast(beast_file)
ggtree(beast_tree, mrsd="2013-01-01") + 
  theme_tree2()
# 
# ---------- *
# https://4va.github.io/biodatasci/r-ggtree.html
# library(tidyverse)
# library(ggtree)
# library(treeio)

tree <- read.tree("tree_newick.nwk")
tree

ggtree(tree) +
  geom_nodepoint()+
  geom_tippoint()+
  geom_tiplab()+
  ggtitle()+
  theme_tree2()

ggtree(tree, branch.length = "none") 
ggtree(tree, branch.length = "none", color = "blue", size = 2, linetype = 3) 

ggtree(tree)+
  geom_text(aes(label = node), hjust = -.3)

ggtree(tree)+
  geom_tiplab()

ggtree(tree) + 
  geom_tiplab() + 
  geom_cladelabel(node=17, label="Some random clade", 
                  color="red2", offset=.8, align=TRUE) + 
  geom_cladelabel(node=21, label="A different clade", 
                  color="blue", offset=.8, align=TRUE) + 
  theme_tree2() + 
  xlim(0, 70) + 
  theme_tree()

ggtree(tree) + 
  geom_tiplab() + 
  geom_hilight(node=17, fill="gold") + 
  geom_hilight(node=21, fill="purple")
# 
# ----- 
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
tree

ggtree(tree)
ggtree(tree) + geom_treescale()
ggtree(tree, branch.length="none")
ggtree(tree, layout="circular") + ggtitle("(Phylogram) circular layout")
# 
# With the groupClade and groupOTU methods you can cluster clades or related OTUs, and assign them different colors for example
# 
tree <- groupClade(tree, node=c(21, 17))
ggtree(tree, aes(color=group, linetype=group)) + geom_tiplab(aes(subset=(group==2)))
# 
## ----- 2009 H1N1.china -------------------------
# at http://beast.bio.ed.ac.uk/
# Read the data
tree <- read.beast("flu_tree_beast.tree")

# supply a most recent sampling date so you get the dates
# and add a scale bar
ggtree(tree, mrsd="2013-01-01") + 
  theme_tree2() 

# Finally, add tip labels and adjust axis
ggtree(tree, mrsd="2013-01-01") + 
  theme_tree2() + 
  geom_tiplab(align=TRUE, linesize=.5) + 
  xlim(1990, 2020)
##
# ======== pheatmap ================
?pheatmap
# Create test matrix
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")
summary(test)

# Draw heatmaps
pheatmap(test)
pheatmap(test, kmeans_k = 2)
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(test, cluster_row = FALSE)
pheatmap(test, legend = FALSE)

# Show text within cells
pheatmap(test, display_numbers = TRUE)

# Gaps in heatmaps
pheatmap(test, cluster_rows = FALSE, gaps_row = c(10, 14))
pheatmap(test, cluster_rows = FALSE, gaps_row = c(10, 14), cutree_col = 2)

# Specifying clustering from distance matrix
drows = dist(test, method = "minkowski")
dcols = dist(t(test), method = "minkowski")
pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
#
#
# ======== lda ================
?lda
?predict

Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
train <- sample(1:150, 75)
table(Iris$Sp[train])

z <- lda(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train)
predict(z, Iris[-train, ])$class
#
#
# ======== kmeans ================
?kmeans
library(graphics)
# a 2-dimensional example
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(x) <- c("x", "y")
(cl <- kmeans(x, 2))
plot(x, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex = 2)

## random starts do help here with too many clusters
## (and are often recommended anyway!):
(cl <- kmeans(x, 5, nstart = 25))
plot(x, col = cl$cluster)
points(cl$centers, col = 1:5, pch = 8)
#
#
# ======== prcomp ================
?prcomp

## the variances of the variables in the
## USArrests data vary by orders of magnitude, so scaling is appropriate
prcomp(USArrests)  # inappropriate
prcomp(USArrests, scale = TRUE)
prcomp(~ Murder + Assault + Rape, data = USArrests, scale = TRUE)
plot(prcomp(USArrests))
summary(prcomp(USArrests, scale = TRUE))
biplot(prcomp(USArrests, scale = TRUE))
#
#
# ======== ggbiplot ================
?ggbiplot
data(wine)
wine.pca <- prcomp(wine, scale. = TRUE)
print(ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, groups = wine.class, ellipse = TRUE, circle = TRUE))
#
#
# ======== Rtsne ================
?Rtsne
iris_unique <- unique(iris) # Remove duplicates
iris_matrix <- as.matrix(iris_unique[,1:4])

# Set a seed if you want reproducible results
set.seed(42)
tsne_out <- Rtsne(iris_matrix,pca = FALSE,perplexity = 30,theta = 0.0) # Run TSNE

# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col = iris_unique$Species, asp = 1)
#
#
#
###################################################
## data-set
## ################################################
#
# ==== tagseq ================
data <- read_excel("TagSeqExample.tab.xlsx")
data <- column_to_rownames(data, var = "gene")

## make a matrix of only highly expressed genes
data_subset <- as.matrix(data[rowSums(data) > 50000,])

# ==== NBA2008 ================
nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
nba$Name <- with(nba, reorder(Name, PTS))
nba.m <- melt(nba)
nba.m <- ddply(nba.m, .(variable), transform, rescale = rescale(value))
head(nba.m)

# ==== weight_height ================

wh <- read_csv("weight-height.csv")

## add to the levels of Gender, so you can change labels to M and F
str(wh)
levels(wh$Gender) <- c(levels(wh$Gender),"M","F")
wh[which(wh[,"Gender"] == "Male"),"Gender"] = "M"
wh[which(wh[,"Gender"] == "Female"),"Gender"] = "F"
head(wh)
tail(wh)

# ==== concentric ring of points (df + labels) ================
x <- generate.point.in.sphere(100,2)
x[,1] <- x[,1] + rnorm(100)*0.1
x[,2] <- x[,2] + rnorm(100)*0.1
plot(x[,1],x[,2])

y <- 10*generate.point.in.sphere(100,2)
y[,1] <- y[,1] + rnorm(100)
y[,2] <- y[,2] + rnorm(100)
plot(y[,1],y[,2],col = "red")
points(x[,1],x[,2],col = "blue")

df <- rbind(x,y)
colnames(df) <- c("X","Y")
labels <- c(rep("2",100),rep("1",100))
#
#
#
#################################################
##Clustering
#################################################
### dataset tagseq
### 
## ===== Hierarchical Clustering 1 (Euclidean) =====
# distance matrix based on euclidean distance
## 
dm <- dist(data_subset)
my_hclust_gene <- hclust(dm, method = "complete")
plot(my_hclust_gene, hang = -1)

par(mar = c(5,5,5,12))
nPar <- list(lab.cex = 0.6, pch = c(NA, 19),cex = 0.7, col = "blue")
ePar = list(col = 2:3, lwd = 2:1)
plot(as.dendrogram(my_hclust_gene),nodePar = nPar, edgePar = ePar,horiz = TRUE)
#
#
#
## ===== Hierarchical Clustering 2 (Pearson correlation) =====
## distance matrix by converting pearson correlation to a distance
##
dm <- as.dist((1 - cor(t(data_subset), method = c("pearson")))/2)
my_hclust_gene <- hclust(dm, method = "complete")
plot(my_hclust_gene, hang = -1)
plot(as.dendrogram(my_hclust_gene))

par(mar = c(5,5,5,12))
nPar <- list(lab.cex = 0.6, pch = c(NA, 19),cex = 0.7, col = "blue")
ePar = list(col = 2:3, lwd = 2:1)
plot(as.dendrogram(my_hclust_gene),nodePar = nPar, edgePar = ePar,horiz = TRUE)


#===== Hierarchical Clustering 3 (Spearman correlation) =====
## distance matrix from spearman correlation
## 
dm <- as.dist((1 - cor(t(data_subset), method = c("spearman")))/2)
my_hclust_gene <- hclust(dm, method = "complete")
plot(my_hclust_gene, hang = -1)
plot(as.dendrogram(my_hclust_gene))

par(mar = c(5,5,5,12))
nPar <- list(lab.cex = 0.6, pch = c(NA, 19),cex = 0.7, col = "blue")
ePar = list(col = 2:3, lwd = 2:1)
plot(as.dendrogram(my_hclust_gene),nodePar = nPar, edgePar = ePar,horiz = TRUE)


#====== Hierarchical Clustering 4 (distance correlation) =====
## dcor
##
dm <- dim(data_subset)

mm <- matrix(0,dm[1],dm[1])
for (i in 1:dm[1]) { 
  for (j in 1:dm[1]) { 
    mm[i,j] = dcor(data_subset[i,],data_subset[j,],1.5);
  }
}

rownames(mm) <- rownames(data_subset)
colnames(mm) <- rownames(data_subset)
# View(mm)

dst <- as.dist(1 - mm)

my_hclust_gene <- hclust(dst, method = "complete")
plot(my_hclust_gene, hang = -1)


par(mar = c(5,5,5,12))
nPar <- list(lab.cex = 0.6, pch = c(NA, 19),cex = 0.7, col = "blue")
ePar = list(col = 2:3, lwd = 2:1)
plot(as.dendrogram(my_hclust_gene),nodePar = nPar, edgePar = ePar,horiz = TRUE)
#
#
#
##====== Clustering 5 (kmeans) =====
# data-set: weight_height 
#
#
wh.km <- kmeans(wh[,c(2:3)],2)
aggregate(wh[,c(2:3)],by = list(cluster = wh.km$cluster), mean)
wh.km$centers

#plot cluster
plot(wh[,c(2:3)], col = wh.km$cluster)
points(wh.km$centers, col = 1:2, pch = 8)


## how many mistakes ?
miss1 <- which(wh.km$cluster[c(1:5000)] == 2)
miss2 <- which(wh.km$cluster[c(5001:10000)] == 1)
length(miss1) 
length(miss2) 

## watch animation of kmeans process using animation & fpc pkgs

kmeans.ani(wh[,c(2:3)],2)
#
#
#
##################################################
## Heat maps
##################################################
##  
# ====== heapmaps 1 ===== ggplot
## dataset NBA2008
p <- ggplot(nba.m, aes(variable, Name)) + 
  geom_tile(aes(fill = rescale), colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue")
p

# ====== heatmaps 2 ===== pheatmap
## dataset tagseq

pheatmap(data_subset)

# ===== heatmaps 3 (scaled normalized data set) === z-cal === pheatmap

##scale rows and normalize

## function to scale values in a list
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

##cal z-scale each row, to create normalized data
data_subset_norm <- t(apply(data_subset, 1, cal_z_score))

# heap map
pheatmap(data_subset_norm)
#
#
#
# ===== heatmaps 4 (clustering tree cut) =====

my_hclust_gene <- hclust(dist(data_subset), method = "complete")
plot(my_hclust_gene, hang = -1)

#split tree into two clusters 
my_gene_col <- cutree(my_hclust_gene, k = 2)
head(my_gene_col)

#split the rows and columns by the tree
pheatmap(data_subset,cutree_rows = 4,cutree_cols = 2)
#
#
#
#################################################
## Dimension reduction
#################################################
# data-set weight_height
# df+label concentric ring of points
# 
# ===== LDA - supervised Linear discriminant analysis =====

r <- lda(formula = Gender ~ ., data = wh)
r

test <- data.frame(Gender = "Male", Height = 69., Weight = 187)
predict(r,newdata = test)


test2 <- wh[c(1:4),]
predict(r,test2)


test3 <- wh[c(5001:5004),]
predict(r,test3)
#
#
#
# ===== PCA - unsupervised Principal component analysis =====
# ===== prcomp =====
# data-set: df+labels concentric ring of points

df.pca <- prcomp(df,center = TRUE, scale = TRUE)
print(df.pca) ## shows the coordinates
plot(df.pca)  ## PCA FAILS to see structure
biplot(df.pca)
summary(df.pca)
#
#
#
# ===== PCA again (regular plot) =====
# ===== princomp =====

pca2 = princomp(df)$scores[,1:2]
head(pca2)

plot(pca2, t = 'n', main = "pca", "cex.main" = 2, "cex.lab" = 1.5)   
# empty plot with no points 

text(pca2, labels = labels,col = labels)   
# PCA document in download gives more insight

# ===== PCA visualization ===== ggbiplot ===

g <- ggbiplot(df.pca, 
              obs.scale = 1, 
              var.scale = 1, 
              groups = labels,
              ellipse = TRUE,
              circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')
g

#
#
#
# ====== ## height-weight PCA =====

wh.pca <- prcomp(wh[,c(2:3)], center = TRUE, scale. = TRUE)
labels <- wh$Gender
print(wh.pca)
plot(wh.pca)  # show variance per PC
biplot(wh.pca)
summary(wh.pca)

tail(wh,2)

predict(wh.pca,newdata = tail(wh[,c(2:3)],2)) 

ggbiplot(wh.pca, obs.scale = 1, var.scale = 1, groups = wh$Gender, 
         ellipse = TRUE,circle = TRUE)

ggbiplot(wh.pca, obs.scale = 1, var.scale = 1,groups = wh$Gender,
         ellipse = TRUE, circle = TRUE, cex.main = 2., cex = 2.) + 
  scale_color_discrete(name = '') + 
  geom_point(aes(colour = wh$Gender), size = 0.5) + 
  theme(legend.direction = 'horizontal',legend.position = 'top')

ggbiplot(wh.pca,labels = labels,ellipse = TRUE, groups = wh$Gender)

#
#
# ===== tSNE - t distributed stochastic neighbor embedding =====

tsne <- Rtsne(df, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

# display results of t-SNE
colors = rainbow(length(unique(labels)))
names(colors) = unique(labels)
par(mgp = c(2.5,1,0))

plot(tsne$Y, t = 'n', main = "tSNE",
     xlab = "tSNE dimension 1",
     ylab = "tSNE dimension 2",
     "cex.main" = 2, "cex.lab" = 1.5)
text(tsne$Y, labels = labels,col = colors[labels])
#S stands for Stochastic, so running it twice will not give same results

tsne2 <- Rtsne(df, pca = FALSE, perplexity = 30, theta = 0.5, dims = 2)

# display results of t-SNE
cols <- rainbow(10)
plot(tsne2$Y, t = 'n')
text(tsne2$Y, labels = labels, col = cols[labels])

# display results of t-SNE in 3D 
# set the dims parameter to 3 instead of 2
tsne3 <- Rtsne(df, pca = FALSE, perplexity = 30, theta = 0.5, dims = 3)

# display results of t-SNE using rgl pkg
cols <- rainbow(2)
plot3d(tsne3$Y, col = cols[labels])
legend3d("topright", legend = '1':'2', pch = 16, col = rainbow(2))
#
#
#
# === Weight Height tSNE =====

set.seed(1) # for reproducibility
tsne <- Rtsne(wh[,c(2:3)],
              dims = 2,
              perplexity = 30,
              pca = FALSE, 
              theta = 0.5)

lbls <- as.character(labels)

# display the result of t-sne

colors <- rainbow(length(unique(lbls)))
names(colors) = unique(lbls)
par(mgp = c(2.5,1,0))

plot(tsne$Y, t='n', main="tSNE",
     xlab="tSNE dimension 1",
     ylab="tSNE dimension 2",
     "cex.main"=2, "cex.lab"=1.5)
text(tsne$Y, labels = lbls, col = colors[lbls])

# - or -
par(mfrow = c(1,2))
plot(tsne$Y, col = "blue", pch = 19, cex = 1)
plot(tsne$Y, col = colors[lbls], pch = 21, cex = 1.5)
