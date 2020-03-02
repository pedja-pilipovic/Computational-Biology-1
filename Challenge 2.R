# Download the data from the following URL:
#
#               http://membres-timc.imag.fr/Olivier.Francois/matrix_2019.txt
#
# The data consists of a matrix of 499 rows and 6421 columns with entries in 0, 1, 2. The objective of
# the challenge is to evaluate the number of clusters (for rows) in the data set and to assign a cluster
# label to each row. The result is a list of 499 cluster labels, one for each row. The output file must
# contain the resulting list formatted as a sequence of integer values separated by space characters as
# follows
#
#               12 12 1 6 7 6 6 3 11 11 ...
#
# A README file describing the options used when analyzing the data is required. The results will be
# evaluated based on the confusion matrix and the number of wrongly classified rows.
# 
# Important comments: Use a dimension reduction algorithm such as principal component analysis
# or multidimensional scaling to reduce the dimension of the data set before applying model-based
# clustering algorithms. Then, prefer using the Mclust algorithm rather than reprogramming your
# EM method.

# Turning off warnings
options(warn = -1)

# Libraries used for this challenge 
library(devtools)
library(mclust)
library(ggplot2)
library(dplyr)
library(Rtsne)
library(ggrepel)

# Importing data set
matrix_2019 <- read.table("C:/Pedja/Skola/Fakultet/M2/Computational Biology/matrix_2019.txt", quote="\"", comment.char="")

# Performing PCA
matrix.pca <- prcomp(matrix_2019)

summary(matrix.pca)
str(matrix.pca)

# Computing standard deviation of each principal component
std_dev <- matrix.pca$sdev

# Compute variance
pr_var <- std_dev^2

# Proportion of variance explained
prop_varex <- pr_var/sum(pr_var)

# Ploting two graphs to see how many components do we need
# First one is for proportion of explained variance and the
# second one is cumulative distribution function 
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")

# We can see that around 450 components explained ~95% of total variance
# Let' check that
(cumsum(prop_varex))[450]
# It is 94.4%

# Conclusion is that we need 450 components, which means our new data 
# will be of dimension 499x450

# We reduced our dimension, but it is still high, which means Mclust will
# need a lot of time to do classification, so instead of PCA we will use
# tSNE

# Instead of tsne() function we will use Rtsne() which is the same thing
# just implemented in cpp, so it is faster. For tSNE we need to choose dimension,
# max iteration and perplexity. 

fit.tsne2 <- Rtsne(matrix_2019, max_iter = 5000, perplexity = 2)
ggplot(data.frame(fit.tsne2$Y), aes(x = data.frame(fit.tsne2$Y)[,1], y = data.frame(fit.tsne2$Y)[,2])) + geom_point(alpha = 0.1) + theme_bw()

fit.tsne5 <- Rtsne(matrix_2019, max_iter = 5000, perplexity = 5)
ggplot(data.frame(fit.tsne5$Y), aes(x = data.frame(fit.tsne5$Y)[,1], y = data.frame(fit.tsne5$Y)[,2])) + geom_point(alpha = 0.1) + theme_bw()

fit.tsne10 <- Rtsne(matrix_2019, max_iter = 5000, perplexity = 10)
ggplot(data.frame(fit.tsne10$Y), aes(x = data.frame(fit.tsne10$Y)[,1], y = data.frame(fit.tsne10$Y)[,2])) + geom_point(alpha = 0.1) + theme_bw()

fit.tsne30 <- Rtsne(matrix_2019, max_iter = 5000, perplexity = 30)
ggplot(data.frame(fit.tsne30$Y), aes(x = data.frame(fit.tsne30$Y)[,1], y = data.frame(fit.tsne30$Y)[,2])) + geom_point(alpha = 0.1) + theme_bw()

fit.tsne50 <- Rtsne(matrix_2019, max_iter = 5000, perplexity = 50)
ggplot(data.frame(fit.tsne50$Y), aes(x = data.frame(fit.tsne50$Y)[,1], y = data.frame(fit.tsne50$Y)[,2])) + geom_point(alpha = 0.1) + theme_bw()

fit.tsne100 <- Rtsne(matrix_2019, max_iter = 5000, perplexity = 100)
ggplot(data.frame(fit.tsne100$Y), aes(x = data.frame(fit.tsne100$Y)[,1], y = data.frame(fit.tsne100$Y)[,2])) + geom_point(alpha = 0.1) + theme_bw()


# But after trying different variation, it is okey
# to let default values (dim = 2 and perplexity = 30) 
fit.tsne <- Rtsne(matrix_2019, max_iter = 5000)

tsne.data <- data.frame(fit.tsne$Y)
colnames(tsne.data) <- c("tsne1", "tsne2")

# Variable for data where we will add new column with class for each row
info.data <- tsne.data

# Ploting columns we got from tSNE
ggplot(tsne.data, aes(x = tsne1, y = tsne2)) + geom_point(alpha = 0.1) + theme_bw()

# Applying Mclust on data from tSNE
mc <- Mclust(tsne.data, G = 1:25)
# We let G to be any number from 1 to 25, and the optimal one is 15, so there are 15 classes

# Adding new variable for classification
info.data$mclust <- factor(mc$classification)

# Data frame having position on plane for every class.
# We need this for ploting a graph with label of every class
mc.cent <- info.data %>% group_by(mclust) %>% select(tsne1, tsne2) %>% summarize_all(mean)

# Ploting 2D plane with classes and label for each class
ggplot(info.data, aes(x = tsne1, y = tsne2, colour = mclust)) + geom_point(alpha = 0.3) + 
        theme_bw() + geom_label_repel(aes(label = mclust), data = mc.cent) + guides(colour = FALSE)

