# 9. Download the data from the following URL:
# 
#       http://membres-timc.imag.fr/Olivier.Francois/data2.txt
#
#    Apply the EM algorithm to the data. Evaluate the convergence of the algorithm
#    by testing several initial values of θ^0. Report estimates for m0 and m1, and display
#    p(z | y) by using the barplot command to visualize the probability matrix of size
#    n x 2.
#
# 10. Install the R package mclust from the CRAN web site. Look at the different options
#     of the Mclust function (models and outputs), and run the Mclust command on the
#     data for G = 1 to 5.

# Turning off warnings
options(warn = -1)

# Libraries
library(mclust)

# Reading data
data2 <- read.table("C:/Pedja/Skola/Fakultet/M2/Computational Biology/data2.txt", quote="\"", comment.char="")
data2 <- data2[, 2]
n <- length(data2)

# Plotting histogram and density function
hist(data2, probability = T, col = 'lightblue', main = "", xlab = "Data")
lines(density(data2), add = T, lwd = 2, col = 'coral')

# We have 5 different parameters to estimate and it will be in the next form
#       theta = (m0, m1, sigma0^2, sigma1^2, p)

# Calculating probabilities p(zi = 1 | yi, θ)
prob <- function(theta)
{
     theta[5]/sqrt(theta[4])*exp(-1/(2*theta[4])*(data2 - theta[2])^2)/((1-theta[5])/sqrt(theta[3])*exp(-1/(2*theta[3])*(data2 - theta[1])^2) + theta[5]/sqrt(theta[4])*exp(-1/(2*theta[4])*(data2 - theta[2])^2))
}

# Initial values and counter for number of iterations
theta_new <- c(-1, 4, 1, 1, 0.1)
count <- 0

# EM algorithm
repeat{
        theta_old <- theta_new
        print(theta_old)
        p <- prob(theta_old)
        n1 <- sum(p)
        n0 <- n - n1
        theta_new1 <- 1/n0 * sum(data2*(1 - p))
        theta_new2 <- 1/n1 * sum(data2*p)
        theta_new3 <- 1/n0 * sum((1-p)*(data2 - theta_new1)^2)
        theta_new4 <- 1/n1 * sum(p*(data2 - theta_new2)^2)
        theta_new5 <- n1/n
        theta_new <- c(theta_new1, theta_new2, theta_new3, theta_new4, theta_new5)
        print(theta_new)
        count <- count + 1 
        if(prod(round(theta_old, 8) == round(theta_new, 8)))
        {
                print(count)
                print(theta_new)
                break
        }
}

# Estimated parameters
theta_EM <- theta_new

# Finding the posterior distribution for Z
z_dist <- prob(theta_EM)

# Plotting barplot
z_prob <- as.matrix(data.frame(z_dist, 1- z_dist))
barplot(t(z_prob))

# Let us see what will Mclust say about our data 
Mclust(data2)
# we have the model ("E", 2), which means that modelNames is "E", i.e. equal variance 
# in one dimension and G is 2, i.e. we have two classes.

Mclust(data2)$BIC
# we have three best models based on BIC criterion, which are ("E", 2), ("V", 2) and ("E", 3).

Mclust(data2)$df 
# We estimated 4 parameters m0, m1, sigma0^2 and p. Note sigma0^2 = sigma1^2 in this case (modelNames = "E")

# Here we can see estimated parameters
Mclust(data2)$parameters
# pro will tell us probabilities that data is in class 0 and class 1, respectively
# mean will give us m0 and m1
# variance$sigmasq will give us sigma0^2

# Let us compare Mclust with our EM algorithm
Mclust(data2, G = 2, modelNames = "V")$parameters
# We got this
# (m0, m1, sigma0^2, sigma1^2, p) = (-0.8417746, 2.3495063, 0.807073, 1.903795, 0.6250671);

# We can compere barplots. In order to do that we need to trnspose z given from Mclust
# and also to subtract from 1 in order to get probabilities for class 0 and 1.
barplot(1-t(Mclust(data2, G = 2, modelNames = "V")$z))

# Let us change G
Mclust(data2, G = 1)
Mclust(data2, G = 3)
Mclust(data2, G = 4)
Mclust(data2, G = 5)

# For G = 4 and 5 we can plot density functions
par(mfrow = c(1,2))
plot(Mclust(data2, G = 4), what = "density")
plot(Mclust(data2, G = 5), what = "density")
