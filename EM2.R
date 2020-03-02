# 2. Write a computer program for drawing samples from p(yi | θ) of size n = 200
#    (for fixed values of (p; m0; m1)). Check if your program is correct by drawing a
#    histogram of the simulated data.

# Turning off the warnings
options(warn = -1)

# Function for generating data from mixture normal distribution
rnormix <- function(n, pi, mu, sigma)
{
     k <- length(pi)
     Z <- sample.int(k, n, replace = TRUE, prob = pi)
     Y <- rnorm(n, mu[Z], sigma[Z])
     return(list(Y, Z))
}

# Fixed values
n <- 200
m0 <- 3
m1 <- 7
mu <- c(m0, m1)
p <- 0.27
pi <- c(p, 1-p)
sigma <- c(1,1)

Mixture <- rnormix(n, pi, mu, sigma)
Data <- c(Mixture[[1]])
# R will assign data into classes 1 and 2 so we need to subtract 1 
Class <- Mixture[[2]] - 1

# Plotting histogram and dansity function
hist(Data, probability = T, col = 'lightblue', main = "", xlab = "Data")
lines(density(Data), add = T, lwd = 2, col = 'coral')

# Calculating probabilities p(zi = 1 | yi, θ)
prob <- function(theta)
{
     exp(-1/2*(Data - theta[2])^2)/(exp(-1/2*(Data - theta[1])^2) + exp(-1/2*(Data - theta[2])^2))
}

# Inital values and counter for number of iterations
theta_new <- c(1, 7)
count <- 0

# EM algorithm
repeat{
          theta_old <- theta_new
          print(theta_old)
          p <- prob(theta_old)
          n1 <- sum(p)
          n0 <- n - n1
          theta_new1 <- 1/n0 * sum(Data*(1 - p))
          theta_new2 <- 1/n1 * sum(Data*p)
          theta_new <- c(theta_new1, theta_new2)
          print(theta_new)
          count <- count + 1 
          if(prod(round(theta_old, 8) == round(theta_new, 8)))
          {
               print(count)
               print(theta_new)
               break
          }
}


