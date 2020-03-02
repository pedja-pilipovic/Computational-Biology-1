# 7. Generate simulated data for known values of θ and z. Apply the EM algorithm
#    to the simulated data, and evaluate the convergence of the algorithm by testing
#    several values of θ^0. Plot histograms for p(z|y) using plot(. ,type =’h’).

# Turning off the warnings
options(warn = -1)

# Ganerating data
set.seed(1)
n <- 500
theta1 <- 1/3
theta2 <- 3/4 
z <- 134
theta <- c(theta1, theta2)
y1 <- rbinom(z-1, 1, theta1)
y2 <- rbinom(n-z+1, 1, theta2)
#y1 <- sample(c(0,1), z, replace = T, prob = c(1 - theta1, theta1))
#y2 <- sample(c(0,1), n-z, replace = T, prob = c(1 - theta2, theta2))
y <- c(y1, y2)

# Function for calucalting R
R <- function(theta, y)
{
        n <- length(y)
        Rz <- 1:n
        for (z in 2:n)
        {
                Rz[z] <- Rz[z-1] * (theta[1]/theta[2])^y[z-1]*((1 - theta[1])/(1 - theta[2]))^(1 - y[z-1])
        }
        return(Rz)
}

# FUnction for calculating Phi
Phi <- function(y)
{
        n <- length(y)
        phi1 <- rep(0, n)
        for (z in 2:n) 
        {
                phi1[z] <- phi1[z-1] + y[z]
        }
        return(phi1)
}

# Calculating expacted values
Exp <- function(theta, y)
{
     n <- length(y)
     E1 <- E2 <- E3 <- E4 <- 0
     r <- R(theta, y)
     Rsum <- sum(r)
     phi <- Phi(y)
     phi1 <- phi
     phi3 <- sum(y) - phi1
     for(z in 1:n)
     {          
          E1 <- E1 + r[z]*phi1[z]
          E2 <- E2 + r[z]*(z - 1 - phi1[z])
          E3 <- E3 + r[z]*phi3[z]
          E4 <- E4 + r[z]*(n - z + 1 - phi3[z])
     }
     E <- c(E1, E2, E3, E4)
     E/Rsum
}

# Initial values for EM algorithm and counter for iterations
theta_new <- c(0.99, 0.75)
count <- 0

# EM algorithm
repeat{
        theta_old <- theta_new
        print(theta_old)
        E <- Exp(theta_old, y)
        theta_new1 <- E[1]/(E[1] + E[2])
        theta_new2 <- E[3]/(E[3] + E[4])
        theta_new <- c(theta_new1, theta_new2)
        count <- count + 1
        if(prod(round(theta_old, 8) == round(theta_new, 8)))
        {
                print(count)
                print(theta_new)
                break
        }
}

# Estimated theta
theta_EM <- theta_new

# Plotting distribution for change point Z
prob <- 0
r <- R(theta_EM, y)
Rsum <- sum(r)
prob <- r/Rsum
plot(prob, type = 'h', xlim = c(100, 150), xlab = "z")

