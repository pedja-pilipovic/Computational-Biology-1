# Download the data from the following URL:
#
#         http://membres-timc.imag.fr/Olivier.Francois/sequence_2019.txt
#
# The data consists of a sequence of 399 binary items. The objective of the challenge
# is to provide a list of (one or more) change points with the following information
#
#     • Most likely change point position, z, in the range [1, 399].
#     • Lower and upper values zl and zu, such that p(z ∈ [zl, zu]) = 0.75.
#     • Estimates of frequencies θ1 and θ2 before and after z.
#     • Number of iterations of the EM algorithm.

# Turning off warnings
options(warn = -1)

# Needed libreries
library(bazar)
library(ggplot2)

# Reading data
sequence_2019 <- read.table("C:/Pedja/Skola/Fakultet/M2/Computational Biology/sequence_2019.txt", quote="\"", comment.char="")
sequence_2019 <- as.vector(t(as.matrix(sequence_2019)))
sequence_2019 <- sequence_2019[-400]
n <- length(sequence_2019)

# Function for calculating R
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

# Function for calculating Phi
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

# Function for calculating expacted values
Exp <- function(theta, y)
{
        n <- length(y)
        E1 <- E2 <- E3 <- E4 <- 0
        r <- R(theta, y)
        Rsum <- sum(r)
        phi1 <- Phi(y)
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


theta_seq_new <- c(0.1,0.2)

# EM algorithm
EM <- function(theta_seq_new, data)
{
     count <- 0
     repeat{
          theta_seq_old <- theta_seq_new
          E <- Exp(theta_seq_old, data)
          theta_seq_new1 <- E[1]/(E[1] + E[2])
          theta_seq_new2 <- E[3]/(E[3] + E[4])
          theta_seq_new <- c(theta_seq_new1, theta_seq_new2)
          #print(theta_seq_new)
          count <- count + 1
          if(prod(round(theta_seq_old, 8) == round(theta_seq_new, 8)))
          {
               #print(count)
               #print(theta_seq_new)
               break
          }
     }
     return(c(theta_seq_new, count))
}

# Distribution of Z and ploting a histogram
z <- function(theta, y)
{
     n <- length(y)
     prob <- 0
     r <- R(theta, y)
     Rsum <- sum(r)
     prob <- r/Rsum
     plot(prob, type = 'h', xlim = c(0, 400), xlab = "z")
     return(prob)
}

theta_in1 <- c(0.1, 0.8)

# First call of EM algorithm with original data
EM_alg <- EM(theta_in1, sequence_2019)
theta_EM <- EM_alg[1:2]
iter <- EM_alg[3]
p <- z(theta_EM, sequence_2019)
# Most probable change point Z
z_max <- which.max(p)
# z_lower and z_upper such that P(Z in [z_lower, z_upper]) = 0.75
zl <- max(which(cumsum(p)<0.125))
zu <- min(which(cumsum(p)>0.875))

# Puting result in the data frame in wanted format
res <- data.frame(V1 = z_max, V2 = zl, V3 = zu, V4 = theta_EM[1], V5 = theta_EM[2], V6 = iter)
colnames(res) <- c("position", "lower","upper", "theta1", "theta2", "iter")

# Data frame for plotitng CDF
prob <- data.frame(prob = p)
prob

# Plotting CDF
ggplot(prob, aes(x = 1:399, y=cumsum(prob))) + geom_line() + 
        geom_point() + geom_hline(yintercept = 0.125) + geom_hline(yintercept = 0.875) +
        labs(x = "Z", y = "CDF") 

# Creating subsequence such that we will keep data with probabilites greater than zero
subsequence <- sequence_2019[!almost.zero(p)]

# Performing the same proces for creating subsequences and getting the result 
result <- function(theta_in, data)
{
     for(i in 2:6)
     {
          EM_alg <- EM(theta_in, data)
          theta_EM <- EM_alg[1:2]
          iter <- EM_alg[3]
          p <- z(theta_EM, data)
          z_max <- which.max(p)
          zl <- max(which(cumsum(p)<0.125))
          zu <- min(which(cumsum(p)>0.875))
          res <- rbind(res, c(z_max, zl,  zu, theta_EM[1], theta_EM[2], iter))
          data <- data[!almost.zero(p)]
     }
     return(res)
}


# Result
res <- result(theta_in1, subsequence)

theta_in2 <- c(0.1, 0.2)

# Adding a new solution from the EM algorithm 
EM_alg2 <- EM(theta_in2, sequence_2019)
theta_EM2 <- EM_alg2[1:2]
iter2 <- EM_alg2[3]
p2 <- z(theta_EM2, sequence_2019)
z_max2 <- which.max(p2)
zl2 <- max(which(cumsum(p2)<0.125))
zu2 <- min(which(cumsum(p2)>0.875))
res <- rbind(res, c(z_max2, zl2,  zu2, theta_EM2[1], theta_EM2[2], iter2))
subsequence2 <- sequence_2019[!almost.zero(p2)]

res <- result(theta_in2, subsequence2)

result <- cbind(1:6, unique(res))
colnames(result)[1] <- "number"

