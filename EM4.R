# 5. Write an EM algorithm for estimating all allele frequencies.
# 6. Application: We observe n = 521 cases of peptic ulcer disease, for which nA = 186,
#    nB = 38, nAB = 13 and nO = 284. Find estimates for pA, pB and pO.

# Data
n <- 521
nA <- 186
nB <- 38
nAB <- 13
nO <- 284

# Function for calculating nAA, nAO, nBB and nBO
find_n <- function(theta)
{
     nA_A <- nA*theta[1]^2/(theta[1]^2 + 2*theta[1]*theta[3])
     nA_O <- nA - nA_A
     nB_B <- nB*theta[2]^2/(theta[2]^2 + 2*theta[2]*theta[3])
     nB_O <- nB - nB_B 
     n_new <- c(nA_A, nA_O, nB_B, nB_O)
}

# Initial values and counter for number of iterations
theta_new <- c(1, 1, 1)
count <- 0

#EM Algorithm
repeat{
        theta_old <- theta_new
        print(theta_old)
        n_new <- find_n(theta_old)
        theta_new1 <- 1/(2*n) * (2*n_new[1] + n_new[2] + nAB)
        theta_new2 <- 1/(2*n) * (2*n_new[3] + n_new[4] + nAB)
        theta_new3 <- 1/(2*n) * (n_new[2] + n_new[4] + 2*nO)
        theta_new <- c(theta_new1, theta_new2, theta_new3)
        print(theta_new)
        count <- count + 1 
        if(prod(round(theta_old, 8) == round(theta_new, 8)))
        {
                 print(count)
                print(theta_new)
                break
        }
}

