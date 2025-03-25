#set.seed(42)  # Ensure reproducibility

# Function to compute \hat{T}_n
compute_T_hat_n = function(Z, Xi, mu_0) {
  n = nrow(Z)
  p = ncol(Z)
  
  # Compute T1
  T1 = sum(sapply(1:n, function(i) {
    sum(sapply((1:n)[-i], function(j) {
      sum(Z[i, ] * Z[j, ])
    }))
  })) / (n * (n - 1))
  
  # Compute T2
  T2 = sum(sapply(1:n, function(i) {
    sum(sapply((1:n)[-i], function(j) {
      sum(Z[i, ] * (Xi[j, ] * mu_0))
    }))
  })) / (n * (n - 1))
  
  # Compute T3
  T3 = sum(sapply(1:n, function(i) {
    sum(sapply((1:n)[-i], function(j) {
      sum(mu_0 * (Xi[i, ] * Xi[j, ] * mu_0))
    }))
  })) / (n * (n - 1))
  
  # Compute T_n
  T_n = T1 - 2 * T2 + T3
  
  # Compute sigma_n,0^2
  sigma_n_0_sq = sum(sapply(1:n, function(i) {
    sum(sapply((1:n)[-i], function(j) {
      z_i <- Z[i, ] - (Xi[i, ] * mu_0)
      z_j <- Z[j, ] - (Xi[j, ] * mu_0)
      sum(z_i * z_j)^2
    }))
  })) * (2 /(n * n * n * (n - 1)))
  
  # Compute {T}_n hat
  T_hat_n = T_n / sqrt(sigma_n_0_sq)
  
  return(T_hat_n)
}

# Function to estimate size and power
estimate_size_power = function(n, p, alpha, beta, num_simulations, significance_level) {
  count_H0 = 0  # Counter for rejected null under H0
  count_H1 = 0  # Counter for rejected null under H1
  
  mu_0 = c(rep(2, p))  # Null mean
  d = c(runif(p, -0.5, 0.5) ) # Generate alternative mean differences
 # mu_alt = mu_0 + alpha * d  # Alternative mean
  
  for (sim in 1:num_simulations) {
    # Generate missing pattern indicators
    Xi = matrix(rbinom(n * p, 1, 1 - beta), n, p)  # Missing pattern
    
    # Generate data under H0 (Null)
    Z_H0 = matrix(rgamma(n * p, shape = 4, scale = 0.5), n, p) * Xi
    T_H0 = compute_T_hat_n(Z_H0, Xi, mu_0)
    
    # Generate data under H1 (Alternative)
    Z_H1 = matrix(rgamma(n * p, shape = 4, scale = 0.5), n, p) * Xi
    for (j in 1:p) {
      Z_H1[, j] = Z_H1[, j] + (alpha/2) * d[j]  # Shift mean under H1
    }
   T_H1 = compute_T_hat_n(Z_H1, Xi, mu_0)
    
    # Compute critical value for rejection (assuming normality)
    critical_value = qnorm(1 - significance_level/2)
    
    # Count rejections
    if (abs(T_H0) > critical_value) 
      count_H0 = count_H0 + 1
   if (abs(T_H1) > critical_value) 
     count_H1 =  count_H1 + 1
  }
  
  # Compute empirical size and power
  size = count_H0/num_simulations 
  power = count_H1/num_simulations

  
  return(list( size = size , power = power, count_H0= count_H0, count_H1 = count_H1, T_H0 = T_H0, T_H1 = T_H1))
}

# Parameters
n = 100  # Sample size
p = 50  # Dimension
alpha = 0.1  # Mean difference factor
beta = runif(p)   # Missing probability
num_simulations = 100  # Number of simulations
significance_level = 0.05  # 5% significance level


# Run the test
result = estimate_size_power(n, p, alpha, beta, num_simulations, significance_level)
#cat(sprintf("Value of test statitics : %.3f\n", result$T_H0))
#cat(sprintf("Value of test statitics : %.3f\n", result$T_H1))
#cat(sprintf("No. of count when H0 are rejected :%.0f\n", result$count_H0))
#cat(sprintf("No. of count when H0 are rejected :%.0f\n", result$count_H1))
#cat(sprintf("Empirical Size: %.3f\n", result$size))
cat(sprintf("Empirical Power: %.3f\n", result$power))