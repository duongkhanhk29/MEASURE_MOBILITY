
Dt <- read.csv("GDIM_2023_03.csv", stringsAsFactors=TRUE)


# Creating a new data frame from Dt
Q_df <- data.frame(
  C1P1 = Dt$tm11 * Dt$P1, C2P1 = Dt$tm12 * Dt$P1, C3P1 = Dt$tm13 * Dt$P1, C4P1 = Dt$tm14 * Dt$P1, C5P1 = Dt$tm15 * Dt$P1,
  C1P2 = Dt$tm21 * Dt$P2, C2P2 = Dt$tm22 * Dt$P2, C3P2 = Dt$tm23 * Dt$P2, C4P2 = Dt$tm24 * Dt$P2, C5P2 = Dt$tm25 * Dt$P2,
  C1P3 = Dt$tm31 * Dt$P3, C2P3 = Dt$tm32 * Dt$P3, C3P3 = Dt$tm33 * Dt$P3, C4P3 = Dt$tm34 * Dt$P3, C5P3 = Dt$tm35 * Dt$P3,
  C1P4 = Dt$tm41 * Dt$P4, C2P4 = Dt$tm42 * Dt$P4, C3P4 = Dt$tm43 * Dt$P4, C4P4 = Dt$tm44 * Dt$P4, C5P4 = Dt$tm45 * Dt$P4,
  C1P5 = Dt$tm51 * Dt$P5, C2P5 = Dt$tm52 * Dt$P5, C3P5 = Dt$tm53 * Dt$P5, C4P5 = Dt$tm54 * Dt$P5, C5P5 = Dt$tm55 * Dt$P5
)

# Creating a new data frame from Dt
P_df <- data.frame(
  C1P1 = Dt$P1,    C2P1 = Dt$P2/2,  C3P1 = Dt$P3/3,  C4P1 = Dt$P4/4,  C5P1 = Dt$P5/5,
  C1P2 = 0,        C2P2 = Dt$P2/2,  C3P2 = Dt$P3/3,  C4P2 = Dt$P4/4,  C5P2 = Dt$P5/5,
  C1P3 = 0,        C2P3 = 0,        C3P3 = Dt$P3/3,  C4P3 = Dt$P4/4,  C5P3 = Dt$P5/5,
  C1P4 = 0,        C2P4 = 0,        C3P4 = 0,        C4P4 = Dt$P4/4,  C5P4 = Dt$P5/5,
  C1P5 = 0,        C2P5 = 0,        C3P5 = 0,        C4P5 = 0,        C5P5 = Dt$P5/5
)

all(abs(rowSums(P_df) - 1) < 1e-6)
all(abs(rowSums(Q_df) - 1) < 1e-6)


calculate_TVD_variance <- function(Q_df, P_df, N) {
  # Initialize an empty data frame to store results
  results <- data.frame(TVD = numeric(nrow(Q_df)), 
                        TVD_variance = numeric(nrow(Q_df)),
                        non_contribution = numeric(nrow(Q_df)),
                        non_contribution_variance = numeric(nrow(Q_df)))
  
  # Identify the non_upward columns
  non_upward <- c("C1P1", "C1P2", "C2P2", "C1P3", "C2P3", "C3P3", 
                  "C1P4", "C2P4", "C3P4", "C4P4", "C1P5", "C2P5", 
                  "C3P5", "C4P5", "C5P5")
  
  # Loop through each row of Q_df and P_df
  for (i in 1:nrow(Q_df)) {
    # Extract the i-th row from both Q_df and P_df
    Q <- as.numeric(Q_df[i, ])  
    P <- as.numeric(P_df[i, ])  
    N_val <- N[i]  
    
    # Normalize Q and P by the sum of their elements
    Q_normalized <- Q / sum(Q)  
    P_normalized <- P / sum(P)  
    
    # Compute absolute differences
    abs_diff <- abs(P_normalized - Q_normalized)
    
    # Calculate Total Variation Distance (TVD)
    TVD <- 0.5 * sum(abs_diff)  
    
    # Calculate variance of TVD
    TVD_variance <- sum((P_normalized - Q_normalized)^2) / N_val  
    
    # Extract indices for non_upward elements
    non_upward_indices <- match(non_upward, names(Q_df))  
    
    # Extract absolute differences for non_upward elements
    abs_diff_non <- abs_diff[non_upward_indices]  
    
    # Compute diagonal contribution to TVD
    if (sum(abs_diff) == 0) {
      non_contribution <- 0  
      non_contribution_variance <- 0
    } else {
      non_contribution <- 0.5 * sum(abs_diff_non)  
      
      # Compute variance for non_upward contribution
      non_contribution_variance <- sum((P_normalized[non_upward_indices] - Q_normalized[non_upward_indices])^2) / N_val  
    }
    
    # Store results in the results data frame
    results$TVD[i] <- TVD
    results$TVD_variance[i] <- TVD_variance
    results$non_contribution[i] <- non_contribution
    results$non_contribution_variance[i] <- non_contribution_variance
  }
  
  # Return the results data frame
  return(results)
}


# Example usage:
N <- as.numeric(Dt$obs)  # Assuming the 'obs' column is numeric

# Call the function
TVD_results_df <- calculate_TVD_variance(Q_df, P_df, N)


meta_df <-cbind(TVD_results_df, Dt[, 1:15])
meta_df$year <- as.factor(meta_df$year)
meta_df$cohort <- as.factor(meta_df$cohort)



# Fit the model with weights and robust standard errors using lm()
library(sandwich)

mf <- lm(Social_Mobility ~ country + survey + year + status + cohort + parent + child, 
         data = meta_df, 
         weights = 1 / sqrt(Social_Mobility_Variance))

# Apply robust standard errors using vcovHC from sandwich
summary(mf, robust = TRUE)


# Load necessary libraries
library(visreg)
library(ggplot2)

visreg(mf, "country", gg=TRUE, fill.par=list(fill="#008DFF33")) + theme_minimal() + ylim(0,1) + coord_flip()



