import pandas as pd
import numpy as np
from scipy.stats import pearsonr

# Function to read scores from Excel file
def read_scores_from_excel(file_path):
    df = pd.read_excel(file_path, sheet_name='bootstrapping', index_col=0)  # Assuming the first column has row labels (e.g., sequence positions)
    return df

# Function to calculate Pearson correlation coefficient
def calculate_correlation(scores, lifespans):
    correlation, _ = pearsonr(scores, lifespans)
    return correlation

# Function to perform bootstrap resampling for correlation
def bootstrap_correlation(scores, lifespans, num_samples=1000):
    n = len(scores)
    bootstrap_indices = np.random.choice(np.arange(n), size=(num_samples, n), replace=True)
    bootstrap_correlations = []
    for indices in bootstrap_indices:
        bootstrap_scores = scores[indices]
        bootstrap_lifespans = lifespans[indices]
        correlation = calculate_correlation(bootstrap_scores, bootstrap_lifespans)
        bootstrap_correlations.append(correlation)
    return np.array(bootstrap_correlations)

# Main function for correlation analysis
def correlation_analysis(scores_file):
    # Read scores from Excel file
    df_scores = read_scores_from_excel(scores_file)
    
    # Extract sequence positions (assuming they are in the first column)
    sequence_positions = df_scores.index
    
    # Initialize a dictionary to store correlations for each lifespan column
    correlations = {}
    
    # Loop over each lifespan column (each column of scores)
    for lifespan_column in df_scores.columns:
        # Extract scores and lifespans
        scores = df_scores[lifespan_column].values
        lifespans = df_scores.index.astype(float)  # Convert sequence positions to float (adjust if necessary)
        
        # Calculate observed correlation
        observed_correlation = calculate_correlation(scores, lifespans)
        
        # Perform bootstrap resampling for correlation
        bootstrap_correlations = bootstrap_correlation(scores, lifespans)
        
        # Compute bootstrap statistics
        bootstrap_mean = np.mean(bootstrap_correlations)
        bootstrap_std = np.std(bootstrap_correlations)
        
        # Store results in the dictionary
        correlations[lifespan_column] = {
            'observed_correlation': observed_correlation,
            'bootstrap_mean': bootstrap_mean,
            'bootstrap_std': bootstrap_std
        }
    
    return correlations

# Example usage
if __name__ == "__main__":
    scores_file = '/Users/romaniosbourne/Desktop/p53_cross-org_study/large seqs/~$primates_new.xlsx'
    
    correlations = correlation_analysis(scores_file)
    
    # Print results for each lifespan column
    for lifespan_column, stats in correlations.items():
        print(f"Lifespan column: {lifespan_column}")
        print(f"  Observed correlation coefficient: {stats['observed_correlation']:.2f}")
        print(f"  Bootstrap mean correlation coefficient: {stats['bootstrap_mean']:.2f}")
        print(f"  Bootstrap standard deviation: {stats['bootstrap_std']:.2f}")
