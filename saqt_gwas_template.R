
# Script for running the saqt_gwas_v1.R code. This describes the process for using the code to explore and visualize the results from the SAQT-GWAS Pipeline workflow.

# Load the saqt_gwas_v1.R code. Be sure to change the path to reflect the location of your copy of saqt_gwas_v1.R
source("/path/to/saqt_gwas_v1.R", print.eval = TRUE)

# Specify the location of the association results directory. Change results_name to name that describes the data in a useful way. 
results_name <- "/path/to/results_directory"

# Initialize the required libraries. 
initialize()


# Analysis function, for example:
linear_assoc_manhattan_plot(results_name, "Measure-Name")

