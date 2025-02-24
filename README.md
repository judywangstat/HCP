# Hierarchical Conformal Prediction for Clustered Data with Missing Responses


## Code
Our code investigates the performance of the HCP method across various aspects, including **marginal coverage** (`HCP_marginal.R`), **conditional coverage** (`HCP_conditional.R`), **local coverage**  (`HCP_local.R`.), and **simultaneous prediction** (`HCP_simul.R`). These codes are located in different folders for tables and figures.
We provide the code for all tables and figures in the paper’s simulation and real data analysis, each with a corresponding folder containing a `main.R` file that calls all necessary scripts and can be run directly to reproduce the results.



## Instructions for Use
First, make sure to install the required R packages by running the following command:
```R
install.packages(c("MASS", "stats", "grf", "quantreg", "doParallel", "doRNG", "lme4", "merTools", "randomForest", "rstudioapi"))
```

To reproduce the program results, simply **download** the entire `simulation` and `real data` folders.
Our program can automatically set the **paths**, so no manual setup is required.
However, to ensure the paths are recognized correctly, please note the following:
1. The code should be executed in the **RStudio** environment.
2. The original **folder structure** must remain unchanged.


If the above two conditions are not met, the file **paths** in each `.R` file must be manually configured, including specifying paths for function calls, loading input data, and saving results.


Once the paths are correctly set, simply click on the `main.R` file in each folder to run the program directly. The results will be automatically saved in the corresponding folder.

## Example
To reproduce the results in **Table 1** of the paper, navigate to [this repository](https://github.com/judywangstat/HCP.git), click the green **"Code"** button in the top right corner, and select **"Download ZIP"**.
After downloading the ZIP file, extract it, navigate to the `simulation` folder, and enter the `Table 1` directory.
Table 1 can be reproduced in either of two ways:
1. **Generate the entire Table 1** – A straightforward and convenient approach that requires no manual setup. The steps are as follows:
     - Open the `main.R` file in the `Table 1` folder.  
     - Run the script.  
     - The complete results for Table 1 will be saved as `"final_results_reproduce.csv"` in the `Table 1` folder.
     
   

## Example
For a single master code file implementing the proposed method, simply run 
