# Hierarchical Conformal Prediction (HCP) for Clustered Data with Missing Responses

## Folder Overview
- **real data/** — unchanged
- **simulation/** — updated (added oracle method comparison)


## Abstract
Our code investigates the performance of the HCP method across various aspects, including **marginal coverage** (`"HCP_marginal.R"`), **conditional coverage** (`"HCP_conditional.R"`), **local coverage**  (`"HCP_local.R"`.), and **simultaneous prediction** (`"HCP_simul.R"`). These codes are located in different folders for tables and figures.
We provide the code for all tables and figures in the paper’s simulation and real data analysis, each with a corresponding folder containing a `"main.R"` file that calls all necessary scripts and can be run directly to reproduce the results.



## Instructions for Use
First, make sure to install the required R packages by running the following command:
```R
install.packages(c("MASS", "stats", "grf", "quantreg", "doParallel", "doRNG", "lme4", "merTools", "randomForest", "rstudioapi", "ggplot2", "gridExtra", "viridis"))
```

Then, to reproduce the program results, simply **download** the entire `"simulation"` and `"real data"` folders.
Our program can automatically set the **paths**, so no manual setup is required.
However, to ensure the paths are recognized correctly, please note the following:
1. The code should be executed in the **RStudio** environment.
2. The original **folder structure** must remain unchanged.


If the above two conditions are not met, the file **paths** in each `.R` file must be manually configured, including specifying paths for function calls, loading input data, and saving results; 
see the **Note** at the end for details. 

Once the paths are correctly set, simply click on the `"main.R"` file in each folder to run the program directly. 
The results will be automatically saved in the corresponding folder as `"final_results_reproduce.csv"`. 
Our original results, stored as `"final_results.csv"`, are in the same folder and can be used for comparison before and after reproduction.
See the **Example**  below for details.  

## Example
To reproduce the results in **Table 1** of the paper, navigate to [this repository](https://github.com/judywangstat/HCP.git), click the green **"Code"** button in the top right corner, and select **"Download ZIP"**.
After downloading the ZIP file, extract it, navigate to the `"simulation"` folder, and enter the `"Table 1"` directory.
Table 1 can be reproduced in either of two ways:
1. **Generate the entire Table 1** – A straightforward and convenient approach that requires no manual setup. The steps are as follows:
     - Open the `"main.R"` file in the `"Table 1"` folder.  
     - Run the script.  
     - The complete results for Table 1 will be saved as `"final_results_reproduce.csv"` in the `"Table 1"` folder.
     
2. **Compute each value in Table 1 individually** – A faster approach, useful for a quick verification.  The steps are as follows:
   - Open a single master code file in the `"Table 1"` folder, such as the HCP method (`"HCP_marginal.R"`), or other methods like `"DWR.R"`, `"LC.R"`, or `"LMEM.R"`.
   - Use the following code to clear the workspace and free memory:
     ```R
     rm(list = ls())  
     gc()
     ```  
   - Set the **sample size** (n) and **scenario** (scenario), for example:  
     ```R
     n <- 100  
     scenario <- "Homo"
     ```  
   - Run the script.  
   - The results will be displayed directly in the console.



## Note 
All scripts can automatically set **paths**, but if any issues occur, the paths need to be configured manually.  For example, in **Table 1**, the paths need to be set manually in the following locations:  

1.  In `"main.R"`, modify the paths for calling methods and saving results as follows:

     - Replace relative paths in `"source()"` with specified paths.  

       Modify the following code:  
        ```R
        source(file.path(main_dir, "HCP_marginal.R"), local = TRUE, chdir = TRUE)
        source(file.path(main_dir, "DWR.R"), local = TRUE, chdir = TRUE)
        source(file.path(main_dir, "LC.R"), local = TRUE, chdir = TRUE)
        source(file.path(main_dir, "LMEM.R"), local = TRUE, chdir = TRUE)
        ```
       to a specific path, for example:  
        ```R
        source("/Users/judy/Desktop/code/Table 1/HCP_marginal.R", local = TRUE, chdir = TRUE)
        source("/Users/judy/Desktop/code/Table 1/DWR.R", local = TRUE, chdir = TRUE)
        source("/Users/judy/Desktop/code/Table 1/LC.R", local = TRUE, chdir = TRUE)
        source("/Users/judy/Desktop/code/Table 1/LMEM.R", local = TRUE, chdir = TRUE)
        ```

      -  Modify the path for saving results:
   
         Modify the following code:  
         ```R
         write.csv(final_results, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = FALSE)
         ```
         to a specific path, for example:  
         ```R
         write.csv(final_results, file = "/Users/judy/Desktop/code/Table 1/final_results_reproduce.csv", row.names = TRUE)
         ```
2. In `"HCP_marginal.R"`, `"DWR.R"`, `"LC.R"`, and `"LMEM.R"`, modify the paths for calling simple functions and passing them into parallel computation as follows:  

    - Replace relative paths in `"source()"` with specified paths.  

      Modify the following code:  
        ```R
        source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)
        clusterExport(cl, "parent_dir")
          clusterEvalQ(cl, {
            library(MASS)
            library(stats)
            library(grf)
            library(quantreg)
            source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)
          })
        ```
      to a specific path, for example:  
        ```R
        source("/Users/judy/Desktop/code/my_functions.R")
        clusterEvalQ(cl, {
            library(MASS)
            library(stats)
            library(grf)
            library(quantreg)
            source("/Users/judy/Desktop/code/my_functions.R")
          })
        ```

