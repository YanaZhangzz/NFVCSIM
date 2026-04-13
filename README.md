# MNFSM Statistical Analysis

This repository provides a collection of R scripts for statistical modeling, focusing on **MNFSM** methodology. It includes core functions for parameter estimation and predictive modeling.

## 📁 Repository Structure

### Core Methodology
- **`JOE_estimation.R`**: Implementation of core functions for parameter estimation.
- **`JOE_grad_hessian.R`**: Functions for computing Gradients and Hessian matrices required for optimization and inference.
- **`JOE_inference.R`**: Scripts for statistical inference.
- **`MNFSM_functions.R`**: The main functional library for the MNFSM framework.
- **`MNFSM_pred_functions.R`**: Dedicated functions for out-of-sample prediction and validation.
- **`net_type_functions.R`**: Auxiliary functions for handling network-structured data or specific net-type configurations.

### Demos & Examples
- **`MNFSM_esti_demo.R`**: A script to execute the complete parameter estimation process.
- **`MNFSM_pred_demo.R`**: An example script showing how to perform predictions using trained models.
