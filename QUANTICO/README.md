# QUANTICO Reproducibility Instructions

1. Simulation study
* The QUANTICO simulation study is performed for two sample size scenarios, n = 100 and n = 200. Within each scenario, multiple sub-scenarios are considered. To reproduce the results provided in the Table 1, please go to the folders ‘Sample Size 100 Codes’ and ‘Sample Size 200 Codes’ and follow the instructions in the corresponding README files.
* In order to reproduce the plots in Figure 3, please go to the folder ‘SIMULATION_PLOTS’ and follow the instructions in the corresponding README file.
* In order to reproduce the results on the coverage of uniform credible interval in QUANTICO (noted in Table S1), please go to the folder ‘Coverage’ and follow the instructions in the corresponding README file.
* In order to recalculate computation time required for QUANTICO (reported in Table S1), please go to the folder ‘COMPUTATION_TIMES’ and follow the instructions in the corresponding README file.

2. Real data analysis
* We provide the QUANTICO analysis of a dataset which is produced to emulate the real dataset. We produce similar plots based on QUANTICO analysis of this synthetic dataset.
* To produce similar plots as in Figure 4, 5(a), S1 and S2, please go to the ‘SYNTHETIC REAL DATA ANALYSIS’ folder, and run ‘QUANTICO_SYNTHETIC_DATA_ANALYSIS.m’ followed by ‘QUANTICO_SYNTHETIC_DATA_ANALYSIS_PLOTS.R’.
* To produce similar plots as shown in Figure 5(b-d), please go to ‘SYNTHETIC REAL DATA ANALYSIS’ and run ‘QUANTICO_SYNTHETIC_outlier_patient_mutations.R’.
