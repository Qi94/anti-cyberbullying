# Game-theoretic modeling and analysis of cyberbullying spreading on OSNs

## 01_raw 
This folder includes the raw data used in this paper (1k nodes, 5k nodes, and 10k nodes), downloaded from 'https://networkrepository.com/'. 
Please feel free to use 'changeFileFormat.py' to convert the source files' format for MATLAB (version: 2023Rb) for your convenience. 
The source data at 100 nodes scale are provided by Yang Qin in '.mat' format initially.

## 02_code 
The subfolder '00_SCR Optimization' includes the implementation of Algorithm 1, SCR Optimization ('SCR_Optimization.m'), along with the corresponding source data in '.mat' format, located in 'asset'. These files are essential for generating optimal strategy profiles. 

To achieve convergence and obtain the relevant strategy profiles for further analysis in the '01_Comparative Experiments' and '02_Sensitivity_Analysis' subfolders, please ensure that the source input file path has been updated to reflect your  network before executing the scripts to get the result. The total execution time depends on your hardware specifications and the network scale. For your reference, it will take approximately 30 minutes for a MacBook Pro with an Apple Silicon M2 Max and 32 GB RAM to reach convergence in a network of 10,000 nodes scale.

In the '01_Comparative Experiments' subfolder, you will find all the scripts necessary for conducting the comparative experiments detailed in Section 5 of the research documentation. The '02_Sensitivity_Analysis' subfolder contains the code used in Section 6. These scripts should be compatible with all network scales mentioned in the paper, provided that the correct source data file path has been configured by the user.

This folder primarily uses a network of 10,000 nodes as the sample input for research results. Before running the script 'SCR_Optimization.m', please ensure that the source input file path has been updated to reflect your data.