# anti-cyberbullying

## 01_raw 
This folder includes the raw data used in this paper (1k nodes, 5k nodes, and 10k nodes), downloaded from 'https://networkrepository.com/'. 
Please feel free to use 'changeFileFormat.py' to convert the source files' format for MATLAB (version: 2023Rb) for your convenience. 
The source data at 100 nodes scale are provided by Yang Qin in '.mat' format initially.

## 02_code 
The subfolder '00_SCR Optimization' contains the implementation of Algorithm 1 SCR Optimization ('SCR_Optimization.m') and the corresponding source data in '.mat' format in 'asset'. These files are used to generate the optimal strategy profiles, please run the scripts under different network scales to reach convergence and get the corresponding strategy profiles for your further investigations in '01_Comparative Experiments' and '02_Sensitivity_Analysis' subfolder. The total execution duration is dependent on your hardware spec and network scale (around 30 mins for Apple Silicon M2 Max to reach the convergence under 10,000 nodes network).

In the '01_Comparative Experiments' subfolder, all the necessary scripts required for conducting comparative experiments as outlined in Section 5 of the research documentation are attached, the '02_Sensitivity_Analysis' subfolder contains the code used in Section 6. These scripts are adapted to all the network scales mentioned in the paper.

This folder primarily utilizes a network of 10,000 nodes as the instance input for research results, please change the source data file to adjust the input scale to get the related research result.