# This work is based on the published work "Dislocation–grain boundary interactions in Ta: numerical, molecular dynamics, and machine learning approaches". Cite: https://doi.org/10.1007/s10853-023-09167-y

# The dislocation-GB interaction was studied using numerical calculations of slip transmission parameters (STP) and molecular dynamics simulations (MDS). The parameters from both STP and MDS are fed as inputs to machine learning (ML) models.

# Here, you can find three folders with each containing STP, MDS, and ML. 
1. Open 'STP.m' using matlab with other *.m files in the same folder. The code Stabix (https://stabix.readthedocs.io/en/latest/) was used. MTEX is required. 
2. The input script for GB generation is 'gb.in' and shearing is 'in.shear_controlled' (credits to Wurong Jian, wurong@ucsb.edu)
3. Open the file 'ML.py' for ML with data_ML.xlsx 

# Mail id: kedharnath1992@gmail.com 
# If you use these scripts, kindly cite "Kedharnath, A., Kapoor, R. & Sarkar, A. Dislocation–grain boundary interactions in Ta: numerical, molecular dynamics, and machine learning approaches. J Mater Sci 59, 243–257 (2024). https://doi.org/10.1007/s10853-023-09167-y".
