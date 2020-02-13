# post-processing-tools

This repository has tools for analysing Molecular Dynamics data. 

# Step to implement fortran 90 codes as library for python code:

To implement fortran 90 code, one should follow the bellow commands line in terminal:


# Radial distribution function: 

In python3: 
```
f2py3 -c lib/radial_distribution_function.f90 -m radial_distribution_evaluation
```