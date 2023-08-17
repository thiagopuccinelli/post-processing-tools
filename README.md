# post-processing-tools

This repository has tools for analysing Molecular Dynamics data. Some of the tools are written in Fortran 90 language
and can be accessed by f2py3 from the class post_processing.py found in the python directory. 

# Step to implement fortran 90 codes as library for python code:

To implement fortran 90 code, one should follow the bellow commands line in terminal inside the python directory:

```
f2py3 -c ../lib/* -m tools
```

# To access the post-processing-tools from its python class: 

In order to access the tools from its class, one should insert in a python code the following lines: 

```
import sys
pathto = sys.argv[1]
sys.path.append(pathto)
from post_processing import *
```

With that said, to execute a python code that will take into account some of the computations, the user needs to 
point out the location of the python directory,

```
python3 PYTHONSCRIPT.py location-to-this-github-directory/python/
```

This way the user can, for instance, read in the trajectory from a MD simulation: 

```
trajectory = post_processing.read_trajectory("TRAJECTORY_FILE",NPART,NSTEPS)
```

In the above command, NPART states the number of particles in the system, NSTEPS is the number of steps in the simulation and "TRAJECTORY_FILE" is the name of the trajectory file. 

