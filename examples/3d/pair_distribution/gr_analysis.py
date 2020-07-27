import numpy as np 
import matplotlib.pyplot as plt 
from post_processing import * 

rdf, r = radial_distribution_3d(500,"snap_LAMMPS_3d.xyz",1000,10)


plt.figure()


plt.plot(r,rdf, ls='-')


plt.show()