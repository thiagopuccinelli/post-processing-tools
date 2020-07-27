import numpy as np 
import matplotlib.pyplot as plt 
from post_processing import * 

rdf, r = radial_distribution_2d(200,"snap_LAMMPS_2d.xyz",1000,20)



plt.figure()


plt.plot(r,rdf, ls='-')


plt.show()