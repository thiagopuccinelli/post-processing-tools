import numpy as np 
import matplotlib.pyplot as plt 
from post_processing import * 

msd, t = msd_evaluation(200,"snap_LAMMPS_2d.xyz",1000,20,10)


data=np.loadtxt("lammps_msd.dat")

plt.figure()


plt.plot(t,msd,ls='-', label="from snapshot")
plt.plot(data[:,0]*0.01,data[:,1], ls ='-', label="from LAMMPS")

plt.legend(loc="best")

plt.show()