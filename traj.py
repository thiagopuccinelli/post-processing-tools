import numpy as np 
import matplotlib.pyplot as plt 
from tools import read_in_trajectory, excess_entropy_from_rdf,rdf_monomeric2d, rdf_monomer_to_monomer
import time
start_time = time.time()

print(read_in_trajectory.__doc__)
filein = '../core-softened-phase-diagram/CSA/temp-0.01/lambda_0.60/run/lammps_temp_0.01_press_8.5.snap'


r = read_in_trajectory(filein, 4000, 100)
# print(r)
# print(np.shape(r)) 

L1 = 1.1*np.absolute(np.max(r[0,:,0]) - np.min(r[0,:,0])) 
L2 = 1.1*np.absolute(np.max(r[0,:,1]) - np.min(r[0,:,1])) 
L3 = 0.0
rdf, r = rdf_monomer_to_monomer(r[:,:,0],r[:,:,1],L1,L2)

sex, csum = excess_entropy_from_rdf(rdf,r)

rho = (0.5 * 4000) / (L1 * L2)

print(sex * rho )

print("--- %s seconds ---" % (time.time() - start_time))

plt.figure()
#plt.ylim(0,12)
data = np.genfromtxt("../core-softened-phase-diagram/CSA/temp-0.01/lambda_0.60/ex_entropy/press_8.5_cumsex.dat")
plt.loglog(r,csum)
plt.loglog(data[:,0],data[:,1],color="red")
plt.show()

