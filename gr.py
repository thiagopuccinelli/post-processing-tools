import numpy as np 
import math 
import matplotlib.pyplot as plt 
from tools import excess_entropy_from_rdf,rdf_monomeric2d, rdf_monomer_to_monomer
import time
start_time = time.time()

##### read in lammps snap
N = int(np.genfromtxt('../core-softened-phase-diagram/CSA/temp-0.01/lambda_0.60/run/lammps_temp_0.01_press_8.5.snap', max_rows=1))
traj = np.genfromtxt(
    '../core-softened-phase-diagram/CSA/temp-0.01/lambda_0.60/run/lammps_temp_0.01_press_8.5.snap', skip_header=2,
    invalid_raise=False)[:, 1:4].reshape(-1, N, 3)

#print(rdf_monomer_to_monomer.__doc__)

snaps = len(traj[:,:,:])
L1 = 1.1*np.absolute(np.max(traj[0,:,0]) - np.min(traj[0,:,0])) 
L2 = 1.1*np.absolute(np.max(traj[0,:,1]) - np.min(traj[0,:,1])) 
L3 = 0.0
print(np.shape(traj[:,:,0]))
rdf, r = rdf_monomer_to_monomer(traj[:100,:,0],traj[:100,:,1],L1,L2)

rho = (0.5 * N) / (L1 * L2)

print(excess_entropy_from_rdf.__doc__)

sex, csum = excess_entropy_from_rdf(rdf,r)

#print(r)
#csum = csum[csum>0]
#print(csum)
print("--- %s seconds ---" % (time.time() - start_time))

plt.figure()
#plt.ylim(0,12)
data = np.genfromtxt("../core-softened-phase-diagram/CSA/temp-0.01/lambda_0.60/ex_entropy/press_8.5_cumsex.dat")
plt.loglog(r,csum)
plt.loglog(data[:,0],data[:,1],color="red")
plt.show()

#np.savetxt("teste.dat", np.column_stack([r,rdf]))




