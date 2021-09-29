import numpy as np 
import matplotlib.pyplot as plt 
import tools as tools 
import time
start_time = time.time()

filein = '../core-softened-phase-diagram/CSA/temp-0.01/lambda_0.25/run/lammps_temp_0.01_press_6.5.snap'

r = tools.read_in_trajectory(filein, 4000, 100)

L1 = 1.1*np.absolute(np.max(r[:,:,0]) - np.min(r[:,:,0])) 
L2 = 1.1*np.absolute(np.max(r[:,:,1]) - np.min(r[:,:,1])) 
L3 = 0.0
rdf, r = tools.rdf_monomer_to_monomer(r[:,:,0],r[:,:,1],L1,L2)

sex, csum, r_s2 = tools.excess_entropy_from_rdf(rdf,r)
r_s2 = r_s2[r_s2>0]
csum = csum[csum>0]

rho = (0.5 * 4000) / (L1 * L2)

print(csum)
print(r_s2)
print("--- %s seconds ---" % (time.time() - start_time))

plt.figure()
#plt.ylim(0,12)
data = np.genfromtxt("../core-softened-phase-diagram/CSA/temp-0.01/lambda_0.25/ex_entropy/press_6.5_cumsex.dat")
plt.plot(r_s2,csum)
plt.plot(data[:,0],data[:,1],color="red")
#plt.yscale('log')
plt.show()

