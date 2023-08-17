import numpy as np 
from post_processing import *
import matplotlib.pyplot as plt 

data = post_processing.read_trajectory("dump.lj",200,26)

print(np.shape(data))

Lx = 1.1*np.absolute(np.max(data[0,:,0]) - np.min(data[0,:,0]))
Ly = 1.1*np.absolute(np.max(data[0,:,1]) - np.min(data[0,:,1]))
rdf,r = post_processing.rdf_compute(data[:,:,0],data[:,:,1],Lx,Ly,256)

print(post_processing.compute_s2_from_rdf2d(rdf,r,200,Lx,Ly,1.0))

plt.figure()

plt.plot(r,rdf)

plt.show()