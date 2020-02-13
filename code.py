import radial_distribution_evaluation 
import numpy as np 
import matplotlib.pyplot as plt 

print(radial_distribution_evaluation.__doc__)

avgc,rc,avgf,rf = radial_distribution_evaluation.radial_distribution_close_to_particle(10000,"dump.simulation_eps_5_0.4",100,21.878)

# LATEX:
plt.rc('text', usetex=True)

fig,ax = plt.subplots()

ax.plot(rc,avgc,ls="-",color="black")
ax.plot(rf,avgf,ls="-",color="red")
ax.set_ylabel(r'$g(r)$',fontsize=28)
ax.set_xlabel(r'$r$',fontsize=28)
#ax.set_xlim(0,3.5)
#ax.set_ylim(0,4.0)
#plt.yticks(np.arange(0.,4.1, step=0.5),fontsize=22)
#plt.xticks(np.arange(0., 3.51, step=0.5),fontsize=22)
ax.grid()
plt.tight_layout()

plt.show()
