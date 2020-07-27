import radial_distribution_evaluation 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm

print(radial_distribution_evaluation.__doc__)

avgc1,rc1,avgf1,rf1 = radial_distribution_evaluation.radial_distribution_close_to_particle(5000,"dump.simulation_eps_5_0_phif_0.15",100,21.878,"0.15")
avgc2,rc2,avgf2,rf2 = radial_distribution_evaluation.radial_distribution_close_to_particle(6000,"dump.simulation_eps_5_0_phif_0.2",100,21.878,"0.20")
avgc3,rc3,avgf3,rf3 = radial_distribution_evaluation.radial_distribution_close_to_particle(7000,"dump.simulation_eps_5_0_phif_0.25",100,21.878,"0.25")
avgc4,rc4,avgf4,rf4 = radial_distribution_evaluation.radial_distribution_close_to_particle(8000,"dump.simulation_eps_5_0_0.3",100,21.878,"0.30")

# LATEX:
plt.rc('text', usetex=True)
mark_size = 10
NC        = 9
cls = cm.rainbow(np.linspace(0, 1, NC))

fig,ax = plt.subplots()

ax.plot(rc1,avgc1,ls="-",lw=2,mec=cls[0],color=cls[0],label=r"WCA, $\phi_f = 0.15$")
ax.plot(rc2,avgc2,ls="-",lw=2,mec=cls[1],color=cls[3],label=r"WCA, $\phi_f = 0.20$")
ax.plot(rc3,avgc3,ls="-",lw=2,mec=cls[2],color=cls[6],label=r"WCA, $\phi_f = 0.25$")
ax.plot(rc4,avgc4,ls="-",lw=2,mec=cls[3],color=cls[8],label=r"WCA, $\phi_f = 0.30$")

ax.set_ylabel(r'$g_{BB}(r)$',fontsize=28)
ax.set_xlabel(r'$r$',fontsize=28)
ax.set_xlim(0,8.0)
ax.set_ylim(0,3.0)
plt.yticks(np.arange(0.,3.1, step=0.5),fontsize=22)
plt.xticks(np.arange(0., 8.1, step=2),fontsize=22)
ax.text(5.0, 0.5, r"\epsilon_{AB} = 5.0", fontsize=26)
plt.legend(loc=1,ncol=1,frameon=False,numpoints = 1,fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("fig7c.pdf")
plt.show()

fig,ax = plt.subplots()

ax.plot(rf1,avgf1,ls="-",lw=2,mec=cls[0],color=cls[0],label=r"WCA, $\phi_f = 0.15$")
ax.plot(rf2,avgf2,ls="-",lw=2,mec=cls[1],color=cls[3],label=r"WCA, $\phi_f = 0.20$")
ax.plot(rf3,avgf3,ls="-",lw=2,mec=cls[2],color=cls[6],label=r"WCA, $\phi_f = 0.25$")
ax.plot(rf4,avgf4,ls="-",lw=2,mec=cls[3],color=cls[8],label=r"WCA, $\phi_f = 0.30$")

ax.set_ylabel(r'$g(r)$',fontsize=28)
ax.set_xlabel(r'$r$',fontsize=28)
ax.set_xlim(0,8.0)
ax.set_ylim(0,3.0)
plt.yticks(np.arange(0.,3.1, step=0.5),fontsize=22)
plt.xticks(np.arange(0., 8.1, step=2),fontsize=22)
ax.text(3.5, 0.5, r"\epsilon_{AB} = 5.0", fontsize=26)
plt.legend(loc=1,ncol=1,frameon=False,numpoints = 1,fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("fig7d.pdf")
plt.show()
