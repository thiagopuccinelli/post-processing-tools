from re import A
import freud 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

temps = [0.1,0.18]#np.arange(0.01, 0.21, 0.01)
temps = np.round(temps, decimals=2)
pressures = np.arange(0.5,15.1,0.5)
pressures = np.sort(np.around(pressures, decimals=1))

for temp in temps:
    for press in pressures:
        print("Processing Temperature: {:.2f}, Pressure: {:.1f}".format(temp,press))
        nsteps = 45
        nmon = 2000
        nmol = 4000
        ndim = 3 
        nbins = 300
        rmax = 25.0
    
        filename = "run/lammps_temp_"+str(temp)+"_press_"+str(press)+".snap"

        N = int(np.genfromtxt(filename, max_rows=1))
        traj = np.genfromtxt(
            filename, skip_header=2,
            invalid_raise=False)[:, 1:4].reshape(-1, N, 3)

        L1 = 1.1*np.absolute(np.max(traj[0,:,0]) - np.min(traj[0,:,0])) 
        L2 = 1.1*np.absolute(np.max(traj[0,:,1]) - np.min(traj[0,:,1])) 
        L3 = 0.0
        box = freud.box.Box(L1,L2,L3)


        type_ids =  np.zeros(nmol) # np.zeros(len(traj))
        for i in range(nmol):
            if i % 2 == 0:
                type_ids[i] = 0.
            else:
                type_ids[i] = 1.

        apos = np.zeros((nsteps,nmon,ndim))
        bpos = np.zeros((nsteps,nmon,ndim))

        for step in range(nsteps):
            j = 0 
            k = 0 
            for pos,i in enumerate(type_ids):
                if i == 0:
                    apos[step][j][0] = traj[step][pos][0]
                    apos[step][j][1] = traj[step][pos][1]
                    apos[step][j][2] = traj[step][pos][2]
                    j += 1 
                if i == 1:
                    bpos[step][k][0] = traj[step][pos][0]
                    bpos[step][k][1] = traj[step][pos][1]
                    bpos[step][k][2] = traj[step][pos][2]
                    k += 1 

        box.center(apos[0][:][:])

        # rdf = freud.density.RDF(bins=nbins,r_max=10.)
        # for step in range(nsteps):
        #     rdf.compute(system=(box,apos[step,:,:]),reset=False)

        newpsi2 = np.zeros(2000)
        newpsi6 = np.zeros(2000)
        for step in range(nsteps):
            ###### hexatic order l=6
            hex_order = freud.order.Hexatic(k=6).compute((box,apos[step,:,:]))
            psi6 = hex_order.particle_order
            mod_psi6 = np.sqrt((psi6[:].real)**2 + (psi6[:].imag)**2)
            for i in range(nmon):
                newpsi6[i] += mod_psi6[i] / nsteps 
            
            ###### l = 2 
            hex_order = freud.order.Hexatic(k = 2)
            hex_order.compute(system=(box, apos[step,:,:]))
            psi2 = hex_order.particle_order
            mod_psi2 = np.sqrt((psi2[:].real)**2 + (psi2[:].imag)**2)
            system_psi2 = mod_psi2.sum()/len(mod_psi2)
            for i in range(nmon):
                newpsi2[i] += mod_psi2[i] / nsteps 
        
        data = np.array([newpsi2,newpsi6])
        
        np.savetxt("results/newpsi62_T_{}_p_{}.dat".format(temp,press),data.T)

