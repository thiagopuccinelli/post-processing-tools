from re import A
import freud 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

temps = np.arange(0.01, 0.21, 0.01)
temps = np.round(temps, decimals=2)
pressures = np.arange(0.5,15.1,0.5)
pressures = np.sort(np.around(pressures, decimals=1))

for temp in temps:
    psi2_res = []
    psi4_res = []
    psi6_res = []
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

        g6_avg = [0 for lenght in range(nbins)] 
        g4_avg = [0 for lenght in range(nbins)]
        g2_avg = [0 for lenght in range(nbins)]
        psi6_avg = 0 
        psi4_avg = 0
        psi2_avg = 0

        for step in range(nsteps):
            ###### hexatic order l=6
            hex_order = freud.order.Hexatic(k=6).compute((box,apos[step,:,:]))
            psi6 = hex_order.particle_order
            mod_psi6 = np.sqrt((psi6[:].real)**2 + (psi6[:].imag)**2)
            system_psi6 = mod_psi6.sum()/len(mod_psi6)
            psi6_avg += system_psi6
            cf6 = freud.density.CorrelationFunction(bins=nbins, r_max=rmax)
            cf6.compute(system=(box, apos[step,:,:]), values=psi6,
                    query_points=apos[step,:,:], query_values=psi6)
            g6 = cf6.correlation
            mod_g6 = np.sqrt((g6[:].real)**2 + (g6[:].imag)**2)
            g6_avg += mod_g6[:]

            ###### l = 4 
            hex_order = freud.order.Hexatic(k = 4)
            hex_order.compute(system=(box, apos[step,:,:]))
            psi4 = hex_order.particle_order
            mod_psi4 = np.sqrt((psi4[:].real)**2 + (psi4[:].imag)**2)
            system_psi4 = mod_psi4.sum()/len(mod_psi4)
            psi4_avg += system_psi4
            cf4 = freud.density.CorrelationFunction(bins=nbins, r_max=rmax)
            cf4.compute(system=(box, apos[step,:,:]), values=psi4,query_points=apos[step,:,:], query_values=psi4)
            g4 = cf4.correlation
            mod_g4 = np.sqrt((g4[:].real)**2 + (g4[:].imag)**2)
            g4_avg += mod_g4[:]

            ###### l = 2 
            hex_order = freud.order.Hexatic(k = 2)
            hex_order.compute(system=(box, apos[step,:,:]))
            psi2 = hex_order.particle_order
            mod_psi2 = np.sqrt((psi2[:].real)**2 + (psi2[:].imag)**2)
            system_psi2 = mod_psi2.sum()/len(mod_psi2)
            psi2_avg += system_psi2
            cf2 = freud.density.CorrelationFunction(bins=nbins, r_max=rmax)
            cf2.compute(system=(box, apos[step,:,:]), values=psi2,query_points=apos[step,:,:], query_values=psi2)
            g2 = cf2.correlation
            mod_g2 = np.sqrt((g2[:].real)**2 + (g2[:].imag)**2)
            g2_avg += mod_g2[:]
            
            
            

        g6_avg = g6_avg[:]/nsteps
        g4_avg = g4_avg[:]/nsteps 
        g2_avg = g2_avg[:]/nsteps 
        psi6_avg /= nsteps
        psi4_avg /= nsteps
        psi2_avg /= nsteps
        psi2_res.append(psi2_avg)
        psi4_res.append(psi4_avg)
        psi6_res.append(psi6_avg)

        np.savetxt('results/g642_temp_'+str(temp)+'_press_'+str(press)+'.dat', np.column_stack((cf6.bin_centers, g6_avg, g4_avg, g2_avg)), header='rbin g6 g4 g2')

    np.savetxt('results/psis_temp_'+str(temp)+'.dat', np.column_stack((pressures, psi2_res, psi4_res, psi6_res)), header='pressure psi2 psi4 psi6')


# fig, ax = plt.subplots(figsize=(5, 3), dpi=200)
# plt.plot(rdf.bin_centers, rdf.rdf, label='freud')
# plt.plot(data[:,0],data[:,1],color="red")
# plt.legend(loc="best")
# plt.show()