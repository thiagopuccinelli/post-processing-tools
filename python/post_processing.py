import numpy as np 
import tools 

class post_processing:

    def __init__(self):
        self.read_trajectory()
        self.rdf_compute2d()
        self.rdf_compute3d()
        self.compute_s2_from_rdf3d()
        self.compute_s2_from_rdf2d()

    def read_trajectory(filename,Npart,snaps):
       trajectory = tools.read_in_trajectory(filename,Npart,snaps)
       return trajectory 

    def rdf_compute3d(rx,ry,rz,Lx,Ly,Lz,nbins,rmax):
        return tools.rdf_3d(rx,ry,rz,Lx,Ly,Lz,nbins,rmax)
    
    def rdf_compute2d(rx,ry,Lx,Ly,nbins,rmax):
        return tools.rdf_monomeric2d(rx,ry,Lx,Ly,nbins,rmax)

    def compute_s2_from_rdf3d(rdf,r,Npart,Lx,Ly,Lz):
        sex = 0.0 
        for i,val in enumerate(rdf):
            if val > 0:
                sex += (rdf[i] * np.log(rdf[i]) - (rdf[i] - 1))*r[i]**2.
        rho = ( Npart / (Lx * Ly * Lz ))
        sex = 2. * np.pi * rho * sex 
        return sex  
    
    def compute_s2_from_rdf2d(rdf,r,Npart,Lx,Ly):
        sex = 0.0 
        for i,val in enumerate(rdf):
            if val > 0:
                sex += (rdf[i] * np.log(rdf[i]) - (rdf[i] - 1))*r[i]
        rho = ( Npart / (Lx * Ly))
        sex =  np.pi * rho * sex 
        return sex  

    def excess_entropy_cumulant_from_rdf(rdf,r,nbins):
        cs2, rcs2 = tools.excess_entropy_cumulant_from_rdf(rdf,r,nbins)
        return cs2[cs2>0],rcs2[rcs2>0]