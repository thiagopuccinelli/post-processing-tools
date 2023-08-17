import numpy as np 
import tools 

class post_processing:

    def __init__(self):
        self.read_trajectory()
        self.rdf_compute()

    def read_trajectory(filename,Npart,snaps):
       trajectory = tools.read_in_trajectory(filename,Npart,snaps)
       return trajectory 

    def rdf_compute(rx,ry,Lx,Ly,nbins):
        return tools.rdf_monomeric2d(rx,ry,Lx,Ly,nbins)
