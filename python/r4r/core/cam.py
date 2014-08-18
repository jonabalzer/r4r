######################################################################################
#
# Copyright (c) 2013, Jonathan Balzer
#
# All rights reserved.
#
# This file is part of the R4R library.
#
# The R4R library is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The R4R library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R4R library. If not, see <http://www.gnu.org/licenses/>.
#
######################################################################################

import numpy as np
import string

class PinholeCamera:
    """
    Pinhole camera model.
    """
    def __init__(self, s=(0,0), f=(1,1), c=(0,0), k=(0,0,0,0,0), alpha=0, F=np.eye(4,dtype='float64')):
        self.s = s
        self.f = f
        self.c = c
        self.k = k
        self.alpha = alpha
        self.F = F
        self.Finv = np.linalg.inv(self.F)

    def __str__(self):
        return (str(self.s)+"\n"+str(self.f)+"\n"+str(self.c)+"\n"+str(self.k)+
                "\n"+str(self.F)+"\n"+str(self.Finv))                
        
    def from_projection_matrix(self,K):
        """
        Sets parameters from a given projection matrix.
        """
        self.f[0] = K[0,0]
        self.f[1] = K[1,1]
        self.c[0] = K[0,2]
        self.c[1] = K[1,2]
        
    def read_from_r4r_file(self,fn):
        """
        Reads parameters from file in R4R format.
        """
        fid = open(fn,'r')
        
        # img size
        fid.readline()            
        line = fid.readline()
        ls = string.split(line)
        self.s = (int(ls[0]),int(ls[1]))       
        
        # focal lengths
        fid.readline()            
        line = fid.readline()
        ls = string.split(line)
        self.f = (float(ls[0]),float(ls[1]))      
        
        # principal point
        fid.readline()            
        line = fid.readline()
        ls = string.split(line)
        self.c = (float(ls[0]),float(ls[1]))        
    
        # distortion parameters
        fid.readline()            
        line = fid.readline()
        ls = string.split(line)
        self.k = (float(ls[0]),float(ls[1]),float(ls[2]),float(ls[3]),float(ls[4]))

        # skew parameter
        fid.readline()            
        line = fid.readline()
        ls = string.split(line)
        self.alpha = float(ls[0])
        
        # F
        fid.readline()            

        for i in range(0,4):
            line = fid.readline()
            ls = string.split(line)
            for j in range(0,4):
                self.F[i,j] = float(ls[j])
        
        self.Finv = np.linalg.inv(self.F)
        
        fid.close()
    
    def write_to_r4r_file(self,fn):
        """
        Exports intrinsic parameters to file in R4R format.
        """
        fid = open(fn,'w')
        fid.write("# size\n")
        fid.write(str(self.s[0])+" "+str(self.s[1])+"\n")
        fid.write("# focal length\n")
        fid.write(str(self.f[0])+" "+str(self.f[1])+"\n")
        fid.write("# principal point\n")
        fid.write(str(self.c[0])+" "+str(self.c[1])+"\n")
        fid.write("# distortion\n")
        fid.write(str(self.k[0])+" "+str(self.k[1])+" "+str(self.k[2])+" "+str(self.k[3])+" "+str(self.k[4])+"\n")
        fid.write("# skew coefficient\n")
        fid.write(str(self.alpha)+"\n")
        fid.write("# frame world -> cam\n")
        
        Fstr = str(self.F)
        Fstr = Fstr.replace('[[ ','')
        Fstr = Fstr.replace(' [ ','')
        Fstr = Fstr.replace(']]','')
        Fstr = Fstr.replace(']','')
        fid.write(Fstr)
        
        fid.close()

    def projection_matrix(self):
        """
        Puts the parameters into a matrix.
        """
        K = np.zeros((3,3),dtype='float64')
        K[0,0] = self.f[0]        
        K[1,1] = self.f[1]        
        K[0,2] = self.c[0]
        K[1,2] = self.c[1]
        K[2,2] = 1.0
        return K
        
    def pixel_grid(self):
        """
        Creates a mesh grid according to image dimensions.
        """
        return np.meshgrid(np.arange(0,self.s[0]),np.arange(0,self.s[1]))
        
    def back_project(self,u,v,z):
        """
        Back-projects a set of pixels and depths to 3-d. CAVEAT: It is assumed
        that there is no radial distortion.
        """
        x = z*(u-self.c[0])/self.f[0]
        y = z*(v-self.c[1])/self.f[1]
                      
        return x,y,z