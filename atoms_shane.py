#basic functionality adapted from Max Amsler. I've added additional routines as I've needed them including:
##Swapping Coordinate Systems##
##Moving All Atoms to the Unit Cell##
##Removing Atoms of the Same Species within a tolerance, ie neighboring atoms##
## The most important idea is that atoms are stored in a dictionary, with keys as each element and values as the coordinate in either Direct or Cartesian coordinates
from numpy import linalg
from helper_func import *
import itertools as it
import numpy as np 
class Atoms:
    def __init__(self):
        self.at={}
        self.bemoved=[]
        self.qat=[]
        self.nat=-1
        self.nmolecule=-1
        self.cellvec=[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
        self.boundcond="unknown"
        self.epot=1.E100
        self.qtot=0.0
        self.coordinates="Direct"
        self.pattern=-1
    def swap_coordinate_system(self):
        if self.coordinates=='Direct': #convert to Cartesian
            for i in self.at.keys():
                new_coord = []
                for j in self.at[i]:
                    xPos = j[0]*self.cellvec[0][0] + j[1]*self.cellvec[1][0] + j[2]*self.cellvec[2][0]
                    yPos = j[0]*self.cellvec[0][1] + j[1]*self.cellvec[1][1] + j[2]*self.cellvec[2][1]
                    zPos = j[0]*self.cellvec[0][2] + j[1]*self.cellvec[1][2] + j[2]*self.cellvec[2][2]
                    new_coord.append([xPos, yPos, zPos])
                self.at[i] = new_coord
            self.coordinates='Cartesian'
        elif self.coordinates == 'Cartesian':
            latCnt = [x[:] for x in [[None]*3]*3]
            for a in range(3):
                for b in range(3):
                    latCnt[a][b] = self.cellvec[b][a]
            detLatCnt = linalg.det(np.array(latCnt))
            for i in self.at.keys():
                new_coord = []
                for j in self.at[i]:
                    aPos = (linalg.det(np.array([[j[0], latCnt[0][1], latCnt[0][2]], [j[1], latCnt[1][1], latCnt[1][2]], [j[2], latCnt[2][1], latCnt[2][2]]]))) / detLatCnt
                    bPos = (linalg.det(np.array([[latCnt[0][0], j[0], latCnt[0][2]], [latCnt[1][0], j[1], latCnt[1][2]], [latCnt[2][0], j[2], latCnt[2][2]]]))) / detLatCnt
                    cPos = (linalg.det(np.array([[latCnt[0][0], latCnt[0][1], j[0]], [latCnt[1][0], latCnt[1][1], j[1]], [latCnt[2][0], latCnt[2][1], j[2]]]))) / detLatCnt
                    new_coord.append([aPos, bPos, cPos])
                self.at[i] = new_coord
            self.coordinates='Direct'
    def translate_pos(self,vector): #routine for translating all atoms. Vector should be in fractional coordinates
        if self.coordinates =='Direct':
            for i in self.at.keys():
                new_coord = []
                for j in self.at[i]:
                    new_coord.append((np.array(j)+np.array(vector)).tolist())
                self.at[i] = new_coord
        elif self.coordinates=='Cartesian':
            self.swap_coordinate_system()
            for i in self.at.keys():
                new_coord = []
                for j in self.at[i]:
                    new_coord.append((np.array(j)+np.array(vector)).tolist())
                self.at[i] = new_coord
            self.swap_coordinate_system()
    def translate_pos_slab(self,vector,axis=2,val=0.5,cond='g'): #routine for translating all atoms. Vector should be in fractional coordinates. Axis is 0,1,2. Val is the cutoff value along that axis, cond is 'l' for less than, 'g' is for greater than val
        init_coord = self.coordinates
        if self.coordinates == 'Cartesian':
            self.swap_coordinate_system()
        for i in self.at.keys():
            new_coord = []
            for j in self.at[i]:
                if cond == 'l':
                    if j[axis] < val:
                        new_coord.append((np.array(j)+np.array(vector)).tolist())
                    else:
                        new_coord.append(j)
                elif cond == 'g':
                    if j[axis] > val:
                        new_coord.append((np.array(j)+np.array(vector)).tolist())
                    else:
                        new_coord.append(j)
            self.at[i] = new_coord
        if init_coord=='Cartesian':
            self.swap_coordinate_system()
    def delete_sec_1_cond(self,axis=2,val=0.5,cond='g'): #routine for deleting all atoms above certain fraction. Axis is 0,1,2. Val is the cutoff value along that axis, cond is 'l' for less than, 'g' is for greater than val
        init_coord = self.coordinates
        if self.coordinates == 'Cartesian':
            self.swap_coordinate_system()
        for i in self.at.keys():
            new_coord = []
            for j in self.at[i]:
                if cond == 'l':
                    if j[axis] > val:
                        new_coord.append(j)
                elif cond == 'g':
                    if j[axis] < val:
                        new_coord.append(j)
            self.at[i] = new_coord
        if init_coord=='Cartesian':
            self.swap_coordinate_system()
    def return_ats_unit(self):
        if self.coordinates=='Direct':
            for i in self.at.keys():
                new_coord = []
                for j in self.at[i]:
                    new_coord.append([k % 1 for k in j])
                self.at[i]=new_coord
        elif self.coordinates=='Cartesian':
            self.swap_coordinate_system()
            for i in self.at.keys():
                new_coord = []
                for j in self.at[i]:
                    new_coord.append([k % 1 for k in j])
                self.at[i]=new_coord
            self.swap_coordinate_system()
    def remove_overlap(self, tol):
        if self.coordinates =='Cartesian':
            print "Please use a direct coordinate system"
            return []
        dist_matrix_array=[]#array of arrays containing distances between point and the closest image
        ind_matrix_array=[]#array of arrays containing distances between point and the closest image
        sat = sorted(self.at.keys())
        for i in sat:
            dist_mat,ind_mat = distance_matrix(self.at[i],self.at[i],self.cellvec)
            dist_matrix_array.append(dist_mat)
            ind_matrix_array.append(ind_mat)
        for i in xrange(len(sat)):
            mat = np.array(dist_matrix_array[i])
            diag = ~np.eye(mat.shape[0],dtype=bool)
            boolmat = np.logical_and(mat < tol,diag)
            arr_1, arr_2 = np.where(boolmat)
            new_arr = []
            for j in reversed(xrange(len(arr_1))):
                if arr_1[j] > arr_2[j]:
                    del self.at[sat[i]][arr_1[j]]
    def swap_axis(self,ax1,ax2):
        sat = sorted(self.at.keys())
        for i in sat:
            for item in self.at[i]:
                item[ax1],item[ax2] =item[ax2],item[ax1] 
        for j in self.cellvec:
            j[ax1],j[ax2] = j[ax2],j[ax1]
        self.cellvec[ax1],self.cellvec[ax2] = self.cellvec[ax2],self.cellvec[ax1]
