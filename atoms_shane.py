#basic functionality adapted from Max Amsler. I've added additional routines as I've needed them including:
##Swapping Coordinate Systems##
##Moving All Atoms to the Unit Cell##
##Removing Atoms of the Same Species within a tolerance, ie neighboring atoms##
## The most important idea is that atoms are stored in a dictionary, with keys as each element and values as the coordinate in either Direct or Cartesian coordinates
from numpy import linalg
from helper_func import *
import itertools as it
import numpy as np 
import copy
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
    def get_sec_1_cond(self,axis=2,val=0.5,cond='g'): #routine for returning all atoms above certain fraction. Axis is 0,1,2. Val is the cutoff value along that axis, cond is 'l' for less than, 'g' is for greater than val
        #if self.coordinates == 'Cartesian':
        #    self.swap_coordinate_system()
        at_dict={}
        for i in self.at.keys():
            new_coord = []
            for j in self.at[i]:
                if cond == 'l':
                    if j[axis] < val:
                        new_coord.append(j)
                elif cond == 'g':
                    if j[axis] > val:
                        new_coord.append(j)
            at_dict[i] = new_coord
        return(at_dict)
    def get_at_plane(self,axis=2,val=0.5,tol=0.01): #routine for returning all atoms above certain fraction. Axis is 0,1,2. Val is the cutoff value along that axis, cond is 'l' for less than, 'g' is for greater than val
        #if self.coordinates == 'Cartesian':
        #    self.swap_coordinate_system()
        at_dict={}
        for i in self.at.keys():
            new_coord = []
            for j in self.at[i]:
                if abs(j[axis] - val) < tol:
                    new_coord.append(j)
            at_dict[i] = new_coord
        return(at_dict)
  #  def get_sec_1_cond(self,axis=2,val=0.5,cond='g'): #routine for returning all atoms above certain fraction. Axis is 0,1,2. Val is the cutoff value along that axis, cond is 'l' for less than, 'g' is for greater than val
  #      #if self.coordinates == 'Cartesian':
  #      #    self.swap_coordinate_system()
  #      for i in self.at.keys():
  #          new_coord = []
  #          for j in self.at[i]:
  #              if cond == 'l':
  #                  if j[axis] < val:
  #                      new_coord.append(j)
  #              elif cond == 'g':
  #                  print j[axis]
  #                  if j[axis] > val:
  #                      new_coord.append(j)
  #          self.at[i] = new_coord
  #      return(self.at)
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
            self.swap_coordinate_system()
            #print "Please use a direct coordinate system"
            #return []
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
    def add_vac(self,vac_thickness,axis = 2, opt = 2):#default is to add to c, pad symetrically. Vac thickness in angstroms
        self.return_ats_unit()
        swap = True
        if self.coordinates == 'Cartesian':
            self.swap_coordinate_system()
            swap = False
        vac_atoms_vector = np.zeros(3)
        vac_atoms_vector[axis] = 1
        vac_atoms_vector = vac_atoms_vector * vac_thickness
        low = 1
        for i in self.at.values():
            for j in i:
                if j[2] > 0.98:
                    print j[2]
                    if j[2] < low:
                        low = j[2]
        trans_val = 1-low+0.0001
        trans_vec = np.zeros(3)
        trans_vec[axis] = 1
        trans_vec = trans_vec * trans_val
        print self.at
        self.translate_pos(trans_vec)
        self.return_ats_unit()
        print self.at
        self.swap_coordinate_system() # swap to Cartesian
        if opt==2: #symmetric vacuum
            self.cellvec[axis][axis] = self.cellvec[axis][axis] + 2 * vac_thickness
        for i in self.at.keys():
            for j in xrange(len(self.at[i])):
                    self.at[i][j] = self.at[i][j] + vac_atoms_vector
        if swap:
            self.swap_coordinate_system()
    def sort_atoms(self,axis):
        print self.at
        sat = sorted(self.at.keys())
        for i in sat:
            #print sorted(self.at[sat][i], key = lambda x: x[axis])
            self.at[i]= sorted(self.at[i], key = lambda x: x[axis])
    def get_max_min(self,axis=2): #get max and minimum coordinate along a certain axis
        swapped = False
        if self.coordinates == 'Cartesian':
            self.swap_coordinate_system()
            swapped = True
        max_val = 0
        min_val = 1
        for i in self.at.values():
            for j in i:
                if j[axis] < min_val:
                    min_val = j[axis]
                if j[axis] > max_val:
                    max_val = j[axis]
        if swapped:
            self.swap_coordinate_system()
            max_val = max_val * self.cellvec[axis][axis]
            min_val = min_val * self.cellvec[axis][axis]
        return max_val, min_val
    def default_vac(self,axis=2,vac_thicc=2): #finds maximum distance between atoms along an axis, and then creates a vacuum that's 
        swapped = False
        if self.coordinates == 'Direct':
            self.swap_coordinate_system()
            swapped = True
        max_val, min_val = self.get_max_min(axis=axis)
        mid_val = (max_val+min_val)/2.
        cell_height = (max_val - min_val) * vac_thicc
        target_center = 0.5 * cell_height
        print mid_val
        trans_vector = np.zeros(3)
        trans_vector[axis] = 1.
        trans_vector = trans_vector * -(mid_val-target_center)/self.cellvec[axis][axis]
        self.translate_pos(trans_vector)
        print trans_vector
        print trans_vector * self.cellvec[axis][axis]
        self.cellvec[axis][axis] = cell_height
        max_val, min_val = self.get_max_min(axis=axis)
        mid_val = (max_val+min_val)/2.
        print mid_val
        print target_center
        if swapped:
            self.swap_coordinate_system()
    def sub_atom(self,in_type,index,out_type):
        #atoms
        pos = self.at[in_type].pop(index)
        if out_type in self.at.keys():
            self.at[out_type].append(pos)
        else:
            self.at[out_type]=[pos]
    def get_surface_area(self, axis=2):
        vecs = [0,1,2]
        vecs.remove(axis)
        cross = np.cross(self.cellvec[vecs[0]],self.cellvec[vecs[1]])
        return(np.linalg.norm(cross))
    def add_noise_atom(self,amp,spec,num,neg=False):
        if self.coordinates=='Direct':
            self.swap_coordinate_system()
        pos=self.at[spec][num]
        if neg==True:
            g=it.product(range(-1,2,1),repeat=3)
        else:
            g=it.combinations_with_replacement(range(1,6,1),3)
        new_structs=[]
        iters=[]
        for i in g:
            new_ats = copy.deepcopy(self)
            delta = np.array(i)*amp
            new_pos = list(np.array(pos)+delta)
            #print(new_pos)
            new_ats.at[spec][num]=new_pos
            new_structs.append(new_ats)
            iters.append(i)
            #print str(i)
            #print i
        return(new_structs,iters)
