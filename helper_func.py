from numpy import linalg
import itertools as it
import numpy as np
def distance_matrix(vec1,vec2,cellvec):
        outmat = []
        indmat = []
        for i in xrange(len(vec1)):
            distvec,indvec = pt_vec_distance(vec1[i],vec2,cellvec)
            outmat.append(distvec)
            indmat.append(indvec)
        return outmat,indmat
def pt_vec_distance(x0,x1,cellvec): #x0 is a single 3-D point, x1 is a vector of points 
    g=list(it.product([-1,0,1], repeat=3))
    x1_hold=[] #will contain distance between x0 and all points in x1 for all images
    for i in g:
        img = np.array(i)
        x_img = np.abs((x1 + img) -x0) # this is an array with the distances between all of the points in the image and x0. Needs to be converted to a cartesian distance
        direct_img_dist = np.dot(x_img,np.transpose(cellvec)) 
        dist_hold = np.linalg.norm(direct_img_dist,axis=1)
        x1_hold.append(dist_hold)
    x1_hold = np.array(x1_hold)
    min_arr = np.min(x1_hold,axis = 0)
    arg_arr = np.argmin(x1_hold,axis = 0)
    return min_arr, arg_arr
