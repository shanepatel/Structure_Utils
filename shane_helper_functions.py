## Secondary stuff so the main atoms routine doesn't become too cluttered ##
import numpy as np
import matplotlib.pyplot as plt
import vasp
from scipy import optimize
import itertools as it
from atoms import *
g=list(it.product([-1,0,1], repeat=3))
def pt_vec_distance(x0,x1,cellvec,system): #x0 is a single point, x1 is a vector of points
    cellpt = [-1,0,1]
    x1 = np.array(x1)
    x1_hold=[]
    for i in g:
        img = np.array(i)
        x_img = np.abs((x1 + img)-x0)
        delta_hold = []
        for j in xrange(len(x1)):
            ## function for finding distance hold up delta_hold.append(x_img[j,0] * dim[0] + x_img[j,1] * dim[1] + x_img[j,2] * dim[2])
        dist_hold = np.linalg.norm(delta_hold,axis=1)
        x1_hold.append(dist_hold)
    x1_hold = np.array(x1_hold)
    min_arr = np.min(x1_hold,axis = 0)
    arg_arr = np.argmin(x1_hold,axis = 0)
    return min_arr, arg_arr
def get_delta(dx,cellvec,system):
    if system == 'Direct':

