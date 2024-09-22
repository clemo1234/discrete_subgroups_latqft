#imports

import numpy as np
import numba
import random
import math
import sys
import matplotlib.pylab as plt
from multiprocessing import Pool
#from pyinstrument import Profiler
from guppy import hpy
from memory_profiler import profile
import gc

#To-Do: Create Class Structure (group class , lattice class, observables), minimize use of global scope variables in local scope sections, Use of pytorch tensor operations, Define params and meaning of functions and classes, Unit tests, HPC versility, Easier modification of action

#profiler = Profiler()
#profiler.start()

# Constants
K = 1000  # Decorrelation time
N = 100      # Number of samples
#PMAX = 648

# Global variables
#P = 648


hit = int()
acc = float()

# Python equivalents of C++ functions


def load_group(filename):
    """
    Add function definition
    """
    #global P, ReTr, ImTr, mult, inv, id, a

    

    with open(filename, 'r') as f:
        #P = int(f.readline().strip())
        
        P = int(f.readline())
        #print("Pval ",P)
        PMAX = P
        ReTr = np.zeros(PMAX)
        ImTr = np.zeros(PMAX)
        mult = np.zeros((PMAX, PMAX), dtype=int)
        inv = np.zeros(PMAX, dtype=int)
        #nn = np.zeros(PMAX, dtype=int)
        #smallgroup = np.zeros(PMAX)
        a = np.zeros(PMAX)
        a = np.full((V * D), 0) 
        if P > PMAX:
            raise ValueError(f"Order of group too large: {P} > {PMAX}")
        
        # Read ReTr values
 
        ReTrtemp = f.readline().split()
        ReTr = np.array(ReTrtemp).astype(float)
        #print(ReTr)
        # Read ImTr values
 
        ImTrtemp = f.readline().split()
        ImTr = np.array(ImTrtemp).astype(float)

        # Read multiplication table
        for n in range(P):
                mult[n,:] = f.readline().split()
        mult = mult.astype(int)
        
        
    

    # Find the identity element
    id_found = False
    for n in range(P):
        for m in range(P):
            if mult[n, m] == m and mult[m, n] == m:
                id = n
                #print(f"id print: {id} \n")
                id_found = True
                break
    if not id_found:
        raise ValueError("Group does not have identity element")

    # Find inverses
    for n in range(P):
        for m in range(P):
            if mult[n, m] == id and mult[m, n] == id:  #mask
                inv[n] = m
                break
        else:
            raise ValueError("Group does not have inverses")
    return P, ReTr, ImTr, mult, inv, id, a

@numba.njit
def step(i, d, s, nn):
    """
    Add function definition
    """
    #print(i,d,s)
    if np.abs(s) != 1:
        raise ValueError("Not implemented ... exiting")
    #if nn.all() == 0:
    #    init_nn()
    return nn[i * 2 * D + d * 2 + (1 + s) // 2]

@numba.njit
def init_nn(nn):
    """
    Add function definition
    """
    #possible grouping into lattice class
    pos = np.zeros(D, dtype=np.uint32)
    #print("Initializing nearest neighbours... ", f"[{Nt}, {Nx}] \n")
    for i in range(V):
        pos,idx = getpos(i, pos)
        #print(pos)
        for d in range(D):
            pos[d] -= 1
            nn[(2*D*i)+(2*d+0)] = getidx(pos)
            pos[d] += 2
            nn[(2*D*i)+((2*d)+1)] = getidx(pos)
            pos[d] -= 1
    return pos, nn, idx
        
        #print(f"{i:4d}: ", end='')
        #for j in range(2 * D):
            #print(f"{nn[2*D*i+j]:4d} ", end='')
        #print()
@numba.njit
def getidx(pos):
    """
    Add function definition
    """
    idx = 0
    for d in range(D - 1, 0, -1):
        idx = ((pos[d] + Nx) % Nx) + (Nx * idx)
    return ((pos[0] + Nt) % Nt) + (Nt * idx)

@numba.njit
def getpos(idx, pos):
    """
    Add function definition
    """
    #pos = [0] * D
    pos[0] = idx % Nt
    idx //= Nt
    for d in range(1, D):
        pos[d] = idx % Nx
        idx //= Nx
    #print(pos)
    return pos,idx

@numba.njit
def plaquette(i, d1, d2, a, mult, inv, id, nn):
    """
    Add function definition
    """
    g = id
    g = mult[g, a[i * D + d1]]
    g = mult[g, a[step(i, d1, 1, nn) * D + d2]]
    g = mult[g, inv[a[step(i, d2, 1, nn) * D + d1]]]
    g = mult[g, inv[a[i * D + d2]]]
    return g

@numba.njit
def polyakov(i, a, id, nn, mult):
    """
    Add function definition
    """
    r = id
    i0 = i
    #print(mult[4])
    while True:
        r = mult[r, a[i * D]]
        #print(f"r val: {a[36]}")
        i = step(i, 0, 1, nn)
        if i0 == i:
            break
    return r

@numba.njit
def getpoly(re,im, pos, a, ReTr, ImTr, id, nn, mult):
    """
    Add function definition
    """
    lre = 0
    lim = 0

    for i in range(V):
        pos, idx = getpos(i,pos)
        if pos[0] != 0:
            continue
        #print(i)
        p = polyakov(i, a, id, nn, mult)
        lre += ReTr[p]
        lim += ImTr[p]
    re = lre / (V / Nt)
    im = lim / (V / Nt)
    return re, im

@numba.njit
def wilson(i, d1, d2, nd1, nd2, a, mult, id, nn, inv):
    """
    Add function definition
    """
    g = id
    for _ in range(nd1):
        g = mult[g, a[i * D + d1]]
        i = step(i, d1, 1, nn)
    for _ in range(nd2):
        g = mult[g, a[i * D + d2]]
        i = step(i, d2, 1, nn)
    for _ in range(nd1):
        i = step(i, d1, -1, nn)
        g = mult[g, inv[a[i * D + d1]]]
    for _ in range(nd2):
        i = step(i, d2, -1, nn)
        g = mult[g, inv[a[i * D + d2]]]
    return g

@numba.njit
def get_staples(i, d, st, a, mult, inv, nn):
    """
    Add function definition
    """
    i1 = step(i, d, 1, nn)
    k = 0
    for d1 in range(D):
        if d1 != d:
            g = a[i1 * D + d1]
            g = mult[g, inv[a[step(i, d1, 1, nn) * D + d]]]
            g = mult[g, inv[a[i * D + d1]]]
            #print(f"this is g: {mult[502][200]}")
            st[k] = g
            k += 1
            g = inv[a[step(i1, d1, -1, nn) * D + d1]]
            i2 = step(i, d1, -1, nn)
            g = mult[g, inv[a[i2 * D + d]]]
            g = mult[g, a[i2 * D + d1]]
            st[k] = g
            k += 1
    return st


@numba.njit
def update(beta1,pos,idx, smallgroup, hit, acc, a, mult, ReTr, ImTr, inv, nn):
    #randsmallgrp = random.randint(0, len(smallgroup) - 1)
    #rand01 = random.random  #random num generator [0,1)
    rand01 = lambda: random.uniform(0., 1.)
    smallgroup_len = len(smallgroup)
    two_D_minus_1 = 2 * (D - 1)
    for d in range(D):
        for parity in range(2):
            lhit = 0
            lacc = 0
            for i in range(0,V):
                pos,idx = getpos(i,pos)
                if (pos[0] + pos[1] + pos[2] + pos[3]) % 2 != parity:
                    continue
                
                #print(staples)
                staples = np.full((2 * (D-1)), 0)
                staples = get_staples(i, d, staples, a, mult, inv, nn)
                #print(staples)
                nhit = 20
                for _ in range(nhit):
                    r1 = 0
                    r2 = 0
                    r1new = 0
                    r2new = 0
                    b = a[i * D + d]
                    #print(f"small: {(smallgroup[int(len(smallgroup) * temp_rand)])}")
                    bnew = mult[b, smallgroup[int(smallgroup_len * rand01())]]
                    for j in range(two_D_minus_1):
                        p1 = mult[b, staples[j]]
                        #print(f"p1 val: {p1}")
                        r1 += ReTr[p1]
                        r2 += ReTr[p1] * ReTr[p1] + ImTr[p1] * ImTr[p1]

                        p1 = mult[bnew, staples[j]]
                        r1new += ReTr[p1]
                        r2new += ReTr[p1] * ReTr[p1] + ImTr[p1] * ImTr[p1]

                    oldact = beta0 + (beta1 * r1 + beta2 * r2)
                    newact = beta0 + (beta1 * r1new) + (beta2 * r2new)
                    probrat = np.exp(-(newact - oldact))
                    lhit += 1
                    #print(f"prob : {temp_rand} and {probrat}")
                    if probrat > rand01():
                        #print(f"bnew {bnew}")
                        a[i * D + d] = bnew
                        lacc += 1

            hit += lhit
            acc += lacc

    return hit, acc, a
            #print(hit,acc)


groupfilename = "myBI.txt"
D = 4
Nt = 2
Nx = 2
beta0 = 0.0
#beta1 = int()
beta2 = 0.0
iseed = random.randint(0, 10000)
#plaqb = [[],[]]
V =Nt*pow(Nx,D-1)
#with open('output_test_python.txt', 'w') as file:


@profile
def simulate(b_vals):
        #global beta1
        
        hit = 0
        acc = 0
        temp =[]
        beta1 = b_vals
        random.seed(2)
        a = np.full((V * D), 0, dtype=int)  # Initialize gauge field array with zeros
        nn = np.full((V * D * 2), 0, dtype=int)
        
        # Load the multiplication table
        P, ReTr, ImTr, mult, inv, id, a = load_group(groupfilename)
        randgrp = lambda: random.randint(0, P)
        #print(ReTr)
        #print(f"check stuff: (a*b)*c: {mult[mult[1][2]][3]} a*(b*c): {mult[1][mult[2][3]]} \n")
        # Initialize small group
        smallgroup = []
        min_retr = min(ReTr)
        nn_retr = min(r for r in ReTr if r > min_retr)
        #print(nn_retr)
        #temp_nn = []
        #for r in ReTr:
        #    if r > min_retr:
        #        temp_nn.append(r)
        #nn_retr = min(temp_nn)
        smallgroup = [i for i, retr in enumerate(ReTr) if min_retr < retr and retr < nn_retr + 1e-6]
        #print(f"Retr: {smallgroup}")
        #print(f"min_retr: {min_retr} nn_retr: {nn_retr}")
        #print(f"small group size: {len(smallgroup)}")

        # Initialize the gauge field
        for m in range(V*D):
            #a[m] = id
            a[m] = randgrp()
        #print(a)
        #step(0, 0, 1)  # Prime the nn table
        pos, nn, idx = init_nn(nn)
        #print("init done")

        # Thermalization
        for _ in range(K*0):
            hit, acc, a = update(beta1, pos,idx,smallgroup, hit, acc, a, mult, ReTr, ImTr, inv,nn)
        #print("thermo done")
    
        for n in range(N):
            for _ in range(K):
                hit, acc, a = update(beta1, pos,idx,smallgroup, hit, acc, a, mult, ReTr, ImTr, inv,nn)
            rep = float()
            imp = float()
            rep, imp = getpoly(rep, imp, pos, a, ReTr, ImTr, id, nn, mult)
            #print(rep,imp)
            simpleplaq = 0
            wloop = np.zeros((Nx, Nt))

            for i in range(V):
                for d1 in range(D):
                    for d2 in range(d1 + 1, D):
                        simpleplaq += ReTr[plaquette(i, d1, d2, a, mult, inv, id, nn)]
            
            for k in range(Nx):
                for l in range(Nt):
                    lwloop = 0
                    for i in range(V):
                        for d1 in range(1, D):
                            lwloop += ReTr[wilson(i, d1, 0, k + 1, l + 1,a, mult, id,nn,inv)]
                    wloop[k, l] = lwloop / (V * D * (D - 1) / 2)

            simpleplaq /= ((V * D * (D - 1)) / 2)
            #print(f"GMES: {999.0} {rep} {imp} {simpleplaq}", end=" ")
            #for k in range(Nx):
            #    print(" ".join(f"{wloop[k, l]:e}" for l in range(Nt)), end=" ")
            #print()

            #print(f"ACC: {acc / hit}")

            #print("CONFIGS: ", end="")
            #print(" ".join(f"{ai:03d}" for ai in a))
            #test = " ".join(f"{wloop[k, l]:e}" for l in range(Nt))
            temp.append(simpleplaq)

            #print(temp)
            #plaqb[0].append(beta1)
            gc.collect()
        return np.mean(temp), beta1
        #file.write(f"Beta 1 value: {beta1} \n")
        #file.write(f"GMES: {999.0} {rep} {imp} {simpleplaq} {test} \n")


#h= hpy()
#print(h.heap())
#print(h.heap()[0].referrers.byvia)
#print(h.heap()[0].rp)
#print(h.heap()[0].sp)

def convert_to_2d_xy_array(coords):
    """
    Convert an array of coordinate tuples into a 2D x, y array.

    Args:
    coords (list of tuples): An array of coordinate tuples, e.g., [(1, 2), (3, 4), (5, 6)]

    Returns:
    np.array: A 2D array where the first row contains x-coordinates and the second row contains y-coordinates
    """
    # Convert the list of tuples to a NumPy array
    coords_array = np.array(coords)
    
    # Transpose the array to switch rows and columns
    xy_array = coords_array.T
    
    return xy_array


def parallel_loop():
    data = np.linspace(2.9/3,4.5/3,20)

    with Pool() as pool:
        results = pool.map(simulate,data)
    return convert_to_2d_xy_array(results)
    

if __name__ == "__main__":
    x,y=parallel_loop()
    #profiler.stop()
    #profiler.print()
    with open('array.txt', 'w') as f:
        f.write("["+','.join(map(str, (1+(x/3)))) + ']' + "\n")
        f.write("["+','.join(map(str, y*3)) + ']')
    plt.scatter(y*3,(1+(x/3)))
    plt.xlabel(r"$\beta_1$")
    plt.ylabel("Plaq. Expectation")
    plt.show()

