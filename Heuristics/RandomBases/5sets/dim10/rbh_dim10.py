
# importing headers

import numpy as np
from matplotlib import pyplot as plt
# from sklearn.decomposition import TruncatedSVD
# import pandas as pd
from scipy.linalg import dft
from scipy.linalg import expm
import cmath
import sympy
from sympy import Matrix, Symbol
from sympy.physics.quantum import TensorProduct
import math
from sympy import symbols
from sympy import expand, factor
from sympy import *
import time
import pandas as pd
import multiprocessing
from multiprocessing import Pool
# import numba
# from numba import jit



#  performs dimension reduction 
def dimension_reduction_svd(Min):
  U,S,V = np.linalg.svd(Min)
  n_component = Min.shape[1] - 1
  VT = np.transpose(V)

  Unew = U
  Sx = np.identity(Min.shape[1])
  for i in range(Min.shape[1]):
    Sx[i][i] = S[i]
  Snew = np.zeros(Min.shape)
  Snew[:Min.shape[0], :Min.shape[1]] = np.diag(Sx)
  Snew = Snew[:, :n_component]

  VTnew = VT[:n_component, :n_component]
  
  Mout = Unew.dot(Snew.dot(VTnew))

  return Mout

# orthogonalise a matrix 
def orthogonalise(Min):
  Q,R = np.linalg.qr(Min)
  return Q

# calculate the distance between 2 bases
def distance(A,B):
  n = A.shape[0]
  Sum = 0
  for i in range(n):
    V1 = np.transpose(np.array(A[i]))
    for j in range(n):
      V2 = np.transpose(np.conj(np.array(B[j])))
      dP = V1 * V2
      dP = pow(abs(sum(dP)),2)
      #print(str(i) + '|' + str(j), dP)
      Sum += (dP * (1-dP))


  Distance = Sum / (n-1)    


  return Distance

# calculates the average squared distance between a set of bases
def ASD(Bases):
  Bases = np.array(Bases)
  #number of bases
  n = Bases.shape[0]
  sum = 0
  for i in range(n):
    for j in range(i+1,n):
      distance_pair = distance(Bases[i],Bases[j])
      sum += float(distance_pair)

  sum = sum * 2
  sum = sum / (n*(n-1))

  return sum    

# performs tensorproduct
def tensorproduct(A,B):
  A = np.array(A)
  B = np.array(B)
  C = TensorProduct(A,B)
  return C

# checks if a matrix is a orthonormal basis
def checkBasis(A):
  A = np.array(A)
  n = A.shape[0] #dimension
  # check 1 : whether the abs value of vectors = 1
  print('CHECK 1 :')
  for i in range(n):
    sum = abs(A[i].dot(np.conj(A[i])))
    #for j in range(n):
     # sum += pow(ab(A[i][j]),2)

    print('abs |' + str(i) + '> = ' + str(sum) )



  # check 2 : the vectors are mutually orthogonal
  print('CHECK 2 :')
  for i in range(n):
    for j in range(n):
      dotProd = abs(A[i].dot(np.conj(A[j])))
      print('<A' + str(i) + '|A' + str(j) + '> = ' + str(dotProd) )


  return



def makeBasis(M,s):
  for i in range(len(M)):
    for j in range(len(M[0])):
      M[i][j] *= s

  return np.array(M)      



def getEigVec(M):
  w,v = np.linalg.eig(M)
  return v




def  calG(A,Bases):
    # calculate the grads
    
    nB = len(Bases)
    nr = len(A)
    G = np.zeros((nr,nr),dtype=complex)
    for i in range(nB): #choose base
        B = Bases[i]
        for a in range(nr): # choose 1st vector
            v1 = A[a]
            v1c = np.conj(v1)
            for b in range(nr): # choose 2nd vector
                v2 = B[b]
                v2c = np.conj(v2)

                ##################
                cross1 = np.outer(v1,v1c)
                cross2 = np.outer(v2,v2c)
                # print(f"v1 : {v1} \nv2 : {v2}")
                # print(f"Cross Prod : {cross}")
                # print(f"Dot product : {dot}") 
                res =  np.matmul(cross1,cross2)
                res_sq = np.matmul(res,res)
                G += res_sq


    G = np.imag(G)
    Gret = (8/((nB * (nB-1)) * (nr-1))) * G
    # print(Gret)

    # print("Checking orthonormality")
    # checkBasis(Gret)

    return Gret



def calV(Grads, step):

    Sigmas = []
    for i in range(len(Grads)):
        Sigmas.append(step * Grads[i])


    Vs = []
    for i in range(len(Sigmas)):
        Vs.append((1j) * Sigmas[i])

    

    return Vs
    # calculate the v values

def unitChange(B,V):

    nr = len(B)
    Bnew = []
    for i in range(nr):
        vb = B[i]
        vbnew = np.matmul(V,vb)
        Bnew.append(np.asarray(vbnew))

    return Bnew


def descent(Bases,step):
    nB = len(Bases)
    nr = len(Bases[0])


    Gradients = []
    for i in range(nB):
        Gradients.append(calG(Bases[i],Bases))
    
    Vmatrices = calV(Gradients,step)

    gradBases = []
    for i in range(nB):
        gradBases.append(unitChange(Bases[i],Vmatrices[i]))

    newBases = []
    for i in range(nB):
        newBases.append(orthogonalise(gradBases[i] + Bases[i]))
    

    return newBases

import random


def measure2(Bases):
    nB = len(Bases)
    d = 1/math.sqrt(10)
    s = []
    for i in range(nB):
        A = Bases[i]
        for j in range(i+1,nB):
            B = Bases[j]
            for x in range(len(A)):
                for y in range(len(B)):
                    v1 = np.array(A[x])
                    v2 = np.array(B[y])
                    p = v1 * v2
                    q = np.abs(np.sum(p))
                    s.append(np.abs(q - d))

    # print(s)
    # print(f"Max : {max(s)}")
    maxs = -1
    for i in range(len(s)):
        if s[i] > maxs:
            maxs = s[i]
    return maxs

if __name__ == "__main__":


    TARGET_MERGED = 0.928
    TARGET_NONMERGED = 0.916

    ITER_Merged = []
    ITER_NonMerged = []

    MaxAsd = []
    

    for run in range(50):
        # print(f"----RUN No. {run}----")
        asdvals = []
    
        R1 = np.random.rand(10,10) - (1j * np.random.rand(10,10))
        R2 = np.random.rand(10,10) + (1j * np.random.rand(10,10))
        R3 = np.random.rand(10,10) - (1j * np.random.rand(10,10))
        R4 = np.random.rand(10,10) + (1j * np.random.rand(10,10))
        R5 = np.random.rand(10,10) + (1j * np.random.rand(10,10))

        R1o = orthogonalise(R1)
        R2o = orthogonalise(R2)
        R3o = orthogonalise(R3)
        R4o = orthogonalise(R4)
        R5o = orthogonalise(R5)

        StartBases = [R1o,R2o,R3o,R4o,R5o]

        sb = ASD(StartBases)
        NewBases = descent(StartBases,0.5)
        NB = NewBases
        max = ASD(NewBases)
        asdvals.append(max)
        max_meas2 = measure2(NewBases)
        itno = -1
        mergepass = -1
        nonmergepass = -1
        # asdvals = []
        # xax = []
        for i in range(10000):
            NB = descent(NB,0.5)
            asd = (ASD(NB))
            itno = i
            if mergepass == -1:
              if asd > TARGET_MERGED:
                ITER_Merged.append(i+1)
                mergepass = 1

            if nonmergepass == -1:
              if asd > TARGET_NONMERGED:
                ITER_NonMerged.append(i+1)
                nonmergepass = 1
            # print(f"iterno -> {i} : ASD -> {asd}")
            # xax.append(i)
            asdvals.append(asd)
            # print(asd)
            if asd > max:
                max = asd
                max_meas2 = measure2(NB)
        
        MaxAsd.append(asdvals)
        print(f"iter{run}. Max ASD : {max} --- Drift : {max_meas2}")
        # MaxAsd.append((max,itno,sb))

    # fo = open("randomStarts.txt","w")
    # fo.write(str(MaxAsd))
    # fo.close()

    print(f"Average Iterations to Cross our Starting Bases -- MERGED : {sum(ITER_Merged)/len(ITER_Merged)} --- NonMerged : {sum(ITER_NonMerged)/len(ITER_NonMerged)}")


    df = pd.DataFrame(MaxAsd)
    df.to_csv('rbh_heur_dim10.csv')




        


