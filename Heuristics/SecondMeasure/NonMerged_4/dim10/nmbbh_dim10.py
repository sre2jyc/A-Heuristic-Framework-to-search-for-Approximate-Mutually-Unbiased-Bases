
# importing headers

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
# from sklearn.decomposition import TruncatedSVD
from scipy.linalg import dft
import cmath
import sympy
from sympy import Matrix, Symbol
from sympy.physics.quantum import TensorProduct
import math
from sympy import symbols
from sympy import expand, factor
from sympy import *
import time

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

# function to merge multiple bases into a large matrix on which dimension reduction is performed
def mergeBases(B,n,d):
  M = []
  nb = n
  for i in range(nb):
    for j in range(d):
      M.append(B[i][j])
      #np.append(M,B[i][j],axis = 0)

  #M = np.array(M)
  return M

# all possible bases
def allPossibleBases(M,d):
  B = []
  n = M.shape[0]
  for i in range(n):
    Bi = []
    for j in range(i,i+d):
      Bi.append(M[j%n])

    Bi = orthogonalise(Bi)  
    B.append(Bi)
  B = np.array(B)  
  return B

def makeBasis(M):
  for i in range(len(M)):
    for j in range(len(M[0])):
      M[i][j] *= s

  return np.array(M)      


def searchFourBoxes(B):
  nG = B.shape[0]  # number of base groups, here 8
  nB = B.shape[1]   # number of bases in each group
  #ASDs = []
  #points = []
  Results = []
  for i in range(nG):
    for j in range(i+1,nG):
      for k in range(j+1,nG):
        for l in range(k+1,nG):
          BG1 = B[i]
          BG2 = B[j]
          BG3 = B[k]
          BG4 = B[l]

          for a in range(nB):
            for b in range(nB):
              for c in range(nB):
                for d in range(nB):
                  B1 = BG1[a]
                  B2 = BG2[b]
                  B3 = BG3[c]
                  B4 = BG4[d]

                  asd = ASD([B1,B2,B3,B4])
                  # print('ASD is : ' + str(asd))
                  #ASDs.append(asd)
                  point = [(a),(b),(c),(d)]
                  Results.append((asd,point))

  return Results                  


def chooseFourBoxes(B):
  AllRes = []
  nB = 4
  for i in range(nB):
    for j in range(i+1,nB):
      for k in range(j+1,nB):
        for l in range(k+1,nB):
          choosedBases = np.array([B[i],B[j],B[k],B[l]])
          Results = searchFourBoxes(choosedBases)
          SortedRes = Results
          SortedRes.sort(reverse = True)
          AllRes.append([B[4],SortedRes[:10]])

  return AllRes          


def getEigVec(M):
  w,v = np.linalg.eig(M)
  return v

def creatingStartingBases():
    pi = math.pi
    o = math.cos(2*pi/11) + (math.sin(2*pi/11)*(1j))
    Z11 = [[1,0,0,0,0,0,0,0,0,0,0],[0,o**1,0,0,0,0,0,0,0,0,0],[0,0,o**2,0,0,0,0,0,0,0,0],[0,0,0,o**3,0,0,0,0,0,0,0],[0,0,0,0,o**4,0,0,0,0,0,0],[0,0,0,0,0,o**5,0,0,0,0,0],[0,0,0,0,0,0,o**6,0,0,0,0],[0,0,0,0,0,0,0,o**7,0,0,0],[0,0,0,0,0,0,0,0,o**8,0,0],[0,0,0,0,0,0,0,0,0,o**9,0],[0,0,0,0,0,0,0,0,0,0,o**10]] 
    Z11 = np.array(Z11)
    X11 = [[0,0,0,0,0,0j,0,0,0,0j,1+0j],[1,0,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,1,0]]
    X11 = np.array(X11)
    M1 = Z11
    M2 = X11
    M3 = np.dot(M2,M1)
    M4 = np.dot(M3,M1)
    M5 = np.dot(M4,M1)
    M6 = np.dot(M5,M1)
    M7 = np.dot(M6,M1)
    M8 = np.dot(M7,M1)
    M9 = np.dot(M8,M1)
    M10 = np.dot(M9,M1)
    M11 = np.dot(M10,M1)
    M12 = np.dot(M11,M1)
    B01 = np.transpose(getEigVec(M1))
    B02 = np.transpose(getEigVec(M2))
    B03 = np.transpose(getEigVec(M3))
    B04 = np.transpose(getEigVec(M4))
    B05 = np.transpose(getEigVec(M5))
    B06 = np.transpose(getEigVec(M6))
    B07 = np.transpose(getEigVec(M7))
    B08 = np.transpose(getEigVec(M8))
    B09 = np.transpose(getEigVec(M9))
    B10 = np.transpose(getEigVec(M10))
    B11 = np.transpose(getEigVec(M11))
    B12 = np.transpose(getEigVec(M12))

    Bases = [B01,B02,B03,B04,B05,B06,B07,B08,B09,B10,B11,B12]

    return Bases


def constructingBoxes(StartBases):

    # MergedBases = mergeBases(StartBases,12,11)
    # MergedBases = np.array(MergedBases)

    # ReducedBases = dimension_reduction_svd(MergedBases)
    # ReducedBases = np.array(ReducedBases)

    RB1 = dimension_reduction_svd(StartBases[0])
    RB2 = dimension_reduction_svd(StartBases[1]) 
    RB3 = dimension_reduction_svd(StartBases[2])
    RB4 = dimension_reduction_svd(StartBases[3]) 
    RB5 = dimension_reduction_svd(StartBases[4])
    RB6 = dimension_reduction_svd(StartBases[5]) 
    RB7 = dimension_reduction_svd(StartBases[6])
    RB8 = dimension_reduction_svd(StartBases[7]) 
    RB9 = dimension_reduction_svd(StartBases[8])
    RB10 = dimension_reduction_svd(StartBases[9]) 
    RB11 = dimension_reduction_svd(StartBases[10])
    RB12 = dimension_reduction_svd(StartBases[11])

    APB1 = allPossibleBases(RB1,10)
    APB2 = allPossibleBases(RB2,10)
    APB3 = allPossibleBases(RB3,10)
    APB4 = allPossibleBases(RB4,10)
    APB5 = allPossibleBases(RB5,10)
    APB6 = allPossibleBases(RB6,10)
    APB7 = allPossibleBases(RB7,10)
    APB8 = allPossibleBases(RB8,10)
    APB9 = allPossibleBases(RB9,10)
    APB10 = allPossibleBases(RB10,10)
    APB11 = allPossibleBases(RB11,10)
    APB12 = allPossibleBases(RB12,10)


    Boxes = [APB1,APB2,APB3,APB4,APB5,APB6,APB7,APB8,APB9,APB10,APB11,APB12]
    return Boxes

def insideProcess(BS): 
  nBoxSets = len(BS)
  results = []
  for i in range(nBoxSets):
    results.append(chooseFourBoxes(BS[i]))
    print(f"Result appended {i}")
  return results

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

def calG(A,Bases):
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



if __name__ == "__main__":
    # importFiles()
    StartBases = creatingStartingBases()
    # print(ASD(StartBases))

    Boxes = constructingBoxes(StartBases)
    Boxes = np.array(Boxes)


    print(Boxes.shape)

    Boxsets = []

    for i in range(Boxes.shape[0]):
        for j in range(i+1,Boxes.shape[0]):
            for k in range(j+1,Boxes.shape[0]):
                for l in range(k+1,Boxes.shape[0]):
                    Boxsets.append(np.array([Boxes[i],Boxes[j],Boxes[k],Boxes[l],[i,j,k,l]],dtype=object))
    



    bb1 = Boxes[3][7]
    bb2 = Boxes[4][5]
    bb3 = Boxes[5][5]
    bb4 = Boxes[6][6]

    BestBases = [bb1,bb2,bb3,bb4]
    # checkBasis(bb1)
    # checkBasis(bb2)
    # checkBasis(bb3)
    # checkBasis(bb4)






    print("NONMERGED 4SETS : DIM10")

    print(f"Average Sq Dist of Best Set of 4 Bases : {ASD(BestBases)}")


    meas2 = measure2(BestBases)
    # print((meas2))
    print(f"Max Value of Measure2 : {(meas2)}")

    NewBases = descent(BestBases,0.5)
    NB = NewBases
    max = ASD(NewBases)
    max_meas2 = measure2(NewBases)
    # max_meas2 = max(max_meas2)
    asdvals = []
    m2vals = []
    asdvals.append(max)
    m2vals.append(max_meas2)
    # xax = []
    for i in range(10000):
      NB = descent(NB,0.5)
      asd = (ASD(NB))
      m2 = measure2(NB)
      print(f"Iter no : {i} --> ASD : {asd}  Drift : {m2}")
    #   xax.append(i)
      asdvals.append(asd)
      m2vals.append(m2)
      if m2 < max_meas2:
        max = asd
        max_meas2 = m2
    

    df = pd.DataFrame(m2vals)
    df.to_csv('nmbbh_heur_dim10.csv')

    # print(checkBasis(NB[0]))
    print(f"Max ASD : {max} --- Drift : {max_meas2}")
    print(f"Corresponding Alpha : {math.sqrt(10) * max_meas2}")

