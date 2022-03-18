
# importing headers

import numpy as np
from matplotlib import pyplot as plt
# from sklearn.decomposition import TruncatedSVD
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
#   print(type(A))
  n = A.shape[0]
  Sum = 0
  for i in range(n):
    V1 = np.transpose(np.array(A[i]))
    for j in range(n):
      V2 = np.transpose(np.conj(np.array(B[j])))
      # print(type(V2))
      dP = V1 * V2
      # print(type(dP))
      dP = pow(abs(sum(dP)),2)
      #print(str(i) + '|' + str(j), dP)
      Sum += (dP * (1-dP))


  Distance = Sum / (n-1)    


  return Distance

# calculates the average squared distance between a set of bases
def ASD(Bases):
  Bases = np.array(Bases)
  # print(f"Shape : {Bases.shape}")
  # print(Bases)
  #number of bases
  n = 4
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
                  #ASDs.append(asd)
                  point = [(a),(b),(c),(d)]
                  
                  print(f"Point : {point} --> ASD : {asd}")
                  Results.append((asd,point))

  return Results                  


def chooseFourBoxes(B):
  AllRes = []
  nB = 4
  for i in range(nB):
    for j in range(i+1,nB):
      for k in range(j+1,nB):
        for l in range(k+1,nB):
          choosedBoxes = np.array([B[i],B[j],B[k],B[l]])
          Results = searchFourBoxes(choosedBoxes)
          SortedRes = Results
          SortedRes.sort(reverse = True)
          AllRes.append([(B[4],SortedRes[:10])])

  return AllRes          

def searchFiveBoxes(B):
  nG = B.shape[0]
  nB = B.shape[1]
  print("Inside sFb")
  Results =  []
  BG0 = B[0]
  BG1 = B[1]
  BG2 = B[2]
  BG3 = B[3]
  BG4 = B[4]
  for a in range(nB):
    for b in range(nB):
      for c in range(nB):
        for d in range(nB):
          for e in range(nB):
            B0 = BG0[a]
            B1 = BG1[b]
            B2 = BG2[c]
            B3 = BG3[d]
            B4 = BG4[e]

            asd = ASD([B0,B1,B2,B3,B4])
            point = [a,b,c,d,e]
            # print("Res sFb")
            Results.append((asd,point))
  print("Res sFb")
  return Results


def chooseFiveBoxes(B):
  AllRes = []
  nB = 5

  choosedBoxes = np.array([B[0],B[1],B[2],B[3],B[4]])
  print("Inside cFB")
  Results = searchFiveBoxes(choosedBoxes)
  print("Res cFB")
  SortedResults = Results
  SortedResults.sort(reverse=True)
  AllRes.append([(B[5],SortedResults[:10])])

  return AllRes


def getEigVec(M):
  w,v = np.linalg.eig(M)
  return v

def creatingStartingBases():
    pi = math.pi
    o = math.cos(2*pi/47) + (math.sin(2*pi/47)*(1j))
    Z47 = [
            [
            o**0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, o**1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, o**2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, o**3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, o**4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, o**5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, o**6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, o**7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, o**8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, o**9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**37, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**38, 0, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**39, 0, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**40, 0, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**41, 0, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**42, 0, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**43, 0, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**44, 0, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**45, 0
            ],
            [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, o**46
            ]
        ]    
    Z47 = np.array(Z47)
    X47 =  [
        [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j
            ],
            [
            1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j, 0j
            ],
            [
            0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 0j, 1 + 0j, 0j
            ]
        ]
    
    X47 = np.array(X47)
    M1 = Z47
    M2 = X47
    # M3 = np.dot(M2,M1)
    # M4 = np.dot(M3,M1)
    # M5 = np.dot(M4,M1)
    # M6 = np.dot(M5,M1)
    # M7 = np.dot(M6,M1)
    # M8 = np.dot(M7,M1)
    # M9 = np.dot(M8,M1)
    # M10 = np.dot(M9,M1)
    # M11 = np.dot(M10,M1)
    # M12 = np.dot(M11,M1)
    Matrices = []

    Matrices.append(M1)
    Matrices.append(M2)
    temp = M2
    for i in range(2,49):
        temp2 = np.dot(M1,temp)
        Matrices.append(temp2)
        temp = temp2
    

    Bases = []

    for i in range(0,48):
        Bases.append(np.transpose((getEigVec(Matrices[i]))))
    

    # B01 = np.transpose(getEigVec(M1))
    # B02 = np.transpose(getEigVec(M2))
    # B03 = np.transpose(getEigVec(M3))
    # B04 = np.transpose(getEigVec(M4))
    # B05 = np.transpose(getEigVec(M5))
    # B06 = np.transpose(getEigVec(M6))
    # B07 = np.transpose(getEigVec(M7))
    # B08 = np.transpose(getEigVec(M8))
    # B09 = np.transpose(getEigVec(M9))
    # B10 = np.transpose(getEigVec(M10))
    # B11 = np.transpose(getEigVec(M11))
    # B12 = np.transpose(getEigVec(M12))

    # Bases = [B01,B02,B03,B04,B05,B06,B07,B08,B09,B10,B11,B12]

    return Bases


def constructingBoxes(StartBases):

    # MergedBases = mergeBases(StartBases,48,47)
    # MergedBases = np.array(MergedBases)

    # ReducedBases = dimension_reduction_svd(MergedBases)
    # ReducedBases = np.array(ReducedBases)
    RedBases = []
    for i in range(len(StartBases)):
      RedBases.append(dimension_reduction_svd(StartBases[i]))

    Boxes = []
    for i in range(len(RedBases)):
        Boxes.append(allPossibleBases(RedBases[i],46))

    # Boxes = [B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11]
    return Boxes

def insideProcess(BS):
  n = 4
  nBoxSets = len(BS)
  results = []
  for i in range(nBoxSets):
    if n == 4 :
      results.append(chooseFourBoxes(BS[i]))
      print(f"Result appended {i}")

    elif n == 5:
      print(f"Inside {i}")

      results.append(chooseFiveBoxes(BS[i]))
      print(f"Result appended {i}")
  
  return results



# def Grad(B1,B):
#   Msum = 0
#   for b2 in range(4):
#   #B1 = B[b1]
#     B2 = B[b2]
#     for i in range(10):
#       for j in range(10):
#         m1 = np.matrix(B1[i])
#         m2 = np.matrix(B2[j])
#         m1c = np.matrix(np.conjugate(B1[i]))
#         m2c = np.matrix(np.conjugate(B2[j]))
#         M = np.transpose(m1) * m1c * np.transpose(m2) * m2c 
        
#         M = np.dot(M,M)
#         Msum += M


#   iM = np.imag(Msum)
#   G0 = (8/(4 * 3 * 9)) * iM   
#   return np.array(G0)  

# def calV(sig):
#   expo = sig * 1j
#   val = 1 + expo + (expo**2/math.factorial(2)) + (expo**3/math.factorial(3)) + (expo**4/math.factorial(4)) + (expo**5/math.factorial(5))
#   return val

# def descent(B, step):
#   G0 = Grad(B[0],B)
#   G1 = Grad(B[1],B)
#   G2 = Grad(B[2],B)
#   G3 = Grad(B[3],B)
   
#   G = np.array([G0,G1,G2,G3])
#   # print('Gradients are : ')
#   # print(G)
#   Sig = G
#   for i in range(4):
#     Sig[i] = step * Sig[i]


#   # print('Sigma is : ' + Sig)
#   # print("Sigma Value : ")
#   # print(Sig)

#   V = []
  
#   for i in range(4):
#     # II = np.identity(6) 
#     II = expm(Sig[i] * 1j)
#     V.append(II)

#   V = np.array(V)  

#   # print("Value of V : ")
#   # print(V)    

#   uB = []

#   for b in range(4):
#     tB = B[b]
#     nB = []
#     for i in range(10):
#       m = np.matrix(tB[i])
#       nv = np.dot((V[b]), np.transpose(m))
#       nB.append(nv)

#     uB.append((nB))

#   #uB.append(B[3])
#   uB = np.array(uB)
#   # print("Updated Bases ")
#   # print(uB)
  
#   return uB



# def toBasis(UpdatedBases):
#   UB = []
#   for i in range(4):
#     bb = []
#     for j in range(10):
#       nv = []
#       for k in range(10):
#         nv.append(UpdatedBases[i][j][k][0])
#       bb.append(nv)
#     UB.append(bb)

#   UB = np.array(UB)
#   return UB 


# def RUN(SB,iter,step):
#   StartingBases = SB
#   maxDis = -1
#   bestBases = SB
#   allAsdVal = []
#   # print('Start Base ASD : ')
#   # print(ASD(SB))
#   for i in range(iter):
#     UpdatedBases = descent(StartingBases,step)
#     UB = toBasis(UpdatedBases)
#     # print('ASD : ')
#     asdVal = ASD(UB)
#     allAsdVal.append(asdVal)
#     print(f"Updates {i} ASD : {asdVal}")
#     if maxDis < asdVal :
#         maxDis = asdVal
#         bestBases = UB
#     StartingBases = UB

#   return [allAsdVal,StartingBases]

def distancePrinter(A,B):
  n = A.shape[0]
  Sum = 0
  for i in range(n):
    V1 = np.transpose(np.array(A[i]))
    for j in range(n):
      V2 = np.transpose(np.conj(np.array(B[j])))
      dP = V1 * V2
      dP = pow(abs(sum(dP)),2)
      print(str(i) + '|' + str(j), dP)
      Sum += (dP * (1-dP))


  Distance = Sum / (n-1)    


  return Distance

def modifiedMUBTest(Bases):
  Bases = np.array(Bases)
  #number of bases
  n = Bases.shape[0]
  sum = 0
  for i in range(n):
    for j in range(i+1,n):
      print(f"For Base Pair {i} and {j}, the vector products are : ")
      distance_pair = distancePrinter(Bases[i],Bases[j])
      sum += float(distance_pair)

  sum = sum * 2
  sum = sum / (n*(n-1))

  return sum


def measure2(Bases):
    nB = len(Bases)
    d = 1/math.sqrt(46)
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

    return s

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


    print(f"Shape of Boxes : {Boxes.shape}")

    # Boxsets = []

    # for i in range(Boxes.shape[0]):
    #     for j in range(i+1,Boxes.shape[0]):
    #         for k in range(j+1,Boxes.shape[0]):
    #             for l in range(k+1,Boxes.shape[0]):
    #                Boxsets.append(np.array([Boxes[i],Boxes[j],Boxes[k],Boxes[l],[i,j,k,l]],dtype=object))
    


    # # creating Sets of Boxsets for each Process
    # Boxsets = Boxsets[:16]
    # BS = []
    # nBs = math.ceil(len(Boxsets)/16)
    # for i in range(nBs):
    #   if i == nBs-1:
    #     BS.append(Boxsets[i*16:])
    #   else:  
    #     BS.append(Boxsets[i*16:(i+1)*16])



    # print(len(BS))
    # print(len(BS[0]))
    # print(len(BS[-1]))

    # start = time.time()
    # # Boxsets = Boxsets[:16]
    # nProcess = 16
    # print(nProcess)
    # processPool = Pool(nProcess)
    # results = processPool.map(insideProcess, BS)
    # processPool.close()
    # processPool.join()
    # end = time.time()
    # print(f"Time taken {end-start}secs")
    # print(results)

    bb1 = Boxes[0][0]
    bb2 = Boxes[1][34]
    bb3 = Boxes[2][33]
    bb4 = Boxes[3][5]

  

    BestBases = [bb1,bb2,bb3,bb4]

    print("NONMERGED 4SETS : DIM46")

    print(f"Average Sq Dist of Best Set of 4 Bases : {ASD(BestBases)}")


    meas2 = measure2(BestBases)
    # print((meas2))
    print(f"Max Value of Measure2 : {max(meas2)}")



    # [0, 34, 33, 5]



    
