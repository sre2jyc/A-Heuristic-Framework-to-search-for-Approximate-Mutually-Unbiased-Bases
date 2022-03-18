
# importing headers

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
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

def makeBasis(M,s):
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
          AllRes.append([(B[4],SortedRes[:10])])

  return AllRes          


def getEigVec(M):
  w,v = np.linalg.eig(M)
  return v

def creatingStartingBases():
    pi = math.pi
    s = (1/math.sqrt(7))
    a = (math.cos(2*pi/7) + math.sin(2*pi/7)*1j) 
    ac = np.conj(a)
    b = (math.cos(4*pi/7) + math.sin(4*pi/7)*1j) 
    bc = np.conj(b)
    c = (math.cos(6*pi/7) + math.sin(6*pi/7)*1j) 
    cc = np.conj(c)
    M0 = np.identity(7) 
    M1 = dft(7)
    M1 = makeBasis(M1,s)
    M2 = [[1,a,cc,b,b,cc,a],[1,b,ac,bc,ac,b,1],[1,c,a,a,c,1,ac],[1,cc,c,cc,1,bc,bc],[1,bc,bc,1,cc,c,cc],[1,ac,1,c,a,a,c],[1,1,b,ac,bc,ac,b]]
    #M2 = np.array(M2) * (1/math.sqrt(7))
    #M2 = np.transpose(M2)
    M2 = makeBasis(M2,s)
    M3 = [[1,b,a,cc,cc,a,b],[1,c,c,1,a,ac,a],[1,cc,bc,c,bc,cc,1],[1,bc,1,ac,b,b,ac],[1,ac,b,b,ac,1,bc],[1,1,cc,bc,c,bc,cc],[1,a,ac,a,1,c,c]]
    #M3 = np.array(M3) * (1/math.sqrt(7))
    #M3 = np.transpose(M3)
    M3 = makeBasis(M3,s)
    M4 = [[1,c,bc,ac,ac,bc,c],[1,cc,1,b,c,c,b],[1,bc,b,bc,1,a,a],[1,ac,cc,a,cc,ac,1],[1,1,ac,cc,a,cc,ac],[1,a,a,1,bc,b,bc],[1,b,c,c,b,1,cc]]
    #M4 = np.array(M4) * (1/math.sqrt(7))
    #M4 = np.transpose(M4)
    M4 = makeBasis(M4,s)
    M5 = [[1,cc,b,a,a,b,cc],[1,bc,cc,cc,bc,1,c],[1,ac,ac,1,b,bc,b],[1,1,a,c,ac,c,a],[1,a,c,ac,c,a,1],[1,b,bc,b,1,ac,ac],[1,c,1,bc,cc,cc,bc]]
    #M5 = np.array(M5) * (1/math.sqrt(7))
    #M5 = np.transpose(M5)
    M5 = makeBasis(M5,s)
    M6 = [[1,bc,ac,c,c,ac,bc],[1,ac,a,ac,1,cc,cc],[1,1,c,b,cc,b,c],[1,a,bc,bc,a,1,b],[1,b,1,a,bc,bc,a],[1,c,b,cc,b,c,1],[1,cc,cc,1,ac,a,ac]]
    #M6 = np.array(M6) * (1/math.sqrt(7))
    #M6 = np.transpose(M6)
    M6 = makeBasis(M6,s)
    M7 = [[1,ac,c,bc,bc,c,ac],[1,1,bc,a,b,a,bc],[1,a,1,cc,ac,ac,cc],[1,b,b,1,c,cc,c],[1,c,cc,c,1,b,b],[1,cc,ac,ac,cc,1,a],[1,bc,a,b,a,bc,1]]
    #M7 = np.array(M7) * (1/math.sqrt(7))
    #M7 = np.transpose(M7)
    M7 = makeBasis(M7,s)

    M1 = np.array(M1)
    M2 = np.array(M2)
    M3 = np.array(M3)
    M4 = np.array(M4)
    M5 = np.array(M5)
    M6 = np.array(M6)
    M7 = np.array(M7)

    Bases = [M0,M1,M2,M3,M4,M5,M6,M7]

    return Bases


def constructingBoxes(StartBases):

    RB1 = dimension_reduction_svd(StartBases[0])
    RB2 = dimension_reduction_svd(StartBases[1]) 
    RB3 = dimension_reduction_svd(StartBases[2])
    RB4 = dimension_reduction_svd(StartBases[3]) 
    RB5 = dimension_reduction_svd(StartBases[4])
    RB6 = dimension_reduction_svd(StartBases[5]) 
    RB7 = dimension_reduction_svd(StartBases[6])
    RB8 = dimension_reduction_svd(StartBases[7]) 

    APB1 = allPossibleBases(RB1,6)
    APB2 = allPossibleBases(RB2,6)
    APB3 = allPossibleBases(RB3,6)
    APB4 = allPossibleBases(RB4,6)
    APB5 = allPossibleBases(RB5,6)
    APB6 = allPossibleBases(RB6,6)
    APB7 = allPossibleBases(RB7,6)
    APB8 = allPossibleBases(RB8,6)



    Boxes = [APB1,APB2,APB3,APB4,APB5,APB6,APB7,APB8]
    return Boxes



def toBasis(UpdatedBases):
  UB = []
  for i in range(4):
    bb = []
    for j in range(10):
      nv = []
      for k in range(10):
        nv.append(UpdatedBases[i][j][k][0])
      bb.append(nv)
    UB.append(bb)

  UB = np.array(UB)
  return UB 

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
    d = 1/math.sqrt(6)
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


def insideProcess(BS):
  
  nBoxSets = len(BS)
  results = []
  for i in range(nBoxSets):
    print(f"Inside {i}")
    results.append(chooseFiveBoxes(BS[i]))
    print(f"Result appended {i}")
  return results

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

    # Boxsets = []

    # for i in range(Boxes.shape[0]):
    #     for j in range(i+1,Boxes.shape[0]):
    #         for k in range(j+1,Boxes.shape[0]):
    #             for l in range(k+1,Boxes.shape[0]):
    #               for m in range(l+1,Boxes.shape[0]):
    #                 Boxsets.append(np.array([Boxes[i],Boxes[j],Boxes[k],Boxes[l],Boxes[m],[i,j,k,l,m]],dtype=object))
    


    # creating Sets of Boxsets for each Process
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

    # [1, 2, 3, 5,7], [(0.9019378204994893, [1, 4, 4, 5, 3]]

    bb1 = Boxes[1][1]
    bb2 = Boxes[2][4]
    bb3 = Boxes[3][4]
    bb4 = Boxes[5][5]
    bb5 = Boxes[7][3]

    BestBases = [bb1, bb2, bb3, bb4, bb5]

    print("NONMERGED 5SETS : DIM6")

    print(f"Average Sq Dist of Best Set of 5 Bases : {ASD(BestBases)}")


    meas2 = measure2(BestBases)
    # print((meas2))
    print(f"Max Value of Measure2 : {(meas2)}")

    NewBases = descent(BestBases,0.5)
    NB = NewBases
    max = ASD(NewBases)
    max_meas2 = measure2(NewBases)
    # max_meas2 = max(max_meas2)
    asdvals = []
    asdvals.append(max)
    # xax = []
    for i in range(10000):
      NB = descent(NB,0.5)
      asd = (ASD(NB))
      m2 = measure2(NB)
      print(f"Iter no : {i} --> ASD : {asd}  Drift : {m2}")
    #   xax.append(i)
      asdvals.append(asd)
      if asd > max:
        max = asd
        max_meas2 = m2
    

    df = pd.DataFrame(asdvals)
    df.to_csv('nmbbh_heur_dim6_(0.5).csv')

    # print(checkBasis(NB[0]))
    print(f"Max ASD : {max} --- Drift : {max_meas2}")