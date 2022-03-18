import numpy as np
import sympy
from sympy import Matrix, Symbol
from sympy.physics.quantum import TensorProduct

def tensorproduct(A,B):
  A = np.array(A)
  B = np.array(B)
  C = TensorProduct(A,B)
  return C

A = np.array([[1,2],[3,4]])

B = np.array([[1,1],[1,1]])

T = tensorproduct(A,B)

# print(T)

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

def getEigVec(M):
  w,v = np.linalg.eig(M)
  return v



px = np.array([[0,1],[1,0]])
py = np.array([[0,-1j],[1j,0]])
pz = np.array([[1,0],[0,-1]])
I = np.identity(2)


m1 = TensorProduct(I,pz)
m2 = TensorProduct(I,px)
m3 = TensorProduct(I,py)
m4 = TensorProduct(pz,px)
m5 = TensorProduct(pz,py)

bases = [m1,m2,m3,m4,m5]
print(bases[3])
# print(ASD(bases))

# print(m1)
# print(distance(np.ones((4,4)),np.identity(4)))


b1 = np.transpose(getEigVec(m1))
b2 = np.transpose(getEigVec(m2))
b3 = np.transpose(getEigVec(m3))
b4 = np.transpose(getEigVec(m4))
b5 = np.transpose(getEigVec(m5))

print(ASD([b1,b2]))
