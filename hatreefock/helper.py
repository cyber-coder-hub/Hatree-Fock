#Read the geometry
import numpy as np
from math import sqrt
#########################################################
#Function to define no of Electron
#On Entry 
#element--> name of the element
#On Exit
#u--> atomic no
########################################################

def no_of_e(element):
  symbol = [
            'H','He',
            'Li','Be','B','C','N','O','F','Ne',
            'Na','Mg','Al','Si','P','S','Cl','Ar',
            'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
            'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
            'Rh', 'Pd', 'Ag', 'Cd',
            'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',  'Eu',
            'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
            'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl','Pb','Bi','Po','At','Rn']
  u=symbol.index(element)+1
  return u
##########################################
# Calculate the distance between two atoms
#On Entry 
#a--> position of first atom, in list containing x,y and z
#b--> position of second atom
#On return
#R--> distance between two atoms
##########################################################
def find_distance(a,b):
 x_1=float(a[0])
 x_2=float(b[0])
 y_1=float(a[1])
 y_2=float(b[1])
 z_1=float(a[2])
 z_2=float(b[2])
 R_square =(x_1-x_2)**2+(y_1-y_2)**2+(z_1-z_2)**2
 R=sqrt(R_square)
 
 return R
##########################################################
#This function reads the 1 electron Hamiltonian from the disk
#
#On Entry
#file_name--> Name of the file to be read from
#basis--> Dimension of the basis set
#
#On Exit
#A-->Numpy array with the 1 electron Hamiltonian elements(nbasis,nbasis)
##########################################################
def file_read_1e(file_name,nbasis):

 #open the file
 input_file=open(file_name)   #open the file
 #read the file using readline to file_content
 file_content=input_file.readlines()# read the content
 #close the file
 input_file.close() 
 A=np.zeros([nbasis,nbasis])

 for line in file_content:
    V_line=line.rstrip()
    V_line=V_line.split() 
    i=int(V_line[0])-1
    j=int(V_line[1])-1
    A[i][j]=float(V_line[2])
    A[j][i]=float(V_line[2])
 return A
 #########################################################
#This function reads the 2 electron Hamiltonian from the disk
#
#On Entry
#file_name--> Name of the file to be read from
#basis--> Dimension of the basis set
#
#On Exit
#twoe-->Numpy array with the 2 electron integrals (nbasis,nbasis,nbasis,nbasis)
##########################################################
def read_2_e(file,nbasis):

#read the file
 input_file=open(file)   #open the file

 file_content=input_file.readlines()# read the content

 input_file.close()            # close the file

#print(file_content)
 twoe_index=[]
 twoe_value=[]

 for line in file_content:
    V_line=line.rstrip()
    V_line=line.split()
    i=int(V_line[0])
    j=int(V_line[1])# change from Dirac to Muliken Ordering
    k=int(V_line[2])
    l=int(V_line[3])
    ijkl=compound_index(i,j,k,l)
    twoe_index.append(ijkl)
    twoe_value.append(V_line[4])
    
 #define the 4 d array
 twoe=np.zeros([nbasis,nbasis,nbasis,nbasis])

 for i in range(nbasis):
   for j in range(nbasis):
     for k in range(nbasis):
       for l in range(nbasis):
           ijkl=compound_index(i+1,j+1,k+1,l+1)
           if ijkl in twoe_index:
                ind=twoe_index.index(ijkl)
                twoe[i,j,k,l]=float(twoe_value[ind])
                
           
 return twoe
 
#####################################################
# The compound index code 
#
#on Input 
#i,j,k,l --> four integers
#
#On Return 
#ijkl-->gives the compound index of I,J,K,L
#
#####################################################
def compound_index(i,j,k,l):

  if i>j:
   ij=i*(i+1)/2+j
  else:
   ij=j*(j+1)/2+i

  if k>l:
   kl=k*(k+1)/2+l
  else:
   kl=l*(l+1)/2+k

  if ij>kl:
   ijkl=ij*(ij+1)/2+kl
  else:
   ijkl=kl*(kl+1)/2+ij

  return ijkl

#########################################################
#This function calculate X=S^-1/2
#
#On input
#S--> overlap matrix dimension (nbasis,nbasis)
#
#On output
#X--> the transformation matrix dimension (nbasis,nbasis)
##########################################################
def get_X(matrix,nbasis):
  x=np.zeros([nbasis,nbasis])
  xTemp=np.zeros([nbasis,nbasis])
  temp=np.zeros([nbasis,nbasis])
  lambda_b,L=np.linalg.eigh(matrix) #genarating eigen values and eigen vector
  #print(lambda)
  #print(L)
  for i in range(nbasis):
    temp[i][i]=(lambda_b[i])**(-0.5)
  xTemp=np.matmul(L,temp)
  x=np.matmul(xTemp,L.transpose())
  return x   
        
    
 
