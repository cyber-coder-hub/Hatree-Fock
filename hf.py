#THIS IS A PROGRAM TO CALCULATE H-F ENERGYOF A CLOSED SHELL MOLECULE
import numpy as np
def atomicNumber(atomSymbol):
     #print(atomSymbol)
     symbol = [
            'NONE','H','He',
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
     z=symbol.index(atomSymbol)
     #print(z)
     return(z)

#READ THE GEOMETRY
inputFile=open('geom.dat')   #opening file
fileContent=inputFile.readlines()  #Reading content
inputFile.close()        #closing file
#print(fileContent)

#Storing the input in the list
tempGeom=[]
for line in fileContent:
    vLine=line.rstrip()     # to remove the trailing spaces
    if len(vLine)>0:   #to check the white space line and remove them
     tempGeom.append(vLine.split())
#print (tempGeom)
#print (tempGeom[1][0])

#No of atoms
nAtom=int(tempGeom[0][0])
print ('Total number of atoms present in the system: '+str(nAtom))
atomSymbol=[]     #this list will contain atom symbols
geom=[]    #this will contain cartesian cordinate
for i in range(1,nAtom+1):    
    atomSymbol.append(tempGeom[i][0])   #storing atom symbols in new list
    geom.append(tempGeom[i][1:])
print (atomSymbol)
print('cartesian co-ordiante in a.u.')
print (geom)

# storing total no of electrons in system
Z=[]     #storing no of electron in each atom
for i in range(nAtom):
    #print(i,atomSymbol[i])
    Z.append(atomicNumber(atomSymbol[i]))
#print (Z)
totalNoOfElectrons=np.sum(Z)
print('total no of elcetrons in the system :'+str(totalNoOfElectrons))




#calculate the nuclear repusion energy
from distance import find_distance
eNuclear=0.0
for i in range(nAtom):
    for j in range(0,i):
        zA=atomicNumber(atomSymbol[i])
        zB=atomicNumber(atomSymbol[j])
        rAB=find_distance(geom[i],geom[j])
        eNuclear+=((zA*zB)/rAB)
print('nuclear repulsion energy: ' +str(eNuclear)+ ' a.u.')



#defininig dimension of basis set
nBasis=7
#read one electron integrals
from readFile import readFile
#read s
s=readFile('s.dat',nBasis)
print(s)
#read kinetic energy
t=readFile('t.dat',nBasis)
#print(t)
#read potential energy
v=readFile('v.dat',nBasis)
#print(v)

#h core
import numpy as np
hCore=np.zeros([nBasis,nBasis])
hCore=np.add(t,v)
print('#####hCore######')
print(hCore)

