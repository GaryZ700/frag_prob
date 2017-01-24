# Saswata Roy writes this.

import networkx as nx
import numpy as np
import os
from collections import Counter
import argparse as ap
import math

parser = ap.ArgumentParser(description="Analysis of the fragmentation in MD from mdlog.x")
parser.add_argument('-f',"--filehead", metavar='filehead',type=str, default="b4500m")
parser.add_argument('-s',"--start", metavar='start',type=int, default=1)
parser.add_argument('-l',"--last", metavar="last",type=int,default=10)

args = parser.parse_args()
filehead = args.filehead
start = args.start
last = args.last

# The necessary functions:

#returns list of atom symbols ex. 'c', 'h', 'li'
def getatoms(atoms):
    #open first mdlog
    f = file('./mdlog.1')
    
    #get number of atoms from md.out
    #may be unstable due to location of first 'atoms'
    atom_num = int(os.popen('grep "atoms" "md.out"').read().split('atoms')[0])
    
    #read all lines contaning atom symbol
    for i in range(0,atom_num+3):
        line = f.readline().split(' ')
        if(i>2):
            atoms.append(line[len(line)-1].split('\n')[0])
    return atoms   
#****************************************************************

#get velocity of specific atom at specific struct
def getv(struct_num,atoms,struct_size,log,atom_num):
    f = file('work','w')
    f.write(str(readstruct(struct_num,struct_size,log)))
    f.close()

    lines = os.popen('sed -n '+str(len(atoms)+3)+','+str(2*len(atoms)+2)+'p ./work').read().split('\n')
    line = lines[atom_num].split(' ')
    line = list(filter(None,line))
    
    x = float(line[1])
    y = float(line[2])
    z = float(line[3])

    return [x,y,z]
#****************************************************************

def dist(X1,X2):
    return ((X1[0]-X2[0])**2 + (X1[1]-X2[1])**2 + (X1[2]-X2[2])**2)**0.5
#***************************************************************
def readstruct(i,struct_size,log):
    return str(os.popen('sed -n '+str(i*struct_size+3)+','+str((i+1)*struct_size+2)+'p ./mdlog.'+str(log)).read())

#get coordinates of specific atom at specific struct
def getx(struct_num,atoms,struct_size,log,atom_num):
    f = file('work','w')
    f.write(str(readstruct(struct_num,struct_size,log)))
    f.close()

    lines = os.popen('sed -n 2'','+str(len(atoms)+1)+'p ./work').read().split('\n')
    line = lines[atom_num].split(' ')
    line = list(filter(None,line))
     
    x = float(line[0])
    y = float(line[1])
    z = float(line[2])

    return [x,y,z]
#****************************************************************

def totalMass(atoms,masses):
    tmass = 0.0
    for atom in atoms:
    
        tmass = tmass + masses[atoms[atom]]    
    return tmass            
#****************************************************************

def fragMass(atoms,masses,frag):
    tmass = 0.0
    for atom in frag:
        tmass = tmass + masses[atoms[atom]]    
    return tmass            
#****************************************************************

def createGraph(longX,atoms,cutoff):
    G = nx.Graph()
    for i in range(len(atoms)):
        G.add_node(i)
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            if(dist(longX[i],longX[j])<cutoff):
                G.add_edge(i,j)
    return G 
#****************************************************************

#get number of structs in mdlog file
def getStructs(md):
    struct = []
    
    for c in range(1,md+1):
        struct.append([])
        struct[c-1] = int(os.popen(str('grep "t=" "mdlog.'+str(c)+'" | wc -l')).read())
    
    return struct

#****************************************************************

#user interface for accessing frag database

def gui(dB,fragTypes,atoms):
    
    opt = [1]    
    sp = "    "


    print("\n")
    print("Molecule: " + str(atoms))
    print("Fragment Versions: "+ str(len(fragTypes)))
    print("\n")
    
    for n in range(0,len(fragTypes)):
        print(sp + numtoatmname(atoms,fragTypes[n]) + " Frequency: " + str(fragTypes[n][len(fragTypes[n])-1])) 

    print("1. Get KE")
    print("2. Get RE")
    print("3. Get VE")

    inp = raw_input("Please choose an option.")
    
    if(inp == 1):
        for n in range(0,len(fragTypes)):
            print(sp + numtoatmname(atoms,fragTypes[n]) + " Frequency: " + str(fragTypes[n][len(fragTypes[n])-1])) 

    
        





#****************************************************************

def numtoatmname(atoms,obj):
    
    finstr = ""
    
    for frag in obj:
        holder = ""
            
        if(not isinstance(frag,int)):
            for atom in frag:        
                holder = holder + " " + atoms[atom]
        finstr = finstr + "   " + holder 

    return finstr


#****************************************************************

def constructI(frag,atoms,X):
    I = [[0.0 for y in range(3)] for x in range(3)]
    for i in range(3):
        for atom in frag:
            I[i][i] += masses[atoms[atom]]*(X[atom][(i+1)%3]**2 + X[atom][(i+2)%3]**2)
        for j in range(3):
            if(i==j):
                continue
            I[i][j] += -masses[atoms[atom]]*(X[atom][i]*X[atom][j])
    return I

 
         
#****************************************************************
#calculates the total KE for each frag
def TKEfrag(frag,longV,atoms):
    KE = 0.0
    for atom in frag:
        for dim in range(3):
            KE += .5*masses[atoms[atom]]*longV[atom][dim]*longV[atom][dim]
    return KE          







#get number of fragmentation types in db
def getFT(dB):
    
    flist = []
    repeat = []
    
    
    flist.append(dB[0][1])
    repeat.append([1])
    
    for b in range(1, len(dB)):
        
        unique = True
        
        for c in range(0,len(flist)):
            if(flist[c] == dB[b][1]):
                unique = False
                repeat[c].append(1)
                
        if(unique):
            flist.append(dB[b][1])
            repeat.append([1])

    for h in range(0, len(repeat)):
        flist[h].append(len(repeat[h]))

    
    return flist    




#****************************************************************

#Calculates the rotational energy for each frag
def RKEfrag(frag,longX,longV,atoms,tolerance):
    newV,newX =  newXV(longX,longV,frag,atoms)
    
    I = constructI(frag,atoms,newX)
    Id,U = np.linalg.eig(I)
    rke = 0.0
    Vpp = np.zeros(3)
    Xpp = np.zeros(3)
  
    for atom in frag:
        Vpp = np.zeros(3)
        Xpp = np.zeros(3)
  

        for dim in range(3):
            Vpp[dim] = np.linalg.norm(np.cross(U[dim],newV[atom]))
            Xpp[dim] = np.linalg.norm(newX[atom] - (np.dot(U[dim],newX[atom]))*U[dim])
            if(Xpp[dim] > tolerance ):
                rke += .5*Id[dim]*((Vpp[dim]/Xpp[dim])**2)

#****************************************************************
def newXV(longX,longV,frag,atoms):
    Xcom = np.zeros(3)
    Vcom = np.zeros(3)
    tmass = fragMass(atoms,masses,frag)
    for atom in frag:
        for dim in range(3):
            Vcom[dim] += (longV[atom][dim]*masses[atoms[atom]])/tmass
            Xcom[dim] += (longX[atom][dim]*masses[atoms[atom]])/tmass
    
    
    newX = np.empty([len(longV),3])
    newV = np.empty([len(longV),3])
    
    for i in range(len(longX)):
        newX[i] = longX[i] - Xcom
        newV[i] = longV[i] - Vcom   
    
    return newV,newX 


#***************************************************************
masses = {'c':21894.2, 'h':1837.29, 'o':21894.2}
cutoff = 3.0
fragments = []
os.chdir(filehead + str(start))
atoms = getatoms([])
struct_size = 3*len(atoms)+3
os.chdir("..")

struct = []
KEnum = 0
KEsum = []
KEd = []
KEdc = 0
dB = [] #the database. We shall be playing with this a bit.
tolerance = 0.00001
#***************************************************************

for siml in range(start,last+1):
    
    #output = open('output','w')
    os.chdir(filehead + str(siml))
    if os.path.isfile("GEO_OPT_FAILED"):
        os.chdir("..")
        continue
    
    #get # of mdlogs
    md = int(os.popen('ls -l mdlog* | wc -l').read())
    #get number of structs in each md find
    nstruct = getStructs(md)
    #This is where the main loop of reading the coordinates begin:
    for m in range(md,md+1):
        #struct loop
        for s in range(1,3):#int(struct[m])-1):
            X = []
            V = []
            for a in range(0,len(atoms)):
                X.append(list(getx(s,atoms,struct_size,m,a)))
                V.append(list(getv(s,atoms,struct_size,m,a)))
            G = createGraph(X,atoms,cutoff)
            #H is the list of subgraphs
            H = [list(yy) for yy in nx.connected_components(G)]

    nn = len(H)
    Vcom = [[0.0 for y in range(3)] for x in range(nn)]
    KEcom = [0.0 for y in range(nn)]
    cou = 0
    RKEd = []
    TKEd = []
    for frag in H:
        tmass = fragMass(atoms,masses,frag)
        RKEd.append(RKEfrag(frag,X,V,atoms,tolerance)
        TKEd.append(TKEfrag(frag,V,atoms))
  
        #This part calculates the KEcom for the frag 
        for atom in frag:
            for dim in range(3):
                Vcom[cou][dim] += (V[atom][dim]*masses[atoms[atom]])/tmass
            for dim in range(3):
                KEtemp = 0.0
                KEtemp += 0.5*(V[atom][dim]*V[atom][dim])*tmass
                KEcom[cou] += KEtemp
        cou += 1
    listemp = []
    listemp.append(nn)
    listemp.append(H)
    listemp.append(Vcom)
    listemp.append(KEcom)
    listemp.append(RKEd)
    listemp.append(TKEd)
    dB.append(listemp)

    print("End of a siml")
    os.chdir('..')

    ha1 = []
    for h1 in H:
        ha2 = []
        for h2 in h1:
            ha2.append(atoms[h2])
        ha1.append(ha2)
    fragments = fragments + [str([list(Hh) for Hh in ha1])]

#for line in dB:
#    print(line)

fragTypes = getFT(dB)
#print str(fragTypes)


gui(dB,fragTypes,atoms)




