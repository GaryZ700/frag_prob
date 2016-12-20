#Written by Gary Zeri 
#Last Edit 12/11/16

#Editing by Saswata

import networkx as nx
import numpy as np
import os
from collections import Counter
#to get number of md files use ls -l mdlog* | wc -l
#IMPORTANT FUNCTIONS#######################################################


#function to get all atoms in molecule
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

#read individual structure from mdlog
#args: struct number, starts counting at 0
#      number of atoms
#      log file number
def readstruct(i,struct_size,log):
    return str(os.popen('sed -n '+str(i*struct_size+3)+','+str((i+1)*struct_size+2)+'p ./mdlog.'+str(log)).read())
    

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


def dist(X1,X2):
    return ((X1[0]-X2[0])**2 + (X1[1]-X2[1])**2 + (X1[2]-X2[2])**2)**0.5

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



def totalMass(atoms,masses):
    tmass = 0.0
    for g in range(0,len(atoms)):
        tmass = tmass + masses[atoms[g]]    
    return tmass            

def createGraph(longX,atoms,cutoff):
    G = nx.Graph()
    for i in range(len(atoms)):
        G.add_node(i)
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            if(dist(longX[i],longX[j])<cutoff):
                G.add_edge(i,j)
    return G 

#END OF IMORTANT FUNCTIONS####################################################


#INIT IMPORTANT THINGS#######################################


#create output file

#data required
masses = {'c':21894.2, 'h':1837.29, 'o':21894.2}
cutoff = 3.0

fragments = []
for siml in range(1,120):
    #output = open('output','w')
    #location of folder containing mdlogs
    loc = './b7000m' + str(siml)
    #move into mdlog directory
    os.chdir(loc)
    if os.path.isfile("GEO_OPT_FAILED"):
        os.chdir("..")
        continue
    #get # of mdlogs
    md = int(os.popen('ls -l mdlog* | wc -l').read())
    
    #get list of atoms
    atoms = getatoms([])
    
    #get amount of lines in single struct
    struct_size = 3*len(atoms)+3
    
    #get number of structs in each md file
    struct = []
    for c in range(0,md):
        struct.append([])
        struct[c] = os.popen(str('grep "t=" "mdlog.'+str(c)+'" | wc -l')).read()
    	
    #atom velocity holder
    V = []
    for h in range(0,len(atoms)):
        V.append([])
    #V=np.array(V)
    
    #coordinates
    # X
    #Vcom = np.array([0,0,0])
    #tmass = totalMass(atoms,masses)
    
# This is where the main loop of reading the coordinates begin:
    for m in range(md,md+1):
        #struct loop
        for s in range(1,3):#int(struct[m])-1):
            X = []
            for a in range(0,len(atoms)):
                X.append(list(getx(s,atoms,struct_size,m,a)))
            G = createGraph(X,atoms,cutoff)
#           H is the list of subgraphs
            H = [list(yy) for yy in nx.connected_components(G)]
#    print H
#   ha1/2 are for converting it into atoms (so that degenerate ones are
#   taken into the same class. (In case that was not clear: atom number 3 and 4
#       are both Hydrogens. So, if the atoms are not distinguishable, except their
#       identity in atom number, then the fragments may not be distinguished.

    ha1 = []
    for h1 in H:
        ha2 = []
        for h2 in h1:
            ha2.append(atoms[h2])
        ha1.append(ha2)
    fragments = fragments + [str([list(Hh) for Hh in ha1])]
    #output.write(str(len(H)))
    os.chdir('..')
    #output.close()




freqfrag = Counter(fragments)
print freqfrag

#print fragments




