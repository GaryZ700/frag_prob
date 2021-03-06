#Written by Gary Zeri 
#Last Edit 12/11/16

#Editing by Saswata

import networkx as nx
import numpy as np
import os
from collections import Counter
import argparse as ap
import math

#starting arguments
#parser = ap.ArgumentParser(description="Analysis of the fragmentation in MD from mdlog.x")

#parser.add_argument('-f',"--filehead", metavar='filehead',type=str, default="b4500m")
#parser.add_argument('-s',"--start", metavar='start',type=int, default=1)
#parser.add_argument('-l',"--last", metavar="last",type=int,default=10)


#args = parser.parse_args()

#filehead = args.filehead
#start = args.start
#last = args.last



#IMPORTANT FUNCTIONS#######################################################


#function to get all atoms in molecule
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
#args: structure number, list of atoms, size of each struct, , number of atom to find
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

    
    #returns list of XYZ coords of particular atom
    return [x,y,z]



def totalMass(atoms,masses):
    tmass = 0.0
    for atom in atoms:
    
        tmass = tmass + masses[atoms[atom]]    
    return tmass            



def fragMass(atoms,masses,frag):
    tmass = 0.0
    for atom in frag:
        tmass = tmass + masses[atoms[atom]]    
    return tmass            


#create nx graph of atoms at respective pos
#args: pos, list of atoms, dist at which bond is "broken"
def createGraph(longX,atoms,cutoff):
    G = nx.Graph()
    for i in range(len(atoms)):
        G.add_node(i)
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            if(dist(longX[i],longX[j])<cutoff):
                G.add_edge(i,j)
    return G 

#get number of structs in mdlog file
def getStructs(md):
    struct = []
    
    #for md files in each siml
    for c in range(1,md+1):
        
        struct.append([])
        struct[c-1] = int(os.popen(str('grep "t=" "mdlog.'+str(c)+'" | wc -l')).read())
    
    #returns array that has #of structs in each md file
    return struct
  


def translate(atoms,frag):
    
    holder = ""

    for atom in frag:
        holder = holder + " " + atoms[atom]
   
    return holder


#return standard deviation of data
def STDev(data,avg):
    print(str(data))   
    Dsum = 0.0

    for val in data:
        Dsum = Dsum + (val-avg)*(val-avg) 
    
    return math.sqrt((Dsum/len(data)))

#get number of MD files in folder
def getMD:
    md = int(os.popen('ls -l mdlog* | wc -l').read())

def getFrags(md,struct,atoms,struct_size,cutoff):
    for m in range(1,md+1):
    #struct loop
        for s in range(0,int(struct[m])):
            X = []
            for a in range(0,len(atoms)):
                X.append(list(getx(s,atoms,struct_size,m,a)))
            G = createGraph(X,atoms,cutoff)
            #H is the list of subgraphs
            H = [list(yy) for yy in nx.connected_components(G)]
 

def getVcom(atoms,masses,H):
    for frag in H:
        tmass = fragMass(atoms,masses,frag)
        vtemp = np.array([0.0,0.0,0.0])
        
        for atom in frag:
            vtemp =  vtemp + ((masses[atoms[atom]] * np.array(getv((struct[md-1]-10),atoms,struct_size,md,atom))/tmass))

        Vcom.append(vtemp)
    	
        

    #for prnt in range(0,len(H)):
        #print("Frag:" + str(H[prnt]) + "  Center of Mass Velocity" + str(Vcom[prnt]))
        print("Frag:" + str(frag) + "  Center of Mass Velocity" + str(vtemp))
    
    
        print("\n")
    #get COM KE for frags
    
    return Vcom
 
def getKEcom(H,Vcom):    
    
    KEcom = []

    for frag in range(0,len(H)):
        tmass = fragMass(atoms,masses,frag)
        
        vsq = 0.0
        for a in range(0,3):
            #v = v + (Vcom[c][a]*Vcom[c][a])
            vsq = vsq + (Vcom[frag][a]*Vcom[frag][a])
            # Not sure if the equation is correct
            # I just put the v=0.0 INSIDE the loop of frag.
            # This is why you were getting some error.
    
        KEcom.append([])
        KEcom[c] = .5*tmass*vsq
    
    return KEcom
 




#



#END OF IMORTANT FUNCTIONS#####################################################


#INIT IMPORTANT THINGS#######################################


#data required
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

#START DOING STUFF HERE#############################################################

for siml in range(start,last+1):
    
    #output = open('output','w')
    
    os.chdir(filehead + str(siml))
    if os.path.isfile("GEO_OPT_FAILED"):
        os.chdir("..")
        continue
    
    
    #get # of mdlogs
    md = int(os.popen('ls -l mdlog* | wc -l').read())
    #get number of structs in each md find
    struct = getStructs(md)

    
    #This is where the main loop of reading the coordinates begin:
    #for md file in siml
    for m in range(1,md+1):
        
        #for struct in md file
        for s in range(0,int(struct[m])):
         
            X = []
            
            #for atom in struct
            for a in range(0,len(atoms)):
                #create list of atom positions
                X.append(list(getx(s,atoms,struct_size,m,a)))
            G = createGraph(X,atoms,cutoff)
            #H is the list of subgraphs
            H = [list(yy) for yy in nx.connected_components(G)]
    #print loc, H
    #   ha1/2 are for converting it into atoms (so that degenerate ones are
    #   taken into the same class. (In case that was not clear: atom number 3 and 4
    #   are both Hydrogens. So, if the atoms are not distinguishable, except their
    #   identity in atom number, then the fragments may not be distinguished.

    ha1 = []
    for h1 in H:
        ha2 = []
        for h2 in h1:
            ha2.append(atoms[h2])
        ha1.append(ha2)
    fragments = fragments + [str([list(Hh) for Hh in ha1])]
    #output.write(str(len(H)))
    #vcom [frag] is center of mass velocity for specific frag
    Vcom = []					
    
    #move into last molecule folder
    #os.chdir("./"+str(filehead)+str(last))
    
    #get number of smd files in last folder
    #md = int(os.popen('ls -l mdlog* | wc -l').read())
    
    #get number of structs in each md file
    #struct = getStructs(md) 
    
    KEcom = []
    c = 0
    
    for frag in H:
        tmass = fragMass(atoms,masses,frag)
        vtemp = np.array([0.0,0.0,0.0])
        
        for atom in frag:
            vtemp =  vtemp + ((masses[atoms[atom]] * np.array(getv((struct[md-1]-10),atoms,struct_size,md,atom))/tmass))

        Vcom.append(vtemp)
    
    #for prnt in range(0,len(H)):
        #print("Frag:" + str(H[prnt]) + "  Center of Mass Velocity" + str(Vcom[prnt]))
        print("Frag:" + str(frag) + "  Center of Mass Velocity" + str(vtemp))
    
    
        print("\n")
    #get COM KE for frags
    
    
    #for frag in H:
        tmass = fragMass(atoms,masses,frag)
        
        vsq = 0.0
        for a in range(0,3):
            #v = v + (Vcom[c][a]*Vcom[c][a])
            vsq = vsq + (vtemp[a]*vtemp[a])
        #Not sure if the equation is correct
        # I just put the v=0.0 INSIDE the loop of frag.
        # This is why you were getting some error.
    
        KEcom.append([])
        KEcom[c] = .5*tmass*vsq
        c = c + 1
    

    #prepare length of KEsum
    if(len(KEsum) != len(KEcom)):
        for l in range(0,len(KEcom)):
            KEsum.append([])
            KEsum[l] = 0.0

            
        
            print("///////////////////////////////")
        KEd = [[0.0]*len(KEcom),[0.0]*((last-start)+1)]
        print str(len(KEd))
        
      # for c in range(0,len(KEd)):
       #    KEd[c] = [0.0]*((last-start)+1)



    for prnt in range(0,len(H)):
        print("Frag:" + str(H[prnt]) + "  Center of Mass KE " + str(KEcom[prnt]))
        
        #calculate KEsum 
        KEsum[prnt] = KEsum[prnt] + KEcom[prnt]
        
        
        
        KEd[prnt][KEdc] = KEcom[prnt]
    

    

    print "\n"
    
    KEtot = 0.0
    for s in range(0,len(KEcom)):
        KEtot = KEtot + KEcom[s]
    
    print("Sum of frag KE: " + str(KEtot))    
        #output.close()
    os.chdir('..')


    #increment KE data counter
    KEdc = KEdc + 1

    print("**************************************\n")

freqfrag = Counter(fragments)
print freqfrag
print "\n"

#print KE average for each fragement & Standard deviation

KEavg = []

#calculate KEavg
for z in range(0, len(KEsum)):
    KEavg.append([])

    KEavg[z] = KEsum[z]/((last-start)+1)
    
    print ("KEAVG  " + str(KEavg[z]))

print translate(atoms,H[0])
print translate(atoms,H[1])

print "\n \n"


sp = "    "

#create summary file write to file
f = file("./Summary.txt","w")
f.write("Summary\n")
f.write("  \n")


for frag in range(0, len(H)):
    
    f.write("Frag '"+ translate(atoms,H[frag]) + sp +"KE Average: " + str(KEavg[frag]) + sp + "KE Deviation: " + str(getSTD(KEd[frag],KEavg[frag])) + "\n")



os.popen("cat ./Summary.txt")












