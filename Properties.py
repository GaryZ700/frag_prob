import numpy as np
import networkx as nx
import constants as const

#contains functions and data structures used to calculate energies, and fragmentation of molecule
class Properties:

###################################################################
    #find mass of molecule
    def mass(self,molecule,atoms):
       
        #init mass 
        mass = 0.0
        
        #sum mass of each atom in molecule
        for atom in molecule:
            mass += const.atomMass[atoms[atom]]

        return mass
###################################################################
    #turn list of atom numbers into atom symbol
    def toAtoms(self,molecule,atoms):
        
        #init list of atom symbols
        atomList = []

        #for each atom number in the molecule, 
        #append an atom character to the atom list
        for atom in molecule:
            atomList.append(atoms[atom])

        return atomList
###################################################################
    #construct moment of inertia tensor for a molecule about the center of mass
    def momentOfInertia(self,molecule,atoms,newX):
       
        #generate list of three empty lists
        I = [[0.0 for y in range(3)] for x in range(3)]
        
        #for each dimension and for each atom in molecule
        #construction of diagonals of moment of inertia tensor
        for dimension in range(3):
            for atom in molecule:
                I[dimension][dimension] += const.atomMass[atoms[atom]]*(newX[atom][(dimension+1)%3]**2 + newX[atom][(dimension+2)%3]**2)

        #calculate off diagnoals of tensor by ....
        for dimension in range(3):
            for dimension2 in range(3):
                for atom in molecule:
                    if(dimension==dimension2):
                        continue
                    I[dimension][dimension2] += -const.atomMass[atoms[atom]]*(newX[atom][dimension]*newX[atom][dimension2])
        
        return I
#################################################################
    #return center of mass properties for molecule 
    def comXV(self,molecule,atoms,position,velocity):
    
        #init center of mass lists 
        Xcom = np.zeros(3)
        Vcom = np.zeros(3)

        #calculate mass of molecule
        moleculeMass = self.mass(molecule,atoms)

        #for each atom in molecule,
        #and for each dimension
        #calculate center of mass velocity and position for each dimension
        for atom in molecule:
            for dimension in range(3):
                
                Vcom[dimension] += (velocity[atom][dimension]*const.atomMass[atoms[atom]]/moleculeMass)
                Xcom[dimension] += (position[atom][dimension]*const.atomMass[atoms[atom]]/moleculeMass)
    
       
        newX = np.empty([len(position),3])
        newV = np.empty([len(velocity),3])
    
        for i in range(len(position)):
            newX[i] = position[i] - Xcom
            newV[i] = velocity[i] - Vcom   
    
        return newX, newV 
#################################################################
    #find rotational kinetic energy of atom group
    def RKE(self,molecule,atoms,position,velocity):
       
        #get center of mass position and velocity for molecule
        Xcom, Vcom = self.comXV(molecule,atoms,position,velocity)

        #calculate moment of inertia tensor using center of mass position
        I = self.momentOfInertia(molecule,atoms,Xcom)

        #find eigen values and vectors of moment of Inertia
        eigenValues, eigenVectors = np.linalg.eig(I)

        #init list for angular momentum
        L = np.zeros(3)

        #calculate angular momentum for each atom in molecule
        for atom in range(len(molecule)):
            L += const.atomMass[atoms[molecule[atom]]] * np.cross(Xcom[atom],Vcom[atom])

        #init rke
        rke = 0.0

        #for each dimension, sum dot product of angular momentum and egien vector,
        #divided by twice the eigen values
        for dimension in range(3):
            rke += (np.dot(L,eigenVectors[dimension]))**2 / (2*eigenValues[dimension])

        return rke
######################################################################
    #find cartesian distances between two points in 3D space
    def distance(self,position1,position2):
        return ((position1[0]-position2[0])**2 + (position1[1]-position2[1])**2 + (position1[2]-position2[2])**2)**0.5
######################################################################        
    #get structure of molecule
    def findFragments(self,atoms,position,cutoff):
        
        #
        atomGraph = nx.Graph()

        #generate a graph with one point per atom        
        for atom in range(len(atoms)):
            atomGraph.add_node(atom)

        print("Position")
        print(position)
        #picks an atom, and a second atom,
        #and then compares the distance between the two atoms,
        #if the distance is less than the cutoff, "bond" the atoms
        for atom in range(len(atoms)):
            for nextAtom in range(atom + 1, len(atoms)):
                if(self.distance(position[atom],position[nextAtom]) < cutoff):
                    atomGraph.add_edge(atom,nextAtom)
        
        #return list of each bonded atom group, or fragment
        fragments =  [list(atoms) for atoms in nx.connected_components(atomGraph)]

         
        return fragments


