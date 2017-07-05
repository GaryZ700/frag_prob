#list of classes used to analyze Turbomole simulations
import os
import json
import networkx as nx
import constants as const
import numpy as np

#stores molecule information and configuration of atoms in MD files
class Molecule:
    
    #generate molecular database folder and molecule json file
    def createMolecule(self,simulation,atoms,data):

        #check if analysis folder exists, if not, create it
        if(os.path.exists("./analysis") == False):
            os.popen("mkdir analysis")

        #move into analysis directory
        os.chdir("./analysis")

        #create simulation folder
        os.popen("mkdir " + simulation)

        #move into simulaion folder
        os.chdir("./" + simulation)

        #create json data file for molecule
        molecule = open("molecule.json",'w')
        json.dump(data,molecule)

        #return to parent directory
        os.chdir("../..")

    def read(self,simulation,timestep,**kargs):
        
        jsonData = open("./analysis/molecule.json",'r')
        data = json.load(jsonData)
        

        print("\n Atoms, Positions, Velocities")
            
        for atom in range(len(data["atoms"])):
           print(data["atoms"][atom])
           print(data[simulation][timestep]["position"][atom])
           print(data[simulation][timestep]["velocity"][atom])
            
        


#contains functions used in parsing data from MD log files
class IO:
    
    #return list of starting and ending simulations or mdlogs in a list
    def simulationRange(self,rangeType):
         
        Range = os.popen("ls " + rangeType +"* -d").read().split(rangeType)

        #remove blank element from start of list
        Range.pop(0)
        
        #get int list of simulation ending numbers
        Range = map(int,Range)
   
        return [min(Range),max(Range)]


    def allParseTimestep(self,atoms,mdFile):
        return 0

    #return list of atoms with positions and velocities at specified timesteps
    def parseTimestep(self,atoms,mdFile,timestep):
        
        data = []

        if(timestep == "last"):
            
            timesteps = os.popen("grep -n 't=' mdlog." + str(mdFile)).read().split("\n")
            timesteps.pop()
            #print(timesteps)

            #gets last timestep in md
            timesteps = timesteps[len(timesteps)-1].split(":")

            timestepLine = int(timesteps[0])
            timestepNumber = float(timesteps[1].split("=")[1])

            data.append({})
            data[0]["timestep"] = timestepNumber
            data[0]["position"] = []
            data[0]["velocity"] = []

        for atom in range((len(atoms)*2) + 1):
            
            #print(str(atom) + " atom number")

            coordinates = filter(None,os.popen("sed '"+ str((timestepLine + atom + 1)) + "!d' mdlog." + str(mdFile)).read().split(" ")) 
         
            #print(coordinates)

            if(atom < len(atoms)):
                coordinates.pop()
                data[0]["position"].append(map(float,coordinates))
            elif(atom != len(atoms)):
                coordinates.pop(0)
                data[0]["velocity"].append(map(float,coordinates))

        return data
        

    #return list of atoms in molecule    
    def parseAtoms(self,folderHead):

        #move into first simulation
        os.chdir(folderHead + str(self.simulationRange(folderHead)[0]))
            
        #open first mdlog
        mdFile = open("mdlog."+str(self.simulationRange("mdlog.")[0]),'r')

        #skip to structure data in file
        for i in range(3):
            next(mdFile)

        
        atoms = []

        for line in mdFile:
                                         
            #parse atom from line
            line = filter(None,line.split(" "))         
            atom = line[len(line)-1].split("\n")[0]

            #if atom in line, append to list
            if(len(atom) == 1):
                atoms.append(atom) 

            else: 
                os.chdir("..")
                return atoms
        

    #used to parse requested data from specifed structs and logs
    #parameters:
    #header = string, starting name of simulations
    #simulationsToParse = string or array, range of simulations to parse, "startInt-EndInt", "all", [list of simulations to parse]
    #MDsToParse = string or array, range of MD files to parse, "startInt-EndInt", "all", [list of MD logs to parse]
    #logsToParse = string or array, range of logs in simulation to parse, "startInt-EndInt", "all", [list of logs to parse]
    def parse(self, atoms, folderHead, timestepsToParse, simulationsToParse =  []):
        
        data = {}
        #for each simulation to parse
        for simulation in range(simulationsToParse[0] - 1, simulationsToParse[1]):
            
            #move simulation folder    
            os.chdir(folderHead+str(simulation + 1))

            print(os.popen("pwd").read())

            data[folderHead+str(simulation + 1)] = {} 

        
    
            MDsToParse = self.simulationRange("mdlog.")
                
            for md in range(MDsToParse[0],(MDsToParse[1] + 1)):
            
                print("mdlog."+str(md))
                data[folderHead+str(simulation + 1)].update(self.parseTimestep(atoms,md,timestepsToParse))             
                
            
            


            #return back to parent folder of simulations
            os.chdir("..")

        
        

        return data 
        




#contains functions and data structures used to calculate energies, and fragmentation of molecule
class Properties:

    def mass(self,molecule,atoms):
        mass = 0.0
        print("Molecule " + str(len(atoms)) )
        print(molecule)
        for atom in molecule:
            mass += const.atomMass[atoms[atom]]

        return mass

    def toAtoms(self,molecule,atoms):
        
        atomList = []
        for atom in molecule:
            atomList.append(atoms[atom])

        return atomList


    #construct moment of inertia tensor for a molecule, using center of mass position
    def momentOfInertia(self,molecule,atoms,newX):
        
        print("Newx")
        print(newX)
        I = [[0.0 for y in range(3)] for x in range(3)]

        for dimension in range(3):
            for atom in molecule:
                print(atom)

                I[dimension][dimension] += const.atomMass[atoms[atom]]*(newX[atom][(dimension+1)%3]**2 + newX[atom][(dimension+2)%3]**2)

        for dimension in range(3):
            for dimension2 in range(3):
                for atom in molecule:
                    if(dimension==dimension2):
                        continue
                    I[dimension][dimension2] += -const.atomMass[atoms[atom]]*(newX[atom][dimension]*newX[atom][dimension2])
        return I

    #return center of mass properties for molecule 
    def comXV(self,molecule,atoms,position,velocity):
    
        #init center of mass lists 
        Xcom = np.zeros(3)
        Vcom = np.zeros(3)
        print("XV")
        print(self.toAtoms(molecule,atoms))
        #calculate mass of molecule
        moleculeMass = self.mass(molecule,atoms)

        print(moleculeMass)
        print(velocity)
        print(position)

        #for each atom in molecule
        for atom in molecule:
            for dimension in range(3):
                Vcom[dimension] += (velocity[atom][dimension]*const.atomMass[atoms[atom]]/moleculeMass)
    
                Xcom[dimension] += (position[atom][dimension]*const.atomMass[atoms[atom]]/moleculeMass)
    
        newX = np.empty([len(position),3])
        newV = np.empty([len(velocity),3])
    
        for i in range(len(position)):
            newX[i] = position[i] - Xcom
            newV[i] = velocity[i] - Vcom   
    
        return newV,newX 


    #find rotational kinetic energy of atom group
    def RKE(self,molecule,atoms,position,velocity):
        
        Xcom, Vcom = self.comXV(molecule,atoms,position,velocity)

        I = self.momentOfInertia(molecule,atoms,Xcom)


        eigenValues, eigenVectors = np.linalg.eig(I)

        L = np.zeros(3)

        for atom in range(len(molecule)):
            L += const.atomMass[atoms[molecule[atom]]] * np.cross(Xcom[atom],Vcom[atom])

        rke = 0.0
        for dimension in range(3):
            rke += (np.dot(L,eigenVectors[dimension]))**2 / (2*eigenValues[dimension])

        print(rke)

        return rke

    




    #find cartesian distances between two points in 3D space
    def distance(self,position1,position2):
        return ((position1[0]-position2[0])**2 + (position1[1]-position2[1])**2 + (position1[2]-position2[2])**2)**0.5
        
    #find occurance of fragmentation in molecule
    def findFragments(self,atoms,position,cutoff):
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


