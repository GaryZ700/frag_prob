#stores molecule information and configuration of atoms in MD files

import os
import json

class Molecule:
   
########################################################################
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
############################################################################################
    #read information from molecule json
    def read(self,simulation,timestep,**kargs):
        
        jsonData = open("./analysis/molecule.json",'r')
        data = json.load(jsonData)
        

        print("\n Atoms, Positions, Velocities")
            
        for atom in range(len(data["atoms"])):
           print(data["atoms"][atom])
           print(data[simulation][timestep]["position"][atom])
           print(data[simulation][timestep]["velocity"][atom])
            
