#list of classes used to analyze turbomole simulations

import os


#stores molecule information and configuration of atoms in MD files
class molecule:
    h = 0

#contains functions used in parsing data from MD log files
class IO:

    
    
    #return list of starting simulation, and ending simulation
    def simulationRange(self,rangeType):
         

        Range = os.popen("ls " + rangeType +"* -d").read().split(rangeType)

        #remove blank element from start of list
        Range.pop(0)
        
        #get int list of simulation ending numbers
        Range = map(int,Range)
        

        return [min(Range),max(Range)]

    #return list of atoms with positions and velocities at timestep
    def parseTimestep(self,mdFile,timestep):
        
        if(timestep[0] == "last"):
            
            timesteps = os.popen("grep -n 't=' mdlog." + mdNum).read().split("\n")
            lineNum = int(timesteps[len(timesteps)-1])
            



        

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
    def parse(self, atoms, folderHead, simulationsToParse =  [], MDsToParse = [], timestepsToParse = []):
        

        #for each simulation to be parse
        for simulation in range(simulationsToParse[0],simulationsToParse[1]):
            #move simulation folder    
            os.chdir(folderHead+str(simulation + 1))

            print(os.popen("pwd").read())
            
            
            if(MDsToParse[0] == "all"):
                MDsToParse = self.simulationRange("mdlog.")
                
            for md in range(MDsToParse[0],(MDsToParse[1] + 1)):
                parseTimestep(md,timeStepstoParse)             
                
                


            #return back to parent folder of simulations
            os.chdir("..")

        
        

    
        




#contains functions and data structures used to calculate energies, and fragmentation of molecule
class properties:
    h = 0
