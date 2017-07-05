#contains functions used in parsing data from MD log files

import os


class IO:
    
###############################################################################    
    #return list of starting and ending simulations or mdlogs in a list
    def simulationRange(self,rangeType):
        
        #get list of simulations or mdlogs in directory
        Range = os.popen("ls " + rangeType +"* -d").read().split(rangeType)

        #remove blank element from start of list
        Range.pop(0)
        
        #get int list of simulation ending numbers
        Range = map(int,Range)

        #return minumum and maximum value of simulation/md log
        return [min(Range),max(Range)]
###############################################################################    
    #parses all timesSteps from md logs
    def allParseTimestep(self,atoms,mdFile):
        return 0
###############################################################################
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
##############################################################################
    #return list of atoms in simulation    
    def parseAtoms(self,folderHead):

        #move into first simulation
        os.chdir(folderHead + str(self.simulationRange(folderHead)[0]))
            
        #open first mdlog
        mdFile = open("mdlog."+str(self.simulationRange("mdlog.")[0]),'r')

        #skip to structure data in file
        for i in range(3):
            next(mdFile)

        #init atoms list
        atoms = []
        

        for line in mdFile:
                                         
            #parse atom from file line
            line = filter(None,line.split(" "))         
            atom = line[len(line)-1].split("\n")[0]

            #if atom in line, append to list
            if(len(atom) == 1):
                atoms.append(atom) 
            
            #if atom not in line, return to parent directory and return list of atoms
            else: 
                os.chdir("..")
                return atoms
######################################################################################

        


