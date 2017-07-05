#master anlysis script for Turbomole Simulations

from IO import IO
from Molecule import Molecule
from Properties import Properties
import argparse as ap
from documentation import analysis as doc
import os
import constants as consts

#Script arguemnts##################################################################

# By default, script will parse all simulations, and mdlogs, and last timestep in each mdlog
# will perform full fragment analysis
# and will present results to user

# -h: help, provide documetation on script

# -ss: int, starting simulation, sumulation to start from, defualt = 0
# -es: int, ending simulation, simulation to end from, defualt = max
# -sm: int, starting md log, md log to start from, defualt = 0
# -em: int, ending md log, md log to end form, defualt = max
# -st: int, all, first, last
# -et: int, 

# -ip: list of lists, individual parsing, list of type [[s,m,t],[s1,m1,t1],...[sn,mn,tn]], where s = simulation int, m = md log int, t = struct int
# -ip: continued, s,m, and t can also be lists in themselves, ex. [ [ [1,2], [4,5,6], [1,2,3] ] ], would parse the first three data structures from the 4th, 5th, and 6th md logs of the first two simulations
# -ip: continued, default = not used

# -nf: t/f, new folder, save anlysis data to new folder, defualt = true
# -o: string of anlysis folder, overwrite, save anlysis data to existing data file, defualt = false

# -fn: string, folder name, name of anylysis folder, defualt = anlysisN, where n is the number of previous analyzations that have occured


# init global variables################################################################################

# main tool classes
io = IO()
molecule = Molecule()


#init script arguments
parser = ap.ArgumentParser(description = doc.overview)

parser.add_argument('-ss',metavar="startSimulation", type=int,help=doc.ss)
parser.add_argument('-es',metavar="endSimulation", type=int,help=doc.es)
#parser.add_argument('-sm',metavar="startMD",type=int, help=doc.sm)
#parser.add_argument('-em',metavar="endMD", type=int, help=doc.em)
parser.add_argument('-ts',metavar="timestep",default="last",type=str, help=doc.ts)
#parser.add_argument('-st',metavar="startStructure", help=doc.st)
#parser.add_argument('-et',metavar="endStructure",default="last", help=doc.et)
#parser.add_argument('-ip',metavar="indvidualParsing",default="[void]", help=doc.ip)
parser.add_argument('-f',metavar ="folderHead",default="b4500m", type=str, help=doc.f)
#parser.add_argument('-nf',metavar="nf", help=doc.nf)
#parser.add_argument('-ss',metavar="ss",default="0", help=doc.ss)
#parser.add_argument('-ss',metavar="ss",default="0", help=doc.ss)

args = parser.parse_args()

print(args.ss)

#init global functions###############################################################################
#create analysis folder and populate with specified MD log data
def firstRun(args):

    simulationsToParse = []
    MDsToParse = []
    timestepsToParse = []
    

    simulationRange = io.simulationRange(args.f)

    if(args.ss is None):
        simulationsToParse.append(simulationRange[0])
    else:
        simulationsToParse.append(args.ss)
    if(args.es is None):
        simulationsToParse.append(simulationRange[1])
    else:
        simulationsToParse.append(args.es)

    atoms = io.parseAtoms(args.f)
    generateAnalysis(atoms,args.f,args.ts,simulationsToParse)



#calculate specified molecular properties
def calculateProperties(timestepData,atoms):

    properties = Properties()
    
    timestepData["fragments"] = {}
    timestepData["fragments"]["fragmentList"] = properties.findFragments(atoms, timestepData["position"],consts.cutoff)

    if(timestepData["fragments"]["fragmentList"] != False):
      
        timestepData["fragments"]["RKE"] = []
        timestepData["fragments"]["KE"] = []

        for fragment in timestepData["fragments"]["fragmentList"]:
    
            timestepData["fragments"]["RKE"].append(properties.RKE(fragment,atoms,timestepData["position"],timestepData["velocity"]))
       
       #timestepData["fragments"]["fragmentList"]["KE"] = properties.KE(timestepData["fragments"]["fragmentList"]) 



    return timestepData

#present specified data to user
def viewResults():
    return 0

#parses requested data from MD logs
#and returns list of 
def generateAnalysis(atoms, folderHead, timestepsToParse, simulationsToParse =  []):
    
    #for each simulation to parse
    for simulation in range(simulationsToParse[0] - 1, simulationsToParse[1]):
        
        #init data list
        data = []

        #move simulation folder  
        os.chdir(folderHead+str(simulation + 1))

        print(os.popen("pwd").read())

        #find range of mdlogs in simulation
        MDsToParse = io.simulationRange("mdlog.")
        
        #for ecah mdlog in simulation, parse specified struct
        #and append timestep to data list
        for md in range(MDsToParse[0],(MDsToParse[1] + 1)):
        
            print("mdlog."+str(md))
            data += io.parseTimestep(atoms,md,timestepsToParse) 
        
        for timestepData in data:
            timestepData = calculateProperties(timestepData,atoms)
            

        #return back to parent folder of simulations
        os.chdir("..")
        

        molecule.createMolecule((folderHead + str(simulation)),atoms,data)


 


#main code goes here#####################################################################################


#check if analysis folder already exists 
if(os.path.exists("./analysis") == False):
    #if not, create analysis folder and populate with simulation data
    firstRun(args)

