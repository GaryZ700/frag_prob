#Script to perform in depth analysis on Turbomole Parser script data

import os
import argparse as ap
from IO import IO
from Molecule import Molecule

molecule = Molecule()
io = IO()

parser = ap.ArgumentParser()

parser.add_argument('-f',metavar="folderHead",default="b4500m",type=str)

args = parser.parse_args()

####################################################################################

def generateReport(folderHead):
    #generate report on md simulations

    #move into analysis folder
    os.chdir("./analysis")

    simulations = io.simulationRange(folderHead)

    report = open("report",'w')

    report.write("Turbomole Analysis Report")
    
    report.write("\n")

    report.write("Simulation: " + folderHead + str(simulations[0]) + " to " + folderHead + str(simulations[1]))
    
    report.write(fragmentPathways(folderHead,simulations))

    report.close()


    #print report on screen
    report = open("report",'r')
    print(report.read())

    #return to parent directory
    os.chdir("..")


###############################################################################
def fragmentPathways(folderHead,simulations):
    
    fragments = []

    for simulation in range(simulations[0],simulations[1] + 1):
        
        timesteps = molecule.load(folderHead + str(simulation))

        for timestep in timesteps:
            fragments += timestep["molecule"]["moleculeList"]
    
    return Counter(fragments)


###############################################################################
#main code goes here


#check if analysis folder exists
if(os.path.exists("./analysis")):
    generateReport(args.f)
else:
    print("Can not perform analysis, analysis folder does not exist in current directory. Try to run mdparser.py before running analysis.py.")
    
    #data = molecule.load(str(folderHead + IO.simulationRange(folderHead)[1]))



