#script to generate molecular data for analysis.py

from tools import IO
from tools import molecule
import argparse as ap
from documentation import analysis as doc

#parse script arguemnts##################################################################

#-h: help, provide documetation on script

#-ss: int, starting simulation, sumulation to start from, defualt = 0
#-es: int, ending simulation, simulation to end from, defualt = max
#-sm: int, starting md log, md log to start from, defualt = 0
#-em: int, ending md log, md log to end form, defualt = max
#-st: int, all, first, last
#-et: int, 

#-ip: list of lists, individual parsing, list of type [[s,m,t],[s1,m1,t1],...[sn,mn,tn]], where s = simulation int, m = md log int, t = struct int
#-ip: continued, s,m, and t can also be lists in themselves, ex. [ [ [1,2], [4,5,6], [1,2,3] ] ], would parse the first three data structures from the 4th, 5th, and 6th md logs of the first two simulations
#-ip: continued, default = not used

#-nf: t/f, new folder, save anlysis data to new folder, defualt = true
#-o: string of anlysis folder, overwrite, save anlysis data to existing data file, defualt = false

#-fn: string, folder name, name of anylysis folder, defualt = anlysisN, where n is the number of previous analyzations that have occured

io = IO()
parser = ap.ArgumentParser(description = doc.overview)

parser.add_argument('-ss',metavar="startSimulation",default="all", help=doc.ss)
parser.add_argument('-es',metavar="endSimulation",default="all", help=doc.es)
parser.add_argument('-sm',metavar="startMD",default="all", help=doc.sm)
parser.add_argument('-em',metavar="endMD",default="all", help=doc.em)
parser.add_argument('-st',metavar="startStructure",default="last", help=doc.st)
parser.add_argument('-et',metavar="endStructure",default="last", help=doc.et)
parser.add_argument('-ip',metavar="indvidualParsing",default="[void]", help=doc.ip)
parser.add_argument('-f',metavar ="folderHead",default="b4500m",help=doc.f)
#parser.add_argument('-nf',metavar="nf", help=doc.nf)
#parser.add_argument('-ss',metavar="ss",default="0", help=doc.ss)
#parser.add_argument('-ss',metavar="ss",default="0", help=doc.ss)

args = parser.parse_args()

print(args)

simulationRange = io.simulationRange(args.f)

if(args.ss == "all"):
    args.ss = simulationRange[0]
if(args.es == "all"):
    args.es = simulationRange[1]




atoms = io.parseAtoms(args.f)



io.parse(atoms,args.f,[args.ss,args.es],[args.sm,args.em],[args.st,args.et])


