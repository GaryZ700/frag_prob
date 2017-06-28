#script to run analysis on turbomole simulations

from tools.py import IO
import argparse as ap

#parse script arguemnts################################################################################################################################################################

#-h: help, provide documetation on script

#-ss: int, starting simulation, sumulation to start from, defualt = 0
#-es: int, ending simulation, simulation to end from, defualt = max
#-sm: int, starting md log, md log to start from, defualt = 0
#-em: int, ending md log, md log to end form, defualt = max
#-st: int, starting structure, structure to start from in md log, defualt - 0
#-et: int, ending structure, structure to end on in md log, defualt = max

#-ip: list of lists, individual parsing, list of type [[s,m,t],[s1,m1,t1],...[sn,mn,tn]], where s = simulation int, m = md log int, t = struct int
#-ip: continued, s,m, and t can also be lists in themselves, ex. [ [ [1,2], [4,5,6], [1,2,3] ] ], would parse the first three data structures from the 4th, 5th, and 6th md logs of the first two simulations
#-ip: continued, defualt = not used

#-nf: t/f, new folder, save anlysis data to new folder, defualt = true
#-o: string of anlysis folder, overwrite, save anlysis data to existing data file, defualt = false

#-fn: string, folder name, name of anylysis folder, defualt = anlysisN, where n is the number of previous analyzations that have occured

parser = ap.ArgumentParser()



#initialize main IO for script
io = IO()
