#hidden documentation file for Turbomole anlysis scripts


#each class contains documentation for appropriate anlysis script

class analysis:

   overview = "Python analysis script for Turbomole Simulations. Script should be run in folder containing simulation folders, not inside of the simulationis themselves. This script will generate an anlysis folder containg molecular data for performing further analysis, along with analysis data. Analysis data can be viewed using the -va option." 

   f = "Simulation folder header name. Default name is set to b4500m."

   ss = "Number of simulation to begin parsing from. If no ending simulation is specified, parsing wil continue until last simulation is reached."

   es = "Number of last simulation that should be parsed. If no starting simulation is specified, parsing will begin at first simulation."

   sm = "MD log where data parsing should begin. If ending log is not specified, parsing will end at last MD log."

   em = "MD log where data parsing should end. If starting log is not specified, parsing will begin from first MD log."

   st = "Starting data structure in MD log. If ending structure is not specified, parsing will end at last structure."

   et = "Ending data structure in MD log. If starting structure is not specified, parsing will begin at first structure. Set et and st to same log if only one structure should be parsed."

   ip = "Individual parsing input, Pass list of lists in form [[s,m,t],[s1,m2,t2...[sn,mn,tn]], where s is the simulation, m is the MD log, and t is the data structure. s, m, and t can also be lists if parsing of more than one item is needed."
   
   ts = "Timestep Parsing Option. Pass a string of 'all', 'first', or 'last'. Will parse respective timesteps from Md logs. Default option is 'last.'"
