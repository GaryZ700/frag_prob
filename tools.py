#list of classes used to analyze turbomole simulations


#stores molecule information and configuration of atoms in MD files
class molecule:

#contains functions used in parsing data from MD log files
class IO:

    #used to parse requested data from specifed structs and logs
    #parameters:
    #header = string, starting name of simulations
    #simulationsToParse = string or array, range of simulations to parse, "startInt-EndInt", "all", [list of simulations to parse]
    #MDsToParse = string or array, range of MD files to parse, "startInt-EndInt", "all", [list of MD logs to parse]
    #logsToParse = string or array, range of logs in simulation to parse, "startInt-EndInt", "all", [list of logs to parse]
    def parse(header, simulationsToParse, MDstoParse, logsToParse):
        :



#contains functions and data structures used to calculate energies, and fragmentation of molecule
class properties:

