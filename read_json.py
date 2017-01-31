import json
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats




f = open("data.json","r")

data_f = json.load(f)
n = len(data_f)
pathways = []
fragments = []
KE = []

KE.append("KEcom")
KE.append("RKE")
KE.append("TKE")
for i in range(n):
    for j in range(i,n):
        if (str(data_f[i]["fragments"]) > str(data_f[j]["fragments"])):
            temp  = data_f[i]
            data_f[i] = data_f[j]
            data_f[j] = temp
    pathway = data_f[i]["fragments"]
    if str(pathway) not in pathways:
        pathways.append(str(pathway))
    for frag in pathway:
        if frag not in fragments:
            fragments.append(frag)

exit = False

listofI = []

while(exit == False):    


    print pathways

    print("""
    Choose fragment:
    """)
    for frags in fragments:
        print fragments.index(frags), frags
     
    print(str(len(fragments))+" Exit")
    chfrag = int(raw_input("Enter number"))
    
    if(chfrag == len(fragments)):
       break

    print("""
    Choose pathway:
    """)
    for path in pathways:
        print pathways.index(path), path
    chpath = int(raw_input("Enter number"))

    print("""
    Choose KE:
    """)
   
    for ke in KE:
         print KE.index(ke), ke
    chke = int(raw_input("Enter number"))

    

    for i in range(n):
        if(str(data_f[i]["fragments"])==pathways[chpath]):
            nl = data_f[i]["fragments"].index(fragments[chfrag])

            listofI.append(data_f[i][KE[chke]][nl])

    print("Would you like to enter another fragment? y/n")
    yn = str(raw_input())
    print("\n")
    if(yn == "y"):
	continue    

    print listofI
   
 
    
    print(np.std(listofI))
    mu = np.mean(listofI)
    stdev = np.std(listofI)
    # gauss = np.random.normal(mu, stdev,size = len(listofI))
    plotlist = []
    xlength = np.linspace(min(listofI),max(listofI),20)
    for x in xlength:
        plotlist.append(scipy.stats.norm(mu,stdev).pdf(x))
#    print(str(gauss,))
    print("\n")
    #print("Standard Deviation: "+str(stdev)+"    Mu: "+str(mu))
    raw_input("Hit Enter to View Graph.")
    
   # plt.hist(gauss, bins = 10, histtype='stepfilled', normed=True, color='b', label='Hist')
 
    plt.hist(listofI,bins = 10,histtype='stepfilled', normed=True, color='r', alpha=0.5, label='Uniform')
    plt.scatter(xlength,plotlist)
    plt.title("Gaussian/Uniform Histogram")
    plt.xlabel("Value")
    plt.ylabel("Probability")
    plt.legend()
    plt.show()
    

    listofI = []
    

 
f2 = open("temp","w")
json.dump(data_f,f2)
f2.close()
