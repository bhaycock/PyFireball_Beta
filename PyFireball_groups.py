from PyFireball import *

""" PyFireball needs to be split up for readability. I also need a full new set of tools to take care of a study
of ~ 1000 or more supercells. I'm calling the group of, e.g. a single percentage, um... I need a name. I need the 
ability to sort by energy, I need to produce simple things like the clustering and the varience output. I need to pop
out the DOS calcs, and to automate a full summary of the directory, etc.

Required functions: 
Make list of directories required.
Find a way to dump off failed runs.
DOS of LE structs with average and all.
Varience  X
boxplot varience.


"""
def listRunsInDirectory(superDirectory):
  import os
  """ Lists all the directories in a directory so as I can run the whole lot into memory"""
  # Create a "listOfRuns" file.
  runsFileLocation = superDirectory + '/listOfRuns'
  command = 'ls -p ' + superDirectory + ' | grep "/" > ' + runsFileLocation
  os.system(command)
  
  # Read in the listOfRuns file and remove Fdata directory... if it's there!
  runsFile = open(runsFileLocation, 'r')
  
  filteredRunFile =[]
  for runDirectory in runsFile.splitlines():
    if runDirectory != 'Fdata/' or runDirectory !='Fdata':
      filteredRunFile.append(runDirectory)  
  runsFile.close()

  command = 'rm ' + runsFileLocation
  os.system(command)
  
  runsFile = open(runsFileLocation, 'w')
  for runDirectory in filteredRunFile:
    runsFile.write(runDirectory + "/n")

  runsFile.close()
  return filteredRunFile
  
def readBunchOfDirectories(superDirectory):
  listOfDirectories = listRunsInDirectory(superDirectory)
  listOfSupercells = []
  for directory in listOfDirectories:
    listOfSupercells.append(directory2SupercellLiteVersion(directory))
  return listOfSupercells
  
def sortSupercellList(listOfSupercells, key = 'etot'):
  """ Need a routine to grab the supercell's CohE, ETOT and CohE per atom, methinks. Damn."""
  if key.lower() == 'etot':
    return (sorted(listOfSupercells, key = lambda cell: cell.etot))
  else if key.lower() == 'cohe':
    return (sorted(listOfSupercells, key = lambda cell: cell.cohE))
  else:
    print "Invalid key in sort Supercell List"
    return None
    
def generateAndPlotClustering(supercellList, speciesOfInterest, rinput, elements, label = "", etype = "etot", relative = True, normalized = True):
  """ idea here is to make a clustering plot for a single list of supercells. Type of energy is choosable, and we can build here to make a list-of-lists for the final push on PyFireball. Plots similar to this can be found in my most recent paper."""
  #Generate a clustering / energy list.
  (energy, clustering) = generateClustering(supercellList, speciesOfInterest, rinput, elements, etype = "etot", relative, normalized)
  #Send to the plotter
  plotClustering(energy, clustering, label, relative, normalized)  
      
def generateClustering(supercellList, speciesOfInterest, rinput, elements, etype = "etot", relative = True, normalized = True):
  #Generate a clustering / energy list.
  energy = []
  clustering = []
  for icell in supercellList:
    if etype.lower() == "etot" : energy.append (icell.etot)
    if etype.lower() == "cohe" : energy.append (icell.cohE)
    clustering.append(icell.clustering(icell, speciesOfInterest, rinput, elements))
  energy = np.asnumeric(energy)
  if relative : energy = energy - min(energy)
  if normalized : energy = energy / max(energy) 
  return (energy, clustering)
  
def plotClustering(energy, clustering, label, relative, normalized):
  # Plots the clustering.
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  if inLabel != "" : ax.scatter(CuInFe_Clustering, CuInFe_NormRelEng, color = 'blue', marker = '*', s=90)
  else: ax.scatter(CuInFe_Clustering, CuInFe_NormRelEng, color = 'blue', marker = '*', s=90, label = inLabel)

  ax.set_xlabel("Clustering Factor ($\AA$)", fontsize = 18)

  if relative and normalized: ax.set_ylabel("Normalized Relative Energy (eV)", fontsize = 18)
  elif relative and not normalized: ax.set_ylabel("Relative Energy (eV)", fontsize = 18)
  elif normalized and not relative: ax.set_ylabel("Normalized Energy (eV)", fontsize = 18)
  else: ax.set_ylabel("Energy (eV)", fontsize = 18)

  ax.axis('tight')

  ax.grid(True)

  #Stuff to marker below from http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
  #Should put legend above the box, the URL gives below the box
  # Shink current axis's height by 10%  
  box = ax.get_position()
  ax.set_position([box.x0, box.y0 + box.height * 0.0, box.width, box.height * 0.9])
# Put a legend above current axis
  ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.14), fancybox=True, shadow=True, ncol=3)
# End of Stuff from URL


  fig.show()
  fig.savefig("ClusterPlot.png")
  return
  
def calculateVarience(supercellList, baseSpecies, seekSpecies, numSeek, etype = "etot", relative = True, normalized = True):
  """ Similar to calculation of clustering. So that's the  first thing.
  according to http://stackoverflow.com/questions/19391149/numpy-mean-and-variance-from-single-function, this is kinda easy... and
  allows for my box-plot implimentation.
  What I want is a routine that, per supercell, returns the list of the lengths."""
  variencesArrayList = []
  variencesSummaryList = []
  energy = []
  for cell in supercellList:
    if etype.lower() == "etot" : energy.append (icell.etot)
    if etype.lower() == "cohe" : energy.append (icell.cohE)
    varArray = PF.varience_array(cell, baseSpecies, seekSpecies, numSeek)
    variencesArrayList.append(varArray)
    varienceSummaryList.append((varArray.mean(), varArray.var())
  
  energy = np.asnumeric(energy)
  if relative : energy = energy - min(energy)
  if normalized : energy = energy / max(energy) 
  return (energy, variencesArrayList, variencesSummaryList)
  
def plotSummaryVarience(energy, varienceSummaryList, inLabel ="", relative, normalized, xlim = [], ylim = [], equilibrium = 0.00)
# Begin plot object
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  # Unpack the summary varience data.
  distList =[]
  varienceList = []
  for ipoint in varienceSummaryList:
    (dist, varience) = ipoint
    distList.append(dist - equilibrium)
    varienceList.append(varience)
  
  #Add this data to the plot:
  if inLabel == "" : ax.errorbar(distList, energy, xerr=varienceList, color = 'red', marker = '^', fmt = ' ')
  else: ax.errorbar(distList, energy, xerr=varienceList, color = 'red', marker = '^', fmt = ' ', label = inLabel)

  #Assemble plot
  if xlim: ax.set_xlim(xlim[0],xlim[1])
  if ylim: ax.set_ylim(ylim[0], ylim[1])
  if not xlim and not ylim: ax.axis('tight')
  if equalibrium != 0.00 : ax.set_xlabel("Average Distance from Equilibrium ($10^{-3} \AA$)", fontsize = 18)
  else ax.set_xlabel("Average Distance between Sites ($10^{-3} \AA$)", fontsize = 18)

  if relative and normalized: ax.set_ylabel("Normalized Relative Energy (eV)", fontsize = 18)
  elif relative and not normalized: ax.set_ylabel("Relative Energy (eV)", fontsize = 18)
  elif normalized and not relative: ax.set_ylabel("Normalized Energy (eV)", fontsize = 18)
  else: ax.set_ylabel("Energy (eV)", fontsize = 18)

 # Vertical line if we have offset for equilibrium
  if equilibrium != 0.00 : ax.plot([0, 0],[min(energy) - 0.05,max(energy)] , color="black", lw=2, linestyle='dashed') 

  ax.grid(True)
  fig.show()
  
  #Stuff to marker below from http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
  #Should put legend above the box, the URL gives below the box
  # Shink current axis's height by 10%  
  box = ax.get_position()
  ax.set_position([box.x0, box.y0 + box.height * 0.0, box.width, box.height * 0.9])

# Put a legend above current axis
  ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22), fancybox=True, shadow=True, ncol=3)
# End of Stuff from URL

  fig.savefig("variancePlot.png")  
  return

  
  
  
  
    