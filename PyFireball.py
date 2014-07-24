import sys
import math
import copy
import matplotlib.pyplot as plt 
import numpy as np

# Add clustering factor and also a clustering in terms of nearest neighbors.
# Want to add ssh capabilities to upload / download files AND to run jobs on (e.g.) MOUNTAINEER
# Need an "Output DOScalc" method.
# Need an "Output BandPlotCalc" method
# Need "Output to XcrysDen (with script)"
# Need boxPlot varience (after sorting multi-inputs)


""" Collection of small but useful functions to assist in accessing FIREBALL data using Python
This set currently contains the following routines:

ReadBas(basfile) - takes the .bas filename and returns a List of lists of atom(Symbol, X, Y, Z), indexed over iatom.

basList2Object(atomList) - TEMPORARY- Accepts the list from ReadBas, returns a list of atom objects.

chargesList2atoms(atomList, chargesList) - TEMPORARY - Accepts the list from basList2Object and the List from readCharges and puts charges into the shells in each atom object.

readDOS(atomList, directory = '.') - reads in the DOS information and appends to the atomList list of atom objects. One can therefore look at the DOS of any group of atoms and shells.

density(atom) - Just returns the terminal the dens_ file for that atom.

computeZtot(atomList) - calculates Ztot (total no. of electrons) in atomList -TODO- Make this a supercell class method.

plotDOS(atomList)(, fermiEnergy, emin = fermiEnergy -2, emax = fermiEnergy + 4, relativeEnergy = TRUE, byShell = FALSE) - NEEDS IMPLEMENTATION - Returns a plotDOS object, which contains methods to plot to screen, output to file, plot to file, etc.

neighList2atomObjects(atomList, neighList) - TEMPORARY - Accepts the list of atom objects and the neighList from readneigh, appends a list of neighbors to each atom object, where the list of neighbors is the atom object of that atom. Also adds the xlm vectors to the atom object so as mbeta can be automagically added to returned neighbours using the .neighbor method
							
symbol(nZ) - accepts an integer, the atomic number and returns the atomic symbol.

ReadLVS(lvsFile) - takes the .lvs filename and returns a List of Lists of LVS(X, Y, Z), indexed over dimension.
			N.B.- use numpys asnumeric() to convert list of lists into Arrays, if needed. 

ReadNeigh(directory = '.') - reads in the NEIGHBORS file from a TG directory (Lightning implementation in future work), returns NeighList #Need to specify working directory for this

ReadCharges(directory = '.') - Returns the CHARGES file in a list-of-lists format, indexed over iatom.
	
readFdata(fdataLocation) PyFIREBALL implementation of the readinfo.f90. Reads in Fdata info and creates species objects. Returns a directory of species indexed by the element names.

NeighborCount(Nlist, basList, iatom, neighspecies, depth, ignorelist = []) - Takes in a Neighbor List of the form output from ReadNeigh and a basis list of the form output by ReadBas, along with the atom concernd and the species we want and depth, returns the number of neeighbors of that depth. This is recursive: the ignorlist means any already-counted atoms are ignored.

NearestSpecList(Nlist, iatom, basList, neighSpecies, numneigh) - IN PROGRESS - Returns an ordered list (distance low-high) of neighbors. This method again uses the neighbor map, sweeps through to find the neighbors that are of species neighSpecies, works out their distance from the iatom, then sorts that list and returns the requested numneigh. Unlike Neighborcount, it does not look to neighbors-of-neighbors

atomDistance(List, List) - returns float of distance between two atoms. The List representing the two atoms is in the form (AtomicNo(int), X(float), Y(float), Z(float).

xlmbeta(LVS) - rerums the xl vectors for cell images.

NNList(atomArray, LVSArray, atom#, depth) - NOT YET IMPLEMENTED- returns the nearestneighbor list for atom, depth gives you how many NNs you want (1st NN, 2nd NN, etc) If you want ONLY the 2nd NN, then run a list comprenesion, "for b not in a" where a is 1st NN List and b is 2nd NN List.
		
NeighborCount(Nlist, iatom, basList, neighspecies, depth) - retruns the number of neighbors of the identified species. Depth to Nth nearest neighbors.

NearestSpecList(Nlist, iatom, basList, neighSpecies, numneigh) - IN PROGRESS -returns nearest neighbor list, sorted from lowest to highest. numneigh is the number of returned neighbors.

plotDOSList(densList, name = "multiDOSplot", legend = True, plotContributions = True, plotHOMO = True, plotOffset = True) - Returns a 3x1 plot of 3 and only 3 DOS objects. I want to add a limit to this, but first I gotta sort offset.

catDOSTots(DOSdict, passHOMO = True):  # Accepts a dictionary of DOS objects and returns a DOS object with the Totals of each dictionary object as a contribution.

These are part of a short project to calculate and present dipoles within the Delafossite system, but can be applied to anything.

There are TG versions and LIGHTNING versions for each of the read sections.

This module also includes implementations of the New-FIREBALL-Framework types:
class species
class atom			# atom should simply point at the species, not be a subclass of it- I will recase this later on, and make getters for the species pointer.
  Contains a distance function, called by the '-' operator, such that I can find the distance easily between two atoms.
  atom.x, atom.y and atom.z are a quick way to access the atom's coordinates.
  I've chosen to override the '+' operator to apply mbeta, a reflection based on how Fireball does it, returning an identical atom, but with the new coords. This means that we can apply the '-' operator on reflections very simply. 
  IN PLANNING- Override the '==' operator to confirm two atoms are the same atom (to avoid mbeta issues) 
  .ratom - returns a list of [atom.x, atom.y, atom.z]
  .charge
  .neutralCharge
  .chargeChange
  .neighbour(neighbor_number) - returns this atom's neighbor by that number. When mbeta for this neighbor is == 0, it returns the neighbor object, when mbeta != it returns a shallow copy of the neighbour's atom object with the coordinates of the image.
  .neighbourCount(element, depth = 1) - Returns the number of neighbours of species element. Default depth (as in 1 is first nearest neighbors, 2 2nd NN, etc) is 1 and ignorelist is for recursion. If you want only second nearest neighbours, for example, use a list comprehension.
  
class supercell - this works as a container for each of multiple supercells within a study (therefore I can compare features).
  .plotPDOS   # I should implement a separate "show PDOS" and "outputPDOS(filename)", I can make plotPDOS return a pyplot.figure argument to implement either quite simply.
  .genDOS       #Done
  .genDOSspec	
class shell
class shellPP
class orbital
class DOS
  __sum__ DOS object
  .offset(self, offsetVal) - Returns a new DOS object with the energies offset by - the value of offsetVal
  .plotme(self, HOMO, name = "DOSplot", legend = True, plotContributions = True, plotHOMO = True, plotOffset = False, xlim = []) - Plots to screen with option to plot to .png file of name [name]
  .writeToFile(filename = 'dens_out.dat') - Output to dens_file (for, e.g. dens_spec)
  PENDING average plotDOS objects (I can do this with a __divide__ magic method)
  PENDING Create PlotDOS object with average and contibuters.
  PENDING output to gnuPlot


TO be Created:

Lightning <-> TG converter (Not to be documented, but all calcs are done in Lightning-mode)
Dipole Calculator.
Angles between 3-atoms.
xsf tools.
readFireballIn(directory = '.') - PENDING - read in fireball.in and fill variables as needed, this can be augmented with atomList.
readFireballDirectory (directory = '.') 
Tool to sum DOS over specific atoms. Maybe plot them?
Tool to plot BS from the kylee-input files and the ek.dat?
readXYZ - Gotta be the same as readBas() with the final frame of the xyz file.

"""
P_abohr = float(0.529177)   # Bohr radii to Angstrom
# Atom and species classes are analgous to LIGHTNING, I'm probably going to alter them soon.

#-----------------------------------------------------------------------------
# Module classes:
# species, atom, shell, shellPP, orbital, DOS, supercell
#-----------------------------------------------------------------------------
class species():
  def __init__(self, mass, nZ, nssh, nssh_PP, norb_max, norb_PP_max, atomicE, rcutoff, rcutoff_PP, xmass, Zval, shells=[], shellsPP=[]):
    #I really don't think I need all of these.
    self.symbol = symbol(nZ)
    self.mass = mass
    self.nZ = nZ
    self.nssh = nssh
    self.nssh_PP = nssh_PP
    self.norb_max = norb_max
    self.norb_PP_max = norb_PP_max
    self.atomicE = atomicE
    self.rcutoff = rcutoff
    self.rcutoff_PP = rcutoff_PP 
    self.xmass = xmass
    self.Zval = Zval
    self.rcutoff_max = max(rcutoff)

    self.shells = shells
    self.shellsPP = shellsPP
#    self.orbitals = orbitals  Now in shells.
    self.rcutoffA = [cutoff * P_abohr for cutoff in rcutoff]
    self.rcutoff_max = max(rcutoff)
    self.rcutoffA_max = max(rcutoff) * P_abohr

    return
    
  def __str__(self):
    return self.symbol + " class object"
    
class atom():
  def __init__(self, atomNumber,element, ratom, Q = 999):
    self.atomNumber = atomNumber
    self.x = ratom[0]
    self.y = ratom[1]
    self.z = ratom[2]    
    # I'm was doing something wrong here, I hadn't yet worked out how to pass inherited values from an instance. I realised that I need to just make some getters and keep this away.
    # Look at: http://en.wikipedia.org/wiki/Creational_pattern, http://en.wikipedia.org/wiki/SOLID_(object-oriented_design), http://legacy.python.org/dev/peps/pep-0008/
    self.__element = element
    self.shells = []
    for ishell in element.shells:
      self.shells.append(shell(ishell.lssh, ishell.Qneutral, ishell.rcutoff, ishell.nafile, ishell.wffile, ishell.orbitals)) 
    
    self.shellsPP = []
    for ishellPP in element.shellsPP:
      self.shellsPP.append(shellPP(ishellPP.lssh))
    
#    self.orbitals = [] These are in the shell, need to confirm I've not made an error with some tests.
#    for iorb in element.orbitals:
#      self.orbitals.append(orbital(iorb.issh, iorb.l, iorb.m))
    
    return
    
  def neighbour(self, neighNum):
    #Returns the neighbour's atom object (by neighbour number), in the case of mbeta != 0, it returns a copy with modified xyz coordinates.
    try: 
      self.neighList
    except NameError:
      print "Call for a neighbour has been carried out, but this atom has no neighbour list."
      return NULL
  
    if self.neighList[neighNum][1] == 0 : return self.neighList[neighNum][0]
    else:
      # Here's where we have some fun augmenting a neighbour with the new position vectors.
      # I think now might be a good time to add the "+" modifier.
      neighbourImage = copy.copy(self.neighList[neighNum][0])

      neighbourImage + [self.xlm[1][self.neighList[neighNum][1]], self.xlm[2][self.neighList[neighNum][1]], self.xlm[3][self.neighList[neighNum][1]]]
      return neighbourImage 

  def neighbourCount(self, element, depth, ignorelist = []):
  # Class method implementation of the global method "NeighbourCount"
  # Difference here is it can be called by atom.neighboutCount(element, depth)
  # The ignorelist[] optional input allows for recursion when counting passed the 1st NN
  
  # First confirm that this instance of atom has got a neighbour map.
    try: 
      self.neighList
    except NameError:
      print "Call for a neighbour count has been carried out, but this atom has no neighbour list."
      return NULL
  # If I get an integer for the species of interest, convert to symbol.
    if (type(element) == 'int') : element = symbol(element)  
  # loop through the neighbours and count those of the species we want.
    
  # Count how many of those are the neighspecies.
    numberOfGuys = 0
    otherGuys = 0
    for jatom in self.neighList:
      if (jatom[0].symbol == element ) & ((jatom) not in ignorelist): 
        numberOfGuys = numberOfGuys + 1
        ignorelist.append(jatom)    
    
    if depth > 1: 
     for jatom in self.neighList:			#have to run loop again, otherwise ignorelist isn't right on first few sends. 
        newGuys = 0
        newGuys = jatom.neighbourCount(element, depth-1, ignorelist)
        otherGuys += newGuys
      
    return numberOfGuys + otherGuys
      
  def __str__(self):
    return "atom #" + str(self.atomNumber) + " -->  " + self.symbol + " @ (" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"
  
  def __sub__(self, other):
    # Some people will cringe at this- but I've superceded the subtraction operator with a distance operator between atoms, because '-' is kinda like a bond, right?
    #check that I'm not looking for the distance between a neighbor (and therefore need to use mbeta) or not:
    test = ("tu", "ple")
    if type(other) == type (test) : 
      print "Ah! we need a modification for mbeta in here!"
      print "use the neighborAsAtom method to convert neighbour to work with mbeta modified co-ords"
    return math.sqrt((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)
    
  def __add__(self, other):
    newRatom = []
    newRatom.append(self.x + other[0])
    newRatom.append(self.y + other[1])
    newRatom.append(self.z + other[2])
    return atom(self.atomNumber,self.__element, newRatom)

  def info(self):					#This is teporary, but means I can create a createDOS method that can take a supercell or an atomList
    return self.__class__.__name__

  #Rather than passing properties from the element class, I'm setting a bunch of pointers and getters.
  @property  
  def ratom(self):
    # Returns a list (ratom in Fireball) of the atoms xyz coordinates.
    return [self.x, self.y, self.z]
  @property
  def charge(self):
  # Returns the charge on the atom by summing the charge (Qcurrent) on each shell
    sumCharges = 0
    for ishell in self.shells:
      sumCharges += ishell.Qcurrent
    return sumCharges
  @property
  def chargeChange(self):
    return self.neutralCharge - self.charge
  @property
  def neutralCharge(self):
    sumCharges = 0
    for ishell in self.shells:
      sumCharges += ishell.Qneutral
    return sumCharges
  @property
  def symbol(self):
    return self.__element.symbol
  @property
  def mass(self):
    return self.__element.mass
  @property  
  def nZ(self):
    return self.__element.nZ
  @property  
  def nssh(self):
    return self.__element.nssh
  @property  
  def nssh_PP(self):
    return self.__element.nssh_PP
  @property  
  def norb_max(self):
    return self.__element.norb_max
  @property  
  def norb_PP_max(self):
    return self.__element.norb_PP_max
  @property  
  def atomicE(self):
    return self.__element.atomicE
  @property  
  def rcutoff(self):
    return self.__element.rcutoff
  @property  
  def rcutoff_PP(self):
    return self.__element.rcutoff_PP
  @property  
  def xmass(self):
    return self.__element.xmass
  @property  
  def Zval(self):
    return self.__element.Zval
  @property  
  def rcutoff_max(self):
    return self.__element.rcutoff_max
  @property  
  def rcutoffA(self):
    return self.__element.rcutoffA
  @property  
  def rcutoffA_max(self):
    return self.__element.rcutoffA_max
  @property
  def whatClass(self):
    return self.__class__.__name__
  
class shell():
  def __init__(self, lssh, Qneutral, rcutoff, nafile, wffile, orbitals=[], Qcurrent = 0):
  # Similar to the atom class, I should make the immutables here into pointers to a shell class that contains them. Qcurrent is the only mutable object in here.
    self.lssh = lssh
    self.Qneutral = Qneutral
    self.Qcurrent = Qcurrent
    self.rcutoff = rcutoff
    self.rcutoffA = rcutoff * P_abohr
    self.nafile = nafile
    self.wffile = wffile
    self.orbitals = []
    for orb in orbitals:
      self.orbitals.append(orbital(orb.l, orb.m))

    return

class shellPP():
  def __init__(self, lssh, cl=9999):
    self.lssh = lssh
    self.cl = cl
    return

class orbital():
  def __init__(self, l, m):
    self.l = l
    self.m = m
    return

class DOS():
  def __init__(self, energies, contributions, HOMO = 999, sumLine = True):
  #  print "Creating DOS object"
    # I'm having a bother working out what this should be.
    # Do I make a table, so as I have it in the same format as a dens_ file,
    # or do I point at the DOS contributions from each shell/orbital pair in the atom.
    # I can change the implementation later.
    for key in contributions:
      if len(energies) != len(contributions[key]):
        print "Error in creation of the DOS object- misalignment of energies."
        print "len(energies) =", len(energies), "len(contributions[key]) =", len(contributions[key]), "key=", key
        return
    self.energies = energies
    self.contributions = contributions
    if HOMO != 999 : self.HOMO = HOMO
    if sumLine:
      self.Totals = []
      for iEnergy in range(len(self.energies)):
        sum = 0
        for key, value in self.contributions.iteritems(): #why didnt I just use value here?
          sum += value[iEnergy]
        self.Totals.append(sum)
    return
    
  def __add__(self, other):
    # Check that energies align.
    if self.energies != other.energies:
      print "You cannot sum DOS info of different energies"
      return
    # Check that number of shells align.
    if set(self.contributions.keys()) != set(other.contributions.iterkeys()):  #set(d_1.keys()) == set(d_2.keys()) from http://stackoverflow.com/questions/3210832/pythonic-way-to-check-if-two-dictionaries-have-the-identical-set-of-keys
      print "You cannot sum DOS info of different orbitals"
      print "self.contributions.keys() != other.contributions.keys()"
      print self.contributions.keys()
      print other.contributions.keys()
      return 
    summedDOS = []
    if self.HOMO != other.HOMO:
      print "HOMO levels are different in each DOS being summed, returning a nulled HOMO"
      HOMO = 999
    else:
      HOMO = self.HOMO
    newcontributions = {}
    for key, value in self.contributions.iteritems():
      newcontributions[key] = []
      for iline in range(len(value)):
        newcontributions[key].append(self.contributions[key][iline] + other.contributions[key][iline])

#    for iLine in range(len(self.contributions)):
#      summedLine = []
#      for item in range(len(self.contributions[iLine])):
#        summedLine.append(self.contributions[iLine][item] + other.contributions[iLine][item])
#      summedDOS.append(summedLine)
    return DOS(self.energies, newcontributions, HOMO) #, sumLine = False)
   
  def __iadd__(self, other):
    return self + other
    
  def offset(self, offsetVal):
    # Returns a new DOS object with the energies offset by - the value of offsetVal
    newEnergies = [energy - offsetVal for energy in self.energies]
    try:
      return DOS(newEnergies, self.contributions, self.HOMO - offsetVal) 
    except AttributeError:   
      return DOS(newEnergies, self.contributions)
    
  def plotme(self, name = "DOSplot", printToFile = True, legend = True, plotContributions = True, plotHOMO = True, plotOffset = False, xlim = [], plotTOT = True, multiHOMO = {}, contriblw = 1, totlw = 2, legendLoc = 'best', returnAxes = False):
    if multiHOMO: colorList = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white'] # Allows me to make HOMOline the same color as the contributer
    if not self.HOMO: plotHOMO = False
    if self.HOMO == 999: plotHOMO = False
    if plotOffset: 
      newHOMO = 0
      print "You have chosen to offset the input DOS objects by their HOMOs... just so you know."
      if plotHOMO and not multiHOMO: offset = -self.HOMO
  #    plotEnergies = [energy - HOMO for energy in self.energies]
  #    newHOMO = 0.0
    else:
      plotEnergies = self.energies
      newHOMO = self.HOMO
    # I want a seperate TOTALS column so as I can highlight it.
    DOSTotals = self.Totals
    fig = plt.figure()
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
    if plotContributions:
      for key, value in sorted(self.contributions.iteritems()):  # Offset each contribution
        if plotOffset and multiHOMO: plotEnergies = [energy - multiHOMO[key] for energy in self.energies] 
        axes.plot(plotEnergies, value, lw = contriblw, label = key) # , label = ["i", "i2", "i3", "i4"]
    plotEnergies = self.energies # Return energies to "normal" once offset for the contributions is completed. 
    if plotTOT: axes.plot(plotEnergies, DOSTotals, lw =totlw, color = 'black', label = "Total")
    if legend :
      if legendLoc == 'Top':
      #Stuff to marker below from http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
      #Should put legend above the box, the URL gives below the box
      # Shink current axis's height by 10%  
        box = axes.get_position()
        axes.set_position([box.x0, box.y0 + box.height * 0.0, box.width, box.height * 0.9])

     # Put a legend above current axis
        axes.legend(loc='upper center', bbox_to_anchor=(0.5, 1.19), fancybox=True, shadow=True, ncol=3)
     # End of Stuff from URL

      else: axes.legend(loc=legendLoc, shadow = True)
    axes.grid(True)
    axes.axis('tight')
#  axes.set_title('PDOS')
    axes.set_xlabel('energy (eV)', fontsize=18)
    axes.set_ylabel('DOS (arb.)', fontsize=18)
    axes.get_yaxis().set_ticklabels([])
    if multiHOMO and not plotOffset:
      colorNumber = 0
      for thisHOMO in sorted(multiHOMO.itervalues()):
        axes.plot([thisHOMO, thisHOMO],[0,max(DOSTotals) * 1.2] , lw=2, color = colorList[colorNumber],  linestyle='dashed') # I need a smart way to deal with colors of the lines. 
        if colorNumber == 7: colorNumber = 0
        else: colorNumber += 1 
    elif plotHOMO : axes.plot([newHOMO, newHOMO],[0,max(DOSTotals) * 1.2] , color="green", lw=2, linestyle='dashed') #HOMO line
    if xlim: axes.set_xlim(xlim) 
    fig.show()
    if printToFile: plt.savefig(name + "DOS")
    if returnAxes: return axes # I THINK this will let me plot a list easier by constructing the axes object using plotme()... cool, huh?
    return

  def plotmeListVersion(self, HOMO, name = "DOSplot", printToFile = True, legend = True, plotContributions = True, plotHOMO = True, plotOffset = False, xlim = [], plotTOT = True):
    if plotOffset:
      print "You have chosen to offset the input DOS objects by their HOMOs. This creates temporary Offset DOS objects and does NOT overwrite the original data... just so you know."
      plotEnergies = [energy - HOMO for energy in self.energies]
      newHOMO = 0.0
    else:
      plotEnergies = self.energies
      newHOMO = HOMO
    # I want a seperate TOTALS column so as I can highlight it.
    DOSTotals = []
    for lines in self.contributions:
      DOSTotals.append(lines[-1])
    fig = plt.figure()
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # left, bottom, width, height (range 0 to 1)
    if plotContributions: axes.plot(plotEnergies,self.contributions, lw = 1) # , label = ["i", "i2", "i3", "i4"]
    if plotTOT: axes.plot(plotEnergies, DOSTotals, lw =2, color = 'black', label = "Total")
    if legend : axes.legend(loc=2, shadow = True)
    axes.grid(True)
    axes.axis('tight')
#  axes.set_title('PDOS')
    axes.set_xlabel('energy (eV)', fontsize=18)
    axes.set_ylabel('DOS (arb.)', fontsize=18)
    if plotHOMO : axes.plot([newHOMO, newHOMO],[0,max(DOSTotals) * 1.2] , color="green", lw=2, linestyle='dashed') #HOMO line
    if xlim: axes.set_xlim(xlim) 
    fig.show()
    if printToFile: plt.savefig(name + "DOS")
    return

  def subPlot3(self, position, fig, legend = True, plotContributions = True, plotHOMO = True, plotOffset = False, xlim = [], plotTOT = True, multiHOMO = {}, contriblw = 1, totlw = 2, legendLoc = 'best', returnAxes = False):
    # Part of plot3DOS, basically returns a pyplot.subplot(31Position) (3 rows, 1 column, plot#) for this DOS object, otherwise this is very similar to the DOS.plotme() method. If I get time to delve into the documentation a bit deeper, I'll work out how to do this with just the plotme() call.
    if multiHOMO: colorList = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white'] # Allows me to make HOMOline the same color as the contributer
    if not self.HOMO: plotHOMO = False
    if self.HOMO == 999: plotHOMO = False
    if plotOffset: 
      newHOMO = 0
      print "You have chosen to offset the input DOS objects by their HOMOs... just so you know."
      if plotHOMO and not multiHOMO: offset = -self.HOMO
  #    plotEnergies = [energy - HOMO for energy in self.energies]
  #    newHOMO = 0.0
    else:
      plotEnergies = self.energies
      newHOMO = self.HOMO
    # I want a seperate TOTALS column so as I can highlight it.
    DOSTotals = self.Totals
    threeonePos = 310 + position
    axes = plt.subplot(threeonePos)
    if plotContributions:
      for key, value in sorted(self.contributions.iteritems()):  # Offset each contribution
        if plotOffset and multiHOMO: plotEnergies = [energy - multiHOMO[key] for energy in self.energies] 
        axes.plot(plotEnergies, value, lw = contriblw, label = key) # , label = ["i", "i2", "i3", "i4"]
    plotEnergies = self.energies # Return energies to "normal" once offset for the contributions is completed. 
    if plotTOT: axes.plot(plotEnergies, DOSTotals, lw =totlw, color = 'black', label = "Total")
    if legend : axes.legend(loc=legendLoc, shadow = True, fontsize = 10)
    axes.grid(True)
    axes.axis('tight')
#  axes.set_title('PDOS')
    axes.set_xlabel('Energy (eV)', fontsize=18)
    if position == 2: axes.set_ylabel('Density of States (arb.)', fontsize=18)
    axes.get_yaxis().set_ticklabels([])
    if multiHOMO and not plotOffset:
      colorNumber = 0
      for key,thisHOMO in sorted(multiHOMO.iteritems()):
        print key, thisHOMO
        axes.plot([thisHOMO, thisHOMO],[0,max(DOSTotals) * 1.2] , lw=2, color = colorList[colorNumber],  linestyle='dashed') # I need a smart way to deal with colors of the lines. 
        if colorNumber == 7: colorNumber = 0
        else: colorNumber += 1 
    elif plotHOMO : axes.plot([newHOMO, newHOMO],[0,max(DOSTotals) * 1.2] , color="green", lw=2, linestyle='dashed') #HOMO line
    if xlim: axes.set_xlim(xlim) 
    #fig.show()
   
    return axes # I THINK this will let me plot a list easier by constructing the axes object using plotme()... cool, huh?    
       
  def writeToFileListVersion(self, filename = 'dens_out.dat'):
    densfile = open(filename, 'w')
    for iline in range(len(self.energies)):
      #densfile.write(str(self.energies[iline]) + '    ' + '    '.join(str(contrib) for contrib in self.contributions[iline]) + '\n')
      densfile.write('%f    ' % self.energies[iline] + '    '.join('%f' % contrib for contrib in self.contributionsList[iline]) + '\n')
    densfile.close() 
    return   
   
  def writeToFile(self, filename = 'dens_out.dat', header = True, columns = True, asList = False):
    if asList:
      self.contributionsList = self.makeContributionsList()
      self.writeToFileListVersion(filename)
    densfile = open(filename, 'w')
    if columns:
      if header: densfile.write('Energy    ' + '    '.join(key for key in sorted(self.contributions.iterkeys())) + '    Totals' + '\n')
      for iLine in range(len(self.energies)):
        densfile.write('%f    ' % self.energies[iLine] + '    '.join('%f' % value[iLine] for (key, value) in sorted(self.contributions.items()))  + '    %f' % self.Totals[iLine] + '\n')  #Sorting by key name (becuse it's shell/ orbital) from http://www.pythoncentral.io/how-to-sort-python-dictionaries-by-key-or-value/
      densfile.close()
      return
    else:
      if header: 
        densfile.write('Energy ' + '    '.join('%f' % contrib for contrib in self.energies) + '\n')
        for key, value in sorted(self.contributions.iteritems()):
          densfile.write(key + '    '.join('%f' % dens for dens in values) + '\n')
        return
      else:
        densfile.write('    '.join('%f' % contrib for contrib in self.energies) + '\n')
        for key, value in sorted(self.contributions.iteritems()):
          densfile.write('    '.join('%f' % dens for dens in values) + '\n')
        return

  def makeContributionsList(self):
    #This is taxing and may require some thought... but maybe not. I can construct a shell, orbital arrangement and I can see if there's a key for that, then reverse the append done in density.
    contribList = []
    for i in range(len(self.energies)):
      contribLine = []
      contribLine.append(self.energies[i])
      for ishell in range(4): #assumung a maximum of 4 shells 
        for iorb in range(2*ishell +1):
          #construct name as was done in atom2shell (and therefore this works for summed DOSes(but not catonated DOSes)
          contributionKey = "shell=" + str(ishell) + "_orbital=" + iorbital
          if contributionKey in self.contributions.keys(): contribLine.append(self.contributions[contributionKey][i])
      contribList.append(contribLine)
    return contribList      
          
            
  def sumtots(self, other):
    # Returns a sum of the Totals column of a list of DOS objects.
    # Check energies align.
    # sum the total columns of both.
    newTotal = []
    totalsByElement = {}
    for itotal in range(len(self.energies)):
      newTotal.append(self.contributions[itotal][-1] + other.contributions[itotal][-1])
    # if self contains previous element totals, pass that data through and append the list of elements.
    # if other contains previous element totals, pass that data through and append the list of elements.

  def cat(self, other):
    #Appends DOS object to the calling object, updates the TOTALS column.... this is a precurser to something to just do it by lists.
    if self.energies != other.energies:
      print "You cannot sum DOS info of different energies"
      return
    self.contributions.update(other.contributions)
    self.Total = []
    for iEnergy in len(self.energies):
      sum = 0
      for key, value in self.contributions.iteritems():
        sum += value[iEnergy]
      self.Total.append(sum)   
    return
      
class supercell(): # Basically, everything is centered around the atom class and atomList, this is more a "container" for that and related details.
    #def __init__(self, basisfile, natoms, nfragments, norbitals, norbitals_new, ztot, atom=[], nkpoints, kpoints=[], vector_lattice=[], xl=[], iquench, T_initial, T_final, T_instantaneous, iconstraint_rcm, iconstraint_vcm, iconstraint_L, iconstraint_KE, ifix_neighbors, rcm=[], rcm_old=[], vcm=[]):
  def __init__(self, directory, basisfile, atomList, lvs=[], xl=[], eigens = []): #ztot = 0.00, nkpoints = 0, kpoints = [], iquench = False, T_initial= 999.9, T_final = 999.9, rcm = [], eigens = []):
  #Another direct copy from the New Fireball Framework. Might need some stuff removed.
    self.directory = directory				# Path to working directory for this specific supercell
    self.basisfile = basisfile                  #name of the basis file    self.lvs = lvs                  #lattice vectors    self.xl = xl                  #the cell lattice vectors      
    self.atomList = atomList   
    self.eigens = eigens 			  # Energy eigenvalues
    return   
  
  @property
  def ztot(self):
  # Computes Ztot- the number of electrons in the system - form atomList.
  # This is a temporary routine until I implement the supercell class.
    ztot = 0
    for atom in self.atomList:
      for shell in atom.shells:
        ztot += shell.Qneutral      
    return ztot 
  
  @property
  def HOMO(self):
    return (self.eigens[(int(self.ztot/2)-1)])

  def clusteringFactor(self, element, rinput):
    return clustering(self, element, rinput)
    
  def varience(self, baseSpecies, seekSpecies, numSeek):
    return varience(self, baseSpecies, seekSpecies, numSeek)
  
  def calculateVdwNeighbours(self):
    neighVdW(self)
    return    
        
  def plotPDOS(self, elements, filename = 'supercellDOS'): 
    # Generate dens_TOT DOS object with dens_spec contributions.
    plotMeDOS = self.genDOS(elements)
    plotMeDOS.plotme(filename)
    return

  def genDOS(self, elements):
    # Generates a DOS object of contributions over the constituent elements in the cell.
    dosHere = {}
    for symbol in elements.keys():
      if [atom for atom in self.atomList if atom.symbol == symbol]:  #Skip if the element isn't in the supercell
        dosHere[symbol] = self.genDOSspec(symbol)
    returnDOS = catDOSTots(dosHere)   
    return returnDOS 
    
  def genDOSspec(self, element):
    # Generate list of atoms of type element.
    atomsOfInterest = [atom for atom in self.atomList if atom.symbol == element]
    if  not atomsOfInterest: 
      print "No atom of that element in this supercell"
      return
    quickList = [atom.atomNumber for atom in atomsOfInterest]
    DOSesOfInterest = [atom2DOS(atom, HOMO = self.HOMO) for atom in atomsOfInterest]
    DOSspec = DOSesOfInterest[0]
#    DOSspec = atom2DOS(atomsOfInterest[0], HOMO = self.HOMO)
    for idos in range(1, len(DOSesOfInterest)):
      DOSspec += DOSesOfInterest[idos]
    return DOSspec
    
  def draw(self, andIsoSurface = False, numNN = 8):
    from mayavi import mlab
    fig1 = mlab.figure(1, bgcolor=(0.0, 0.0, 0.0), size=(450, 500))
    mlab.clf()
    fig1.scene.disable_render = True
    for iatom in self.atomList:
      thisAtom = mlab.points3d([iatom.x], [iatom.y], [iatom.z], color=symbolColour(iatom.symbol))	# Add color for species here, going to use jmols colours
      sortedNeighbours = neighSortByDistance(self, iatom)
      for (ineigh, mbeta) in sortedNeighbours[0:numNN]:
        if mbeta == 0:
      	  thisBondX = [iatom.x, ineigh.x] 
      	  thisBondY = [iatom.y, ineigh.y] 
      	  thisBondZ = [iatom.z, ineigh.z]       	
          thisAtomsBonds = mlab.plot3d(thisBondX, thisBondY, thisBondZ ,tube_radius=0.1)
    #if andIsoSurface:
      # To be implemented
    fig1.scene_disable_render = False
    mlab.show()
    return
  
  @property  
  def etot(self, logfile = "runlog.log"):
    import os
    """ returns etot of the supercell, if it's already been saved, then returns that... otherwise looks for a log file
    and then grabs the final ETOT out of answer.xyz"""
    try:
      return self.ETOT
    except AttributeError:
      # okay, so we don't have an ETOT, lets open the answer.xyz
      command = "grep ETOT " + self.directory + logfile + " | tail -1 > " + self.directory + "/ETOT_tail.dat"
      os.system(command)
      ETOTfile = open(self.directory + "/ETOT_tail.dat", "r")
      self.ETOT = float(ETOTfile.readline().split()[-1])
      return self.ETOT
  
  @property
  def cohE(self, logfile = "runlog.log"):
    import os
    """ returns the cohesive energy of the supercell by checking for the cohE in the supercell object. If the cohE isn't there, 
    then it looks for a runlog file, firstly by checking for a file named in the logfile variable. Failing that, it grabs the 
    final entry of cohE for any file with that in the supercell directory"""
    try:
      return self.COHE
    except AttributeError:
      # okay, so we don't have an ETOT, lets open the answer.xyz
      command = "grep CohesiveE " + self.directory + logfile + " | tail -1 > " + self.directory + "/CohesiveE_tail.dat"
      os.system(command)
      CohEfile = open(self.directory + "/CohesiveE_tail.dat", "r")
      self.COHE = float(CohEfile.readline().split()[-1])
      return self.COHE

         
#-----------------------------------------------------------------------------
# Module methods.
# ReadBas, ReadLVS, ReadNeigh, ReadCharges, ReadDOS, density, atom2DOS, atomDistance,
# NeighborCount, NearestSpecList, symbol, basList2Object, chargesList2atoms,
# neighList2atomObjects, xlmbeta, computeZtot, readEigen, readInfo,
# plotDOSList, catDOSTots 
#-----------------------------------------------------------------------------
def ReadBas(basFileName):
  basFile = open(basFileName, 'r')
  atomList =[]
  #Discard Header
  print ""
  print "There are %d atoms in this supercell" %int(basFile.readline())

  #Loop over the rest of the file and read in species no., X, Y, Z
  for atomInfo in basFile.readlines():
    atomList.append(atomInfo.split())
  
  # Convert the strings to integers and floats
  for atom in atomList:
    atom[0] = int(atom[0])
    atom[1:] = [float(i) for i in atom[1:]]
  
  #Close File
  basFile.close()
  
  #return array
  return atomList
    
def ReadLVS(LVSFileName):
  LVSFile = open(LVSFileName, 'r')
  LVSList = []

  #Loop over the rest of the file and read in species no., X, Y, Z
  for dimension in LVSFile.readlines():
    LVSList.append(dimension.split())
  
  # Convert the strings to floats
  for dimension in LVSList:
    dimension[:] = [float(i) for i in dimension]
  
  # Close File
  LVSFile.close()
  
  #return list
  return LVSList
  
def ReadNeigh(directory = '.'):
  # Reads in the NEIGHBORS file and returns a list of lists of neighbors indexed over iatom.
  # TAKE CARE with this function, the traditional iatom was indexed from 1, whereas python 
  # indexes from 0, the "fix" I have applied is to subract 1 from each atom index on input. 
  # this means that you can find the xyz of a neighbor of an atom by: atomList[NeighborList[atom][neighNumber]][1:]
  # so, for atom number 13 (which is 12 in the array, so we're now calling it 12), and given the implementation of 
  # other subroutines in this module, we can simply address the atom number of neighbor 6 by:
  # NeighborList[12][5] and it's XYZ coordinates by atomList[NeighborList[12][5]][1:] or its species number by
  # atomList[NeighborList[12][5]][0]

  #Open NEIGHBORS file
  if directory == '.' : NeighFile = open ('NEIGHBORS', 'r')
  else: NeighFile = open (directory + '/NEIGHBORS', 'r')
  # The issue with the NEIGHBORS file is that it's very human-readable, it basically has an entry for every single atom which
  # begins with the number of neighbors, then there's four colums of data, corresponding to atom #, Neighbor #, mbeta and atom # of neighbor.
  # mbeta is the image cell.
  
  # To try and make this a little more clear, here is the section from FIREBALL-TG that carries out this readin:
  
#"""  
#! If there is a neighbor file from a previous run, or the user added one,
#! then initialize the input neighbor map accordingly.
#        if (nstepi .gt. 1) then
#
#        open (unit = 20, file = 'NEIGHBORS', status = 'unknown')
#        read (20,*)
#        do iatom = 1, natoms
#         read (20,*) num_neigh
#         neighn(iatom) = num_neigh
#         do ineigh = 1, num_neigh
#          read (20,*) katom, kneigh, mbeta, jatom
#          if (katom .ne. iatom .or. kneigh .ne. ineigh) then
#           write (*,*) ' iatom, katom = ', iatom, katom
#           write (*,*) ' ineigh, kneigh = ', ineigh, kneigh
#           write (*,*) ' Problem in NEIGHBORS, atoms not lined up correctly. '
#           write (*,*) ' Fix and restart! '
#           stop
#          end if
#          neigh_b(ineigh,iatom) = mbeta
#          neigh_j(ineigh,iatom) = jatom
#         end do
#        end do
#        close (unit = 20)
#"""
  # So, you can see, there's a basic set of "things" going on, that are easy enough to replicate.
  # For my current purposes, I only need the mbeta = 0 cases, but all are coming in this way for the time being.
  # I'm not sure if I want an option to return the "traditional" neigh_j and neigh_b arrays here also. Ideally I want
  # to output a list, indexed over iatom, of lists of neighbors. I can easily make this a list of lists of lists indexed
  # over iatom, mbeta, which might be the best thing to do.
  
  
  # Read in the first line, which is the number of atoms.
  natoms = int(NeighFile.readline().split()[0])
  
  #Create new lists
  NeighList = []
  mbeta = []
  jatom = []
  
  # read in info to the "old style" mbeta, jatom configuration and build lists for the lists of lists of lists outputs.
  for iatom in range(natoms):
    num_neigh = int(NeighFile.readline())
    jatom_temp = []
    mbeta_temp = []
 
    for ineigh in range(num_neigh):
      NeighString = NeighFile.readline()
      NeighSplit = NeighString.split()
      
      # We're only interested in NeighSplit[2] (mbeta) and NeighSplit[3] (jatom)
      jatom_temp.append(int(NeighSplit[3]))
      mbeta_temp.append(int(NeighSplit[2]))

    jatom.append(jatom_temp)
    mbeta.append(mbeta_temp)

   # Finally map the jatom/ mbeta set to be in the form of NeighList[atom][mbeta][Neighbor]
    
  for iatom in range(natoms):
    NeighList.append([])
    mbeta_max = max(mbeta[iatom])

    for mbeta_count in range(mbeta_max+1):
      NeighList[iatom].append([])
      num_neigh = len(jatom[iatom])
 
      for ineigh in range(num_neigh):
        if mbeta[iatom][ineigh] == mbeta_count:
          NeighList[iatom][mbeta_count].append(jatom[iatom][ineigh])


  NeighFile.close()
  return NeighList      

def ReadCharges(directory = '.'):
  # This reads in the CHARGES file and returns a List-of-Lists, again indexed over iatom.
  # Big advantage here over fortran is that it can total on the fly
  if directory == '.' : CHARGESFile = open ('CHARGES', 'r')
  else: CHARGESFile = open (directory + '/CHARGES', 'r')
  CHARGESList = []
  
  for CHARGELine in CHARGESFile.readlines()[1:]:
    CHARGESList.append(CHARGELine.split())

  # Convert the strings to floats and add a Total Charge column
  for CHARGELine in CHARGESList:
    CHARGELine[:] = [float(i) for i in CHARGELine]
  
  return CHARGESList

def readDOS(atomList, directory = '.'):
  # atomList is the result of readBas, and therefore a list of atoms. My current issue is if I create a whole new DOS object or append DOS information to the atoms. Dammit. Gotta think about this for a minute. Nope! atom class HAS to contain each atom / shell's DOS contribution, that way I can pull specific DOS metrics.
  # DOS info is arranged with the shells, each will contain a list of tuples containing E and DOS, the atom object will contain a DOS "getter", which will sum-on-the-fly it's shells. supercell object can then sum-on-the-fly with it's atoms.
      # Check that the 3 digit dens atom 1 exists:
  try:
    fmt = "%03d"
    filename = directory + "/dens_"+ fmt % 1 + ".dat"
    filecheck = open (filename, 'r')
    filecheck.close()
    print "3format"
  except IOError:    # 3 digit doesnt, check 4 digit
    try:			#There has got to be something more elegant than nested "try" statements.
      fmt = "%04d"
      filename = directory + "/dens_"+ fmt % 1 + ".dat"
      filecheck = open (filename, 'r')
      filecheck.close()
      print "4format"
    except IOError:
      print "Call for DOS read-in, but there's no dens_XXX.dat files in ", directory, fmt
      return
  print "opened a checkfile, fmt = ", fmt
  

  for atom in atomList:
    # Create the Energy list.
    atom.DOSEnergies = []
    atom.DOScontribution = []
    for ishell in atom.shells:
      for orbital in ishell.orbitals:
        orbital.DOScontribution = []
        
    #construct dens_ file name.
    densfile = open (directory + "/dens_"+ fmt % (atom.atomNumber +1) + ".dat", 'r')
    densList = densfile.readlines()
    densfile.close()
    for densLine in densList:
      E_Dens = densLine.split()
      icolumn = 0
      atom.DOSEnergies.append(float(E_Dens[0]))
      for ishell in atom.shells:
        for orbital in ishell.orbitals:
          icolumn += 1
          orbital.DOScontribution.append(float(E_Dens[icolumn]))
      icolumn += 1
      atom.DOScontribution.append(float(E_Dens[icolumn]))      
  return  

def density(atom):   #returns the density of states for a particular atom, returns a list-of-lists of Energy, DOS_orbital, Total DOS
  reconstructingDens = []
  for i in range(len(atom.DOSEnergies)):
    densLine = []
    densLine.append(atom.DOSEnergies[i])
    for ishell in atom.shells:
      for orbital in ishell.orbitals:
        densLine.append(orbital.DOScontribution[i])  
    densLine.append(atom.DOScontribution[i])
    reconstructingDens.append(densLine)
  return reconstructingDens

def atom2DOS_list(atom, HOMO = 999): # returns object of DOS class with an input of an atom class # Being retired for a version that places a dictionary with names of the shells in there instead.
  thisdensity = density(atom)
  thisEnergy = []
  thisContributions = []
  for line in thisdensity:
    thisEnergy.append(line[0])
    thisContributions.append(line[1:-1])
  return DOS(thisEnergy, thisContributions, HOMO)

def atom2DOS(atom, HOMO = 999):
# Returns a DOS object with a dictionary of lists (instead of the previous list of lists) for the contributions.
  #First thing, we need names of the orbitals.
#  thisdensity = density(atom) #Density if a list of lists whereby the list is a list of Energy, contrib1, contrib2, etc.
  thisEnergy = copy.copy(atom.DOSEnergies)
    
  thisContributions = {}
  for ishell in range(len(atom.shells)):
      for iorbital in range(len(atom.shells[ishell].orbitals)):
        contributionKey = "shell=" + str(ishell) + "_orbital=" + str(iorbital) 
        thisContributions[contributionKey] = copy.copy(atom.shells[ishell].orbitals[iorbital].DOScontribution)
  return DOS(thisEnergy, thisContributions, HOMO)
         
def atomDistance(atom1, atom2):
  return math.sqrt((atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2 + (atom1[3] - atom2[3])**2 )

def NeighborCount(Nlist, basList, iatom, neighspecies, depth, ignorelist = []):
  # Takes in a Neighbor List of the form output from ReadNeigh and a basis list of the form output by ReadBas, 
  # along with the atom concernd and the species we want and depth, returns the number of neighbors of that depth.
  # For recursion, the ignorlist means any already-counted atoms are ignored.
  if not ignorelist : ignorelist = []
  
  # Confirm iatom from the call.
  print "iatom is %d, which is %d", iatom, basList[iatom][0]
  # Grab list of neighbors from Nlist for that iatom
  ListOfAtoms = Nlist[iatom][0]

  # Count how many of those are the neighspecies.
  numberOfGuys = 0
  otherGuys = 0
  for jatom in ListOfAtoms:
    if (basList[jatom-1][0] == neighspecies ) & ((jatom - 1) not in ignorelist): 
      numberOfGuys = numberOfGuys + 1
      ignorelist.append(jatom - 1)    
    
  if depth > 1: 
   for jatom in ListOfAtoms:			#have to run loop again, otherwise ignorelist isn't right on first few sends. 
      newGuys = 0
      newGuys = NeighborCount(Nlist, basList, jatom-1, neighspecies, depth-1, ignorelist)
      otherGuys += newGuys
      
  return numberOfGuys + otherGuys

def NearestSpecList(Nlist, iatom, basList, neighSpecies, numneigh):
  # Returns an ordered list (distance low-high) of neighbors.
  # This method again uses the neighbor map, sweeps through to find
  # the neighbors that are of species neighSpecies, works out their 
  # distance from the iatom, then sorts that list and returns the 
  # requested numneigh. Unlike Neighborcount, it does not look to neighbors-of-neighbors
  ignorelist = []
  returnlist = []
  
  # I will need to rebuild an mbeta array of which vectors are applied to make the neighbors list.
  numMbetas = len(Nlist[iatom])
  print numMbetas
  
  return

def symbol(nZ):
  periodic = ['H','He','Li','Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar','K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er','Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra','Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf','Es', 'Fm', 'Md', 'No', 'Lw']
  if nZ == 0:
    print "You have just tried to access the name of an element with atomic no. 0"
    print "Consult a periodic table and start again"
    sys.exit(1)
  return periodic[nZ-1]
def symbolColour(elementName):
  #Returns a 3-tuple of the colours (from JMOL) for each element name. Default is 0.99, 0.99, 0.99
  coloursByElement = {"H" : (0.99609375,0.99609375,0.99609375), "He" : (0.84765625,0.99609375,0.99609375), "Li" : (0.796875,0.5,0.99609375), "Be" : (0.7578125,0.99609375,0), "B" : (0.99609375,0.70703125,0.70703125), "C" : (0.5625,0.5625,0.5625), "N" : (0.1875,0.3125,0.96875), "O" : (0.99609375,0.05078125,0.05078125), "F" : (0.5625,0.875,0.3125), "Ne" : (0.69921875,0.88671875,0.95703125), "Na" : (0.66796875,0.359375,0.9453125), "Mg" : (0.5390625,0.99609375,0), "Al" : (0.74609375,0.6484375,0.6484375), "Si" : (0.9375,0.78125,0.625), "P" : (0.99609375,0.5,0), "S" : (0.99609375,0.99609375,0.1875), "Cl" : (0.12109375,0.9375,0.12109375), "Ar" : (0.5,0.81640625,0.88671875), "K" : (0.55859375,0.25,0.828125), "Ca" : (0.23828125,0.99609375,0), "Sc" : (0.8984375,0.8984375,0.8984375), "Ti" : (0.74609375,0.7578125,0.77734375), "V" : (0.6484375,0.6484375,0.66796875), "Cr" : (0.5390625,0.59765625,0.77734375), "Mn" : (0.609375,0.4765625,0.77734375), "Fe" : (0.875,0.3984375,0.19921875), "Co" : (0.9375,0.5625,0.625), "Ni" : (0.3125,0.8125,0.3125), "Cu" : (0.78125,0.5,0.19921875), "Zn" : (0.48828125,0.5,0.6875), "Ga" : (0.7578125,0.55859375,0.55859375), "Ge" : (0.3984375,0.55859375,0.55859375), "As" : (0.73828125,0.5,0.88671875), "Se" : (0.99609375,0.62890625,0), "Br" : (0.6484375,0.16015625,0.16015625), "Kr" : (0.359375,0.71875,0.81640625), "Rb" : (0.4375,0.1796875,0.6875), "Sr" : (0,0.99609375,0), "Y" : (0.578125,0.99609375,0.99609375), "Zr" : (0.578125,0.875,0.875), "Nb" : (0.44921875,0.7578125,0.78515625), "Mo" : (0.328125,0.70703125,0.70703125), "Tc" : (0.23046875,0.6171875,0.6171875), "Ru" : (0.140625,0.55859375,0.55859375), "Rh" : (0.0390625,0.48828125,0.546875), "Pd" : (0,0.41015625,0.51953125), "Ag" : (0.75,0.75,0.75), "Cd" : (0.99609375,0.84765625,0.55859375), "In" : (0.6484375,0.45703125,0.44921875), "Sn" : (0.3984375,0.5,0.5), "Sb" : (0.6171875,0.38671875,0.70703125), "Te" : (0.828125,0.4765625,0), "I" : (0.578125,0,0.578125), "Xe" : (0.2578125,0.6171875,0.6875), "Cs" : (0.33984375,0.08984375,0.55859375), "Ba" : (0,0.78515625,0), "La" : (0.4375,0.828125,0.99609375), "Ce" : (0.99609375,0.99609375,0.77734375), "Pr" : (0.84765625,0.99609375,0.77734375), "Nd" : (0.77734375,0.99609375,0.77734375), "Pm" : (0.63671875,0.99609375,0.77734375), "Sm" : (0.55859375,0.99609375,0.77734375), "Eu" : (0.37890625,0.99609375,0.77734375), "Gd" : (0.26953125,0.99609375,0.77734375), "Tb" : (0.1875,0.99609375,0.77734375), "Dy" : (0.12109375,0.99609375,0.77734375), "Ho" : (0,0.99609375,0.609375), "Er" : (0,0.8984375,0.45703125), "Tm" : (0,0.828125,0.3203125), "Yb" : (0,0.74609375,0.21875), "Lu" : (0,0.66796875,0.140625), "Hf" : (0.30078125,0.7578125,0.99609375), "Ta" : (0.30078125,0.6484375,0.99609375), "W" : (0.12890625,0.578125,0.8359375), "Re" : (0.1484375,0.48828125,0.66796875), "Os" : (0.1484375,0.3984375,0.5859375), "Ir" : (0.08984375,0.328125,0.52734375), "Pt" : (0.8125,0.8125,0.875), "Au" : (0.99609375,0.81640625,0.13671875), "Hg" : (0.71875,0.71875,0.8125), "Tl" : (0.6484375,0.328125,0.30078125), "Pb" : (0.33984375,0.34765625,0.37890625), "Bi" : (0.6171875,0.30859375,0.70703125), "Po" : (0.66796875,0.359375,0), "At" : (0.45703125,0.30859375,0.26953125), "Rn" : (0.2578125,0.5078125,0.5859375), "Fr" : (0.2578125,0,0.3984375), "Ra" : (0,0.48828125,0), "Ac" : (0.4375,0.66796875,0.9765625), "Th" : (0,0.7265625,0.99609375), "Pa" : (0,0.62890625,0.99609375), "U" : (0,0.55859375,0.99609375), "Np" : (0,0.5,0.99609375), "Pu" : (0,0.41796875,0.99609375), "Am" : (0.328125,0.359375,0.9453125), "Cm" : (0.46875,0.359375,0.88671875), "Bk" : (0.5390625,0.30859375,0.88671875), "Cf" : (0.62890625,0.2109375,0.828125), "Es" : (0.69921875,0.12109375,0.828125), "Fm" : (0.69921875,0.12109375,0.7265625), "Md" : (0.69921875,0.05078125,0.6484375), "No" : (0.73828125,0.05078125,0.52734375), "Lr" : (0.77734375,0,0.3984375), "Rf" : (0.796875,0,0.34765625), "Db" : (0.81640625,0,0.30859375), "Sg" : (0.84765625,0,0.26953125), "Bh" : (0.875,0,0.21875), "Hs" : (0.8984375,0,0.1796875), "Mt" : (0.91796875,0,0.1484375)}
  if elementName in coloursByElement.keys():
    return coloursByElement[elementName]
  else: return (0.99, 0.99, 0.99)

def basList2Object(atomList, elements):
  ## Accepts the atomList argument, which is the list of atoms from ReadBas(), returns a list of atom objects. May add charges, Q, here.
  print "Converting atomList to atomList"
  objectList = []
  atomNumber = 0
  for iatom in atomList:
    elementName = symbol(iatom[0])
    objectList.append(atom(atomNumber, elements[elementName], iatom[1:]))
    atomNumber += 1
    
  return objectList       


def chargesList2atoms(atomList, charges):
  #accepts the list of atoms from basList2Object and adds the charges to each atom object.
  #Obviously, this is silly, because the atom instantiator can accepts a charges argument.
  # I'm going to iterate over a counter, I think. This way I can avoid embarassing addressing problems.
  
  #Make sure nothing's gone horribly wrong:
  if len(atomList) != len(charges): print "Something's gone horribly wrong"   
  
  for iatom in range (len(atomList)):
    for ishell in range(len(atomList[iatom].shells)):
      atomList[iatom].shells[ishell].Qcurrent = charges[iatom][ishell] 
  return

def neighList2atomObjects(atomList, neighList, xlm):
  # Accepts a list of atom objects and the neighbors list from readNeigh, adds a list of atoms neighbors to each atom.
  # This method also adds the xlm vectors to the atom object. This way I can neatly call the "atom.neighbour()" method later.
  # I'm going to iterate over a counter, I think. This way I can avoid embarassing addressing problems.
  # atom.neighList[0] is the atom object, atom.neighlist[1] is the mbeta value for these guys.
  	
  # Make sure nothings gone horribly horribly wrong
  if len(atomList) != len(neighList): print "Something's gone horribly horribly wrong"
  for iatom in atomList:
    iatom.neighList = []
    iatom.xlm = xlm
    #First things first, how many mbetas have we here:
    mbeta_max = len(neighList[iatom.atomNumber])
    for mbeta in range(mbeta_max):
      if neighList[iatom.atomNumber][mbeta]:
        for jatom in neighList[iatom.atomNumber][mbeta]:
          iatom.neighList.append((atomList[jatom -1], mbeta)) 
  return 


def xlmbeta(LVS):
  mbox = 4
  xxl = [0.0]*4 
  for i in range(0,4):
    xxl[i] = [0.0] * 729
    
# Initialize 27 "neighboring" cells at lattice vector xl(3,0:124).
# The units are angstrom, and NOT units of atomic lattice.
# The central cell is 0,0,0.
# We set up xxl in units of atomic lattice,
# then put xxl ------> xl in angstom.

  xxl [1][1] = 1.0

  xxl [1][2] = -1.0

  xxl [2][3] = 1.0

  xxl [2][4] = -1.0

  xxl [3][5] = 1.0

  xxl [3][6] = -1.0

  xxl [1][7] = 1.0
  xxl [2][7] = 1.0
  
  xxl [1][8] = -1.0
  xxl [2][8] = 1.0

  xxl [1][9] = -1.0
  xxl [2][9] = -1.0

  xxl [1][10] = 1.0
  xxl [2][10] = -1.0

  xxl [1][11] = 1.0
  xxl [2][11] = 1.0
  xxl [3][11] = 1.0

  xxl [1][12] = -1.0
  xxl [2][12] = 1.0
  xxl [3][12] = 1.0

  xxl [1][13] = -1.0
  xxl [2][13] = -1.0
  xxl [3][13] = 1.0

  xxl [1][14] = 1.0
  xxl [3][14] = 1.0  
  xxl [2][14] = -1.0

  xxl [1][15] = 1.0
  xxl [2][15] = 1.0    
  xxl [3][15] = -1.0

  xxl [1][16] = -1.0
  xxl [3][16] = -1.0
  xxl [2][16] = 1.0

  xxl [1][17] = -1.0
  xxl [2][17] = -1.0
  xxl [3][17] = -1.0

  xxl [1][18] = 1.0
  xxl [2][18] = -1.0
  xxl [3][18] = -1.0

  xxl [1][19] = 1.0
  xxl [3][19] = 1.0

  xxl [1][20] = -1.0
  xxl [3][20] = 1.0

  xxl [2][21] = 1.0
  xxl [3][21] = 1.0

  xxl [2][22] = -1.0
  xxl [3][22] = 1.0

  xxl [1][23] = 1.0
  xxl [3][23] = -1.0

  xxl [1][24] = -1.0
  xxl [3][24] = -1.0

  xxl [2][25] = 1.0
  xxl [3][25] = -1.0

  xxl [2][26] = -1.0
  xxl [3][26] = -1.0

# Set up the rest of xxl. Set up the bottom layer, the three middle
# layers and then the top layer. The inside 3X3X3 has already been set up
# by the data statement. Set up:
# (i) bottom layer at z=-2, and is 5X5.
# (ii) the three middle layers which enclose the inner 3X3.
# (iii) the top layer at z=2, and is 5X5.
# We have 26 lbeta's so far, so let's start at 26 and increase the number.
  lbeta = 26

# -4 layer for the cube 9x9x9
  for ix in range(-mbox, mbox +1):
    for iy in range (-mbox, mbox +1):
      lbeta = lbeta + 1
      xxl[1][lbeta] = float(ix)
      xxl[2][lbeta] = float(iy)
      xxl[3][lbeta] = -4.0

# -3 layer for the cube 7x7x7
  for ix in range(-mbox, mbox +1):
    for iy in range (-mbox, mbox +1):
      lbeta = lbeta + 1
      xxl[1][lbeta] = float(ix)
      xxl[2][lbeta] = float(iy)
      xxl[3][lbeta] = -3.0

# -2 layer for the cube 5x5x5
  for ix in range(-mbox, mbox +1):
    for iy in range (-mbox, mbox +1):
      lbeta = lbeta + 1
      xxl[1][lbeta] = float(ix)
      xxl[2][lbeta] = float(iy)
      xxl[3][lbeta] = -2.0

# Middle layers
  for midl in range (-1, 1 +1):
    for ix in range(-mbox, mbox +1):
      for iy in range (-mbox, mbox +1):
      
# Skip the "core" 3X3X3 cube, and only include the "skin"
         if (abs(ix) > 1) | (abs(iy) > 1) :
            lbeta = lbeta + 1
            xxl[1][lbeta] = float(ix)
            xxl[2][lbeta] = float(iy)
            xxl[3][lbeta] = float(midl)
     

# Now the 2.top layer, the cube 5x5x5
  for ix in range(-mbox, mbox +1):
    for iy in range (-mbox, mbox +1):
      lbeta = lbeta + 1
      xxl[1][lbeta] = float(ix)
      xxl[2][lbeta] = float(iy)
      xxl[3][lbeta] = 2.0

# Now the 3.top layer, the cube 7x7x7
  for ix in range(-mbox, mbox +1):
    for iy in range (-mbox, mbox +1):
      lbeta = lbeta + 1
      xxl[1][lbeta] = float(ix)
      xxl[2][lbeta] = float(iy)
      xxl[3][lbeta] = 3.0

# Now the 4.top layer, the cube 9x9x9
  for ix in range(-mbox, mbox +1):
    for iy in range (-mbox, mbox +1):
      lbeta = lbeta + 1
      xxl[1][lbeta] = float(ix)
      xxl[2][lbeta] = float(iy)
      xxl[3][lbeta] = 4.0

  mbeta_max = lbeta
  print 'mbeta_max = ', mbeta_max
# allocate xl
  xl = [0.0]*4 
  for i in range(0,4): 
    xl[i] = [0.0] * 729
    
# put xxl into xl with real angstrom units.
  for mbeta in range(0,mbeta_max):
   xl[1][mbeta] = xxl[1][mbeta]*LVS[0][0] + xxl[1][mbeta]*LVS[0][1] + xxl[3][mbeta]*LVS[0][2]		#I'm reversing this representation
   xl[2][mbeta] = xxl[1][mbeta]*LVS[1][0] + xxl[2][mbeta]*LVS[1][1] + xxl[3][mbeta]*LVS[1][2]   
   xl[3][mbeta] = xxl[1][mbeta]*LVS[2][0] + xxl[2][mbeta]*LVS[2][1] + xxl[3][mbeta]*LVS[2][2]
   
  return xl   

def makeCells(lvs):
  # I implemented the above method, xlmbeta, based on what I found in the FIREBALL source, I think that routine was in error from when I attempted to make the neighbours_vdW array. This is the PyFireball implementation of the Lightning make_cells function, which I think is more what I need.
  mbox = 4 
  mbetaMax = (2 * mbox + 1)**3 - 1

  xl = np.zeros([mbetaMax+1, 3])
  mbeta = 1
  for iz in range(-mbox, mbox +1):
    for iy in range(-mbox, mbox +1):
      for ix in range(-mbox, mbox +1):
        if (not ( ix == 0 and iy == 0 and iz == 0)):
          xl[mbeta] = float(ix) * lvs[0] +  float(iy) * lvs[1] + float(iz) * lvs[2]
          mbeta += 1
        
  return xl
  
def computeZtot(atomList):
  # Computes Ztot- the number of electrons in the system - form atomList.
  # This is a temporary routine until I implement the supercell class.
  ztot = 0
  for atom in atomList:
    for shell in atom.shells:
      ztot += shell.Qneutral
      
  return ztot 

def readEigen(directory = '.'):
  if directory == '.' : eigenFile = open ('eigen.dat', 'r')
  else: eigenFile = open (directory + '/eigen.dat', 'r')
  eigenList = []
  
  for eigenLine in eigenFile.readlines()[2:]:
    for eigenValue in eigenLine.split():
      eigenList.append(float(eigenValue))
  
  return eigenList
  
  
def readInfo(fdataLocation):
  """ PyFireball implementation of the readinfo.f90 and associated subroutines.
  Reads in the info.dat file from the Fdata and returns a directory of the 
  species, indexed by the element symbol"""
  # Declare species info dictionary.
  elements = {}
  
  infoLocation = fdataLocation + '/info.dat'
  try:
    infoFile = open(infoLocation, 'r')
  except IOError:
    print "  No info.dat file found at", fdataLocation
    print "  Please specify Fdata Location"
    sys.exit(1)

  # Say Hello
  print ' '
  print ' Reading from the info.dat file '
    
  info = infoFile.readlines()
  
  print ' This FData was created by: ', info[0]
  infoFile.close()
  nspecies = int(info[1].split()[0])
  print ' '
  print ' Number of species in database = ', nspecies

  # At this point, readinfo.f90 likes to not read in unused basis. I want to read it all in for a number of reasons.
  # The general form of the info.dat file is two lines of a header followed by each species information in this form:
  """
  ======================================================================
     1          - Information for this species 
    O           - Element 
      8         - Nuclear Z 
     15.999     - Atomic Mass 
     2          - Number of shells; L for each shell 
      0  1
     2          - Number of shells; L for each shell  (Pseudopotential) 
      0  1
     0.71 - Radial cutoffs PP 
    2.00   4.00
     3.40   3.80
    O/3/0.15/008_340.wf1       O/3/0.15/008_380.wf2     
    O/3/0.15/008_380.na0       O/3/0.15/008_340.na1       O/3/0.15/008_380.na2     
    -426.99209   - Atomic energy 
  ======================================================================
  """
  # In order to simplfy the read in, we are going to use relative line numbers. Therefore, there's an offset ( = 2 + 14*(ispecies - 1))
  
  for ispecies in range(nspecies):
    startLine = 2 + 16* ispecies
    elementLine = info[startLine + 2].split()
    NZLine = info[startLine + 3].split()
    massLine = info[startLine + 4].split()
    shellLine = info[startLine + 5].split()
    Lshells = info[startLine + 6].split()
    L_PPshellLine = info[startLine + 7].split()
    L_PPshells = info[startLine + 8].split()
    RcutoffPP = info[startLine + 9].split()
    Occupations = info[startLine + 10].split()
    cutoffs = map(float, info[startLine + 11].split())
    wfFiles = info[startLine + 12].split()
    naFiles = info[startLine + 13].split()
    energyLine = info[startLine + 14].split()    
    
    # Break down the info file into what we need for the species object.    
    Element = elementLine[0]
    NZ = int(NZLine[0])
    mass = float(massLine[0])
    nshells = int(shellLine[0])
    lssh = []
    norb_max = 0
    for ishell in range (nshells):
      lssh.append(int(Lshells[ishell]))
      norb_max += 2*(int(Lshells[ishell])) + 1
      
    nshellsPP = int(L_PPshellLine[0])
    lsshPP = []
    norbPP_max = 0
    for ishell in range (nshellsPP):
      lsshPP.append(int(L_PPshells[ishell]))
      norbPP_max += 2*(int(L_PPshells[ishell])) + 1
    rcutoffPP = float(RcutoffPP[0])
    Qneutral = []
    for ishell in range (nshells):
      Qneutral.append(float(Occupations[ishell]))

    atomicE = float(energyLine[0])
    Zval = int(sum(lssh))

    # Create the shells for this species. Thank goodness this need only be done once.
    # Of course, I have all the data, so this should be fairly simple.(later: ... famous last words- I've messed something up! orbitals?)
    shells = []
    for i in range(0, nshells):
      orbitals = []    
      for iorb in range(2*lssh[i] + 1):
        orbitals.append(orbital(lssh[i], iorb - lssh[i]))
      shells.append(shell(lssh[i], Qneutral[i], cutoffs[i], naFiles[i], wfFiles[i], orbitals))
    shellsPP =[]
    for i in range(0, nshellsPP):
      shellsPP.append(shellPP(lsshPP[i]))

    # Create species object.
    elements[Element] = species(mass, NZ, nshells, nshellsPP, norb_max, norbPP_max, atomicE, cutoffs, RcutoffPP, mass, Zval, shells, shellsPP)
    
  print elements.keys()

  return elements

def plotDOSList(densList, name = "multiDOSplot", legend = True, plotContributions = True, plotHOMO = True, plotOffset = True, xlim = [], saveFile = True): #Need a xlim in here, too. 
  if len(densList) != 3 :
    print "plotDOSList is supposed to only deal with Lists of 3 density objects"
    return
  if plotOffset:
    print "You have chosen to offset the input DOS objects by their HOMOs. This creates temporary Offset DOS objects and does NOT overwrite the original data... just so you know."
    # Loop and call offset.
    newList =[]
    for DOSobject in densList:
      newList.append(DOSobject.offset(DOSobject.HOMO))
    densList = newList
  f = plt.figure()
  plt.subplots_adjust(hspace = 0.001)
  
  ax1 = plt.subplot(311)
  if plotContributions:
    for key, value in densList[0].contributions.iteritems():
      ax1.plot(densList[0].energies, value, lw = 1, label = key)
  #Totals line
  ax1.plot(densList[0].energies, densList[0].Totals, lw = 2, color = 'black')
  ax1.grid(True)
  ax1.get_yaxis().set_ticklabels([])
  ax1.set_xlabel('energy (eV)', fontsize = 18)
  ax1.axis('tight')
  if plotHOMO: ax1.plot([densList[0].HOMO, densList[0].HOMO],[0,max(densList[0].Totals) * 1.2] , color="green", lw=2, linestyle='dashed')
  if xlim: ax1.set_xlim(xlim) 
  
  ax2 = plt.subplot(312, sharex = ax1)
  if plotContributions:
    for key, value in densList[1].contributions.iteritems(): 
      ax2.plot(densList[1].energies, value, lw = 1, label = key)
  #Totals line
  ax2.plot(densList[1].energies, densList[1].Totals, lw = 2, color = 'black')  
  ax2.grid(True)
  ax2.get_yaxis().set_ticklabels([])
  ax2.set_xlabel('energy (eV)', fontsize = 18)
  ax2.set_ylabel('DOS (arb.)', fontsize=18)
  ax2.axis('tight')
  if plotHOMO: ax2.plot([densList[1].HOMO, densList[1].HOMO],[0,max(densList[1].Totals) * 1.2] , color="green", lw=2, linestyle='dashed')
  if xlim: ax2.set_xlim(xlim)
  
  ax3 = plt.subplot(313, sharex = ax1)
  if plotContributions:
    for key, value in densList[2].contributions.iteritems(): 
      ax3.plot(densList[2].energies, value, lw = 1, label = key)
  #Totals line
  ax3.plot(densList[2].energies, densList[2].Totals, lw = 2, color = 'black')
  ax3.grid(True)
  ax3.get_yaxis().set_ticklabels([])
  if plotHOMO: ax3.plot([densList[2].HOMO, densList[2].HOMO],[0,max(densList[2].Totals) * 1.2] , color="green", lw=2, linestyle='dashed')
  if xlim: ax3.set_xlim(xlim)
  ax3.set_xlabel('energy (eV)', fontsize = 18)  
  #ax3.axis('tight')  #Uncomment if desired.
  
  xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
  plt.setp(xticklabels, visible=False)
  f.show()
  if saveFile: f.savefig(name + "multiDOS")
  return

def plot3DOS(densList, name = "multiDOSplot", legend = True, plotContributions = True, plotHOMO = True, xlim = [], saveFile = True, subLegend = False, plotOffset = False, plotTOT = False, multiHOMO = [], contriblw = 1, totlw = 2, legendLoc = 'best', subLegendLoc = 'best'): #Need a xlim in here, too. 
# TO Do: Sort out legend size, sort out HOMO lines for passing multiHOMO, sort 
  # Takes a list of three DOS objects. Creates a nice stacked plot of them. This had better work.
  if len(densList) != 3 :
    print "plotDOSList is supposed to only deal with Lists of 3 density objects"
    return
  f = plt.figure()
  plt.subplots_adjust(hspace = 0.001)
  if multiHOMO: 
    multiHOMO1 = multiHOMO[0]
    multiHOMO2 = multiHOMO[1]
    multiHOMO3 = multiHOMO[2]
  else: 
    multiHOMO1 = {}
    multiHOMO2 = {}
    multiHOMO3 = {}
  ax1 = densList[0].subPlot3(1, f, legend = subLegend, plotContributions = plotContributions, plotHOMO = plotHOMO, plotOffset = plotOffset, xlim = xlim, plotTOT = plotTOT, multiHOMO = multiHOMO1, contriblw = contriblw, totlw = totlw, legendLoc = subLegendLoc)
  ax2 = densList[1].subPlot3(2, f, legend = subLegend, plotContributions = plotContributions, plotHOMO = plotHOMO, plotOffset = plotOffset, xlim = xlim, plotTOT = plotTOT, multiHOMO = multiHOMO2, contriblw = contriblw, totlw = totlw, legendLoc = subLegendLoc)
  ax3 = densList[2].subPlot3(3, f, legend = subLegend, plotContributions = plotContributions, plotHOMO = plotHOMO, plotOffset = plotOffset, xlim = xlim, plotTOT = plotTOT, multiHOMO = multiHOMO3, contriblw = contriblw, totlw = totlw, legendLoc = subLegendLoc)
  
  # That was easy, next is to change the few small things for the figure to work:
  xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
  plt.setp(xticklabels, visible=False)
  f.show()
  if saveFile: f.savefig(name + "3DOS")  
  
  # Hmm. I think I want a method to return the axes objects for 
  return


def directory2Supercell(elements, directory = '.', basFile = 'answer.bas', lvsFile = 'input.lvs'):
  #Create empty supercell (this is cleaner for generation, IMHO)
  directoryCell = supercell(directory,'',[],[],[],[])
  
  #Fill in all the blanks as required from the directory.
  directoryCell.BasList = ReadBas(directory + '/' + basFile)
  directoryCell.LVSList = ReadLVS(directory + '/' + lvsFile)
  directoryCell.lvs = np.asarray(directoryCell.LVSList)
  directoryCell.neighList = ReadNeigh(directory)
  directoryCell.CHARGES = ReadCharges(directory)
  #directoryCell.xlm = xlmbeta(directoryCell.LVSList)
  directoryCell.xl = makeCells(directoryCell.lvs)  
  directoryCell.atomList = basList2Object(directoryCell.BasList, elements) #I'm assuming elements gets passed through via the scope
  neighList2atomObjects(directoryCell.atomList, directoryCell.neighList, directoryCell.xl)
  chargesList2atoms(directoryCell.atomList, directoryCell.CHARGES)
  readDOS(directoryCell.atomList, directory)
  directoryCell.eigens = readEigen(directory)
  readInNeighVdw(directoryCell)
  return directoryCell

def directory2SupercellLiteVersion(elements, directory = '.', basFile = 'answer.bas', lvsFile = 'input.lvs', xl = None):
  """ Very similar to the above, but doesn't read in DOS or anything I might not yet need."""
  #Create empty supercell (this is cleaner for generation, IMHO)
  directoryCell = supercell(directory,'',[],[],[],[])
  
  #Fill in all the blanks as required from the directory.
  directoryCell.BasList = ReadBas(directory + '/' + basFile)
  directoryCell.LVSList = ReadLVS(directory + '/' + lvsFile)
  directoryCell.lvs = np.asarray(directoryCell.LVSList)
  directoryCell.CHARGES = ReadCharges(directory)
  if xl: directoryCell.xl = xl
  else: directoryCell.xl = makeCells(directoryCell.lvs)  
  directoryCell.atomList = basList2Object(directoryCell.BasList, elements) #I'm assuming elements gets passed through via the scope
  chargesList2atoms(directoryCell.atomList, directoryCell.CHARGES)
  directoryCell.eigens = readEigen(directory)
  readInNeighVdw(directoryCell)
  return directoryCell

def catDOSTots(DOSdict, passHOMO = True):
  # Accepts a dictionary of DOS objects and returns a DOS object with the Totals of each dictionary object as a contribution.
  # To avoid dumb mistakes, first we loop through all key, value tuples and confirm that the energies are the same (because this is the addressing)
  for value1 in DOSdict.itervalues():
    for value2 in DOSdict.itervalues():
      if value1.energies != value2.energies:
        print "Misalignment of energies in catDOSTots"
        print "This can lead to strange errors."
        print "Either re-run FIREBALL with aligned DOS"
        print "Or bin your data in the DOS objects"
        return
  returnContributions = {}
  for key in DOSdict.keys():
    returnContributions[key] = DOSdict[key].Totals
    returnEnergies = DOSdict[key].energies		#These two lines ares not necessarily the most efficient, but they work.
    returnHOMO = DOSdict[key].HOMO
  if passHOMO: return DOS(returnEnergies, returnContributions, returnHOMO)
  else: return DOS(returnEnergies, returnContributions)
  
def editFireballIn(cell, changes, iterateNsteps = True):
#"""Edits supercell's directories fireball.in, passes through all flags
#   except the declared flags in the changes dictionary, which are either edited or added"""
# Lists of parameters in each section in fireball.in
  import os
  # move fireball.in to fireball.in.old
  command = "mv " + cell.directory + "/fireball.in " + cell.directory +"/fireball.in.old"
  
  os.system(command) 
  # open fireball.in.old
  fireballIn = open(cell.directory + "/fireball.in.old", "r")

  # check for lines that need to be added.
  appendedLines = []
  for fireballLine in fireballIn.readlines():
    fireballLine = fireballLine.split()
    if fireballLine[0] in changes:
      fireballLine[2] = str(changes[fireballLine[0]])
      del changes[fireballLine[0]]
    if (fireballLine[0] == 'nstepi' or fireballLine[0] == 'nstepf') and iterateNsteps == True:
      fireballLine[2] = str(int(fireballLine[2]) +1)
    appendedLines.append(fireballLine)
    
   # writeout to file, add in any dictionary items remaining as you do.
  fireballOutFile = open(cell.directory + "/fireball.in", "w")
  for fireballLine in appendedLines: #fireballIn:
    deleteList = []
    if fireballLine[0] == "&OPTION":
      fireballOutFile.write(" ".join(fireballLine) + "\n")
      for key, value in changes.iteritems():
        if fireballInSection(key) == "OPTIONS":
          fireballOutFile.write(key + "=" + str(value) +"\n")
          deleteList.append(key)
    elif fireballLine[0] == "&OUTPUT":
      fireballOutFile.write(" ".join(fireballLine) + "\n")
      for key, value in changes.iteritems():
        if fireballInSection(key) == "OUTPUTS":
          fireballOutFile.write(key + "=" + str(value) +"\n")
          deleteList.append(key)
    elif fireballLine[0] == "&QUENCH":
      fireballOutFile.write(" ".join(fireballLine) + "\n")
      for key, value in changes.iteritems():
        if fireballInSection(key) == "QUENCH":
          fireballOutFile.write(key + "=" + str(value) +"\n")
          deleteList.append(key)
    elif fireballLine[0] == "&MESH":
      fireballOutFile.write(" ".join(fireballLine) + "\n")
      for key, value in changes.iteritems():
        if fireballInSection(key) == "MESH":
          fireballOutFile.write(key + "=" + str(value) +"\n")
          deleteList.append(key)
    elif fireballLine[0] == "&TDS":
      fireballOutFile.write(" ".join(fireballLine) + "\n")
      for key, value in changes.iteritems():
        if fireballInSection(key) == "TDS":
          fireballOutFile.write(key + "=" + str(value) +"\n")
          deleteList.append(key)
    else: fireballOutFile.write(" ".join(fireballLine) + "\n")
    for delMe in deleteList: del changes[delMe]
    
  if changes: 
    print "ERROR: label not in fireball Labels"
    print changes
  
  fireballOutFile.close()
  fireballIn.close()
  return
  
def fireballInSection(label):
#""" Returns the SECTION name for the label for the fireball.in file"""
  OPTIONS = ["iharris", "idogs", "ihubbard", "ihorsfield", "imcweda", "iks", "igsn", "iqout", "qstate", "iquench", "icluster", "iensemble", "ifixcharge", "ifixneigh", "iumbrella", "ibarrier", "ivdw", "iimage", "idynmat", "iharmonic", "iconstraints(1)", "iconstraints(2)", "iconstraints(3)", "iconstraints(4)", "iendtemp", "ineb", "itrans", "basisfile", "lvsfile", "kptpreference", "acfile", "xvfile", "nstepi", "nstepf", "dt", "T_initial", "T_final", "max_scf_iterations", "bmix", "ialgmix", "sigmatol", "tempfe", "rescal", "itdse", "ibias", "rescal", "xyz2line", "imdet", "nddt"]
  OUTPUTS = ["iwrtcdcoefs", "iwrtcharges", "iwrtdensity", "iwrteigen", "iwrtefermi", "iwrtfpieces", "iwrthampiece", "iwrtcomponents", "iwrtneigh", "iwrtneigh_com", "iwrtxyz", "iwrtdos", "iwrthop", "iwrtatom", "iwrtpop", "iwrtHS", "iwrtvel", "iwrtden", "iwrtewf", "iwrtxsf", "idensimport", "iwrtpsit", "iwrtqt"]
  QUENCH = ["energytol", "forcetol", "T_want", "taurelax"]
  MESH = ["Ecut", "iewform", "npbands", "ewfewin_min", "ewfewin_max", "ifixg0", "g0(1:3)"]
  TDS = ["netime", "nexcite", "idelec", "hoccup", "eband", "np2es", "isp2es"]
  if label in OPTIONS: return "OPTIONS"
  elif label in OUTPUTS: return "OUTPUTS"
  elif label in QUENCH: return "QUENCH"
  elif label in MESH: return "MESH"
  elif label in TDE: return "TDS"
  else : return None

def supercellDirectory2DOScalc(elements, cell):  
  # Instead of supercell2DOScalc, this just creates the dos.optional file and edits the fireball.in file.
  # create the dos.optional file
  writeDOSinput(cell%directory)
  
  # edit fireball.in
  changes_to_fireballIn = {"iwrtdos" : 1, "iwrtcharges" : 1, "wrtcdcoefs" : 1, "iwrteigen" : 1, "iwrtpop" : 1}  
  editFireballIn(cell, changes_to_fireballIn, iterateNsteps = True)
  
  # Thats it... you need to rewrite the job script or call fireball.x now
  return 
  
def writeDOSinput(outDirectory, cell, atomStart = 1, atomEnd = 0, initEnergy=-4.00, energyStep=0.02, nEnergySteps=400):
  """ Writes out a Fireball-stlye dos.input file to the outDirectory, this is of the form:
  1.0                   ! scale factor (leave 1.0)
  1       108           ! list of atoms to analyze DOS
  400                   ! number of energy steps
  -2.0   0.01          ! initial energy, energy step 
  0                     ! leave untouched
  0.0     0.0           ! leave untouched
  0.05                  ! imaginary part of Green function (controls energy level smearing)
  """
  if atomEnd == 0 : atomEnd = len(cell.atomList)

  dosInputFileName = outDirectory + "/dos.optional"
  dosInputFile = open(dosInputFileName, 'w')
  dosInputFile.write("0                   ! scale factor (leave 1.0)")
  dosInputFile.write(str(atomStart) + "      " + str(atomEnd) +"          ! list of atoms to analyze DOS")
  dosInputFile.write(str(nEnergySteps) + "                  ! number of energy steps")
  dosInputFile.write(str(initEnergy) + "      " + str(energyStep) + "         ! initial energy, energy step ")
  dosInputFile.write("0                     ! leave untouched")
  dosInputFile.write("0.0     0.0           ! leave untouched")
  dosInputFile.write("0.05                  ! imaginary part of Green function (controls energy level smearing)")
  dosInputFile.close()
  return


def neighVdW(supercell):
  # Equivilant of FIREBALL's neighbors_vdW routine, these neighbours are very different from the standard neighbors.
  # This routine takes >> 1 hour on my MacBook Air 2011, so I'm going to make a write out method so as the calculated neighbors need only be
  # done once for any supercell.   # This method needs to become a PyPy or equivilant module, it's waaaaayyy to slow uncompiled (about an hour per supercell) Hmm- Maybe just make a fortran prog and call it from system
  import time
  rangeMax = 20.0
  # Loop through all atoms.
  for iatom in supercell.atomList:
    iatom.neighVdw = []   
  # Loop through all possible xlms
    for mbeta in range(len(supercell.xl)):
    
  # Loop through all atoms, except the current one.
      for candidateNeighAtom in supercell.atomList:
        startTime = time.time()
        print "In loop, start time"
        if candidateNeighAtom == iatom and mbeta == 0: continue
             
        #jatom = candidateNeighAtom + supercell.xl[mbeta]		# Rather than making a whole atom(takes ages), lets just take the positions.
        jatomPosition_x = candidateNeighAtom.x + supercell.xl[mbeta][0]
        jatomPosition_y = candidateNeighAtom.y + supercell.xl[mbeta][1]
        jatomPosition_z = candidateNeighAtom.z + supercell.xl[mbeta][2] 
        
        print "Time to do jatom position =", preTime - postTime
        z = math.sqrt((iatom.x - jatomPosition_x)**2 + (iatom.y - jatomPosition_y)**2 + (iatom.z - jatomPosition_z)**2)
        #z = iatom - jatom
  # add list to this atom.
        checktime = time.time()
        if z < rangeMax:
          iatom.neighVdw.append((candidateNeighAtom, mbeta))
        print "ifTime = ",  checktime -time.time() 
        print "Looptime =", startTime - time.time()
        #return
  return

def writeoutNeighVdw(supercell):
  # Routine to write out the neighbours_vdW information, the find_neighbors is a very time consuming process and therefore this is a shortcut to
  # restarting analysis as needed. Also, it should be simple enough, the neighbours are per atom as a list of tuples of atom, mbeta. I therefore 
  # only need those tuples per atom for this to be able to be read in per atom.
  neighFileName = supercell.directory + "/neighFilePF.dat"
  print neighFileName
  neighFile = open(neighFileName, 'w')
  for iatom in supercell.atomList:
    neighFile.write("atom: "+ str(iatom.atomNumber) + "\n") 
    for neighAtom in iatom.neighVdw:
      neighFile.write(str(neighAtom[0].atomNumber) + "  " + str(neighAtom[1]) + "\n")
  return 

def readInNeighVdw(supercell):
  # Routine to write out the neighbours_vdW information, the find_neighbors is a very time consuming process and therefore this is a shortcut to
  # restarting analysis as needed. Also, it should be simple enough, the neighbours are per atom as a list of tuples of atom, mbeta. I therefore 
  # only need those tuples per atom for this to be able to be read in per atom.
  neighFileName = supercell.directory + "/neighFilePF.dat"
  neighFile = open(neighFileName, 'r')
  for neighInfo in neighFile.readlines():
    neighLine = neighInfo.split()
    if neighLine[0] == "atom:": 
      iatom = supercell.atomList[int(neighLine[1])]
      iatom.neighVdw = []
    else:
      iatom.neighVdw.append((supercell.atomList[int(neighLine[0])], int(neighLine[1])))
  return

def neighDist(supercell, iatom, neighAtom):
  jatom = neighAtom[0] + supercell.xl[neighAtom[1]]
  return iatom - jatom
  
def neighSortByDistance(supercell, iatom):
  SortDistList = sorted(iatom.neighVdw, key = lambda jatom: neighDist(supercell, iatom, jatom))
  return SortDistList

def clustering (supercell, speciesOfInterest, rinput):		# Not sure if I need elements here, but anyways, and I want rinput out somehow.
  # This Bloody clustering algorithm uses the "neighbors_vdw" neighbor map, rather than the actual neighbor map.
  # I have to implement that, now. Going to call it neighVDW or something, it's the same, but with a cutoff of 20Ang, it seems.
  
  # so-called mean-squate list
  msList = []
  # Loop through all atoms.
  
  # If one is a speciesOfInterest: Loop through that atoms neighbors
  speciesOfInterestList = [atom for atom in supercell.atomList if atom.symbol == speciesOfInterest]
  maxClusterSize = len(speciesOfInterestList)
  for iatom in speciesOfInterestList:
    neighboursOfInterestList = [jatom for jatom in iatom.neighVdw if jatom[0].symbol == speciesOfInterest]
    for neighAtom in neighboursOfInterestList:
      mbeta = neighAtom[1]
      #print mbeta    
      jatom = neighAtom[0] + supercell.xl[mbeta]	# Be sure to apply new position if jatom is in a reflection cell
       
  # If neighbour also == speciesOfInterest:
  # Find distance between the two sites
  
# From FIREBALL:
#                    ms_test = (Rinput - z)**2
#                    ms_max = maxval(ms(1:max_cluster_size))
#                    ! if ms_test is smaller than any of the components in ms,
#                    ! then replace the maximum component by this number
#                    if (ms_test .lt. ms_max) then
#                      index_ms_max = maxloc(ms(1:max_cluster_size),1)
#                      ms(index_ms_max) = ms_test

# Max_Cluster_Size = # of SpeciesOfInterest -1 in FIREBALL code.

# Pythonic method:

# find distance between these neighbor pairs (check that iatom != jatom and mbeta != 0)


# append (Rinput - distance(iatom, jatom))^2 to msList
      #print rinput, iatom, jatom, iatom - jatom
      msList.append((rinput - (iatom - jatom))**2)
#end loops over atoms (do I need that?)
  
# final ms is the first max_cluster_size of atoms in the sorted msList (sort low-high), sum these, divide by max_cluster_size, then SQRT
  msList.sort() 
  clustering = math.sqrt((sum(msList[0:maxClusterSize]))/maxClusterSize)
  
  return clustering

def sortAllNeighboursByDistance(cell):
  """ Resorts all neighbourList's of all atoms in the supercell by distance """
  for iatom in cell.atomList:
    iatom.neighVdw = neighSortByDistance(cell, iatom)

def varience_array(supercell, baseSpecies, seekSpecies, numSeek):
  """ returns the array of distances from baseSpecies to all atoms in seekSpecies (as a list of species) for all baseSpecies
  (single argument) atoms in the supercell. numSeek defines how many of the closesest neighbour atoms of the seekSpecies kinds to use."""
  ### Quick list comprehension to get a list of the baseSpecies
  baseSpeciesList = [iatom for iatom in supercell.atomList if iatom.symbol == baseSpecies]
  distList = []
  for iatom in baseSpeciesList:
  # Generate list of neighbours that are part of the seekSpecies:
    seekNeighbours = [ineigh for ineigh in iatom.neighVdw if ineigh[0].symbol in seekSpecies]
    # append to list of the required atoms:
    tempList = []
    for ineigh in seekNeighbours:
      tempList.append(neighDist(supercell, iatom, ineigh))  
    distList += sorted(tempList)[0:numSeek]
  distArray = np.asarray(distList)
  return distArray
  
def varience(supercell, baseSpecies, seekSpecies, numSeek):
  """ returns the mean and varience of distance from baseSpecies to all atoms in seekSpecies (as a list of species) for all baseSpecies
  (single argument) atoms in the supercell. numSeek defines how many of the closesest neighbour atoms of the seekSpecies kinds to use. """
  distList = varience_array(supercell, baseSpecies, seekSpecies, numSeek)
  return (distList.mean(), distList.var())        

def drawAtoms(supercell, atomList, andIsoSurface = False, numNN = 6, implicitBonds = False, chargeCloud = False, chargeLabel = False, nameLabel = False, speciesLabel = False, numLabel = False, atomsAlso = True, secretZero = False):
  from mayavi import mlab
  fig1 = mlab.figure(1, bgcolor=(0.0, 0.0, 0.0), size=(450, 500))
  mlab.clf()
  fig1.scene.disable_render = True
  previousAtomDrawn = False
  if chargeCloud:
    #initialize the charges cloud:
    chargeCloud = []
    x_charge = []
    y_charge = []
    z_charge = []
    if secretZero:
      x_charge.append(0.00)
      y_charge.append(0.00)
      z_charge.append(0.00)
      chargeCloud.append(0.00)

  for iatom in atomList:
    if atomsAlso: thisAtom = mlab.points3d([iatom.x], [iatom.y], [iatom.z], color=symbolColour(iatom.symbol))	# Using jmols colours
    if chargeCloud:
      x_charge.append(iatom.x)
      y_charge.append(iatom.y)
      z_charge.append(iatom.z)
      chargeCloud.append(iatom.charge)
      
    if implicitBonds:
      if previousAtomDrawn:
        thisBondX = [iatom.x, prevAtom[0]] 
      	thisBondY = [iatom.y, prevAtom[1]] 
      	thisBondZ = [iatom.z, prevAtom[2]]       	
        thisAtomsBonds = mlab.plot3d(thisBondX, thisBondY, thisBondZ, tube_radius=0.1)
        prevAtom[0] = iatom.x
        prevAtom[1] = iatom.y
        prevAtom[2] = iatom.z
      else:
        previousAtomDrawn = True
        prevAtom =[]
        prevAtom.append(iatom.x)
        prevAtom.append(iatom.y)
        prevAtom.append(iatom.z)
        
    if numNN != 0 : 
      sortedNeighbours = neighSortByDistance(supercell, iatom)
      for (ineigh, mbeta) in sortedNeighbours[0:numNN]:
        if mbeta == 0:
          thisBondX = [iatom.x, ineigh.x] 
      	  thisBondY = [iatom.y, ineigh.y] 
      	  thisBondZ = [iatom.z, ineigh.z]       	
          if implicitBonds: thisAtomsBonds = mlab.plot3d(thisBondX, thisBondY, thisBondZ ,tube_radius=0.1)
          thisNeigh = mlab.points3d([ineigh.x], [ineigh.y], [ineigh.z], color=symbolColour(ineigh.symbol))
          if chargeCloud:
            x_charge.append(ineigh.x)
            y_charge.append(ineigh.y)
            z_charge.append(ineigh.z)
            chargeCloud.append(ineigh.charge)
        else:
          jatomX = ineigh.x + supercell.xl[mbeta][0]
          jatomY = ineigh.y + supercell.xl[mbeta][1]
          jatomZ = ineigh.z + supercell.xl[mbeta][2]
          thisBondX = [iatom.x, ineigh.x + supercell.xl[mbeta][0]]
          thisBondY = [iatom.y, ineigh.y + supercell.xl[mbeta][1]]
          thisBondZ = [iatom.z, ineigh.z + supercell.xl[mbeta][2]]
          if implicitBonds: thisAtomsBonds = mlab.plot3d(thisBondX, thisBondY, thisBondZ ,tube_radius=0.1)
          thisNeigh = mlab.points3d([jatomX], [jatomY], [jatomZ], color=symbolColour(ineigh.symbol))
          if chargeCloud:
            x_charge.append(jatomX)
            y_charge.append(jatomY)
            z_charge.append(jatomZ)
            chargeCloud.append(ineigh.charge)
          
    #if andIsoSurface:
      # To be implemented
  if chargeCloud:
    if secretZero: # because the charge "cloud" size is normalized, secretzero helps make it plotable when all atom charges are basically the same.
      x_charge[0] = sum(x_charge) / (len(x_charge) -1)
      y_charge[0] = sum(y_charge) / (len(y_charge) -1)
      z_charge[0] = sum(z_charge) / (len(z_charge) -1)
      
    chargeAtom = mlab.points3d(x_charge, y_charge, z_charge, chargeCloud, scale_mode='scalar', scale_factor = 1, opacity = 0.50, color=symbolColour(iatom.symbol))    #[iatom.charge]
  # Add labels as required... for some reason this wants to be at the end:
  for iatom in atomList:
    if chargeLabel: chargeLabels = mlab.text3d(iatom.x+1.0, iatom.y, iatom.z, str(iatom.charge), scale=(0.5, 0.5, 0.5))
    if nameLabel: nameLabels = mlab.text3d(iatom.x+1.0, iatom.y, iatom.z, str(iatom), scale=(0.5, 0.5, 0.5))
    if speciesLabel: speciesLabels = mlab.text3d(iatom.x+1.0, iatom.y, iatom.z, iatom.symbol,scale=(0.5, 0.5, 0.5))
    if numLabel: numLabels = mlab.text3d(iatom.x+1.0, iatom.y, iatom.z, str(iatom.atomNumber),scale=(0.5, 0.5, 0.5)) 
    if numNN != 0 : 
      for (ineigh, mbeta) in sortedNeighbours[0:numNN]:
        jatomX = ineigh.x + supercell.xl[mbeta][0]
        jatomY = ineigh.y + supercell.xl[mbeta][1]
        jatomZ = ineigh.z + supercell.xl[mbeta][2]
        if chargeLabel: chargeLabels = mlab.text3d(jatomX+1.0, jatomY, jatomZ, str(ineigh.charge), scale=(0.5, 0.5, 0.5))
        if nameLabel: nameLabels = mlab.text3d(jatomX+1.0, jatomY, jatomZ, str(ineigh), scale=(0.5, 0.5, 0.5))
        if speciesLabel: speciesLabels = mlab.text3d(jatomX+1.0, jatomY, jatomZ, ineigh.symbol,scale=(0.5, 0.5, 0.5))
        if numLabel: numLabels = mlab.text3d(jatomX+1.0, jatomY, jatomZ, str(ineigh.atomNumber),scale=(0.5, 0.5, 0.5)) 

        
        
        
  fig1.scene.disable_render = True

  #mlab.show()
  
  return fig1     

def fart():
  basFile ="DOSrun/answer.bas"
  lvsFile = "DOSrun/input.lvs"
  run(basFile, lvsFile)
  return

def run(basFile, lvsFile):
  elements = readInfo('./Fdata')
  BasList = ReadBas(basFile)
  LVSList = ReadLVS(lvsFile)
  
  NeighList = ReadNeigh()
  
  CHARGES = ReadCharges()
    
#  elements = readInfo('./Fdata')
  atomList = basList2Object(BasList, elements)
  chargesList2atoms(atomList, CHARGES)

  #Combine position and CHARGES data into a table
  xlm = xlmbeta(LVSList)
  mbeta = 0  
  neighList2atomObjects(atomList, NeighList, xlm)
  
  newAtom = atomList[130] 
  print "newAtom.charge", newAtom.charge
  print "newAtom.neutralCharge", newAtom.neutralCharge
  print "newAtom.chargeChange", newAtom.chargeChange  
  readDOS(atomList, directory = './DOSrun')
  print "len(newAtom.shells)", len(newAtom.shells)
  print "newAtom.shells[0].orbitals", newAtom.shells[0].orbitals
  print "len(newAtom.shells[0].orbitals", len(newAtom.shells[0].orbitals)
  print "len(newAtom.shells[1].orbitals", len(newAtom.shells[1].orbitals)
#  print "newAtom.shells[0].orbitals[0].l", newAtom.shells[0].orbitals[0].l
#  print "newAtom.orbitals.l", newAtom.orbitals.l
  print "newAtom.shells[0].lssh", newAtom.shells[0].lssh
  print "newAtom.shells[1].lssh", newAtom.shells[1].lssh
  print "elements['O'].shells", len(elements['O'].shells)
  print "elements['Fe'].shells", len(elements['Fe'].shells)
  print "elements['Al'].shells", len(elements['Al'].shells)
  print "elements['Cu'].shells", len(elements['Cu'].shells)
  print " "       
  
  dens_new = atom2DOS(atomList[9])
  dens_new.writeToFile('./Fart/testOut.dat')
  print dens_new.contributions[0]
  dens_new.plotme(-2.5, plotContributions = True)
#  dens_new2 = dens_new.offset(-2.5)
#  dens_new2.plotme(-2.5, plotContributions = True)
  densList = []
  #for i in range(8, 11):
#  densList.append(atom2DOS(atomList[80]))
#  densList.append(atom2DOS(atomList[180]))
#  densList.append(atom2DOS(atomList[280]))
  
#  plotDOSList(densList, plotHOMO = False)  
  ztot = computeZtot(atomList)
  print "Ztot = ", computeZtot(atomList)
  #read eigen and HOMO now required.
  
  eigens = readEigen("./DOSrun")
  print "len(eigens)", len(eigens)
  print "HOMO = ", eigens[int(ztot)/2]
  print "Eigen High = ", eigens[-1]
  print "Eigen Low =", eigens[0]
  
  densList.append(atom2DOS(atomList[80], eigens[int(ztot)/2]))
  densList.append(atom2DOS(atomList[180], eigens[int(ztot)/2]))
  densList.append(atom2DOS(atomList[280], eigens[int(ztot)/2]))
  plotDOSList(densList, plotOffset = False, xlim = [-1,1])
  plotDOSList(densList, plotOffset = True, xlim =[-1,1])
  plotDOSList(densList, plotOffset = False)
  
  # Lets see if we FINALLY make a supercell, it needs to have a directory, a basis file name, (norbitals), atomList, xl, eigenvals, ztot and a TOTAL DOS, to say the least.
  firstCellDirectory = "./DOSrun" 
#  firstCell = supercell(firstCellDirectory, "./DOSrun/input.bas", atomList, LVSList, xlm, eigens)
  
  readCell = directory2Supercell(elements, 'DOSrun')

  return 

def test(testCell, elements):
  #Combine position and CHARGES data into a table
  newAtom = testCell.atomList[130] 
  print "newAtom.charge", newAtom.charge
  print "newAtom.neutralCharge", newAtom.neutralCharge
  print "newAtom.chargeChange", newAtom.chargeChange  
  print "len(newAtom.shells)", len(newAtom.shells)
  print "newAtom.shells[0].orbitals", newAtom.shells[0].orbitals
  print "len(newAtom.shells[0].orbitals", len(newAtom.shells[0].orbitals)
  print "len(newAtom.shells[1].orbitals", len(newAtom.shells[1].orbitals)
#  print "newAtom.shells[0].orbitals[0].l", newAtom.shells[0].orbitals[0].l
#  print "newAtom.orbitals.l", newAtom.orbitals.l
  print "newAtom.shells[0].lssh", newAtom.shells[0].lssh
  print "newAtom.shells[1].lssh", newAtom.shells[1].lssh
  print "elements['O'].shells", len(elements['O'].shells)
  print "elements['Fe'].shells", len(elements['Fe'].shells)
  print "elements['Al'].shells", len(elements['Al'].shells)
  print "elements['Cu'].shells", len(elements['Cu'].shells)
  print " "       
  
  dens_new = atom2DOS(testCell.atomList[9])
  print dens_new.contributions.keys()
  dens_new.plotme(-2.5, plotContributions = True)
#  dens_new2 = dens_new.offset(-2.5)
#  dens_new2.plotme(-2.5, plotContributions = True)
  densList = []
  #for i in range(8, 11):
#  densList.append(atom2DOS(atomList[80]))
#  densList.append(atom2DOS(atomList[180]))
#  densList.append(atom2DOS(atomList[280]))
  
#  plotDOSList(densList, plotHOMO = False)  
  ztot = computeZtot(testCell.atomList)
  print "Ztot = ", computeZtot(testCell.atomList)

  #read eigen and HOMO now required.
  print "len(eigens)", len(testCell.eigens)
  print "HOMO = ", testCell.eigens[int(ztot)/2]
  print "Eigen High = ", testCell.eigens[-1]
  print "Eigen Low =", testCell.eigens[0]
  
  densList.append(atom2DOS(testCell.atomList[80], testCell.eigens[int(ztot)/2]))
  densList.append(atom2DOS(testCell.atomList[180], testCell.eigens[int(ztot)/2]))
  densList.append(atom2DOS(testCell.atomList[280], testCell.eigens[int(ztot)/2]))
  plotDOSList(densList, plotOffset = False, xlim = [-1,1])
  plotDOSList(densList, plotOffset = True, xlim =[-1,1])
  plotDOSList(densList, plotOffset = False)
  
  return 
  
def testDenswrite(cell, elements):
  dens_new = atom2DOS(cell.atomList[9])    
  dens_new.writeToFile('./Test/testOut.dat')

def quickStart():
  print "Reading in Fdata"
  elements = readInfo('./MultiCupratePaper/Al_Fdata/')
  cell1 = directory2Supercell(elements, './MultiCupratePaper/CuAl0.98Fe0.02/90/')
  return (elements, cell1)
  
def Oldmain():
  #Here, we will parse command-line arguments for the .bas and .lvs files. CHARGES is assumed to stay the same. 
  # Make a list of command line arguments, omitting the [0] element
  # which is the script itself.
  args = sys.argv[1:]

  if not args:
    print 'usage: basFile lvsFile'
    sys.exit(1)
  if len(args) > 2:
    print 'usage: basFile lvsFile'
    sys.exit(1)

  basFile = args[0]
  lvsFile = args[1]
  
  run(basFile, lvsFile)
  
  return   

def main():
  # From the command line, we expect an argument for the directory to load.
  # Make a list of command line arguments, omitting the [0] element
  # which is the script itself.
  args = sys.argv[1:]

  if not args:
    print 'usage: PyFireball [directoryName]'
    print 'assuming directory is .'
    directory = '.'
  elif len(args) > 1:
    print 'usage: PyFireball [directoryName]'
    sys.exit(1)
  else:
    directory = args[0]
  # Read FData  
  elements = readInfo('./Fdata')
  # Read supercell
  readCell = directory2Supercell(elements, directory)  
  
  # Run test
  #test(readCell, elements)
  #fart()
#  testDenswrite(readCell, elements)
#  fart()
  dens_spec13 = readCell.genDOSspec('Al')
  dens_spec13.writeToFile('./Test_dens_spec_13.dat')
  dens_spec13.plotme('AlDOS')
  readCell.plotPDOS(elements)
  return
def test3Plot():
  InDOS_combined= catDOSTots(InDOSes)
  GaDOS_combined= catDOSTots(GaDOSes)
  listOfDens = [InDOS_combined, GaDOS_combined, InDOS_combined]
  plot3DOS(listOfDens)
  return  
if __name__ == '__main__':
  #readinfoMain()
  main()