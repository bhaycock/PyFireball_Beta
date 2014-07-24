# FUCK! Can't I just do this with proper clustering analysis via scikit?
# I'd still have to generate the list of atom pairs distances, wouldn't I?
# and then, what would I do? I'd pop the list (by number of pairs) into a clustering algorithm
# giving me... well, in the case of, e.g.,  6Fe... nothing... unless I get a fecking average.

# BUUUUTTT, I can do this with an evalution of the variance- clustering analysis on the 
# number of atom pairs and vicinity of Fe to that pair, I think.

# Hmm.

import sys
import math
import copy
import matplotlib.pyplot as plt 
import PyFireball as PF

# I require rinput, supercell, species of interest.
# That appears to be all.
def clustering (supercell, speciesOfInterest, rinput, elements):		# Not sure if I need elements here, but anyways, and I want rinput out somehow.
  # Loop through all atoms.
  
  # If one is a speciesOfInterest: Loop through that atoms neighbors
  speciesOfInterestList = atom for atom in supercell.atomList where atom.symbol == speciesOfInterest
  maxClusterSize = len(speciesOfInterestList) - 1
  for iatom in speciesOfInterestList:
    neighbourOfInterestList = jatom for jatom in iatom.neighbourList where jatom.symbol == speciesOfInterest
    for jatom in neighboursOfInterestList:
       
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
      msList.append((Rinput - (iatom - jatom))^2)
#end loops over atoms (do I need that?)
  
# final ms is the first max_cluster_size of atoms in the sorted msList (sort low-high), sum these, divide by max_cluster_size, then SQRT
  msList.sort() 
  clustering = sqrt((sum(msList[0:maxClusterSize +1]))/maxClusterSize)
  
 return clustering


