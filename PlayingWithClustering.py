import PyFireball as PF
(elements, CuAl) = PF.quickStart()
PF.neighVdW(CuAl) or PF.readInNeighVdw(CuAl)
PF.clustering(CuAl, 'Fe', 3.0351)

sortList = PF.neighSortByDistance(CuAl, CuAl.atomList[0])

for neighAtom in SortDistList[1:30]:
  jatom = neighAtom[0] + CuAl.xl[neighAtom[1]]
  print iatom - jatom, jatom