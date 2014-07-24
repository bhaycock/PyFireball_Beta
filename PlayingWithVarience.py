import PyFireball as PF
(elements, CuAlTest) = PF.quickStart()
PF.neighVdW(CuAlTest) or PF.readInNeighVdw(CuAlTest)
PF.clustering(CuAlTest, 'Fe', 3.0351)

seekSpecies = ['Fe', 'Cu']
PF.sortAllNeighboursByDistance(CuAlTest)
seekSpecies = ['Fe', 'Al']
PF.varience(CuAlTest, 'Cu', seekSpecies, 6 )
