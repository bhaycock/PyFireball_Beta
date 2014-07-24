import PyFireball as PF
print "Reading in Fdata"
elements = PF.readInfo('./MultiCupratePaper/Al_Fdata/')
elements.update(PF.readInfo('./MultiCupratePaper/Ga_Fdata/'))
elements.update(PF.readInfo('./MultiCupratePaper/In_Fdata/'))

print "Reading in Supercells"
CuAl098Fe002 = PF.directory2Supercell(elements, './MultiCupratePaper/CuAl0.98Fe0.02/90/')
CuAl099Fe001 = PF.directory2Supercell(elements, './MultiCupratePaper/CuAl0.99Fe0.01/DOSrun/')
CuAl100Fe000 = PF.directory2Supercell(elements, './MultiCupratePaper/CuAl1.00Fe0.00/CuAl1.00Fe0.00_NewDOS_for_Even_energies/001/')
CuGa098Fe002 = PF.directory2Supercell(elements, './MultiCupratePaper/CuGa0.98Fe0.02/57/')
CuGa099Fe001 = PF.directory2Supercell(elements, './MultiCupratePaper/CuGa0.99Fe0.01/')
CuGa100Fe000 = PF.directory2Supercell(elements, './MultiCupratePaper/CuGa1.00Fe0.00/')
CuIn098Fe002 = PF.directory2Supercell(elements, './MultiCupratePaper/CuIn0.98Fe0.02/3/')
CuIn099Fe001 = PF.directory2Supercell(elements, './MultiCupratePaper/CuIn0.99Fe0.01/001/')
CuIn100Fe000 = PF.directory2Supercell(elements, './MultiCupratePaper/CuIn1.00Fe0.00/001/')

print "Creating HOMO list"
InMultiHOMO = {'CuIn$_{0.99}$Fe$_{0.01}$O$_2$' : CuIn099Fe001.HOMO, 'CuIn$_{1.00}$Fe$_{0.00}$O$_2$' : CuIn100Fe000.HOMO, 'CuIn$_{0.98}$Fe$_{0.02}$O$_2$' : CuIn098Fe002.HOMO}
GaMultiHOMO = {'CuGa$_{0.99}$Fe$_{0.01}$O$_2$' : CuGa099Fe001.HOMO, 'CuGa$_{1.00}$Fe$_{0.00}$O$_2$' : CuGa100Fe000.HOMO, 'CuGa$_{0.98}$Fe$_{0.02}$O$_2$' : CuGa098Fe002.HOMO}
AlMultiHOMO = {'CuAl$_{0.99}$Fe$_{0.01}$O$_2$' : CuAl099Fe001.HOMO, 'CuAl$_{1.00}$Fe$_{0.00}$O$_2$' : CuAl100Fe000.HOMO, 'CuAl$_{0.98}$Fe$_{0.02}$O$_2$' : CuAl098Fe002.HOMO}
multiHOMOList = [AlMultiHOMO, GaMultiHOMO, InMultiHOMO]


print "Generating TOTAL DOS image"

InDOSes = {'CuIn$_{0.98}$Fe$_{0.02}$O$_2$' : CuIn098Fe002.genDOS(elements), 'CuIn$_{0.99}$Fe$_{0.01}$O$_2$' :CuIn099Fe001.genDOS(elements), 'CuIn$_{1.00}$Fe$_{0.00}$O$_2$' :CuIn100Fe000.genDOS(elements)}
AlDOSes = {'CuAl$_{0.98}$Fe$_{0.02}$O$_2$' : CuAl098Fe002.genDOS(elements), 'CuAl$_{0.99}$Fe$_{0.01}$O$_2$' :CuAl099Fe001.genDOS(elements), 'CuAl$_{1.00}$Fe$_{0.00}$O$_2$' :CuAl100Fe000.genDOS(elements)}
GaDOSes = {'CuGa$_{0.98}$Fe$_{0.02}$O$_2$' : CuGa098Fe002.genDOS(elements), 'CuGa$_{0.99}$Fe$_{0.01}$O$_2$' :CuGa099Fe001.genDOS(elements), 'CuGa$_{1.00}$Fe$_{0.00}$O$_2$' :CuGa100Fe000.genDOS(elements)}

InDOS_combined= PF.catDOSTots(InDOSes)
GaDOS_combined = PF.catDOSTots(GaDOSes)
AlDOS_combined = PF.catDOSTots(AlDOSes)

listOfDens = [AlDOS_combined, GaDOS_combined, InDOS_combined]

PF.plot3DOS(listOfDens, multiHOMO = multiHOMOList, subLegend = True, contriblw =2, xlim = [-4, 4], name = 'TotalDOS_Fig5')

# Figure 6 is the Cu-only contributions:
print "Generating Cu DOS image"

CuOnlyAlDOSes = {'CuAl$_{0.98}$Fe$_{0.02}$O$_2$' : CuAl098Fe002.genDOSspec('Cu'), 'CuAl$_{0.99}$Fe$_{0.01}$O$_2$' :CuAl099Fe001.genDOSspec('Cu'), 'CuAl$_{1.00}$Fe$_{0.00}$O$_2$' :CuAl100Fe000.genDOSspec('Cu')}
CuOnlyInDOSes = {'CuIn$_{0.98}$Fe$_{0.02}$O$_2$' : CuIn098Fe002.genDOSspec('Cu'), 'CuIn$_{0.99}$Fe$_{0.01}$O$_2$' :CuIn099Fe001.genDOSspec('Cu'), 'CuIn$_{1.00}$Fe$_{0.00}$O$_2$' :CuIn100Fe000.genDOSspec('Cu')}
CuOnlyGaDOSes = {'CuGa$_{0.98}$Fe$_{0.02}$O$_2$' : CuGa098Fe002.genDOSspec('Cu'), 'CuGa$_{0.99}$Fe$_{0.01}$O$_2$' :CuGa099Fe001.genDOSspec('Cu'), 'CuGa$_{1.00}$Fe$_{0.00}$O$_2$' :CuGa100Fe000.genDOSspec('Cu')}

CuOnlyInDOS_combined= PF.catDOSTots(CuOnlyInDOSes)
CuOnlyGaDOS_combined = PF.catDOSTots(CuOnlyGaDOSes)
CuOnlyAlDOS_combined = PF.catDOSTots(CuOnlyAlDOSes)


listOfDens = [CuOnlyAlDOS_combined, CuOnlyGaDOS_combined, CuOnlyInDOS_combined]

PF.plot3DOS(listOfDens, multiHOMO = multiHOMOList, subLegend = True, contriblw =2, xlim = [-4, 4], name = 'CuDOS_Fig6')


# Figure 8 is the O-only contributions:
print "Generating O DOS image"

OOnlyAlDOSes = {'CuAl$_{0.98}$Fe$_{0.02}$O$_2$' : CuAl098Fe002.genDOSspec('O'), 'CuAl$_{0.99}$Fe$_{0.01}$O$_2$' :CuAl099Fe001.genDOSspec('O'), 'CuAl$_{1.00}$Fe$_{0.00}$O$_2$' :CuAl100Fe000.genDOSspec('O')}
OOnlyInDOSes = {'CuIn$_{0.98}$Fe$_{0.02}$O$_2$' : CuIn098Fe002.genDOSspec('O'), 'CuIn$_{0.99}$Fe$_{0.01}$O$_2$' :CuIn099Fe001.genDOSspec('O'), 'CuIn$_{1.00}$Fe$_{0.00}$O$_2$' :CuIn100Fe000.genDOSspec('O')}
OOnlyGaDOSes = {'CuGa$_{0.98}$Fe$_{0.02}$O$_2$' : CuGa098Fe002.genDOSspec('O'), 'CuGa$_{0.99}$Fe$_{0.01}$O$_2$' :CuGa099Fe001.genDOSspec('O'), 'CuGa$_{1.00}$Fe$_{0.00}$O$_2$' :CuGa100Fe000.genDOSspec('O')}

OOnlyInDOS_combined= PF.catDOSTots(OOnlyInDOSes)
OOnlyGaDOS_combined = PF.catDOSTots(OOnlyGaDOSes)
OOnlyAlDOS_combined = PF.catDOSTots(OOnlyAlDOSes)


listOfDens = [OOnlyAlDOS_combined, OOnlyGaDOS_combined, OOnlyInDOS_combined]

PF.plot3DOS(listOfDens, multiHOMO = multiHOMOList, subLegend = True, contriblw =2, xlim = [-4, 4], name = 'ODOS_Fig8')

# Figure 9 is the Bsite-only contributions:
print "Generating BSite DOS image"

BsiteOnlyAlDOSes = {'CuAl$_{0.98}$Fe$_{0.02}$O$_2$' : CuAl098Fe002.genDOSspec('Al'), 'CuAl$_{0.99}$Fe$_{0.01}$O$_2$' :CuAl099Fe001.genDOSspec('Al'), 'CuAl$_{1.00}$Fe$_{0.00}$O$_2$' :CuAl100Fe000.genDOSspec('Al')}
BsiteOnlyInDOSes = {'CuIn$_{0.98}$Fe$_{0.02}$O$_2$' : CuIn098Fe002.genDOSspec('In'), 'CuIn$_{0.99}$Fe$_{0.01}$O$_2$' :CuIn099Fe001.genDOSspec('In'), 'CuIn$_{1.00}$Fe$_{0.00}$O$_2$' :CuIn100Fe000.genDOSspec('In')}
BsiteOnlyGaDOSes = {'CuGa$_{0.98}$Fe$_{0.02}$O$_2$' : CuGa098Fe002.genDOSspec('Ga'), 'CuGa$_{0.99}$Fe$_{0.01}$O$_2$' :CuGa099Fe001.genDOSspec('Ga'), 'CuGa$_{1.00}$Fe$_{0.00}$O$_2$' :CuGa100Fe000.genDOSspec('Ga')}

BsiteOnlyInDOS_combined= PF.catDOSTots(BsiteOnlyInDOSes)
BsiteOnlyGaDOS_combined = PF.catDOSTots(BsiteOnlyGaDOSes)
BsiteOnlyAlDOS_combined = PF.catDOSTots(BsiteOnlyAlDOSes)


listOfDens = [BsiteOnlyAlDOS_combined, BsiteOnlyGaDOS_combined, BsiteOnlyInDOS_combined]

PF.plot3DOS(listOfDens, multiHOMO = multiHOMOList, subLegend = True, contriblw =2, xlim = [-4, 4], name = 'BsiteDOS_Fig9')

# Figure 7 is the Fe-only contributions:
print "Generating Fe DOS image"

FeOnlyAlDOSes = {'CuAl$_{0.98}$Fe$_{0.02}$O$_2$' : CuAl098Fe002.genDOSspec('Fe'), 'CuAl$_{0.99}$Fe$_{0.01}$O$_2$' :CuAl099Fe001.genDOSspec('Fe')}
FeOnlyInDOSes = {'CuIn$_{0.98}$Fe$_{0.02}$O$_2$' : CuIn098Fe002.genDOSspec('Fe'), 'CuIn$_{0.99}$Fe$_{0.01}$O$_2$' :CuIn099Fe001.genDOSspec('Fe')}
FeOnlyGaDOSes = {'CuGa$_{0.98}$Fe$_{0.02}$O$_2$' : CuGa098Fe002.genDOSspec('Fe'), 'CuGa$_{0.99}$Fe$_{0.01}$O$_2$' :CuGa099Fe001.genDOSspec('Fe')}

FeOnlyInDOS_combined= PF.catDOSTots(FeOnlyInDOSes)
FeOnlyGaDOS_combined = PF.catDOSTots(FeOnlyGaDOSes)
FeOnlyAlDOS_combined = PF.catDOSTots(FeOnlyAlDOSes)

listOfDens = [FeOnlyAlDOS_combined, FeOnlyGaDOS_combined, FeOnlyInDOS_combined]

PF.plot3DOS(listOfDens, multiHOMO = multiHOMOList, subLegend = True, contriblw =2, xlim = [-4, 4], name = 'FeDOS_Fig7')

for key, value in GaMultiHOMO.iteritems():
  print key, value    

print "BAM!"
