PyFireball_Beta
===============

Collection of Tools for handling High-Throughput calculations producted by the FIREBALL ab-initio MD package, written in Python.

PyFireball-Beta
===============

Collection of tools for processing High-Throughput FIREBALL calculations.

The PyFIREBALL toolbox is a collection of tools for analysis and generation of FIREBALL jobs. Principally, the toolbox exists to assist in the analysis of High Throughput studies.
PyFIREBALL is as “pythonic” as possible and therefore has the following features:
•	an atom class, which contains all information about a specific atom in any cell.
o	therefore every atom object is unique and can interact with other atom objects:
•	e.g. atom1 – atom2 = distance between the atoms.
•	any string function (such as print) on the atom tells me about the atom’s specific features, including it’s position, species, etc.
•	atom.Neighbors(‘species’, N) returns a list of Nth neighbours of a specific species… therefore an atom object can tell you about how close it is to an impurity.
•	…. everything else FIREBALL keeps about an atom is in the atom object, for every single atom.
•	a supercell class, which acts as a container for atoms. It also contains the supercell specific data, and a bunch of methods such as supercell.plotDOS(), which will carry out all the calculations and generate a rather pretty DOS plot immediately.
•	Ability to directly edit the fireball.in (and related files), so I can observe 1000+ supercells, and then send 10 to carry out DOS or Optical calcuations.
•	Draw supercells, atom groups, charge densities, etc. directly from the command line.
•	Directly probe neighbours-of-neighbours, etc.
The philosophy behind PyFIREBALL Beta is that you can load the module into iPython Notebook, IDLE, etc, and immediately begin working on your dataset. Supercell objects can be created with the readDirectory() method, and Fdata can be directly read in via the readInfo() method. Then you can quickly process, analyse and report on your findings in an ETL fashion. The roadmap for PyFIREBALL includes further visualization abailities and xsf file handling, allowing for the multi step processes in generation and analysis of HT calcuations to be simplified with your own scripts or within an interactive python console.
