#! /usr/bin/env python
 
# author: Pierre Poulain
# contributors: Justine Guegan, Edithe Selwa
# last update: 20110414
#
# this Python script compute principal axes from a PDB file
# it also produces a .pml script for a nice rendering with Pymol
 
#==========================================================================
# import required modules
#==========================================================================
import sys
import os.path
import numpy
 
#==========================================================================
# define data or hard-coded parameters
#==========================================================================
# scale factor to enhance the length of axis in Pymol
scale_factor = 20
 
#==========================================================================
# define functions
#==========================================================================
 
#==========================================================================
def read_pdb_xyz(pdb_name):
    """
read xyz from pdb
return:
[[x1 y1 z1]
[x2 y2 z2]
[.. .. ..]
[xn yn zn]]
   """
    xyz = []
    pdb_file = open(pdb_name, 'r')
    for line in pdb_file:
        if line.startswith("ATOM"):
            # extract x, y, z coordinates for carbon alpha atoms
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            if line[12:16].strip() == "CA":
                                xyz.append([x, y, z])
    pdb_file.close()
    return xyz
 
#==========================================================================
# start program
#==========================================================================
 
# check if argument is there
if len(sys.argv) == 2:
        pdb_name = sys.argv[1]
else:
    message = """
ERROR: missing pdb filename as argument
usage: %s file.pdb""" %(sys.argv[0])
    sys.exit(message)
 
# check if argument is an existing file
if not os.path.exists(pdb_name):
    sys.exit("ERROR: file %s does not seem to exist" %(pdb_name))
 
 
#--------------------------------------------------------------------------
# compute principal axes
#--------------------------------------------------------------------------
# read pdb
xyz = read_pdb_xyz(pdb_name)
print "%d CA atomes found if %s" %(len(xyz), pdb_name)
 
#create coordinates array
coord = numpy.array(xyz, float)
 
# compute geometric center
center = numpy.mean(coord, 0)
print "geometric center coordinates:\n", center
 
# center with geometric center
coord = coord - center
 
# compute principal axis matrix
inertia = numpy.dot(coord.transpose(), coord)
e_values, e_vectors = numpy.linalg.eig(inertia)
# warning eigen values are not necessary ordered!
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html
print "(unordered) eigen values:"
print e_values
print "(unordered) eigen vectors:"
print e_vectors
 
#--------------------------------------------------------------------------
# order eigen values (and eigen vectors)
#
# axis1 is the principal axis with the biggest eigen value (eval1)
# axis2 is the principal axis with the second biggest eigen value (eval2)
# axis3 is the principal axis with the smallest eigen value (eval3)
#--------------------------------------------------------------------------
for i in xrange(len(e_values)):
        # find biggest eigen value
        if e_values[i] == max(e_values):
                eval1 = e_values[i]
                axis1 = e_vectors[:,i]
        # find smallest eigen value
        elif e_values[i] == min(e_values):
                eval3 = e_values[i]
                axis3 = e_vectors[:,i]
        # middle eigen value
        else:
                eval2 = e_values[i]
                axis2 = e_vectors[:,i]
 
print "inertia axis are now ordered !"

print "the first principal axis is in red"
print "coordinates: ", axis1
print "eigen value: ", eval1
print
print "the second principal axis is in green"
print "coordinates:", axis2
print "eigen value:", eval2
print
print "the third principal axis is in blue"
print "coordinates:", axis3
print "eigen value:", eval3
