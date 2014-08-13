from scipy.spatial import Delaunay
from fvcomClass import FVCOM
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import tri


filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'

ind = [-66.3419, -66.3324, 44.2755, 44.2815]

test = FVCOM(filename, ax=ind)


z = np.vstack((test.lon[test.region_e], test.lat[test.region_e])).T
points = map(tuple,z)
delTri = Delaunay(points)
delTri.simplices
delTri.simplices.shape

npoints = np.asarray(points)
plt.triplot(npoints[:,0], npoints[:,1], delTri.simplices.copy())
plt.show()

z = np.vstack((test.lon, test.lat)).T
points = map(tuple,z)
delTri = Delaunay(points)
delTri.simplices
delTri.simplices.shape

npoints = np.asarray(points)
plt.triplot(npoints[:,0], npoints[:,1], delTri.simplices.copy())
plt.show()


Tri = tri.Triangulation(test.lon, test.lat, test.trinodes)
Tri2 = tri.Triangulation(test.lon, test.lat, test.trinodes[test.region_e])

plt.triplot(Tri)
plt.triplot(Tri2, 'r--')
plt.show()

Tri3 = tri.Triangulation(test.lon[test.region_n], test.lat[test.region_n])
plt.triplot(Tri3)
plt.show()

z = np.vstack((test.lon[test.region_n], test.lat[test.region_n])).T
points = map(tuple,z)
delTri = Delaunay(points)
delTri.simplices
delTri.simplices.shape

npoints = np.asarray(points)
plt.triplot(npoints[:,0], npoints[:,1], delTri.simplices.copy())
plt.show()

