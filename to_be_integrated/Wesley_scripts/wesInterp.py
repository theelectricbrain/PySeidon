from fvcomClass import FVCOM
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker

filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
test = FVCOM(filename)
trinodes = test.trinodes

# Wes Comment: ONly works with 1 point right now
X = -66.3385
Y = 44.277
index = test.closest_point([-66.3385],[44.277])[0]

newtri = Tri.Triangulation(test.lon[trinodes[index]],
                           test.lat[trinodes[index]], np.array([[0,1,2]]))

trif = newtri.get_trifinder()
trif.__call__(-66.3385, 44.277)

averEl = np.mean(test.el[:, trinodes[index]], axis=0)
inter = Tri.LinearTriInterpolator(newtri, averEl)
zi = inter(-66.3385, 44.277)



fig = plt.figure(figsize=(18,10))

ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(test.lat[trinodes[index]])*np.pi/180.0)))

plt.tricontourf(newtri, averEl, cmap=plt.cm.jet)
plt.triplot(newtri, 'ko-')

plt.plot(X,Y, 'ko-')

xy = (X,Y)
#ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='offset points')
ax.annotate('(%s, %s)' % xy, xy=xy, xycoords='data')

plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.title('Interpolating TriGrid')
cbar = plt.colorbar()
cbar.set_label('Time Averaged Elev', rotation=-90,labelpad=30)

scale = 1
ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
ax.xaxis.set_major_formatter(ticks)
ax.yaxis.set_major_formatter(ticks)


plt.show()

print '\nLon values at nodes'
print test.lon[trinodes[index]]

print '\nLat values at nodes'
print test.lat[trinodes[index]]

print '\nAverage Elevation at nodes'
print averEl

print '\nInterpolated value at point'
print zi
