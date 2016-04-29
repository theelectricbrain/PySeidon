#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import seaborn
import pandas as pd
try:
    from osgeo import ogr
    from osgeo import osr
    havegdal=True
except ImportError:
    print 'Gdal is not installed, osgeo cannot be imported. Saving to shape files will not be possible.'
    havegdal=False
import zipfile
# Local import
from pyseidon_dvt.utilities.windrose import WindroseAxes
from pyseidon_dvt.utilities.interpolation_utils import *
#from miscellaneous import depth_at_FVCOM_element as depth_at_ind

# Kml header
###
kml_groundoverlay = '''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.0">
<Document>
<GroundOverlay>
  <name>__NAME__</name>
  <color>__COLOR__</color>
  <visibility>__VISIBILITY__</visibility>
  <Icon>
    <href>overlay.png</href>
  </Icon>
  <LatLonBox>
    <south>__SOUTH__</south>
    <north>__NORTH__</north>
    <west>__WEST__</west>
    <east>__EAST__</east>
  </LatLonBox>
</GroundOverlay>
<ScreenOverlay>
    <name>Legend</name>
    <Icon>
        <href>legend.png</href>
    </Icon>
    <overlayXY x="0" y="0" xunits="fraction" yunits="fraction"/>
    <screenXY x="0.015" y="0.075" xunits="fraction" yunits="fraction"/>
    <rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>
    <size x="0" y="0" xunits="pixels" yunits="pixels"/>
</ScreenOverlay>
</Document>
</kml>
'''
###

class PlotsFvcom:
    """
    **'Plots' subset of FVCOM class gathers plotting functions**
    """
    def __init__(self, variable, grid, debug):
        self._debug = debug
        #Back pointer
        setattr(self, '_var', variable)
        setattr(self, '_grid', grid)

        return

    def _def_fig(self):
        """Defines figure window"""
        self._fig = plt.figure(figsize=(18,10))
        plt.rc('font',size='22')


    def colormap_var(self, var, title=' ', cmin=[], cmax=[], cmap=[], bb = None,
                     degree=True, mesh=False, isoline = 'bathy',
                     dump=False, png=False, shapefile=False, kmz=False, debug=False, **kwargs):
        """
        2D xy colormap plot of any given variable and mesh.

        Input:
          - var = gridded variable, 1 D numpy array (nele or nnode)

        Options:
          - title = plot title, string
          - cmin = minimum limit colorbar
          - cmax = maximum limit colorbar
          - cmap = matplolib colormap
          - bb = bounding box, ex.:
                 bb = [minimun longitude, maximun longitude, minimun latitude, maximum latitude]
                 bb = [minimun x, maximun x, minimun y, maximum y]
          - mesh = True, with mesh; False, without mesh
          - degree = boolean, coordinates in degrees (True) or meters (False)
          - isloline = 'bathy': bathymetric isolines, 'var': variable isolines, 'none': no isolines
          - dump = boolean, dump profile data in csv file
          - png = boolean, save map as png
          - shapefile = boolean, save map as shapefile
          - kmz = boolean, save map as kmz
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
        debug = debug or self._debug
        if debug:
            print 'Plotting grid...'
        # Figure if var had nele or nnode dimensions
        if var.shape[0] == self._grid.nele:
            dim = self._grid.nele
        elif var.shape[0] == self._grid.nnode:
            dim = self._grid.nnode
        else:
            print "Var has the wrong dimension, var.shape[0]= Grid.nele or nnode"
            return

        # Bounding box nodes, elements and variable
        if degree:
            lon = self._grid.lon[:]
            lat = self._grid.lat[:]
            if debug:
                print "Computing bounding box..."
            if bb is None:
                if self._grid._ax == []:
                    self._grid._ax = [lon.min(), lon.max(),
                                     lat.min(), lat.max()]
                bb = self._grid._ax

            if not hasattr(self._grid, 'triangleLL'):        
                # Mesh triangle
                if debug:
                    print "Computing triangulation..."
                trinodes = self._grid.trinodes[:] 
                tri = Tri.Triangulation(lon, lat, triangles=trinodes)
                self._grid.triangleLL = tri
            else:
                tri = self._grid.triangleLL

        else:
            x = self._grid.x[:]
            y = self._grid.y[:]
            if debug:
                print "Computing bounding box..."
            if bb is None:
                bb = [x.min(), x.max(), y.min(), y.max()]

            if not hasattr(self._grid, 'triangleXY'):        
                # Mesh triangle
                if debug:
                    print "Computing triangulation..."
                trinodes = self._grid.trinodes[:] 
                tri = Tri.Triangulation(x, y, triangles=trinodes)
                self._grid.triangleXY = tri
            else:
                tri = self._grid.triangleXY

        #setting limits and levels of colormap
        if cmin==[]:
            if debug:
                print "Computing cmin..."
            cmin=var[:].min()
        if cmax==[]:
            if debug:
                print "Computing cmax..."
            cmax=var[:].max()
        step = (cmax-cmin) / 50.0
        levels=np.arange(cmin, (cmax+step), step)   # depth contours to plot

        #Figure window params
        self._def_fig()
        if degree:
            self._ax = self._fig.add_subplot(111,
                       aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
        else:
            self._ax = self._fig.add_subplot(111, aspect=1.0)

        #Plotting functions
        if debug:
            print "Computing colormap..."
        if cmap==[]:
            f = self._ax.tripcolor(tri, var[:],vmax=cmax,vmin=cmin,cmap=plt.cm.gist_earth)
        else:
            f = self._ax.tripcolor(tri, var[:],vmax=cmax,vmin=cmin,cmap=cmap) 
        if mesh:
            plt.triplot(tri, color='white', linewidth=0.5)

        #Label and axis parameters
        if degree:
            self._ax.set_ylabel('Latitude')
            self._ax.set_xlabel('Longitude')
        else:
            self._ax.set_ylabel('Distance (m)')
            self._ax.set_xlabel('Distance (m)')
        self._ax.patch.set_facecolor('0.5')
        cbar=self._fig.colorbar(f, ax=self._ax)
        cbar.set_label(title, rotation=-90,labelpad=30)
        scale = 1

        if degree:
            ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
            self._ax.xaxis.set_major_formatter(ticks)
            self._ax.yaxis.set_major_formatter(ticks)

        self._ax.set_xlim([bb[0],bb[1]])
        self._ax.set_ylim([bb[2],bb[3]])

        # Isolines
        if isoline == 'bathy':
            cs = self._ax.tricontour(tri, self._grid.h, colors='w', linewidths=1.0)
            plt.clabel(cs, fontsize=11, inline=1)
            plt.figtext(.12, .95, "Notes: white lines = bathymetric isolines", size='x-small')
        elif isoline == 'var':
            if var.shape[0] == self._grid.nele:
                vari = interp_linear_to_nodes(var, self._grid.xc, self._grid.yc, self._grid.x, self._grid.y)
            else:
                vari = var
            bounds=np.linspace(cmin,cmax,11)
            cs = self._ax.tricontour(tri, vari[:], vmin=cmin, vmax=cmax, colors='w', linewidths=1.0, levels=bounds)
            plt.clabel(cs, fontsize=11, inline=1)
            plt.figtext(.12, .95, "Notes: white lines = isolines", size='x-small')

        # Show plot
        self._ax.grid()
        self._fig.show()

        # Saving
        title = title.replace(" ", "_")
        title = title.replace("(", "_")
        title = title.replace(")", "_")
        title = title.replace("-", "_")
        title = title.replace("/", "_")
        title = title.replace(".", "_")
        savename=title.lower()
        if png:
            if kmz:
                self._fig.savefig("overlay.png", bbox_inches='tight', transparent=True)
            self._fig.savefig(savename+".png", bbox_inches='tight')

        if dump:
            if degree:
                self._dump_map_data_as_csv(var, self._grid.lonc, self._grid.latc,
                                           title=savename, varLabel='map',
                                           xLabel=' ', yLabel=' ', **kwargs)
            else:
                self._dump_map_data_as_csv(var, self._grid.xc, self._grid.yc, title=savename,
                                           varLabel='map', xLabel=' ', yLabel=' ', **kwargs)
        if debug or self._debug:
            print '...Passed'

        if shapefile and havegdal:
            if var.shape[0] == self._grid.nnode:
                if debug: print "Interpolating var..."
                var = interpN(var, self._grid.trinodes, self._grid.aw0, debug=debug)
            if degree:
                self._save_map_as_shapefile(var, self._grid.lon, self._grid.lat,
                                            title=savename, varLabel='map', debug=debug)
            else:
                self._save_map_as_shapefile(var, self._grid.x, self._grid.y,
                                            title=savename, varLabel='map', debug=debug)
        elif shapefile and not havegdal:
            print 'Shape file cannot be saved. Missing gdal.'

        if kmz:
            name = kwargs.pop('name', title)
            color = kwargs.pop('color', '9effffff')
            visibility = str( kwargs.pop('visibility', 1) )
            kmzfile = kwargs.pop('kmzfile', savename+'.kmz')
            pixels = kwargs.pop('pixels', 2048)  # pixels of the max. dimension
            units = kwargs.pop('units', '-')

            # Deformation dur to projection
            # if degree:
            #     geo_aspect = np.cos(lat.mean()*np.pi/180.0)
            #     xsize = lon.ptp()*geo_aspect
            #     ysize = lat.ptp()
            #     xmax = lon.max()
            #     ymax = lat.max()
            #     xmin = lon.min()
            #     ymin = lat.min()
            # else:
            #     geo_aspect = np.cos(self._grid.lat[:].mean()*np.pi/180.0)
            #     xsize = x.ptp()*geo_aspect
            #     ysize = y.ptp()
            #     xmax = x.max()
            #     ymax = y.max()
            #     xmin = x.min()
            #     ymin = y.min()
            geo_aspect = np.cos(self._grid.lat[:].mean()*np.pi/180.0)
            xsize = np.abs(bb[1] - bb[0]) * geo_aspect
            ysize = np.abs(bb[3] - bb[2])
            aspect = ysize/xsize
            if aspect > 1.0:
                figsize = (30.0/aspect, 30.0)
            else:
                figsize = (30.0, 30.0*aspect)

            plt.ioff()
            fig = plt.figure(figsize=figsize, facecolor=None, frameon=False, dpi=pixels // 5)
            # fig = figure(facecolor=None, frameon=False, dpi=pixels//10)
            ax = fig.add_axes([0, 0, 1, 1])
            if cmap == []:
                pc = ax.tripcolor(tri, var[:], vmax=cmax, vmin=cmin, cmap=plt.cm.gist_earth)
            else:
                pc = ax.tripcolor(tri, var[:], vmax=cmax, vmin=cmin, cmap=cmap)
            #ax.set_xlim(xmin, xmax)
            #ax.set_ylim(ymin, ymax)
            ax.set_xlim([bb[0], bb[1]])
            ax.set_ylim([bb[2], bb[3]])
            ax.set_axis_off()
            plt.savefig('overlay.png', dpi=200, transparent=True)

            # Write kmz
            fz = zipfile.ZipFile(kmzfile, 'w')
            fz.writestr(savename+'.kml', kml_groundoverlay.replace('__NAME__', name)\
                                                          .replace('__COLOR__', color)\
                                                          .replace('__VISIBILITY__', visibility)\
                                                          .replace('__SOUTH__', str(bb[2]))\
                                                          .replace('__NORTH__', str(bb[3]))\
                                                          .replace('__EAST__', str(bb[1]))\
                                                          .replace('__WEST__', str(bb[0])))
            fz.write('overlay.png')
            os.remove('overlay.png')

            # colorbar png
            fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
            ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
            cb = fig.colorbar(pc, cax=ax)
            cb.set_label(units, color='0.9')
            for lab in cb.ax.get_yticklabels():
                plt.setp(lab, 'color', '0.9')

            plt.savefig('legend.png', transparent=True)
            fz.write('legend.png')
            os.remove('legend.png')
            fz.close()

    def rose_diagram(self, direction, norm, png=False, title="rose_diagram"):

        """
        Plots rose diagram

        Inputs:
          - direction = 1D array
          - norm = 1D array

        Options:
          - title = plot title, string
          - png = boolean, saves rose diagram as png
        """
        #Convertion
        #TR: not quite sure here, seems to change from location to location
        #    express principal axis in compass
        direction = np.mod(90.0 - direction, 360.0)

        #Create new figure
        #fig = plt.figure(figsize=(18,10))
        #plt.rc('font',size='22')
        self._def_fig()      
        rect = [0.1, 0.1, 0.8, 0.8]
        ax = WindroseAxes(self._fig, rect)#, axisbg='w')
        self._fig.add_axes(ax)
        #Rose
        ax.bar(direction, norm , normed=True, opening=0.8, edgecolor='white')
        #adjust legend
        l = ax.legend(shadow=True, bbox_to_anchor=[-0.1, 0], loc='lower left')
        plt.setp(l.get_texts(), fontsize=10)
        plt.xlabel('Rose diagram in % of occurrences - Colormap of norms')
        self._fig.show()

        # Saving
        title = title.replace(" ", "_")
        title = title.replace("(", "_")
        title = title.replace(")", "_")
        title = title.replace("-", "_")
        title = title.replace("/", "_")
        title = title.replace(".", "_")
        savename=title.lower()
        if png:
            self._fig.savefig(savename+".png", bbox_inches='tight')

    def plot_xy(self, x, y, xerror=[], yerror=[],
                title='xy_plot', xLabel=' ', yLabel=' ',
                png=False, dump=False, **kwargs):
        """
        Simple X vs Y plot

        Inputs:
          - x = 1D array
          - y = 1D array

        Options:
          - xerror = error on 'x', 1D array
          - yerror = error on 'y', 1D array
          - title = plot title, string
          - xLabel = title of the x-axis, string
          - yLabel = title of the y-axis, string
          - png = boolean, saves map as png
          - dump = boolean, dump profile data in csv file
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
        #fig = plt.figure(figsize=(18,10))
        #plt.rc('font',size='22')
        self._def_fig()
        self._ax = self._fig.add_subplot(111)         
        self._ax.plot(x, y, label=title)
        scale = 1
        self._ax.set_ylabel(yLabel)
        self._ax.set_xlabel(xLabel)
        self._ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        self._ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        self._ax.grid(b=True, which='major', color='w', linewidth=1.5)
        self._ax.grid(b=True, which='minor', color='w', linewidth=0.5)
        if not yerror==[]:
            self._ax.fill_between(x, y-yerror, y+yerror,
            alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', antialiased=True)        
        if not xerror==[]:
            self._ax.fill_betweenx(y, x-xerror, x+xerror,
            alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', antialiased=True)
        if (not xerror==[]) or (not yerror==[]): 
            blue_patch = mpatches.Patch(color='#089FFF',
                         label='Standard deviation',alpha=0.2)
            plt.legend(handles=[blue_patch],loc=1, fontsize=12)
            #plt.legend([blue_patch],loc=1, fontsize=12)

        self._fig.show()
        # Saving
        title = title.replace(" ", "_")
        title = title.replace("(", "_")
        title = title.replace(")", "_")
        title = title.replace("-", "_")
        title = title.replace("/", "_")
        title = title.replace(".", "_")
        savename=title.lower().replace(" ","_")
        if png:
            self._fig.savefig(savename+".png", bbox_inches='tight')

        if dump: self._dump_profile_data_as_csv(x, y,xerror=xerror, yerror=yerror,
                                                title=savename, xLabel=xLabel,
                                                yLabel=yLabel, **kwargs)

    def Histogram(self, y, title='Histogram', xLabel=' ', yLabel=' ', bins=50,
                  png=False, dump=False, **kwargs):
        """
        Histogram plot

        Inputs:
          - bins = list of bin edges
          - y = 1D array

        Options:
          - title = plot title, string
          - xLabel = title of the x-axis, string
          - yLabel = title of the y-axis, string
          - bins = number of bins, integer
          - png = boolean, saves histogram as png
          - dump = boolean, dump profile data in csv file
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
        ## the histogram of the data
        #fig = plt.figure(figsize=(18,10))
        self._def_fig()
        self._ax = self._fig.add_subplot(111)        
        density, bins = np.histogram(y, bins=bins, normed=True, density=True)
        unity_density = density / density.sum()
        widths = bins[:-1] - bins[1:]
        # To plot correct percentages in the y axis 
        self._ax.bar(bins[1:], unity_density, width=widths)
        formatter = ticker.FuncFormatter(lambda v, pos: str(v * 100))
        self._ax.yaxis.set_major_formatter(formatter)
        self._ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        self._ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        self._ax.grid(b=True, which='major', color='w', linewidth=1.5)
        self._ax.grid(b=True, which='minor', color='w', linewidth=0.5)

        plt.ylabel(yLabel)
        plt.xlabel(xLabel)

        self._fig.show() 

        # Saving
        title = title.replace(" ", "_")
        title = title.replace("(", "_")
        title = title.replace(")", "_")
        title = title.replace("-", "_")
        title = title.replace("/", "_")
        title = title.replace(".", "_")
        savename=title.lower()
        if png:
            self._fig.savefig(savename+".png", bbox_inches='tight')
        if dump: self._dump_profile_data_as_csv(bins[1:], unity_density,
                                                title=savename, xLabel=xLabel,
                                                yLabel=yLabel, **kwargs)

    def add_points(self, x, y, label=' ', color='black'):
        """
        Adds scattered points (x,y) on current figure,
        where x and y are 1D arrays of the same lengths.

        Inputs:
          - x = float number or list of float numbers
          - y = float number or list of float numbers

        Options:
          - Label = a string
          - Color = a string, 'red', 'green', etc. or gray shades like '0.5'
        """
        plt.scatter(x, y, s=50, color=color)
        #TR : annotate does not work on my machine !?
        plt.annotate(label, xy=(x, y), xycoords='data', xytext=(-20, 20),
                     textcoords='offset points', ha='right',
                     arrowprops=dict(arrowstyle="->", shrinkA=0),
                     fontsize=12)

    def _dump_profile_data_as_csv(self, x, y, xerror=[], yerror=[],
                                 title=' ', xLabel=' ', yLabel=' ', **kwargs):
        """
        Dumps profile data in csv file

        Inputs:
          - x = 1D array
          - y = 1D array

        Options:
          - xerror = error on 'x', 1D array
          - yerror = error on 'y', 1D array
          - title = file name, string
          - xLabel = name of the x-data, string
          - yLabel = name of the y-data, string
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
        if title == ' ': title = 'dump_profile_data'
        filename=title + '.csv'
        if xLabel == ' ': xLabel = 'X'
        if yLabel == ' ': yLabel = 'Y'
        if not xerror == []:
            df = pd.DataFrame({xLabel:x[:], yLabel:y[:], 'error': xerror[:]})
        elif not yerror == []:
            df = pd.DataFrame({xLabel:x[:], yLabel:y[:], 'error': yerror[:]})
        else:
            df = pd.DataFrame({xLabel:x[:], yLabel:y[:]})
        df.to_csv(filename, encoding='utf-8', **kwargs)

    def _dump_map_data_as_csv(self, var, x, y, title=' ',
                              varLabel=' ', xLabel=' ', yLabel=' ', **kwargs):
        """
        Dumps map data in csv file

        Inputs:
          - var = gridded variable, 1 D numpy array (nele or nnode)
          - x = coordinates, 1 D numpy array (nele or nnode)
          - y = coordinates, 1 D numpy array (nele or nnode)

        Options:
          - title = file name, string
          - xLabel = name of the x-data, string
          - yLabel = name of the y-data, string
          - kwargs = keyword options associated with pandas.DataFrame.to_csv, such as:
                     sep, header, na_rep, index,...etc
                     Check doc. of "to_csv" for complete list of options
        """
        if title == ' ': title = 'dump_map_data'
        filename=title + '.csv'
        if varLabel == ' ': varLabel = 'Z'
        if xLabel == ' ': xLabel = 'X'
        if yLabel == ' ': yLabel = 'Y'
        df = pd.DataFrame({xLabel:x[:], yLabel:y[:], varLabel: var[:]})
        df.to_csv(filename, encoding='utf-8', **kwargs)

    def _save_map_as_shapefile(self, var, x, y, title=' ', varLabel=' ', debug=False):
        """
        Saves map as shapefile

        Inputs:
          - var = gridded variable, 1 D numpy array (nele or nnode)
          - x = coordinates, 1 D numpy array (nele or nnode)
          - y = coordinates, 1 D numpy array (nele or nnode)

        Options:
          - title = file name, string
          - kwargs = keyword options associated with ???
        """

        debug = debug or self._debug
        if debug:
            print 'Converting map to shapefile...'
        if title == ' ':
            title = 'save_map_data'
        else:  # reformat file name
            title = title.replace(" ", "_")
            title = title.replace("(", "_")
            title = title.replace(")", "_")
            title = title.replace("-", "_")
            title = title.replace("/", "_")
            title = title.replace(".", "_")

        filename=title + '.shp'

        # Projection
        #epsg_in=4326
        epsg_in = 3857  # Google Projection

        # give alternative file name is already exists
        driver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(filename):
            filename = filename[:-4] + "_bis.shp"

        shapeData = driver.CreateDataSource(filename)

        spatialRefi = osr.SpatialReference()
        spatialRefi.ImportFromEPSG(epsg_in)

        lyr = shapeData.CreateLayer("poly_layer", spatialRefi, ogr.wkbPolygon)

        #Features
        if varLabel==' ': valLabel = 'var'
        lyr.CreateField(ogr.FieldDefn(varLabel, ogr.OFTReal))

        if debug: print "Writing ESRI Shapefile %s..." % filename
        lon = x[:]
        lat = y[:]
        trinodes = self._grid.trinodes[:]

        if debug: print "Writing Node Array"
        cnt = 0
        for row in trinodes:
            val1 = -999
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for val in row:
                if val1 == -999:
                    val1 = val
                ring.AddPoint(lon[val], lat[val])
            #Add 1st point to close ring
            ring.AddPoint(lon[val1], lat[val1])

            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            #Now add field values from array
            feat = ogr.Feature(lyr.GetLayerDefn())
            feat.SetGeometry(poly)
            feat.SetField(varLabel, var[cnt])

            lyr.CreateFeature(feat)
            feat.Destroy()
            poly.Destroy()

            val1 = -999
            cnt += 1

        shapeData.Destroy()
        if debug: print "Finished writing Shapefile Mesh. [Total Nodes: %d]" % cnt           
