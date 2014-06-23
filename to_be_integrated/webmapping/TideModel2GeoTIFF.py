#! /usr/bin/env python

import csv
import sys, getopt
import os
import subprocess as sub
import shlex

"""
This program is a simple CSV parser that reads a multi-column csv 
tidal model file and converts each column to a indivual GeoTIFF
EPSG 900913

Args: 	input csv file
	output prefix for GeoTIFF's
	region extent (wesn)
	region resolution
	ASCII polygon mask [optional] -- (see gmt grdmask)
	
Requirements: Requires GMT v. >5 with full resolution coastlines
and GDAL installed and in the PATH
http://gmt.soest.hawaii.edu/
http://www.gdal.org/	

TODO:	- option to get region from data
	- optional projection [default is EPSG 900913]
	- lots more?

"""



def main(argv):
	inputfile = ''
	outputfile = ''
	reg = ''
	resolution = ''
	poly_file = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:r:e:-p:",["ifile=","ofile="])
	except getopt.GetoptError:
		print 'TideModel2GeoTIFF.py -i <inputfile> -o <outfileprefix> -r <west/east/south/north> -e <resolution> -p <polygon_mask>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'TideModel2GeoTIFF.py -i <inputfile> -o <outfileprefix> -r <west/east/south/north> -e <resolution> -p <polygon_mask>'
			print 'where <west/east/south/north> must decimal degree geographic coordinates and '
			print '<resolution> is the output frid resolution in meters and'
			print '<polygon_mask> is a polygon mask file used by grdmask'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
		elif opt == '-r':
			reg = arg
		elif opt == '-e':
			resolution = arg
		elif opt == '-p':
			poly_file = arg
	if inputfile == '':
		print 'Error: no input file'
		print 'Usage: TideModel2GeoTIFF.py -i <inputfile> -o <outfileprefix> -r <west/east/south/north> -e <resolution> -p <polygon_mask>'
		sys.exit(2)
	if outputfile == '':
		print 'Error: no outputfile prefix'
		print 'Usage: TideModel2GeoTIFF.py -i <inputfile> -o <outfileprefix> -r <west/east/south/north> -e <resolution> -p <polygon_mask>'
		sys.exit(2)
	if reg == '':
		print 'Error: no grid region specified'
		print 'Usage: TideModel2GeoTIFF.py -i <inputfile> -o <outfileprefix> -r <west/east/south/north> -e <resolution> -p <polygon_mask>'
		sys.exit(2)
	if resolution == '':
		print 'Error: no grid resolution specified'
		print 'Usage: TideModel2GeoTIFF.py -i <inputfile> -o <outfileprefix> -r <west/east/south/north> -e <resolution> -p <polygon_mask>'
		sys.exit(2)
	
	region = '-R' + reg
	print 'Input file is ' + inputfile
	print 'Output file prefix is ' + outputfile
	print "Region = " + region
	print "Resolution = " + resolution
	if poly_file != '':
		print "Custom mask file is " + poly_file
	
	#Define column array 
	my_cols = ["MaxSpeed", "PowerDen", "PowerAsym", "FloodDev", "EbbDev", "DirDev"]
	cnt = 2
	for name in my_cols:
		print "Building " + name + " files...."
		#First parse csv file and output 3 column file for gridding
		with open(inputfile,"rb") as source:
			rdr= csv.reader( source )
			with open("tmp_file.csv" ,"wb") as result:
				wtr= csv.writer( result, delimiter=' ' )
				for r in rdr:
					if (r[cnt] == "NaN" or r[cnt] == "Inf"):
						print "Skipping NaN value"
					elif (name == 'PowerAsym' and float(r[cnt]) > 10000.0):
						#Skip really big values
						print "Skipping Power Asymmetry values (" + r[cnt] + ") greater than 20.0"
					else:
						wtr.writerow( (r[0], r[1], r[cnt]) )
		
		print "Finished writing CSV."
	
		#1st make chunk size set for gdal compat
		args_run = "gmtset IO_NC4_CHUNK_SIZE c"
		args = shlex.split(args_run)
		proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
		for line in proc.stderr:
			print(">>>> " + line)
	
		#Do the gridding
		print "Running sphinterpolate to build map..."
		args_run = "sphinterpolate tmp_file.csv " + region + " -I" + resolution + "e -Gmap1.cdf -Q1 -T -V"
		args = shlex.split(args_run)
		proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
		for line in proc.stderr:
			print(">>>> " + line)
		
		#reset negative values to zero
		#May only need to some fields??
		args_run = "grdclip map1.cdf -Sb0.0/0.0 -Gmap2.cdf"
		args = shlex.split(args_run)
		proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
		for line in proc.stderr:
			print(">>>> " + line)
	
		#Extract the mask	
		print "Extracting coastline mask, this may take a while..."
		args_run = "grdlandmask " + region + " -I" + resolution + "e -Gmask_map.cdf -Df -N1/NaN/NaN/NaN/NaN -V"
		args = shlex.split(args_run)
		proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
		for line in proc.stderr:
			print(">>>> " + line)
			
		if (poly_file == ''):
			#Apply the mask to the map
			args_run = "grdmath map2.cdf mask_map.cdf MUL = " + outputfile + "_" + name + ".cdf"
			args = shlex.split(args_run)
			proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
			for line in proc.stderr:
				print(">>>> " + line)
		else:
			#Apply the mask to the map
			args_run = "grdmath map2.cdf mask_map.cdf MUL = tmp_map1.cdf"
			args = shlex.split(args_run)
			proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
			for line in proc.stderr:
				print(">>>> " + line)
	
		if (poly_file != ''):
			#apply additional custom mask such as Block Rock
			print "Apply custom mask to map..."
			args_run = "grdmask " + poly_file + " " + region + " -I" + resolution + "e  -Gmask_map.cdf -N1/NaN/NaN -V"
			args = shlex.split(args_run)
			proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
			for line in proc.stderr:
				print(">>>> " + line)
			#Apply the mask to the map
			args_run = "grdmath tmp_map1.cdf mask_map.cdf MUL = " + outputfile + "_" + name + ".cdf"
			args = shlex.split(args_run)
			proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
			for line in proc.stderr:
				print(">>>> " + line)
					
		print "Cleaning up temp files"
		os.remove("tmp_file.csv")
		os.remove("mask_map.cdf")
		os.remove("map1.cdf")
		os.remove("map2.cdf")
		if (poly_file != ''):
			os.remove("tmp_map1.cdf")
	
		##Now make a GeoTiff
	
		#First add projection info
		print "Adding projection info"
		args_run = "gdal_translate -a_srs EPSG:4326 " + outputfile + "_" + name + ".cdf -a_nodata -999 grid2.tif"
		args = shlex.split(args_run)
		proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
		for line in proc.stderr:
			print(">>>> " + line)
		
		#Now re-project
		print "Re-Projecting to EPSG900913"
		args_run = "gdalwarp -t_srs EPSG:3857 -r bilinear -co \"TILED=YES\" grid2.tif " + outputfile + "_" + name + ".tif"
		args = shlex.split(args_run)
		proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
		for line in proc.stderr:
			print(">>>> " + line)
	
		#Add level of detail
		args_run = "gdaladdo " + outputfile + "_" + name + ".tif -r average  2 4 8 16"
		args = shlex.split(args_run)
		proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
		for line in proc.stderr:
			print(">>>> " + line)
		
		os.remove("grid2.tif")
	
		print "Finished Outputfile: " + outputfile + "_" + name + ".tif"
		print ""
		
		cnt += 1
	
	
if __name__ == "__main__":
   main(sys.argv[1:])



