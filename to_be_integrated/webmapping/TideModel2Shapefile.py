#! /usr/bin/env python

import csv
import sys, getopt
import subprocess as sub
import shlex

def main(argv):
	inputfile = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print 'TideModel2Shapefile.py -i <inputfile> -o <outfileprefix>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'TideModel2Shapefile.py -i <inputfile> -o <outfileprefix>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
	if inputfile == '':
		print 'Error: no input file'
		print 'Usage: TideModel2Shapefile.py -i <inputfile> -o <outfileprefix>'
		sys.exit(2)
	if outputfile == '':
		print 'Error: no outputfile prefix'
		print 'Usage: TideModel2Shapefile.py -i <inputfile> -o <outfileprefix>'
		sys.exit(2)	
	print 'Input file is ' + inputfile
	print 'Output file prefix is ' + outputfile
	
	file_read = csv.DictReader(open(inputfile, 'rb'), ["Longitude", "Latitude", "MaxSpeed", "PowerDen", "PowerAsym", "FloodDev", "EbbDev", "DirDev"], delimiter=',')
	file_write = csv.DictWriter(open(outputfile+'.csv', 'wb'), ["Longitude", "Latitude", "MaxSpeed", "PowerDen", "PowerAsym", "FloodDev", "EbbDev", "DirDev"], delimiter=',')
	file_write.writeheader()
	file_write.writerows(file_read)
	
	print "Finished writing CSV with header"
	
	txt_file = open(outputfile+'.vrt', "w")
	txt_file.write("<OGRVRTDataSource>\n\t")
	txt_file.write("<OGRVRTLayer name=\""+outputfile+"\">\n\t\t")
	txt_file.write("<SrcDataSource>"+outputfile+".csv</SrcDataSource>\n\t\t")
	txt_file.write("<GeometryType>wkbPoint</GeometryType>\n\t\t")
	txt_file.write("<LayerSRS>WGS84</LayerSRS>\n\t\t")
	txt_file.write("<GeometryField encoding=\"PointFromColumns\" x=\"Longitude\" y=\"Latitude\"/>\n\t")
	txt_file.write("</OGRVRTLayer>\n")
	txt_file.write("</OGRVRTDataSource>")
	txt_file.close()
	
	print "Finished writing OGR VRT File"
	
	args_run = "ogr2ogr -t_srs EPSG:900913 " + outputfile + "-EPSG900913.shp " + outputfile + ".vrt"
	args = shlex.split(args_run)
	proc = sub.Popen(args, stderr=sub.PIPE, shell=False)
	for line in proc.stderr:
		print(">>>> " + line)
		
	print "Finished creating new EPSG900913 Shapefile: " + outputfile + "-EPSG900913.shp "
	
if __name__ == "__main__":
   main(sys.argv[1:])



