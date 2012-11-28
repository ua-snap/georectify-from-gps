#!/usr/bin/env python

###############################################################################
# importScout.py 
#
# Project:  Aeryon Scout data processing
# Purpose:  Reads images from the Aeryon scout system and geocodes them to UTM projection 
# Author:   Scott Arko
#
###############################################################################
# Copyright (c) 2012, Scott Arko 
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
# 
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
###############################################################################
# NOTES:
#
#  This software was designed to read a variety of formats that can come from teh Scout.  
#  however, only the tif version currently works.  The dng import should work, however, both GDAL and PIL
#  read the dng files at half resolution.  Once this is resolved, the workflow for tif and dng are 
#  identical. 
#
#
#


# Import necessary modules

# Required modules
import os
from osgeo import gdal
from PIL import Image
from PIL.ExifTags import TAGS
import sys
import math
import numpy as np
from pyproj import transform,Proj

# Optional modules

try:
    import matplotlib.pylab as plt
    from matplotlib import pyplot
    import mpl_toolkits.mplot3d.axes3d as p3
    mpl = True
except ImportError:
    mpl = False



# Various functions

def Usage():
    print "############################"
    print "Usage:"
    print "    importScout.py [-d] [-align] filelist"
    print ""
    print " -d -> Debug mode that will display an image of the data being read for verification.  Requires matplotlib python binding installed"
    print " -align -> Performs FFT-based alignment of adjacent images.  NOT IMPLEMENTED AT CURRENT TIME"
    print ""
    print "This program can currently only read tif images that have been converted from "
    print "Scout dng files along with their corresponding .nfo files.  The .tif and .nfo file"
    print "must have the same name"
    print "############################"

def getExif(fn):
    ret = {}
    i = Image.open(fn)
    info = i._getexif()
    for tag, value in info.items():
        decoded = TAGS.get(tag, tag)
        ret[decoded] = value
    return ret

def getGPSData(exif):
    # Read the exif data returned by getExif and parse out GPS lat, lon, alt and bearing
    lat = [float(x)/float(y) for x, y in exif['GPSInfo'][2]]
    latref = exif['GPSInfo'][1]
    lon = [float(x)/float(y) for x, y in exif['GPSInfo'][4]]
    lonref = exif['GPSInfo'][3]
    alt = float(exif['GPSInfo'][6][0])/float(exif['GPSInfo'][6][1])
    bearing = float(exif['GPSInfo'][24][0])/float(exif['GPSInfo'][24][1])
   
    lat = lat[0] + lat[1]/60 + lat[2]/3600
    lon = lon[0] + lon[1]/60 + lon[2]/3600
    if latref == 'S':
        lat = -lat
    if lonref == 'W':
        lon = -lon
       
    return lat,lon,alt,bearing

def calcImageParams(file,lat,lon,alt,bearing):
    # This function calculates the necessary image parameters to generate a 
    # standard world file.  Requires lat,lon,alt,bearing provided by getGPSData
    fh = gdal.Open(file)
    xsize = fh.RasterXSize
    ysize = fh.RasterYSize
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init='epsg:32611')
    (east,north) = transform(inProj,outProj,lon,lat)
    xres = 0.0064
    yres = 0.0064
    xUL = east-(xres*xsize/2)
    yUL = north + (yres*ysize/2)
    bearingRad = bearing*3.14159/180
    tanb = math.tan(bearingRad)
    line2 = -1*tanb*xres
    line3 = -1*tanb*yres

def createTransform(lat,lon,alt,bearing):
    # This function calculates the necessary image parameters to generate a 
    # standard world file.  Requires lat,lon,alt,bearing provided by getGPSData
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init='epsg:32611')
    print inProj
    lat = lat*180.0/np.pi
    lon = lon*180.0/np.pi
    east,north = transform(inProj,outProj,lon,lat)
    print 'Center Point: '+str(east)+' '+str(north)
    (xres,yres) = calcResolution(alt)
    # hard coded to scout parameters
    xsize = 2592
    ysize = 1944
    xd = xres*xsize/2
    yd = yres*ysize/2
    offset = math.atan(xd/yd)
    radius = math.sqrt(xd**2 + yd**2)
    bearingRad = bearing*3.14159/180
    xUL = east + radius*math.sin(bearingRad + offset)
    yUL = north + radius*math.cos(bearingRad + offset)
    tanb = math.tan(bearingRad)
    if bearing <= 180:
        xrot = -1*tanb*xres
        yrot = 1*tanb*yres
    else: 
        xrot = -tanb*xres
        yrot = 1*tanb*yres
    if bearing > 270 or bearing < 90:
        xres = xres*np.cos(bearingRad)
        yres = yres*np.cos(bearingRad)
    else:
        xres = xres*np.cos(bearingRad-np.pi)
        yres = yres*np.cos(bearingRad-np.pi)
    return xUL,xres,xrot,yUL,yrot,yres

def parseNFO(file):
    try:
        f = open(file)
        j = 0
        lines = 28
        key = []
        value = []
        for line in f:
            (k,v) = line.split('=')
            key.append(k.rstrip())
            value.append(v.rstrip())
            j = j+1
        return key,value
    except ValueError as e:
        raise ValueError("NFO file is malformed")
    
def writeGDALFile(filename,transform,proj,datar,datag,datab):
    (y,x) = datar.shape
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    dst_datatype = gdal.GDT_Byte
    dst_ds = driver.Create(filename,x,y,3,dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(datar)
    dst_ds.GetRasterBand(2).WriteArray(datag)
    dst_ds.GetRasterBand(3).WriteArray(datab)
    dst_ds.SetGeoTransform(transform)
    dst_ds.SetProjection(proj)
    return 1

def read_gdal_file(file,bearing,plot=0):
    filehandle = gdal.Open(file)
    banddata = filehandle.GetRasterBand(1)
    datar = banddata.ReadAsArray()
    banddatag = filehandle.GetRasterBand(2)
    datag = banddatag.ReadAsArray()
    banddatab = filehandle.GetRasterBand(3)
    datab = banddatab.ReadAsArray()
    if bearing > 90 and bearing < 270:
        datar = np.fliplr(np.flipud(datar))
        datag= np.fliplr(np.flipud(datag))
        datab= np.fliplr(np.flipud(datab))
    if plot == 1:
        pyplot.figure()
        pyplot.imshow(datag)
        pyplot.show()
    return filehandle.RasterXSize,filehandle.RasterYSize,datar,datag,datab

def calcResolution(alt):
    # Define parameters for Scout camera
    xsize = 2592
    ysize = 1944
    xsizemm = 5.7024
    ysizemm = 4.2768
    fl = 7.5
    xfov = 2*math.atan(xsizemm/(2*fl))
    yfov = 2*math.atan(ysizemm/(2*fl))
    xfovm = 2*math.tan(xfov/2)*alt
    yfovm = 2*math.tan(yfov/2)*alt
    xrescm = xfovm/xsize
    yrescm = -1*yfovm/ysize
    return xrescm,yrescm




if __name__=='__main__':
    # parse command line options
    
    cl = sys.argv[1:]
    
    i = 1 # i will be start index of file list
    count = 0 # count will be number of files in file list
    plot = 0 
    for item in cl:
        item.rstrip()
        if item == '-d' and mpl == True:
            #debug plot mode
            plot = 1
            i = i+1
        elif item == '-align':
            align = True
            i = i+1
        else:
            count = count+1
    
    if count ==0:
        Usage()
   
    # Define UTM 11N projection with WKT.  This would need to be generalized for 
    # data in other areas.  Test data is in 11N.
    proj = "PROJCS[\"WGS 84 / UTM zone 11N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.01745329251994328,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",-117],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],AUTHORITY[\"EPSG\",\"32611\"],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]"
    
    files = sys.argv[i:]
    salen = count 
    i = 0
    lat = np.zeros([salen])
    lon = np.zeros([salen])
    alt = np.zeros([salen])
    bearing = np.zeros([salen])
    tranform = ('','','','','','')
    
    for file in files:
        print file
        if 'jpg' in file: 
                a = getExif(file)
                (lat[i],lon[i],alt[i],bearing[i]) = getGPSData(a)
                calcImageParams(file,lat[i],lon[i],alt[i],bearing[i])
                print lat[i],' ',lon[i],' ',alt[i],' ',bearing[i]
                i = i + 1
    
        elif 'dng' in file:
            (xsize,ysize,data) = read_gdal_file(file,bearing,2)
            nfoFile = file.replace('dng','nfo')
            key,value = parseNFO(nfoFile)
            lat = float(value[key.index('gps_lat_deg')])
            lon = float(value[key.index('gps_lon_deg')])
            alt = float(value[key.index('alt_agl')])
            bearing = float(value[key.index('yaw_deg')])
            trans = createTransform(lat,lon,alt,bearing)
            tifFile = file.replace('dng','tif')
            writeGDALFile('temp.tif',trans,proj,data)
            cmd = 'gdalwarp -t_srs EPSG:32611 temp.tif '+tifFile
            print cmd
            os.system(cmd)
    
        elif 'tif' in file:
            try:
                nfoFile = file.replace('tif','nfo')
                tiffFile = file.replace('tif','tiff')

                # If the output file exists, skip.  TODO, make it possible to
                # override with --force or --overwrite.
                if os.path.exists(tiffFile):
                    print "File already processed, skipping.\n"
                    continue

                key,value = parseNFO(nfoFile)
                lat = float(value[key.index('gps_lat_deg')])
                lon = float(value[key.index('gps_lon_deg')])
                lat = lat*np.pi/180.0
                lon = lon*np.pi/180.0
                alt = float(value[key.index('alt_agl')])
                bearing = float(value[key.index('yaw_deg')])
                (xsize,ysize,datar,datag,datab) = read_gdal_file(file,bearing,plot)
                # bearing is corrected as data are read in as though they are collected along northern bearing line
                if bearing > 90 and bearing < 270:
                    bearing = bearing - 180
                trans = createTransform(lat,lon,alt,bearing)
                writeGDALFile('temp.tif',trans,proj,datar,datag,datab)
                cmd = 'gdalwarp -dstnodata 0,0,0  -t_srs EPSG:32611 temp.tif '+tiffFile+' > /dev/null'
                #print cmd
                os.system(cmd)
            except IOError as e:
                print "Unable to process file (" + file + "): I/O error({0}): {1}".format(e.errno, e.strerror)
            except ValueError as e:
                print "Unable to process file (" + file + "): " + str(e)

