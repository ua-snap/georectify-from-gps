#!/usr/bin/env python

###############################################################################
# importScout.py 
#
# Project:  Aeryon Scout data processing
# Purpose:  Reads images from the Aeryon scout system and geocodes them to UTM projection 
# Authors:   Scott Arko, Bruce Crevensten
#
###############################################################################
# Copyright (c) 2013, Scott Arko, Bruce Crevensten 
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
from pyproj import transform, Proj, Geod
import pprint

# import local modules
sys.path.append('.')
import transformations

# Optional modules
try:
    import vlfeat as vl
    sift = True
except ImportError:
    sift = False

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
    print "    importScout.py [-d] [-align] [-sift] filelist"
    print ""
    print " -d -> Debug mode that will display an image of the data being read for verification.  Requires matplotlib python binding installed"
    print " -align -> Performs FFT-based alignment of adjacent images.  NOT IMPLEMENTED AT CURRENT TIME"
    print " -sift ->  Performs SIFT-based alignment of adjacent images.  Better for images with rotation error.  NOT IMPLEMENTED AT CURRENT TIME"
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
    lat = [float(x) / float(y) for x, y in exif['GPSInfo'][2]]
    latref = exif['GPSInfo'][1]
    lon = [float(x) / float(y) for x, y in exif['GPSInfo'][4]]
    lonref = exif['GPSInfo'][3]
    alt = float(exif['GPSInfo'][6][0]) / float(exif['GPSInfo'][6][1])
    bearing = float(exif['GPSInfo'][24][0]) / float(exif['GPSInfo'][24][1])

    lat = lat[0] + lat[1] / 60 + lat[2] / 3600
    lon = lon[0] + lon[1] / 60 + lon[2] / 3600
    if latref == 'S':
        lat = -lat
    if lonref == 'W':
        lon = -lon

    return lat, lon, alt, bearing


def calcImageParams(file, lat, lon, alt, bearing):
    # This function calculates the necessary image parameters to generate a
    # standard world file.  Requires lat,lon,alt,bearing provided by getGPSData
    fh = gdal.Open(file)
    #xsize = fh.RasterXSize
    #ysize = fh.RasterYSize
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init='epsg:32611')
    (east, north) = transform(inProj, outProj, lon, lat)
    xres = 0.0064
    yres = 0.0064
    #xUL = east - (xres * xsize / 2)
    #yUL = north + (yres * ysize / 2)
    bearingRad = bearing * 3.14159 / 180
    tanb = math.tan(bearingRad)
    line2 = -1 * tanb * xres
    line3 = -1 * tanb * yres


def createTransform(lat, lon, alt, bearing):
    # This function calculates the necessary image parameters to generate a
    # standard world file.  Requires lat,lon,alt,bearing provided by getGPSData
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init='epsg:4326')
    lat = lat * 180.0 / np.pi
    lon = lon * 180.0 / np.pi
    east, north = transform(inProj, outProj, lon, lat)
    print 'Center Point: ' + str(east) + ' ' + str(north)
    (xres, yres) = calcResolution(alt)

    # hard coded to scout parameters for PhotoS3 camera
    xsize = 640  # for flir tau, 2592=s3
    ysize = 480  # for flir tau, 1944=s3

    # Calculate the distance (in m from center to corner of image)
    xd = xres * xsize / 2
    yd = yres * ysize / 2
    radius = math.sqrt(xd ** 2 + yd ** 2)
    offset = math.atan(xd / yd)
    bearingRad = bearing * np.pi / 180
    xUL = east + radius * math.sin(bearingRad + offset)
    yUL = north + radius * math.cos(bearingRad + offset)

    if bearing > 180:
        b2 = bearing - 360
        bearingRad = b2 * np.pi / 180.0

    xrot = -1 * xres * math.sin(bearingRad)
    yrot = yres * math.sin(bearingRad)
    xres = xres * math.cos(bearingRad)
    yres = yres * math.cos(bearingRad)

    return xUL, xres, xrot, yUL, yrot, yres


def rotate(angle, x, y):
    angle = math.radians(angle)
    xr = (x * math.cos(angle)) - (y * math.sin(angle))
    yr = (x * math.sin(angle)) + (y * math.cos(angle))
    return xr, yr


def createBoxTransform(lat, lon, alt, bearing):
    # This function calculates the necessary image parameters to generate a
    # standard world file.  Requires lat,lon,alt,bearing provided by getGPSData
    inProj = Proj(init='EPSG:4326')
    outProj = Proj(init='EPSG:3338')
    lat = lat * 180.0 / np.pi
    lon = lon * 180.0 / np.pi
    east, north = transform(inProj, outProj, lon, lat)

    print 'Center Point: ' + str(east) + ' ' + str(north)

    (xSizeM, ySizeM) = calcResolutionFlir(alt)

    print "Altitude {} Bearing {}".format(alt, bearing)
    print "Footprint in meters: {} x {} y".format(xSizeM * 2, ySizeM * 2)

    print "xSizeM {} ySizeM {}".format(xSizeM, ySizeM)

    offsetUrX, offsetUrY = rotate(bearing, xSizeM, ySizeM)  # ur = upper right
    offsetUlX, offsetUlY = rotate(bearing, -xSizeM, ySizeM)  # ul
    offsetLlX, offsetLlY = rotate(bearing, -xSizeM, -ySizeM)  # ll
    offsetLrX, offsetLrY = rotate(bearing, xSizeM, -ySizeM)  # lr

    print "Offsets: {} {}\n{} {}\n{} {}\n{} {}".format(offsetUrX, offsetUrY, offsetUlX, offsetUlY, offsetLlX, offsetLlY, offsetLrX, offsetLrY)

    lonUR = east + offsetUrX
    latUR = north + offsetUrY
    lonUL = east + offsetUlX
    latUL = north + offsetUlY
    lonLL = east + offsetLlX
    latLL = north + offsetLlY
    lonLR = east + offsetLrX
    latLR = north + offsetLrY

    print "{},{},0\n{},{},0\n{},{},0\n{},{},0".format(lonUR, latUR, lonUL, latUL, lonLL, latLL, lonLR, latLR)

    return {
        'lonUR': lonUR,
        'latUR': latUR,
        'lonUL': lonUL,
        'latUL': latUL,
        'lonLL': lonLL,
        'latLL': latLL,
        'lonLR': lonLR,
        'latLR': latLR
    }


def parseNFO(file):
    try:
        f = open(file)
        j = 0
        #lines = 28
        key = []
        value = []
        for line in f:
            (k, v) = line.split('=')
            key.append(k.rstrip())
            value.append(v.rstrip())
            j = j + 1
        return key, value
    except ValueError as e:
        raise ValueError("NFO file is malformed (" + e + ")")


def writeGDALFile(filename, transform, proj, datar, datag, datab):
    (y, x) = datar.shape
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    dst_datatype = gdal.GDT_UInt16
    dst_ds = driver.Create(filename, x, y, 3, dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(datar)
    dst_ds.GetRasterBand(2).WriteArray(datag)
    dst_ds.GetRasterBand(3).WriteArray(datab)
    dst_ds.SetGeoTransform(transform)
    dst_ds.SetProjection(proj)
    return 1


def writeSimpleFile(filename, data):
    (y, x) = data.shape
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    dst_datatype = gdal.GDT_Byte
    dst_ds = driver.Create(filename, x, y, 1, dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(data)
    return 1


def read_gdal_file(file, bearing, plot=0):
    filehandle = gdal.Open(file)
    banddata = filehandle.GetRasterBand(1)
    datar = banddata.ReadAsArray() + 1
    banddatag = filehandle.GetRasterBand(2)
    datag = banddatag.ReadAsArray() + 1
    banddatab = filehandle.GetRasterBand(3)
    datab = banddatab.ReadAsArray() + 1
    if bearing > 90 and bearing < 270:
        datar = np.fliplr(np.flipud(datar))
        datag = np.fliplr(np.flipud(datag))
        datab = np.fliplr(np.flipud(datab))
    if plot == 1:
        pyplot.figure()
        pyplot.imshow(datag)
        pyplot.show()
    return filehandle.RasterXSize, filehandle.RasterYSize, datar, datag, datab


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
    print xrescm
    print yrescm
    return xrescm,yrescm


def calcResolutionFlir(altM):

    xfovRad = 32 * (np.pi / 180)
    yfovRad = 26 * (np.pi / 180)

    xfovm = altM * math.tan(xfovRad / 2)
    yfovm = altM * math.tan(yfovRad / 2)

    print "xFovM {} yFovM {}".format(xfovm, yfovm)
    return xfovm, yfovm


def alignImages(master,slave,plot=0,aligntype=0):
    # Open the files
    mh = gdal.Open(master)
    sh = gdal.Open(slave)
    # Get band 1 data and read as a numpy array
    md = mh.GetRasterBand(2)
    sd = sh.GetRasterBand(2)
    # Get the geo transoform
    mtrans = mh.GetGeoTransform()
    strans = sh.GetGeoTransform()
    # Need to get corner coords from transform
    mxsize = mh.RasterXSize
    mysize = mh.RasterYSize
    sxsize = sh.RasterXSize
    sysize = sh.RasterYSize
    muleast = mtrans[0]
    mulnorth = mtrans[3]
    suleast = strans[0]
    sulnorth = strans[3]
    mlleast = muleast + mxsize*mtrans[1]
    mllnorth = mulnorth + mysize*mtrans[5]
    slleast = suleast + sxsize*strans[1]
    sllnorth = sulnorth + sysize*strans[5]

    # Determine the extent of overlap and the indices 
    # for master and slave overlap region
    if muleast > suleast:
        uleast = muleast
        muleindex = 0
        suleindex = int((muleast-suleast)/strans[1])
    else:
        uleast = suleast
        suleindex = 0
        muleindex = int((suleast-muleast)/mtrans[1])
    if mlleast < slleast:
        lleast = mlleast
        mlleindex = mxsize-1
        slleindex = sxsize - int((slleast-mlleast)/strans[1])
    else:
        lleast = slleast
        slleindex = sxsize-1
        mlleindex = mxsize - int((mlleast-slleast)/mtrans[1])

    if mulnorth < sulnorth:
        ulnorth = mulnorth
        mulnindex = 0
        sulnindex = int(-1*(sulnorth-mulnorth)/strans[5])
    else:
        ulnorth = sulnorth
        sulnindex = 0
        mulnindex = int(-1*(mulnorth-sulnorth)/mtrans[5])
    if mllnorth > sllnorth:
        llnorth = mllnorth
        mllnindex = mysize-1 
        sllnindex = sysize - int((sllnorth - mllnorth)/strans[5])
    else:
        llnorth = sllnorth
        sllnindex = sysize-1 
        mllnindex = mysize - int((mllnorth - sllnorth)/mtrans[5])

    # Extract overlap region in master and use it to ingest overlap region 
    # from slave at same resolution
    xsize = (mllnindex-mulnindex)/5
    ysize = (mlleindex-muleindex)/5
    print xsize,' ',ysize
    mdata = md.ReadAsArray(muleindex,mulnindex,mlleindex-muleindex,mllnindex-mulnindex,ysize,xsize)
    #mdata = mdtemp[mulnindex:mllnindex,muleindex:mlleindex]
    #xsize,ysize = np.shape(mdata)
    sdata = sd.ReadAsArray(suleindex,sulnindex,slleindex-suleindex,sllnindex-sulnindex,ysize,xsize)

    # Perform a 2-sigma stretch of the non-zero elements 
    mdmean = np.mean(mdata.ravel()[np.flatnonzero(mdata)])
    mdstd = np.std(mdata.ravel()[np.flatnonzero(mdata)])
    sdmean = np.mean(sdata.ravel()[np.flatnonzero(sdata)])
    sdstd = np.std(sdata.ravel()[np.flatnonzero(sdata)])

    md2 = (mdata - mdmean)/(2*mdstd)
    np.putmask(md2,md2>1,1)
    np.putmask(md2,md2<-1,-1)
    md2 = md2*127 + 127

    mdata = md2

    sd2 = (sdata - sdmean)/(2*sdstd)
    np.putmask(sd2,sd2>1,1)
    np.putmask(sd2,sd2<-1,-1)
    sd2 = sd2*127 + 127

    sdata = sd2

    # calculate scaling factors between two images (this should only be the result of 
    # pixel size difference if all is correct
    yscale = float((slleindex-suleindex))/float(ysize)
    xscale = float((sllnindex-sulnindex))/float(xsize)

    if plot == 1:
        pyplot.figure()
        pyplot.imshow(mdata,cmap='gray',interpolation='nearest')
        pyplot.figure()
        pyplot.imshow(sdata,cmap='gray',interpolation='nearest')
        pyplot.show()

    ref_patch = mdata
    test_patch = sdata

    if aligntype == 0:

        # This is SIFT align.  It is different than FFT-based and requires that the images being aligned
        # are oriented to their centers.  
        xsize = (int(mh.RasterXSize)-400)/4
        ysize = (int(mh.RasterYSize)-400)/4

        ref_patch = md.ReadAsArray(200,200,int(mh.RasterXSize)-200,int(mh.RasterYSize)-200,xsize,ysize)
        test_patch = sd.ReadAsArray(200,200,int(sh.RasterXSize)-200,int(sh.RasterYSize)-200,xsize,ysize)

        writeSimpleFile('file1.tif',ref_patch) 
        writeSimpleFile('file2.tif',test_patch) 
        vl.process_image('file1.tif','file1.sift')
        l1,d1 = vl.read_features_from_file('file1.sift')
        vl.process_image('file2.tif','file2.sift')
        l2,d2 = vl.read_features_from_file('file2.sift')
        m = vl.match_twosided(d1,d2)

        im = np.array(Image.open('file1.tif'))
        im2 = np.array(Image.open('file2.tif'))
        pyplot.figure()
        vl.plot_features(im,l1,True)
        pyplot.gray()

        pyplot.figure()
        pyplot.imshow(im)
        
        # draw lines for matches
        for i in range(len(m)):
            if m[i] > 0:
                pyplot.plot([l1[i,0], l2[m[i,0],0]], [l1[i,1], l2[m[i,0],1]], 'c')



        pyplot.figure()
        vl.plot_matches(im,im2,l1,l2,m)
        pyplot.gray()

        pyplot.figure()
        for i in range(len(m)):
         if m[i] > 0:
             pyplot.scatter(l1[i,0], l2[m[i,0],0]) #, [locs1[i,1], locs2[matchscores[i,0],1]], 'c')

        pyplot.figure()
        for i in range(len(m)):
         if m[i] > 0:
             pyplot.scatter(l1[i,1], l2[m[i,0],1])

        pyplot.figure()
        for i in range(len(m)):
         if m[i] > 0:
             pyplot.scatter(l1[i,1], math.sqrt((l2[m[i,0],1]-l1[i,1])**2 + (l2[m[i,0],0]-l1[i,0])**2))

        pyplot.show()

    else:
    
        fft1 = np.fft.fft2(ref_patch)
        fft2 = np.fft.fft2(test_patch)
        fft2_conj = np.conjugate(fft2)
        
        result = fft1 * fft2_conj/abs(fft1 * fft2_conj)
        shift = np.fft.ifft2(result)

        realshift = abs(shift)
        shiftmax = realshift.max()
        shiftmin = realshift.min()
        shiftmean = realshift.mean()
        shiftsnr = shiftmax / shiftmean

        maxloc = np.argmax(realshift)
        dims = realshift.shape
        maxindex = np.unravel_index(maxloc, dims)
        xloc = maxindex[0]
        yloc = maxindex[1]
        snr = shiftsnr
        
        if xloc > xsize / 2:
            xloc = xloc - xsize
        if yloc > ysize / 2:
            yloc = yloc - ysize
        
        radloc = np.sqrt(xloc * xloc + yloc * yloc)
    
        xloc = float(xloc) * mtrans[1] * 5.0
        yloc = float(yloc) * mtrans[5] * 5.0
    
        print "%.2f, %.2f, %.2f, %.2f" % (xloc, yloc, radloc, snr)
    
    return 1
    


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
        elif item == '-sift' and sift==True:
            siftalign = True
            i = i+1
        else:
            count = count+1
    
    if count ==0:
        Usage()
   
    # Define UTM 11N projection with WKT.  This would need to be generalized for 
    # data in other areas.  Test data is in 16N.
    proj = "PROJCS[\"WGS 84 / UTM zone 16N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.01745329251994328,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",-117],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],AUTHORITY[\"EPSG\",\"32611\"],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]"

    files = sys.argv[i:]
    salen = count
    i = 0
    lat = np.zeros([salen])
    lon = np.zeros([salen])
    alt = np.zeros([salen])
    bearing = np.zeros([salen])

    for file in files:
        print file
        if 'jpg' in file:
            try:
                exif = getExif(file)
                (lat, lon, alt, bearing) = getGPSData(exif)
                calcImageParams(file, lat, lon, alt, bearing)
                lat = lat * np.pi / 180.0
                lon = lon * np.pi / 180.0
                # bearing is corrected as data are read in as though they are collected along northern bearing line
                if bearing > 90 and bearing < 270:
                    bearing = bearing - 180
                trans = createBoxTransform(lat, lon, alt, bearing)
                trans['sourceFile'] = file
                trans['outputFile'] = os.path.splitext(file)[0] + '.tif'
                gdalTransformCommand = 'gdal_translate -of GTiff -a_srs EPSG:3338 -gcp 0 0 {lonUL} {latUL} -gcp 640 0 {lonUR} {latUR} -gcp 0 480 {lonLL} {latLL} -gcp 640 480 {lonLR} {latLR} {sourceFile} {outputFile}'.format(**trans)
                print gdalTransformCommand
                #os.system(gdalTransformCommand)
                gdalWarpCommand = 'gdalwarp -dstnodata 0,0,0,0 -t_srs EPSG:3786 {outputFile} {outputFile}'.format(**trans)
                print gdalWarpCommand
                #os.system(gdalWarpCommand)

            except IOError as e:
                print "Unable to process file (" + file + "): I/O error({0}): {1}".format(e.errno, e.strerror)
            except ValueError as e:
                print "Unable to process file (" + file + "): " + str(e)

        elif 'dng' in file:
            (xsize, ysize, data) = read_gdal_file(file, bearing, 2)
            nfoFile = file.replace('dng', 'nfo')
            key, value = parseNFO(nfoFile)
            lat = float(value[key.index('gps_lat_deg')])
            lon = float(value[key.index('gps_lon_deg')])
            alt = float(value[key.index('alt_agl')])
            bearing = float(value[key.index('yaw_deg')])
            trans = createTransform(lat, lon, alt, bearing)
            tifFile = file.replace('dng', 'tif')
            writeGDALFile('temp.tif', trans, proj, data)
            cmd = 'gdalwarp -t_srs EPSG:32611 temp.tif ' + tifFile
            print cmd
            os.system(cmd)

        elif 'tif' in file:
            try:
                nfoFile = file.replace('tif', 'nfo')
                tiffFile = file.replace('tif', 'tiff')

                # If the output file exists, skip.  TODO, make it possible to
                # override with --force or --overwrite.
                if False and os.path.exists(tiffFile):
                    print "File already processed, skipping.\n"
                    continue

                key, value = parseNFO(nfoFile)
                lat = float(value[key.index('gps_lat_deg')])
                lon = float(value[key.index('gps_lon_deg')])
                lat = lat * np.pi / 180.0
                lon = lon * np.pi / 180.0
                alt = float(value[key.index('alt_agl')])
                bearing = float(value[key.index('yaw_deg')])
                (xsize, ysize, datar, datag, datab) = read_gdal_file(file, bearing, plot)
                # bearing is corrected as data are read in as though they are collected along northern bearing line
                if bearing > 90 and bearing < 270:
                    bearing = bearing - 180
                trans = createTransform(lat, lon, alt, bearing)
                writeGDALFile('temp.tif', trans, proj, datar, datag, datab)
                cmd = 'gdalwarp -dstnodata 0,0,0  -t_srs EPSG:32611 temp.tif ' + tiffFile
                os.system(cmd)
            except IOError as e:
                print "Unable to process file (" + file + "): I/O error({0}): {1}".format(e.errno, e.strerror)
            except ValueError as e:
                print "Unable to process file (" + file + "): " + str(e)

