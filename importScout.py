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
from pyproj import transform, Proj
import pprint
import argparse
import logging

# For debugging.  Example: pp.pprint(whatever)
pp = pprint.PrettyPrinter(indent=4)

# Define the template for the KML overlay.  Trailing slash on the first line
# keeps the XML declaration valid.
kmlTemplate = """\
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
    <name>KML Preview</name>
    <open>1</open>
    <description>Coverage preview of assets</description>
 <Style id="examplePolyStyle">
      <LineStyle>
        <width>1.5</width>
      </LineStyle>
      <PolyStyle>
        <color>7d00ffff</color>
      </PolyStyle>
  </Style>
  <Placemark>
    <name>KML Coverage Preview</name>
    <styleUrl>#examplePolyStyle</styleUrl>
    <MultiGeometry>
        {}
    </MultiGeometry>
  </Placemark>
  </Document>
</kml>
"""

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

    logging.debug('getGPSData: lat {} lon {} alt {} bearing {}'.format(lat, lon, alt, bearing))
    return lat, lon, alt, bearing


def createTransform(lat, lon, alt, bearing):
    # This function calculates the necessary image parameters to generate a
    # standard world file.  Requires lat,lon,alt,bearing provided by getGPSData
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init='epsg:4326')
    lat = lat * 180.0 / np.pi
    lon = lon * 180.0 / np.pi
    east, north = transform(inProj, outProj, lon, lat)
    logging.debug('Original lon {} lat {} Transformed lon {} lat {}'.format(lon, lat, east, north))
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

    inProj = Proj(init='EPSG:4326')
    outProj = Proj(init='EPSG:3857')
    #bearing += 90  # correct for scout northing

    east, north = transform(inProj, outProj, lon, lat)
    logging.debug('Original lon {} lat {} Transformed lon {} lat {}'.format(lon, lat, east, north))

    logging.debug('Center Point: ' + str(east) + ' ' + str(north))

    (xSizeM, ySizeM) = calcResolutionFlir(alt)

    logging.debug("Altitude {} Bearing {}".format(alt, bearing))
    logging.debug("Footprint in meters: {} x {} y".format(xSizeM * 2, ySizeM * 2))
    logging.debug("xSizeM {} ySizeM {}".format(xSizeM, ySizeM))

    offsetUrX, offsetUrY = rotate(bearing, xSizeM, ySizeM)  # ur = upper right
    offsetUlX, offsetUlY = rotate(bearing, -xSizeM, ySizeM)  # ul
    offsetLlX, offsetLlY = rotate(bearing, -xSizeM, -ySizeM)  # ll
    offsetLrX, offsetLrY = rotate(bearing, xSizeM, -ySizeM)  # lr

    logging.debug("Offsets: {} {}\n{} {}\n{} {}\n{} {}".format(offsetUrX, offsetUrY, offsetUlX, offsetUlY, offsetLlX, offsetLlY, offsetLrX, offsetLrY))

    lonUR = east + offsetUrX
    latUR = north + offsetUrY
    lonUL = east + offsetUlX
    latUL = north + offsetUlY
    lonLL = east + offsetLlX
    latLL = north + offsetLlY
    lonLR = east + offsetLrX
    latLR = north + offsetLrY

    logging.debug("{},{},0\n{},{},0\n{},{},0\n{},{},0".format(lonUR, latUR, lonUL, latUL, lonLL, latLL, lonLR, latLR))

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


def getKmlPoly(trans):

    outProj = Proj(init='EPSG:4326')
    inProj = Proj(init='EPSG:3857')

    kml = dict()
    kml['lonUL'], kml['latUL'] = transform(inProj, outProj, trans['lonUL'], trans['latUL'])
    kml['lonUR'], kml['latUR'] = transform(inProj, outProj, trans['lonUR'], trans['latUR'])
    kml['lonLL'], kml['latLL'] = transform(inProj, outProj, trans['lonLL'], trans['latLL'])
    kml['lonLR'], kml['latLR'] = transform(inProj, outProj, trans['lonLR'], trans['latLR'])

    return """
    <Polygon>
      <extrude>1</extrude>
      <altitudeMode>clampToGround</altitudeMode>
      <outerBoundaryIs>
        <LinearRing>
          <coordinates>
            {lonUL},{latUL},0
            {lonUR},{latUR},0
            {lonLR},{latLR},0
            {lonLL},{latLL},0
            {lonUL},{latUL},0
          </coordinates>
        </LinearRing>
      </outerBoundaryIs>
    </Polygon>
    """.format(**kml)


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

    logging.debug("xFovM {} yFovM {}".format(xfovm, yfovm))
    return xfovm, yfovm


if __name__ == '__main__':

    # Set up default values, these should be made configurable via command line in the future

    # Tempfile location for scratch work when generating geotifs
    tempFile = '/tmp/temp.tif'

    # Define UTM 11N projection with WKT.  This would need to be generalized for
    # data in other areas.  Test data is in 16N.
    proj = 'PROJCS["NAD83 / Alaska Albers",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",55],PARAMETER["standard_parallel_2",65],PARAMETER["latitude_of_center",50],PARAMETER["longitude_of_center",-154],PARAMETER["false_easting",0],PARAMETER["false_northing",0],AUTHORITY["EPSG","3338"],AXIS["X",EAST],AXIS["Y",NORTH]]'

    # Handle command line options
    parser = argparse.ArgumentParser(description='Georeference aerial images by GPS coordinates, and create coverage previews.')
    parser.add_argument('--kml', '-k', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('filelist', nargs=argparse.REMAINDER, type=file)
    args = parser.parse_args()
    pp.pprint(args)

    if(args.debug == True):
        logLevel = logging.DEBUG
    else:
        logLevel = logging.INFO

    # Set up logging
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logLevel)

    # Variable for accumulating KML information, if requested
    kml = ''

    for file in args.filelist:
        logging.info('Processing {}...'.format(file.name))

        # Assuming for now that JPG === infrared/FLIR.
        if 'jpg' in file.name:
            try:
                exif = getExif(file)
                (lat, lon, alt, bearing) = getGPSData(exif)

                # bearing is corrected as data are read in as though they are collected along northern bearing line
                if bearing > 90 and bearing < 270:
                    bearing = bearing - 180

                trans = createBoxTransform(lat, lon, alt, bearing)

                trans['sourceFile'] = file.name
                trans['tempFile'] = tempFile
                trans['outputFile'] = os.path.splitext(file.name)[0] + '.tif'

                if(True == args.kml):
                    kml += getKmlPoly(trans)
                else:
                    if(os.path.exists(trans['tempFile'])):
                        os.unlink(trans['tempFile'])

                    if(os.path.exists(trans['outputFile'])):
                        os.unlink(trans['outputFile'])

                    gdalTransformCommand = 'gdal_translate -a_srs EPSG:3857 -of GTiff -gcp 0 0 {lonUL} {latUL} -gcp 640 0 {lonUR} {latUR} -gcp 0 480 {lonLL} {latLL} -gcp 640 480 {lonLR} {latLR} {sourceFile} {tempFile}'.format(**trans)
                    logging.debug(gdalTransformCommand)

                    os.system(gdalTransformCommand)
                    gdalWarpCommand = 'gdalwarp -s_srs EPSG:3857 -t_srs EPSG:4326 -dstalpha -dstnodata 0,0,0,0 {tempFile} {outputFile}'.format(**trans)
                    logging.debug(gdalWarpCommand)
                    os.system(gdalWarpCommand)

            except IOError as e:
                logging.warning("Unable to process file (" + file.name + "): I/O error({0}): {1}".format(e.errno, e.strerror))
            except ValueError as e:
                logging.warning("Unable to process file (" + file.name + "): " + str(e))

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
            logging.debug(cmd)
            os.system(cmd)

        elif 'tif' in file:
            try:
                nfoFile = file.replace('tif', 'nfo')
                tiffFile = file.replace('tif', 'tiff')

                # If the output file exists, skip.  TODO, make it possible to
                # override with --force or --overwrite.
                if False and os.path.exists(tiffFile):
                    logging.info("File already processed, skipping.\n")
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

    if(True == args.kml):
        if(args.kml):
            if(os.path.exists('coverage.kml')):
                os.unlink('coverage.kml')
        kmlFile = open('coverage.kml', 'w')
        kmlFile.write(kmlTemplate.format(kml))
        kmlFile.close()
