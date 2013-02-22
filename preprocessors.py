import math
from PIL import Image
from PIL.ExifTags import TAGS
import pprint


def parse_nfo(file):
        nfo = {}
        try:
            for line in file:
                (k, v) = line.split('=')
                k = k.rstrip()
                v = v.rstrip()
                nfo[k] = v
            return nfo
        except ValueError as e:
            raise ValueError("NFO file is malformed (" + e + ")")


def parse_exif(file):
        exif = {}
        i = Image.open(file)
        info = i._getexif()
        pprint.pprint(info)
        for tag, value in info.items():
            decoded = TAGS.get(tag, tag)
            exif[decoded] = value
        return exif


# Responsible for figuring out image and camera properties from a file
class ImagePreprocessor:
    def __init__(self, file):
        self.file = file

        # Image properties
        self.lat = None
        self.lon = None
        self.alt = None
        self.bearing = None

        # Camera properties
        self.xRes = None  # xResolution of physical sensor (pixels)
        self.yRes = None  # yResolution of physical sensor (pixels)
        self.xFovRad = None  # xFOV in Radians
        self.yFovRad = None  # yFov in Radians
        self.xSizeMm = None  # x sensor size (mm)
        self.ySizeMm = None  # y sensor size (mm)
        self.focalLength = None  # focal length (mm)

    # To be implemented by child classes
    def parse_camera_params(self):
        pass

    # to be implemented by child classes
    def parse_image_params(self):
        pass

    def process(self):
        self.parse_camera_params()
        self.parse_image_params()


class NullPreprocessor(ImagePreprocessor):
    def process(self):
        raise 'Could not find preprocessor for file {}'.format(self.file.name)


class ScoutFlirJpgPreprocessor(ImagePreprocessor):
    def parse_camera_params(self):
        self.xFovRad = math.radians(32)
        self.yFovRad = math.radians(26)

    def parse_image_params(self):
        exif = parse_exif()
        # Read the exif data returned by getExif and
        # parse out GPS lat, lon, alt and bearing
        lat = [float(x) / float(y) for x, y in exif['GPSInfo'][2]]
        latref = exif['GPSInfo'][1]
        lon = [float(x) / float(y) for x, y in exif['GPSInfo'][4]]
        lonref = exif['GPSInfo'][3]
        alt = float(exif['GPSInfo'][6][0]) / float(exif['GPSInfo'][6][1])
        bearing = float(exif['GPSInfo'][24][0]) / float(exif['GPSInfo'][24][1])

        # Convert bearing/lat/lon to usable values
        if bearing > 90 and bearing < 270:
            bearing = bearing - 180

        lat = lat[0] + lat[1] / 60 + lat[2] / 3600
        lon = lon[0] + lon[1] / 60 + lon[2] / 3600
        if latref == 'S':
            lat = -lat
        if lonref == 'W':
            lon = -lon

        logging.debug('parseFromExif: lat {} lon {} alt {} bearing {}'.format(lat, lon, alt, bearing))
        (self.lat, self.lon, self.alt, self.bearing) = (lat, lon, alt, bearing)


class ScoutPhotoS3NfoPreprocessor(ImagePreprocessor):

    def parse_camera_params(self):
        self.xRes = 2592
        self.yRes = 1944
        self.xSizeMm = 5.7024
        self.ySizeMm = 4.2768
        self.focalLength = 7.5
        self.xFovRad = 0.557711
        self.yFovRad = 0.423049

    def parse_image_params(self):
        # TODO: make this check if it should use the NFO or other means.
        nfo = parse_nfo()
        self.lat = float(nfo['gps_lat_deg'])
        self.lon = float(nfo['gps_lon_deg'])
        self.bearing = float(nfo['yaw_deg'])
        self.alt = float(nfo['alt_agl'])
        logging.debug('parsed from NFO: lat {} lon {} alt {} bearing {}'.format(self.lat, self.lon, self.alt, self.bearing))


"""
Swamp for reintegrating



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


"""
