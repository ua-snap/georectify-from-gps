from pyproj import transform, Proj
import math


# From a GeoImage, perform the georectification calculations
# Starting by hardcoding for this to be just the FLIR, will make it a superclass later
class GeoReferencer:
    def __init__(self, geoImage):
        self.geoImage = geoImage

    def transform(self):

        inProj = Proj(init='EPSG:4326')
        outProj = Proj(init='EPSG:3857')

        logging.debug(self.geoImage.lon)
        logging.debug(self.geoImage.lat)
        logging.debug(pprint.pformat(self.geoImage))

        (east, north) = transform(inProj, outProj, self.geoImage.lon, self.geoImage.lat)

        logging.debug('Original lon {} lat {} Transformed lon {} lat {}'.format(self.geoImage.lon, self.geoImage.lat, east, north))
        logging.debug('self.geoImage.xFovRad={}'.format(self.geoImage.xFovRad))
        logging.debug('self.geoImage.yFovRad={}'.format(self.geoImage.xFovRad))

        self.xFovM = self.geoImage.alt * math.tan(self.geoImage.xFovRad / 2)
        self.yFovM = self.geoImage.alt * math.tan(self.geoImage.yFovRad / 2)

        logging.debug("xFovM {} yFovM {}".format(self.xFovM, self.yFovM))
        logging.debug("Footprint in meters: {} x {} y".format(self.xFovM * 2, self.yFovM * 2))

        offsetUrX, offsetUrY = rotate(self.geoImage.bearing, self.xFovM, self.yFovM)  # ur = upper right
        offsetUlX, offsetUlY = rotate(self.geoImage.bearing, -self.xFovM, self.yFovM)  # ul
        offsetLlX, offsetLlY = rotate(self.geoImage.bearing, -self.xFovM, -self.yFovM)  # ll
        offsetLrX, offsetLrY = rotate(self.geoImage.bearing, self.xFovM, -self.yFovM)  # lr

        logging.debug("Offsets: {} {}\n{} {}\n{} {}\n{} {}".format(offsetUrX, offsetUrY, offsetUlX, offsetUlY, offsetLlX, offsetLlY, offsetLrX, offsetLrY))

        self.lonUR = east + offsetUrX
        self.latUR = north + offsetUrY
        self.lonUL = east + offsetUlX
        self.latUL = north + offsetUlY
        self.lonLL = east + offsetLlX
        self.latLL = north + offsetLlY
        self.lonLR = east + offsetLrX
        self.latLR = north + offsetLrY

        return {
            'lonUR': self.lonUR,
            'latUR': self.latUR,
            'lonUL': self.lonUL,
            'latUL': self.latUL,
            'lonLL': self.lonLL,
            'latLL': self.latLL,
            'lonLR': self.lonLR,
            'latLR': self.latLR
        }


# Convenience function to rotate a point in space.
def rotate(angle, x, y):
    angle = math.radians(angle)
    xr = (x * math.cos(angle)) - (y * math.sin(angle))
    yr = (x * math.sin(angle)) + (y * math.cos(angle))
    return xr, yr


"""
swamp for reintegrating




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

"""