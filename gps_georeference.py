#!/usr/bin/env python

# Copyright (c) 2013, Scott Arko, Bruce Crevensten
import pprint
import argparse
import logging
from georeferencer import *
from preprocessors import *
from postprocessors import *


# For debugging.  Example: pp.pprint(whatever)
pp = pprint.PrettyPrinter(indent=4)


def get_input_processor_from_file(file):
    """
    Figure out what camera is being used, and return an appropriate
    decoder / preprocessor object.  Possible cases that aren't implemented
    yet return a NullPreprocessor.

    Incoming file is assumed to be a file handle.
    """
    preprocessor = None

    # If the source file is a .jpg, check the EXIF tags
    if 'jpg' in file.name:
        exif = parse_exif()
        if('BOARD_FLIR_TAU_640' == exif['Model']):
            preprocessor = ScoutFlirJpgPreprocessor(file)

    elif 'dng' in file.name:
    # If the source file is a .dng read the associated .nfo to figure the camera
        nfoFilename = None  # magic business to get the right name needed here
        nfo = parse_nfo(nfoFilename)
        if('BOARD_MT9P031_SUNEX' == nfo['camera_model']):
            preprocessor = NullPreprocessor(file)  # This is a high priority path.

    elif 'nfo' in file.name:
    # If the source file is .nfo, read that
        nfo = parse_nfo(file)
        if('BOARD_MT9P031_SUNEX' == nfo['camera_model']):
            preprocessor = ScoutPhotoS3NfoPreprocessor(file)

    else:
        preprocessor = NullPreprocessor(file)

    return preprocessor


if __name__ == '__main__':

    # Set up default values, these should be made configurable via command line in the future

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

    postprocessors = []

    ## TBD: parallelize this after watching to see if that'd impact some of the
    ## gdal / ufraw-batch work.
    for file in args.filelist:
        logging.info('Processing {}...'.format(file.name))

        try:
            pre = get_input_processor_from_file(file)
            pre.process()
            r = GeoReferencer(pre)
            trans = r.transform()

            for postprocessor in postprocessors:
                postprocessor.process(file, trans)

        except IOError as e:
            logging.warning("Unable to process file (" + file.name + "): I/O error({0}): {1}".format(e.errno, e.strerror))
        except ValueError as e:
            logging.warning("Unable to process file (" + file.name + "): " + str(e))

    logging.info('Finished.')
