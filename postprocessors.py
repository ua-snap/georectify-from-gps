from pyproj import transform, Proj

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

# Variable for accumulating KML information, if requested
kml = ''


def get_kml_poly(trans):
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


"""
Swamp for reintegrating additional postprocessing actions

if(True == args.kml):
            if(args.kml):
                if(os.path.exists('coverage.kml')):
                    os.unlink('coverage.kml')
            kmlFile = open('coverage.kml', 'w')
            kmlFile.write(kmlTemplate.format(kml))
            kmlFile.close()
    

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

            elif 'nfo' in file.name:
                i = ScoutPhotoS3GeoImage(file)
                i.parseCameraParams()
                i.parseImageParams()

                r = GeoReferencer(i)
                trans = r.transform()
                if(True == args.kml):
                    kml += getKmlPoly(trans)



# Tempfile location for scratch work when generating geotifs
tempFile = '/tmp/temp.tif'

# Define UTM 11N projection with WKT.  This would need to be generalized for
# data in other areas.  Test data is in 16N.
proj = 'PROJCS["NAD83 / Alaska Albers",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",55],PARAMETER["standard_parallel_2",65],PARAMETER["latitude_of_center",50],PARAMETER["longitude_of_center",-154],PARAMETER["false_easting",0],PARAMETER["false_northing",0],AUTHORITY["EPSG","3338"],AXIS["X",EAST],AXIS["Y",NORTH]]'
"""
