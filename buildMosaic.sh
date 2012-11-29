#! /usr/bin/env bash
for i in `seq 1 26`;
  do
    if [ $i -lt 10 ]
    then
      i="0"$i
    fi
    ls | grep "2012-11-06-flight_0$i.*tiff" | xargs gdal_merge.py -v -n 0 -o 2012-11-06-flight_0$i.mosaic.tiff
  done