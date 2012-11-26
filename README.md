# georectify-from-gps

Project:  Aeryon Scout data processing

Purpose:  Reads images from the Aeryon scout system and geocodes them to UTM projection

Author:   Scott Arko

This software was designed to read a variety of formats that can come from the Scout.  However, only the tif version currently works.  The dng import should work, however, both GDAL and PIL read the dng files at half resolution.  Once this is resolved, the workflow for tif and dng are identical. 

## Usage

1.  You need to have the osgeo, PIL, pyproj and numpy python modules installed.  matplotlib is required if you want to use the debug mode to inspect the images during processing. 

2.  Neither gdal (osgeo) or PIL read the Scout dng images properly.  They read as half resolution and I can't figure out why.  Therefore, I convert the dng files to tifs before processing.  I use ufraw-batch and ImageMagick to accomplish this.  Basically,

```ufraw-batch *.dng```

This will convert all the dng files to ppm, then

```mogrify --auto-levels -format tif *.ppm```

This will create tif files with the same basename as the dng and nfo files.  After this, I create a new directory and move the tif files into it.  Then copy the nfo files into the same directory as the tifs.  

Then you can run the script.  

```./importScout.py [-d] file file file ... file```

I would recommend running it with one file to begin just to see if it works.  It will make a system call to gdalwarp to do the final warping, but if you have osgeo python installed, you also should have gdalwarp.

## License

Copyright (c) 2012, Scott Arko 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.

You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the
Free Software Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.