ADSimDetector Releases
======================

The latest untagged master branch can be obtained at
https://github.com/areaDetector/ADSimDetector.

Tagged source code and pre-built binary releases prior to R2-0 are included
in the areaDetector releases available via links at
https://cars.uchicago.edu/software/epics/areaDetector.html.

Prior to R2-3 the files in this repository were contained in the ADExample
repository.  ADSimDetector was split into its own repository on Nov. 22, 2016
and the first release of the new repository was called R2-3 to continue the
release numbering from ADExample.

Tagged source code releases from R2-3 onward can be obtained at 
https://github.com/areaDetector/ADSimDetector/releases.

Tagged source code releases from R2-0 to R2-2 can be obtained at 
https://github.com/areaDetector/ADExample/releases.

Tagged prebuilt binaries from R2-3 onward can be obtained at
https://cars.uchicago.edu/software/pub/ADSimDetector.

Tagged prebuilt binaries from R2-0 to R2-2 onward can be obtained at
https://cars.uchicago.edu/software/pub/ADExample.

The versions of EPICS base, asyn, and other synApps modules used for each release can be obtained from 
the EXAMPLE_RELEASE_PATHS.local, EXAMPLE_RELEASE_LIBS.local, and EXAMPLE_RELEASE_PRODS.local
files respectively, in the configure/ directory of the appropriate release of the 
[top-level areaDetector](https://github.com/areaDetector/areaDetector) repository.


Release Notes
=============

R2-11 (April 7, 2024)
=========================
* Moved documentation to https://areadetector.github.io/areaDetector/ADSimDetector/simDetector.html
* Added .bob OPI files for Phoebus.
* Use new asyn NDArrayDriver::updateTimeStamps method.
* Initialize pointer to NULL.
* Improve validation of Min/Size/Bin.


R2-10 (October 22, 2019)
=========================
* Added support for NDArray datatypes NDInt64 and NDUInt64
  * This requires asyn R4-37 and ADCore R3-8 
* Added iocSimDetector/testHDF5Compression.py to test HDF5 compression.
* Added Time64 attribute to iocSimDetector/simDetectorAttributes.xml.  This tests attributes of type NDAttrUInt64.


R2-9 (May 29, 2019)
===================
* Requires asyn R4-35 because of changes to the asynPortClient class used in simDetectorNoIOC (see below).
* Changes in the simulation modes and parameters to improve utility and performance.
  * The Noise parameter now applies to all simulation modes, not just to Peaks mode.
  * The Noise computation only calls rand() when the image is reset, i.e. the data type,
    dimensions, etc. are changed.  In normal cases it just adds the precomputed noise
    to each pixel starting at a random pixel.  This means the noise spatial pattern is the same
    in each successive image, but starting at a different pixel.  This makes it look quite random
    while being very fast.
  * A new Offset parameter now applies to all simulation modes.
  * The sineNoise and sineOffset parameters were removed, replaced by Noise and Offset.
  * PeakVariation was changed from longout record to ao record so it is not constrained to be an integer.
  * Removed the dependence of the image data on the AcquireTime.  This was confusing and meant the data
    changed significantly when the AcquireTime changed.
  * Improved performance of Peaks mode by more than 10X.  Previously it computed the Gaussian peak profile
    for each peak in the image.  Now it only does the Gausssian calculation once into a scratch buffer
    and then copies that peak multiple times into the output image at different locations.
  * A new simulation mode Offset&Noise was added.  This mode only does offset and noise, and
    is the fastest simulation mode.  It can generate >700 frames/s using 1024x1024 images and asynFloat64
    data type, which is >5GB/s.  In asynUInt8 mode it can generate >3000 frames/s which is >3GB/s.
* simDetectorNoIOC
  * This directory contains a pure C++ application that instantiates a simDetector, a statistics plugin,
    and an HDF5 file plugin outside the context of an EPICS IOC.  
    It is intended in part to demonstrate that areaDetector drivers and plugins can be used with other
    control systems such as Tango. 
  * This application has been considerably improved and simplified by using the new versions of the
    asynPortClient and asynParamClient classes in asyn R4-35.
* iocHDF5Test
  * Fixed the startup scripts.


R2-8 (July 1, 2018)
===================
* Changed Makefile to use addprefix to add -I to user-defined include file directory paths
* Fixed so it builds OK if WITH_HDF5=NO.
* Changed simDetector.adl to have the new records from ADCore R3-3.
* Updated the edl, ui, and opi autoconvert directories to contain the conversions
  from the most recent adl files and improved converters.
* Fixed comments in the files in the iocBoot directories.  Thanks to Lewis Muir for this.


R2-7 (January 28, 2018)
========================
* Fixed simDetectorNoIOCApp/src/Makefile so it handles external locations of HDF, SZIP, and XML2 correctly.
* Add simDetectorApp/op/Makefile to automatically convert adl files to edl, ui, and opi.
* Fixed medm adl files to improve the autoconversion to other display manager files.
* Added op/Makefile to automatically convert adl files to edl, ui, and opi files.
* Updated the edl, ui, and opi autoconvert directories to contain the conversions
  from the most recent adl files.


R2-6 (July 4, 2017)
========================
* Change startup scripts in iocSimDetector.  st.cmd now loads envPaths and st_base.cmd, so it can be
  run with no copying of files.  st.cmd.linux and st.cmd.windows can still be used for conveniently
  switching host architectures in the same directory by copying envPaths to envPaths.linux or envPaths.windows
  as appropriate.


R2-5 (May 5, 2017)
========================
* Removed parameter counting from driver, not needed with asyn R4-31.
* Added room for new PVs for attribute file status and macros in medm screen.
* Moved NDPluginFFT loading from this IOC to EXAMPLE_commonPlugins.cmd.
* Added new required maxBuffers argument when creating NDPluginStdArrays in simDetectorNoIOC.cpp.


R2-4 (February 21, 2017)
========================
* Add new parameters ADSerialNumber, ADFirmwareVersion, ADSDKVersion. Minor change to
  driver and medm screen.

  
R2-3 (November 22, 2016)
========================
* ADSimDetector was split into its own repository from ADExample.  The release numbering
  starts where ADExample left off.


R2-2 (October 28, 2016)
========================
* Created a new driver in ADExample/exampleApp/ADCSimDetectorSrc.
  This driver generates times-series data for 8 signals as a 2-D array [8, NumTimePoints].
  The signals are common waveforms (sine, cosine, square wave, sawtooth, random noise, etc.)
  
  The driver was written in part to test the new NDPluginTimeSeries and NDPluginFFT plugins 
  in ADCore R2-5. But it can also serve as a starting point for writing drivers that are used 
  to collect time-series data like ADCs, waveform digitizers, etc.  

  There is a new iocs/iocADCSimDetector directory that creates an application for testing
  the driver.
* Added a new sine-wave simulation mode to simDetector.  This simulation mode constructs images
  as the sum or products of two sine waves in the X and Y directions.  It was designed in part
  to test the new NDPluginFFT plugin in ADCore R2-5.
* Added new iocPvaDriver example IOC.  It implements the new pvaDriver from ADCore R2-5 which
  receives EPICS V4 NTNDArrays over EPICS V4 PVAccess, converts them to NDArrays and calls
  areaDetector plugins.


R2-1 (January 9, 2016)
========================
* Changed SchemaLocation in XML files to relative path to ADCore.
* Fixed bug in simDetector::writeInt32.  ImageMode could be initialized.
  Thanks to Phil Atkin for fixing this.
* Changed st.cmd to use FTVL=SHORT rather than USHORT for 16-bit images.
  This is more efficient, and works fine with ImageJ which always treats 16-bit images as unsigned.


R2-0-1 (September 23, 2015)
========================
Changed iocs/iocSimDetector*/configure/RELEASE, replacing ADEXAMPLE_TOP with
ADEXAMPLE.  The _TOP is not needed, and it broke the st.cmd IOC startup script.


R2-0 (September 18, 2015)
========================
This is the first release of this repository.  It contains the simDetector driver and
example IOCS.

The files in this this repository were previously located in the ADCore repository.


R1-9-1 and earlier
==================
Release notes are part of the
[areaDetector Release Notes](https://cars.uchicago.edu/software/epics/areaDetectorReleaseNotes.html).
