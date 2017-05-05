ADSimDetector Releases
======================

The latest untagged master branch can be obtained at
https://github.com/areaDetector/ADSimDetector.

Tagged source code and pre-built binary releases prior to R2-0 are included
in the areaDetector releases available via links at
http://cars.uchicago.edu/software/epics/areaDetector.html.

Prior to R2-3 the files in this repository were contained in the ADExample
repository.  ADSimDetector was split into its own repository on Nov. 22, 2016
and the first release of the new repository was called R2-3 to continue the
release numbering from ADExample.

Tagged source code releases from R2-3 onward can be obtained at 
https://github.com/areaDetector/ADSimDetector/releases.

Tagged source code releases from R2-0 to R2-2 can be obtained at 
https://github.com/areaDetector/ADExample/releases.

Tagged prebuilt binaries from R2-3 onward can be obtained at
http://cars.uchicago.edu/software/pub/ADSimDetector.

Tagged prebuilt binaries from R2-0 to R2-2 onward can be obtained at
http://cars.uchicago.edu/software/pub/ADExample.

The versions of EPICS base, asyn, and other synApps modules used for each release can be obtained from 
the EXAMPLE_RELEASE_PATHS.local, EXAMPLE_RELEASE_LIBS.local, and EXAMPLE_RELEASE_PRODS.local
files respectively, in the configure/ directory of the appropriate release of the 
[top-level areaDetector](https://github.com/areaDetector/areaDetector) repository.


Release Notes
=============

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
[areaDetector Release Notes](http://cars.uchicago.edu/software/epics/areaDetectorReleaseNotes.html).
