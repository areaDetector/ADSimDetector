ADExample Releases
===============

The latest untagged master branch can be obtained at
https://github.com/areaDetector/ADExample.

Tagged source code and pre-built binary releases prior to R2-0 are included
in the areaDetector releases available via links at
http://cars.uchicago.edu/software/epics/areaDetector.html.

Tagged source code releases from R2-0 onward can be obtained at 
https://github.com/areaDetector/ADExample/releases.

Tagged prebuilt binaries from R2-0 onward can be obtained at
http://cars.uchicago.edu/software/pub/ADExample.

The versions of EPICS base, asyn, and other synApps modules used for each release can be obtained from 
the EXAMPLE_RELEASE_PATHS.local, EXAMPLE_RELEASE_LIBS.local, and EXAMPLE_RELEASE_PRODS.local
files respectively, in the configure/ directory of the appropriate release of the 
[top-level areaDetector](https://github.com/areaDetector/areaDetector) repository.


Release Notes
=============

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
