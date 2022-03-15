TOP=../..
include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE

#======== DETECTOR LIBRARY ==============

# The following gets rid of the -fno-implicit-templates flag on vxWorks, 
# so we get automatic template instantiation.
# This is what we want for miscellaneous/asynPortDriver.cpp
ifeq (vxWorks,$(findstring vxWorks, $(T_A)))
CODE_CXXFLAGS=
endif

# Regenerate param set files from PVI YAML files if RULES_PVI exists
PVI_PARAMSETS=YES
PVI_YAML_PATHS = $(ADCORE)/pvi $(TOP)/pvi
-include $(PVI)/configure/RULES_PVI

INC += simDetectorParamSet.h
INC += simDetector.h

LIBRARY_IOC = simDetector
LIB_SRCS += simDetector.cpp

DBD += simDetectorSupport.dbd

include $(ADCORE)/ADApp/commonLibraryMakefile

#=============================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE

