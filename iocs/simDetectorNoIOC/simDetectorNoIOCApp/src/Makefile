TOP=../..
include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE

# The following gets rid of the -fno-implicit-templates flag on vxWorks, 
# so we get automatic template instantiation.
# This is what we want for miscellaneous/asynPortDriver.cpp
ifeq (vxWorks,$(findstring vxWorks, $(T_A)))
CODE_CXXFLAGS=
endif

ifeq ($(EPICS_LIBCOM_ONLY),YES)
  USR_CXXFLAGS += -DEPICS_LIBCOM_ONLY
endif

PROD_IOC_Linux  += simDetectorNoIOCApp
PROD_IOC_WIN32  += simDetectorNoIOCApp
PROD_IOC_Darwin += simDetectorNoIOCApp
PROD_SRCS += simDetectorNoIOC.cpp

PROD_LIBS += simDetector

include $(ADCORE)/ADApp/commonDriverMakefile

ifdef HDF5_INCLUDE
  USR_INCLUDES += $(addprefix -I $(HDF5_INCLUDE))
endif
ifdef SZIP_INCLUDE
  USR_INCLUDES += $(addprefix -I $(SZIP_INCLUDE))
endif
ifdef XML2_INCLUDE
  USR_INCLUDES += $(addprefix -I $(XML2_INCLUDE))
endif

#=============================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE

