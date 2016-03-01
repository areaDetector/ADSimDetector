< envPaths
errlogInit(20000)

dbLoadDatabase("$(TOP)/dbd/ADCSimDetectorApp.dbd")
ADCSimDetectorApp_registerRecordDeviceDriver(pdbbase) 

# Prefix for all records
epicsEnvSet("PREFIX", "13ADCSIM1:")
# The port name for the detector
epicsEnvSet("PORT",   "SIM1")
# The queue size for all plugins
epicsEnvSet("QSIZE",  "20")
# The maximim image width; used for row profiles in the NDPluginStats plugin
epicsEnvSet("XSIZE",  "8")
# The maximim image height; used for column profiles in the NDPluginStats plugin
epicsEnvSet("YSIZE",  "2000")
# The maximum number of time series points in the NDPluginStats plugin
epicsEnvSet("NCHANS", "2048")
# The maximum number of time series points in the NDPluginTimeSeries plugin
epicsEnvSet("TSPOINTS", "2048")
# The maximum number of frames buffered in the NDPluginCircularBuff plugin
epicsEnvSet("CBUFFS", "500")
# The search path for database files
epicsEnvSet("EPICS_DB_INCLUDE_PATH", "$(ADCORE)/db")

asynSetMinTimerPeriod(0.001)

# The EPICS environment variable EPICS_CA_MAX_ARRAY_BYTES needs to be set to a value at least as large
# as the largest image that the standard arrays plugin will send.
# That vlaue is $(XSIZE) * $(YSIZE) * sizeof(FTVL data type) for the FTVL used when loading the NDStdArrays.template file.
# The variable can be set in the environment before running the IOC or it can be set here.
# It is often convenient to set it in the environment outside the IOC to the largest array any client 
# or server will need.  For example 10000000 (ten million) bytes may be enough.
# If it is set here then remember to also set it outside the IOC for any CA clients that need to access the waveform record.  
# Do not set EPICS_CA_MAX_ARRAY_BYTES to a value much larger than that required, because EPICS Channel Access actually
# allocates arrays of this size every time it needs a buffer larger than 16K.
# Uncomment the following line to set it in the IOC.
#epicsEnvSet("EPICS_CA_MAX_ARRAY_BYTES", "10000000")

# Create an ADCimDetector driver
# ADCSimDetectorConfig(const char *portName, int numTimePoints, int dataType,
#                      int maxBuffers, int maxMemory, int priority, int stackSize)
ADCSimDetectorConfig("$(PORT)", $(YSIZE), 7, 0, 0)
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetector.template",  "P=$(PREFIX),R=det1:,  PORT=$(PORT),ADDR=0,TIMEOUT=1")
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetectorN.template", "P=$(PREFIX),R=det1:1:,PORT=$(PORT),ADDR=0,TIMEOUT=1")
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetectorN.template", "P=$(PREFIX),R=det1:2:,PORT=$(PORT),ADDR=1,TIMEOUT=1")
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetectorN.template", "P=$(PREFIX),R=det1:3:,PORT=$(PORT),ADDR=2,TIMEOUT=1")
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetectorN.template", "P=$(PREFIX),R=det1:4:,PORT=$(PORT),ADDR=3,TIMEOUT=1")
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetectorN.template", "P=$(PREFIX),R=det1:5:,PORT=$(PORT),ADDR=4,TIMEOUT=1")
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetectorN.template", "P=$(PREFIX),R=det1:6:,PORT=$(PORT),ADDR=5,TIMEOUT=1")
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetectorN.template", "P=$(PREFIX),R=det1:7:,PORT=$(PORT),ADDR=6,TIMEOUT=1")
dbLoadRecords("$(ADEXAMPLE)/db/ADCSimDetectorN.template", "P=$(PREFIX),R=det1:8:,PORT=$(PORT),ADDR=7,TIMEOUT=1")

# Create a standard arrays plugin, set it to get data from ADCSDetector driver.
NDStdArraysConfigure("Image1", 3, 0, "$(PORT)", 0)
# This creates a waveform large enough for 100000x8 arrays.
dbLoadRecords("NDStdArrays.template", "P=$(PREFIX),R=image1:,PORT=Image1,ADDR=0,TIMEOUT=1,NDARRAY_PORT=$(PORT),TYPE=Float64,FTVL=DOUBLE,NELEMENTS=800000")

# Time series plugin
NDTimeSeriesConfigure("TS1", $(QSIZE), 0, "$(PORT)", 0, 8, SIM_TIME_STEP)
dbLoadRecords("$(ADCORE)/db/NDTimeSeries.template",  "P=$(PREFIX),R=TS:,   PORT=TS1,ADDR=0,TIMEOUT=1,NDARRAY_PORT=$(PORT),NDARRAY_ADDR=0,NCHANS=$(TSPOINTS),ENABLED=1")
dbLoadRecords("$(ADCORE)/db/NDTimeSeriesN.template", "P=$(PREFIX),R=TS:1:, PORT=TS1,ADDR=0,TIMEOUT=1,NCHANS=$(TSPOINTS)")
dbLoadRecords("$(ADCORE)/db/NDTimeSeriesN.template", "P=$(PREFIX),R=TS:2:, PORT=TS1,ADDR=1,TIMEOUT=1,NCHANS=$(TSPOINTS)")
dbLoadRecords("$(ADCORE)/db/NDTimeSeriesN.template", "P=$(PREFIX),R=TS:3:, PORT=TS1,ADDR=2,TIMEOUT=1,NCHANS=$(TSPOINTS)")
dbLoadRecords("$(ADCORE)/db/NDTimeSeriesN.template", "P=$(PREFIX),R=TS:4:, PORT=TS1,ADDR=3,TIMEOUT=1,NCHANS=$(TSPOINTS)")
dbLoadRecords("$(ADCORE)/db/NDTimeSeriesN.template", "P=$(PREFIX),R=TS:5:, PORT=TS1,ADDR=4,TIMEOUT=1,NCHANS=$(TSPOINTS)")
dbLoadRecords("$(ADCORE)/db/NDTimeSeriesN.template", "P=$(PREFIX),R=TS:6:, PORT=TS1,ADDR=5,TIMEOUT=1,NCHANS=$(TSPOINTS)")
dbLoadRecords("$(ADCORE)/db/NDTimeSeriesN.template", "P=$(PREFIX),R=TS:7:, PORT=TS1,ADDR=6,TIMEOUT=1,NCHANS=$(TSPOINTS)")
dbLoadRecords("$(ADCORE)/db/NDTimeSeriesN.template", "P=$(PREFIX),R=TS:8:, PORT=TS1,ADDR=7,TIMEOUT=1,NCHANS=$(TSPOINTS)")

## Load all other plugins using commonPlugins.cmd
< $(ADCORE)/iocBoot/commonPlugins.cmd
set_requestfile_path("$(ADEXAMPLE)/exampleApp/Db")

#asynSetTraceIOMask("$(PORT)",0,2)
#asynSetTraceMask("$(PORT)",0,255)

iocInit()

# save things every thirty seconds
create_monitor_set("auto_settings.req", 30, "P=$(PREFIX)")
