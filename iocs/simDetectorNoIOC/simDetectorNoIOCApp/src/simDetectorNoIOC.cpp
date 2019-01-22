/* simDetectorNoIOC.cpp
 *
 * This is an example of creating a simDetector and controlling it from outside an IOC
 *
 * Author: Mark Rivers
 *         University of Chicago
 *
 * Created:  October 27, 2014
 *
 */

#include <epicsThread.h>
#include <asynPortClient.h>
#include <NDPluginStats.h>
#include <NDPluginStdArrays.h>
#include <NDFileHDF5.h>
#include <simDetector.h>

#define NUM_FRAMES 10
#define FILE_PATH "/tmp/"
#define REPORT_LEVEL 10

#ifndef EPICS_LIBCOM_ONLY
  #include <dbAccess.h>
#endif

class callbackStruct {
public:
  callbackStruct(const char *paramStringIn, asynParamClient *pClientIn, class simDetectorDemo *pSimDetectorDemoIn)
  : paramString(paramStringIn), pClient(pClientIn), pSimDetectorDemo(pSimDetectorDemoIn), pasynUser(0) {};
  const char *paramString;
  asynParamClient *pClient;
  class simDetectorDemo *pSimDetectorDemo;
  asynUser *pasynUser;
};


class epicsShareClass simDetectorDemo
{
public:
  simDetectorDemo();
  ~simDetectorDemo();
  void testAcquire();
  void int32Callback(callbackStruct *pCallback, epicsInt32 data);
  void float64Callback(callbackStruct *pCallback, epicsFloat64 data);
  void NDArrayCallback(callbackStruct *pCallback, NDArray *pArray);

private:
  simDetector    *pSimDetector_;
  asynPortClient *pSimClient_;
  NDPluginStats  *pStatsPlugin_;
  asynPortClient *pStatsClient_;
  NDFileHDF5     *pHDF5Plugin_;
  asynPortClient *pHDF5Client_;
  bool isAcquiring_;
};

void int32CallbackC(void *drvPvt, asynUser *pasynUser, epicsInt32 data)
{
  callbackStruct *pCallback = (callbackStruct*)drvPvt;
  pCallback->pasynUser = pasynUser;
  pCallback->pSimDetectorDemo->int32Callback(pCallback, data);
}

void simDetectorDemo::int32Callback(callbackStruct *pCallback, epicsInt32 data)
{
  if (strcmp(pCallback->paramString, ADAcquireString) == 0) {
    isAcquiring_ = data ? true : false;
  }
}

void float64CallbackC(void *drvPvt, asynUser *pasynUser, epicsFloat64 data)
{
  callbackStruct *pCallback = (callbackStruct*)drvPvt;
  pCallback->pasynUser = pasynUser;
  pCallback->pSimDetectorDemo->float64Callback(pCallback, data);
}


void simDetectorDemo::float64Callback(callbackStruct *pCallback, epicsFloat64 data)
{
  if (strcmp(pCallback->paramString, NDPluginStatsMeanValueString) == 0) {
    printf("Statistics: mean value=%f\n", data);
  }
}

void NDArrayCallbackC(void *drvPvt, asynUser *pasynUser, void* pData)
{
  callbackStruct *pCallback = (callbackStruct*)drvPvt;
  pCallback->pasynUser = pasynUser;
  pCallback->pSimDetectorDemo->NDArrayCallback(pCallback, (NDArray *)pData);
}


void simDetectorDemo::NDArrayCallback(callbackStruct *pCallback, NDArray *pArray)
{
  if (strcmp(pCallback->paramString, NDArrayDataString) == 0) {
    printf("NDArray.uniqueId=%d\n", pArray->uniqueId);
  }
}

simDetectorDemo::simDetectorDemo()
{

  // Create a simDetector driver
  pSimDetector_ =  new simDetector("SIM1", 1024, 1024, NDUInt8, 0, 0, 0, 0);
  // Create an asynPortClient for the simDetector
  pSimClient_   =  new asynPortClient("SIM1");
  pSimClient_->write(NDArrayCallbacksString, 1);           // Enable NDArray callbacks
  pSimClient_->write(ADGainString, 1.0);                   // Set the Gain to 1.0
  pSimClient_->write(ADImageModeString, ADImageMultiple);  // Set the image mode to Multiple
  pSimClient_->write(SimNoiseString, 50.0);                // Set random noise to 50 counts
  pSimClient_->write(ADAcquireTimeString, 0.2);            // 0.2 seconds per image
  pSimClient_->write(ADNumImagesString, NUM_FRAMES);       // Collect 10 frames
  pSimClient_->write(ADNumImagesString, NUM_FRAMES);       // Collect 10 frames

  // Create a statistics plugin getting its data from the simDetector
  pStatsPlugin_ =  new NDPluginStats("STATS1", 20, 0, "SIM1", 0, 0, 0, 0, 0);
  pStatsClient_ =  new asynPortClient("STATS1"); 
  pStatsPlugin_->start();  // Start the plugin
  pStatsClient_->write(NDPluginDriverEnableCallbacksString, 1);  // Enable callbacks to this plugin
  pStatsClient_->write(NDPluginStatsComputeStatisticsString, 1); // Enable computing basic statistics

  // Create an HDF5 plugin getting its data from the simDetector
  pHDF5Plugin_  = new NDFileHDF5("HDF5",   20, 0, "SIM1", 0, 0, 0);
  pHDF5Client_  = new asynPortClient("HDF5"); 
  pHDF5Plugin_->start();  // Start the plugin
  pHDF5Client_->write(NDPluginDriverEnableCallbacksString, 1);  // Enable callbacks to this plugin
  pHDF5Client_->write(NDFileNameString, "test");                // Set the file name
  pHDF5Client_->write(NDFilePathString, FILE_PATH);             // Set the file path
  pHDF5Client_->write(NDFileNumberString, 1);                   // Set the file number
  pHDF5Client_->write(NDAutoIncrementString, 1);                // Enable file number auto-increment
  pHDF5Client_->write(NDFileTemplateString, "%s%s_%3.3d.h5");   // Set the file name format string (C-style)
  pHDF5Client_->write(NDFileWriteModeString, NDFileModeStream); // Set the file mode to stream (multiple arrays per file)
  pHDF5Client_->write(NDFileNumCaptureString, NUM_FRAMES);      // Number of arrays to stream before closing file 
  pHDF5Client_->write(NDFileLazyOpenString, 1);                 // Wait to open file till first frame arrives 

  // Enable callbacks for some parameters
  asynInt32Client *pAcquire = (asynInt32Client*)pSimClient_->getParamClient(ADAcquireString);
  callbackStruct *pAcquireCallback = new callbackStruct(ADAcquireString, pAcquire, this);
  pAcquire->registerInterruptUser(int32CallbackC, pAcquireCallback);

  asynFloat64Client *pMean = (asynFloat64Client*)pStatsClient_->getParamClient(NDPluginStatsMeanValueString);
  callbackStruct *pMeanCallback = new callbackStruct(NDPluginStatsMeanValueString, pMean, this);
  pMean->registerInterruptUser(float64CallbackC, pMeanCallback);

  asynGenericPointerClient *pNDArray = (asynGenericPointerClient*)pSimClient_->getParamClient(NDArrayDataString);
  callbackStruct *pNDArrayCallback = new callbackStruct(NDArrayDataString, pNDArray, this);
  pNDArray->registerInterruptUser(NDArrayCallbackC, pNDArrayCallback);

}

simDetectorDemo::~simDetectorDemo()
{}

void simDetectorDemo::testAcquire()
{
  // Start the HDF5 plugin streaming
  pHDF5Client_->write(NDFileCaptureString, 1);

  // Start the simDetector acquiring
  pSimClient_->write(ADAcquireString, 1);

  // Wait for acquisition to complete. isAcquiring_ is cleared in the int32Callback function
  while (isAcquiring_) {
    epicsThreadSleep(0.1);
  }
  int numHDF5Captured;
  pHDF5Client_->read(NDFileNumCapturedString, &numHDF5Captured);
  printf("HDF5 numCaptured = %d\n", numHDF5Captured);
}

int main(int argc, char **argv)
{
#ifndef EPICS_LIBCOM_ONLY
  // Must set this for callbacks to work if EPICS_LIBCOM_ONLY is not defined
  interruptAccept = 1;  
#endif
  // Create the object and acquire 3 times
  simDetectorDemo demo;
  demo.testAcquire();
  demo.testAcquire();
  demo.testAcquire();
}



