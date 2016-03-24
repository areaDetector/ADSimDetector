/* ADCSimDetector.cpp
 *
 * This is a driver for a simulated ADC detector.
 *
 * Author: Mark Rivers
 *         University of Chicago
 *
 * Created:  February 28, 2016
 *
 */

#include <stdlib.h>
#include <math.h>

#include <epicsTime.h>
#include <epicsThread.h>
#include <epicsEvent.h>
#include <iocsh.h>

#include "asynNDArrayDriver.h"
#include <epicsExport.h>
#include "ADCSimDetector.h"

static const char *driverName = "ADCSimDetector";

// Some systems don't define M_PI in math.h
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

static void simTaskC(void *drvPvt)
{
    ADCSimDetector *pPvt = (ADCSimDetector *)drvPvt;

    pPvt->simTask();
}

/** Constructor for ADCSimDetector; most parameters are simply passed to ADDriver::ADDriver.
  * After calling the base class constructor this method creates a thread to compute the simulated detector data,
  * and sets reasonable default values for parameters defined in this class, asynNDArrayDriver and ADDriver.
  * \param[in] portName The name of the asyn port driver to be created.
  * \param[in] numTimePoints The initial number of time points.
  * \param[in] dataType The initial data type (NDDataType_t) of the arrays that this driver will create.
  * \param[in] maxBuffers The maximum number of NDArray buffers that the NDArrayPool for this driver is
  *            allowed to allocate. Set this to -1 to allow an unlimited number of buffers.
  * \param[in] maxMemory The maximum amount of memory that the NDArrayPool for this driver is
  *            allowed to allocate. Set this to -1 to allow an unlimited amount of memory.
  * \param[in] priority The thread priority for the asyn port driver thread if ASYN_CANBLOCK is set in asynFlags.
  * \param[in] stackSize The stack size for the asyn port driver thread if ASYN_CANBLOCK is set in asynFlags.
  */
ADCSimDetector::ADCSimDetector(const char *portName, int numTimePoints, NDDataType_t dataType,
                               int maxBuffers, size_t maxMemory, int priority, int stackSize)

    : asynNDArrayDriver(portName, MAX_SIGNALS, NUM_SIM_DETECTOR_PARAMS, maxBuffers, maxMemory,
               0, 0, /* No interfaces beyond those set in ADDriver.cpp */
               ASYN_CANBLOCK | ASYN_MULTIDEVICE, /* asyn flags*/
               1,                                /* autoConnect=1 */
               priority, stackSize),
    uniqueId_(0), acquiring_(0)

{
    int status = asynSuccess;
    const char *functionName = "ADCSimDetector";

    /* Create the epicsEvents for signaling to the simulate task when acquisition starts and stops */
    this->startEventId_ = epicsEventCreate(epicsEventEmpty);
    if (!this->startEventId_) {
        printf("%s:%s epicsEventCreate failure for start event\n",
            driverName, functionName);
        return;
    }
    this->stopEventId_ = epicsEventCreate(epicsEventEmpty);
    if (!this->stopEventId_) {
        printf("%s:%s epicsEventCreate failure for stop event\n",
            driverName, functionName);
        return;
    }

    createParam(SimAcquireString,         asynParamInt32, &P_Acquire);
    createParam(SimAcquireTimeString,   asynParamFloat64, &P_AcquireTime);
    createParam(SimElapsedTimeString,   asynParamFloat64, &P_ElapsedTime);
    createParam(SimTimeStepString,      asynParamFloat64, &P_TimeStep);
    createParam(SimNumTimePointsString,   asynParamInt32, &P_NumTimePoints);
    createParam(SimPeriodString,        asynParamFloat64, &P_Period);
    createParam(SimAmplitudeString,     asynParamFloat64, &P_Amplitude);
    createParam(SimOffsetString,        asynParamFloat64, &P_Offset);
    createParam(SimFrequencyString,     asynParamFloat64, &P_Frequency);
    createParam(SimPhaseString,         asynParamFloat64, &P_Phase);
    createParam(SimNoiseString,         asynParamFloat64, &P_Noise);

    status |= setIntegerParam(P_NumTimePoints, numTimePoints);
    status |= setIntegerParam(NDDataType, dataType);
    status |= setDoubleParam(P_TimeStep, 0.001);
    status |= setDoubleParam(P_Amplitude, 1.0);
    status |= setDoubleParam(P_Offset,    0.0);
    status |= setDoubleParam(P_Period,    1.0);
    status |= setDoubleParam(P_Phase,     0.0);
    status |= setDoubleParam(P_Noise,     0.0);

    if (status) {
        printf("%s: unable to set parameters\n", functionName);
        return;
    }

    /* Create the thread that updates the images */
    status = (epicsThreadCreate("SimDetTask",
                                epicsThreadPriorityMedium,
                                epicsThreadGetStackSize(epicsThreadStackMedium),
                                (EPICSTHREADFUNC)simTaskC,
                                this) == NULL);
    if (status) {
        printf("%s:%s epicsThreadCreate failure for simulation task\n",
            driverName, functionName);
        return;
    }
}

/** Template function to compute the simulated detector data for any data type */
template <typename epicsType> void ADCSimDetector::computeArraysT()
{
    size_t dims[2];
    int numTimePoints;
    int i, j;
    NDDataType_t dataType;
    epicsType *pData;
    double acquireTime;
    double timeStep;
    double rndm;
    double amplitude[MAX_SIGNALS];
    double period[MAX_SIGNALS];
    double frequency[MAX_SIGNALS];
    double phase[MAX_SIGNALS];
    double noise[MAX_SIGNALS];
    double offset[MAX_SIGNALS];
    
    getIntegerParam(NDDataType, (int *)&dataType);
    getIntegerParam(P_NumTimePoints, &numTimePoints);
    getDoubleParam(P_TimeStep, &timeStep);
    getDoubleParam(P_AcquireTime, &acquireTime);

    dims[0] = MAX_SIGNALS;
    dims[1] = numTimePoints;

    if (this->pArrays[0]) this->pArrays[0]->release();
    this->pArrays[0] = pNDArrayPool->alloc(2, dims, dataType, 0, 0);
    pData = (epicsType *)this->pArrays[0]->pData;
    memset(pData, 0, MAX_SIGNALS * numTimePoints * sizeof(epicsType));

    for (j=0; j<MAX_SIGNALS; j++) {
        getDoubleParam(j, P_Amplitude, amplitude+j);
        getDoubleParam(j, P_Offset,    offset+j);
        getDoubleParam(j, P_Period,    period+j);
        frequency[j] = 1. / period[j];
        setDoubleParam(j, P_Frequency, frequency[j]);
        getDoubleParam(j, P_Phase,     phase+j);
        getDoubleParam(j, P_Noise,     noise+j);
        phase[j] = phase[j]/360.0;
    }
    for (i=0; i<numTimePoints; i++) {
        rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
        // Signal 0 is a sin wave
        j = 0;
        pData[MAX_SIGNALS*i + j] = (epicsType)(offset[j] + noise[j] * rndm + amplitude[j] * 
                                               sin((elapsedTime_ * frequency[j] + phase[j]) * 2. * M_PI));
        // Signal 1 is a cos wave
        j = 1;
        pData[MAX_SIGNALS*i + j] = (epicsType)(offset[j] + noise[j] * rndm + amplitude[j] * 
                                               cos((elapsedTime_ * frequency[j] + phase[j]) * 2. * M_PI));
        // Signal 2 is a square wave
        j = 2;
        pData[MAX_SIGNALS*i + j] = (epicsType)(offset[j] + noise[j] * rndm + amplitude[j] * 
                                              (sin((elapsedTime_ * frequency[j] + phase[j]) * 2. * M_PI) > 0 ? 1.0 : -1.0));
        // Signal 3 is a sawtooth
        j = 3;
        pData[MAX_SIGNALS*i + j] = (epicsType)(offset[j] + noise[j] * rndm + amplitude[j] * 
                                               -2.0/M_PI * atan(1./tan((elapsedTime_ * frequency[j] + phase[j]) * M_PI)));
        // Signal 4 is white noise
        j = 4;
        pData[MAX_SIGNALS*i + j] = (epicsType)(offset[j] + noise[j] * rndm + amplitude[j] * rndm);

        // Signal 5 is signal 0 + signal 1
        j = 5;
        pData[MAX_SIGNALS*i + j] = (epicsType)(offset[j] + noise[j] * rndm + amplitude[j] * 
                                               pData[MAX_SIGNALS*i + 0] + pData[MAX_SIGNALS*i + 1]) ;

        // Signal 6 is signal 0 * signal 1
        j = 6;
        pData[MAX_SIGNALS*i + j] = (epicsType)(offset[j] + noise[j] * rndm + amplitude[j] * 
                                               pData[MAX_SIGNALS*i + 0] * pData[MAX_SIGNALS*i + 1]) ;

        // Signal 7 is 4 sin waves
        j = 7;
        pData[MAX_SIGNALS*i + j] = (epicsType)(offset[j] + noise[j] * rndm + amplitude[j] *
                                              (sin((elapsedTime_ * 1.*frequency[j] + phase[j]) * 2. * M_PI) +
                                               sin((elapsedTime_ * 2.*frequency[j] + phase[j]) * 2. * M_PI) +
                                               sin((elapsedTime_ * 3.*frequency[j] + phase[j]) * 2. * M_PI) +
                                               sin((elapsedTime_ * 4.*frequency[j] + phase[j]) * 2. * M_PI)));

        elapsedTime_ += timeStep;
        if ((acquireTime > 0) && (elapsedTime_ > acquireTime)) {
            setAcquire(0);
            setIntegerParam(P_Acquire, 0);
            break;
        }
    }
    setDoubleParam(P_ElapsedTime, elapsedTime_);
}

/** Computes the new image data */
void ADCSimDetector::computeArrays()
{
    int dataType;
    getIntegerParam(NDDataType, &dataType); 

    switch (dataType) {
        case NDInt8:
            computeArraysT<epicsInt8>();
            break;
        case NDUInt8:
            computeArraysT<epicsUInt8>();
            break;
        case NDInt16:
            computeArraysT<epicsInt16>();
            break;
        case NDUInt16:
            computeArraysT<epicsUInt16>();
            break;
        case NDInt32:
            computeArraysT<epicsInt32>();
            break;
        case NDUInt32:
            computeArraysT<epicsUInt32>();
            break;
        case NDFloat32:
            computeArraysT<epicsFloat32>();
            break;
        case NDFloat64:
            computeArraysT<epicsFloat64>();
            break;
    }
}

void ADCSimDetector::setAcquire(int value)
{
    if (value && !acquiring_) {
        /* Send an event to wake up the simulation task */
        epicsEventSignal(this->startEventId_); 
    }
    if (!value && acquiring_) {
        /* This was a command to stop acquisition */
        /* Send the stop event */
        epicsEventSignal(this->stopEventId_); 
    }
}

/** This thread calls computeImage to compute new image data and does the callbacks to send it to higher layers.
  * It implements the logic for single, multiple or continuous acquisition. */
void ADCSimDetector::simTask()
{
    int status = asynSuccess;
    NDArray *pImage;
    epicsTimeStamp startTime;
    int numTimePoints;
    int arrayCounter;
    double timeStep;
    int i;
    const char *functionName = "simTask";

    this->lock();
    /* Loop forever */
    while (1) {
        /* Has acquisition been stopped? */
        status = epicsEventTryWait(this->stopEventId_);
        if (status == epicsEventWaitOK) {
            acquiring_ = 0;
        }
       
        /* If we are not acquiring then wait for a semaphore that is given when acquisition is started */
        if (!acquiring_) {
          /* Release the lock while we wait for an event that says acquire has started, then lock again */
            asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                "%s:%s: waiting for acquire to start\n", driverName, functionName);
            this->unlock();
            status = epicsEventWait(this->startEventId_);
            this->lock();
            acquiring_ = 1;
            elapsedTime_ = 0.0;
        }

        /* Update the data */
        computeArrays();

        pImage = this->pArrays[0];

        /* Put the frame number and time stamp into the buffer */
        pImage->uniqueId = uniqueId_++;
        getIntegerParam(NDArrayCounter, &arrayCounter);
        arrayCounter++;
        setIntegerParam(NDArrayCounter, arrayCounter);
        epicsTimeGetCurrent(&startTime);
        pImage->timeStamp = startTime.secPastEpoch + startTime.nsec / 1.e9;
        updateTimeStamp(&pImage->epicsTS);

        /* Get any attributes that have been defined for this driver */
        this->getAttributes(pImage->pAttributeList);

        /* Call the NDArray callback */
        /* Must release the lock here, or we can get into a deadlock, because we can
         * block on the plugin lock, and the plugin can be calling us */
        this->unlock();
        doCallbacksGenericPointer(pImage, NDArrayData, 0);
        this->lock();

        /* Call the callbacks to update any changes */
        for (i=0; i<MAX_SIGNALS; i++) {
            callParamCallbacks(i);
        }

        /* Sleep for numTimePoint * timeStep seconds */
        getIntegerParam(P_NumTimePoints, &numTimePoints);
        getDoubleParam(P_TimeStep, &timeStep);
        unlock();
        epicsThreadSleep(numTimePoints * timeStep);
        lock();
    }
}


/** Called when asyn clients call pasynInt32->write().
  * This function performs actions for some parameters, including ADAcquire, ADColorMode, etc.
  * For all parameters it sets the value in the parameter library and calls any registered callbacks..
  * \param[in] pasynUser pasynUser structure that encodes the reason and address.
  * \param[in] value Value to write. */
asynStatus ADCSimDetector::writeInt32(asynUser *pasynUser, epicsInt32 value)
{
    int function = pasynUser->reason;
    int addr;
    asynStatus status = asynSuccess;

    getAddress(pasynUser, &addr);
 
    /* Set the parameter and readback in the parameter library.  This may be overwritten when we read back the
     * status at the end, but that's OK */
    status = setIntegerParam(addr, function, value);

    if (function == P_Acquire) {
        setAcquire(value);
    } else {
        /* If this parameter belongs to a base class call its method */
        if (function < FIRST_SIM_DETECTOR_PARAM) status = asynNDArrayDriver::writeInt32(pasynUser, value);
    }

    /* Do callbacks so higher layers see any changes */
    callParamCallbacks(addr);

    if (status)
        asynPrint(pasynUser, ASYN_TRACE_ERROR,
              "%s:writeInt32 error, status=%d function=%d, value=%d\n",
              driverName, status, function, value);
    else
        asynPrint(pasynUser, ASYN_TRACEIO_DRIVER,
              "%s:writeInt32: function=%d, value=%d\n",
              driverName, function, value);
    return status;
}



/** Report status of the driver.
  * Prints details about the driver if details>0.
  * It then calls the ADDriver::report() method.
  * \param[in] fp File pointed passed by caller where the output is written to.
  * \param[in] details If >0 then driver details are printed.
  */
void ADCSimDetector::report(FILE *fp, int details)
{

    fprintf(fp, "ADC simulation detector %s\n", this->portName);
    if (details > 0) {
        int numTimePoints, dataType;
        getIntegerParam(P_NumTimePoints, &numTimePoints);
        getIntegerParam(NDDataType, &dataType);
        fprintf(fp, "  # time points:   %d\n", numTimePoints);
        fprintf(fp, "      Data type:   %d\n", dataType);
    }
    /* Invoke the base class method */
    asynNDArrayDriver::report(fp, details);
}


/** Configuration command, called directly or from iocsh */
extern "C" int ADCSimDetectorConfig(const char *portName, int numTimePoints, int dataType,
                                 int maxBuffers, int maxMemory, int priority, int stackSize)
{
    new ADCSimDetector(portName, numTimePoints, (NDDataType_t)dataType,
                    (maxBuffers < 0) ? 0 : maxBuffers,
                    (maxMemory < 0) ? 0 : maxMemory, 
                    priority, stackSize);
    return(asynSuccess);
}

/** Code for iocsh registration */
static const iocshArg ADCSimDetectorConfigArg0 = {"Port name",     iocshArgString};
static const iocshArg ADCSimDetectorConfigArg1 = {"# time points", iocshArgInt};
static const iocshArg ADCSimDetectorConfigArg2 = {"Data type",     iocshArgInt};
static const iocshArg ADCSimDetectorConfigArg3 = {"maxBuffers",    iocshArgInt};
static const iocshArg ADCSimDetectorConfigArg4 = {"maxMemory",     iocshArgInt};
static const iocshArg ADCSimDetectorConfigArg5 = {"priority",      iocshArgInt};
static const iocshArg ADCSimDetectorConfigArg6 = {"stackSize",     iocshArgInt};
static const iocshArg * const ADCSimDetectorConfigArgs[] = {&ADCSimDetectorConfigArg0,
                                                            &ADCSimDetectorConfigArg1,
                                                            &ADCSimDetectorConfigArg2,
                                                            &ADCSimDetectorConfigArg3,
                                                            &ADCSimDetectorConfigArg4,
                                                            &ADCSimDetectorConfigArg5,
                                                            &ADCSimDetectorConfigArg6};
static const iocshFuncDef configADCSimDetector = {"ADCSimDetectorConfig", 7, ADCSimDetectorConfigArgs};
static void configADCSimDetectorCallFunc(const iocshArgBuf *args)
{
    ADCSimDetectorConfig(args[0].sval, args[1].ival, args[2].ival, args[3].ival,
                         args[4].ival, args[5].ival, args[6].ival);
}


static void ADCSimDetectorRegister(void)
{

    iocshRegister(&configADCSimDetector, configADCSimDetectorCallFunc);
}

extern "C" {
epicsExportRegistrar(ADCSimDetectorRegister);
}
