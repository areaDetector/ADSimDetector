/* ADCSimDetector.cpp
 *
 * This is a driver for a simulated area detector.
 *
 * Author: Mark Rivers
 *         University of Chicago
 *
 * Created:  March 20, 2008
 *
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include <epicsTime.h>
#include <epicsThread.h>
#include <epicsEvent.h>
#include <iocsh.h>

#include "asynNDArrayDriver.h"
#include <epicsExport.h>
#include "ADCSimDetector.h"

static const char *driverName = "ADCSimDetector";

/** Template function to compute the simulated detector data for any data type */
template <typename epicsType> void ADCSimDetector::computeArraysT()
{
    //int status = asynSuccess;
    //(epicsType *)this->pRaw->pData;
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

static void simTaskC(void *drvPvt)
{
    ADCSimDetector *pPvt = (ADCSimDetector *)drvPvt;

    pPvt->simTask();
}

/** This thread calls computeImage to compute new image data and does the callbacks to send it to higher layers.
  * It implements the logic for single, multiple or continuous acquisition. */
void ADCSimDetector::simTask()
{
    int status = asynSuccess;
    int acquire=0;
    NDArray *pImage;
    epicsTimeStamp startTime;
    const char *functionName = "simTask";

    this->lock();
    /* Loop forever */
    while (1) {
        /* Has acquisition been stopped? */
        status = epicsEventTryWait(this->stopEventId_);
        if (status == epicsEventWaitOK) {
            acquire = 0;
        }
       
        /* If we are not acquiring then wait for a semaphore that is given when acquisition is started */
        if (!acquire) {
          /* Release the lock while we wait for an event that says acquire has started, then lock again */
            asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                "%s:%s: waiting for acquire to start\n", driverName, functionName);
            this->unlock();
            status = epicsEventWait(this->startEventId_);
            this->lock();
            acquire = 1;
            currentPoint_ = 0;
        }

        /* Get the current time */
        epicsTimeGetCurrent(&startTime);

        /* Update the data */
        computeArrays();

        pImage = this->pArrays[0];

        /* Put the frame number and time stamp into the buffer */
        pImage->uniqueId = arrayCounter_++;
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
        callParamCallbacks();
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
    int acquiring;
    int addr;
    asynStatus status = asynSuccess;

    callParamCallbacks();
    getAddress(pasynUser, &addr);
 
    getIntegerParam(P_Acquire, &acquiring);
    /* Set the parameter and readback in the parameter library.  This may be overwritten when we read back the
     * status at the end, but that's OK */
    status = setIntegerParam(function, value, addr);

    if (function == P_Acquire) {
        if (value && !acquiring) {
            /* Send an event to wake up the simulation task.
             * It won't actually start generating new images until we release the lock below */
            epicsEventSignal(this->startEventId_); 
        }
        if (!value && acquiring) {
            /* This was a command to stop acquisition */
            /* Send the stop event */
            epicsEventSignal(this->stopEventId_); 
        }
    } else {
        /* If this parameter belongs to a base class call its method */
        if (function < FIRST_SIM_DETECTOR_PARAM) status = asynNDArrayDriver::writeInt32(pasynUser, value);
    }

    /* Do callbacks so higher layers see any changes */
    callParamCallbacks();

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

/** Constructor for ADCSimDetector; most parameters are simply passed to ADDriver::ADDriver.
  * After calling the base class constructor this method creates a thread to compute the simulated detector data,
  * and sets reasonable default values for parameters defined in this class, asynNDArrayDriver and ADDriver.
  * \param[in] portName The name of the asyn port driver to be created.
  * \param[in] numPointPoinst The initial number of time points.
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
               0, 1, /* ASYN_CANBLOCK=0, ASYN_MULTIDEVICE=0, autoConnect=1 */
               priority, stackSize)

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
    createParam(SimNumTimePointsString,   asynParamInt32, &P_NumTimePoints);
    createParam(SimAmplitudeString,     asynParamFloat64, &P_Amplitude);
    createParam(SimFrequencyString,     asynParamFloat64, &P_Frequency);
    createParam(SimPhaseString,         asynParamFloat64, &P_Phase);
    createParam(SimNoiseString,         asynParamFloat64, &P_Noise);

    status |= setIntegerParam(NDArraySizeX,    numTimePoints);
    status |= setIntegerParam(P_NumTimePoints, numTimePoints);
    status |= setIntegerParam(NDArraySize, 0);
    status |= setIntegerParam(NDDataType, dataType);
    status |= setDoubleParam(P_Amplitude, 1.0);
    status |= setDoubleParam(P_Frequency, 1.0);
    status |= setDoubleParam(P_Phase, 0.0);
    status |= setDoubleParam(P_Noise, 0.0);

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
