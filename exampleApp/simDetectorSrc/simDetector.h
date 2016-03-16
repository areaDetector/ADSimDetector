#include <epicsEvent.h>
#include "ADDriver.h"

/** Simulation detector driver; demonstrates most of the features that areaDetector drivers can support. */
class epicsShareClass simDetector : public ADDriver {
public:
    simDetector(const char *portName, int maxSizeX, int maxSizeY, NDDataType_t dataType,
                int maxBuffers, size_t maxMemory,
                int priority, int stackSize);

    /* These are the methods that we override from ADDriver */
    virtual asynStatus writeInt32(asynUser *pasynUser, epicsInt32 value);
    virtual asynStatus writeFloat64(asynUser *pasynUser, epicsFloat64 value);
    virtual void setShutter(int open);
    virtual void report(FILE *fp, int details);
    void simTask(); /**< Should be private, but gets called from C, so must be public */

protected:
    int SimGainX;
    #define FIRST_SIM_DETECTOR_PARAM SimGainX
    int SimGainY;
    int SimGainRed;
    int SimGainGreen;
    int SimGainBlue;
    int SimNoise;
    int SimResetImage;
    int SimMode;
    int SimPeakStartX;
    int SimPeakStartY;
    int SimPeakWidthX;
    int SimPeakWidthY;
    int SimPeakNumX;
    int SimPeakNumY;
    int SimPeakStepX;
    int SimPeakStepY;
    int SimPeakHeightVariation;
    int SimXSine1Amplitude;
    int SimXSine1Frequency;
    int SimXSine1Phase;
    int SimXSine1Offset;
    int SimXSine1Noise;
    int SimXSine2Amplitude;
    int SimXSine2Frequency;
    int SimXSine2Phase;
    int SimXSine2Offset;
    int SimXSine2Noise;
    int SimYSine1Amplitude;
    int SimYSine1Frequency;
    int SimYSine1Phase;
    int SimYSine1Offset;
    int SimYSine1Noise;
    int SimYSine2Amplitude;
    int SimYSine2Frequency;
    int SimYSine2Phase;
    int SimYSine2Offset;
    int SimYSine2Noise;

    #define LAST_SIM_DETECTOR_PARAM SimYSine2Noise

private:
    /* These are the methods that are new to this class */
    template <typename epicsType> int computeArray(int sizeX, int sizeY);
    template <typename epicsType> int computeLinearRampArray(int sizeX, int sizeY);
    template <typename epicsType> int computePeaksArray(int sizeX, int sizeY);
    template <typename epicsType> int computeSineArray(int sizeX, int sizeY);
    int computeImage();

    /* Our data */
    epicsEventId startEventId_;
    epicsEventId stopEventId_;
    NDArray *pRaw_;
    double *xSine1_;
    double *xSine2_;
    double *ySine1_;
    double *ySine2_;
    double xSineCounter_;
    double ySineCounter_;
};

typedef enum {
    SimModeLinearRamp,
    SimModePeaks,
    SimModeSine
}SimModes_t;

#define SimGainXString                "SIM_GAIN_X"
#define SimGainYString                "SIM_GAIN_Y"
#define SimGainRedString              "SIM_GAIN_RED"
#define SimGainGreenString            "SIM_GAIN_GREEN"
#define SimGainBlueString             "SIM_GAIN_BLUE"
#define SimNoiseString                "SIM_NOISE"
#define SimResetImageString           "RESET_IMAGE"
#define SimModeString                 "SIM_MODE"
#define SimPeakStartXString           "SIM_PEAK_START_X"
#define SimPeakStartYString           "SIM_PEAK_START_Y"
#define SimPeakWidthXString           "SIM_PEAK_WIDTH_X"
#define SimPeakWidthYString           "SIM_PEAK_WIDTH_Y"
#define SimPeakNumXString             "SIM_PEAK_NUM_X"
#define SimPeakNumYString             "SIM_PEAK_NUM_Y"
#define SimPeakStepXString            "SIM_PEAK_STEP_X"
#define SimPeakStepYString            "SIM_PEAK_STEP_Y"
#define SimPeakHeightVariationString  "SIM_PEAK_HEIGHT_VARIATION"
#define SimXSine1AmplitudeString      "SIM_XSINE1_AMPLITUDE"
#define SimXSine1FrequencyString      "SIM_XSINE1_FREQUENCY"
#define SimXSine1PhaseString          "SIM_XSINE1_PHASE"
#define SimXSine1OffsetString         "SIM_XSINE1_OFFSET"
#define SimXSine1NoiseString          "SIM_XSINE1_NOISE"
#define SimXSine2AmplitudeString      "SIM_XSINE2_AMPLITUDE"
#define SimXSine2FrequencyString      "SIM_XSINE2_FREQUENCY"
#define SimXSine2PhaseString          "SIM_XSINE2_PHASE"
#define SimXSine2OffsetString         "SIM_XSINE2_OFFSET"
#define SimXSine2NoiseString          "SIM_XSINE2_NOISE"
#define SimYSine1AmplitudeString      "SIM_YSINE1_AMPLITUDE"
#define SimYSine1FrequencyString      "SIM_YSINE1_FREQUENCY"
#define SimYSine1PhaseString          "SIM_YSINE1_PHASE"
#define SimYSine1OffsetString         "SIM_YSINE1_OFFSET"
#define SimYSine1NoiseString          "SIM_YSINE1_NOISE"
#define SimYSine2AmplitudeString      "SIM_YSINE2_AMPLITUDE"
#define SimYSine2FrequencyString      "SIM_YSINE2_FREQUENCY"
#define SimYSine2PhaseString          "SIM_YSINE2_PHASE"
#define SimYSine2OffsetString         "SIM_YSINE2_OFFSET"
#define SimYSine2NoiseString          "SIM_YSINE2_NOISE"

#define NUM_SIM_DETECTOR_PARAMS ((int)(&LAST_SIM_DETECTOR_PARAM - &FIRST_SIM_DETECTOR_PARAM + 1))

