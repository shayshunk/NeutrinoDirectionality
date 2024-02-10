#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"

#define pi 3.14159265358979323846

// Invariables
int const bins = 301;  // Number of bins for our histograms
int const totalDataLines = 4000;
int const totalSimLines = 500;
float const histogramMax = 150.5;  // Maximum value of bin for histograms
float const segmentWidth = 145.7;  // Distance between segment centers in mm
float const atmosphericScaling = 1.00025443769309;  // Atmosphering scaling coefficient for background subtraction
char const* dataPath = "/home/shay/Documents/PROSPECTData/IBD_Data/SEER_DS_period_%s/Period_%s_files.txt";
char const* dataFileName = "/home/shay/Documents/PROSPECTData/IBD_Data/SEER_DS_period_%s/%s/AD1_IBD_2022_DS_SEER.root";
char const* simPath = "/home/shay/Documents/PROSPECTData/MC_Data/DS_SEER_MC/period_%s/period_%s.txt";
char const* simFileName = "/home/shay/Documents/PROSPECTData/MC_Data/DS_SEER_MC/period_%s/%s/AD1_IBD_2022_DS_SEER.root";

// Variables
int lineCounter = 0;
double livetimeOff = 0, livetimeOn = 0;

// Print flags
bool DETECTOR_VERBOSITY = 1;
bool LIVETIME_VERBOSITY = 1;
bool IBDCOUNT_VERBOSITY = 1;
bool MEAN_VERBOSITY = 1;
bool SYSTEMATIC_MEAN_VERBOSITY = 1;
bool ANGLES_STATISTICS = 1;
bool COVARIANCE_VERBOSITY = 1;

// Utilities to refer to parameters easily
enum Directions
{
    X = 0,
    Y,
    Z,
    DirectionSize
};

enum Datasets
{
    Data = 0,
    DataUnbiased,
    Sim,
    SimUnbiased,
    DatasetSize
};

enum Signals
{
    CorrelatedReactorOn = 0,
    CorrelatedReactorOff,
    AccidentalReactorOn,
    AccidentalReactorOff,
    TotalDifference,
    SignalSize
};

struct TreeValues
{
    double Esmear;
    double nCaptTime;
    double xRx;
    int promptSegment;
    int delayedSegment;
    int dataSet;
    int direction;
    int period;
    double promptPosition;
    double delayedPosition;
    bool reactorOn;
};

struct IBDValues
{
    std::array<std::array<float, DirectionSize>, DatasetSize> effectiveIBD;
    std::array<std::array<float, DirectionSize>, DatasetSize> totalIBD;
    std::array<std::array<float, DirectionSize>, DatasetSize> totalIBDError;
    std::array<std::array<float, DirectionSize>, DatasetSize> mean;
    std::array<std::array<float, DirectionSize>, DatasetSize> sigma;
    std::array<std::array<float, DirectionSize>, DatasetSize> sigmaSystematics;
};

struct AngleValues
{
    std::array<double, DatasetSize> phi;
    std::array<double, DatasetSize> phiError;
    std::array<double, DatasetSize> phiErrorSystematics;
    std::array<double, DatasetSize> theta;
    std::array<double, DatasetSize> thetaError;
    std::array<double, DatasetSize> thetaErrorSystematics;
    double phiTrue, phiTrueError, thetaTrue, thetaTrueError;
};

struct CovarianceValues
{
    std::array<double, DatasetSize> phiError;
    std::array<double, DatasetSize> phiErrorSystematics;
    std::array<double, DatasetSize> thetaError;
    std::array<double, DatasetSize> thetaErrorSystematics;
    std::array<double, DatasetSize> tilt;
    std::array<double, DatasetSize> tiltSystematics;
};

std::string DatasetToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0:
            name = "Data";
            break;
        case 1:
            name = "Data Unbiased";
            break;
        case 2:
            name = "Simulation";
            break;
        case 3:
            name = "Simulation Unbiased";
            break;
        case 4:
            name = "True";
            break;
        default:
            name = "Data";
    }

    return name;
}

std::string SignalToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0:
            name = "Correlated Reactor On";
            break;
        case 1:
            name = "Correlated Reactor Off";
            break;
        case 2:
            name = "Accidental Reactor On";
            break;
        case 3:
            name = "Accidental Reactor Off";
            break;
        case 4:
            name = "Total Difference";
            break;
        default:
            name = "Total Difference";
    }

    return name;
}

std::string AxisToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0:
            name = "X";
            break;
        case 1:
            name = "Y";
            break;
        case 2:
            name = "Z";
            break;
        default:
            name = "X";
    }

    return name;
}

class Directionality
{
  public:
    Directionality();
    void FillHistogramUnbiased(int signalSet);
    void FillHistogram();
    void SetUpHistograms(int Data, int period = 0);
    void CalculateUnbiasing();
    void SubtractBackgrounds();
    void AddSystematics();
    void CalculateAngles();
    void CalculateCovariances();
    void OffsetTheta();
    void PrintAngles();
    void FillOutputFile();

  private:
    // Histogram to count IBDs
    array<array<array<TH1F, DirectionSize>, SignalSize>, DatasetSize> histogram;

    // Values grabbed from ROOT tree
    double Esmear;
    double nCaptTime;
    double xRx;
    double promptPosition, delayedPosition;
    int promptSegment, delayedSegment;
    int dataSet;
    int direction;
    int period;
    bool reactorOn;

    // Storing final counts
    std::array<std::array<float, DirectionSize>, DatasetSize> effectiveIBD;
    std::array<std::array<float, DirectionSize>, DatasetSize> totalIBD;
    std::array<std::array<float, DirectionSize>, DatasetSize> totalIBDError;
    std::array<std::array<float, DirectionSize>, DatasetSize> mean;
    std::array<std::array<float, DirectionSize>, DatasetSize> sigma;
    std::array<std::array<float, DirectionSize>, DatasetSize> sigmaSystematics;

    // Storing final angles
    std::array<double, DatasetSize> phi;
    std::array<double, DatasetSize> phiError;
    std::array<double, DatasetSize> phiErrorSystematics;
    std::array<double, DatasetSize> theta;
    std::array<double, DatasetSize> thetaError;
    std::array<double, DatasetSize> thetaErrorSystematics;
    double phiTrue, phiTrueError, thetaTrue, thetaTrueError;

    // Storing final ellipse values
    std::array<double, DatasetSize> phiEllipseError;
    std::array<double, DatasetSize> phiEllipseErrorSystematics;
    std::array<double, DatasetSize> thetaEllipseError;
    std::array<double, DatasetSize> thetaEllipseErrorSystematics;
    std::array<double, DatasetSize> tilt;
    std::array<double, DatasetSize> tiltSystematics;
};
