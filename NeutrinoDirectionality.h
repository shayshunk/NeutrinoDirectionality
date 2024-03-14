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

// Print flags
bool DETECTOR_VERBOSITY = 0;
bool LIVETIME_VERBOSITY = 0;
bool IBDCOUNT_VERBOSITY = 0;
bool NCOUNT_VERBOSITY = 0;
bool MEAN_VERBOSITY = 0;
bool SYSTEMATIC_MEAN_VERBOSITY = 0;
bool ANGLES_STATISTICS = 0;
bool COVARIANCE_VERBOSITY = 0;

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
    // Variables
    double livetimeOff = 0, livetimeOn = 0;

    // Main functions
    Directionality();
    void ReadFileList(int dataset, int periodNo);
    void SetUpHistograms(int dataset, int periodNo = 0);
    void FillHistogramUnbiased(int signalSet);
    void FillHistogram();
    void CalculateUnbiasing();
    void SubtractBackgrounds();
    void AddSystematics();
    void CalculateAngles();
    void CalculateCovariances();
    void OffsetTheta();
    void PrintAngles();
    void FillOutputFile();

    // Inline functions
    inline void ResetLineNumber() { lineNumber = 0; }
    inline void ResetLineCounter() { lineCounter = 0; }
    inline void ResetIndex() { index = 0; }

  private:
    // Histogram to count IBDs
    std::array<std::array<std::array<TH1F, DirectionSize>, SignalSize>, DatasetSize> histogram;

    // File list
    std::array<std::string, 4045> files;

    // Values grabbed from ROOT tree
    double Esmear;
    double nCaptTime;
    double xRx;
    double promptPosition, delayedPosition;
    int promptSegment, delayedSegment;
    int dataSet;
    int direction;
    int period;
    int lineNumber = 0, lineCounter = 0;
    std::size_t index = 0;
    bool reactorOn;

    // Invariables
    int const bins = 301;  // Number of bins for our histograms
    int const totalDataLines = 4040;
    int const totalSimLines = 500;
    float const histogramMax = 150.5;  // Maximum value of bin for histograms
    float const segmentWidth = 145.7;  // Distance between segment centers in mm
    float const atmosphericScaling = 1.00025443769309;  // Atmosphering scaling coefficient for background subtraction
    char const* dataPath = "/home/shay/Documents/PROSPECTData/IBD_Data/SEER_DS_period_%s/Period_%s_files.txt";
    char const* dataFileName = "/home/shay/Documents/PROSPECTData/IBD_Data/SEER_DS_period_%s/%s/AD1_IBD_2022_DS_SEER.root";
    char const* simPath = "/home/shay/Documents/PROSPECTData/MC_Data/DS_SEER_MC/period_%s/period_%s.txt";
    char const* simFileName = "/home/shay/Documents/PROSPECTData/MC_Data/DS_SEER_MC/period_%s/%s/AD1_IBD_2022_DS_SEER.root";

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
