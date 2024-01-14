#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include "TError.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TVectorD.h"

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

enum Directions
{
    x = 0,
    y,
    z,
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

struct IBDvalues
{
    std::array<std::array<double, DirectionSize>, DatasetSize> effectiveIBD;
    std::array<std::array<double, DirectionSize>, DatasetSize> mean;
    std::array<std::array<double, DirectionSize>, DatasetSize> sigma;
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
            name = "x";
            break;
        case 1:
            name = "y";
            break;
        case 2:
            name = "z";
            break;
        default:
            name = "x";
    }

    return name;
}
