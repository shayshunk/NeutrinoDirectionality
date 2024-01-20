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
#include "TVector2.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TF1.h"

#define pi 3.14159265358979323846

// Invariables
int const xBins = 12;
int const yBins = 9;
int const zBins = 301;  // Number of bins for our histograms
int const totalDataLines = 4000;
int const totalSimLines = 500;
float const zHistogramMax = 150.5;  // Maximum value of bin for histogram
float const segmentWidth = 145.7;  // Distance between segment centers in mm
float const atmosphericScaling = 1.00025443769309;  // Atmosphering scaling coefficient for background subtraction
char const* dataPath = "/home/shay/Documents/PROSPECTData/MC_Data/PerfectSimulation/SimFileGoodRunsList.txt";
char const* dataFileName = "/home/shay/Documents/PROSPECTData/MC_Data/PerfectSimulation/%s/Jun_Perfect_IBD_2020.root";

// Variables
int lineCounter = 0;
double livetimeOff = 0, livetimeOn = 0;

enum Directions
{
    X = 0,
    Y,
    Z,
    DirectionSize
};

enum Signals
{
    CorrelatedReactorOn = 0,
    AccidentalReactorOn,
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
    int direction;
    int period;
    double promptPosition;
    double delayedPosition;
    bool reactorOn;
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

