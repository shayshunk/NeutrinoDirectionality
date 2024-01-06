#include "TError.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TLeaf.h"

#include <iostream>
#include <string>
#include <array>
#include <fstream>

const int bins = 301; // Number of bins for our histograms
const float histogramMax = 150.5; // Maximum value of bin for histograms
const float segmentWidth = 145.7; // Distance between segment centers in mm
const char* dataPath = "/project/prospect/tmp/Analyzed/Analyzed_2022_IBD_DS_SEER_Data/DS_SEER_Data_NewCal/SEER_DS_period_%s/Period_%s_files.txt";
const char* fileName = "/project/prospect/tmp/Analyzed/Analyzed_2022_IBD_DS_SEER_Data/DS_SEER_Data_NewCal/SEER_DS_period_%s/%s/AD1_IBD_2022_DS_SEER.root";

enum Directions
{
    x = 0, y, z
};

enum Datasets
{
    Data = 0, DataUnbiased, Sim, SimUnbiased, True
};

enum Signals
{
    CorrelatedReactorOn = 0, CorrelatedReactorOff, AccidentalReactorOn, AccidentalReactorOff, TotalDifference
};

struct TreeValues
{
    double Esmear;
    double nCaptTime;
    int promptSegment;
    int neutronSegment;
    int dataSet;
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
        case 0: { name = "Data"; } break;
        case 1: { name = "DataUnbiased"; } break;
        case 2: { name = "Sim"; } break;
        case 3: { name = "SimUnbiased"; } break;
        case 4: { name = "True"; } break;
    }

    return name;
}

std::string SignalToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0: { name = "CorrelatedReactorOn"; } break;
        case 1: { name = "CorrelatedReactorOff"; } break;
        case 2: { name = "AccidentalReactorOn"; } break;
        case 3: { name = "AccidentalReactorOff"; } break;
        case 4: { name = "TotalDifference"; } break;
    }

    return name;
}

std::string AxisToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0: { name = "x"; } break;
        case 1: { name = "y"; } break;
        case 2: { name = "z"; } break;
    }

    return name;
}

