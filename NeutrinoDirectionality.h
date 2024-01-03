#include "TError.h"
#include "TKey.h"
#include "TFile.h"
#include "TH1D.h"

#include <iostream>
#include <string>
#include <array>

enum Datasets
{
    Data = 0, DataUnbiased, Sim, SimUnbiased, True
};

enum Signals
{
    CorrelatedReactorOn = 0, CorrelatedReactorOff, AccidentalReactorOn, AccidentalReactorOff, TotalDifference
};

const int bins = 301; // Number of bins for our histograms
const float histogramMax = 150.5; // Maximum value of bin for histograms

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