#include "TH1D.h"
#include "TError.h"
#include "TKey.h"
#include "TFile.h"

enum Datasets
{
    Data = 0, DataUnbiased, Sim, SimUnbiased, True
};

enum Signals
{
    CorrelatedOn = 0, CorrelatedOff, AcceleratedOn, AcceleratedOff
};

