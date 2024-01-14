
#include "NeutrinoDirectionality.h"

#include "DetectorConfig.h"

using std::cout, std::string, std::ifstream, std::array, std::getline;

void FillDetectorConfig()
{
    // Filling values based on different periods
    int noSegments = 154;
    int noPeriods = 5;

    for (int i = 0; i < noPeriods; i++)
    {
        int excludeSize = excludeList[i].size();
        int counter = 0, tmp = 0;
        vector<int> period;

        for (int j = 0; j < noSegments; j++)
        {
            // Filling live segments by checking exclude list
            tmp = excludeList[i][counter];

            if (j == tmp)
            {
                period.push_back(0);
                if (counter + 1 < excludeSize)
                {
                    counter += 1;
                }
            }
            else
            {
                period.push_back(1);
            }
        }

        detectorConfig.push_back(period);
    }

    cout << "Below is the detector configuration.\n";
    cout << "--------------------------------------------\n";

    for (int i = 0; i < detectorConfig.size(); i++)
    {
        if (i != 5)
            cout << "Detector configuration for period: " << i + 1 << '\n';
        else
            cout << "Detector configuration for simulation: \n";

        for (int j = 0; j < detectorConfig[i].size(); j++)
        {
            if (detectorConfig[i][j])
            {
                cout << "\u25A0 ";
            }
            else
            {
                cout << "\u25A1 ";
            }
            //cout << detectorConfig[i][j] << " ";
            if ((j + 1) % 14 == 0)
                cout << '\n';
        }
        cout << '\n';
    }
    cout << "--------------------------------------------\n";
}

bool checkNeighbor(int periodNo, int segment, char direction)
{
    // Used for dead segment calculations

    bool neighbor = false;

    periodNo = periodNo - 1;

    switch (direction)
    {
        case 'r':
            neighbor = detectorConfig[periodNo][segment + 1];
            break;
        case 'l':
            neighbor = detectorConfig[periodNo][segment - 1];
            break;
        case 'u':
            neighbor = detectorConfig[periodNo][segment + 14];
            break;
        case 'd':
            neighbor = detectorConfig[periodNo][segment - 14];
            break;
        default:
            cout << "That direction doesn't exist!\n";
            return false;
    }

    return neighbor;
}

void FillHistogramUnbiased(array<array<array<std::shared_ptr<TH1D>, DirectionSize>, SignalSize>, DatasetSize>& histogram,
                           TreeValues& currentEntry,
                           int signalSet)
{
    bool posDirection = false, negDirection = false, success = false;
    double weight = 0;

    // Check for live neighbors in different directions based on which axis
    // we're filling
    if (currentEntry.direction == x)
    {
        posDirection = checkNeighbor(currentEntry.period, currentEntry.promptSegment, 'r');
        negDirection = checkNeighbor(currentEntry.period, currentEntry.promptSegment, 'l');
    }
    else if (currentEntry.direction == y)
    {
        posDirection = checkNeighbor(currentEntry.period, currentEntry.promptSegment, 'u');
        negDirection = checkNeighbor(currentEntry.period, currentEntry.promptSegment, 'd');
    }

    // Need to weight accidental datasets by deadtime correction factor
    weight = (signalSet == AccidentalReactorOff || signalSet == AccidentalReactorOn) ? currentEntry.xRx : 1;

    // Dataset + 1 returns the unbiased version of that dataset
    if (posDirection && !negDirection)
        histogram[currentEntry.dataSet + 1][signalSet][currentEntry.direction]->Fill(segmentWidth, weight);
    else if (!posDirection && negDirection)
        histogram[currentEntry.dataSet + 1][signalSet][currentEntry.direction]->Fill(-segmentWidth, weight);
    else if (posDirection && negDirection)
        histogram[currentEntry.dataSet + 1][signalSet][currentEntry.direction]->Fill(0.0, weight);
}

void FillHistogram(array<array<array<std::shared_ptr<TH1D>, DirectionSize>, SignalSize>, DatasetSize>& histogram,
                   TreeValues& currentEntry)
{
    // Applying energy cut
    if (currentEntry.Esmear < 0.8 || currentEntry.Esmear > 7.4)
    {
        return;
    }

    if (currentEntry.nCaptTime > pow(10, 3) && currentEntry.nCaptTime < 120 * pow(10, 3))  // Correlated Dataset
    {
        // Calculate neutron displacement
        double displacement = currentEntry.delayedPosition - currentEntry.promptPosition;

        // Figure out whether the reactor is on and assign signal index
        int signalSet = currentEntry.reactorOn ? CorrelatedReactorOn : CorrelatedReactorOff;

        // Fill regular dataset with displacement
        histogram[currentEntry.dataSet][signalSet][currentEntry.direction]->Fill(displacement);

        // Fill dead segment correction dataset
        if (currentEntry.promptSegment == currentEntry.delayedSegment)
        {
            FillHistogramUnbiased(histogram, currentEntry, signalSet);
        }
    }
    else if (currentEntry.nCaptTime > pow(10, 6))  // Accidental Dataset
    {
        // Calculate neutron displacement
        double displacement = currentEntry.delayedPosition - currentEntry.promptPosition;

        // Figure out whether the reactor is on and assign signal index
        int signalSet = currentEntry.reactorOn ? AccidentalReactorOn : AccidentalReactorOff;

        // Fill regular dataset with displacement
        histogram[currentEntry.dataSet][signalSet][currentEntry.direction]->Fill(displacement, currentEntry.xRx);

        // Fill dead segment correction dataset
        if (currentEntry.promptSegment == currentEntry.delayedSegment)
        {
            FillHistogramUnbiased(histogram, currentEntry, signalSet);
        }
    }
}

void SetUpHistograms(array<array<array<std::shared_ptr<TH1D>, DirectionSize>, SignalSize>, DatasetSize>& histogram,
                     int dataSet,
                     int period = 0)
{
    // Declaring some variables for use later
    int totalLines = 0;
    bool reactorOn = true;

    // Figuring out dataset
    char const* path;
    char const* fileName;
    if (dataSet == Data)
    {
        path = dataPath;
        fileName = dataFileName;
    }
    else if (dataSet == Sim)
    {
        path = simPath;
        fileName = simFileName;
    }

    // Combining names into file list name
    string fileList = Form(path, std::to_string(period).c_str(), std::to_string(period).c_str());

    // Opening and checking file list
    ifstream file;
    file.open(fileList, ifstream::in);

    if (!(file.is_open() && file.good()))
    {
        cout << "File list not found! Exiting.\n";
        cout << "Trying to find: " << fileList << '\n';
        return;
    }

    while (file.good() && !file.eof())
    {
        lineCounter++;

        if (dataSet == Data || dataSet == DataUnbiased)
        {
            if (lineCounter % 200 == 0)
            {
                cout << "Reading file: " << lineCounter << "/" << totalDataLines << '\r';
                cout.flush();
            }
        }
        else if (dataSet == Sim || dataSet == SimUnbiased)
        {
            if (lineCounter % 50 == 0)
            {
                cout << "Reading file: " << lineCounter << "/" << totalSimLines << '\r';
                cout.flush();
            }
        }

        // Reading file list
        string line;
        getline(file, line);

        // Combining names into root file name
        TString rootFilename = Form(fileName, std::to_string(period).c_str(), line.data());

        if (rootFilename.Contains(" 0"))
        {
            rootFilename.ReplaceAll(" 0", "");
            reactorOn = false;
        }
        else if (rootFilename.Contains(" 1"))
        {
            rootFilename.ReplaceAll(" 1", "");
            reactorOn = true;
        }

        // Open the root file
        auto rootFile = std::make_unique<TFile>(rootFilename);

        // Declare deadtime correction coefficient (veto deadtime correction)
        double xRx;

        // Going into empty scope to let the pointers die out for safety
        {
            TVectorD* runtime = (TVectorD*)rootFile->Get("runtime");
            TVectorD* promptVeto = (TVectorD*)rootFile->Get("accumulated/P2kIBDPlugin.tveto_prompt");  // prompt veto deadtime
            TVectorD* delayedVeto = (TVectorD*)rootFile->Get("accumulated/P2kIBDPlugin.tveto_delayed");  // delayed veto deadtime
            xRx = runtime->Max() / (runtime->Max() - promptVeto->Max()) * runtime->Max() / (runtime->Max() - delayedVeto->Max());

            if (reactorOn && dataSet == Data)
            {
                livetimeOn += runtime->Max() / xRx;
            }
            else if (dataSet == Data)
            {
                livetimeOff += runtime->Max() / xRx;
            }
        }

        // Grab rootTree and cast to unique pointer
        TTree* rootTree = (TTree*)rootFile->Get("P2kIBDPlugin/Tibd");

        long nEntries = rootTree->GetEntries();

        for (long i = 0; i < nEntries; i++)
        {
            rootTree->GetEntry(i);

            for (int direction = x; direction < DirectionSize; direction++)
            {
                // Intializing struct of relevant values
                TreeValues currentEntry;

                // Grabbing relevant values from the rootTree entry
                currentEntry.promptPosition = rootTree->GetLeaf("xyz")->GetValue(direction);
                currentEntry.delayedPosition = rootTree->GetLeaf("n_xyz")->GetValue(direction);
                currentEntry.promptSegment = rootTree->GetLeaf("maxseg")->GetValue(0);
                currentEntry.delayedSegment = rootTree->GetLeaf("n_seg")->GetValue(0);

                // We throw out events where the neutron moves in a direction we're not checking
                if (currentEntry.promptSegment != currentEntry.delayedSegment
                    && currentEntry.promptPosition == currentEntry.delayedPosition)
                    continue;

                // Copying some loop values into current entry
                currentEntry.Esmear = rootTree->GetLeaf("Esmear")->GetValue(0);
                currentEntry.nCaptTime = rootTree->GetLeaf("ncapt_dt")->GetValue(0);
                currentEntry.xRx = xRx;
                currentEntry.dataSet = dataSet;
                currentEntry.period = period;
                currentEntry.direction = direction;
                currentEntry.reactorOn = reactorOn;

                FillHistogram(histogram, currentEntry);
            }
        }
        // Returns the next character in the input sequence, without extracting it: The character is left as the next character
        // to be extracted from the stream
        file.peek();
        rootFile->Close();
    }
}

void CalculateAngles(array<array<array<std::shared_ptr<TH1D>, DirectionSize>, SignalSize>, DatasetSize>& histogram, IBDvalues& neutrinoCounts)
{
}

void CalculateUnbiasing(array<array<array<std::shared_ptr<TH1D>, DirectionSize>, SignalSize>, DatasetSize>& histogram, IBDvalues& neutrinoCounts)
{
    // Defining variables used in calculation. Check the error propagation technote for details on the method
    double rPlus = 0, rMinus = 0;
    double p = 0, pError = 0;
    double nPlus = 0, nPlusPlus = 0, nMinus = 0, nMinusMinus = 0, nPlusMinus = 0;
    double nPlusError = 0, nPlusPlusError = 0, nMinusError = 0, nMinusMinusError = 0, nPlusMinusError = 0;

    for (int dataset = DataUnbiased; dataset < DatasetSize; dataset+=2)
    {   
        for (int direction = x; direction < z; direction++)
        {
            // Dataset - 1 returns the biased dataset counts
            nPlus = histogram[dataset - 1][TotalDifference][direction]->GetBinContent(296);
            nPlusPlus = histogram[dataset][TotalDifference][direction]->GetBinContent(297);
            nMinus = histogram[dataset - 1][TotalDifference][direction]->GetBinContent(6);
            nMinusMinus = histogram[dataset][TotalDifference][direction]->GetBinContent(5);
            nPlusMinus = histogram[dataset][TotalDifference][direction]->GetBinContent(151);

            nPlusError = histogram[dataset - 1][TotalDifference][direction]->GetBinError(296);
            nPlusPlusError = histogram[dataset][TotalDifference][direction]->GetBinError(297);
            nMinusError = histogram[dataset - 1][TotalDifference][direction]->GetBinError(6);
            nMinusMinusError = histogram[dataset][TotalDifference][direction]->GetBinError(5);
            nPlusMinusError = histogram[dataset][TotalDifference][direction]->GetBinError(151);

            rPlus = nPlus / (nPlusPlus + nPlusMinus);
            rMinus = nMinus / (nMinusMinus + nPlusMinus);

            p = segmentWidth * (rPlus - rMinus) / (rPlus + rMinus + 1);

            pError = segmentWidth * pow(1 / ((nMinus * (nPlusMinus + nPlusPlus) + (nMinusMinus + nPlusMinus) * (nPlus + nPlusMinus + nPlusPlus))), 2) * sqrt(pow((nMinusMinus + nPlusMinus) * (nPlusMinus + nPlusPlus), 2) * (pow(nPlusError * (2 * nMinus + nMinusMinus + nPlusMinus), 2)
                                                                                                                                        + pow(nMinusError * (2 * nPlus + nPlusPlus + nPlusMinus), 2))
                                                                                                    + pow((nPlus * (nPlusMinus + nMinusMinus) * (2 * nMinus + nMinusMinus + nPlusMinus) * nPlusPlusError), 2)
                                                                                                    + pow((nPlusMinusError * (nPlus * pow((nMinusMinus + nPlusMinus), 2) + nMinus * (2 * nMinusMinus * nPlus - 2 * nPlus * nPlusPlus - pow((nPlusMinus + nPlusPlus), 2)))), 2)
                                                                                                    + pow((nMinus * (nPlusMinus + nPlusPlus) * (2 * nPlus + nPlusMinus + nPlusPlus) * nMinusMinusError), 2));
            
            neutrinoCounts.mean[dataset][direction] = p;
            neutrinoCounts.sigma[dataset][direction] = pError;
        }
        neutrinoCounts.effectiveIBD[dataset][z] = neutrinoCounts.effectiveIBD[dataset - 1][z]; 
        neutrinoCounts.mean[dataset][z] = neutrinoCounts.mean[dataset - 1][z];
        neutrinoCounts.sigma[dataset][z] = neutrinoCounts.sigma[dataset - 1][z];
    }

    // Printing out values
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        cout << "Mean and sigma values for: " << DatasetToString(dataset) << '\n';
        for (int direction = x; direction < DirectionSize; direction++)
        {
            cout << "p" << AxisToString(direction) << ": " << neutrinoCounts.mean[dataset][direction] << " ± " << neutrinoCounts.sigma[dataset][direction] << '\n';
        }
        cout << "--------------------------------------------\n";
    }
}

void SubtractBackgrounds(array<array<array<std::shared_ptr<TH1D>, DirectionSize>, SignalSize>, DatasetSize>& histogram)
{
    /* IBD events = (Correlated - Accidental/100)_{reactor on} + (-livetimeOn/livetimeOff*Correlated +
    livetimeOn/livetimeOff*Accidental/100)_{reactor off} */

    // Defining variables for IBD background subtraction
    double totalIBDs = 0, totalIBDErr = 0, effIBDs = 0;

    IBDvalues neutrinoCounts;

    cout << "--------------------------------------------\n";
    cout << "Subtracting backgrounds.\n";

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        cout << "Total IBD events for: " << DatasetToString(dataset) << '\n';
        for (int direction = x; direction < DirectionSize; direction++)
        {
            // Have to static cast raw pointer to shared pointer to keep up safety measures
            histogram[dataset][TotalDifference][direction] = std::shared_ptr<TH1D>(
                static_cast<TH1D*>(histogram[dataset][CorrelatedReactorOn][direction]->Clone("Displacements")));
            histogram[dataset][TotalDifference][direction]->Add(histogram[dataset][AccidentalReactorOn][direction].get(),
                                                                -1. / 100.);

            if (dataset == Data || dataset == DataUnbiased)
            {
                histogram[dataset][TotalDifference][direction]->Add(histogram[dataset][CorrelatedReactorOff][direction].get(),
                                                                    -livetimeOn * atmosphericScaling / livetimeOff);
                histogram[dataset][TotalDifference][direction]->Add(histogram[dataset][AccidentalReactorOff][direction].get(),
                                                                    livetimeOn * atmosphericScaling / (100 * livetimeOff));
            }

            totalIBDs = histogram[dataset][TotalDifference][direction]->IntegralAndError(
                0, histogram[dataset][TotalDifference][direction]->GetNbinsX() + 1, totalIBDErr);
            effIBDs = pow(totalIBDs, 2) / pow(totalIBDErr, 2);  // Effective IBD counts. Done by Poisson Distribution
                                                                // N^2/(sqrt(N)^2) = N; Eff. counts = counts^2/counts_err^2
            cout << AxisToString(direction) << ": " << totalIBDs << " ± " << totalIBDErr << ". Effective IBD counts: " << effIBDs
                 << '\n';

            neutrinoCounts.effectiveIBD[dataset][direction] = effIBDs;
            
            if (dataset == Data || dataset == Sim)
            {
                neutrinoCounts.mean[dataset][direction] = histogram[dataset][TotalDifference][direction]->GetMean();
                neutrinoCounts.sigma[dataset][direction] = histogram[dataset][TotalDifference][direction]->GetStdDev() / sqrt(effIBDs);
            }
        }
        cout << "--------------------------------------------\n";
    }

    CalculateUnbiasing(histogram, neutrinoCounts);
}

int main()
{
    // Ignore Warnings
    gErrorIgnoreLevel = kError;

    // Fill detector configuration
    FillDetectorConfig();

    // Check that configurations are properly set up
    cout << "Checking right neighbor for segment 44, period 4: " << checkNeighbor(4, 44, 'r') << '\n';
    cout << "Checking up neighbor for segment 82, period 3: " << checkNeighbor(3, 82, 'u') << '\n';

    // Set up what we're measuring. Check enums in header for what the ints are
    array<float, 5> phi, phiError;
    array<float, 5> theta, thetaError;

    // Need histograms for counting each variable. Check enums in header for
    // what the ints are Don't need an array for the true reactor direction
    array<array<array<std::shared_ptr<TH1D>, DirectionSize>, SignalSize>, DatasetSize> histogram;

    // Set up histograms for all 3 directions
    for (int i = Data; i < DatasetSize; i++)  // Dataset
    {
        for (int j = CorrelatedReactorOn; j < SignalSize; j++)  // Signal set
        {
            // No reactor off for simulations
            if ((i == Sim || i == SimUnbiased) && (j == CorrelatedReactorOff || j == AccidentalReactorOff))
                continue;

            for (int c = x; c < DirectionSize; c++)
            {
                string dataset = DatasetToString(i);
                string signalSet = SignalToString(j);
                string axis = AxisToString(c);
                string histogramName = dataset + "_" + signalSet + "_" + axis;
                histogram[i][j][c]
                    = std::make_shared<TH1D>(histogramName.c_str(), dataset.c_str(), bins, -histogramMax, histogramMax);
            }
        }
    }

    // Filling data histograms
    cout << "Filling data histograms!\n";
    cout << "--------------------------------------------\n";
    for (int period = 1; period < 6; period++)  // 5 periods of PROSPECT Data
    {
        SetUpHistograms(histogram, Data, period);
    }
    lineCounter = 0;

    cout << "Successfully filled data histogram!\n";
    cout << "--------------------------------------------\n";
    cout << "Total livetime for all Reactor Off events: " << livetimeOff << '\n';
    cout << "Total livetime for all Reactor On events: " << livetimeOn << '\n';
    cout << "--------------------------------------------\n";

    // Filling simulation histograms
    cout << "Filling simulation histograms!\n";
    cout << "--------------------------------------------\n";
    for (int period = 1; period < 6; period++)
    {
        SetUpHistograms(histogram, Sim, period);
    }

    cout << "Successfully filled simulation histogram!\n";

    SubtractBackgrounds(histogram);

    // Set up our output file
    /* auto outputFile = std::make_unique<TFile>("Directionality.root",
    "recreate"); outputFile->Close(); */

    return 0;
}