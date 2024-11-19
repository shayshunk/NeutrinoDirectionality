#include "MobilityTiming.h"

#include "DetectorConfig.h"
#include "Formatting.h"
#include "Timer.h"

using std::cout, std::string, std::ifstream, std::array, std::getline;

Directionality::Directionality()
{
    for (int dataset = Data; dataset < DatasetSize; dataset++)  // Dataset
    {
        for (int signalSet = CorrelatedReactorOn; signalSet < TotalDifference; signalSet++)  // Signal set
        {
            for (int direction = X; direction < DirectionSize; direction++)
            {
                // No reactor off for simulations
                if ((dataset == Sim || dataset == SimUnbiased)
                    && (signalSet == CorrelatedReactorOff || signalSet == AccidentalReactorOff))
                    continue;

                string data = DatasetToString(dataset);
                string signal = SignalToString(signalSet);
                string axis = AxisToString(direction);
                string histogramName = data + " " + signal + " " + axis;
                histogram[dataset][signalSet][direction]
                    = TH1F(histogramName.c_str(), data.c_str(), bins, -histogramMax, histogramMax);
            }
        }
    }

    timingHistogram = TH1F("Timing", "timing", 1000, pow(10, 6), pow(10, 8));

    ResetLineNumber();
}

void Directionality::ReadFileList(int dataset, int periodNo)
{
    char const* path;

    switch (dataset)
    {
        case Data:
            path = dataPath;
            break;
        case Sim:
            path = simPath;
            break;
        case PerfSim:
            path = perfSimPath;
            break;
        default:
            path = dataPath;
            break;
    }

    // Combining names into file list name
    string fileList;

    if (periodNo != 6)
        fileList = Form(path, std::to_string(periodNo).c_str(), std::to_string(periodNo).c_str());
    else
        fileList = path;

    // Opening and checking file list
    ifstream file;
    file.open(fileList, ifstream::in);

    if (!(file.is_open() && file.good()))
    {
        cout << "File list not found! Exiting.\n";
        cout << "Trying to find: " << fileList << '\n';
        return;
    }

    while (file.good() && getline(file, files[lineNumber]))
    {
        lineNumber++;
    }

    files[lineNumber] = "Done";
    lineNumber++;
}

void Directionality::SetUpHistograms(int dataset, int periodNo)
{
    // Declaring some variables for use later
    int totalLines = 0;
    reactorOn = true;

    dataSet = dataset;
    period = periodNo;

    char const* fileName;

    switch (dataset)
    {
        case Data:
            fileName = dataFileName;
            break;
        case Sim:
            fileName = simFileName;
            break;
        case PerfSim:
            fileName = perfSimFileName;
            break;
        default:
            fileName = dataFileName;
            break;
    }

    while (index < files.size())
    {
        if (files[index] == "Done")
        {
            index++;
            return;
        }

        if (dataSet == Data || dataSet == DataUnbiased)
        {
            if (lineCounter % 100 == 0)
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
        else if (dataSet == PerfSim || dataSet == PerfSimUnbiased)
        {
            cout << "Reading file: " << lineCounter << "/" << totalPerfSimLines << '\r';
            cout.flush();
        }

        // Combining names into root file name
        TString rootFilename;
        if (period != 6)
            rootFilename = Form(fileName, std::to_string(period).c_str(), files[index].data());
        else
            rootFilename = Form(fileName, files[index].data());

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

            for (direction = X; direction < DirectionSize; direction++)
            {
                // Grabbing relevant values from the rootTree entry
                promptPosition = rootTree->GetLeaf("xyz")->GetValue(direction);
                delayedPosition = rootTree->GetLeaf("n_xyz")->GetValue(direction);
                promptSegment = rootTree->GetLeaf("maxseg")->GetValue(0);
                delayedSegment = rootTree->GetLeaf("n_seg")->GetValue(0);

                // We throw out events where the neutron moves in a direction we're not checking
                if (promptSegment != delayedSegment && promptPosition == delayedPosition)
                    continue;

                if (promptSegment == 32 || delayedSegment == 32)
                    continue;

                // Copying some loop values into current entry
                Esmear = rootTree->GetLeaf("Esmear")->GetValue(0);
                nCaptTime = rootTree->GetLeaf("ncapt_dt")->GetValue(0);

                FillHistogram();
            }
        }
        // Returns the next character in the input sequence, without extracting it: The character is left as the next character
        // to be extracted from the stream
        rootFile->Close();

        lineCounter++;
        index++;
    }
}

void Directionality::FillHistogram()
{
    // Applying energy cut
    if (Esmear < 0.8 || Esmear > 7.4)
    {
        return;
    }

    int signalSet = 0;

    if (nCaptTime > 20 * pow(10, 3) && nCaptTime < 120 * pow(10, 3))  // Correlated Dataset
    {
        // Calculate neutron displacement
        double displacement = delayedPosition - promptPosition;

        // Figure out whether the reactor is on and assign signal index
        signalSet = reactorOn ? CorrelatedReactorOn : CorrelatedReactorOff;

        // Fill regular dataset with displacement but Z only takes same segment
        if (direction != Z)
        {
            // Correcting center to center displacment
            if (displacement > 0)
                displacement = segmentWidth;
            else if (displacement < 0)
                displacement = -segmentWidth;

            histogram[dataSet][signalSet][direction].Fill(displacement);
        }
        else if (direction == Z)
        {
            double displacement = delayedPosition - promptPosition;
            histogram[dataSet][signalSet][direction].Fill(displacement);
        }
    }
    else if (nCaptTime > 7 * pow(10, 6))  // Accidental Dataset
    {
        // Calculate neutron displacement
        double displacement = delayedPosition - promptPosition;

        // Figure out whether the reactor is on and assign signal index
        signalSet = reactorOn ? AccidentalReactorOn : AccidentalReactorOff;

        // Fill regular dataset with displacement but Z only takes same segment
        if (direction != Z)
        {
            // Correcting center to center displacment
            if (displacement > 0)
                displacement = segmentWidth;
            else if (displacement < 0)
                displacement = -segmentWidth;

            histogram[dataSet][signalSet][direction].Fill(displacement, xRx);
        }
        else if (direction == Z)
        {
            double displacement = delayedPosition - promptPosition;
            histogram[dataSet][signalSet][direction].Fill(displacement, xRx);
        }

        if (maxAccidental < nCaptTime)
        {
            maxAccidental = nCaptTime;
        }

        timingHistogram.Fill(nCaptTime);
    }
}

void Directionality::SubtractBackgrounds()
{
    /* IBD events = (Correlated - Accidental/100)_{reactor on} + (-livetimeOn/livetimeOff*Correlated +
    livetimeOn/livetimeOff*Accidental/100)_{reactor off} */

    // Defining variables for IBD background subtraction
    double totalIBDs = 0, totalIBDErr = 0, effIBDs = 0;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        for (int direction = X; direction < DirectionSize; direction++)
        {
            string histogramName;
            string data = DatasetToString(dataset);
            histogramName = data + " Total Difference " + AxisToString(direction);

            // Copying Correlated Reactor On to start
            histogram[dataset][TotalDifference][direction] = TH1F(histogram[dataset][CorrelatedReactorOn][direction]);
            histogram[dataset][TotalDifference][direction].SetNameTitle(histogramName.c_str(), data.c_str());

            histogram[dataset][TotalDifference][direction].Add(&histogram[dataset][AccidentalReactorOn][direction], -1. / 100);

            if (dataset == Data || dataset == DataUnbiased)
            {
                histogram[dataset][TotalDifference][direction].Add(&histogram[dataset][CorrelatedReactorOff][direction],
                                                                   -livetimeOn * atmosphericScaling / livetimeOff);

                histogram[dataset][TotalDifference][direction].Add(&histogram[dataset][AccidentalReactorOff][direction],
                                                                   livetimeOn * atmosphericScaling / (100 * livetimeOff));
            }

            totalIBDs = histogram[dataset][TotalDifference][direction].IntegralAndError(
                0, histogram[dataset][TotalDifference][direction].GetNbinsX() + 1, totalIBDErr);
            effIBDs = pow(totalIBDs, 2) / pow(totalIBDErr, 2);  // Effective IBD counts. Done by Poisson Distribution
                                                                // N^2/(sqrt(N)^2) = N; Eff. counts = counts^2/counts_err^2

            totalIBD[dataset][direction] = totalIBDs;
            totalIBDError[dataset][direction] = totalIBDErr;

            if (direction == Z)
            {
                effectiveIBD[dataset][direction] = effIBDs;
                continue;
            }

            if (dataset == Data || dataset == Sim || dataset == PerfSim)
            {
                mean[dataset][direction] = histogram[dataset][TotalDifference][direction].GetMean();
                sigma[dataset][direction] = histogram[dataset][TotalDifference][direction].GetStdDev() / sqrt(effIBDs);
                effectiveIBD[dataset][direction] = effIBDs;
            }
            else
            {
                effectiveIBD[dataset][direction] = effIBDs;
            }
        }

        if (dataset == DataUnbiased || dataset == SimUnbiased)
            continue;

        // Z is fit to a Guassian and only takes same segment inputs
        // Possible thanks to 1mm resolution in Z
        TF1 gaussian("Fit", "gaus", -140, 140);

        histogram[dataset][TotalDifference][Z].Fit("Fit", "RQ");

        float zMean = gaussian.GetParameter(1);
        float zError = gaussian.GetParError(1);

        mean[dataset][Z] = zMean;
        sigma[dataset][Z] = zError;

        // Deleting fit because we don't want the plot options stuck here
        delete histogram[dataset][TotalDifference][Z].GetListOfFunctions()->FindObject("Fit");
    }

    // Printing out values
    if (IBDCOUNT_VERBOSITY)
    {
        for (int dataset = Data; dataset < DatasetSize; dataset++)
        {
            cout << "Total and Effective IBD Events for: " << boldOn << DatasetToString(dataset) << resetFormats << '\n';

            for (int direction = X; direction < DirectionSize; direction++)
            {
                cout << boldOn << AxisToString(direction) << ": " << resetFormats << totalIBD[dataset][direction] << " ± "
                     << totalIBDError[dataset][direction] << boldOn << ". Effective IBD counts: " << resetFormats
                     << effectiveIBD[dataset][direction] << ".\n";
            }

            cout << "--------------------------------------------\n";
        }
    }

    cout << boldOn << cyanOn << "Subtracted backgrounds.\n" << resetFormats;
    cout << "--------------------------------------------\n";
}

void Directionality::FillOutputFile()
{
    // Set up our output file
    TFile outputFile("Mobility.root", "recreate");

    outputFile.cd();

    // Printing histograms to file just in case
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        for (int signalSet = CorrelatedReactorOn; signalSet < SignalSize; signalSet++)
        {
            // No reactor off for simulations
            if ((dataset == Sim || dataset == SimUnbiased)
                && (signalSet == CorrelatedReactorOff || signalSet == AccidentalReactorOff))
                continue;

            for (int direction = X; direction < DirectionSize; direction++)
            {
                histogram[dataset][signalSet][direction].Write();
            }
        }
    }

    timingHistogram.Write();

    cout << boldOn << cyanOn << "Filled output file: " << resetFormats << blueOn << boldOn << "Mobility.root!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    outputFile.Close();
}

int main(int argc, char* argv[])
{
    // Ignore Warnings (mostly for time honestly)
    gErrorIgnoreLevel = kError;

    // Take ownership of histograms
    TH1::AddDirectory(kFALSE);

    // Setting up timer
    Timer timer;

    // Setting up Directionality class
    Directionality neutrinoDirection;

    // Testing detector config
    /* cout << "Checking right neighbor for segment 25: " << CheckNeighbor(6, 25, 'r') << '\n';
    cout << "Checking left neighbor for segment 25: " << CheckNeighbor(6, 25, 'l') << '\n';
    cout << "Checking up neighbor for segment 25: " << CheckNeighbor(6, 25, 'u') << '\n';
    cout << "Checking down neighbor for segment 25: " << CheckNeighbor(6, 25, 'd') << '\n'; */

    // Filling data histograms
    cout << "--------------------------------------------\n";
    cout << "Filling data histograms!\n";
    cout << "--------------------------------------------\n";
    for (int period = 1; period <= 5; period++)  // 5 periods of PROSPECT Data
    {
        neutrinoDirection.ReadFileList(Data, period);
    }

    for (int period = 1; period <= 5; period++)  // 5 periods of PROSPECT Data
    {
        neutrinoDirection.SetUpHistograms(Data, period);
    }

    neutrinoDirection.ResetLineCounter();
    neutrinoDirection.ResetLineNumber();
    neutrinoDirection.ResetIndex();

    cout << boldOn << cyanOn << "Successfully filled data histogram!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    if (LIVETIME_VERBOSITY)
    {
        cout << "Total livetime for all" << boldOn << " Reactor Off " << resetFormats
             << "events: " << neutrinoDirection.livetimeOff << '\n';
        cout << "Total livetime for all" << boldOn << " Reactor On " << resetFormats
             << "events: " << neutrinoDirection.livetimeOn << '\n';
        cout << "--------------------------------------------\n";
    }

    // Filling simulation histograms
    cout << "Filling simulation histograms!\n" << resetFormats;
    cout << "--------------------------------------------\n";
    for (int period = 1; period <= 5; period++)  // 5 periods of PROSPECT Data
    {
        neutrinoDirection.ReadFileList(Sim, period);
    }
    for (int period = 1; period <= 5; period++)
    {
        neutrinoDirection.SetUpHistograms(Sim, period);
    }

    cout << boldOn << cyanOn << "Successfully filled simulation histogram!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    neutrinoDirection.ResetLineCounter();
    neutrinoDirection.ResetLineNumber();
    neutrinoDirection.ResetIndex();

    // Filling simulation histograms
    cout << "Filling perfect simulation histograms!\n" << resetFormats;
    cout << "--------------------------------------------\n";
    neutrinoDirection.ReadFileList(PerfSim, 6);
    neutrinoDirection.SetUpHistograms(PerfSim, 6);

    cout << boldOn << cyanOn << "Successfully filled perfect simulation histogram!\n" << resetFormats;
    cout << "--------------------------------------------\n";
    // Should be self-explanatory
    neutrinoDirection.SubtractBackgrounds();
    neutrinoDirection.FillOutputFile();

    cout << boldOn << cyanOn << "Maximum accidental time window: " << resetFormats << neutrinoDirection.maxAccidental;
    cout << "--------------------------------------------\n";

    return 0;
}