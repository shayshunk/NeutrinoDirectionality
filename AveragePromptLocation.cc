#include "AveragePromptLocation.h"

#include "Formatting.h"

using std::cout, std::string, std::ifstream, std::array, std::getline;

void FillHistogram(array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>& histogram, TreeValues& currentEntry)
{
    // Applying energy cut
    if (currentEntry.Esmear < 0.8 || currentEntry.Esmear > 7.4)
    {
        return;
    }

    if (currentEntry.nCaptTime > pow(10, 3) && currentEntry.nCaptTime < 120 * pow(10, 3))  // Correlated Dataset
    {
        float promptPosition;

        switch (currentEntry.direction)
        {
            case X:
                promptPosition = (currentEntry.promptSegment % 14);
                correlatedxCounts++;
                break;
            case Y:
                promptPosition = floor(currentEntry.promptSegment / 14);
                correlatedyCounts++;
                break;
            case Z:
                promptPosition = currentEntry.promptPosition;
                correlatedzCounts++;
                break;
            default:
                promptPosition = currentEntry.promptPosition;
        }

        // Fill regular dataset with prompt location
        histogram[CorrelatedReactorOn][currentEntry.direction]->Fill(promptPosition);
    }
    else if (currentEntry.nCaptTime > pow(10, 6))  // Accidental Dataset
    {
        float promptPosition;

        switch (currentEntry.direction)
        {
            case X:
                promptPosition = (currentEntry.promptSegment % 14);
                accidentalxCounts++;
                break;
            case Y:
                promptPosition = floor(currentEntry.promptSegment / 14);
                accidentalyCounts++;
                break;
            case Z:
                promptPosition = currentEntry.promptPosition;
                accidentalzCounts++;
                break;
            default:
                promptPosition = currentEntry.promptPosition;
        }

        // Fill regular dataset with prompt location
        histogram[AccidentalReactorOn][currentEntry.direction]->Fill(promptPosition);
    }
}

void SetUpHistograms(array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>& histogram, int periodNo)
{
    // Declaring some variables for use later
    int totalLines = 0;

    // Figuring out period
    string period = std::to_string(periodNo);

    // Combining names into file list name
    string fileList = Form(dataPath, std::to_string(periodNo).c_str(), std::to_string(periodNo).c_str());

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

        if (lineCounter % 1 == 0)
        {
            cout << "Reading file: " << lineCounter << "/" << totalSimLines << '\r';
            cout.flush();
        }

        // Reading file list
        string line;
        getline(file, line);

        // Combining names into root file name
        TString rootFilename = Form(dataFileName, period.c_str(), line.data());

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
        }

        // Grab rootTree and cast to unique pointer
        TTree* rootTree = (TTree*)rootFile->Get("P2kIBDPlugin/Tibd");

        long nEntries = rootTree->GetEntries();

        for (long i = 0; i < nEntries; i++)
        {
            rootTree->GetEntry(i);

            for (int direction = X; direction < DirectionSize; direction++)
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
                currentEntry.direction = direction;

                if (direction == X && correlatedxCounts == 200000)
                    continue;
                else if (direction == Y && correlatedyCounts == 200000)
                    continue;
                if (direction == Z && correlatedzCounts == 200000)
                    continue;

                FillHistogram(histogram, currentEntry);
            }
        }
        // Returns the next character in the input sequence, without extracting it: The character is left as the next character
        // to be extracted from the stream
        file.peek();
        rootFile->Close();
    }
}

void SubtractBackgrounds(array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>& histogram)
{
    /* IBD events = (Correlated - Accidental/100)_{reactor on} + (-livetimeOn/livetimeOff*Correlated +
    livetimeOn/livetimeOff*Accidental/100)_{reactor off} */

    // Defining variables for IBD background subtraction
    double totalIBDs = 0, totalIBDErr = 0, effIBDs = 0;

    for (int direction = X; direction < DirectionSize; direction++)
    {
        // Have to static cast raw pointer to shared pointer to keep up safety measures
        histogram[TotalDifference][direction]
            = std::shared_ptr<TH1F>(static_cast<TH1F*>(histogram[CorrelatedReactorOn][direction]->Clone("Displacements")));
        histogram[TotalDifference][direction]->Add(histogram[AccidentalReactorOn][direction].get(), -1. / 100.);

        totalIBDs = histogram[TotalDifference][direction]->IntegralAndError(
            0, histogram[TotalDifference][direction]->GetNbinsX() + 1, totalIBDErr);
        effIBDs = pow(totalIBDs, 2) / pow(totalIBDErr, 2);  // Effective IBD counts. Done by Poisson Distribution
                                                            // N^2/(sqrt(N)^2) = N; Eff. counts = counts^2/counts_err^2

        float mean = histogram[TotalDifference][direction]->GetMean();
        float sigma = histogram[TotalDifference][direction]->GetStdDev() / sqrt(effIBDs);

        cout << boldOn << underlineOn << "Values for: " << resetFormats << AxisToString(direction) << '\n';
        cout << "Average prompt location: " << mean << " ± " << sigma << '\n';
        cout << "Total IBDs: " << totalIBDs << " ± " << totalIBDErr << '\n';
        cout << "Effective IBDs: " << effIBDs << '\n';
    }
}

int main()
{
    // Ignore Warnings
    gErrorIgnoreLevel = kError;

    // Need histograms for counting each variable
    array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize> histogram;

    // Set up histograms for all 3 directions
    for (int j = CorrelatedReactorOn; j < SignalSize; j++)  // Signal set
    {
        for (int c = X; c < DirectionSize; c++)
        {
            string signalSet = SignalToString(j);
            string axis = AxisToString(c);
            string histogramName = signalSet + "_" + axis;

            int bins;
            float histogramMin, histogramMax;

            switch (c)
            {
                case X:
                    bins = xBins;
                    histogramMin = 0.5;
                    histogramMax = xBins + 0.5;
                    break;
                case Y:
                    bins = yBins;
                    histogramMin = 0.5;
                    histogramMax = yBins + 0.5;
                    break;
                case Z:
                    bins = zBins;
                    histogramMin = -zHistogramMax;
                    histogramMax = zHistogramMax;
                    break;
                default:
                    cout << "Mistake!\n";
            }

            histogram[j][c] = std::make_shared<TH1F>(histogramName.c_str(), signalSet.c_str(), bins, histogramMin, histogramMax);
        }
    }

    // Filling histograms
    cout << "--------------------------------------------\n";
    cout << "Filling simulation histograms!\n";
    cout << "--------------------------------------------\n";

    for (int period = 1; period <= 5; period++)
    {
        SetUpHistograms(histogram, period);

        cout << "Correlated x counts for period " << period << " : " << correlatedxCounts << '\n';
        cout << "Accidental x counts for period " << period << " : " << accidentalxCounts << '\n';
        cout << "Correlated y counts for period " << period << " : " << correlatedyCounts << '\n';
        cout << "Accidental y counts for period " << period << " : " << accidentalyCounts << '\n';
        cout << "Correlated z counts for period " << period << " : " << correlatedzCounts << '\n';
        cout << "Accidental z counts for period " << period << " : " << accidentalzCounts << '\n';

        correlatedxCounts = 0;
        correlatedyCounts = 0;
        correlatedzCounts = 0;
        accidentalxCounts = 0;
        accidentalyCounts = 0;
        accidentalzCounts = 0;
    }

    cout << boldOn << cyanOn << "Successfully filled simulation histogram!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    SubtractBackgrounds(histogram);

    return 0;
}