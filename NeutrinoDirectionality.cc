#include "NeutrinoDirectionality.h"

#include "DetectorConfig.h"
#include "Formatting.h"

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

    if (DETECTOR_VERBOSITY)
    {
        cout << "--------------------------------------------\n";
        cout << "Below is the detector configuration.\n";
        cout << "--------------------------------------------\n";

        for (int i = 0; i < noPeriods; i++)
        {
            cout << "Detector configuration for period: " << i + 1 << '\n';

            for (int j = 140; j >= 0; j -= 14)
            {
                for (int k = 0; k < 14; k++)
                {
                    if (detectorConfig[i][j + k])
                    {
                        cout << "\u25A0 ";
                    }
                    else
                    {
                        cout << "\u25A1 ";
                    }
                }

                cout << '\n';
            }
            cout << '\n';
        }
    }
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

void FillHistogramUnbiased(array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize>& histogram,
                           TreeValues& currentEntry,
                           int signalSet)
{
    bool posDirection = false, negDirection = false, success = false;

    // Need to weight accidental datasets by deadtime correction factor
    double weight = (signalSet == AccidentalReactorOff || signalSet == AccidentalReactorOn) ? currentEntry.xRx : 1;

    // Check for live neighbors in different directions based on which axis
    // we're filling
    if (currentEntry.direction == X)
    {
        posDirection = checkNeighbor(currentEntry.period, currentEntry.promptSegment, 'r');
        negDirection = checkNeighbor(currentEntry.period, currentEntry.promptSegment, 'l');
    }
    else if (currentEntry.direction == Y)
    {
        posDirection = checkNeighbor(currentEntry.period, currentEntry.promptSegment, 'u');
        negDirection = checkNeighbor(currentEntry.period, currentEntry.promptSegment, 'd');
    }
    else  // Fill Z with same segment events
    {
        double displacement = currentEntry.delayedPosition - currentEntry.promptPosition;
        histogram[currentEntry.dataSet][signalSet][currentEntry.direction]->Fill(displacement, weight);
    }

    // Dataset + 1 returns the unbiased version of that dataset
    if (posDirection && !negDirection)
        histogram[currentEntry.dataSet + 1][signalSet][currentEntry.direction]->Fill(segmentWidth, weight);
    else if (!posDirection && negDirection)
        histogram[currentEntry.dataSet + 1][signalSet][currentEntry.direction]->Fill(-segmentWidth, weight);
    else if (posDirection && negDirection)
        histogram[currentEntry.dataSet + 1][signalSet][currentEntry.direction]->Fill(0.0, weight);
}

void FillHistogram(array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize>& histogram,
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

        // Fill regular dataset with displacement but Z only takes same segment
        if (currentEntry.direction != Z)
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

        // Fill regular dataset with displacement but Z only takes same segment
        if (currentEntry.direction != Z)
            histogram[currentEntry.dataSet][signalSet][currentEntry.direction]->Fill(displacement, currentEntry.xRx);

        // Fill dead segment correction dataset
        if (currentEntry.promptSegment == currentEntry.delayedSegment)
        {
            FillHistogramUnbiased(histogram, currentEntry, signalSet);
        }
    }
}

void SetUpHistograms(array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize>& histogram,
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

                if (direction != Z)  // Cubical distance cuts in Z for X and Y
                {
                    float zPromptPosition = rootTree->GetLeaf("xyz")->GetValue(Z);
                    float zDelayedPosition = rootTree->GetLeaf("n_xyz")->GetValue(Z);

                    float displacement = zPromptPosition - zDelayedPosition;

                    if (displacement > 60)
                        continue;
                }

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

void CalculateUnbiasing(array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize>& histogram,
                        IBDValues& neutrinoCounts)
{
    // Defining variables used in calculation. Check the error propagation technote for details on the method
    double rPlus = 0, rMinus = 0;
    double p = 0, pError = 0;
    double nPlus = 0, nPlusPlus = 0, nMinus = 0, nMinusMinus = 0, nPlusMinus = 0;
    double nPlusError = 0, nPlusPlusError = 0, nMinusError = 0, nMinusMinusError = 0, nPlusMinusError = 0;

    for (int dataset = DataUnbiased; dataset < DatasetSize; dataset += 2)
    {
        for (int direction = X; direction < Z; direction++)
        {
            // Dataset - 1 returns the biased dataset counts
            // Grabbing data from filled bins, rest should be empty
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

            pError
                = segmentWidth
                  * pow(1 / ((nMinus * (nPlusMinus + nPlusPlus) + (nMinusMinus + nPlusMinus) * (nPlus + nPlusMinus + nPlusPlus))),
                        2)
                  * sqrt(
                      pow((nMinusMinus + nPlusMinus) * (nPlusMinus + nPlusPlus), 2)
                          * (pow(nPlusError * (2 * nMinus + nMinusMinus + nPlusMinus), 2)
                             + pow(nMinusError * (2 * nPlus + nPlusPlus + nPlusMinus), 2))
                      + pow((nPlus * (nPlusMinus + nMinusMinus) * (2 * nMinus + nMinusMinus + nPlusMinus) * nPlusPlusError), 2)
                      + pow((nPlusMinusError
                             * (nPlus * pow((nMinusMinus + nPlusMinus), 2)
                                + nMinus * (2 * nMinusMinus * nPlus - 2 * nPlus * nPlusPlus - pow((nPlusMinus + nPlusPlus), 2)))),
                            2)
                      + pow((nMinus * (nPlusMinus + nPlusPlus) * (2 * nPlus + nPlusMinus + nPlusPlus) * nMinusMinusError), 2));

            neutrinoCounts.mean[dataset][direction] = p;
            neutrinoCounts.sigma[dataset][direction] = pError;
        }
        neutrinoCounts.effectiveIBD[dataset][Z] = neutrinoCounts.effectiveIBD[dataset - 1][Z];
        neutrinoCounts.mean[dataset][Z] = neutrinoCounts.mean[dataset - 1][Z];
        neutrinoCounts.sigma[dataset][Z] = neutrinoCounts.sigma[dataset - 1][Z];
    }

    cout << boldOn << cyanOn << "Calculated Means.\n" << resetFormats;
    cout << "--------------------------------------------\n";

    // Printing out values
    if (MEAN_VERBOSITY)
    {
        for (int dataset = Data; dataset < DatasetSize; dataset++)
        {
            cout << "Mean and sigma values for: " << boldOn << DatasetToString(dataset) << resetFormats << '\n';
            for (int direction = X; direction < DirectionSize; direction++)
            {
                cout << boldOn << "p" << AxisToString(direction) << ": " << resetFormats
                     << neutrinoCounts.mean[dataset][direction] << " ± " << neutrinoCounts.sigma[dataset][direction] << '\n';
            }
            cout << "--------------------------------------------\n";
        }
    }
}

IBDValues SubtractBackgrounds(array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize>& histogram)
{
    /* IBD events = (Correlated - Accidental/100)_{reactor on} + (-livetimeOn/livetimeOff*Correlated +
    livetimeOn/livetimeOff*Accidental/100)_{reactor off} */

    // Defining variables for IBD background subtraction
    double totalIBDs = 0, totalIBDErr = 0, effIBDs = 0;
    IBDValues neutrinoCounts;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        for (int direction = X; direction < DirectionSize; direction++)
        {
            string histogramName;
            histogramName = DatasetToString(dataset) + " Total Difference " + AxisToString(direction);

            // Have to static cast raw pointer to shared pointer to keep up safety measures
            histogram[dataset][TotalDifference][direction] = std::shared_ptr<TH1F>(
                static_cast<TH1F*>(histogram[dataset][CorrelatedReactorOn][direction]->Clone(histogramName.c_str())));
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

            neutrinoCounts.totalIBD[dataset][direction] = totalIBDs;
            neutrinoCounts.totalIBDError[dataset][direction] = totalIBDErr;

            if (direction == Z)
            {
                neutrinoCounts.effectiveIBD[dataset][direction] = effIBDs;
                continue;
            }

            if (dataset == Data || dataset == Sim)
            {
                neutrinoCounts.mean[dataset][direction] = histogram[dataset][TotalDifference][direction]->GetMean();
                neutrinoCounts.sigma[dataset][direction] = histogram[dataset][TotalDifference][direction]->GetStdDev()
                                                           / sqrt(effIBDs);
                neutrinoCounts.effectiveIBD[dataset][direction] = effIBDs;
            }
            else
            {
                neutrinoCounts.effectiveIBD[dataset][direction] = neutrinoCounts.effectiveIBD[dataset - 1][direction];
            }
        }

        if (dataset == DataUnbiased || dataset == SimUnbiased)
            continue;

        // Z is fit to a Guassian and only takes same segment inputs
        // Possible thanks to 1mm resolution in Z
        TF1* gaussian = new TF1("Fit", "gaus", -140, 140);

        histogram[dataset][TotalDifference][Z]->Fit("Fit", "RQ");

        float zMean = gaussian->GetParameter(1);
        float zError = gaussian->GetParError(1);

        neutrinoCounts.mean[dataset][Z] = zMean;
        neutrinoCounts.sigma[dataset][Z] = zError;

        // Deleting fit because we don't want the plot options stuck here
        delete histogram[dataset][TotalDifference][Z]->GetListOfFunctions()->FindObject("Fit");
    }

    // Printing out values
    if (IBDCOUNT_VERBOSITY)
    {
        for (int dataset = Data; dataset < DatasetSize; dataset++)
        {
            cout << "Total and Effective IBD Events for: " << boldOn << DatasetToString(dataset) << resetFormats << '\n';

            for (int direction = X; direction < DirectionSize; direction++)
            {
                cout << boldOn << AxisToString(direction) << ": " << resetFormats << neutrinoCounts.totalIBD[dataset][direction]
                     << " ± " << neutrinoCounts.totalIBDError[dataset][direction] << boldOn
                     << ". Effective IBD counts: " << resetFormats << neutrinoCounts.effectiveIBD[dataset][direction] << ".\n";
            }

            cout << "--------------------------------------------\n";
        }
    }

    cout << boldOn << cyanOn << "Subtracted backgrounds.\n" << resetFormats;
    cout << "--------------------------------------------\n";

    CalculateUnbiasing(histogram, neutrinoCounts);

    return neutrinoCounts;
}

void AddSystematics(IBDValues& neutrinoCounts)
{
    // Systematics extracted from BiPo study

    // Defining variables for readability of code
    double sigmaX, sigmaY, sigmaZ;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        sigmaX = neutrinoCounts.sigma[dataset][X];
        sigmaY = neutrinoCounts.sigma[dataset][Y];
        sigmaZ = neutrinoCounts.sigma[dataset][Z];

        neutrinoCounts.sigmaSystematics[dataset][X] = sqrt(pow(sigmaX, 2) + pow(0.25, 2) + pow(0.08, 2));
        neutrinoCounts.sigmaSystematics[dataset][Y] = sqrt(pow(sigmaY, 2) + pow(0.39, 2) + pow(0.08, 2));
        neutrinoCounts.sigmaSystematics[dataset][Z] = sqrt(pow(sigmaZ, 2) + pow(0.05, 2) + pow(0.09, 2));
    }

    cout << boldOn << cyanOn << "Added Systematics.\n" << resetFormats;
    cout << "--------------------------------------------\n";
}

AngleValues CalculateAngles(IBDValues const& neutrinoCounts)
{
    AngleValues finalAngles;

    // Defining variables for readability of code
    double px, py, pz;
    double sigmaX, sigmaY, sigmaZ;
    double sigmaXSystematics, sigmaYSystematics, sigmaZSystematics;
    double effIBDX, effIBDY, effIBDZ;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        // Grabbing values from current dataset
        px = neutrinoCounts.mean[dataset][X];
        py = neutrinoCounts.mean[dataset][Y];
        pz = neutrinoCounts.mean[dataset][Z];

        sigmaX = neutrinoCounts.sigma[dataset][X];
        sigmaY = neutrinoCounts.sigma[dataset][Y];
        sigmaZ = neutrinoCounts.sigma[dataset][Z];

        sigmaXSystematics = neutrinoCounts.sigmaSystematics[dataset][X];
        sigmaYSystematics = neutrinoCounts.sigmaSystematics[dataset][Y];
        sigmaZSystematics = neutrinoCounts.sigmaSystematics[dataset][Z];

        effIBDX = neutrinoCounts.effectiveIBD[dataset][X];
        effIBDY = neutrinoCounts.effectiveIBD[dataset][Y];
        effIBDZ = neutrinoCounts.effectiveIBD[dataset][Z];

        // phi = arctan(y / x)
        double tanPhi = py / px;
        double phi = atan(tanPhi) * 180.0 / pi;
        double tanPhiError = sqrt(pow((sigmaX * py) / pow(px, 2), 2) + pow(sigmaY / px, 2));
        double phiError = (tanPhiError / (1 + pow(tanPhi, 2))) * 180.0 / pi;
        double tanPhiErrorSystematics = sqrt(pow((sigmaXSystematics * py) / pow(px, 2), 2) + pow(sigmaYSystematics / px, 2));
        double phiErrorSystematics = (tanPhiErrorSystematics / (1 + pow(tanPhi, 2))) * 180.0 / pi;

        // theta = arctan(z / sqrt(x^2 + y^2))
        double tanTheta = pz / sqrt(pow(px, 2) + pow(py, 2));
        double theta = atan(tanTheta) * 180.0 / pi;
        double tanThetaError = sqrt((1 / (px * px + py * py))
                                    * (pow((px * pz * sigmaX / (px * px + py * py)), 2)
                                       + pow((py * pz * sigmaY / (px * px + py * py)), 2) + pow(sigmaZ, 2)));
        double thetaError = (tanThetaError / (1 + pow(tanTheta, 2))) * 180.0 / pi;
        double tanThetaErrorSystematics
            = sqrt((1 / (px * px + py * py))
                   * (pow((px * pz * sigmaXSystematics / (px * px + py * py)), 2)
                      + pow((py * pz * sigmaYSystematics / (px * px + py * py)), 2) + pow(sigmaZSystematics, 2)));
        double thetaErrorSystematics = (tanThetaErrorSystematics / (1 + pow(tanTheta, 2))) * 180.0 / pi;

        // Storing values
        finalAngles.phi[dataset] = phi;
        finalAngles.phiError[dataset] = phiError;
        finalAngles.phiErrorSystematics[dataset] = phiErrorSystematics;
        finalAngles.theta[dataset] = theta;
        finalAngles.thetaError[dataset] = thetaError;
        finalAngles.thetaErrorSystematics[dataset] = thetaErrorSystematics;
    }

    // Calculating "true" neutrino direction
    // Based on Figure 1, https://doi.org/10.1103/PhysRevD.103.032001
    float xTrue = -5.97, yTrue = -5.09, zTrue = 1.19;
    float xTrueError = 0.1, yTrueError = 0.1, zTrueError = 0.1;

    // Scaling by average prompt location
    xTrue += 64.61 / 1000;
    yTrue += 13.11 / 1000;
    zTrue += 0.73 / 1000;

    // Same angle calculation as above
    float tanPhiTrue = yTrue / xTrue;
    float phiTrue = atan(tanPhiTrue) * 180.0 / pi;
    float tanPhiTrueError = sqrt(pow((yTrue * xTrueError) / (xTrue * xTrue), 2) + pow(yTrueError / xTrue, 2));
    float phiTrueError = tanPhiTrueError / (1 + pow(tanPhiTrue, 2)) * 180.0 / pi;

    float tanThetaTrue = zTrue / sqrt(pow(xTrue, 2) + pow(yTrue, 2));
    float thetaTrue = atan(tanThetaTrue) * 180.0 / pi;
    float tanThetaTrueError
        = sqrt(pow(1 / sqrt(xTrue * xTrue + yTrue * yTrue), 2)
               * (pow(xTrue * zTrue / (xTrue * xTrue + yTrue * yTrue) * xTrueError, 2)
                  + pow(yTrue * zTrue / (xTrue * xTrue + yTrue * yTrue) * yTrueError, 2) + pow(zTrueError, 2)));
    float thetaTrueError = tanThetaTrueError / (1 + pow(tanPhiTrue, 2)) * 180.0 / pi;

    finalAngles.phiTrue = phiTrue;
    finalAngles.phiTrueError = phiTrueError;
    finalAngles.thetaTrue = thetaTrue;
    finalAngles.thetaTrueError = thetaTrueError;

    return finalAngles;
}

CovarianceValues CalculateCovariances(IBDValues const& neutrinoCounts, AngleValues const& finalAngles)
{
    CovarianceValues oneSigmaEllipse;

    // Calculating covariances
    array<array<float, 2>, 2> covarianceMatrix;
    array<array<float, 2>, 2> covarianceMatrixSystematics;

    // Defining variables for readability
    float px, py, pz;
    float sigmaX, sigmaY, sigmaZ;
    float sigmaXSystematics, sigmaYSystematics, sigmaZSystematics;
    float phi, theta;

    // From error propagation document
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        // Grabbing values of current dataset from struct
        px = neutrinoCounts.mean[dataset][X];
        py = neutrinoCounts.mean[dataset][Y];
        pz = neutrinoCounts.mean[dataset][Z];

        sigmaX = neutrinoCounts.sigma[dataset][X];
        sigmaY = neutrinoCounts.sigma[dataset][Y];
        sigmaZ = neutrinoCounts.sigma[dataset][Z];

        sigmaXSystematics = neutrinoCounts.sigmaSystematics[dataset][X];
        sigmaYSystematics = neutrinoCounts.sigmaSystematics[dataset][Y];
        sigmaZSystematics = neutrinoCounts.sigmaSystematics[dataset][Z];

        phi = finalAngles.phi[dataset];
        theta = finalAngles.theta[dataset];

        // Filling out first matrix
        covarianceMatrix[0][0] = (pow(sigmaX, 2) * pow(py, 2) / (pow(px, 4))) + pow(sigmaY, 2) / pow(px, 2);
        covarianceMatrix[0][1] = (pow(pz, 2) / (px * pow(pow(px, 2) + pow(py, 2), (3. / 2.))))
                                 * (pow(sigmaY, 2) - pow(sigmaX, 2));
        covarianceMatrix[1][0] = ((py * pz) / (px * pow(pow(px, 2) + pow(py, 2), (3. / 2.)))) * (pow(sigmaY, 2) - pow(sigmaX, 2));
        covarianceMatrix[1][1] = ((pow(px, 2) * pow(pz, 2) * pow(sigmaX, 2) + pow(py, 2) * pow(pz, 2) * pow(sigmaY, 2))
                                  / (pow(pow(px, 2) + pow(py, 2), 3)))
                                 + (pow(sigmaZ, 2)) / (pow(px, 2) + pow(py, 2));

        covarianceMatrixSystematics[0][0] = (pow(sigmaXSystematics, 2) * pow(py, 2) / (pow(px, 4)))
                                            + pow(sigmaYSystematics, 2) / pow(px, 2);
        covarianceMatrixSystematics[0][1] = (pow(pz, 2) / (px * pow(pow(px, 2) + pow(py, 2), (3. / 2.))))
                                            * (pow(sigmaYSystematics, 2) - pow(sigmaXSystematics, 2));
        covarianceMatrixSystematics[1][0] = ((py * pz) / (px * pow(pow(px, 2) + pow(py, 2), (3. / 2.))))
                                            * (pow(sigmaYSystematics, 2) - pow(sigmaXSystematics, 2));
        covarianceMatrixSystematics[1][1]
            = ((pow(px, 2) * pow(pz, 2) * pow(sigmaXSystematics, 2) + pow(py, 2) * pow(pz, 2) * pow(sigmaYSystematics, 2))
               / (pow(pow(px, 2) + pow(py, 2), 3)))
              + (pow(sigmaZSystematics, 2)) / (pow(px, 2) + pow(py, 2));

        // Including angles
        covarianceMatrix[0][0] = covarianceMatrix[0][0] * pow(cos(phi * pi / 180), 4);
        covarianceMatrix[0][1] = covarianceMatrix[0][1] * pow(cos(phi * pi / 180), 2) * pow(cos(theta * pi / 180), 2);
        covarianceMatrix[1][0] = covarianceMatrix[1][0] * pow(cos(phi * pi / 180), 2) * pow(cos(theta * pi / 180), 2);
        covarianceMatrix[1][1] = covarianceMatrix[1][1] * pow(cos(theta * pi / 180), 4);

        covarianceMatrixSystematics[0][0] = covarianceMatrixSystematics[0][0] * pow(cos(phi * pi / 180), 4);
        covarianceMatrixSystematics[0][1] = covarianceMatrixSystematics[0][1] * pow(cos(phi * pi / 180), 2)
                                            * pow(cos(theta * pi / 180), 2);
        covarianceMatrixSystematics[1][0] = covarianceMatrixSystematics[1][0] * pow(cos(phi * pi / 180), 2)
                                            * pow(cos(theta * pi / 180), 2);
        covarianceMatrixSystematics[1][1] = covarianceMatrixSystematics[1][1] * pow(cos(theta * pi / 180), 4);

        // Eigenvalues
        float a = covarianceMatrix[0][0], aSystematics = covarianceMatrixSystematics[0][0];
        float b = covarianceMatrix[0][1], bSystematics = covarianceMatrixSystematics[0][1];
        float c = covarianceMatrix[1][0], cSystematics = covarianceMatrixSystematics[1][0];
        float d = covarianceMatrix[1][1], dSystematics = covarianceMatrixSystematics[1][1];

        float lambda1 = ((a + d) + sqrt(pow(a, 2) - 2 * a * d + 4 * b * c + pow(d, 2))) / 2;
        float lambda2 = ((a + d) - sqrt(pow(a, 2) - 2 * a * d + 4 * b * c + pow(d, 2))) / 2;

        float lambda1Sytematics = ((aSystematics + dSystematics)
                                   + sqrt(pow(aSystematics, 2) - 2 * aSystematics * dSystematics
                                          + 4 * bSystematics * cSystematics + pow(dSystematics, 2)))
                                  / 2;
        float lambda2Sytematics = ((aSystematics + dSystematics)
                                   - sqrt(pow(aSystematics, 2) - 2 * aSystematics * dSystematics
                                          + 4 * bSystematics * cSystematics + pow(dSystematics, 2)))
                                  / 2;

        // Eigenvector to calculate tilt
        float vector1_1 = lambda1 - d;
        float vector1_2 = c;

        float normalizer = 1.0 / c;
        vector1_1 *= normalizer;
        vector1_2 *= normalizer;

        float tilt = atan(vector1_1) * 180.0 / pi;
        tilt = 90 + tilt;

        float vector1_1Systematics = lambda1Sytematics - dSystematics;
        float vector1_2Systematics = cSystematics;

        float normalizerSystematics = 1.0 / cSystematics;
        vector1_1Systematics *= normalizerSystematics;
        vector1_2Systematics *= normalizerSystematics;

        float tiltSystematics = atan(vector1_1Systematics) * 180.0 / pi;
        tiltSystematics = 90 + tiltSystematics;

        // Calculating final angle errors
        float phiError = sqrt(2.291 * lambda1) * 180.0 / pi;
        float thetaError = sqrt(2.291 * lambda2) * 180.0 / pi;

        float phiErrorSystematics = sqrt(2.291 * lambda1Sytematics) * 180.0 / pi;
        float thetaErrorSystematics = sqrt(2.291 * lambda2Sytematics) * 180.0 / pi;

        // Filling struct
        oneSigmaEllipse.phiError[dataset] = phiError;
        oneSigmaEllipse.phiErrorSystematics[dataset] = phiErrorSystematics;
        oneSigmaEllipse.thetaError[dataset] = thetaError;
        oneSigmaEllipse.thetaErrorSystematics[dataset] = thetaErrorSystematics;
        oneSigmaEllipse.tilt[dataset] = tilt;
        oneSigmaEllipse.tiltSystematics[dataset] = tiltSystematics;
    }

    // Prints out the 1 sigma values if COVARIANCE_VERBOSITY is set to 1
    // Change at top of file
    if (COVARIANCE_VERBOSITY)
    {
        for (int dataset = Data; dataset < DatasetSize; dataset++)
        {
            cout << "The 1 sigma ellipse for: " << boldOn << DatasetToString(dataset) << resetFormats << " with systematics.\n";
            cout << greenOn;
            cout << boldOn << underlineOn << "Phi:" << resetFormats << greenOn << " " << finalAngles.phi[dataset] << "\u00B0 ± "
                 << oneSigmaEllipse.phiErrorSystematics[dataset] << "\u00B0.\n";
            cout << boldOn << underlineOn << "Theta:" << resetFormats << greenOn << " " << finalAngles.theta[dataset]
                 << "\u00B0 ± " << oneSigmaEllipse.thetaErrorSystematics[dataset] << "\u00B0.\n";
            cout << boldOn << underlineOn << "Tilt:" << resetFormats << greenOn << " "
                 << oneSigmaEllipse.tiltSystematics[dataset] << "\u00B0.\n" << resetFormats;
            cout << "--------------------------------------------\n";

            cout << "The 1 sigma ellipse for: " << boldOn << DatasetToString(dataset) << resetFormats
                 << " without systematics.\n";
            cout << greenOn;
            cout << boldOn << underlineOn << "Phi:" << resetFormats << greenOn << " " << finalAngles.phi[dataset] << "\u00B0 ± "
                 << oneSigmaEllipse.phiError[dataset] << "\u00B0.\n";
            cout << boldOn << underlineOn << "Theta:" << resetFormats << greenOn << " " << finalAngles.theta[dataset]
                 << "\u00B0 ± " << oneSigmaEllipse.thetaError[dataset] << "\u00B0.\n";
            cout << boldOn << underlineOn << "Tilt:" << resetFormats << greenOn << " " << oneSigmaEllipse.tilt[dataset]
                 << "\u00B0.\n" << resetFormats;
            cout << "--------------------------------------------\n";
        }
    }

    // Final cone of uncertainty
    float a = oneSigmaEllipse.phiErrorSystematics[DataUnbiased];
    float b = oneSigmaEllipse.thetaErrorSystematics[DataUnbiased];
    theta = 90 - finalAngles.theta[DataUnbiased];
    float solidAngle = pi * a * b * sin(theta * pi / 180.0);
    float solidAngleRadians = solidAngle * pow((pi / 180.0), 2);
    float coneAngle = acos(1 - (solidAngleRadians / (2 * pi))) * 180.0 / pi;

    cout << boldOn << greenOn << underlineOn << "Cone of Uncertainty:" << resetFormats << greenOn << " " << coneAngle << "\u00B0"
         << '\n' << resetFormats; 
    cout << "--------------------------------------------\n";

    return oneSigmaEllipse;
}

void OffsetTheta(AngleValues& finalAngles)
{
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        finalAngles.theta[dataset] = 90 - finalAngles.theta[dataset];
    }

    finalAngles.thetaTrue = 90 + finalAngles.thetaTrue;
}

void PrintAngles(AngleValues const& finalAngles)
{
    cout << boldOn << cyanOn << "Final Angles!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        float phiError, thetaError;
        if (dataset == Sim || dataset == SimUnbiased)
        {
            phiError = finalAngles.phiError[dataset];
            thetaError = finalAngles.thetaError[dataset];
        }
        else if (ANGLES_STATISTICS)
        {
            phiError = finalAngles.phiError[dataset];
            thetaError = finalAngles.thetaError[dataset];
        }
        else
        {
            phiError = finalAngles.phiErrorSystematics[dataset];
            thetaError = finalAngles.thetaErrorSystematics[dataset];
        }

        cout << "Angle values for: " << boldOn << DatasetToString(dataset) << resetFormats << '\n';
        cout << greenOn;
        cout << boldOn << underlineOn << "ϕ:" << resetFormats << greenOn << " " << finalAngles.phi[dataset] << "\u00B0 ± "
             << phiError << "\u00B0.\n";
        cout << boldOn << underlineOn << "θ:" << resetFormats << greenOn << " " << finalAngles.theta[dataset] << "\u00B0 ± "
             << thetaError << "\u00B0.\n"
             << resetFormats;
        cout << "--------------------------------------------\n";
    }

    cout << "Angle values for: " << boldOn << "True Neutrino Direction" << resetFormats << '\n';
    cout << greenOn;
    cout << boldOn << underlineOn << "ϕ:" << resetFormats << greenOn << " " << finalAngles.phiTrue << "\u00B0 ± "
         << finalAngles.phiTrueError << "\u00B0.\n";
    cout << boldOn << underlineOn << "θ:" << resetFormats << greenOn << " " << finalAngles.thetaTrue << "\u00B0 ± "
         << finalAngles.thetaTrueError << "\u00B0.\n" << resetFormats;
    cout << "--------------------------------------------\n";
}

void FillOutputFile(array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize> const& histogram,
                    AngleValues const& finalAngles,
                    CovarianceValues const& oneSigmaEllipse)
{
    // Set up our output file
    auto outputFile = std::make_unique<TFile>("Directionality.root", "recreate");

    outputFile->cd();

    TVector3* angleOutput;
    TVector2* ellipseOutput;
    string outputName;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        angleOutput
            = new TVector3(finalAngles.phi[dataset], finalAngles.phiError[dataset], finalAngles.phiErrorSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Phi";
        outputFile->WriteTObject(angleOutput, outputName.c_str());

        ellipseOutput = new TVector2(oneSigmaEllipse.phiError[dataset], oneSigmaEllipse.phiErrorSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Phi Ellipse";
        outputFile->WriteTObject(ellipseOutput, outputName.c_str());

        angleOutput = new TVector3(
            finalAngles.theta[dataset], finalAngles.thetaError[dataset], finalAngles.thetaErrorSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Theta";
        outputFile->WriteTObject(angleOutput, outputName.c_str());

        ellipseOutput = new TVector2(oneSigmaEllipse.thetaError[dataset], oneSigmaEllipse.thetaErrorSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Theta Ellipse";
        outputFile->WriteTObject(ellipseOutput, outputName.c_str());

        ellipseOutput = new TVector2(oneSigmaEllipse.tilt[dataset], oneSigmaEllipse.tiltSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Tilt Ellipse";
        outputFile->WriteTObject(ellipseOutput, outputName.c_str());
    }

    ellipseOutput = new TVector2(finalAngles.phiTrue, finalAngles.phiTrueError);
    outputName = "True Phi";
    outputFile->WriteTObject(ellipseOutput, outputName.c_str());

    ellipseOutput = new TVector2(finalAngles.thetaTrue, finalAngles.thetaTrueError);
    outputName = "True Theta";
    outputFile->WriteTObject(ellipseOutput, outputName.c_str());

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
                histogram[dataset][signalSet][direction]->Write();
            }
        }
    }

    cout << boldOn << cyanOn << "Filled output file: " << resetFormats << blueOn << boldOn << "Directionality.root!\n"
         << resetFormats;
    cout << "--------------------------------------------\n";

    outputFile->Close();
}

int main(int argc, char* argv[])
{
    // Ignore Warnings (mostly for time honestly)
    gErrorIgnoreLevel = kError;

    // Using command line arguments for verbosity control
    for (int i = 1; i < argc; i++)
    {
        if (string(argv[i]) == "-D")
            DETECTOR_VERBOSITY = 1;
        else if (string(argv[i]) == "-L")
            LIVETIME_VERBOSITY = 1;
        else if (string(argv[i]) == "-I")
            IBDCOUNT_VERBOSITY = 1;
        else if (string(argv[i]) == "-M")
            MEAN_VERBOSITY = 1;
        else if (string(argv[i]) == "-S")
            ANGLES_STATISTICS = 1;
        else if (string(argv[i]) == "-C")
            COVARIANCE_VERBOSITY = 1;
    }

    // Take ownership of histograms
    TH1::AddDirectory(kFALSE);

    // Fill detector configuration
    FillDetectorConfig();

    // Set up what we're measuring
    IBDValues neutrinoCounts;
    AngleValues finalAngles;
    CovarianceValues oneSigmaEllipse;

    // Need histograms for counting each variable
    array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize> histogram;

    // Set up histograms for all 3 directions
    for (int dataset = Data; dataset < DatasetSize; dataset++)  // Dataset
    {
        for (int signalSet = CorrelatedReactorOn; signalSet < SignalSize; signalSet++)  // Signal set
        {
            // No reactor off for simulations
            if ((dataset == Sim || dataset == SimUnbiased)
                && (signalSet == CorrelatedReactorOff || signalSet == AccidentalReactorOff))
                continue;

            for (int direction = X; direction < DirectionSize; direction++)
            {
                string data = DatasetToString(dataset);
                string signal = SignalToString(signalSet);
                string axis = AxisToString(direction);
                string histogramName = data + " " + signal + " " + axis;
                histogram[dataset][signalSet][direction]
                    = std::make_shared<TH1F>(histogramName.c_str(), data.c_str(), bins, -histogramMax, histogramMax);
            }
        }
    }

    // Filling data histograms
    cout << "--------------------------------------------\n";
    cout << "Filling data histograms!\n";
    cout << "--------------------------------------------\n";
    for (int period = 1; period < 6; period++)  // 5 periods of PROSPECT Data
    {
        SetUpHistograms(histogram, Data, period);
    }
    lineCounter = 0;

    cout << boldOn << cyanOn << "Successfully filled data histogram!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    if (LIVETIME_VERBOSITY)
    {
        cout << "Total livetime for all" << boldOn << " Reactor Off " << resetFormats << "events: " << livetimeOff << '\n';
        cout << "Total livetime for all" << boldOn << " Reactor On " << resetFormats << "events: " << livetimeOn << '\n';
        cout << "--------------------------------------------\n";
    }

    // Filling simulation histograms
    cout << "Filling simulation histograms!\n" << resetFormats;
    cout << "--------------------------------------------\n";
    for (int period = 1; period < 6; period++)
    {
        SetUpHistograms(histogram, Sim, period);
    }

    cout << boldOn << cyanOn << "Successfully filled simulation histogram!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    neutrinoCounts = SubtractBackgrounds(histogram);
    AddSystematics(neutrinoCounts);
    finalAngles = CalculateAngles(neutrinoCounts);
    oneSigmaEllipse = CalculateCovariances(neutrinoCounts, finalAngles);
    OffsetTheta(finalAngles);
    PrintAngles(finalAngles);
    FillOutputFile(histogram, finalAngles, oneSigmaEllipse);

    return 0;
}