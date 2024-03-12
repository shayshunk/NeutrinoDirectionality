#include "NeutrinoDirectionality.h"

#include "DetectorConfig.h"
#include "Formatting.h"
#include "Timer.h"

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

bool CheckNeighbor(int periodNo, int segment, char direction)
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

    ResetLineNumber();
}

void Directionality::ReadFileList(int dataset, int periodNo)
{
    char const* path;

    if (dataset == Data)
        path = dataPath;

    else if (dataset == Sim)
        path = simPath;

    // Combining names into file list name
    string fileList = Form(path, std::to_string(periodNo).c_str(), std::to_string(periodNo).c_str());

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
    if (dataset == Data)
        fileName = dataFileName;

    else if (dataset == Sim)
        fileName = simFileName;

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

        // Combining names into root file name
        TString rootFilename = Form(fileName, std::to_string(period).c_str(), files[index].data());

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

                if (direction != Z)  // Cubical distance cuts in Z for X and Y
                {
                    float zPromptPosition = rootTree->GetLeaf("xyz")->GetValue(Z);
                    float zDelayedPosition = rootTree->GetLeaf("n_xyz")->GetValue(Z);

                    float displacement = zPromptPosition - zDelayedPosition;

                    if (displacement > 60)
                        continue;
                }

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

void Directionality::FillHistogramUnbiased(int signalSet)
{
    bool posDirection = false, negDirection = false, success = false;

    // Need to weight accidental datasets by deadtime correction factor
    double weight = (signalSet == AccidentalReactorOff || signalSet == AccidentalReactorOn) ? xRx : 1;

    // Check for live neighbors in different directions based on which axis
    // we're filling
    if (direction == X)
    {
        posDirection = CheckNeighbor(period, promptSegment, 'r');
        negDirection = CheckNeighbor(period, promptSegment, 'l');
    }
    else if (direction == Y)
    {
        posDirection = CheckNeighbor(period, promptSegment, 'u');
        negDirection = CheckNeighbor(period, promptSegment, 'd');
    }
    else  // Fill Z with same segment events
    {
        double displacement = delayedPosition - promptPosition;
        histogram[dataSet][signalSet][direction].Fill(displacement, weight);
    }

    // Dataset + 1 returns the unbiased version of that dataset
    if (posDirection && !negDirection)
        histogram[dataSet + 1][signalSet][direction].Fill(segmentWidth, weight);
    else if (!posDirection && negDirection)
        histogram[dataSet + 1][signalSet][direction].Fill(-segmentWidth, weight);
    else if (posDirection && negDirection)
        histogram[dataSet + 1][signalSet][direction].Fill(0.0, weight);
}

void Directionality::FillHistogram()
{
    // Applying energy cut
    if (Esmear < 0.8 || Esmear > 7.4)
    {
        return;
    }

    if (nCaptTime > pow(10, 3) && nCaptTime < 120 * pow(10, 3))  // Correlated Dataset
    {
        // Calculate neutron displacement
        double displacement = delayedPosition - promptPosition;

        // Figure out whether the reactor is on and assign signal index
        int signalSet = reactorOn ? CorrelatedReactorOn : CorrelatedReactorOff;

        // Fill regular dataset with displacement but Z only takes same segment
        if (direction != Z)
            histogram[dataSet][signalSet][direction].Fill(displacement);

        // Fill dead segment correction dataset
        if (promptSegment == delayedSegment)
        {
            FillHistogramUnbiased(signalSet);
        }
    }
    else if (nCaptTime > pow(10, 6))  // Accidental Dataset
    {
        // Calculate neutron displacement
        double displacement = delayedPosition - promptPosition;

        // Figure out whether the reactor is on and assign signal index
        int signalSet = reactorOn ? AccidentalReactorOn : AccidentalReactorOff;

        // Fill regular dataset with displacement but Z only takes same segment
        if (direction != Z)
            histogram[dataSet][signalSet][direction].Fill(displacement, xRx);

        // Fill dead segment correction dataset
        if (promptSegment == delayedSegment)
        {
            FillHistogramUnbiased(signalSet);
        }
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

            if (dataset == Data || dataset == Sim)
            {
                mean[dataset][direction] = histogram[dataset][TotalDifference][direction].GetMean();
                sigma[dataset][direction] = histogram[dataset][TotalDifference][direction].GetStdDev() / sqrt(effIBDs);
                effectiveIBD[dataset][direction] = effIBDs;
            }
            else
            {
                effectiveIBD[dataset][direction] = effectiveIBD[dataset - 1][direction];
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

    CalculateUnbiasing();
}

void Directionality::CalculateUnbiasing()
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
            nPlus = histogram[dataset - 1][TotalDifference][direction].GetBinContent(296);
            nPlusPlus = histogram[dataset][TotalDifference][direction].GetBinContent(297);
            nMinus = histogram[dataset - 1][TotalDifference][direction].GetBinContent(6);
            nMinusMinus = histogram[dataset][TotalDifference][direction].GetBinContent(5);
            nPlusMinus = histogram[dataset][TotalDifference][direction].GetBinContent(151);

            nPlusError = histogram[dataset - 1][TotalDifference][direction].GetBinError(296);
            nPlusPlusError = histogram[dataset][TotalDifference][direction].GetBinError(297);
            nMinusError = histogram[dataset - 1][TotalDifference][direction].GetBinError(6);
            nMinusMinusError = histogram[dataset][TotalDifference][direction].GetBinError(5);
            nPlusMinusError = histogram[dataset][TotalDifference][direction].GetBinError(151);

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

            mean[dataset][direction] = p;
            sigma[dataset][direction] = pError;
        }
        effectiveIBD[dataset][Z] = effectiveIBD[dataset - 1][Z];
        mean[dataset][Z] = mean[dataset - 1][Z];
        sigma[dataset][Z] = sigma[dataset - 1][Z];
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
                cout << boldOn << "p" << AxisToString(direction) << ": " << resetFormats << mean[dataset][direction] << " ± "
                     << sigma[dataset][direction] << '\n';
            }
            cout << "--------------------------------------------\n";
        }
    }
}

void Directionality::AddSystematics()
{
    // Systematics extracted from BiPo study

    // Defining variables for readability of code
    double sigmaX, sigmaY, sigmaZ;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        sigmaX = sigma[dataset][X];
        sigmaY = sigma[dataset][Y];
        sigmaZ = sigma[dataset][Z];

        sigmaSystematics[dataset][X] = sqrt(pow(sigmaX, 2) + pow(0.25, 2) + pow(0.08, 2));
        sigmaSystematics[dataset][Y] = sqrt(pow(sigmaY, 2) + pow(0.39, 2) + pow(0.08, 2));
        sigmaSystematics[dataset][Z] = sqrt(pow(sigmaZ, 2) + pow(0.06, 2) + pow(0.06, 2));
    }

    cout << boldOn << cyanOn << "Added Systematics.\n" << resetFormats;
    cout << "--------------------------------------------\n";

    // Printing out values
    if (SYSTEMATIC_MEAN_VERBOSITY)
    {
        for (int dataset = Data; dataset < DatasetSize; dataset++)
        {
            cout << "Mean and sigma values for: " << boldOn << DatasetToString(dataset) << resetFormats << '\n';
            for (int direction = X; direction < DirectionSize; direction++)
            {
                cout << boldOn << "p" << AxisToString(direction) << ": " << resetFormats << mean[dataset][direction] << " ± "
                     << sigmaSystematics[dataset][direction] << '\n';
            }
            cout << "--------------------------------------------\n";
        }
    }
}

void Directionality::CalculateAngles()
{
    // Defining variables for readability of code
    double px, py, pz;
    double sigmaX, sigmaY, sigmaZ;
    double sigmaXSystematics, sigmaYSystematics, sigmaZSystematics;
    double effIBDX, effIBDY, effIBDZ;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        // Grabbing values from current dataset
        px = mean[dataset][X];
        py = mean[dataset][Y];
        pz = mean[dataset][Z];

        sigmaX = sigma[dataset][X];
        sigmaY = sigma[dataset][Y];
        sigmaZ = sigma[dataset][Z];

        sigmaXSystematics = sigmaSystematics[dataset][X];
        sigmaYSystematics = sigmaSystematics[dataset][Y];
        sigmaZSystematics = sigmaSystematics[dataset][Z];

        effIBDX = effectiveIBD[dataset][X];
        effIBDY = effectiveIBD[dataset][Y];
        effIBDZ = effectiveIBD[dataset][Z];

        // phi = arctan(y / x)
        double tanPhi = py / px;
        double phiTemp = atan(tanPhi) * 180.0 / pi;
        double tanPhiError = sqrt(pow((sigmaX * py) / pow(px, 2), 2) + pow(sigmaY / px, 2));
        double phiErrorTemp = (tanPhiError / (1 + pow(tanPhi, 2))) * 180.0 / pi;
        double tanPhiErrorSystematics = sqrt(pow((sigmaXSystematics * py) / pow(px, 2), 2) + pow(sigmaYSystematics / px, 2));
        double phiErrorSystematicsTemp = (tanPhiErrorSystematics / (1 + pow(tanPhi, 2))) * 180.0 / pi;

        // theta = arctan(z / sqrt(x^2 + y^2))
        double tanTheta = pz / sqrt(pow(px, 2) + pow(py, 2));
        double thetaTemp = atan(tanTheta) * 180.0 / pi;
        double tanThetaError = sqrt((1 / (px * px + py * py))
                                    * (pow((px * pz * sigmaX / (px * px + py * py)), 2)
                                       + pow((py * pz * sigmaY / (px * px + py * py)), 2) + pow(sigmaZ, 2)));
        double thetaErrorTemp = (tanThetaError / (1 + pow(tanTheta, 2))) * 180.0 / pi;
        double tanThetaErrorSystematics
            = sqrt((1 / (px * px + py * py))
                   * (pow((px * pz * sigmaXSystematics / (px * px + py * py)), 2)
                      + pow((py * pz * sigmaYSystematics / (px * px + py * py)), 2) + pow(sigmaZSystematics, 2)));
        double thetaErrorSystematicsTemp = (tanThetaErrorSystematics / (1 + pow(tanTheta, 2))) * 180.0 / pi;

        // Storing values
        phi[dataset] = phiTemp;
        phiError[dataset] = phiErrorTemp;
        phiErrorSystematics[dataset] = phiErrorSystematicsTemp;
        theta[dataset] = thetaTemp;
        thetaError[dataset] = thetaErrorTemp;
        thetaErrorSystematics[dataset] = thetaErrorSystematicsTemp;
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
    phiTrue = atan(tanPhiTrue) * 180.0 / pi;
    float tanPhiTrueError = sqrt(pow((yTrue * xTrueError) / (xTrue * xTrue), 2) + pow(yTrueError / xTrue, 2));
    phiTrueError = tanPhiTrueError / (1 + pow(tanPhiTrue, 2)) * 180.0 / pi;

    float tanThetaTrue = zTrue / sqrt(pow(xTrue, 2) + pow(yTrue, 2));
    thetaTrue = atan(tanThetaTrue) * 180.0 / pi;
    float tanThetaTrueError
        = sqrt(pow(1 / sqrt(xTrue * xTrue + yTrue * yTrue), 2)
               * (pow(xTrue * zTrue / (xTrue * xTrue + yTrue * yTrue) * xTrueError, 2)
                  + pow(yTrue * zTrue / (xTrue * xTrue + yTrue * yTrue) * yTrueError, 2) + pow(zTrueError, 2)));
    thetaTrueError = tanThetaTrueError / (1 + pow(tanPhiTrue, 2)) * 180.0 / pi;
}

void Directionality::CalculateCovariances()
{
    // Calculating covariances
    array<array<float, 2>, 2> covarianceMatrix;
    array<array<float, 2>, 2> covarianceMatrixSystematics;

    // Defining variables for readability
    float px, py, pz;
    float sigmaX, sigmaY, sigmaZ;
    float sigmaXSystematics, sigmaYSystematics, sigmaZSystematics;
    float phiTemp, thetaTemp;

    // From error propagation document
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        // Grabbing values of current dataset from struct
        px = mean[dataset][X];
        py = mean[dataset][Y];
        pz = mean[dataset][Z];

        sigmaX = sigma[dataset][X];
        sigmaY = sigma[dataset][Y];
        sigmaZ = sigma[dataset][Z];

        sigmaXSystematics = sigmaSystematics[dataset][X];
        sigmaYSystematics = sigmaSystematics[dataset][Y];
        sigmaZSystematics = sigmaSystematics[dataset][Z];

        phiTemp = phi[dataset];
        thetaTemp = theta[dataset];

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
        covarianceMatrix[0][0] = covarianceMatrix[0][0] * pow(cos(phiTemp * pi / 180), 4);
        covarianceMatrix[0][1] = covarianceMatrix[0][1] * pow(cos(phiTemp * pi / 180), 2) * pow(cos(thetaTemp * pi / 180), 2);
        covarianceMatrix[1][0] = covarianceMatrix[1][0] * pow(cos(phiTemp * pi / 180), 2) * pow(cos(thetaTemp * pi / 180), 2);
        covarianceMatrix[1][1] = covarianceMatrix[1][1] * pow(cos(thetaTemp * pi / 180), 4);

        covarianceMatrixSystematics[0][0] = covarianceMatrixSystematics[0][0] * pow(cos(phiTemp * pi / 180), 4);
        covarianceMatrixSystematics[0][1] = covarianceMatrixSystematics[0][1] * pow(cos(phiTemp * pi / 180), 2)
                                            * pow(cos(thetaTemp * pi / 180), 2);
        covarianceMatrixSystematics[1][0] = covarianceMatrixSystematics[1][0] * pow(cos(phiTemp * pi / 180), 2)
                                            * pow(cos(thetaTemp * pi / 180), 2);
        covarianceMatrixSystematics[1][1] = covarianceMatrixSystematics[1][1] * pow(cos(thetaTemp * pi / 180), 4);

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

        float tiltTemp = atan(vector1_1) * 180.0 / pi;
        tiltTemp += 90;

        float vector1_1Systematics = lambda1Sytematics - dSystematics;
        float vector1_2Systematics = cSystematics;

        float normalizerSystematics = 1.0 / cSystematics;
        vector1_1Systematics *= normalizerSystematics;
        vector1_2Systematics *= normalizerSystematics;

        float tiltSystematicsTemp = atan(vector1_1Systematics) * 180.0 / pi;
        tiltSystematicsTemp += 90;

        // Calculating final angle errors
        float phiErrorTemp = sqrt(2.291 * lambda1) * 180.0 / pi;
        float thetaErrorTemp = sqrt(2.291 * lambda2) * 180.0 / pi;

        float phiErrorSystematicsTemp = sqrt(2.291 * lambda1Sytematics) * 180.0 / pi;
        float thetaErrorSystematicsTemp = sqrt(2.291 * lambda2Sytematics) * 180.0 / pi;

        // Filling array
        phiEllipseError[dataset] = phiErrorTemp;
        phiEllipseErrorSystematics[dataset] = phiErrorSystematicsTemp;
        thetaEllipseError[dataset] = thetaErrorTemp;
        thetaEllipseErrorSystematics[dataset] = thetaErrorSystematicsTemp;
        tilt[dataset] = tiltTemp;
        tiltSystematics[dataset] = tiltSystematicsTemp;
    }

    // Prints out the 1 sigma values if COVARIANCE_VERBOSITY is set to 1
    // Change at top of file
    if (COVARIANCE_VERBOSITY)
    {
        for (int dataset = Data; dataset < DatasetSize; dataset++)
        {
            cout << "The 1 sigma ellipse for: " << boldOn << DatasetToString(dataset) << resetFormats << " with systematics.\n";
            cout << greenOn;
            cout << boldOn << underlineOn << "Phi:" << resetFormats << greenOn << " " << phi[dataset] << "\u00B0 ± "
                 << phiEllipseErrorSystematics[dataset] << "\u00B0.\n";
            cout << boldOn << underlineOn << "Theta:" << resetFormats << greenOn << " " << theta[dataset] << "\u00B0 ± "
                 << thetaEllipseErrorSystematics[dataset] << "\u00B0.\n";
            cout << boldOn << underlineOn << "Tilt:" << resetFormats << greenOn << " " << tiltSystematics[dataset] << "\u00B0.\n"
                 << resetFormats;
            cout << "--------------------------------------------\n";

            cout << "The 1 sigma ellipse for: " << boldOn << DatasetToString(dataset) << resetFormats
                 << " without systematics.\n";
            cout << greenOn;
            cout << boldOn << underlineOn << "Phi:" << resetFormats << greenOn << " " << phi[dataset] << "\u00B0 ± "
                 << phiEllipseError[dataset] << "\u00B0.\n";
            cout << boldOn << underlineOn << "Theta:" << resetFormats << greenOn << " " << theta[dataset] << "\u00B0 ± "
                 << thetaEllipseError[dataset] << "\u00B0.\n";
            cout << boldOn << underlineOn << "Tilt:" << resetFormats << greenOn << " " << tilt[dataset] << "\u00B0.\n"
                 << resetFormats;
            cout << "--------------------------------------------\n";
        }
    }

    // Final cone of uncertainty
    float a = phiEllipseError[DataUnbiased];
    float b = thetaEllipseError[DataUnbiased];
    thetaTemp = 90 - theta[DataUnbiased];
    float solidAngle = pi * a * b * sin(thetaTemp * pi / 180.0);
    float solidAngleRadians = solidAngle * pow((pi / 180.0), 2);
    float coneAngle = acos(1 - (solidAngleRadians / (2 * pi))) * 180.0 / pi;

    cout << boldOn << greenOn << underlineOn << "Cone of Uncertainty:" << resetFormats << greenOn << " " << coneAngle << "\u00B0"
         << '\n'
         << resetFormats;
    cout << "--------------------------------------------\n";
}

void Directionality::OffsetTheta()
{
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        theta[dataset] = 90 - theta[dataset];
    }

    thetaTrue = 90 + thetaTrue;
}

void Directionality::PrintAngles()
{
    cout << boldOn << cyanOn << "Final Angles!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        float phiErrorTemp, thetaErrorTemp;

        if (dataset == Sim || dataset == SimUnbiased || ANGLES_STATISTICS)
        {
            phiErrorTemp = phiError[dataset];
            thetaErrorTemp = thetaError[dataset];
        }
        else
        {
            phiErrorTemp = phiErrorSystematics[dataset];
            thetaErrorTemp = thetaErrorSystematics[dataset];
        }

        cout << "Angle values for: " << boldOn << DatasetToString(dataset) << resetFormats << '\n';
        cout << greenOn;
        cout << boldOn << underlineOn << "ϕ:" << resetFormats << greenOn << " " << phi[dataset] << "\u00B0 ± " << phiErrorTemp
             << "\u00B0.\n";
        cout << boldOn << underlineOn << "θ:" << resetFormats << greenOn << " " << theta[dataset] << "\u00B0 ± "
             << thetaErrorTemp << "\u00B0.\n"
             << resetFormats;
        cout << "--------------------------------------------\n";
    }

    cout << "Angle values for: " << boldOn << "True Neutrino Direction" << resetFormats << '\n';
    cout << greenOn;
    cout << boldOn << underlineOn << "ϕ:" << resetFormats << greenOn << " " << phiTrue << "\u00B0 ± " << phiTrueError
         << "\u00B0.\n";
    cout << boldOn << underlineOn << "θ:" << resetFormats << greenOn << " " << thetaTrue << "\u00B0 ± " << thetaTrueError
         << "\u00B0.\n"
         << resetFormats;
    cout << "--------------------------------------------\n";

    float phiOffset = fabs(phi[Data] - phi[DataUnbiased]);
    float phiOffsetPercent = phiOffset / phi[DataUnbiased] * 100;
    float thetaOffset = fabs(theta[Data] - theta[DataUnbiased]);
    float thetaOffsetPercent = thetaOffset / theta[DataUnbiased] * 100;

    cout << "Offsets between " << boldOn << "Data" << resetFormats << " and " << boldOn << "Data Unbiased:\n";
    cout << resetFormats << greenOn << "ϕ:" << phiOffset << " or " << phiOffsetPercent << "%.\n";
    cout << "θ:" << thetaOffset << " or " << thetaOffsetPercent << "%.\n" << resetFormats;
    cout << "--------------------------------------------\n";
    
    phiOffset = fabs(phi[Sim] - phi[SimUnbiased]);
    phiOffsetPercent = phiOffset / phi[SimUnbiased] * 100;
    thetaOffset = fabs(theta[Sim] - theta[SimUnbiased]);
    thetaOffsetPercent = thetaOffset / theta[SimUnbiased] * 100;
    cout << "Offset between " << boldOn << "Sim" << resetFormats << " and " << boldOn << "Sim Unbiased:\n";
    cout << resetFormats << greenOn << "ϕ: " << phiOffset << "\u00B0 or " << phiOffsetPercent << "%.\n";
    cout << "θ: " << thetaOffset << "\u00B0 or " << thetaOffsetPercent << "%.\n" << resetFormats;
    cout << "--------------------------------------------\n";
}

void Directionality::FillOutputFile()
{
    // Set up our output file
    TFile outputFile("Directionality.root", "recreate");

    outputFile.cd();

    TVector3* angleOutput;
    TVector2* ellipseOutput;
    string outputName;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        angleOutput = new TVector3(phi[dataset], phiError[dataset], phiErrorSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Phi";
        outputFile.WriteTObject(angleOutput, outputName.c_str());

        ellipseOutput = new TVector2(phiError[dataset], phiErrorSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Phi Ellipse";
        outputFile.WriteTObject(ellipseOutput, outputName.c_str());

        angleOutput = new TVector3(theta[dataset], thetaError[dataset], thetaErrorSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Theta";
        outputFile.WriteTObject(angleOutput, outputName.c_str());

        ellipseOutput = new TVector2(thetaError[dataset], thetaErrorSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Theta Ellipse";
        outputFile.WriteTObject(ellipseOutput, outputName.c_str());

        ellipseOutput = new TVector2(tilt[dataset], tiltSystematics[dataset]);
        outputName = DatasetToString(dataset) + " Tilt Ellipse";
        outputFile.WriteTObject(ellipseOutput, outputName.c_str());
    }

    ellipseOutput = new TVector2(phiTrue, phiTrueError);
    outputName = "True Phi";
    outputFile.WriteTObject(ellipseOutput, outputName.c_str());

    ellipseOutput = new TVector2(thetaTrue, thetaTrueError);
    outputName = "True Theta";
    outputFile.WriteTObject(ellipseOutput, outputName.c_str());

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

    cout << boldOn << cyanOn << "Filled output file: " << resetFormats << blueOn << boldOn << "Directionality.root!\n"
         << resetFormats;
    cout << "--------------------------------------------\n";

    outputFile.Close();
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
        else if (string(argv[i]) == "-SM")
            SYSTEMATIC_MEAN_VERBOSITY = 1;
        else if (string(argv[i]) == "-S")
            ANGLES_STATISTICS = 1;
        else if (string(argv[i]) == "-C")
            COVARIANCE_VERBOSITY = 1;
    }

    // Take ownership of histograms
    TH1::AddDirectory(kFALSE);

    // Setting up timer
    Timer timer;

    // Setting up Directionality class
    Directionality neutrinoDirection;

    // Fill detector configuration
    FillDetectorConfig();

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

    neutrinoDirection.SubtractBackgrounds();
    neutrinoDirection.AddSystematics();
    neutrinoDirection.CalculateAngles();
    neutrinoDirection.CalculateCovariances();
    neutrinoDirection.OffsetTheta();
    neutrinoDirection.PrintAngles();
    neutrinoDirection.FillOutputFile();

    return 0;
}