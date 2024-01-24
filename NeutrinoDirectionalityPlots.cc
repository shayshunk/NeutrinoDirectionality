#include "Formatting.h"
#include "NeutrinoDirectionality.h"
#include "Plotter.h"

using std::string, std::cout, std::array;

void MakeDisplacementPlots(array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize> const& histogram)
{
    // Taking ownership of histograms
    TH1::AddDirectory(kFALSE);

    // Setting up some basic padding and removing statistics boxes
    SetBasicPlotStyle();

    // Where we want our plots saved
    string plotDirectory = "DirectionalityPlots";

    // Checking if the directory already exists, making it if not
    SetPlotDirectory(plotDirectory);

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        for (int direction = X; direction < Z; direction++)
        {
            string plotName = DatasetToString(dataset) + " Displacement " + AxisToString(direction);
            string histName = DatasetToString(dataset);
            TCanvas canvas(plotName.c_str(), plotName.c_str(), 2000, 1700);

            // Temporary histogram for x and y because we want a 3 bin plot
            TH1F tempHist(histName.c_str(), histName.c_str(), 3, 0, 3);

            // Setting up where to get the values from
            int lowBin, highBin, mediumBin = 151;

            // These come from the main analysis code
            if (dataset == Data || dataset == Sim)
            {
                lowBin = 6;
                highBin = 296;
            }
            else
            {
                lowBin = 5;
                highBin = 297;
            }

            // Grabbing values from the main histograms
            tempHist.SetBinContent(1, histogram[dataset][TotalDifference][direction]->GetBinContent(lowBin));
            tempHist.SetBinError(1, histogram[dataset][TotalDifference][direction]->GetBinError(lowBin));

            tempHist.SetBinContent(2, histogram[dataset][TotalDifference][direction]->GetBinContent(mediumBin));
            tempHist.SetBinError(2, histogram[dataset][TotalDifference][direction]->GetBinError(mediumBin));

            tempHist.SetBinContent(3, histogram[dataset][TotalDifference][direction]->GetBinContent(highBin));
            tempHist.SetBinError(3, histogram[dataset][TotalDifference][direction]->GetBinError(highBin));

            // Setting up titles and draw styles
            string histTitle = "Delayed - Prompt Positions: ";
            histTitle = histTitle + AxisToString(direction) + " (mm)";
            tempHist.GetXaxis()->SetTitle(histTitle.c_str());
            tempHist.GetYaxis()->SetTitle("IBD Events");

            tempHist.SetLineColor(kRed);
            tempHist.SetLineWidth(4);
            tempHist.SetLineStyle(1);
            tempHist.SetMarkerStyle(20);
            tempHist.SetMarkerSize(1.3);
            tempHist.SetMarkerColor(kBlack);

            // Applying header function
            SetTitleStyles(&tempHist);

            // Custom labels and size
            // For some reason the label size that's normal on other plots is tiny for the X labels here
            tempHist.GetXaxis()->SetBinLabel(1, "-D");
            tempHist.GetXaxis()->SetBinLabel(2, "0");
            tempHist.GetXaxis()->SetBinLabel(3, "D");
            tempHist.GetXaxis()->SetLabelSize(0.07);

            tempHist.Draw("HISTS ][");

            string fullPath = plotDirectory + "/" + plotName + ".png";
            canvas.SaveAs(fullPath.c_str());
        }

        // Different plot style for Z
        string plotName = DatasetToString(dataset) + " Displacement Z";
        string histName = DatasetToString(dataset);
        TCanvas canvas(plotName.c_str(), plotName.c_str(), 2000, 1700);

        auto tempHist = std::unique_ptr<TH1F>(static_cast<TH1F*>(histogram[dataset][TotalDifference][Z]->Clone()));

        // Setting up titles and draw styles
        string histTitle = "Delayed - Prompt Positions: ";
        histTitle = histTitle + "Z (mm)";
        tempHist->GetXaxis()->SetTitle(histTitle.c_str());
        tempHist->GetYaxis()->SetTitle("IBD Events");

        tempHist->SetLineColor(kRed);
        tempHist->SetLineWidth(2);
        tempHist->SetLineStyle(1);
        tempHist->SetMarkerStyle(20);
        tempHist->SetMarkerSize(1);
        tempHist->SetMarkerColor(kRed);

        // Applying header function
        SetTitleStyles(tempHist.get());

        // Setting up fit
        TF1 gaussian("Fit", "gaus", -140, 140);
        gaussian.SetLineColor(kBlue);
        gaussian.SetLineWidth(2);
        gaussian.SetLineStyle(1);

        tempHist->Fit("Fit", "RQ");
        string zMean = Form("%.2f", gaussian.GetParameter(1));
        string zError = Form("%.2f", gaussian.GetParError(1));
        string zSigma = Form("%.2f", gaussian.GetParameter(2));
        float zChiSquare = gaussian.GetChisquare();
        float zNDF = gaussian.GetNDF();
        float zChiSquareNDF = zChiSquare / zNDF;
        string zChi2NDF = Form("%.2f", zChiSquareNDF);

        tempHist->Draw("P E1 X0");

        // Setting up legend
        float x_legend = 0.6;
        float y_legend = 0.7;

        TLegend legend(x_legend, y_legend, x_legend + 0.25, y_legend + 0.2);
        legend.SetFillColor(0);
        legend.SetFillStyle(0);
        legend.SetTextFont(82);
        legend.SetTextSize(0.03);

        //legend->SetHeader("Fit Parameters", "C");
        string meanText = "#mu: " + zMean + " #pm " + zError;
        legend.AddEntry(&gaussian, meanText.c_str(), "l");

        string sigmaText = "#sigma: " + zSigma;
        legend.AddEntry(&gaussian, sigmaText.c_str(), "l");

        string chiSquareText = "#frac{#chi^{2}}{ndf}: " + zChi2NDF;
        legend.AddEntry(&gaussian, chiSquareText.c_str(), "l");

        legend.Draw();

        string fullPath = plotDirectory + "/" + plotName + ".png";
        canvas.SaveAs(fullPath.c_str());
    }
}

int NeutrinoDirectionalityPlots()
{
    // Taking ownership of TH1 pointers
    TH1::AddDirectory(kFALSE);

    // Setting up what we need to plot
    array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize> histogram;
    float phiTrue, phiTrueError, thetaTrue, thetaTrueError;
    array<float, DatasetSize> phi;
    array<float, DatasetSize> phiError;
    array<float, DatasetSize> phiErrorSystematics;
    array<float, DatasetSize> theta;
    array<float, DatasetSize> thetaError;
    array<float, DatasetSize> thetaErrorSystematics;
    array<float, DatasetSize> tilt;
    array<float, DatasetSize> tiltSystematics;

    // Opening the root file
    auto rootFile = std::make_unique<TFile>("Directionality.root", "read");

    // Grabbing angles and errors
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        // Grabbing phi
        string vectorName = DatasetToString(dataset) + " Phi";
        string ellipseName = vectorName + " Ellipse";

        TVector3* input = (TVector3*)rootFile->Get(vectorName.c_str());
        phi[dataset] = input->X();
        TVector2* inputErrors = (TVector2*)rootFile->Get(ellipseName.c_str());
        phiError[dataset] = inputErrors->X();
        phiErrorSystematics[dataset] = inputErrors->Y();

        // Grabbing theta
        vectorName = DatasetToString(dataset) + " Theta";
        ellipseName = vectorName + " Ellipse";

        input = (TVector3*)rootFile->Get(vectorName.c_str());
        theta[dataset] = input->X();
        inputErrors = (TVector2*)rootFile->Get(ellipseName.c_str());
        thetaError[dataset] = inputErrors->X();
        thetaErrorSystematics[dataset] = inputErrors->Y();

        // Grabbing tilt
        vectorName = DatasetToString(dataset) + " Tilt Ellipse";

        input = (TVector3*)rootFile->Get(vectorName.c_str());
        tilt[dataset] = input->X();
        tiltSystematics[dataset] = input->Y();
    }

    // Grabbing true angles
    TVector2* input = (TVector2*)rootFile->Get("True Phi");
    phiTrue = input->X();
    phiTrueError = input->Y();

    input = (TVector2*)rootFile->Get("True Theta");
    thetaTrue = input->X();
    thetaTrueError = input->Y();

    // Grabbing histograms
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        for (int signalSet = CorrelatedReactorOn; signalSet < SignalSize; signalSet++)
        {
            for (int direction = X; direction < DirectionSize; direction++)
            {
                // Get histogram and static cast to shared pointer for safety
                string histogramName = DatasetToString(dataset) + " " + SignalToString(signalSet) + " " + AxisToString(direction);

                histogram[dataset][signalSet][direction]
                    = std::shared_ptr<TH1F>(static_cast<TH1F*>(rootFile->Get(histogramName.c_str())));
            }
        }
    }

    MakeDisplacementPlots(histogram);

    return 0;
}