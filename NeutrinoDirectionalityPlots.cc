#include "DirectionalityPlots.h"
#include "Formatting.h"
#include "Plotter.h"

using std::string, std::cout, std::array;

void MakeDisplacementPlots(array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize> const& histogram)
{
    // Taking ownership of histograms
    TH1::AddDirectory(kFALSE);

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        for (int direction = X; direction < Z; direction++)
        {
            string plotName = DatasetToString(dataset) + " Displacement " + AxisToString(direction);
            string histName = DatasetToString(dataset);

            // Setting up canvas, no top margin because no title
            TCanvas canvas(plotName.c_str(), plotName.c_str());
            canvas.SetCanvasSize(2000, 1700);
            canvas.SetTopMargin(0.05);

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

        // Setting up canvas, no top margin because no title
        TCanvas canvas(plotName.c_str(), plotName.c_str(), 2000, 1700);
        canvas.SetCanvasSize(2000, 1700);
        canvas.SetTopMargin(0.05);

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

        TLegend legend(0.65, 0.7, 0.95, 0.95);
        legend.SetFillColor(0);
        legend.SetFillStyle(0);
        legend.SetTextFont(82);
        legend.SetTextSize(0.03);

        // legend->SetHeader("Fit Parameters", "C");
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

void MakeAnglePlots(FinalValues const& finalValues)
{
    // Setting up canvas, no top margin because no title
    TCanvas canvas("Angles", "Angles");

    canvas.SetCanvasSize(1600, 1550);
    canvas.SetTopMargin(0.05);
    canvas.SetBottomMargin(0.125);

    // Need to draw a frame to set bounds
    auto frame = canvas.DrawFrame(93.5, 33.5, 108.5, 48.5);

    // Setting up axis labels
    frame->GetXaxis()->SetTitle("#theta (deg)");
    frame->GetYaxis()->SetTitle("#phi (deg)");
    frame->GetXaxis()->CenterTitle(kTRUE);
    frame->GetYaxis()->CenterTitle(kTRUE);

    // Making sure the font is legible
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);

    // Error ellipses for each data point
    TEllipse Data(finalValues.theta[DataUnbiased],
                  finalValues.phi[DataUnbiased],
                  finalValues.thetaErrorSystematics[DataUnbiased],
                  finalValues.phiErrorSystematics[DataUnbiased],
                  0,
                  360,
                  finalValues.tiltSystematics[DataUnbiased]);
    Data.SetFillColor(kAzure - 3);
    Data.SetFillStyle(3001);
    Data.Draw();

    TEllipse DataStatistics(finalValues.theta[DataUnbiased],
                            finalValues.phi[DataUnbiased],
                            finalValues.thetaError[DataUnbiased],
                            finalValues.phiError[DataUnbiased],
                            0,
                            360,
                            finalValues.tilt[DataUnbiased]);
    DataStatistics.SetFillColor(kAzure + 3);
    DataStatistics.SetFillStyle(3001);
    DataStatistics.Draw();

    TEllipse Simulation(finalValues.theta[SimUnbiased],
                        finalValues.phi[SimUnbiased],
                        finalValues.thetaError[SimUnbiased],
                        finalValues.phiError[SimUnbiased],
                        0,
                        360,
                        finalValues.tilt[SimUnbiased]);
    Simulation.SetFillColor(kRed + 2);
    Simulation.SetFillStyle(3001);
    Simulation.Draw();

    TEllipse SimulationUncorrected(finalValues.theta[Sim],
                                   finalValues.phi[Sim],
                                   finalValues.thetaError[Sim],
                                   finalValues.phiError[Sim],
                                   0,
                                   360,
                                   finalValues.tilt[Sim]);
    SimulationUncorrected.SetFillColor(kMagenta + 3);
    SimulationUncorrected.SetFillStyle(3001);
    SimulationUncorrected.Draw();

    TEllipse True(finalValues.thetaTrue, finalValues.phiTrue, finalValues.thetaTrueError, finalValues.phiTrueError);
    True.SetFillColor(kGreen + 2);
    True.SetFillStyle(3001);
    True.Draw();

    // Setting up data markers
    TMarker dataPoint(finalValues.theta[DataUnbiased], finalValues.phi[DataUnbiased], 43);
    dataPoint.SetMarkerSize(5.5);
    dataPoint.SetMarkerColor(kAzure - 2.25);
    dataPoint.Draw();

    TMarker dataSystematicsPoint(finalValues.theta[DataUnbiased], finalValues.phi[DataUnbiased], 43);
    dataSystematicsPoint.SetMarkerSize(5.5);
    dataSystematicsPoint.SetMarkerColor(kAzure + 3);

    TMarker simPoint(finalValues.theta[SimUnbiased], finalValues.phi[SimUnbiased], 29);
    simPoint.SetMarkerSize(4.5);
    simPoint.SetMarkerColor(kRed + 2.25);
    simPoint.Draw();

    TMarker simUncorrectedPoint(finalValues.theta[Sim], finalValues.phi[Sim], 43);
    simUncorrectedPoint.SetMarkerSize(4.5);
    simUncorrectedPoint.SetMarkerColor(kMagenta + 3.25);
    simUncorrectedPoint.Draw();

    TMarker truePoint(finalValues.thetaTrue, finalValues.phiTrue, 33);
    truePoint.SetMarkerSize(5);
    truePoint.SetMarkerColor(kGreen + 2.25);
    truePoint.Draw();

    TLegend legend(0.54, 0.75, 0.95, 0.95);
    legend.SetTextFont(62);
    legend.SetTextSize(0.03);
    legend.AddEntry(&dataPoint, "Data", "p");
    legend.AddEntry(&dataSystematicsPoint, "Data - No Systematics", "p");
    legend.AddEntry(&simPoint, "Simulation", "p");
    legend.AddEntry(&simUncorrectedPoint, "Uncorrected Sim", "p");
    legend.AddEntry(&truePoint, "True Neutrino Direction", "p");
    legend.Draw();

    string fullPath = plotDirectory + "/FinalPlot.png";
    canvas.SaveAs(fullPath.c_str());
}

int NeutrinoDirectionalityPlots()
{
    // Taking ownership of TH1 pointers
    TH1::AddDirectory(kFALSE);

    // Don't open the canvases and then destroy them
    gROOT->SetBatch(kTRUE);

    // Setting up some basic padding and removing statistics boxes
    SetBasicPlotStyle();

    // Where we want our plots saved
    plotDirectory = "DirectionalityPlots";

    // Checking if the directory already exists, making it if not
    SetPlotDirectory(plotDirectory);

    // Setting up what we need to plot
    array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize> histogram;
    FinalValues finalValues;

    // Opening the root file
    auto rootFile = std::make_unique<TFile>("Directionality.root", "read");

    // Grabbing angles and errors
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        // Grabbing phi
        string vectorName = DatasetToString(dataset) + " Phi";
        string ellipseName = vectorName + " Ellipse";

        TVector3* input = (TVector3*)rootFile->Get(vectorName.c_str());
        finalValues.phi[dataset] = input->X();
        TVector2* inputErrors = (TVector2*)rootFile->Get(ellipseName.c_str());
        finalValues.phiError[dataset] = inputErrors->X();
        finalValues.phiErrorSystematics[dataset] = inputErrors->Y();

        // Grabbing theta
        vectorName = DatasetToString(dataset) + " Theta";
        ellipseName = vectorName + " Ellipse";

        input = (TVector3*)rootFile->Get(vectorName.c_str());
        finalValues.theta[dataset] = input->X();
        inputErrors = (TVector2*)rootFile->Get(ellipseName.c_str());
        finalValues.thetaError[dataset] = inputErrors->X();
        finalValues.thetaErrorSystematics[dataset] = inputErrors->Y();

        // Grabbing tilt
        vectorName = DatasetToString(dataset) + " Tilt Ellipse";

        input = (TVector3*)rootFile->Get(vectorName.c_str());
        finalValues.tilt[dataset] = input->X();
        finalValues.tiltSystematics[dataset] = input->Y();
    }

    // Grabbing true angles
    TVector2* input = (TVector2*)rootFile->Get("True Phi");
    finalValues.phiTrue = input->X();
    finalValues.phiTrueError = input->Y();

    input = (TVector2*)rootFile->Get("True Theta");
    finalValues.thetaTrue = input->X();
    finalValues.thetaTrueError = input->Y();

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

    rootFile->Close();

    MakeDisplacementPlots(histogram);
    MakeAnglePlots(finalValues);

    return 0;
}