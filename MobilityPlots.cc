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
            canvas.Update();

            // Temporary histogram for x and y because we want a 3 bin plot
            TH1F tempHist(histName.c_str(), histName.c_str(), 3, 0, 3);

            // Setting up where to get the values from
            int lowBin, highBin, mediumBin = 151;

            // These come from the main analysis code
            lowBin = 5;
            highBin = 297;

            // Grabbing values from the main histograms
            tempHist.SetBinContent(1, histogram[dataset][TotalDifference][direction]->GetBinContent(lowBin));
            tempHist.SetBinError(1, histogram[dataset][TotalDifference][direction]->GetBinError(lowBin));

            tempHist.SetBinContent(2, histogram[dataset][TotalDifference][direction]->GetBinContent(mediumBin));
            tempHist.SetBinError(2, histogram[dataset][TotalDifference][direction]->GetBinError(mediumBin));

            tempHist.SetBinContent(3, histogram[dataset][TotalDifference][direction]->GetBinContent(highBin));
            tempHist.SetBinError(3, histogram[dataset][TotalDifference][direction]->GetBinError(highBin));

            // Setting up titles and draw styles
            if (direction == X)
                tempHist.GetXaxis()->SetTitle("#bf{#it{s}}_{#it{x}}: Delayed - Prompt");
            else
                tempHist.GetXaxis()->SetTitle("#bf{#it{s}}_{#it{y}}: Delayed - Prompt");

            tempHist.GetYaxis()->SetTitle("IBD Events");

            tempHist.SetLineColor(kRed);
            tempHist.SetLineWidth(9);
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

            tempHist.GetXaxis()->SetTitleSize(0.062);
            tempHist.GetYaxis()->SetTitleSize(0.062);

            string fullPath = plotDirectory + "/" + plotName + ".png";
            canvas.SaveAs(fullPath.c_str());
        }

        // Different plot style for Z
        string plotName = DatasetToString(dataset) + " Displacement Z";
        string histName = DatasetToString(dataset);

        // Setting up canvas, no top margin because no title
        TCanvas canvas(plotName.c_str(), plotName.c_str(), 2000, 1700);
        canvas.SetCanvasSize(2000, 1700);
        canvas.SetTopMargin(0.075);

        auto tempHist = std::unique_ptr<TH1F>(static_cast<TH1F*>(histogram[dataset][TotalDifference][Z]->Clone()));

        // Setting up titles and draw styles
        tempHist->SetTitle("Z Displacement (0 - 20 #mus)");
        tempHist->GetXaxis()->SetTitle("#bf{#it{s}}_{#it{z}}: Delayed - Prompt (mm)");
        tempHist->GetYaxis()->SetTitle("IBD Events");

        tempHist->SetLineColor(kRed);
        tempHist->SetLineWidth(4);
        tempHist->SetLineStyle(1);
        tempHist->SetMarkerStyle(20);
        tempHist->SetMarkerSize(1);
        tempHist->SetMarkerColor(kRed);

        // Applying header function
        SetTitleStyles(tempHist.get());
        gStyle->SetOptTitle(1);

        // Getting total counts
        int counts = tempHist->GetEntries();

        // Setting up fit
        TF1 gaussian("Fit", "gaus", -140, 140);
        gaussian.SetLineColor(kBlue);
        gaussian.SetLineWidth(4);
        gaussian.SetLineStyle(1);

        tempHist->Fit("Fit", "RQ");
        string zMean = Form("%.2f", gaussian.GetParameter(1));
        string zError = Form("%.2f", gaussian.GetParError(1));
        string zSigma = Form("%.2f", gaussian.GetParameter(2));
        string zSigmaError = Form("%.2f", gaussian.GetParError(2));
        int zChiSquare = round(gaussian.GetChisquare());
        int zNDF = gaussian.GetNDF();
        string zChi2 = Form("%i", zChiSquare);
        string zNDFText = Form("%i", zNDF);
        string zChi2NDF = "#frac{" + zChi2 + "}{" + zNDFText + "}";

        tempHist->Draw("P E1 X0");

        tempHist->GetXaxis()->SetTitleSize(0.062);
        tempHist->GetYaxis()->SetTitleSize(0.062);
        // tempHist->GetYaxis()->SetRangeUser(0, 450);

        // Setting up legend
        float x_legend = 0.6;
        float y_legend = 0.7;

        TLegend legend(0.625, 0.75, 0.95, 0.925);
        legend.SetFillColor(0);
        legend.SetFillStyle(0);
        // legend.SetTextFont(82);
        legend.SetTextSize(0.04);

        // legend->SetHeader("Fit Parameters", "C");
        /* string meanText = "Gaussian fit, #it{#mu} = " + zMean + " #pm " + zError + " mm";
        legend.AddEntry(&gaussian, meanText.c_str(), "l"); */

        string sigmaText = "#sigma: " + zSigma + " mm";
        legend.AddEntry(&gaussian, sigmaText.c_str(), "l");

        string countsText = "Entries: " + std::to_string(counts);
        legend.AddEntry(tempHist.get(), countsText.c_str(), "p");

        /* string chiSquareText = "#frac{#chi^{2}}{ndf}:  " + zChi2NDF;
        legend.AddEntry(&gaussian, chiSquareText.c_str(), "l"); */

        legend.Draw();

        string fullPath = plotDirectory + "/" + plotName + ".png";
        canvas.SaveAs(fullPath.c_str());

        cout << "Z Mean for " << boldOn << DatasetToString(dataset) << " = " << zMean << " ± " << zError << resetFormats << '\n';
        cout << "Z Sigma for " << boldOn << DatasetToString(dataset) << " = " << zSigma << " ± " << zSigmaError << resetFormats
             << '\n';
    }
}

int MobilityPlots()
{
    // Taking ownership of TH1 pointers
    TH1::AddDirectory(kFALSE);

    // Don't open the canvases and then destroy them
    gROOT->SetBatch(kTRUE);

    // Setting up some basic padding and removing statistics boxes
    SetBasicPlotStyle();

    // Where we want our plots saved
    plotDirectory = "MobilityPlots";

    // Checking if the directory already exists, making it if not
    SetPlotDirectory(plotDirectory);

    // Setting up what we need to plot
    array<array<array<std::shared_ptr<TH1F>, DirectionSize>, SignalSize>, DatasetSize> histogram;

    // Opening the root file
    auto rootFile = std::make_unique<TFile>("Mobility_20.root", "read");

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

    return 0;
}