#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"

void SetBasicPlotStyle()
{
    // Canvas Style

    // Turning off statistics box
    gStyle->SetOptStat(0);
    // Turning off automatic histogram titles
    gStyle->SetOptTitle(0);

    // Padding and font sizes
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetTitleFont(22, "t");
    gStyle->SetTitleSize(0.06, "t");
}

void SetPlotDirectory(std::string plotDirectory)
{
    if (0 != system(Form("test -d %s",plotDirectory.c_str())))
        system(Form("mkdir -m 755 -p %s", plotDirectory.c_str()));
}

void SetTitleStyles(TH1* histogram)
{
    // X axis
    histogram->GetXaxis()->CenterTitle(true);
    histogram->GetXaxis()->SetTitleFont(22);
    histogram->GetXaxis()->SetTitleSize(0.05);
    histogram->GetXaxis()->SetLabelFont(82);
    histogram->GetXaxis()->SetLabelSize(0.04);

    // Y axis
    histogram->GetYaxis()->CenterTitle(true);
    histogram->GetYaxis()->SetTitleFont(22);
    histogram->GetYaxis()->SetTitleSize(0.05);
    histogram->GetYaxis()->SetLabelFont(82);
    histogram->GetYaxis()->SetLabelSize(0.04);
    histogram->GetYaxis()->SetMaxDigits(4);

    // Z axis
    histogram->GetZaxis()->CenterTitle(true);
    histogram->GetZaxis()->SetTitleFont(22);
    histogram->GetZaxis()->SetTitleSize(0.05);
    histogram->GetZaxis()->SetLabelFont(82);
    histogram->GetZaxis()->SetLabelSize(0.04);
}