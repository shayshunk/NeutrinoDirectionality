#include <time.h>
#include <vector>
#include <map>
#include <iostream>
#include "TH1D.h"
#include "TDatime.h"
#include "TVectorD.h"
#include "TChain.h"
#include "TSystem.h"
#include "TChainElement.h"
#include "TCollection.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "BP.C"
#include "PROSPECT_Style.cc"

using std::cout;

int BiPoxyz() 
{
    // Defining useful physical constants
    const double n2f = 1/12.0;
    const double tauBiPo = 0.1643/log(2);

    // Initilizaing Data Structure
    BP *bp = new BP();

    TChain *ch = bp->chain;

    // Setting up histograms and styles
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    TH1D *hp_x = new TH1D("hp_x", "Alpha-Beta X Position Difference", 2001, -1000, 1000);
    hp_x->SetLineColor(kBlue);
    hp_x->SetLineWidth(2);

    TH1D *hp_y = new TH1D("hp_y", "Alpha-Beta Y Position Difference", 2001, -1000, 1000);
    hp_y->SetLineColor(kRed);
    hp_y->SetLineWidth(2);

    TH1D *hp_z = new TH1D("hp_z", "Alpha-Beta Z Position Difference", 2001, -1000, 1000);
    hp_z->SetLineColor(kBlack);
    hp_z->SetLineWidth(2);
    hp_z->Sumw2();

    TH1D *hf_x = new TH1D("hf_x", "Alpha-Beta X Position Difference (Accidental)", 2001, -1000, 1000);
    hf_x->SetLineColor(kRed);
    hf_x->SetLineWidth(2);

    TH1D *hf_y = new TH1D("hf_y", "Alpha-Beta Y Position Difference (Accidental)", 2001, -1000, 1000);
    hf_y->SetLineColor(kBlack);
    hf_y->SetLineWidth(2);

    TH1D *hf_z = new TH1D("hf_z", "Alpha-Beta Z Position Difference (Accidental)", 2001, -1000, 1000);
    hf_z->SetLineColor(kBlue);
    hf_z->SetLineWidth(2);
    hf_z->Sumw2();

    TH1D *hx = new TH1D("hx", "Alpha-Beta X Position Difference (Acc Subtr)", 2001, -1000, 1000);
    hx->SetLineColor(kMagenta);
    hx->SetLineWidth(2);

    TH1D *hy = new TH1D("hy", "Alpha-Beta Y Position Difference (Acc Subtr)", 2001, -1000, 1000);
    hy->SetLineColor(kMagenta);
    hy->SetLineWidth(2);

    TH1D *hz = new TH1D("hz", "Alpha-Beta Z Position Difference (Acc Subtr)", 2001, -1000, 1000);
    hz->SetLineColor(kMagenta);
    hz->SetLineWidth(2);

    //------------------------------
    // Setting boundary cuts
    double hAE = 0.98, lAE = 0.73, hApsd = 0.34, lApsd = 0.17; // alpha
    double highBE = 4.0, lowBE = 0, hPpsd = 0.22, lPpsd = 0.05; // beta
    double t_start = 7.0e-4, t_end = 3 * tauBiPo;

    // Beginning data reading loop
    int noEntries = ch->GetEntries();

    for (int i = 0; i < noEntries; ++i)
    {
        bp->GetEntry(i);
        if (i % 500000 == 0)
        {
            cout << "Entry " << i << " of " << noEntries << "\n";
            double prog = (double(i))/double(noEntries);
            cout << prog*100 << "% complete.\n";
        }

        auto multPrompt = bp->mult_prompt;
        auto multFar = bp->mult_far;

        // Apply alpha cuts

        double alphaE = bp->aE;
        double alphaPSD = bp->aPSD;
        double alphaZ = bp->az;

        if (!(abs(alphaZ) < 1000))
        {
            continue;
        }
        if (alphaE < lAE || alphaE > hAE)
        {
            continue;
        }
        if (alphaPSD < lApsd || alphaPSD > hApsd)
        {
            continue;
        }

        // Filling Histograms

        for (int j = 0; j < multPrompt; ++j)
        {
            if (bp->pmult_clust->at(j) != bp->pmult_clust_ioni->at(j))
            {
                continue; // Throwing out clusters with recoils mixed in
            }
            if (bp->pEtot->at(j) < lowBE || bp->pEtot->at(j) > highBE)
            {
                continue; // Optional beta energy cut
            }
            if (multPrompt == 0)
            {
                continue;
            }

            double alphaT = bp->at;
            double betaT = bp->pt->at(j);
            double dt = alphaT - betaT;

            if (dt > t_start && dt < t_end)
            {
                int alphaX = bp->aseg % 14;
                int alphaY = bp->aseg / 14;
                int betaX = bp->pseg->at(j) % 14;
                int betaY = bp->pseg->at(j) / 14;
                double betaZ = bp->pseg->at(j);

                double dx = 145.7 * (alphaX - betaX);
                double dy = 145.7 * (alphaY - betaY);
                double dz = 145.7 * (alphaZ - betaZ);

                hp_x->Fill(dx);
                hp_y->Fill(dy);
                hp_z->Fill(dz);
            }
        }

        for (int j = 0; j < multFar; ++j)
        {
            if (bp->fmult_clust->at(j) != bp->fmult_clust_ioni->at(j))
            {
                continue; // Throwing out clusters with recoils mixed in 
            }
            if (bp->fEtot->at(j) < lowBE || bp->fEtot->at(j) > highBE)
            {
                continue; // Optional beta energy cut
            }
            if (multFar == 0)
            {
                continue;
            }

            double alphaT = bp->at;
            double betaT = bp->ft->at(j);
            double dt = (betaT - alphaT) - 10.0*tauBiPo;
            dt *= n2f;

            if (dt > 0 && dt < (t_end - t_start))
            {
                int alphaX = bp->aseg % 14;
                int alphaY = bp->aseg / 14;
                int betaX = bp->fseg->at(j) % 14;
                int betaY = bp->fseg->at(j) / 14;
                double betaZ = bp->fseg->at(j);

                double dx = 145.7 * (alphaX - betaX);
                double dy = 145.7 * (alphaY - betaY);
                double dz = 145.7 * (alphaZ - betaZ);

                hf_x->Fill(dx, n2f);
                hf_y->Fill(dy, n2f);
                hf_z->Fill(dz, n2f);
            }
        }
    }

    double x_legend = 0.6;
    double y_legend = 0.77;

    TCanvas *c = new TCanvas("c", "c", 0, 0, 1600, 1000);
    c->Divide(2, 1);

    c->cd(1);

    hp_x->Draw("hist");
    hf_x->Draw("hist same");
    hp_x->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hp_x->GetYaxis()->SetTitle("Counts");

    c->cd(2);

    hx = (TH1D*)hp_x->Clone("hx");
    hx->Add(hf_x, -1);
    hx->Draw("hist");
    hx->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hx->GetYaxis()->SetTitle("Counts");
    
    TLegend *leg = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.02);
    
    TLegendEntry* l11 = leg->AddEntry( (TObject*)0, Form("Mean = %.2f", hx->GetMean()), "" );
    l11->SetTextColor(kBlack);
    TLegendEntry* l12 = leg->AddEntry( (TObject*)0, Form("Std Dev = %.2f", hx->GetStdDev()), "" );
    l12->SetTextColor(kBlack);

    leg->Draw();

    c->SaveAs("X_AlphaBeta.png");        

    TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 1600, 1000);
    c1->Divide(2, 1);

    c1->cd(1);

    hp_y->Draw("hist");
    hf_y->Draw("hist same");
    hp_y->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hp_y->GetYaxis()->SetTitle("Counts");

    c1->cd(2);

    hy = (TH1D*)hp_y->Clone("hy");
    hy->Add(hf_y, -1);
    hy->Draw("hist");
    hy->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hy->GetYaxis()->SetTitle("Counts");
    
    TLegend *leg1 = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
    leg1->SetBorderSize(1);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextFont(62);
    leg1->SetTextSize(0.02);
    
    TLegendEntry* l111 = leg1->AddEntry( (TObject*)0, Form("Mean = %.2f", hy->GetMean()), "" );
    l111->SetTextColor(kBlack);
    TLegendEntry* l121 = leg1->AddEntry( (TObject*)0, Form("Std Dev = %.2f", hy->GetStdDev()), "" );
    l121->SetTextColor(kBlack);

    leg1->Draw();

    c1->SaveAs("Y_AlphaBeta.png");           

    TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 1600, 1000);
    c2->Divide(2, 1);

    c2->cd(1);

    hp_z->Draw("hist");
    hf_z->Draw("hist same");
    hp_z->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hp_z->GetYaxis()->SetTitle("Counts");

    c2->cd(2);

    hz = (TH1D*)hp_z->Clone("hz");
    hz->Add(hf_z, -1);
    hz->Draw("hist");
    hz->GetXaxis()->SetTitle("Beta/Alpha Displacement (mm)");
    hz->GetYaxis()->SetTitle("Counts");

    TLegend *leg2 = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
    leg2->SetBorderSize(1);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(62);
    leg2->SetTextSize(0.02);
    
    TLegendEntry* l112 = leg2->AddEntry( (TObject*)0, Form("Mean = %.2f", hz->GetMean()), "" );
    l112->SetTextColor(kBlack);
    TLegendEntry* l122 = leg2->AddEntry( (TObject*)0, Form("Std Dev = %.2f", hz->GetStdDev()), "" );
    l122->SetTextColor(kBlack);

    leg2->Draw();        

    c2->SaveAs("Z_AlphaBeta.png");   

    return 0;
}