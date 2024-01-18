#include "NeutrinoDirectionality.h"
#include "Formatting.h"
#include "Plotter.h"


int main()
{
    TH1::AddDirectory(kFALSE);

    //set the environment to apply styles
    SetRootEnv();

    vector<string> vars;
    vars.push_back("X");
    vars.push_back("Y");
    vars.push_back("Z");

    // Open the root file
    TFile *f_input = new TFile("OutputDirectionality.root", "read");

    // Angles
    double phiDeg_data, phiDegErr_data, thetaDeg_data, thetaDegErr_data;
    double phiDeg_realsim, phiDegErr_realsim, thetaDeg_realsim, thetaDegErr_realsim;
    double phiDeg_perfsim, phiDegErr_perfsim, thetaDeg_perfsim, thetaDegErr_perfsim;
    double phiDeg_true, phiDegErr_true, thetaDeg_true, thetaDegErr_true;

    // Make a histogram for each variable
    vector<TH1D*> h_data_Diff, h_realsim_Diff, h_perfsim_Diff;
    for( vector<string>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
         cout << " var " << var->c_str() << endl;

         h_data_Diff.push_back   ( (TH1D*)f_input->Get(Form("Data_Diff_%s",    var->c_str())) );
         h_realsim_Diff.push_back( (TH1D*)f_input->Get(Form("RealSim_Diff_%s", var->c_str())) );
         h_perfsim_Diff.push_back( (TH1D*)f_input->Get(Form("PerfSim_Diff_%s", var->c_str())) );
    }

    // Get angles
    TVector2 *read_phi_data = (TVector3*)f_input->Get("Data Phi");
    double phi_data = read_phi_data->X();

    TVector2 *read_theta_data = (TVector2*)f_input->Get("Data Theta");
    double theta_data = read_theta_data->X();

    TVector2 *read_phi_realsimBias = (TVector2*)f_input->Get("Simulation Unbiased Phi");
    double phi_realsimBias = read_phi_realsimBias->X();
    //double phiErr_realsimBias = read_phi_realsimBias->Y();

    TVector2 *read_theta_realsimBias = (TVector2*)f_input->Get("Simulated Unbiased Theta");
    double theta_realsimBias = read_theta_realsimBias->X();
    //double thetaErr_realsimBias = read_theta_realsimBias->Y();

    TVector2 *read_phi_true = (TVector2*)f_input->Get("True Phi");
    double phi_true = read_phi_true->X();
    double phiErr_true = read_phi_true->Y();

    TVector2 *read_theta_true = (TVector2*)f_input->Get("True Theta");
    double theta_true = read_theta_true->X();
    double thetaErr_true = read_theta_true->Y();

    TVector2 *read_phi_realsim = (TVector2*)f_input->Get("Simulation Phi");
    double phi_realsim = read_phi_realsim->X();
    //double phiErr_realsim = read_phi_realsim->Y();

    TVector2 *read_theta_realsim = (TVector2*)f_input->Get("Simulation Theta");
    double theta_realsim = read_theta_realsim->X();
    double thetaErr_realsim = read_theta_realsim->Y();

    //Close the root file
    f_input->Close();

    vector<vector<TH1D*>> h_Diff;
    h_Diff.push_back(h_data_Diff);
    h_Diff.push_back(h_realsim_Diff);
    h_Diff.push_back(h_perfsim_Diff);

    vector<string> hist;
    hist.push_back("Unbiased Data");
    hist.push_back("Unbiased RealSim");
    hist.push_back("PerfSim");

    vector<string> hist_names;
    hist_names.push_back("Data");
    hist_names.push_back("Unbiased Realistic Simulation");
    hist_names.push_back("Perfect Simulation");

    double Datax[2] = {90 - theta_data, thetaErr_data};
    double Datay[2] = {phi_data, phiErr_data};
    double RealSimx[2] = {90 - theta_realsim, thetaErr_realsim};
    double RealSimy[2] = {phi_realsim, phiErr_realsim};
    double RealSimxBias[2] = {90 - theta_realsimBias, thetaErr_realsimBias};
    double RealSimyBias[2] = {phi_realsimBias, phiErr_realsimBias};
    double PerfSimx[2] = {90 - theta_perfsim, thetaErr_perfsim};
    double PerfSimy[2] = {phi_perfsim, phiErr_perfsim};
    double Truex[2] = {90 - theta_true, thetaErr_true};
    double Truey[2] = {phi_true, phiErr_true};


    // ==========================================================
    //     Now let's make plots!
    // ==========================================================
    for( unsigned int h = 0; h < h_Diff.size(); ++h )
    {
         for( unsigned int ivar = 0; ivar < vars.size(); ++ivar )
         {
               if( vars[ivar] == "X" || vars[ivar] == "Y" )
               { 
                    TString cName = Form( "Separation_Prompt_Delayed_3bins_%s_%s", hist[h].c_str(), vars[ivar].c_str() );
                    TCanvas c( cName, cName, 1200, 800 );
                         
                    TH1D * tmp_h_IBD = new TH1D("h", "h", 3, 0.0, 3.0);
                         
                    for( int ibin = 1; ibin < h_Diff[h][ivar]->GetNbinsX(); ibin++ )
                    {
                         //cout << "bin " << ibin << ", content " << h_IBD[h][ivar]->GetBinContent(ibin) << endl;
                              
                         tmp_h_IBD->SetBinContent(1, h_Diff[h][ivar]->GetBinContent(6));
                         tmp_h_IBD->SetBinError(1, h_Diff[h][ivar]->GetBinError(6));
                         
                         tmp_h_IBD->SetBinContent(2, h_Diff[h][ivar]->GetBinContent(151));
                         tmp_h_IBD->SetBinError(2, h_Diff[h][ivar]->GetBinError(151));
                         
                         tmp_h_IBD->SetBinContent(3, h_Diff[h][ivar]->GetBinContent(296));
                         tmp_h_IBD->SetBinError(3, h_Diff[h][ivar]->GetBinError(296));
                    }

                    tmp_h_IBD->GetXaxis()->SetTitle( Form("Delayed - Prompt: %s (mm)", vars[ivar].c_str()) );
                    tmp_h_IBD->GetYaxis()->SetTitle("IBD Events / mm");
                    tmp_h_IBD->GetXaxis()->SetNdivisions(103);
                    tmp_h_IBD->GetYaxis()->SetNdivisions(509);
                    tmp_h_IBD->SetMaximum( 1.35 * tmp_h_IBD->GetMaximum() );
                    ApplyAxisStyle( tmp_h_IBD );
                    tmp_h_IBD->SetLineColor(kRed);
                    tmp_h_IBD->SetLineWidth(4);
                    tmp_h_IBD->SetLineStyle(1);
                    tmp_h_IBD->SetMarkerStyle(20);
                    tmp_h_IBD->SetMarkerSize(1.3);
                    tmp_h_IBD->SetMarkerColor(kBlack);

                    tmp_h_IBD->GetXaxis()->SetBinLabel(1, "-D");
                    tmp_h_IBD->GetXaxis()->SetBinLabel(2, "0");
                    tmp_h_IBD->GetXaxis()->SetBinLabel(3, "D");

                    tmp_h_IBD->GetXaxis()->SetLabelSize(0.075);
                    tmp_h_IBD->Draw("HISTS ][");
                    //tmp_h_IBD->Draw("P E1 X0");

                    double x_legend = 0.6;
                    double y_legend = 0.77;

                    TLegend *leg = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
                    leg->SetBorderSize(0);
                    leg->SetFillColor(0);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(62);
                    leg->SetTextSize(0.03);

		    TLegendEntry* l11 = leg->AddEntry( (TObject*)0, Form("Mean = %.2f", h_Diff[h][ivar]->GetMean()), "" );
                    l11->SetTextColor(kBlack);
                    TLegendEntry* l12 = leg->AddEntry( (TObject*)0, Form(" RMS = %.2f", h_Diff[h][ivar]->GetRMS()), "" );
                    l12->SetTextColor(kBlack);

                    leg->Draw();



                    AddHistoTitle(Form("%s", hist_names[h].c_str()), 0.05, 62);
                    c.SaveAs( Form("%s.png", c.GetName()) );
               }
               else
               {
                    TString cName = Form( "Separation_Prompt_Delayed_%s_%s", hist[h].c_str(), vars[ivar].c_str() );
                    TCanvas c( cName, cName, 1200, 800 );
                         
                    h_Diff[h][ivar]->GetXaxis()->SetTitle( Form("Delayed - Prompt: %s (mm)", vars[ivar].c_str()) );
                    h_Diff[h][ivar]->GetYaxis()->SetTitle( Form("IBD Events / %0.0f mm", h_Diff[h][ivar]->GetBinWidth(1)) );
                    h_Diff[h][ivar]->GetXaxis()->SetNdivisions(509);
                    h_Diff[h][ivar]->GetYaxis()->SetNdivisions(509);
                    h_Diff[h][ivar]->SetMaximum( 1.35 * h_Diff[h][ivar]->GetMaximum() ); 
                    ApplyAxisStyle( h_Diff[h][ivar] );
                    h_Diff[h][ivar]->SetLineColor(kBlack);
                    h_Diff[h][ivar]->SetLineWidth(3);
                    h_Diff[h][ivar]->SetLineStyle(1);
                    h_Diff[h][ivar]->SetMarkerStyle(20);
                    h_Diff[h][ivar]->SetMarkerSize(1.3);
                    h_Diff[h][ivar]->SetMarkerColor(kBlack);
                    h_Diff[h][ivar]->Draw("P E1 X0");

                    double x_legend = 0.6;
                    double y_legend = 0.77;

                    TLegend *leg = new TLegend( x_legend, y_legend, x_legend + 0.25, y_legend + 0.08 );
                    leg->SetBorderSize(0);
                    leg->SetFillColor(0);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(62);
                    leg->SetTextSize(0.03);
                    
                    TLegendEntry* l11 = leg->AddEntry( (TObject*)0, Form("Mean = %.2f", h_Diff[h][ivar]->GetMean()), "" );
                    l11->SetTextColor(kBlack);
                    TLegendEntry* l12 = leg->AddEntry( (TObject*)0, Form(" RMS = %.2f", h_Diff[h][ivar]->GetRMS()), "" );
                    l12->SetTextColor(kBlack);

                    leg->Draw();

                    AddHistoTitle(Form("%s", hist_names[h].c_str()), 0.05, 62);
                    c.SaveAs( Form("%s.png", c.GetName()) );
               }
         }
    }

    {
     TString cName = Form( "Neutrino_Angles_Unbiased" );
     TCanvas c( cName, cName, 2000, 1600 );

     TMultiGraph *mg = new TMultiGraph();
     mg->GetXaxis()->SetLimits(93.5, 108.5);
     mg->GetYaxis()->SetRangeUser(33.5, 48.5);
     mg->GetXaxis()->SetTitle("#theta (deg)");
     mg->GetYaxis()->SetTitle("#phi (deg)");
     mg->GetXaxis()->CenterTitle(kTRUE);
     mg->GetYaxis()->CenterTitle(kTRUE);
     mg->Draw("AP");

     // Uncertainties
     /* TEllipse *Data = new TEllipse(Truex[0], Truey[0], Datax[1], Datay[1]);
     Data->SetFillColor(kCyan + 1); 
     Data->SetFillStyle(3001);
     Data->Draw(); */

     // Uncorrected errors without systematics: phi = 2.89, theta = 2.84, tilt = -89.7032056 degrees
     // Uncorrected errors with systematics, phi = 4.59, theta = 3.00, tilt = -86.0763931
     // Corrected errors with systematics, phi = 4.6, theta = 2.45, tilt = -88.4317274
     // Corrected errors without systematics, phi = 3.65, theta = 2.29, tilt = -89.8796279

     TEllipse *DataFinal = new TEllipse(Datax[0], Datay[0], 4.60, 2.45, 0, 360, 271.57);
     DataFinal->SetFillColor(kAzure - 3); 
     DataFinal->SetFillStyle(3001);
     DataFinal->Draw();

     TEllipse *DataNoSys = new TEllipse(Datax[0], Datay[0], 3.65, 2.29, 0, 360, 270.23);
     DataNoSys->SetFillColor(kAzure + 4); 
     DataNoSys->SetFillStyle(3001);
     DataNoSys->Draw();
     
     TEllipse *RealSimBias = new TEllipse(RealSimxBias[0], RealSimyBias[0], RealSimxBias[1], RealSimyBias[1]);
     RealSimBias->SetFillColor(kRed + 2); 
     RealSimBias->SetFillStyle(3001);
     RealSimBias->Draw();

     TEllipse *RealSim = new TEllipse(RealSimx[0], RealSimy[0], RealSimx[1], RealSimy[1]);
     RealSim->SetFillColor(kRed + 4); 
     RealSim->SetFillStyle(3001);
     RealSim->Draw();

     /* TEllipse *PerfectSim = new TEllipse(PerfSimx[0], PerfSimy[0], PerfSimx[1], PerfSimy[1]);
     PerfectSim->SetFillColor(kRed + 2);
     PerfectSim->SetFillStyle(3001);
     PerfectSim->Draw(); */
     
     TEllipse *True = new TEllipse(Truex[0], Truey[0], Truex[1], Truey[1]);
     True->SetFillColor(kGreen + 2); 
     True->SetFillStyle(3001);
     True->Draw();

     c.Update();

     // Angles
     TGraph* grData = new TGraph(1, Datax, Datay);
     grData->SetMarkerStyle(43);
     grData->SetMarkerSize(5.5);
     grData->SetMarkerColor(kAzure - 3);
     
     TGraph* grDataNoSys = new TGraph(1, Datax, Datay);
     grDataNoSys->SetMarkerStyle(43);
     grDataNoSys->SetMarkerSize(4.5);
     grDataNoSys->SetMarkerColor(kAzure + 4);

     TGraph* grRealSimBias = new TGraph(1, RealSimxBias, RealSimyBias);
     grRealSimBias->SetMarkerStyle(43);
     grRealSimBias->SetMarkerSize(4.5);
     grRealSimBias->SetMarkerColor(kRed + 2);

     TGraph* grRealSim = new TGraph(1, RealSimx, RealSimy);
     grRealSim->SetMarkerStyle(43);
     grRealSim->SetMarkerSize(4.5);
     grRealSim->SetMarkerColor(kRed + 4);

     /* TGraph* grPerfSim = new TGraph(1, PerfSimx, PerfSimy);
     grPerfSim->SetMarkerStyle(43);
     grPerfSim->SetMarkerSize(4.5);
     grPerfSim->SetMarkerColor(kRed + 2); */

     TGraph* grTrue = new TGraph(1, Truex, Truey);
     grTrue->SetMarkerStyle(33);
     grTrue->SetMarkerSize(5);
     grTrue->SetMarkerColor(kGreen + 2);

     mg->Add(grRealSim);
     mg->Add(grRealSimBias);
     mg->Add(grData);
     mg->Add(grTrue);
     //mg->Add(grPerfSim);
     c.Update();
     
     TLegend *leg = new TLegend(0.48, 0.7, 0.72, 0.88);
     leg->SetBorderSize(0);
     leg->SetTextFont(62);
     leg->SetTextSize(0.035);
     leg->AddEntry(grData, "Data", "p");
     leg->AddEntry(grDataNoSys, "Data - No Systematics", "p");
     leg->AddEntry(grRealSimBias, "Simulation", "p");
     leg->AddEntry(grRealSim, "Uncorrected Sim", "p");
     //leg->AddEntry(grPerfSim, "Sim - no inop. segments", "p");
     leg->AddEntry(grTrue, "True Neutrino Direction", "p");
     leg->Draw();

     AddHistoTitle("Average Reconstructed Neutrino Direction", 0.05, 62);
     c.SaveAs( Form("%s_PaperFinalPlot.png", c.GetName()) );
    }


    printf("\nProgram Complete!\n\n");

    return 0;

} // End Program
