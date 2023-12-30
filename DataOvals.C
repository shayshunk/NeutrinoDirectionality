#include <iostream>
#include "TSystem.h"
#include "GeneralHeader.h"
#include <string>
#include "TLegend.h"
#include "TBox.h"
#include "TEllipse.h"
// X is Theta
// Y is Phi
void DataOvals(){
// C is Christian
// J is Jin
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000, 2000, 1000);
	c1->cd(); gStyle->SetOptStat(0);
        // 90 degrees - atan(y/x):
        //double Realx[2] = {-8.59, 0.71}, Realy[2] = {49.29, 0.73};
        // atan2(y,x):
        double Realx[2] = {-8.625, 0.71}, Realy[2] = {40.451, 0.73};
	// Christian's Values
	double CxData[2] = {-10.98, 2.90}, CyData[2] = {49.01, 3.43};
	double CxSim[2] = {-12.31, 0.49}, CySim[2] = {46.65, 0.60};
	double CxPer[2] = {-9.33, 0.31}, CyPer[2] = {49.35, 0.37};

	// My (Jin's) Values
        TFile *fData = new TFile("OutputDirectionality.root");
        TTree *TData = (TTree*)fData->Get("Tnd");
        TFile *fSim = new TFile("Simulations/OutputDirSim.root");
        TTree *TSim = (TTree*)fSim->Get("Tnd");
        TFile *fPer = new TFile("Simulations/OutputDirPer.root");
        TTree *TPer = (TTree*)fPer->Get("Tnd");

        double Datax[2], Datay[2]; TData->GetEntry(0);
        Datax[0] = TData->GetLeaf("thetaDeg")->GetValue(0);
        Datax[1] = TData->GetLeaf("thetaDegErr")->GetValue(0);
        Datay[0] = TData->GetLeaf("phiDeg")->GetValue(0);
        Datay[1] = TData->GetLeaf("phiDegErr")->GetValue(0);

        double Simx[2], Simy[2]; TSim->GetEntry(0);
        // Simx[0] = -13.34; Simx[1] = 0.48;
        Simx[0] = TSim->GetLeaf("thetaDeg")->GetValue(0);
        Simx[1] = TSim->GetLeaf("thetaDegErr")->GetValue(0);
        // Simy[0] = 43.32; Simy[1] = 0.56;
        Simy[0] = TSim->GetLeaf("phiDeg")->GetValue(0);
        Simy[1] = TSim->GetLeaf("phiDegErr")->GetValue(0);

        double Perx[2], Pery[2]; TPer->GetEntry(0);
        // Perx[0] = -9.01; Perx[1] = 0.31;
        Perx[0] = TPer->GetLeaf("thetaDeg")->GetValue(0);
        Perx[1] = TPer->GetLeaf("thetaDegErr")->GetValue(0);
        // Pery[0] = 40.61; Pery[1] = 0.31;
        Pery[0] = TPer->GetLeaf("phiDeg")->GetValue(0);
        Pery[1] = TPer->GetLeaf("phiDegErr")->GetValue(0);

        TGraph* grReal = new TGraph(1, Realx, Realy);
        grReal->SetMarkerStyle(33);
        grReal->SetMarkerSize(4);
        grReal->SetMarkerColor(kGreen + 1);

	TGraph* CgrData = new TGraph(1, CxData, CyData);
	CgrData->SetMarkerStyle(29);
	CgrData->SetMarkerSize(3);
	CgrData->SetMarkerColor(kBlue);

	TGraph* CgrSim = new TGraph(1, CxSim, CySim);
	CgrSim->SetMarkerStyle(29);
	CgrSim->SetMarkerSize(3);
	CgrSim->SetMarkerColor(kRed);

	TGraph* CgrPer = new TGraph(1, CxPer, CyPer);
	CgrPer->SetMarkerStyle(29);
	CgrPer->SetMarkerSize(3);
	CgrPer->SetMarkerColor(kPink + 8);

        TGraph* grData = new TGraph(1, Datax, Datay);
        grData->SetMarkerStyle(39);
        grData->SetMarkerSize(3);
        grData->SetMarkerColor(kViolet - 6);

        TGraph* grSim = new TGraph(1, Simx, Simy);
        grSim->SetMarkerStyle(39);
        grSim->SetMarkerSize(3);
        grSim->SetMarkerColor(kOrange + 1);

        TGraph* grPer = new TGraph(1, Perx, Pery);
        grPer->SetMarkerStyle(39);
        grPer->SetMarkerSize(3);
        grPer->SetMarkerColor(kPink - 9);

	TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle("Angle Measurements of Data & Simulation");
	mg->Add(grReal);
	//mg->Add(CgrData);
	//mg->Add(CgrSim);
	//mg->Add(CgrPer);
	mg->Add(grData);
	mg->Add(grSim);
	mg->Add(grPer);
	mg->Draw("AP");

	mg->GetXaxis()->SetLimits(-27, 0);
	mg->GetYaxis()->SetRangeUser(30, 60.5);
	mg->GetXaxis()->SetTitle("#theta (degrees)");
	mg->GetYaxis()->SetTitle("#phi (degrees)");
	mg->Draw("AP");
	c1->Update();

	// Uncertainties
	TEllipse *CDatao = new TEllipse(CxData[0], CyData[0], CxData[1], CyData[1]);
	CDatao->SetFillColor(kBlue); CDatao->SetFillStyle(3002);
	//CDatao->Draw();

	TEllipse *CSimo = new TEllipse(CxSim[0], CySim[0], CxSim[1], CySim[1]);
	CSimo->SetFillColor(kRed); CSimo->SetFillStyle(3001);
	//CSimo->Draw();

	TEllipse *CPero = new TEllipse(CxPer[0], CyPer[0], CxPer[1], CyPer[1]);
	CPero->SetFillColor(kPink + 8); CPero->SetFillStyle(3003);
	//CPero->Draw();

        TEllipse *Datao = new TEllipse(Datax[0], Datay[0], Datax[1], Datay[1]);
        Datao->SetFillColor(kViolet - 6); Datao->SetFillStyle(3002);
        Datao->Draw();

        TEllipse *Simo = new TEllipse(Simx[0], Simy[0], Simx[1], Simy[1]);
        Simo->SetFillColor(kOrange + 1); Simo->SetFillStyle(3001);
        Simo->Draw();

        TEllipse *Pero = new TEllipse(Perx[0], Pery[0], Perx[1], Pery[1]);
        Pero->SetFillColor(kPink - 9); Pero->SetFillStyle(3003);
        Pero->Draw();

	TEllipse *Realo = new TEllipse(Realx[0], Realy[0], Realx[1], Realy[1]);
	Realo->SetFillColor(kGreen + 1); Realo->SetFillStyle(3002);
	Realo->Draw();

	auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
	//legend->AddEntry(CgrData, "Christian Data", "p");
	//legend->AddEntry(CgrSim, "Christian Sim", "p");
	//legend->AddEntry(CgrPer, "Christian Perfect Sim", "p");
        //legend->SetTextSize(30);
	legend->AddEntry(grData, "Data", "p");
	legend->AddEntry(grSim, "Realistic Sim", "p");
	legend->AddEntry(grPer, "Perfect Sim", "p");
	legend->AddEntry(grReal, "True Reactor Direction", "p");
	legend->Draw();
	c1->SaveAs("DataOvals.png");
	printf("\nProgram Complete \n");
}


