// Tree Structure
#include<TTree.h>
using namespace std;

struct struct_ND{
    // Initializing variables from PromptDelayedDataSeaparation0813.C
    
    Double_t liveTimeOff = 0;
    Double_t liveTimeOn = 0;

    // Determining IBD Count
    Double_t totalIBD = 0;
    Double_t totalIBDerr = 0;
    Double_t effIBD = 0;

    Double_t momX = 0;
    Double_t momY = 0;
    Double_t momZ = 0;

    Double_t sigmaX = 0;
    Double_t sigmaY = 0;
    Double_t sigmaZ = 0;

    // Initializing vars from AngleFormula.C
    Double_t phiDeg;
    Double_t phiErr;
    Double_t thetaDeg;
    Double_t thetaErr;

};

void setupNDTreeFill(TTree& T_ND, struct_ND& s_ND){

    T_ND.Branch("liveTimeOff", &s_ND.liveTimeOff, "liveTimeOff/D");
    T_ND.Branch("liveTimeOn", &s_ND.liveTimeOn, "liveTimeOn/D");
    T_ND.Branch("totalIBD", &s_ND.totalIBD, "totalIBD/D");
    T_ND.Branch("totalIBDErr", &s_ND.totalIBDerr, "totalIBDErr/D");
    T_ND.Branch("effIBD", &s_ND.effIBD, "effIBD/D");
    T_ND.Branch("momX", &s_ND.momX, "momX/D");
    T_ND.Branch("momY", &s_ND.momY, "momY/D");
    T_ND.Branch("momZ", &s_ND.momZ, "momZ/D");
    T_ND.Branch("sigmaX", &s_ND.sigmaX, "sigmaX/D");
    T_ND.Branch("sigmaY", &s_ND.sigmaY, "sigmaY/D");
    T_ND.Branch("sigmaZ", &s_ND.sigmaZ, "sigmaZ/D");
    T_ND.Branch("phiDeg", &s_ND.phiDeg, "phiDeg/D");
    T_ND.Branch("phiErr", &s_ND.phiErr, "phiErr/D");
    T_ND.Branch("thetaDeg", &s_ND.thetaDeg, "thetaDeg/D");
    T_ND.Branch("thetaErr", &s_ND.thetaErr, "thetaErr/D");

}

void setupNDTreeRead(TTree& T_ND, struct_ND& s_ND){
    T_ND.SetBranchAddress("liveTimeOff", &s_ND.liveTimeOff);
    T_ND.SetBranchAddress("liveTimeOn", &s_ND.liveTimeOn);
    T_ND.SetBranchAddress("totalIBD", &s_ND.totalIBD);
    T_ND.SetBranchAddress("totalIBDErr", &s_ND.totalIBDerr);
    T_ND.SetBranchAddress("effIBD", &s_ND.effIBD);
    T_ND.SetBranchAddress("momX", &s_ND.momX);
    T_ND.SetBranchAddress("momY", &s_ND.momY);
    T_ND.SetBranchAddress("momZ", &s_ND.momZ);
    T_ND.SetBranchAddress("sigmaX", &s_ND.sigmaX);
    T_ND.SetBranchAddress("sigmaY", &s_ND.sigmaY);
    T_ND.SetBranchAddress("sigmaZ", &s_ND.sigmaZ);
    T_ND.SetBranchAddress("phiDeg", &s_ND.phiDeg);
    T_ND.SetBranchAddress("phiErr", &s_ND.phiErr);
    T_ND.SetBranchAddress("thetaDeg", &s_ND.thetaDeg);
    T_ND.SetBranchAddress("thetaErr", &s_ND.thetaErr);
}