#include <TSystem.h>
#include <stdio.h>
#include "neutrinoDir.h"
#include "GeneralHeader.h"

using namespace std;
#define use_atan2

// Global Constants
const int h_bins = 301;
const int index_offset = 150;
float h_max = 150.5;
double delayed[3] = {0};
double prompt[3] = {0};

// Defining Structs
struct DirectionalityValues{
    // Deadtime corrected livetime and reactor runtime
    double xyz = 0;
    double liveTimeOff = 0;
    double liveTimeOn = 0;
    double totalIBD = 0;
    double totalIBDErr = 0;
    double effIBD = 0;

    // Correlated and Accidentals
    double corOff[h_bins] = {0};
    double accOff[h_bins] = {0};
    double corOn[h_bins] = {0};
    double accOn[h_bins] = {0};

    double corOffLTF[h_bins] = {0};
    double accOffLTF[h_bins] = {0};
    double corOnLTF[h_bins] = {0};
    double accOnLTF[h_bins] = {0};

};

typedef struct DirectionalityValues Struct;

// Data

Struct FindDirValues(const char *fName, const char *TFileName, const char *release, int xyorz = 0){

    Struct s;
    s.xyz = xyorz;

    cout.precision(10);

    // TTree

    double corOnRx = 0;
    double corOnpp = 0;
    double nCaptDt = 0;
    double Esmear = 0;

    int pSeg = 0;
    int nSeg = 0;
    int roundIndex = 0;
    
    float diffIndex = 0;

    long npOrigCorOff = 0;
    long nmOrigCorOff = 0;
    long n0OrigCorOff = 0;
    long npOrigAccOff = 0;
    long nmOrigAccOff = 0;
    long n0OrigAccOff = 0;

    long npOrigCorOn = 0;
    long nmOrigCorOn = 0;
    long n0OrigCorOn = 0;
    long npOrigAccOn = 0;
    long nmOrigAccOn = 0;
    long n0OrigAccOn = 0;

    long n0DiffCorOff = 0;
    long n0DiffAccOff = 0;
    long n0DiffCorOn = 0;
    long n0DiffAccOn = 0;

    TH1F *hDataCorOff = new TH1F("hDataCorOff", "Correlated Off", h_bins, -h_max, h_max);
    TH1F *hDataAccOff = new TH1F("hDataAccOff", "Accidentals Off", h_bins, -h_max, h_max);
    TH1F *hDataCorOn = new TH1F("hDataCorOn", "Correlated On", h_bins, -h_max, h_max);
    TH1F *hDataAccOn = new TH1F("hDataAccOn", "Accidentals On", h_bins, -h_max, h_max);

    /* Correlated IBDs involve a neutron capture (delayed signal) within 1-120 µs of the position (prompt signal).
       Accidental IBDs have the delayed signal > 1ms after the prompt.*/

    // File manipulation (opening)
    ifstream file;
    file.open(fName, ifstream::in);

    if(!(file.is_open() && file.good())){
        cout << "The file is either bad or not found. Exiting.";
        return s;
    }
    
    while(file.good() && (!file.eof())){
        string line;
        getline(file, line);

        TString str = Form("%s%s%s%s", gSystem->Getenv("VETO_OUTDIR"), release, line.data(), TFileName);
        str.Remove(0,6);
        // Correcting a formatting issue within file directory
        // 0 -> Reactor off, 1 -> Reactor on

        if (str.Contains("0")){
            str.ReplaceAll("0", "");
            TFile *f = new TFile(str);

            TVectorD *rtOff = ((TVectorD*)f->Get("runtime"));
            TVectorD *promptvOff = ((TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_prompt"));
            TVectorD *delayedvOff = ((TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_delayed"));

            double xRtOff = (rtOff->Max() / (rtOff->Max() - promptvOff->Max())) * (rtOff->Max() / (rtOff->Max() - delayedvOff->Max()));

            TTree *Tree = (TTree*)f->Get("P2kIBDPlugin/Tibd");

            long nEntries = Tree->GetEntries();
            cout << "We have " << nEntries << " entries\n";

            Double_t n_xyz, xyz, Esmear;

            Tree->SetBranchAddress("n_xyz", &n_xyz);

            for (long i = 0; i < nEntries; i++){
                
                // Esmear is prompt energy and it has to be within PROSPECT boundary conditions.
                Tree->GetEntry(i);

                delayed[xyorz] = Tree->GetLeaf("n_xyz")->GetValue(xyorz);

            }
        }

    }



}