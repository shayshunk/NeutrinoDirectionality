///////////////////
// Jin Oueslati
///////////////////

// Finish error equations for angles
// New Version - Manipulating Bin Content and Error, and usage of each error instead of gaussian width P
// diff_index == 0 does not equal pseg == nseg b/c the first accounts for the other directional adjacent bins; while doing x: dy == 0 for up and down are included in center bin.
#include <TSystem.h>
#include <stdio.h>
#include "NeutrinoDirectionality.h"
#include "GeneralHeader.h"
#include "SegSignal.h"
//#include "PerSegCrossCheck.h" // Cross check header
using namespace std;
#define use_atan2
#define use_avg_equation
#define use_subtraction
//#define use_complex_errors
// Global Constants
int z_max = 450, z_bins = (2*z_max), x_max = 850, x_bins = (2*z_max), index_offset = 150, round_index = 0;
const int h_bins = 301, nsegments = 154;
float h_max = 150.5, diff_index = 0, d_seg = 145.7; // width of the segment (mm)
double delayed[3] = {0}, prompt[3] = {0}, dataXArray[h_bins] = {0}, dataYArray[h_bins] = {0}, ErrXArray[h_bins] = {0}, ErrYArray[h_bins] = {0}, deltaXSquaredArray[3] = {0}, deltaYSquaredArray[3] = {0};
double rX_pos = 0, rX_neg = 0, rY_pos = 0, rY_neg = 0, n0_center = 0; // Ratios and center count
TH1F *hX = new TH1F("hX", "Prompt & Delayed Captured Events X Direction", x_bins, -x_max, x_max);
TH1F *hY = new TH1F("hY", "Prompt & Delayed Captured Events Y Direction", x_bins, -x_max, x_max);
TH1F *hZ = new TH1F("hZ", "Prompt & Delayed Captured Events Z Direction", z_bins, -z_max, z_max);

//---------------------------------------------------------------------------------------
struct DirectionalityValues{
  // Deadtime Corrected Livetime & Reactor Total Runtime
  double xyz=0, livetimeOff=0, averagexyz=0, stddev=0, livetimeOn=0, totalIBDs=0, totalIBDsErr=0, effIBDs=0;
  // Correlated and Accidentals
  double corOffErrArray[h_bins]={0}, accOffErrArray[h_bins]={0}, corOnErrArray[h_bins]={0}, accOnErrArray[h_bins]={0};
  double corOffLTFArray[h_bins]={0}, accOffLTFArray[h_bins]={0}, corOnLTFArray[h_bins]={0}, accOnLTFArray[h_bins]={0};
};
typedef struct DirectionalityValues Struct;

// My Data Version
Struct FindDirectionalityValues(const char *fname, const char *TFilename, const char *release, int xyorz=0){
  Struct s;
  cout.precision(10);
  // 0 - x, 1 - y, 2 - z
  s.xyz = xyorz;
  // Use arrays to add and histograms to count
  // Setting the bins on the actual values
  TH1F *hData_corOff = new TH1F("hData_corOff", "", h_bins, -h_max, h_max);
  TH1F *hData_accOff = new TH1F("hData_accOff", "", h_bins, -h_max, h_max);
  TH1F *hData_corOn = new TH1F("hData_corOn", "", h_bins, -h_max, h_max);
  TH1F *hData_accOn = new TH1F("hData_accOn", "", h_bins, -h_max, h_max);

  double ncapt_dt = 0, Esmear = 0, xRxOn=0, xRxOff=0;
  int nseg = 0, pseg = 0;
  double delta_XL_squared = 0, delta_XR_squared = 0, delta_X0_squared = 0, delta_YL_squared = 0, delta_YR_squared = 0, delta_Y0_squared = 0; // Errors associated with each direction
  // Necessary for equation calculation
  double np_corOff=0, np0_corOff=0, nm_corOff=0, nm0_corOff=0, n0_corOff=0;
  double np_accOff=0, np0_accOff=0, nm_accOff=0, nm0_accOff=0, n0_accOff=0;
  double np_corOn=0, np0_corOn=0, nm_corOn=0, nm0_corOn=0, n0_corOn=0;
  double np_accOn=0, np0_accOn=0, nm_accOn=0, nm0_accOn=0, n0_accOn=0;

  // Original - No split between X & Y Dimension; Y adjacent bins concantenate with true center bin the X dimension and vice versa
  double np_orig_corOff=0, nm_orig_corOff=0, n0_orig_corOff=0, np_orig_accOff=0, nm_orig_accOff=0, n0_orig_accOff=0, np_orig_corOn=0, nm_orig_corOn=0, n0_orig_corOn=0, np_orig_accOn=0, nm_orig_accOn=0, n0_orig_accOn=0;
  // Adjacent or Same Segments:
  bool adj_up = 0, adj_down = 0, adj_right = 0, adj_left = 0, same_up = 0, same_down = 0, same_right = 0, same_left = 0;
  // struct in SegmentSignal.h: Segment list seg_w and bool list has_seg.
  s_SegmentAdjacentSegmentInfo sasi;
  fill_segment_info(sasi);
  //s_PerfectAdjacentSegmentInfo pasi; fill_perfect_segment_info(pasi);
  cout << "\n\n0 for x, 1 for y, 2 for z: " << xyorz << endl;
  // Opening File
  std::ifstream file;
  file.open(fname, std::ifstream::in);
  if (!(file.is_open()&&file.good())){
    printf("Good runs file not found. Exiting\n");
    //return -1;
    return s;
  }
  while (file.good()&!file.eof()){
    string line;
    getline(file, line);
    TString st = Form("%s/%s/%s/%s",gSystem->Getenv("VETO_OUTDIR"), release, line.data(), TFilename);
    st.Remove(0,6); // Formating issue within file directory to correct it.
    // 0 refers to reactor off while 1 is reactor on.
    if (st.Contains(" 0")){
      st.ReplaceAll(" 0", "");
      TFile *f=new TFile(st);
      TVectorD *rtOff = ((TVectorD*)f->Get("runtime"));
      TVectorD *promptvOff = ((TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_prompt"));
      TVectorD *delayedvOff = ((TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_delayed")); // Singular Elements in each given vector (Max).
      xRxOff = (rtOff->Max() / (rtOff->Max() - promptvOff->Max() ) ) * (rtOff->Max() / (rtOff->Max() - delayedvOff->Max() ) );
      TTree *Th = (TTree*)f->Get("P2kIBDPlugin/Tibd");
      long nentries = Th->GetEntries();
      //cout << "Entries: " << nentries << endl;
      for (long i = 0; i < nentries; i++){
        // Esmear refers to Prompt Energy and it is within Prospect's boundary conditions. Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal.
        Th->GetEntry(i);
        delayed[xyorz] = Th->GetLeaf("n_xyz")->GetValue(xyorz);
        prompt[xyorz] = Th->GetLeaf("xyz")->GetValue(xyorz);
        Esmear = Th->GetLeaf("Esmear")->GetValue(0);
	      ncapt_dt = Th->GetLeaf("ncapt_dt")->GetValue(0);
        nseg = Th->GetLeaf("n_seg")->GetValue(0);
        pseg = Th->GetLeaf("maxseg")->GetValue(0);
        // Reset Booleans to False
        adj_right = 0; adj_left = 0; same_right = 0; same_left = 0; adj_up = 0; adj_down = 0; same_up = 0; same_down = 0;
        // Correlated Off
        if (Esmear > 0.8 && Esmear < 7.2 && ncapt_dt > pow(10, 3) && ncapt_dt < 120 * pow(10,3)){
          diff_index = (delayed[xyorz] - prompt[xyorz]); // 0, -D, +D in X & Y
          round_index = round(delayed[xyorz] - prompt[xyorz]);
	        // hData histos fill the NON-LTF Arrays. DeadTime Correction should be applied here not later for accidentals.
          hData_corOff->Fill(diff_index);
          // LTF - Live Time Fraction
          s.corOffLTFArray[round_index + index_offset] += 1.;
          if (pseg == nseg) n0_orig_corOff++;
          else if (diff_index > 0) np_orig_corOff++;
          else if (diff_index < 0) nm_orig_corOff++;
          #ifdef use_subtraction
            if (pseg == nseg) n0_corOff++;
            if (xyorz == 0){
              hX->Fill(prompt[xyorz]); 
              hX->Fill(delayed[xyorz]);
              if (pseg == nseg){same_right = sasi.has_right_seg[pseg]; same_left = sasi.has_left_seg[pseg];}
              else if (diff_index > 0) adj_right = sasi.has_right_seg[pseg];
              else if (diff_index < 0) adj_left = sasi.has_left_seg[pseg];
              // if in given list and in same (adjacent) segment:
              if(same_right) np0_corOff++; // r(+)
              if(same_left) nm0_corOff++; // r(-)
              if(adj_right) np_corOff++; // r(+)
              if(adj_left) nm_corOff++; // r(-)
            }
          else if (xyorz == 1){
            hY->Fill(prompt[xyorz]); 
            hY->Fill(delayed[xyorz]);
            if (pseg == nseg){same_up = sasi.has_up_seg[pseg]; same_down = sasi.has_down_seg[pseg];}
            else if (diff_index > 0) adj_up = sasi.has_up_seg[pseg];
            else if (diff_index < 0) adj_down = sasi.has_down_seg[pseg];
            if(same_up) np0_corOff++; // r(+)
            if(same_down) nm0_corOff++; // r(-)
            if(adj_up) np_corOff++; // r(+)
            if(adj_down) nm_corOff++; // r(-)
          }
          else if (xyorz == 2){hZ->Fill(prompt[xyorz]); hZ->Fill(delayed[xyorz]);}
#endif
        }
        // Accidentals Off (IBDs with neutrons capturing over 1 ms after the prompt).
        else if (Esmear > 0.8 && Esmear < 7.2 && ncapt_dt > pow(10, 6)){
          diff_index = (delayed[xyorz] - prompt[xyorz]);
          round_index = round(delayed[xyorz] - prompt[xyorz]);
          hData_accOff->Fill(diff_index);
          // LTF - Live Time Fraction
          s.accOffLTFArray[round_index + index_offset] += 1.*(xRxOff);
	  // Original Comparison
	  if (pseg == nseg) n0_orig_accOff += 1.*(xRxOff);
	  else if (diff_index > 0) np_orig_accOff += 1.*(xRxOff);
	  else if (diff_index < 0) nm_orig_accOff += 1.*(xRxOff);
#ifdef use_subtraction
          if (pseg == nseg) n0_accOff += 1.*(xRxOff);
          if (xyorz == 0){
            hX->Fill(prompt[xyorz]); hX->Fill(delayed[xyorz]);
            if (pseg == nseg){same_right = sasi.has_right_seg[pseg]; same_left = sasi.has_left_seg[pseg];}
            else if (diff_index > 0) adj_right = sasi.has_right_seg[pseg];
            else if (diff_index < 0) adj_left = sasi.has_left_seg[pseg];
            if(same_right) np0_accOff += 1.*(xRxOff); // r(+)
            if(same_left) nm0_accOff += 1.*(xRxOff); // r(-)
            if(adj_right) np_accOff += 1.*(xRxOff); // r(+)
            if(adj_left) nm_accOff += 1.*(xRxOff); // r(-)
          }
          else if (xyorz == 1){
            hY->Fill(prompt[xyorz]); hY->Fill(delayed[xyorz]);
            if (pseg == nseg){same_up = sasi.has_up_seg[pseg]; same_down = sasi.has_down_seg[pseg];}
            else if (diff_index > 0) adj_up = sasi.has_up_seg[pseg];
            else if (diff_index < 0) adj_down = sasi.has_down_seg[pseg];
            if(same_up) np0_accOff += 1.*(xRxOff); // r(+)
            if(same_down) nm0_accOff += 1.*(xRxOff); // r(-)
            if(adj_up) np_accOff += 1.*(xRxOff); // r(+)
            if(adj_down) nm_accOff += 1.*(xRxOff); // r(-)
          }
          else if (xyorz == 2){hZ->Fill(prompt[xyorz]); hZ->Fill(delayed[xyorz]);}
#endif
        }
      }
      // Returns the next character in the input sequence, without extracting it: The character is left as the next character to be extracted from the stream.
      file.peek();
      f->Close();
      s.livetimeOff += (rtOff->Max())/(xRxOff);
    }
    else if (st.Contains(" 1")){
      st.ReplaceAll(" 1", "");
      TFile *f=new TFile(st);
      TVectorD *rtOn = ((TVectorD*)f->Get("runtime"));
      TVectorD *promptvOn = ((TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_prompt"));
      TVectorD *delayedvOn = ((TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_delayed"));
      xRxOn = (rtOn->Max() / (rtOn->Max() - promptvOn->Max() ) ) * (rtOn->Max() / (rtOn->Max() - delayedvOn->Max() ) );
      TTree *Th = (TTree*)f->Get("P2kIBDPlugin/Tibd");
      long nentries = Th->GetEntries();
      for (long i = 0; i < nentries; i++){
        Th->GetEntry(i);
        delayed[xyorz] = Th->GetLeaf("n_xyz")->GetValue(xyorz);
        prompt[xyorz] = Th->GetLeaf("xyz")->GetValue(xyorz);
        Esmear = Th->GetLeaf("Esmear")->GetValue(0);
        ncapt_dt = Th->GetLeaf("ncapt_dt")->GetValue(0);
        nseg = Th->GetLeaf("n_seg")->GetValue(0);
        pseg = Th->GetLeaf("maxseg")->GetValue(0);
        adj_right = 0; adj_left = 0; same_right = 0; same_left = 0; adj_up = 0; adj_down = 0; same_up = 0; same_down = 0; // Reset Booleans to False
        // Correlated On
        if (Esmear > 0.8 && Esmear < 7.2 && ncapt_dt > pow(10, 3) && ncapt_dt < 120 * pow(10,3)){
          diff_index = (delayed[xyorz] - prompt[xyorz]);
          round_index = round(delayed[xyorz] - prompt[xyorz]);
          hData_corOn->Fill(diff_index);
          // LTF - Live Time Fraction
          s.corOnLTFArray[round_index + index_offset] += 1.;
	  // Original Comparison
	  if (pseg == nseg) n0_orig_corOn++;
	  else if (diff_index > 0) np_orig_corOn++;
	  else if (diff_index < 0) nm_orig_corOn++;
#ifdef use_subtraction
          if (pseg == nseg) n0_corOn++;
          if (xyorz == 0){
            hX->Fill(prompt[xyorz]); hX->Fill(delayed[xyorz]);
            if (pseg == nseg){same_right = sasi.has_right_seg[pseg]; same_left = sasi.has_left_seg[pseg];}
            else if (diff_index > 0) adj_right = sasi.has_right_seg[pseg];
            else if (diff_index < 0) adj_left = sasi.has_left_seg[pseg];
            if(same_right) np0_corOn++; // r(+)
            if(same_left) nm0_corOn++; // r(-)
            if(adj_right) np_corOn++; // r(+)
            if(adj_left) nm_corOn++; // r(-)
          }
          else if (xyorz == 1){
            hY->Fill(prompt[xyorz]); hY->Fill(delayed[xyorz]);
            if (pseg == nseg){same_up = sasi.has_up_seg[pseg]; same_down = sasi.has_down_seg[pseg];}
            else if (diff_index > 0) adj_up = sasi.has_up_seg[pseg];
            else if (diff_index < 0) adj_down = sasi.has_down_seg[pseg];
            if(same_up) np0_corOn++; // r(+)
            if(same_down) nm0_corOn++; // r(-)
            if(adj_up) np_corOn++; // r(+)
            if(adj_down) nm_corOn++; // r(-)
          }
          else if (xyorz == 2){hZ->Fill(prompt[xyorz]); hZ->Fill(delayed[xyorz]);}
#endif
        }
        // Accidental On
        else if (Esmear > 0.8 && Esmear < 7.2 && ncapt_dt > pow(10, 6)){
          diff_index = (delayed[xyorz] - prompt[xyorz]);
          round_index = round(delayed[xyorz] - prompt[xyorz]);
          hData_accOn->Fill(diff_index);
          // LTF - Live Time Fraction
          s.accOnLTFArray[round_index + index_offset] += 1.*(xRxOn);
	  // Original Comparison
	  if (pseg == nseg) n0_orig_accOn += 1.*(xRxOn);
	  else if (diff_index > 0) np_orig_accOn += 1.*(xRxOn);
	  else if (diff_index < 0) nm_orig_accOn += 1.*(xRxOn);
#ifdef use_subtraction
          if (pseg == nseg) n0_accOn += 1.*(xRxOn);
          if (xyorz == 0){
            hX->Fill(prompt[xyorz]); hX->Fill(delayed[xyorz]);
            if (pseg == nseg){same_right = sasi.has_right_seg[pseg]; same_left = sasi.has_left_seg[pseg];}
            else if (diff_index > 0) adj_right = sasi.has_right_seg[pseg];
            else if (diff_index < 0) adj_left = sasi.has_left_seg[pseg];
            if(same_right) np0_accOn += 1.*(xRxOn); // r(+)
            if(same_left) nm0_accOn += 1.*(xRxOn); // r(-)
            if(adj_right) np_accOn += 1.*(xRxOn); // r(+)
            if(adj_left) nm_accOn += 1.*(xRxOn); // r(-)
          }
          else if (xyorz == 1){
            hY->Fill(prompt[xyorz]); hY->Fill(delayed[xyorz]);
            if (pseg == nseg){same_up = sasi.has_up_seg[pseg]; same_down = sasi.has_down_seg[pseg];}
            else if (diff_index > 0) adj_up = sasi.has_up_seg[pseg];
            else if (diff_index < 0) adj_down = sasi.has_down_seg[pseg];
            if(same_up) np0_accOn += 1.*(xRxOn); // r(+)
            if(same_down) nm0_accOn += 1.*(xRxOn); // r(-)
            if(adj_up) np_accOn += 1.*(xRxOn); // r(+)
            if(adj_down) nm_accOn += 1.*(xRxOn); // r(-)
          }
          else if (xyorz == 2){hZ->Fill(prompt[xyorz]); hZ->Fill(delayed[xyorz]);}
#endif
        }
      }
      file.peek();
      f->Close();
      s.livetimeOn += (rtOn->Max())/(xRxOn);
    }
  }
  // Used for Z Calculations
  // Additional Factors to calculate total IBD Count
  double totalIBDscorOff = 0, totalIBDsaccOff = 0, totalIBDscorOn = 0, totalIBDsaccOn = 0;
  double npAvg = 0, np0Avg = 0, nmAvg = 0, nm0Avg = 0, n0Avg = 0, n_prime = 0;
  // Comment this loop out if you are not doing it via the histogram method
  for (int k = 0; k < h_bins; k++){
    s.corOffErrArray[k] = hData_corOff->GetBinContent(k+1);
    s.accOffErrArray[k] = hData_accOff->GetBinContent(k+1);
    s.corOnErrArray[k] = hData_corOn->GetBinContent(k+1);
    s.accOnErrArray[k] = hData_accOn->GetBinContent(k+1);
  }
  for (unsigned int i = 0; i < h_bins; i++){totalIBDscorOff += s.corOffLTFArray[i]; totalIBDsaccOff += s.accOffLTFArray[i]; totalIBDscorOn += s.corOnLTFArray[i]; totalIBDsaccOn += s.accOnLTFArray[i];}
  // Should be correct if the 10^4us factor is right for the time windows
  s.totalIBDs = (totalIBDscorOn - totalIBDsaccOn/100.0) - (s.livetimeOn/s.livetimeOff) * (totalIBDscorOff - totalIBDsaccOff/100.0);
  s.totalIBDsErr = sqrt(totalIBDscorOn + totalIBDsaccOn/10000.0 + pow((s.livetimeOn/s.livetimeOff), 2) * (totalIBDscorOff + totalIBDsaccOff/10000.0));
  s.effIBDs = pow(s.totalIBDs, 2)/pow(s.totalIBDsErr, 2);

#ifdef use_subtraction
    // Used for X & Y Calculations
  if (xyorz == 0 || xyorz == 1){
    /*
    // Comparison to NeutrinoDirectionality.cc
    cout << "Original Measurements:" << endl;
    double npOrigAvg = ((np_orig_corOn - np_orig_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(np_orig_corOff - np_orig_accOff/100.0));
    double nmOrigAvg = ((nm_orig_corOn - nm_orig_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(nm_orig_corOff - nm_orig_accOff/100.0));
    double n0OrigAvg = ((n0_orig_corOn - n0_orig_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(n0_orig_corOff - n0_orig_accOff/100.0));
    cout << "np Orig Avg: " << npOrigAvg << " nm Orig Avg: " << nmOrigAvg << " n0 OrigAvg: " << n0OrigAvg << "\n" << endl;
    */

    // Dead Time Correction for accidentals is done when looping through files, accidental correction time is 100 microseconds which can be distributed here; along with livetime for Rx time
    npAvg = ((np_corOn - np_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(np_corOff - np_accOff/100.0));
    np0Avg = ((np0_corOn - np0_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(np0_corOff - np0_accOff/100.0));
    nmAvg = ((nm_corOn - nm_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(nm_corOff - nm_accOff/100.0));
    nm0Avg = ((nm0_corOn - nm0_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(nm0_corOff - nm0_accOff/100.0));
    n0Avg = ((n0_corOn - n0_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(n0_corOff - n0_accOff/100.0));
    n0_center = n0Avg;
    cout << "Unbiased Measurements:" << endl;
    //cout << "np corOff " << np_corOff << " nm corOff " << nm_corOff << " np0 corOff " << np0_corOff << " nm0 corOff " << nm0_corOff << " n0 corOff " << n0_corOff << endl;
    //cout << "np accOff " << np_accOff/100.0 << " nm accOff " << nm_accOff/100.0 << " np0 accOff " << np0_accOff/100.0 << " nm0 accOff " << nm0_accOff/100.0 << " n0 accOff " << n0_accOff/100.0 << endl;
    //cout << "np corOn " << np_corOn << " nm corOn " << nm_corOn << " np0 corOn " << np0_corOn << " nm0 corOn " << nm0_corOn << " n0 corOn " << n0_corOn << endl;
    //cout << "np accOn " << np_accOn/100.0 << " nm accOn " << nm_accOn/100.0 << " np0 accOn " << np0_accOn/100.0 << " nm0 accOn " << nm0_accOn/100.0 << " n0 accOn " << n0_accOn/100.0 << endl;
    cout << "np Avg: " << npAvg << " nm Avg: " << nmAvg << " np0 Avg: " << np0Avg << " nm0 Avg: " << nm0Avg << "\nn0 Avg: " << n0Avg << endl;
  }
  if (xyorz == 0){
    // Measuring IBD Count, average and sttdev - independent of eqn: // Unsigned is always positive
    rX_pos = ((double)npAvg/np0Avg), rX_neg = ((double)nmAvg/nm0Avg);
    cout << "r+: " << rX_pos << " r-: " << rX_neg << endl;
    s.averagexyz = d_seg * ((rX_pos - rX_neg)/(1.0 + rX_pos + rX_neg));
    double np_unbiased = (rX_pos * n0Avg), nm_unbiased = (rX_neg * n0Avg);
    s.totalIBDs = np_unbiased + nm_unbiased + n0Avg;
    n_prime = npAvg + nmAvg + n0Avg;
#ifdef use_complex_errors
    /* L - Left, R - Right, C - Compliment
    // N0 = N0L_C + N0L
    // XL = N0 * (NL / N0L)

    // Original Study: deltaSquared_XL = NL
    // Accurate Detailed Version: deltaSquared_XL = (NL*(N0/N0L)^2) + (N0L*((NL*N0L - N0*NL)/N0L^2)^2) + (N0_C*((NL/N0L)^2))
    delta_XL_squared = ((nmAvg * pow((n0_center/nm0Avg),2)) + (nm0Avg * pow((((nmAvg * nm0Avg) - (n0_center * nmAvg))/ pow(nm0Avg,2)),2)) + ((n0_center - nm0Avg)*pow((nmAvg/nm0Avg),2)));
    delta_XR_squared = ((npAvg * pow((n0_center/np0Avg),2)) + (np0Avg * pow((((npAvg * np0Avg) - (n0_center * npAvg))/ pow(np0Avg,2)),2)) + ((n0_center - np0Avg)*pow((npAvg/np0Avg),2)));
    delta_X0_squared = (n0_center); */
    delta_XL_squared = pow((d_seg + s.averagexyz),2) * ((nmAvg * pow((n0_center/nm0Avg),2)) + (nm0Avg * pow((((nmAvg * nm0Avg) - (n0_center * nmAvg))/ pow(nm0Avg,2)),2)) + ((n0_center - nm0Avg)*pow((nmAvg/nm0Avg),2)));
    delta_XR_squared = pow((d_seg - s.averagexyz),2) * ((npAvg * pow((n0_center/np0Avg),2)) + (np0Avg * pow((((npAvg * np0Avg) - (n0_center * npAvg))/ pow(np0Avg,2)),2)) + ((n0_center - np0Avg)*pow((npAvg/np0Avg),2)));
    delta_X0_squared = pow(s.averagexyz,2) * (n0_center);
    deltaXSquaredArray[0] = delta_XL_squared; deltaXSquaredArray[1] = delta_X0_squared; deltaXSquaredArray[2] = delta_XR_squared;
    cout << "Detailed Delta Values: X-^2 " << delta_XL_squared <<  " X0^2 " << delta_X0_squared << " X+^2 " << delta_XR_squared << endl;
    // n_prime * sigma^2 = NA0[avg^2 + (r+)*(D - avg)^2 +  (r-)*(D + avg)^2]
    s.stddev = sqrt((delta_XL_squared + delta_XR_squared + delta_X0_squared)/n_prime); // Normalized by total count

#else
    delta_XL_squared = nm_unbiased * pow((d_seg + s.averagexyz),2);
    delta_XR_squared = np_unbiased * pow((d_seg - s.averagexyz),2);
    delta_X0_squared = n0_center * pow((s.averagexyz),2);
    deltaXSquaredArray[0] = delta_XL_squared; deltaXSquaredArray[1] = delta_X0_squared; deltaXSquaredArray[2] = delta_XR_squared;
    cout << "Original Delta Values: X-^2 " << delta_XL_squared <<  " X0^2 " << delta_X0_squared << " X+^2 " << delta_XR_squared << endl;
    s.stddev = sqrt((delta_XL_squared + delta_XR_squared + delta_X0_squared)/n_prime);
#endif
  }
  else if (xyorz == 1){
    rY_pos = ((double)npAvg/np0Avg), rY_neg = ((double)nmAvg/nm0Avg);
    cout << "r+: " << rY_pos << " r-: " << rY_neg << endl;
    s.averagexyz = d_seg * ((rY_pos - rY_neg)/(1.0 + rY_pos + rY_neg));
    double np_unbiased = (rY_pos * n0Avg), nm_unbiased = (rY_neg * n0Avg);
    s.totalIBDs = np_unbiased + nm_unbiased + n0Avg;
    n_prime = npAvg + nmAvg + n0Avg;

#ifdef use_complex_errors
    delta_YL_squared = pow((d_seg + s.averagexyz),2) * ((nmAvg * pow((n0_center/nm0Avg),2)) + (nm0Avg * pow((((nmAvg * nm0Avg) - (n0_center * nmAvg))/ pow(nm0Avg,2)),2)) + ((n0_center - nm0Avg)*pow((nmAvg/nm0Avg),2)));
    delta_YR_squared = pow((d_seg - s.averagexyz),2) * ((npAvg * pow((n0_center/np0Avg),2)) + (np0Avg * pow((((npAvg * np0Avg) - (n0_center * npAvg))/ pow(np0Avg,2)),2)) + ((n0_center - np0Avg)*pow((npAvg/np0Avg),2)));
    delta_Y0_squared = pow(s.averagexyz,2) * (n0_center);
    deltaYSquaredArray[0] = delta_YL_squared; deltaYSquaredArray[1] = delta_Y0_squared; deltaYSquaredArray[2] = delta_YR_squared;
    cout << "Detailed Delta Values: Y-^2 " << delta_YL_squared <<  " Y0^2 " << delta_Y0_squared << " Y+^2 " << delta_YR_squared << endl;
    s.stddev = sqrt((delta_YL_squared + delta_YR_squared + delta_Y0_squared)/n_prime);
#else
    delta_YL_squared = nm_unbiased * pow((d_seg + s.averagexyz),2);
    delta_YR_squared = np_unbiased * pow((d_seg - s.averagexyz),2);
    delta_Y0_squared = n0_center * pow((s.averagexyz),2);
    deltaYSquaredArray[0] = delta_YL_squared; deltaYSquaredArray[1] = delta_Y0_squared; deltaYSquaredArray[2] = delta_YR_squared;
    cout << "Original Delta Values: Y-^2 " << delta_YL_squared <<  " Y0^2 " << delta_Y0_squared << " Y+^2 " << delta_YR_squared << endl;
    s.stddev = sqrt((delta_YL_squared + delta_YR_squared + delta_Y0_squared)/n_prime);
#endif
  }
#endif
  cout << "IBD Count: " << s.totalIBDs << endl;
  //cout << "Total IBDs Count: " << s.totalIBDs << "\n" << endl;
  // Done by Poisson Distribution N^2/(sqrt(N)^2) = N; Eff. counts = counts^2/counts_err^2
  /*Since we’re doing a background subtraction, those errors are going to add rather than subtract,
  the error in counts is higher than just sqrt(N). The 19k is the counts squared over the error squared
  which you can effectively use as your N value when doing the error analysis of your theta and phi values*/
  return s;
} // End Program
//-------------------------------------------------------------------------------------------------------------------------------------------

// Angle Calculation
int main(){
  // Ignore Warnings
  gErrorIgnoreLevel=kError;
  s_nd snd;
  Struct result;
  TTree *Tnd = new TTree("Tnd","Neutrino Directionality Tree ");
  TFile *fout = new TFile("OutputDataSeg.root","RECREATE");

  setupNDTreeFill(*Tnd,snd);
  char file_path[300], file_list[400], X_title[300], Y_title[300], Z_title[300];
  string root_file = "AD1_IBD_2020.root";
  sprintf(file_path,"home/prospect-collab/converted_data/Analyzed/Analyzed_2020A_IBD_v23.1");
  sprintf(file_list,"/home/mmo58/prospect_bundle/PROSPECT2x_Analysis/Analysis/AnalyzerConfig/2019B_GoodRuns_RxStatus.txt");

  double sigma = 0, effIBDs = 0, IBD_Count_X = 0, IBD_Count_Y = 0, IBD_Count_Z = 0;
  TH1F *hDataXDiff = new TH1F("hDataXDiff", "", h_bins, -h_max, h_max);
  TH1F *hDataYDiff = new TH1F("hDataYDiff", "", h_bins, -h_max, h_max);
  TH1F *hDataZDiff = new TH1F("hDataZDiff", "", h_bins, -h_max, h_max);

  // For X,Y and Z
  for (int i=0;i<3;i++){
    result = FindDirectionalityValues(file_list, root_file.c_str(), file_path, i);
    if (result.xyz == 0){
      // The zero'th bin is an underflow bin and is not on the visible plot. Root starts the first bin in a histogram at index=1.
      //for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinContent(j+1, ((result.corOnLTFArray[j] - result.accOnLTFArray[j]/100.0) - (result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]/100.0))/result.totalIBDs);
      //for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinError(j+1, (sqrt(result.corOnErrArray[j] + result.accOnErrArray[j]/10000.0 + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffErrArray[j] + result.accOffErrArray[j]/10000.0))/result.totalIBDs));
      snd.momX = result.averagexyz;
      snd.sigmaX = result.stddev; // Errors added in quadature
      cout << "Sigma X: " << snd.sigmaX;
      IBD_Count_X = result.totalIBDs;
    }
    else if (result.xyz == 1){
      //for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinContent(j+1, ((result.corOnLTFArray[j] - result.accOnLTFArray[j]/100.0) - (result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]/100.0))/result.totalIBDs);
      //for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinError(j+1, (sqrt(result.corOnErrArray[j] + result.accOnErrArray[j]/10000.0 + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffErrArray[j] + result.accOffErrArray[j]/10000.0))/result.totalIBDs));
      snd.momY = result.averagexyz;
      snd.sigmaY = result.stddev;
      cout << "Sigma Y: " << snd.sigmaY;
      IBD_Count_Y = result.totalIBDs;
    }
    else if (result.xyz == 2){
      for (int j = 0; j < h_bins; j++) hDataZDiff->SetBinContent(j+1, ((result.corOnLTFArray[j] - result.accOnLTFArray[j]/100.0) - (result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]/100.0)));
      for (int j = 0; j < h_bins; j++) hDataZDiff->SetBinError(j+1, (sqrt(result.corOnErrArray[j] + result.accOnErrArray[j]/10000.0 + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffErrArray[j] + result.accOffErrArray[j]/10000.0))));
      snd.momZ = hDataZDiff->GetMean(1);
      snd.sigmaZ = hDataZDiff->GetStdDev(1);
      cout << "Sigma Z: " << snd.sigmaZ;
      IBD_Count_Z = result.totalIBDs;
      sprintf(Z_title, "Z Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
      hDataZDiff->SetTitle(Z_title); hDataZDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataZDiff->GetYaxis()->SetTitle("Total IBD Events");
    }
  }
  // Manipulate Bin Content
  for (int k = 0; k < h_bins; k++){
    // Skips 0 b/c that is the underflow bin so [6] = -145, k = 5
    if (k == 5) {dataXArray[k] = (rX_neg * n0_center); dataYArray[k] = (rY_neg * n0_center); ErrXArray[k] = sqrt(deltaXSquaredArray[0]); ErrYArray[k] = sqrt(deltaYSquaredArray[0]);}
    //else if (k == 150) {dataXArray[k] = ((rY_neg + 1.0 + rY_pos) * n0_center); dataYArray[k] = ((rX_neg + 1.0 + rX_pos) * n0_center); ErrXArray[k] = sqrt(deltaXSquaredArray[1]); ErrYArray[k] = sqrt(deltaYSquaredArray[1]);}
    else if (k == 150) {dataXArray[k] = (n0_center); dataYArray[k] = (n0_center); ErrXArray[k] = sqrt(deltaXSquaredArray[1]); ErrYArray[k] = sqrt(deltaYSquaredArray[1]);}
    else if (k == 295) {dataXArray[k] = (rX_pos * n0_center); dataYArray[k] = (rY_pos * n0_center); ErrXArray[k] = sqrt(deltaXSquaredArray[2]); ErrYArray[k] = sqrt(deltaYSquaredArray[2]);}
    else {dataXArray[k] = 0.; dataYArray[k] = 0.;}
  }
  for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinContent(j+1, dataXArray[j]);
  for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinError(j+1, ErrXArray[j]);
  sprintf(X_title, "X Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
  hDataXDiff->SetTitle(X_title);
  hDataXDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataXDiff->GetYaxis()->SetTitle("Total IBD Events");
  for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinContent(j+1, dataYArray[j]);
  for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinError(j+1, ErrYArray[j]);
  sprintf(Y_title, "Y Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
  hDataYDiff->SetTitle(Y_title);
  hDataYDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataYDiff->GetYaxis()->SetTitle("Total IBD Events");

  // ADD STUFF TO TTREE
  snd.effIBDs = result.effIBDs;
  snd.livetimeOff = result.livetimeOff;
  snd.livetimeOn = result.livetimeOn;
  snd.totalIBDs = result.totalIBDs;
  snd.totalIBDsErr = result.totalIBDsErr;

  // Largest Sigma value that best accounts for statistical error (usually in z direction). (P in equation).
  if ((snd.sigmaX > snd.sigmaY) && (snd.sigmaX > snd.sigmaZ)) sigma = snd.sigmaX;
  else if ((snd.sigmaY > snd.sigmaX) && (snd.sigmaY > snd.sigmaZ)) sigma = snd.sigmaY;
  else if ((snd.sigmaZ > snd.sigmaX) && (snd.sigmaZ > snd.sigmaY)) sigma = snd.sigmaZ;

  // Angle Phi with respect to vertical between detector and reactor in xy plane
#ifdef use_atan2
  double phiRad = atan2(snd.momY, snd.momX);
#else
  double phiRad = pi/2.0 - atan(snd.momY/snd.momX);
#endif

  double phiDeg = phiRad * 180.0/pi;
#ifdef use_complex_errors
  double phiErrTerm1 = ((snd.sigmaX*snd.momY)/(pow(snd.momX,2)*sqrt(IBD_Count_X)));
  double phiErrTerm2 = snd.sigmaY/(snd.momX*sqrt(IBD_Count_Y));
  double phiErr = sqrt(pow(phiErrTerm1,2) + pow(phiErrTerm2,2));
#else
  double phiErr = sqrt( pow( (sigma * snd.momY)/(sqrt(snd.effIBDs) * pow(snd.momX, 2) ), 2) + pow((sigma)/(sqrt(snd.effIBDs) * snd.momX), 2 ) );  // The error with arctan(x) with x being tan(phi) or tan(theta) propagates as dx/(1+x^2)
#endif

  // General Uncertainty Equation: delta tan(phi) / tan(phi).
  double phiRadErr = phiErr/(1 + pow((pi/2.0 - phiRad), 2));
  double phiDegErr = phiRadErr * 180.0/pi;

  // Angle Theta with respect to horizontal between detector and reactor in xz plane
#ifdef use_atan2
  double thetaRad = atan2(snd.momZ, (sqrt(pow(snd.momX, 2) + pow(snd.momY, 2))));
#else
  double thetaRad = atan(snd.momZ/(sqrt(pow(snd.momX, 2) + pow(snd.momY, 2))));
#endif
  double thetaDeg = thetaRad * 180.0/pi;
#ifdef use_complex_errors
  double thetaErr_r = pow(snd.momX,2) + pow(snd.momY,2);
  double thetaErrTerm1 = pow((snd.sigmaZ/sqrt(IBD_Count_Z)),2);
  double thetaErrTerm2 = pow((snd.sigmaY/sqrt(IBD_Count_Y)),2) * pow(((snd.momY*snd.momZ)/thetaErr_r),2);
  double thetaErrTerm3 = pow((snd.sigmaX/sqrt(IBD_Count_X)),2) * pow(((snd.momX*snd.momZ)/thetaErr_r),2);
  double thetaErr = sqrt( pow((1/thetaErr_r),2) * (thetaErrTerm1 + thetaErrTerm2 + thetaErrTerm3) );
#else
  double thetaErr = sqrt( (pow( (sigma)/(sqrt(snd.effIBDs) * sqrt(pow(snd.momX, 2) + pow(snd.momY, 2))) , 2)) * (1 + pow( (snd.momY * snd.momZ)/(pow(snd.momX, 2) + pow(snd.momY, 2) ), 2) + pow( (snd.momX * snd.momZ)/(pow(snd.momX, 2)+ pow(snd.momY, 2) ), 2) ) );
#endif
  double thetaRadErr = thetaErr/(1 + pow(thetaRad, 2));
  double thetaDegErr = thetaRadErr * 180.0/pi;
  snd.phiDeg = phiDeg;
  snd.phiDegErr = phiDegErr;
  snd.thetaDeg = thetaDeg;
  snd.thetaDegErr = thetaDegErr;
  Tnd->Fill();

  cout << "\n\npx: " << snd.momX << " py: " << snd.momY << " pz: " << snd.momZ << endl;
  cout << "Phi in degress: " << phiDeg << endl; cout << "Phi Error in degrees: " << phiDegErr << endl;
  cout << "Theta in degrees: " << thetaDeg << endl; cout << "Theta Error in degrees: " << thetaDegErr << endl;

  gStyle->SetOptStat(1);
  // Set the Errors in X to zero
  gStyle->SetErrorX(0);
  hDataXDiff->SetMarkerStyle(8); hDataXDiff->SetMarkerColor(1); hDataXDiff->SetLineColor(1);
  hDataYDiff->SetMarkerStyle(8); hDataYDiff->SetMarkerColor(2); hDataYDiff->SetLineColor(2);
  hDataZDiff->SetMarkerStyle(8); hDataZDiff->SetMarkerColor(4); hDataZDiff->SetLineColor(4);
  fout->cd();
  hX->Write(); hY->Write();
  hDataXDiff->Write(); hDataYDiff->Write(); hDataZDiff->Write();
  Tnd->Write();
  fout->Close();
  printf("\nProgram Complete.\n\n");
  return 0; // Closes Function
} // End Program