/*/////////////////////
// Jin Oueslati
///////////////////////
UP TO DATE 4/26/2022
For information regarding IBD Selection Cuts refer to PROSPECT2xAnalysis/Analysis/PhysPulse/IBDCutset.hh and IBD Selection Rules.pdf

AccWindow = 100.0 number of accidental windows (matched to OnTime)
AtmScale = 1.0?
IBD Selection = (CorOn - (AccOn^2/AccWindow)) - (AtmScale*(LiveTimeOn/LiveTimeOff)*(CorOff - (AccOff^2/AccWindow)))

From Christian's workCombining AngleFormulae.C with PromptDelayedDataSeparation.C and adding Histograms.
Vectors can be used for all calculations and bins for histograms but I prefer Arrays/Histos for calculations due to the IBD Selection Rules pdf*/

#include <TSystem.h>
#include <stdio.h>
#include "NeutrinoDirectionality.h"
#include "../GeneralHeader.h"
using namespace std;
#define use_atan2
//#define normalize
//#define extra_xf

// Global Constants
const int h_bins = 301, index_offset = 150;
float h_max = 150.5;
double delayed[3] = {0}, prompt[3] = {0};
//---------------------------------------------------------------------------------------
struct DirectionalityValues{
  // Deadtime Corrected Livetime & Reactor Total Runtime
  double xyz=0, livetimeOff=0, livetimeOn=0, totalIBDs=0, totalIBDsErr=0, effIBDs=0;
  // [301] comes fromthe largest separation we’re allowed to have is +/- 144 mm. I took all the histograms to +/- 150 and binned them in 1 mm increments, which is 301 bins including 0.
  // Correlated and Accidentals
  double corOffArray[h_bins]={0}, accOffArray[h_bins]={0}, corOnArray[h_bins]={0}, accOnArray[h_bins]={0};
  double corOffLTFArray[h_bins]={0}, accOffLTFArray[h_bins]={0}, corOnLTFArray[h_bins]={0}, accOnLTFArray[h_bins]={0};
  // Same as Above but I'm having trouble implementing the vectors into the histograms
};
typedef struct DirectionalityValues Struct;
// My Data Version
Struct FindDirectionalityValues(const char *fname, const char *TFilename, const char *release, int xyorz=0){
  Struct s;
  cout.precision(10);
  s.xyz = xyorz;
  // TTree Parameters
  double corOn_rx = 0, corOn_pp = 0;
  double ncapt_dt = 0, Esmear = 0;
  int pseg = 0, nseg = 0, round_index = 0;
  float diff_index = 0;

  long np_orig_corOff = 0, nm_orig_corOff = 0, n0_orig_corOff = 0, np_orig_accOff = 0, nm_orig_accOff = 0, n0_orig_accOff = 0;
  long np_orig_corOn = 0, nm_orig_corOn = 0, n0_orig_corOn = 0, np_orig_accOn = 0, nm_orig_accOn = 0, n0_orig_accOn = 0;
  long n0_diff_corOff = 0, n0_diff_accOff = 0, n0_diff_corOn = 0, n0_diff_accOn = 0;

  TH1F *hData_corOff = new TH1F("hData_corOff", "", h_bins, -h_max, h_max);
  TH1F *hData_accOff = new TH1F("hData_accOff", "", h_bins, -h_max, h_max);
  TH1F *hData_corOn = new TH1F("hData_corOn", "", h_bins, -h_max, h_max);
  TH1F *hData_accOn = new TH1F("hData_accOn", "", h_bins, -h_max, h_max);
  /* Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal.
  Accidentals are IBDs with neutrons capturing over 1 ms after the prompt.
  We want the Correlated IBDs since they’re actually physical, so that’s why we subtract them off.

  accOnX, corOffX, etc. are the vectors that save the dead time factor x for the dead time correction.
  When I use histograms, I fill them with the x factor, which is constant over each file you look at corOffAvg,
  etc. are the variables I use to calculate the average values.
  In histogram terms, it’s the x value times the y value, both of which you need to know where the average is located.*/

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
    st.Remove(0,6);
    // Formating issue within file directory to correct it.
    // 0 refers to reactor off while 1 is reactor on.
    if (st.Contains(" 0")){
      st.ReplaceAll(" 0", "");
      TFile *f=new TFile(st);
      TVectorD *rtOff = ((TVectorD*)f->Get("runtime"));
      TVectorD *promptvOff = ((TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_prompt"));
      TVectorD *delayedvOff = ((TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_delayed"));
      double xRxOff = (rtOff->Max() / (rtOff->Max() - promptvOff->Max() ) ) * (rtOff->Max() / (rtOff->Max() - delayedvOff->Max() ) );

      TTree *Th = (TTree*)f->Get("P2kIBDPlugin/Tibd");
      long nentries = Th->GetEntries();
      //cout << "Entries: " << nentries << endl;
      for (long i = 0; i < nentries; i++){
        // Esmear refers to Prompt Energy and it is within Prospect's boundary conditions.
        // Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal.
        Th->GetEntry(i);
        delayed[xyorz] = Th->GetLeaf("n_xyz")->GetValue(xyorz);
	      prompt[xyorz] = Th->GetLeaf("xyz")->GetValue(xyorz);
        Esmear = Th->GetLeaf("Esmear")->GetValue(0);
        ncapt_dt = Th->GetLeaf("ncapt_dt")->GetValue(0);
        pseg = Th->GetLeaf("maxseg")->GetValue(0);
        nseg = Th->GetLeaf("n_seg")->GetValue(0);
        round_index = round(delayed[xyorz] - prompt[xyorz]);
        diff_index = (delayed[xyorz] - prompt[xyorz]);
        // Correlated Off
        if (Esmear > 0.8 && Esmear < 7.2 && ncapt_dt > pow(10, 3) && ncapt_dt < 120 * pow(10,3)){
          // In the Tibd TTree: n_xyz is the neutron position, and xyz is the prompt position. They are both arrays of dim 3.
          // ->GetValue(index) from the leaf, the index is the index position for the array.
          // index + 151 thing was to align the histogram bin with the array index since you can’t negatively index an array.
          #ifndef extra_xf
              if (diff_index == 0) n0_diff_corOff++;
              if (pseg == nseg) n0_orig_corOff++;
              else if (diff_index > 0) np_orig_corOff++;
              else if (diff_index < 0) nm_orig_corOff++;
              s.corOffLTFArray[round_index + index_offset] += 1;
          #else
              if (diff_index == 0) n0_diff_corOff += xRxOff;
              if (pseg == nseg) n0_orig_corOff += xRxOff;
              else if (diff_index > 0) np_orig_corOff += xRxOff;
              else if (diff_index < 0) nm_orig_corOff += xRxOff;
              s.corOffLTFArray[round_index + index_offset] += xRxOff;
          #endif
          hData_corOff->Fill(diff_index);
        }
        // Accidentals are IBDs with neutrons capturing over 1 ms after the prompt.
        else if (Esmear > 0.8 && Esmear < 7.2 && ncapt_dt > pow(10, 6)){
          #ifndef extra_xf
              if (diff_index == 0) n0_diff_accOff += 1.*xRxOff;
              if (pseg == nseg) n0_orig_accOff += 1.*xRxOff;
              else if (diff_index > 0) np_orig_accOff += 1.*xRxOff;
              else if (diff_index < 0) nm_orig_accOff += 1.*xRxOff;
                    s.accOffLTFArray[round_index + index_offset] += 1.*xRxOff;
          #else
              if (diff_index == 0) n0_diff_accOff += 1.*(xRxOff*xRxOff);
              if (pseg == nseg) n0_orig_accOff += 1.*(xRxOff*xRxOff);
              else if (diff_index > 0) np_orig_accOff += 1.*(xRxOff*xRxOff);
              else if (diff_index < 0) nm_orig_accOff += 1.*(xRxOff*xRxOff);
                    s.accOffLTFArray[round_index + index_offset] += 1.*(xRxOff*xRxOff);
          #endif
          hData_accOff->Fill(diff_index);
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
      double xRxOn = (rtOn->Max() / (rtOn->Max() - promptvOn->Max() ) ) * (rtOn->Max() / (rtOn->Max() - delayedvOn->Max() ) );
      TTree *Th = (TTree*)f->Get("P2kIBDPlugin/Tibd");
      long nentries = Th->GetEntries();
      for (long i = 0; i < nentries; i++){
        Th->GetEntry(i);
        delayed[xyorz] = Th->GetLeaf("n_xyz")->GetValue(xyorz);
        prompt[xyorz] = Th->GetLeaf("xyz")->GetValue(xyorz);
        Esmear = Th->GetLeaf("Esmear")->GetValue(0);
        ncapt_dt = Th->GetLeaf("ncapt_dt")->GetValue(0);
        pseg = Th->GetLeaf("maxseg")->GetValue(0);
        nseg = Th->GetLeaf("n_seg")->GetValue(0);
        round_index = round(delayed[xyorz] - prompt[xyorz]);
        diff_index = (delayed[xyorz] - prompt[xyorz]);
        // Correlated On
        if (Esmear > 0.8 && Esmear < 7.2 && ncapt_dt > pow(10, 3) && ncapt_dt < 120 * pow(10,3)){
          #ifndef extra_xf
              if (diff_index == 0) n0_diff_corOn++;
              if (pseg == nseg) n0_orig_corOn++;
              else if (diff_index > 0) np_orig_corOn++;
              else if (diff_index < 0) nm_orig_corOn++;
              s.corOnLTFArray[round_index + index_offset] += 1;
          #else
              if (diff_index == 0) n0_diff_corOn += 1.*xRxOn;
              if (pseg == nseg) n0_orig_corOn += 1.*xRxOn;
              else if (diff_index > 0) np_orig_corOn += 1.*xRxOn;
              else if (diff_index < 0) nm_orig_corOn += 1.*xRxOn;
                    s.corOnLTFArray[round_index + index_offset] += 1.*xRxOn;
          #endif
          hData_corOn->Fill(diff_index);
        }
        // Accidental On
        else if (Esmear > 0.8 && Esmear < 7.2 && ncapt_dt > pow(10, 6)){
#ifndef extra_xf
	  if (diff_index == 0) n0_diff_accOn += 1.*xRxOn;
	  if (pseg == nseg) n0_orig_accOn += 1.*xRxOn;
	  else if (diff_index > 0) np_orig_accOn += 1.*xRxOn;
	  else if (diff_index < 0) nm_orig_accOn += 1.*xRxOn;
          s.accOnLTFArray[round_index + index_offset] += 1.*xRxOn;
#else
	  if (diff_index == 0) n0_diff_accOn += 1.*(xRxOn*xRxOn);
	  if (pseg == nseg) n0_orig_accOn += 1.*(xRxOn*xRxOn);
	  else if (diff_index > 0) np_orig_accOn += 1.*(xRxOn*xRxOn);
	  else if (diff_index < 0) nm_orig_accOn += 1.*(xRxOn*xRxOn);
          s.accOnLTFArray[round_index + index_offset] += 1.*(xRxOn*xRxOn);
#endif
          hData_corOn->Fill(diff_index);
        }
      }
      file.peek();
      f->Close();
      s.livetimeOn += (rtOn->Max())/(xRxOn);
    }
  }
  // Comment this loop out if you are not doing it via the histogram method
  for (int k = 0; k < h_bins; k++){
    s.corOffArray[k] = hData_corOff->GetBinContent(k+1);
    s.accOffArray[k] = hData_accOff->GetBinContent(k+1);
    s.corOnArray[k] = hData_corOn->GetBinContent(k+1);
    s.accOnArray[k] = hData_accOn->GetBinContent(k+1);
  }
  // Additional Factors to calculate total IBD Count
  double totalIBDscorOff = 0, totalIBDsaccOff = 0, totalIBDscorOn = 0, totalIBDsaccOn = 0;
  for (unsigned int i = 0; i < h_bins; i++){totalIBDscorOff += s.corOffLTFArray[i]; totalIBDsaccOff += s.accOffLTFArray[i]; totalIBDscorOn += s.corOnLTFArray[i]; totalIBDsaccOn += s.accOnLTFArray[i];}

  // Should be correct if the 10^4us factor is right for the time window in err, but 100 for actual count (corrected in LTF Vectors)???
  s.totalIBDs = (totalIBDscorOn - (totalIBDsaccOn/100.0)) - (s.livetimeOn/s.livetimeOff) * (totalIBDscorOff - (totalIBDsaccOff/100.0));
  s.totalIBDsErr = sqrt(totalIBDscorOn + (totalIBDsaccOn/100000.0) + pow((s.livetimeOn/s.livetimeOff), 2) * (totalIBDscorOff + (totalIBDsaccOff/100000.0)));
  // Done by Poisson Distribution N^2/(sqrt(N)^2) = N; Eff. counts = counts^2/counts_err^2
  s.effIBDs = pow(s.totalIBDs, 2) / pow(s.totalIBDsErr, 2);

  // New
  double npOrigAvg = ((np_orig_corOn - np_orig_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(np_orig_corOff - np_orig_accOff/100.0));
  double nmOrigAvg = ((nm_orig_corOn - nm_orig_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(nm_orig_corOff - nm_orig_accOff/100.0));
  double n0OrigAvg = ((n0_orig_corOn - n0_orig_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(n0_orig_corOff - n0_orig_accOff/100.0));
  double n0DiffAvg = ((n0_diff_corOn - n0_diff_accOn/100.0) - (s.livetimeOn/s.livetimeOff)*(n0_diff_corOff - n0_diff_accOff/100.0));

  // reset histograms h->Reset("ICESM");
  cout << "\n0 for x, 1 for y, 2 for z: " << xyorz << endl;
  //cout << "np Avg: " << npOrigAvg << " nm Avg: " << nmOrigAvg << " n0 Avg (pseg == nseg): " << n0OrigAvg << " n0 Avg (Diff == 0): " << n0DiffAvg << "\n" << endl;
  cout << "Total IBDs Count: " << s.totalIBDs << "\nIBD Error Count: " << s.totalIBDsErr << "\nEff IBDs = IBDs^2/Error^2: " << round(s.effIBDs) << endl;
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
  setupNDTreeFill(*Tnd,snd);
  TFile *fout = new TFile("OutputDirectionality.root","RECREATE");
  string root_file = "AD1_IBD_2020.root";
  char file_path[400], file_list[400];
  sprintf(file_path,"home/prospect-collab/converted_data/Analyzed/Analyzed_2020A_IBD_v23.1");
  sprintf(file_list,"/home/mmo58/prospect_bundle/PROSPECT2x_Analysis/Analysis/AnalyzerConfig/2019B_GoodRuns_RxStatus.txt");

  double sigma = 0, effIBDs = 0.;
  //double IBDsX = 0, IBDsY = 0, IBDsZ = 0;
  TH1F *hDataXDiff = new TH1F("hDataXDiff", "", h_bins, -h_max, h_max);
  TH1F *hDataYDiff = new TH1F("hDataYDiff", "", h_bins, -h_max, h_max);
  TH1F *hDataZDiff = new TH1F("hDataZDiff", "", h_bins, -h_max, h_max);
  char X_title[300], Y_title[300], Z_title[300];
  // For X,Y and Z
  for (int i=0;i<3;i++){
    // Aggragete - all in one go
    result = FindDirectionalityValues(file_list, root_file.c_str(), file_path, i);
    // Could make these into Arrays if more data was added.
    if (result.xyz == 0){
      // The zero'th bin is an underflow bin and is not on the visible plot. Root starts the first bin in a histogram at index=1.
      // Changing j to j+1 just changes mean and StdDev values
#ifdef normalize
      for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - result.accOnLTFArray[j]) - (result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]))/result.totalIBDs);
      for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + (result.accOnArray[j]/10000.0) + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + (result.accOffArray[j]/10000.0)))/result.totalIBDs));
      sprintf(X_title, "X Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
      hDataXDiff->SetTitle(X_title);
      hDataXDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataXDiff->GetYaxis()->SetTitle("Fraction of IBD Events");
#else
# ifdef acc_scale
      // Currently Used
      for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - result.accOnLTFArray[j]) - ((result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]))));
      for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + result.accOnArray[j]/10000.0 + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + result.accOffArray[j]/10000.0))));
# else
      for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - (result.accOnLTFArray[j]/100.0)) - ((result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - (result.accOffLTFArray[j]/100.0)))));
      for (int j = 0; j < h_bins; j++) hDataXDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + (result.accOnArray[j]/10000.0) + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + (result.accOffArray[j]/10000.0)))));
# endif
      sprintf(X_title, "X Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
      hDataXDiff->SetTitle(X_title);
      hDataXDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataXDiff->GetYaxis()->SetTitle("Total IBD Events");
#endif
      //TH1::GetMean(int axis) returns the mean value along axis
      //TH1::GetStdDev(int axis)  returns the sigma distribution along axis
      snd.momX = hDataXDiff->GetMean(1);
      snd.sigmaX = hDataXDiff->GetStdDev(1);
      // Testing Different methods to calculate IBD Counts: Works
      //for (int j = 0; j < h_bins; j++) IBDsX += ((result.corOnLTFArray[j] - result.accOnLTFArray[j]) - (result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]));
      //for (int j = 0; j < h_bins; j++) IBDsX += hDataXDiff->GetBinContent(j);
    }
    else if (result.xyz == 1){
#ifdef normalize
      for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - result.accOnLTFArray[j]) - (result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]))/result.totalIBDs);
      for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + result.accOnArray[j]/10000.0 + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + result.accOffArray[j]/10000.0))/result.totalIBDs));
      sprintf(Y_title, "Y Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
      hDataYDiff->SetTitle(Y_title);
      hDataYDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataYDiff->GetYaxis()->SetTitle("Fraction of IBD Events");
#else
# ifdef acc_scale
      // Currently Used
      for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - result.accOnLTFArray[j]) - ((result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]))));
      for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + result.accOnArray[j]/10000.0 + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + result.accOffArray[j]/10000.0))));
# else
      for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - (result.accOnLTFArray[j]/100.0)) - ((result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - (result.accOffLTFArray[j]/100.0)))));
      for (int j = 0; j < h_bins; j++) hDataYDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + (result.accOnArray[j]/10000.0) + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + (result.accOffArray[j]/10000.0)))));
# endif
      sprintf(Y_title, "Y Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
      hDataYDiff->SetTitle(Y_title);
      hDataYDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataYDiff->GetYaxis()->SetTitle("Total IBD Events");
#endif
      snd.momY = hDataYDiff->GetMean(1);
      snd.sigmaY = hDataYDiff->GetStdDev(1);
    }
    else if (result.xyz == 2){
#ifdef normalize
      for (int j = 0; j < h_bins; j++) hDataZDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - result.accOnLTFArray[j]) - (result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]))/result.totalIBDs);
      for (int j = 0; j < h_bins; j++) hDataZDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + result.accOnArray[j]/10000.0 + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + result.accOffArray[j]/10000.0))/result.totalIBDs));
      sprintf(Z_title, "Z Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
      hDataZDiff->SetTitle(Z_title);
      hDataZDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataZDiff->GetYaxis()->SetTitle("Fraction of IBD Events");
#else
# ifdef acc_scale
      // Currently Used
      for (int j = 0; j < h_bins; j++) hDataZDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - result.accOnLTFArray[j]) - ((result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]))));
      for (int j = 0; j < h_bins; j++) hDataZDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + result.accOnArray[j]/10000.0 + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + result.accOffArray[j]/10000.0))));
# else
      for (int j = 0; j < h_bins; j++) hDataZDiff->SetBinContent(j+1,((result.corOnLTFArray[j] - (result.accOnLTFArray[j]/100.0)) - ((result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - (result.accOffLTFArray[j]/100.0)))));
      for (int j = 0; j < h_bins; j++) hDataZDiff->SetBinError(j+1, (sqrt(result.corOnArray[j] + (result.accOnArray[j]/10000.0) + pow((result.livetimeOn/result.livetimeOff), 2) * (result.corOffArray[j] + (result.accOffArray[j]/10000.0)))));
# endif
      sprintf(Z_title, "Z Separation of Prompt and Delayed Events: IBD = %f", round(result.totalIBDs));
      hDataZDiff->SetTitle(Z_title);
      hDataZDiff->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]"); hDataZDiff->GetYaxis()->SetTitle("Total IBD Events");
#endif
      snd.momZ = hDataZDiff->GetMean(1);
      snd.sigmaZ = hDataZDiff->GetStdDev(1);
      //for (int j = 0; j < h_bins; j++) IBDsZ += ((result.corOnLTFArray[j] - result.accOnLTFArray[j]) - (result.livetimeOn/result.livetimeOff) * (result.corOffLTFArray[j] - result.accOffLTFArray[j]));
    }
  }
  // ADD STUFF TO TTREE
  snd.effIBDs = result.effIBDs;
  snd.livetimeOff = result.livetimeOff; snd.livetimeOn = result.livetimeOn;
  snd.totalIBDs = result.totalIBDs; snd.totalIBDsErr = result.totalIBDsErr;

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
  // The error with arctan(x) with x being tan(phi) or tan(theta) propagates as dx/(1+x^2)
  double phierr = sqrt( pow( (sigma * snd.momY)/(sqrt(snd.effIBDs) * pow(snd.momX, 2) ), 2) + pow((sigma)/(sqrt(snd.effIBDs) * snd.momX), 2 ) );
  // General Uncertainty Equation: delta tan(phi) / tan(phi).
  double phiRadErr = phierr/(1 + pow((pi/2.0 - phiRad), 2));
  double phiDegErr = phiRadErr * 180.0/pi;

  // Angle Theta with respect to horizontal between detector and reactor in xz plane
  //double thetaRad = atan(snd.momZ/(sqrt(pow(snd.momX, 2) + pow(snd.momY, 2))));
#ifdef use_atan2
  double thetaRad = atan2(snd.momZ, (sqrt(pow(snd.momX, 2) + pow(snd.momY, 2))));
#else
  double thetaRad = atan(snd.momZ/(sqrt(pow(snd.momX, 2) + pow(snd.momY, 2))));
#endif

  double thetaDeg = thetaRad * 180.0/pi;
  double thetaerr = sqrt( (pow( (sigma)/(sqrt(snd.effIBDs) * sqrt(pow(snd.momX, 2) + pow(snd.momY, 2))) , 2)) * (1 + pow( (snd.momY * snd.momZ)/(pow(snd.momX, 2) + pow(snd.momY, 2) ), 2) + pow( (snd.momX * snd.momZ)/(pow(snd.momX, 2)+ pow(snd.momY, 2) ), 2) ) );
  double thetaRadErr = thetaerr/(1 + pow(thetaRad, 2));
  double thetaDegErr = thetaRadErr * 180.0/pi;

  snd.phiDeg = phiDeg;
  snd.phiDegErr = phiDegErr;
  snd.thetaDeg = thetaDeg;
  snd.thetaDegErr = thetaDegErr;
  Tnd->Fill();

  cout << "\npx: " << snd.momX << " py: " << snd.momY << " pz: " << snd.momZ;
  cout << "\nLargest sigma value: " << sigma << endl;
  cout << "Phi in degress: " << phiDeg << endl;
  cout << "Phi Error in degrees: " << phiDegErr << endl;
  cout << "Theta in degrees: " << thetaDeg << endl;
  cout << "Theta Error in degrees: " << thetaDegErr << endl;

  gStyle->SetOptStat(1);
  // Set the Errors in X to zero
  gStyle->SetErrorX(0);
  hDataXDiff->SetMarkerStyle(8); hDataXDiff->SetMarkerColor(1); hDataXDiff->SetLineColor(1);
  hDataYDiff->SetMarkerStyle(8); hDataYDiff->SetMarkerColor(2); hDataYDiff->SetLineColor(2);
  hDataZDiff->SetMarkerStyle(8); hDataZDiff->SetMarkerColor(4); hDataZDiff->SetLineColor(4);
  fout->cd();
  hDataXDiff->Write(); hDataYDiff->Write(); hDataZDiff->Write();
  Tnd->Write();
  fout->Close();
  printf("\nProgram Complete.\n\n");
  return 0; // Closes Function
} // End Program