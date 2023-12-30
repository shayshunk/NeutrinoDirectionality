
//This code is an adaptation of Manjinder Oueslati's and Christian Nave's codes

/*
For information regarding IBD Selection Cuts refer to PROSPECT2xAnalysis/Analysis/PhysPulse/IBDCutset.hh and IBD Selection Rules.pdf

IBD selection is defined at https://docdb.wlab.yale.edu/prospect/docs/0028/002802/006/Neutron_Mobility_Technote.pdf (DocDB 2802-v6)
IBD Selection = (CorrelatedOn*x - AccidentalOn*x^2/AccWindow) - AtmScale*(LiveTimeOn/LiveTimeOff)*(CorrelatedOff*x - AccidentalOff*x^2/AccWindow)
  Where:
  . x = runtime/(runtime - prompt veto deadtime) * runtime/(runtime - delayed veto deadtime) //deadtime correction factor for each file
  . AccWindow = 100.0 number of accidental windows (matched to OnTime) //accidental scaling factor
  . AtmScale = 1.0? //atmospheric scaling

-> Realistic Simulation & Perfect Simulation:
   . In this simulation there are no backgrounds (Off) and no Accidentals
   . This code does not simulate backgrounds
   . Use three different calibration files to simulate time drifting across data period
   . Accidentals are not that significant but I will leave the code for them
*/

#include "GeneralHeader.h"
#include "DetectorConfig.h"
#include<vector>
#include<string>

using std::ifstream;
using std::cout;
using std::string;
using std::getline;
using std::vector;

void fillDetectorConfig(vector<vector<int>>& detectorConfig){
    // Filling values based on different periods

    int noSegments = 154;
    int noPeriods = 6;

    for (int i = 0; i < noPeriods; i++) {

        int excludeSize = excludeList[i].size();
        int counter = 0, tmp = 0;
        vector<int> period;

        for (int j = 0; j < noSegments; j++) {
            // Filling live segments by checking exclude list
            tmp = excludeList[i][counter];

            if (j == tmp) 
            {
                period.push_back(0);
                if (counter + 1 < excludeSize) 
                {
                    counter += 1;
                }
            }
            else{
                period.push_back(1);
            }
        }

        detectorConfig.push_back(period);
    }

    cout << "Below is the detector configuration.\n";

    for (int i = 0; i < detectorConfig.size(); i++) 
    {
        cout << "Detector configuration for period " << i + 1 << "\n";
        for (int j = 0; j < detectorConfig[i].size(); j++) 
        {   
            cout << detectorConfig[i][j] << " ";
            if ((j + 1) % 14 == 0)
                cout << "\n";
        }
        cout << "\n";
    }
}

bool checkNeighbor(vector<vector<int>>& detectorConfig, int periodNo, int segNo, char dir)
{
    // Used for dead segment calculations

    bool neighbor = false;

    periodNo = periodNo - 1;

    switch(dir){
        case 'r':
           neighbor = detectorConfig[periodNo][segNo + 1];
           break;
        case 'l':
           neighbor = detectorConfig[periodNo][segNo - 1];
           break;
        case 'u':
           neighbor = detectorConfig[periodNo][segNo + 14];
           break;
        case 'd':
           neighbor = detectorConfig[periodNo][segNo - 14];
           break;
        default:
            cout << "Segment number: " << segNo << " has no live neighbor in direction " << dir << "!\n";
            return false;
    }

    return neighbor;

}


int main()
{
    // Ignore Warnings
    gErrorIgnoreLevel = kError;

    int periodNo = 2;
    string period = std::to_string(periodNo);

    vector<vector<int>> detectorConfig;
    fill();
    fillDetectorConfig(detectorConfig);

    cout << checkNeighbor(detectorConfig, 1, 39, 'r') << "\n";
    cout << checkNeighbor(detectorConfig, 1, 43, 'r') << "\n";
    cout << checkNeighbor(detectorConfig, 1, 91, 'r') << "\n";

    vector<string> vars;
    vars.push_back("X");
    vars.push_back("Y");
    vars.push_back("Z");

    // Save stuff into a root file
    TFile *f_output = new TFile("OutputDirectionality.root", "recreate");

    // Angles
    double phiDeg_data, phiDegErr_data, thetaDeg_data, thetaDegErr_data; 
    double phiDeg_dataBias, phiDegErr_dataBias, thetaDeg_dataBias, thetaDegErr_dataBias; 
    double phiDeg_realsim, phiDegErr_realsim, thetaDeg_realsim, thetaDegErr_realsim;
    double phiDeg_realsimBias, phiDegErr_realsimBias, thetaDeg_realsimBias, thetaDegErr_realsimBias; 
    double phiDeg_perfsim, phiDegErr_perfsim, thetaDeg_perfsim, thetaDegErr_perfsim; 
    double phiDeg_true, phiDegErr_true, thetaDeg_true, thetaDegErr_true; 

    // Make a histogram for each variable
    vector<TH1D*> h_data_CorOn, h_data_CorOff, h_data_AccOn, h_data_AccOff, h_data_Diff;
    vector<TH1D*> h_data_CorOnBias, h_data_CorOffBias, h_data_AccOnBias, h_data_AccOffBias, h_data_DiffBias;
    vector<TH1D*> h_realsim_Cor, h_realsim_Acc, h_realsim_Diff;
    vector<TH1D*> h_realsim_CorBias, h_realsim_AccBias, h_realsim_DiffBias;
    vector<TH1D*> h_perfsim_Cor, h_perfsim_Acc, h_perfsim_Diff;

    // [300] comes from the largest separation we're allowed to have is +/- 144 mm. I took all the histograms to +/- 150 and binned them in 1 mm increments, which is 300 bins
    //int h_bins = 300; 
    //double h_max = 150.0;
    // [301] comes from the largest separation we’re allowed to have is +/- 144 mm. I took all the histograms to +/- 150 and binned them in 1 mm increments, which is 301 bins including 0
    int h_bins = 301; 
    double h_max = 150;

    for (vector<string>::iterator var = vars.begin(); var != vars.end(); ++var)
    {

        h_data_CorOn.push_back ( new TH1D(Form("Data_CorrIBD_On_%s",  var->c_str()), "Data", h_bins, -h_max, h_max) );
        h_data_CorOff.push_back( new TH1D(Form("Data_CorrIBD_Off_%s", var->c_str()), "Data", h_bins, -h_max, h_max) );
        h_data_AccOn.push_back ( new TH1D(Form("Data_AccIBD_On_%s",   var->c_str()), "Data", h_bins, -h_max, h_max) );
        h_data_AccOff.push_back( new TH1D(Form("Data_AccIBD_Off_%s",  var->c_str()), "Data", h_bins, -h_max, h_max) );
        h_data_Diff.push_back  ( new TH1D(Form("Data_Diff_%s",        var->c_str()), "Data", h_bins, -h_max, h_max) );

        h_data_CorOnBias.push_back ( new TH1D(Form("Data_CorrIBD_On_Bias%s",  var->c_str()), "UnbiasedData", h_bins, -h_max, h_max) );
        h_data_CorOffBias.push_back( new TH1D(Form("Data_CorrIBD_Off_Bias%s", var->c_str()), "UnbiasedData", h_bins, -h_max, h_max) );
        h_data_AccOnBias.push_back ( new TH1D(Form("Data_AccIBD_On_Bias%s",   var->c_str()), "UnbiasedData", h_bins, -h_max, h_max) );
        h_data_AccOffBias.push_back( new TH1D(Form("Data_AccIBD_Off_Bias%s",  var->c_str()), "UnbiasedData", h_bins, -h_max, h_max) );
        h_data_DiffBias.push_back  ( new TH1D(Form("Data_Diff_Bias%s",        var->c_str()), "UnbiasedData", h_bins, -h_max, h_max) );

        h_realsim_Cor.push_back ( new TH1D(Form("RealSim_CorrIBD_%s", var->c_str()), "RealSim", h_bins, -h_max, h_max) );
        h_realsim_Acc.push_back ( new TH1D(Form("RealSim_AccIBD_%s",  var->c_str()), "RealSim", h_bins, -h_max, h_max) );
        h_realsim_Diff.push_back( new TH1D(Form("RealSim_Diff_%s",    var->c_str()), "RealSim", h_bins, -h_max, h_max) );

        h_realsim_CorBias.push_back ( new TH1D(Form("RealSim_CorrIBD_Bias%s", var->c_str()), "UnbiasedRealSim", h_bins, -h_max, h_max) );
        h_realsim_AccBias.push_back ( new TH1D(Form("RealSim_AccIBD_Bias%s",  var->c_str()), "UnbiasedRealSim", h_bins, -h_max, h_max) );
        h_realsim_DiffBias.push_back( new TH1D(Form("RealSim_Diff_Bias%s",    var->c_str()), "UnbiasedRealSim", h_bins, -h_max, h_max) );

        h_perfsim_Cor.push_back ( new TH1D(Form("PerfSim_CorrIBD_%s", var->c_str()), "PerfSim", h_bins, -h_max, h_max) );
        h_perfsim_Acc.push_back ( new TH1D(Form("PerfSim_AccIBD_%s",  var->c_str()), "PerfSim", h_bins, -h_max, h_max) );
        h_perfsim_Diff.push_back( new TH1D(Form("PerfSim_Diff_%s",    var->c_str()), "PerfSim", h_bins, -h_max, h_max) );
    }

    //----------------------------------------------------
    // Fill the data: histograms & angles
    //----------------------------------------------------
    {
        int countlines = 0;
        bool ReactorOn = true;
        double livetimeOff = 0.0; //Total livetime for all RxOff runs or reactor off total runtime 
        double livetimeOn = 0.0;  //Total livetime for all RxOn runs or reactor on total runtime

        for (int c = 1; c <= 5; c++)
        {
            periodNo = c;
            period = std::to_string(c);

             //string file_list = "/home/al2667/PROSPECT2x_Analysis/Analysis/AnalyzerConfig/2019B_GoodRuns_RxStatus.txt";
            string file_list = Form("/project/prospect/tmp/Analyzed/SEER_DS_period_%s/Period_%s_files.txt", period.c_str(), period.c_str());
            // Opening File
            ifstream file;
            file.open(file_list, ifstream::in);
            if( !(file.is_open() && file.good()) ){
                cout << "Good runs file not found. Exiting\n";
                
                return -1;
            }   

            while (file.good() && !file.eof())
            {

                    countlines = countlines + 1;
                    if( countlines % 100 == 0 )
                        cout << "Looking at Data file: " << countlines << "/4000.\n";

                    string line;
                    getline(file, line);
                    //TString st = Form("/project/prospect/tmp/prospect-collab/converted_data/Analyzed/Analyzed_2020A_IBD_v23.1/%s/AD1_IBD_2020.root", line.data()); // Yale Meitner server
                    TString st = Form("/project/prospect/tmp/Analyzed/SEER_DS_period_%s/%s/AD1_IBD_2022_DS_SEER_fid4.root", period.c_str(), line.data()); // Yale Meitner server

                    //cout << "st = " << st << "\n";
                    // Use above for debugging 

                    // Formating issue within file directory to correct it
                    // 0 refers to reactor off while 1 is reactor on
                    if (st.Contains(" 0")) 
                    {
                        st.ReplaceAll(" 0", "");
                        ReactorOn = false;
                    }

                    if (st.Contains(" 1"))
                    {
                        st.ReplaceAll(" 1", "");
                        ReactorOn = true;
                    }

                    TFile *f = new TFile(st);

                    //cout << "st = " << st << "\n";
                    // Use above for debugging

                    TVectorD *rt = (TVectorD*)f->Get("runtime"); // total duration of the run
                    TVectorD *promptv = (TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_prompt"); // prompt veto deadtime
                    TVectorD *delayedv = (TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_delayed"); // delayed veto deadtime
                    double xRx = rt->Max() / (rt->Max() - promptv->Max()) * rt->Max() / (rt->Max() - delayedv->Max()); // deadtime correction coefficient (veto deadtime correction)
                    if( ReactorOn )
                        livetimeOn += rt->Max() / xRx;
                    else
                        livetimeOff += rt->Max() / xRx;
                
                    TTree *Th = (TTree*)f->Get("P2kIBDPlugin/Tibd");
                    long nentries = Th->GetEntries();

                    for (long i = 0; i < nentries; i++)
                    {
                        Th->GetEntry(i);

                        for (unsigned int ivar = 0; ivar != vars.size(); ++ivar)
                        {
                            double Esmear = Th->GetLeaf("Esmear")->GetValue(0); //prompt energy in MeV. It is within PROSPECT's boundary conditions
                            double ncapt_dt = Th->GetLeaf("ncapt_dt")->GetValue(0); //neutron capture occurs a dt after the prompt signal. It is measured in ns
                            int pseg = Th->GetLeaf("maxseg")->GetValue(0); //prompt segment
                            int nseg = Th->GetLeaf("n_seg")->GetValue(0); //neutron segment
                            double delayed = Th->GetLeaf("n_xyz")->GetValue(ivar); //n_xyz is the neutron position in mm. Array of dim 3 (x, y, z)
                            double prompt = Th->GetLeaf("xyz")->GetValue(ivar); //xyz is the prompt position in mm. Array of dim 3 (x, y, z)
                            double diff_index = delayed - prompt;
                            
                            if (pseg == nseg || diff_index > 0 || diff_index < 0)
                            {
                                // Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal
                                if (Esmear > 0.8 && Esmear < 7.4 && ncapt_dt > pow(10, 3) && ncapt_dt < 120 * pow(10, 3))
                                {
                                    if (ReactorOn)
                                    {
                                        h_data_CorOn[ivar]->Fill(diff_index);

                                        if (pseg == nseg)
                                        {
                                            bool posDir = false;
                                            bool negDir = false;

                                            if (ivar == 0)
                                            {
                                                posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'r');
                                                negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'l');
                                            }
                                            else if (ivar == 1)
                                            {
                                                posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'u');
                                                negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'd');
                                            }

                                            if (posDir && !negDir)
                                                h_data_CorOnBias[ivar]->Fill(diff_index + 145);
                                            else if (!posDir && negDir)
                                                h_data_CorOnBias[ivar]->Fill(diff_index - 145);
                                            else if (posDir && negDir)
                                                h_data_CorOnBias[ivar]->Fill(diff_index);
                                        }
                                    }
                                    else
                                    {
                                        h_data_CorOff[ivar]->Fill(diff_index);

                                        if (pseg == nseg)
                                        {
                                            bool posDir = false;
                                            bool negDir = false;

                                            if (ivar == 0)
                                            {
                                                posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'r');
                                                negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'l');
                                            }
                                            else if (ivar == 1)
                                            {
                                                posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'u');
                                                negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'd');
                                            }

                                            if (posDir && !negDir)
                                                h_data_CorOffBias[ivar]->Fill(diff_index + 145);
                                            else if (!posDir && negDir)
                                                h_data_CorOffBias[ivar]->Fill(diff_index - 145);
                                            else if (posDir && negDir)
                                                h_data_CorOffBias[ivar]->Fill(diff_index);
                                        }
                                    }
                                }
            
                                // Accidentals are IBDs with neutrons capturing over 1 ms after the prompt
                                if (Esmear > 0.8 && Esmear < 7.4 && ncapt_dt > pow(10, 6))
                                {
                                    if (ReactorOn)
                                    {
                                        h_data_AccOn[ivar]->Fill(diff_index, xRx);

                                        if (pseg == nseg)
                                        {
                                            bool posDir = false;
                                            bool negDir = false;

                                            if (ivar == 0)
                                            {
                                                posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'r');
                                                negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'l');
                                            }
                                            else if (ivar == 1)
                                            {
                                                posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'u');
                                                negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'd');
                                            }

                                            if (posDir && !negDir)
                                                h_data_AccOnBias[ivar]->Fill(diff_index + 145, xRx);
                                            else if (!posDir && negDir)
                                                h_data_AccOnBias[ivar]->Fill(diff_index - 145, xRx);
                                            else if (posDir && negDir)
                                                h_data_AccOnBias[ivar]->Fill(diff_index, xRx);
                                        }
                                    }
                                    else
                                    {
                                        h_data_AccOff[ivar]->Fill(diff_index, xRx);

                                        if (pseg == nseg)
                                        {
                                            bool posDir = false;
                                            bool negDir = false;

                                            if (ivar == 0)
                                            {
                                                posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'r');
                                                negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'l');
                                            }
                                            else if (ivar == 1)
                                            {
                                                posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'u');
                                                negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'd');
                                            }

                                            if (posDir && !negDir)
                                                h_data_AccOffBias[ivar]->Fill(diff_index + 145, xRx);
                                            else if (!posDir && negDir)
                                                h_data_AccOffBias[ivar]->Fill(diff_index - 145, xRx);
                                            else if (posDir && negDir)
                                                h_data_AccOffBias[ivar]->Fill(diff_index, xRx);
                                        }
                                    }
                                }
                            }
                        }//end loop over variables
                    }//end loop over the root file

                    // Returns the next character in the input sequence, without extracting it: The character is left as the next character to be extracted from the stream
                    file.peek();
                    f->Close();
                }
            }//end of while

       cout << " Total livetime for all RxOff " << livetimeOff << "\n";
       cout << " Total livetime for all RxOn  " << livetimeOn  << "\n";

       vector<double> mean, sigma, meanBias, sigmaBias, effIBD;
       double totalIBDs = 0.0, totalIBDsErr = 0.0, effIBDs = 0.0, atm_scaling = 1.00025443769309;
       for (unsigned int ivar = 0; ivar != vars.size(); ++ivar)
       {
            //IBD events = (Correlated - Accidental/100)_{reactor on} - livetimeOn/livetimeOff*(Correlated - Accidental/100)_{reactor off}
            //IBD events = (Correlated - Accidental/100)_{reactor on} + (-livetimeOn/livetimeOff*Correlated + livetimeOn/livetimeOff*Accidental/100)_{reactor off}
            h_data_Diff[ivar] = (TH1D*)h_data_CorOn[ivar]->Clone("DataDiff");
            h_data_Diff[ivar]->Add(h_data_AccOn[ivar], -1.0/100.0);
            h_data_Diff[ivar]->Add(h_data_CorOff[ivar], -livetimeOn*atm_scaling/livetimeOff);
            h_data_Diff[ivar]->Add(h_data_AccOff[ivar], livetimeOn*atm_scaling/livetimeOff/100.0);

            h_data_DiffBias[ivar] = (TH1D*)h_data_CorOnBias[ivar]->Clone("DataDiffBias");
            h_data_DiffBias[ivar]->Add(h_data_AccOnBias[ivar], -1.0/100.0);
            h_data_DiffBias[ivar]->Add(h_data_CorOffBias[ivar], -livetimeOn*atm_scaling/livetimeOff);
            h_data_DiffBias[ivar]->Add(h_data_AccOffBias[ivar], livetimeOn*atm_scaling/livetimeOff/100.0);
 
            totalIBDs = h_data_Diff[ivar]->IntegralAndError(0, h_data_Diff[ivar]->GetNbinsX() + 1, totalIBDsErr);
            effIBDs = pow(totalIBDs, 2) / pow(totalIBDsErr, 2); //Effective IBD counts. Done by Poisson Distribution N^2/(sqrt(N)^2) = N; Eff. counts = counts^2/counts_err^2
            cout << Form(" Total IBD events %0.0f ± %0.0f. Effective IBD counts %0.0f", totalIBDs, totalIBDsErr, effIBDs) << "\n";

            effIBD.push_back(effIBDs);

            h_data_Diff[ivar]->SetName( Form("Data_Diff_%s", vars[ivar].c_str()) );   
            h_data_Diff[ivar]->SetTitle( Form("%s Separation of Prompt and Delayed Events: IBD = %0.0f", vars[ivar].c_str(), totalIBDs) );
            h_data_Diff[ivar]->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]");
            h_data_Diff[ivar]->GetYaxis()->SetTitle("Total IBD Events");

            mean.push_back(h_data_Diff[ivar]->GetMean());
            sigma.push_back(h_data_Diff[ivar]->GetStdDev());

            if (ivar == 0 || ivar == 1)
            {
                double rPlus, rMinus, px, pxErr, D = 145.7;
                double np, npp, npm, nm, nmm, n0;
                double x1, x2, x3, x4, x5;
                double npErr, nppErr, npmErr, nmErr, nmmErr, n0Err;
                double x1Err, x2Err, x3Err, x4Err, x5Err;
                double rPlusErr, rMinusErr;

                np = h_data_Diff[ivar]->GetBinContent(297);
                npp = h_data_DiffBias[ivar]->GetBinContent(297);
                nm = h_data_Diff[ivar]->GetBinContent(5);
                nmm = h_data_DiffBias[ivar]->GetBinContent(5);
                n0 = h_data_Diff[ivar]->GetBinContent(151);
                npm = h_data_DiffBias[ivar]->GetBinContent(151);

                npErr = h_data_Diff[ivar]->GetBinError(297);
                nppErr = h_data_DiffBias[ivar]->GetBinError(297);
                nmErr = h_data_Diff[ivar]->GetBinError(5);
                nmmErr = h_data_DiffBias[ivar]->GetBinError(5);
                n0Err = h_data_Diff[ivar]->GetBinError(151);
                npmErr = h_data_DiffBias[ivar]->GetBinError(151);
                
                x1 = np;
                x2 = npp + npm;
                x3 = nm;
                x4 = nmm + npm;
                
                rPlus = x1/x2;
                rMinus = x3/x4;

                x1Err = npErr;
                x2Err = sqrt(pow(nppErr, 2) + pow(npmErr, 2));
                x3Err = nmErr;
                x4Err = sqrt(pow(nmmErr, 2) + pow(npmErr, 2));

                rPlusErr = rPlus * sqrt(pow((x2Err/x2), 2) + pow((x1Err/x1), 2));
                rMinusErr = rMinus * sqrt(pow((x4Err/x4), 2) + pow((x3Err/x3), 2));

                px = D * ( (rPlus - rMinus) / (rPlus + rMinus + 1));

                pxErr = D * pow(1 / ((nm * (npm + npp) + (nmm + npm) * (np + npm + npp))), 2) * sqrt(pow((nmm + npm) * (npm + npp), 2) * (pow(npErr * (2 * nm + nmm + npm), 2)
                                                                                                                                        + pow(nmErr * (2 * np + npp + npm), 2))
                                                                                                    + pow((np * (npm + nmm) * (2 * nm + nmm + npm) * nppErr), 2)
                                                                                                    + pow((npmErr * (np * pow((nmm + npm), 2) + nm * (2 * nmm * np - 2 * np * npp - pow((npm + npp), 2)))), 2)
                                                                                                    + pow((nm * (npm + npp) * (2 * np + npm + npp) * nmmErr), 2));

                meanBias.push_back(px);
                sigmaBias.push_back(pxErr);
                
                /*
                cout << "Errors for variable: " << ivar << "\n";
                cout << "rPlus Error: " << rPlusErr << "\n";
                cout << "rMinus Error: " << rMinusErr << "\n";
                cout << "x1 Error: " << x1Err << "\n";
                cout << "x2 Error: " << x2Err << "\n";
                cout << "x3 Error: " << x3Err << "\n";
                cout << "x4 Error: " << x4Err << "\n";
                cout << "n0 Error: " << n0Err << "\n";
                cout << "nPlus Error: " << npErr << "\n";
                cout << "nMinus Error: " << nmErr << "\n";
                cout << "nPlusPlus Error: " << nppErr << "\n";
                cout << "nMinusMinus Error: " << nmmErr << "\n";
                cout << "rPlus: " << rPlus << "\n";
                cout << "rMinus: " << rMinus << "\n";
                cout << "x1: " << x1 << "\n";
                cout << "x2: " << x2 << "\n";
                cout << "x3: " << x3 << "\n";
                cout << "x4: " << x4 << "\n";
                cout << "n0: " << n0 << "\n";
                cout << "nPlus: " << np << "\n";
                cout << "nMinus: " << nm << "\n";
                cout << "nPlusPlus: " << npp << "\n";
                cout << "nMinusMinus: " << nmm << "\n"; */
            }
            else{
                meanBias.push_back(h_data_Diff[ivar]->GetMean());
                sigmaBias.push_back(h_data_Diff[ivar]->GetStdDev());
            }
       }

       // Largest sigma value that best accounts for statistical error (usually in z direction). P in equation
       double lsigma = 0, lsigmaBias = 0;

       if ((sigma[0] > sigma[1]) && (sigma[0] > sigma[2])) 
        lsigma = sigma[0];
       else if ((sigma[1] > sigma[0]) && (sigma[1] > sigma[2])) 
        lsigma = sigma[1];
       else 
        lsigma = sigma[2];
    
       sigmaBias[0] = sigmaBias[0] * sqrt(effIBD[0]);
       sigmaBias[1] = sigmaBias[1] * sqrt(effIBD[0]);

       if ((sigmaBias[0] > sigmaBias[1]) && (sigmaBias[0] > sigmaBias[2])) 
        lsigmaBias = sigmaBias[0];
       else if ((sigmaBias[1] > sigmaBias[0]) && (sigmaBias[1] > sigmaBias[2])) 
        lsigmaBias = sigmaBias[1];
       else 
        lsigmaBias = sigmaBias[2];

       // Angle Phi with respect to vertical between detector and reactor in xy plane
       double phiRad = atan2(mean[1], mean[0]); //phi = arctan(y/x)
       phiDeg_data = phiRad * 180.0/pi;
       double phiRadErr = sqrt( pow(lsigma/sqrt(effIBDs)*mean[1]/(mean[0]*mean[0]), 2) + pow(lsigma/sqrt(effIBDs)/(mean[0]), 2) ) / (1 + pow(mean[1]/mean[0], 2)); 
       phiDegErr_data = phiRadErr * 180.0/pi;

       // Angle Theta with respect to horizontal between detector and reactor in xz plane
       double thetaRad = atan2(mean[2], sqrt(mean[0]*mean[0] + mean[1]*mean[1])); //theta = arctan(z/sqrt(x^2 + y^2))
       thetaDeg_data = thetaRad * 180.0/pi;
       double thetaRadErr = sqrt( pow(lsigma/sqrt(effIBDs)/sqrt(mean[0]*mean[0] + mean[1]*mean[1]), 2) * (pow(mean[0]*mean[2]/(mean[0]*mean[0] + mean[1]*mean[1]), 2) + pow(mean[1]*mean[2]/(mean[0]*mean[0] + mean[1]*mean[1]), 2) + 1) ) / (1 + pow(mean[2]/sqrt(mean[0]*mean[0] + mean[1]*mean[1]), 2));
       thetaDegErr_data = thetaRadErr * 180.0/pi;

       cout << "  Mean x: "  << mean[0]  << ", Mean y: "  << mean[1]  << ", Mean z: "  << mean[2]  << "\n";
       cout << "  Sigma x: " << sigma[0] << ", Sigma y: " << sigma[1] << ", Sigma z: " << sigma[2] << "\n";
       cout << "  Largest sigma value: "    << lsigma           << "\n";
       cout << "  Phi in degrees: "         << phiDeg_data      << "\n";
       cout << "  Phi Error in degrees: "   << phiDegErr_data   << "\n";
       cout << "  Theta in degrees: "       << thetaDeg_data    << "\n";
       cout << "  Theta Error in degrees: " << thetaDegErr_data << "\n\n";

       // Unbiased Angle Calc
       // Angle Phi with respect to vertical between detector and reactor in xy plane
       double phiRadBias = atan2(meanBias[1], meanBias[0]); //phi = arctan(y/x)
       phiDeg_dataBias = phiRadBias * 180.0/pi;
       double phiRadErrBias = sqrt( pow(sigmaBias[0]/sqrt(effIBD[0])*meanBias[1]/(meanBias[0]*meanBias[0]), 2) + pow(sigmaBias[1]/sqrt(effIBD[1])/(meanBias[0]), 2) ) / (1 + pow(meanBias[1]/meanBias[0], 2)); 
       phiDegErr_dataBias = phiRadErrBias * 180.0/pi;

       // Angle Theta with respect to horizontal between detector and reactor in xz plane
       double thetaRadBias = atan2(meanBias[2], sqrt(meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1])); //theta = arctan(z/sqrt(x^2 + y^2))
       thetaDeg_dataBias = thetaRadBias * 180.0/pi;
       double thetaRadErrBias = sqrt( (1 / (meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1])) * ( pow((meanBias[0]*meanBias[2]*sigmaBias[0]/sqrt(effIBD[0]) / (meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1])), 2) + pow((meanBias[1]*meanBias[2]*sigmaBias[1]/sqrt(effIBD[1]) / (meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1])), 2 ) + pow(sigmaBias[2]/sqrt(effIBD[2]), 2))) / (1 + pow(meanBias[2]/sqrt(meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1]), 2));
       thetaDegErr_dataBias = thetaRadErrBias * 180.0/pi;

       cout << "  Mean Unbiased x: "  << meanBias[0]  << ", Mean Unbiased y: "  << meanBias[1]  << ", Mean Unbiased z: "  << meanBias[2]  << "\n";
       cout << "  Sigma Unbiased x: " << sigmaBias[0] << ", Sigma Unbiased y: " << sigmaBias[1] << ", Sigma Unbiased z: " << sigmaBias[2] << "\n";
       cout << "  Largest Unbiased sigma value: "    << lsigmaBias          << "\n";
       cout << "  Unbiased Phi in degrees: "         << phiDeg_dataBias      << "\n";
       cout << "  Unbiased Phi Error in degrees: "   << phiDegErr_dataBias   << "\n";
       cout << "  Unbiased Theta in degrees: "       << thetaDeg_dataBias    << "\n";
       cout << "  Unbiased Theta Error in degrees: " << thetaDegErr_dataBias << "\n\n";
    }
    //----------------------------------------------------

    //----------------------------------------------------
    // Fill the realistic simulation: histograms & angles
    //----------------------------------------------------
    {
       string file_list = "/project/prospect/tmp/mmo58/prospect_bundle/MyWork/NeutrinoDirectionality/Simulations/SimFileGoodRunsList.txt";
       // Opening File
       std::ifstream file;
       file.open(file_list, std::ifstream::in);
       if( !(file.is_open() && file.good()) )
       {
           printf("Good runs file not found. Exiting\n");

           return -1;
       }
      
       int countlines = 0;

       while( file.good() && !file.eof() )
       {
            periodNo = 6;
            countlines = countlines + 1;
            cout << "Looking at Realistic Simulation file " << countlines << " / 20" << "\n";

            string line;
            getline(file, line);
            TString st = Form("/project/prospect/tmp/mmo58/prospect_bundle/MyWork/NeutrinoDirectionality/Simulations/H5Files/%s/Jun_AD1_IBD_2020.root", line.data()); // Yale Meitner server
            //cout << "st = " << st << "\n";
            TFile *f = new TFile(st);

            TVectorD *rt = (TVectorD*)f->Get("runtime"); //total duration of the run
            TVectorD *promptv = (TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_prompt"); //prompt veto deadtime
            TVectorD *delayedv = (TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_delayed"); //delayed veto deadtime
            double xRx = rt->Max() / (rt->Max() - promptv->Max()) * rt->Max() / (rt->Max() - delayedv->Max()); //deadtime correction coefficient (veto deadtime correction)
        
            TTree *Th = (TTree*)f->Get("P2kIBDPlugin/Tibd");
            long nentries = Th->GetEntries();
            for( long i = 0; i < nentries; i++ )
            {
                Th->GetEntry(i);

                for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
                {
                    double Esmear = Th->GetLeaf("Esmear")->GetValue(0); //prompt energy in MeV. It is within PROSPECT's boundary conditions
                    double ncapt_dt = Th->GetLeaf("ncapt_dt")->GetValue(0); //neutron capture occurs a dt after the prompt signal. It is measured in ns
                    int pseg = Th->GetLeaf("maxseg")->GetValue(0); //prompt segment
                    int nseg = Th->GetLeaf("n_seg")->GetValue(0); //neutron segment
                    double delayed = Th->GetLeaf("n_xyz")->GetValue(ivar); //n_xyz is the neutron position in mm. Array of dim 3 (x, y, z)
                    double prompt = Th->GetLeaf("xyz")->GetValue(ivar); //xyz is the prompt position in mm. Array of dim 3 (x, y, z)
                    double diff_index = delayed - prompt;

                    if (pseg == nseg || diff_index > 0 || diff_index < 0)
                    {
                        // Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal
                        if (Esmear > 0.8 && Esmear < 7.4 && ncapt_dt > pow(10, 3) && ncapt_dt < 120 * pow(10, 3))
                        {
                            h_realsim_Cor[ivar]->Fill(diff_index);

                            if (pseg == nseg)
                            {
                                bool posDir = false;
                                bool negDir = false;

                                if (ivar == 0)
                                {
                                    posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'r');
                                    negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'l');
                                }
                                else if (ivar == 1)
                                {
                                    posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'u');
                                    negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'd');
                                }

                                if (posDir && !negDir)
                                    h_realsim_CorBias[ivar]->Fill(diff_index + 145);
                                else if (!posDir && negDir)
                                    h_realsim_CorBias[ivar]->Fill(diff_index - 145);
                                else if (posDir && negDir)
                                    h_realsim_CorBias[ivar]->Fill(diff_index);
                            }
                        }
    
                        // Accidentals are IBDs with neutrons capturing over 1 ms after the prompt
                        if (Esmear > 0.8 && Esmear < 7.4 && ncapt_dt > pow(10, 6))
                        {
                            h_realsim_Acc[ivar]->Fill(diff_index, xRx);

                            if (pseg == nseg)
                            {
                                bool posDir = false;
                                bool negDir = false;

                                if (ivar == 0)
                                {
                                    posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'r');
                                    negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'l');
                                }
                                else if (ivar == 1)
                                {
                                    posDir = checkNeighbor(detectorConfig, periodNo, pseg, 'u');
                                    negDir = checkNeighbor(detectorConfig, periodNo, pseg, 'd');
                                }

                                if (posDir && !negDir)
                                    h_realsim_AccBias[ivar]->Fill(diff_index + 145, xRx);
                                else if (!posDir && negDir)
                                    h_realsim_AccBias[ivar]->Fill(diff_index - 145, xRx);
                                else if (posDir && negDir)
                                    h_realsim_AccBias[ivar]->Fill(diff_index, xRx);
                            }
                        }
                    }
                }//end loop over variables
            }//end loop over the root file

            // Returns the next character in the input sequence, without extracting it: The character is left as the next character to be extracted from the stream
            file.peek();
            f->Close();
      
       }//end of while

       vector<double> mean, sigma, meanBias, sigmaBias, effIBD;
       double totalIBDs = 0.0, totalIBDsErr = 0.0, effIBDs = 0.0;
       for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
       {
            //IBD events = (Correlated - Accidental/100)
            h_realsim_Diff[ivar] = (TH1D*)h_realsim_Cor[ivar]->Clone("RealSimDiff");
            h_realsim_Diff[ivar]->Add(h_realsim_Acc[ivar], -1.0/100.0);

            h_realsim_DiffBias[ivar] = (TH1D*)h_realsim_CorBias[ivar]->Clone("RealSimDiffBias");
            h_realsim_DiffBias[ivar]->Add(h_realsim_AccBias[ivar], -1.0/100.0);
 
            totalIBDs = h_realsim_Diff[ivar]->IntegralAndError(0, h_realsim_Diff[ivar]->GetNbinsX() + 1, totalIBDsErr);
            effIBDs = pow(totalIBDs, 2) / pow(totalIBDsErr, 2); //Effective IBD counts. Done by Poisson Distribution N^2/(sqrt(N)^2) = N; Eff. counts = counts^2/counts_err^2
            cout << Form(" Total IBD events %0.0f ± %0.0f. Effective IBD counts %0.0f", totalIBDs, totalIBDsErr, effIBDs) << "\n";

            effIBD.push_back(effIBDs);

            h_realsim_Diff[ivar]->SetName( Form("RealSim_Diff_%s", vars[ivar].c_str()) );   
            h_realsim_Diff[ivar]->SetTitle( Form("%s Separation of Prompt and Delayed Events: IBD = %0.0f", vars[ivar].c_str(), totalIBDs) );
            h_realsim_Diff[ivar]->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]");
            h_realsim_Diff[ivar]->GetYaxis()->SetTitle("Total IBD Events");

            mean.push_back(h_realsim_Diff[ivar]->GetMean());
            sigma.push_back(h_realsim_Diff[ivar]->GetStdDev());

            if (ivar == 0 || ivar == 1){

                double rPlus, rMinus, px, pxErr, D = 145.7;
                double np, npp, npm, nm, nmm, n0;
                double x1, x2, x3, x4, x5;
                double npErr, nppErr, npmErr, nmErr, nmmErr, n0Err;
                double x1Err, x2Err, x3Err, x4Err, x5Err;
                double rPlusErr, rMinusErr;

                np = h_realsim_Diff[ivar]->GetBinContent(296);
                npp = h_realsim_DiffBias[ivar]->GetBinContent(296);
                nm = h_realsim_Diff[ivar]->GetBinContent(6);
                nmm = h_realsim_DiffBias[ivar]->GetBinContent(6);
                n0 = h_realsim_Diff[ivar]->GetBinContent(151);
                npm = h_realsim_DiffBias[ivar]->GetBinContent(151);

                npErr = h_realsim_Diff[ivar]->GetBinError(296);
                nppErr = h_realsim_DiffBias[ivar]->GetBinError(296);
                nmErr = h_realsim_Diff[ivar]->GetBinError(6);
                nmmErr = h_realsim_DiffBias[ivar]->GetBinError(6);
                n0Err = h_realsim_Diff[ivar]->GetBinError(151);
                npmErr = h_realsim_DiffBias[ivar]->GetBinError(151);
                
                x1 = np;
                x2 = npp + npm;
                x3 = nm;
                x4 = nmm + npm;
                
                rPlus = x1/x2;
                rMinus = x3/x4;

                x1Err = npErr;
                x2Err = sqrt(pow(nppErr, 2) + pow(npmErr, 2));
                x3Err = nmErr;
                x4Err = sqrt(pow(nmmErr, 2) + pow(npmErr, 2));

                rPlusErr = rPlus * sqrt(pow((x2Err/x2), 2) + pow((x1Err/x1), 2));
                rMinusErr = rMinus * sqrt(pow((x4Err/x4), 2) + pow((x3Err/x3), 2));

                px = D * ( (rPlus - rMinus) / (rPlus + rMinus + 1));

                pxErr = D * pow(1 / ((nm * (npm + npp) + (nmm + npm) * (np + npm + npp))), 2) * sqrt(pow((nmm + npm) * (npm + npp), 2) * (pow(npErr * (2 * nm + nmm + npm), 2)
                                                                                                                                        + pow(nmErr * (2 * np + npp + npm), 2))
                                                                                                    + pow((np * (npm + nmm) * (2 * nm + nmm + npm) * nppErr), 2)
                                                                                                    + pow((npmErr * (np * pow((nmm + npm), 2) + nm * (2 * nmm * np - 2 * np * npp - pow((npm + npp), 2)))), 2)
                                                                                                    + pow((nm * (npm + npp) * (2 * np + npm + npp) * nmmErr), 2));

                meanBias.push_back(px);
                sigmaBias.push_back(pxErr);
            }
            else{
                meanBias.push_back(h_realsim_Diff[ivar]->GetMean());
                sigmaBias.push_back(h_realsim_Diff[ivar]->GetStdDev());
            }
       }

       // Largest sigma value that best accounts for statistical error (usually in z direction). P in equation
       double lsigma = 0.0;
       if( (sigma[0] > sigma[1]) && (sigma[0] > sigma[2]) ) lsigma = sigma[0];
       if( (sigma[1] > sigma[0]) && (sigma[1] > sigma[2]) ) lsigma = sigma[1];
       if( (sigma[2] > sigma[0]) && (sigma[2] > sigma[1]) ) lsigma = sigma[2];

       sigmaBias[0] = sigmaBias[0] * sqrt(effIBD[0]);
       sigmaBias[1] = sigmaBias[1] * sqrt(effIBD[0]);

       // Angle Phi with respect to vertical between detector and reactor in xy plane
       double phiRad = atan2(mean[1], mean[0]); //phi = arctan(y/x)
       phiDeg_realsim = phiRad * 180.0/pi;
       double phiRadErr = sqrt( pow(lsigma/sqrt(effIBDs)*mean[1]/(mean[0]*mean[0]), 2) + pow(lsigma/sqrt(effIBDs)/(mean[0]), 2) ) / (1 + pow(mean[1]/mean[0], 2)); 
       phiDegErr_realsim = phiRadErr * 180.0/pi;

       // Angle Theta with respect to horizontal between detector and reactor in xz plane
       double thetaRad = atan2(mean[2], sqrt(mean[0]*mean[0] + mean[1]*mean[1])); //theta = arctan(z/sqrt(x^2 + y^2))
       thetaDeg_realsim = thetaRad * 180.0/pi;
       double thetaRadErr = sqrt( pow(lsigma/sqrt(effIBDs)/sqrt(mean[0]*mean[0] + mean[1]*mean[1]), 2) * (pow(mean[0]*mean[2]/(mean[0]*mean[0] + mean[1]*mean[1]), 2) + pow(mean[1]*mean[2]/(mean[0]*mean[0] + mean[1]*mean[1]), 2) + 1) ) / (1 + pow(mean[2]/sqrt(mean[0]*mean[0] + mean[1]*mean[1]), 2));
       thetaDegErr_realsim = thetaRadErr * 180.0/pi;

       cout << "  Mean x: "  << mean[0]  << ", Mean y: "  << mean[1]  << ", Mean z: "  << mean[2]  << "\n";
       cout << "  Sigma x: " << sigma[0] << ", Sigma y: " << sigma[1] << ", Sigma z: " << sigma[2] << "\n";
       cout << "  Largest sigma value: "    << lsigma              << "\n";
       cout << "  Phi in degrees: "         << phiDeg_realsim      << "\n";
       cout << "  Phi Error in degrees: "   << phiDegErr_realsim   << "\n";
       cout << "  Theta in degrees: "       << thetaDeg_realsim    << "\n";
       cout << "  Theta Error in degrees: " << thetaDegErr_realsim << "\n\n";

       // Unbiased Angle Calc
       // Angle Phi with respect to vertical between detector and reactor in xy plane
       double phiRadBias = atan2(meanBias[1], meanBias[0]); //phi = arctan(y/x)
       phiDeg_realsimBias = phiRadBias * 180.0/pi;
       double phiRadErrBias = sqrt( pow(sigmaBias[0]/sqrt(effIBD[0])*meanBias[1]/(meanBias[0]*meanBias[0]), 2) + pow(sigmaBias[1]/sqrt(effIBD[1])/(meanBias[0]), 2) ) / (1 + pow(meanBias[1]/meanBias[0], 2)); 
       phiDegErr_realsimBias = phiRadErrBias * 180.0/pi;

       // Angle Theta with respect to horizontal between detector and reactor in xz plane
       double thetaRadBias = atan2(meanBias[2], sqrt(meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1])); //theta = arctan(z/sqrt(x^2 + y^2))
       thetaDeg_realsimBias = thetaRadBias * 180.0/pi;
       double thetaRadErrBias = sqrt( (1 / (meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1])) * ( pow((meanBias[0]*meanBias[2]*sigmaBias[0]/sqrt(effIBD[0]) / (meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1])), 2) + pow((meanBias[1]*meanBias[2]*sigmaBias[1]/sqrt(effIBD[1]) / (meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1])), 2 ) + pow(sigmaBias[2]/sqrt(effIBD[2]), 2))) / (1 + pow(meanBias[2]/sqrt(meanBias[0]*meanBias[0] + meanBias[1]*meanBias[1]), 2));
       thetaDegErr_realsimBias = thetaRadErrBias * 180.0/pi;


       cout << "  Unbiased Mean x: "  << meanBias[0]  << ", Unbiased Mean y: "  << meanBias[1]  << ", Mean z: "  << meanBias[2]  << "\n";
       cout << "  Unbiased Sigma x: " << sigmaBias[0] << ", Unbiased Sigma y: " << sigmaBias[1] << ", Sigma z: " << sigmaBias[2] << "\n";
       cout << "  Largest sigma value: "    << lsigma              << "\n";
       cout << "  Unbiased Phi in degrees: "         << phiDeg_realsimBias      << "\n";
       cout << "  Unbiased Phi Error in degrees: "   << phiDegErr_realsimBias   << "\n";
       cout << "  Unbiased Theta in degrees: "       << thetaDeg_realsimBias    << "\n";
       cout << "  Unbiased Theta Error in degrees: " << thetaDegErr_realsimBias << "\n\n";
    }
    //----------------------------------------------------

    //----------------------------------------------------
    // Fill the perfect simulation: histograms & angles
    //----------------------------------------------------
    {
       string file_list = "/project/prospect/tmp/mmo58/prospect_bundle/MyWork/NeutrinoDirectionality/Simulations/SimFileGoodRunsList.txt";
       // Opening File
       std::ifstream file;
       file.open(file_list, std::ifstream::in);
       if( !(file.is_open() && file.good()) )
       {
           printf("Good runs file not found. Exiting\n");

           return -1;
       }
      
       int countlines = 0;

       while( file.good() && !file.eof() )
       {
           countlines = countlines + 1;
           cout << "Looking at Perfect Simulation file " << countlines << " / 20" << "\n";
 
           string line;
           getline(file, line);
           TString st = Form("/project/prospect/tmp/mmo58/prospect_bundle/MyWork/NeutrinoDirectionality/Simulations/H5Files/%s/Jun_Perfect_IBD_2020.root", line.data()); // Yale Meitner server
           //cout << "st = " << st << "\n";
           TFile *f = new TFile(st);

           TVectorD *rt = (TVectorD*)f->Get("runtime"); //total duration of the run
           TVectorD *promptv = (TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_prompt"); //prompt veto deadtime
           TVectorD *delayedv = (TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_delayed"); //delayed veto deadtime
           double xRx = rt->Max() / (rt->Max() - promptv->Max()) * rt->Max() / (rt->Max() - delayedv->Max()); //deadtime correction coefficient (veto deadtime correction)
      
           TTree *Th = (TTree*)f->Get("P2kIBDPlugin/Tibd");
           long nentries = Th->GetEntries();
           for( long i = 0; i < nentries; i++ )
           {
                Th->GetEntry(i);

                for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
                {
                     double Esmear = Th->GetLeaf("Esmear")->GetValue(0); //prompt energy in MeV. It is within PROSPECT's boundary conditions
                     double ncapt_dt = Th->GetLeaf("ncapt_dt")->GetValue(0); //neutron capture occurs a dt after the prompt signal. It is measured in ns
                     int pseg = Th->GetLeaf("maxseg")->GetValue(0); //prompt segment
                     int nseg = Th->GetLeaf("n_seg")->GetValue(0); //neutron segment
                     double delayed = Th->GetLeaf("n_xyz")->GetValue(ivar); //n_xyz is the neutron position in mm. Array of dim 3 (x, y, z)
                     double prompt = Th->GetLeaf("xyz")->GetValue(ivar); //xyz is the prompt position in mm. Array of dim 3 (x, y, z)
                     double diff_index = delayed - prompt;

                     if( pseg == nseg || diff_index > 0 || diff_index < 0 )
                     {
                         // Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal
                         if( Esmear > 0.8 && Esmear < 7.4 && ncapt_dt > pow(10, 3) && ncapt_dt < 120 * pow(10, 3) )
                             h_perfsim_Cor[ivar]->Fill(diff_index, xRx);
      
                         // Accidentals are IBDs with neutrons capturing over 1 ms after the prompt
                         if( Esmear > 0.8 && Esmear < 7.4 && ncapt_dt > pow(10, 6) )
                             h_perfsim_Acc[ivar]->Fill(diff_index, xRx*xRx);
                     }

                }//end loop over variables

           }//end loop over the root file

           // Returns the next character in the input sequence, without extracting it: The character is left as the next character to be extracted from the stream
           file.peek();
           f->Close();
      
       }//end of while

       vector<double> mean, sigma;
       double totalIBDs = 0.0, totalIBDsErr = 0.0, effIBDs = 0.0;
       for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
       {
            //IBD events = (Correlated - Accidental/100)
            h_perfsim_Diff[ivar] = (TH1D*)h_perfsim_Cor[ivar]->Clone("RealSimDiff");
            h_perfsim_Diff[ivar]->Add(h_perfsim_Acc[ivar], -1.0/100.0);
 
            totalIBDs = h_perfsim_Diff[ivar]->IntegralAndError(0, h_perfsim_Diff[ivar]->GetNbinsX() + 1, totalIBDsErr);
            effIBDs = pow(totalIBDs, 2) / pow(totalIBDsErr, 2); //Effective IBD counts. Done by Poisson Distribution N^2/(sqrt(N)^2) = N; Eff. counts = counts^2/counts_err^2
            cout << Form(" Total IBD events %0.0f ± %0.0f. Effective IBD counts %0.0f", totalIBDs, totalIBDsErr, effIBDs) << "\n";

            h_perfsim_Diff[ivar]->SetName( Form("PerfSim_Diff_%s", vars[ivar].c_str()) );   
            h_perfsim_Diff[ivar]->SetTitle( Form("%s Separation of Prompt and Delayed Events: IBD = %0.0f", vars[ivar].c_str(), totalIBDs) );
            h_perfsim_Diff[ivar]->GetXaxis()->SetTitle("Difference in position of prompt and delayed events [mm]");
            h_perfsim_Diff[ivar]->GetYaxis()->SetTitle("Total IBD Events");

            mean.push_back(h_perfsim_Diff[ivar]->GetMean());
            sigma.push_back(h_perfsim_Diff[ivar]->GetStdDev());
       }

       // Largest sigma value that best accounts for statistical error (usually in z direction). P in equation
       double lsigma = 0.0;
       if( (sigma[0] > sigma[1]) && (sigma[0] > sigma[2]) ) lsigma = sigma[0];
       if( (sigma[1] > sigma[0]) && (sigma[1] > sigma[2]) ) lsigma = sigma[1];
       if( (sigma[2] > sigma[0]) && (sigma[2] > sigma[1]) ) lsigma = sigma[2];

       // Angle Phi with respect to vertical between detector and reactor in xy plane
       double phiRad = atan2(mean[1], mean[0]); //phi = arctan(y/x)
       phiDeg_perfsim = phiRad * 180.0/pi;
       double phiRadErr = sqrt( pow(lsigma/sqrt(effIBDs)*mean[1]/(mean[0]*mean[0]), 2) + pow(lsigma/sqrt(effIBDs)/(mean[0]), 2) ) / (1 + pow(mean[1]/mean[0], 2)); 
       phiDegErr_perfsim = phiRadErr * 180.0/pi;

       // Angle Theta with respect to horizontal between detector and reactor in xz plane
       double thetaRad = atan2(mean[2], sqrt(mean[0]*mean[0] + mean[1]*mean[1])); //theta = arctan(z/sqrt(x^2 + y^2))
       thetaDeg_perfsim = thetaRad * 180.0/pi;
       double thetaRadErr = sqrt( pow(lsigma/sqrt(effIBDs)/sqrt(mean[0]*mean[0] + mean[1]*mean[1]), 2) * (pow(mean[0]*mean[2]/(mean[0]*mean[0] + mean[1]*mean[1]), 2) + pow(mean[1]*mean[2]/(mean[0]*mean[0] + mean[1]*mean[1]), 2) + 1) ) / (1 + pow(mean[2]/sqrt(mean[0]*mean[0] + mean[1]*mean[1]), 2));
       thetaDegErr_perfsim = thetaRadErr * 180.0/pi;

       cout << "  Mean x: "  << mean[0]  << ", Mean y: "  << mean[1]  << ", Mean z: "  << mean[2]  << "\n";
       cout << "  Sigma x: " << sigma[0] << ", Sigma y: " << sigma[1] << ", Sigma z: " << sigma[2] << "\n";
       cout << "  Largest sigma value: "    << lsigma              << "\n";
       cout << "  Phi in degrees: "         << phiDeg_perfsim      << "\n";
       cout << "  Phi Error in degrees: "   << phiDegErr_perfsim   << "\n";
       cout << "  Theta in degrees: "       << thetaDeg_perfsim    << "\n";
       cout << "  Theta Error in degrees: " << thetaDegErr_perfsim << "\n\n";
    }
    //----------------------------------------------------

    //----------------------------------------------------
    // Fill the true angles
    //----------------------------------------------------
    {
       // From figure 1 of the PRD paper, https://doi.org/10.1103/PhysRevD.103.032001 
       double x = 5.97, y = 5.09, z = -1.19; // Distance in [m]
       double xErr = 0.1, yErr = 0.1, zErr = 0.1; // Error in [m]. ±10 cm according to section "Experimental layout" in page 4

       // Angle Phi with respect to vertical between detector and reactor in xy plane
       double phiRad = atan2(y, x); //phi = arctan(y/x)
       phiDeg_true = phiRad * 180.0/pi;
       double phiRadErr = sqrt( pow(y/(x*x)*xErr, 2) + pow(1/x*yErr, 2) ) / (1 + pow(y/x, 2)); 
       phiDegErr_true = phiRadErr * 180.0/pi;

       // Angle Theta with respect to horizontal between detector and reactor in xz plane
       double thetaRad = atan2(z, sqrt(x*x + y*y)); //theta = arctan(z/sqrt(x^2 + y^2))
       thetaDeg_true = thetaRad * 180.0/pi;
       double thetaRadErr = sqrt( pow(1/sqrt(x*x + y*y), 2) * (pow(x*z/(x*x + y*y)*xErr, 2) + pow(y*z/(x*x + y*y)*yErr, 2) + pow(zErr, 2)) ) / (1 + pow(z/sqrt(x*x + y*y), 2));
       thetaDegErr_true = thetaRadErr * 180.0/pi;

       cout << "Looking at True Reactor" << "\n";
       cout << "  Phi in degrees: "         << phiDeg_true      << "\n";
       cout << "  Phi Error in degrees: "   << phiDegErr_true   << "\n";
       cout << "  Theta in degrees: "       << thetaDeg_true    << "\n";
       cout << "  Theta Error in degrees: " << thetaDegErr_true << "\n\n";
    }
    //----------------------------------------------------

    //----------------------------------------------------
    // Save stuff into the root file
    //----------------------------------------------------
    f_output->cd();
    for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
    {
         h_data_Diff[ivar]->Write();
         h_data_DiffBias[ivar]->Write();
         h_data_AccOffBias[ivar]->Write();
         h_data_AccOnBias[ivar]->Write();
         h_data_CorOffBias[ivar]->Write();
         h_data_CorOnBias[ivar]->Write();
         h_realsim_Diff[ivar]->Write();
         h_realsim_DiffBias[ivar]->Write();
         h_perfsim_Diff[ivar]->Write();
    }
    TVector2 *phi_data = new TVector2( phiDeg_data, phiDegErr_data );
    f_output->WriteTObject( phi_data, "phiDeg_data" ); 
    TVector2 *phi_dataBias = new TVector2( phiDeg_dataBias, phiDegErr_dataBias );
    f_output->WriteTObject( phi_dataBias, "phiDeg_dataBias" ); 
    TVector2 *theta_data = new TVector2( thetaDeg_data, thetaDegErr_data );
    f_output->WriteTObject( theta_data, "thetaDeg_data" );
    TVector2 *theta_dataBias = new TVector2( thetaDeg_dataBias, thetaDegErr_dataBias );
    f_output->WriteTObject( theta_dataBias, "thetaDeg_dataBias" ); 
    TVector2 *phi_realsim = new TVector2( phiDeg_realsim, phiDegErr_realsim );
    f_output->WriteTObject( phi_realsim, "phiDeg_realsim" ); 
    TVector2 *phi_realsimBias = new TVector2( phiDeg_realsimBias, phiDegErr_realsimBias );
    f_output->WriteTObject( phi_realsimBias, "phiDeg_realsimBias" );
    TVector2 *theta_realsim = new TVector2( thetaDeg_realsim, thetaDegErr_realsim );
    f_output->WriteTObject( theta_realsim, "thetaDeg_realsim" );
    TVector2 *theta_realsimBias = new TVector2( thetaDeg_realsimBias, thetaDegErr_realsimBias );
    f_output->WriteTObject( theta_realsimBias, "thetaDeg_realsimBias" ); 
    TVector2 *phi_perfsim = new TVector2( phiDeg_perfsim, phiDegErr_perfsim );
    f_output->WriteTObject( phi_perfsim, "phiDeg_perfsim" ); 
    TVector2 *theta_perfsim = new TVector2( thetaDeg_perfsim, thetaDegErr_perfsim );
    f_output->WriteTObject( theta_perfsim, "thetaDeg_perfsim" ); 
    TVector2 *phi_true = new TVector2( phiDeg_true, phiDegErr_true );
    f_output->WriteTObject( phi_true, "phiDeg_true" ); 
    TVector2 *theta_true = new TVector2( thetaDeg_true, thetaDegErr_true );
    f_output->WriteTObject( theta_true, "thetaDeg_true" ); 
    f_output->Close();
    //----------------------------------------------------


    printf("Program Complete!\n\n");

    return 0;

} // End Program
