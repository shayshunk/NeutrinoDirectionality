/*
For information regarding IBD Selection Cuts refer to PROSPECT2xAnalysis/Analysis/PhysPulse/IBDCutset.hh and IBD Selection Rules.pdf

IBD selection is defined at https://docdb.wlab.yale.edu/prospect/docs/0030/003086/004/PRD_technote.pdf (DocDB 3086-v4)
IBD Selection = (CorrelatedOn - AccidentalOn*x/AccWindow) - AtmScale*(LiveTimeOn/LiveTimeOff)*(CorrelatedOff - AccidentalOff*x/AccWindow)
  Where:
  . x = runtime/(runtime - prompt veto deadtime) * runtime/(runtime - delayed veto deadtime) //deadtime correction factor for each file
  . AccWindow = 100.0 number of accidental windows (matched to OnTime) //accidental scaling factor
  . AtmScale = 1.00025443769309 //atmospheric scaling

Realistic Simulation:
. In this simulation there are no backgrounds (Off) and no Accidentals
. This code does not simulate backgrounds
. Use three different calibration files to simulate time drifting across data period
. Accidentals are not that significant but I will leave the code for them

Realistic Simulation per Period:
. The perfect simulation is used as a realistic simulation per period
. The map of live segments was extracted for each period by looking at the segment number where the prompt is detected
. So we know the dead segments, and they are removed in the perfect simulation
. In this way, we do not need to generate a new simulation

Branch names:
. https://docdb.wlab.yale.edu/prospect/docs/0012/001215/001/SimpleIBD.pdf (DocDB 1215-v1)
. https://github.com/PROSPECT-collaboration/PROSPECT2x_Analysis/blob/master/Analysis/PhysPulse/IBDTree.hh
*/

#include "../utils/GeneralHeader.h"
#include "../utils/Plotter.h"
#include "../utils/LiveSegments.h"

double atm_scaling = 1.00025443769309; //atmospheric scaling factor

void AddBins( vector<double>& binsLowEdge, const double binWidth, const double lowBin, const double highBin )
{
    double x = lowBin;
    if( binsLowEdge.size() && lowBin == binsLowEdge.back() )
        x += binWidth;
    while( x <= highBin )
    {
        binsLowEdge.push_back( x );
        x += binWidth;
    }
}

vector<double> GetBins( string var )
{
    vector<double> binsLowEdge;

    if( var == "X" || var == "Y" )
    {
        double tmpBins[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }

    if( var == "Z" )
        AddBins( binsLowEdge, 1.0, -150.5, 150.5 );

    return binsLowEdge;
}


int UnbiasedSegCorrPRLHists( string period )
{
    // Ignore Warnings
    gErrorIgnoreLevel = kError;

    TH1::SetDefaultSumw2();

    vector<string> vars;
    vars.push_back("X");
    vars.push_back("Y");
    vars.push_back("Z");

    string histFileName = HistDir("UnbiasedSegCorr_PRL_Hists");
    //string histFileName = HistDir("UnbiasedSegCorr_PRL_OldCal_Hists");
    histFileName += Form("/Hists_UnbiasedSegCorr_PRL_%s.root", period.c_str());
    cout << "The histogram file name is " << histFileName << endl;

    // Save stuff into a root file
    TFile *f_output = new TFile(histFileName.c_str(), "recreate");

    // Make a histogram for each variable
    vector<TH1D*> h_data_coracc_on, h_data_coracc_off, h_data_acc_on, h_data_acc_off, h_data_IBD;
    vector<TH1D*> h_realsim_coracc, h_realsim_acc, h_realsim_IBD;
    vector<TH1D*> h_perfsim_coracc, h_perfsim_acc, h_perfsim_IBD;

    for( vector<string>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
         cout << " var " << var->c_str() << endl;

         vector<double> binsLowEdge = GetBins(*var);

         h_data_coracc_on.push_back ( new TH1D(Form("Data_CorAcc_RxOn_%s",  var->c_str()), "Data", binsLowEdge.size()-1, &binsLowEdge[0]) );
         h_data_coracc_off.push_back( new TH1D(Form("Data_CorAcc_RxOff_%s", var->c_str()), "Data", binsLowEdge.size()-1, &binsLowEdge[0]) );
         h_data_acc_on.push_back    ( new TH1D(Form("Data_Acc_RxOn_%s",     var->c_str()), "Data", binsLowEdge.size()-1, &binsLowEdge[0]) );
         h_data_acc_off.push_back   ( new TH1D(Form("Data_Acc_RxOff_%s",    var->c_str()), "Data", binsLowEdge.size()-1, &binsLowEdge[0]) );
         h_data_IBD.push_back       ( new TH1D(Form("Data_IBD_%s",          var->c_str()), "Data", binsLowEdge.size()-1, &binsLowEdge[0]) );

         h_realsim_coracc.push_back( new TH1D(Form("RealSim_CorAcc_%s", var->c_str()), "RealSim", binsLowEdge.size()-1, &binsLowEdge[0]) );
         h_realsim_acc.push_back   ( new TH1D(Form("RealSim_Acc_%s",    var->c_str()), "RealSim", binsLowEdge.size()-1, &binsLowEdge[0]) );
         h_realsim_IBD.push_back   ( new TH1D(Form("RealSim_IBD_%s",    var->c_str()), "RealSim", binsLowEdge.size()-1, &binsLowEdge[0]) );

         h_perfsim_coracc.push_back( new TH1D(Form("PerfSim_CorAcc_%s", var->c_str()), "PerfSim", binsLowEdge.size()-1, &binsLowEdge[0]) );
         h_perfsim_acc.push_back   ( new TH1D(Form("PerfSim_Acc_%s",    var->c_str()), "PerfSim", binsLowEdge.size()-1, &binsLowEdge[0]) );
         h_perfsim_IBD.push_back   ( new TH1D(Form("PerfSim_IBD_%s",    var->c_str()), "PerfSim", binsLowEdge.size()-1, &binsLowEdge[0]) );
    }

    double runtime_on = 0.0; //run time for reactor on periods
    double runtime_off = 0.0; //run time for reactor off periods
    double livetime_on = 0.0; //live time for reactor on periods
    double livetime_off = 0.0; //live time for reactor off periods
    double sec_to_days = 1.0/86400.0;

    period.erase(period.begin(), period.end() - 1);
    cout << " Period " << period << endl;    

    //----------------------------------------------------
    // Fill the data: histograms
    //----------------------------------------------------
    {
       //string file_list = Form("/project/prospect/tmp/DS_SEER_Data_OldCal/SEER_DS_period_%s/Period_%s_files.txt", period.c_str(), period.c_str());
       string file_list = Form("/project/prospect/tmp/DS_SEER_Data_NewCal/SEER_DS_period_%s/Period_%s_files.txt", period.c_str(), period.c_str());

       // Counting the number of lines in the text file
       string lin;
       int numLines = 0;
       ifstream in;
       in.open(file_list, std::ifstream::in);
       while( !in.eof() )
       {
              getline(in, lin);
              numLines = numLines + 1;
       }
       in.close();

       // Opening File
       std::ifstream file;
       file.open(file_list, std::ifstream::in);
       if( !(file.is_open() && file.good()) )
       {
           printf("Good runs file not found. Exiting\n");

           return -1;
       }

       int countlines = 0;
       bool ReactorOn = true;

       while( file.good() && !file.eof() )
       {
           countlines = countlines + 1;
           if( countlines % 100 == 0 ) cout << "Looking at Data file " << countlines << " / " << numLines << endl;

           string line;
           getline(file, line);
           //TString st = Form("/project/prospect/tmp/DS_SEER_Data_OldCal/SEER_DS_period_%s/%s/AD1_IBD_2022_DS_SEER_fid4.root", period.c_str(), line.data()); // Yale Meitner server
           TString st = Form("/project/prospect/tmp/DS_SEER_Data_NewCal/SEER_DS_period_%s/%s/AD1_IBD_2022_DS_SEER.root", period.c_str(), line.data()); // Yale Meitner server
           //cout << "st = " << st << endl;

           // Formating issue within file directory to correct it
           // 0 refers to reactor off while 1 is reactor on
           if( st.Contains(" 0") )
           {
               st.ReplaceAll(" 0", "");
               ReactorOn = false;
           }

           if( st.Contains(" 1") )
           {
               st.ReplaceAll(" 1", "");
               ReactorOn = true;
           }

           TFile *f = new TFile(st);
           //cout << "st = " << st << endl;

           TVectorD *rt = (TVectorD*)f->Get("runtime"); //total duration of the run
           TVectorD *promptv = (TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_prompt"); //prompt veto deadtime
           TVectorD *delayedv = (TVectorD*)f->Get("accumulated/P2kIBDPlugin.tveto_delayed"); //delayed veto deadtime
           double xRx = rt->Max() / (rt->Max() - promptv->Max()) * rt->Max() / (rt->Max() - delayedv->Max()); //deadtime correction coefficient (veto deadtime correction)
           if( ReactorOn )
           {
               runtime_on += rt->Max();
               livetime_on += rt->Max() / xRx;
           }
           else
           {
               runtime_off += rt->Max();
               livetime_off += rt->Max() / xRx;
           }

           TTree *Th = (TTree*)f->Get("P2kIBDPlugin/Tibd");
           long nentries = Th->GetEntries();
           for( long i = 0; i < nentries; i++ )
           {
                Th->GetEntry(i);

                for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
                {
                     double Esmear = Th->GetLeaf("Esmear")->GetValue(0); //prompt energy in MeV. It is within PROSPECT's boundary conditions
                     double ncapt_dt = Th->GetLeaf("ncapt_dt")->GetValue(0); //neutron capture occurs a dt after the prompt signal. It is measured in ns
                     int pseg = Th->GetLeaf("maxseg")->GetValue(0); //maximum-energy segment number
                     int nseg = Th->GetLeaf("n_seg")->GetValue(0); //neutron capture segment number
                     double delayed = Th->GetLeaf("n_xyz")->GetValue(ivar); //n_xyz is the neutron position in mm. Array of dim 3 (x, y, z)
                     double prompt = Th->GetLeaf("xyz")->GetValue(ivar); //xyz is the prompt position in mm. Array of dim 3 (x, y, z)
                     double diff_index = delayed - prompt; //separation of prompt and delayed events

                     int m = 0;
                     string direction_pos, direction_neg; // right/left/up/down direction
                     if( vars[ivar] == "X" )
                     {
                         m = 1;
                         direction_pos = "right";
                         direction_neg = "left";
                     }
                     if( vars[ivar] == "Y" )
                     {
                         m = 14;
                         direction_pos = "up";
                         direction_neg = "down";
                     }

                     if( 0.8 < Esmear && Esmear < 7.4 )
                     {
                         // Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal
                         if( 1.0 < ncapt_dt/1000.0 && ncapt_dt/1000.0 < 120.0 ) // IBD candidates
                         {
                             if( ReactorOn )
                             {
                                 if( vars[ivar] == "X" || vars[ivar] == "Y" )
                                 {
                                     if( pseg == nseg - m ) //N+
                                     {
                                         diff_index = 1.0;
                                     }
                                     if( pseg == nseg + m ) //N-
                                     {
                                         diff_index = 2.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0++
                                     {
                                         diff_index = 3.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) ) //N0--
                                     {
                                         diff_index = 4.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0+-
                                     {
                                         diff_index = 5.0;
                                     }
                                 }

                                 h_data_coracc_on[ivar]->Fill(diff_index);
                             }
                             else
                             {
                                 if( vars[ivar] == "X" || vars[ivar] == "Y" )
                                 {
                                     if( pseg == nseg - m ) //N+
                                     {
                                         diff_index = 1.0;
                                     }
                                     if( pseg == nseg + m ) //N-
                                     {
                                         diff_index = 2.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0++
                                     {
                                         diff_index = 3.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) ) //N0--
                                     {
                                         diff_index = 4.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0+-
                                     {
                                         diff_index = 5.0;
                                     }
                                 }

                                 h_data_coracc_off[ivar]->Fill(diff_index);
                             }
                         }
                         else // Accidentals IBDs
                         {
                             if( ReactorOn )
                             {
                                 if( vars[ivar] == "X" || vars[ivar] == "Y" )
                                 {
                                     if( pseg == nseg - m ) //N+
                                     {
                                         diff_index = 1.0;
                                     }
                                     if( pseg == nseg + m ) //N-
                                     {
                                         diff_index = 2.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0++
                                     {
                                         diff_index = 3.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) ) //N0--
                                     {
                                         diff_index = 4.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0+-
                                     {
                                         diff_index = 5.0;
                                     }
                                 }

                                 h_data_acc_on[ivar]->Fill(diff_index, xRx/100.0);
                             }
                             else
                             {
                                 if( vars[ivar] == "X" || vars[ivar] == "Y" )
                                 {
                                     if( pseg == nseg - m ) //N+
                                     {
                                         diff_index = 1.0;
                                     }
                                     if( pseg == nseg + m ) //N-
                                     {
                                         diff_index = 2.0;
                                     }                            
                                     if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0++
                                     {
                                         diff_index = 3.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) ) //N0--
                                     {
                                         diff_index = 4.0;
                                     }
                                     if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0+-
                                     {
                                         diff_index = 5.0;
                                     }
                                 }

                                 h_data_acc_off[ivar]->Fill(diff_index, xRx/100.0);
                             }
                         }

                     }//end of if Esmear

                }//end loop over variables

           }//end loop over the root file

           // Returns the next character in the input sequence, without extracting it: The character is left as the next character to be extracted from the stream
           file.peek();
           f->Close();

       }//end of while

       cout << Form(" Rx-Off run time (days) %0.2f", runtime_off*sec_to_days) << endl;
       cout << Form(" Rx-On  run time (days) %0.2f", runtime_on*sec_to_days)  << "\n\n";
    }
    //----------------------------------------------------

    //----------------------------------------------------
    // Fill the realistic simulation: histograms
    //----------------------------------------------------
    {
       for( int j = 1; j <= 100; j++ )
       {
            if( j % 10 == 0 ) cout << "Looking at Realistic Simulation file " << j << " / 100" << endl;

            TString st = Form("/project/prospect/tmp/DS_SEER_MC/period_%s/p2x_p%s/AD1_IBD_RUN%i/AD1_IBD_2022_DS_SEER.root", period.c_str(), period.c_str(), j); // Yale Meitner server
            //cout << "st = " << st << endl;
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
                      int pseg = Th->GetLeaf("maxseg")->GetValue(0); //maximum-energy segment number
                      int nseg = Th->GetLeaf("n_seg")->GetValue(0); //neutron capture segment number
                      double delayed = Th->GetLeaf("n_xyz")->GetValue(ivar); //n_xyz is the neutron position in mm. Array of dim 3 (x, y, z)
                      double prompt = Th->GetLeaf("xyz")->GetValue(ivar); //xyz is the prompt position in mm. Array of dim 3 (x, y, z)
                      double diff_index = delayed - prompt; //separation of prompt and delayed events
            
                      int m = 0;
                      string direction_pos, direction_neg; // right/left/up/down direction
                      if( vars[ivar] == "X" )
                      {
                          m = 1;
                          direction_pos = "right";
                          direction_neg = "left";
                      }
                      if( vars[ivar] == "Y" )
                      {
                          m = 14;
                          direction_pos = "up";
                          direction_neg = "down";
                      }
            
                      if( 0.8 < Esmear && Esmear < 7.4 )
                      {
                          // Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal
                          if( 1.0 < ncapt_dt/1000.0 && ncapt_dt/1000.0 < 120.0 ) // IBD candidates
                          {
                              if( vars[ivar] == "X" || vars[ivar] == "Y" )
                              {
                                  if( pseg == nseg - m ) //N+
                                  {
                                      diff_index = 1.0;
                                  }
                                  if( pseg == nseg + m ) //N-
                                  {
                                      diff_index = 2.0;
                                  }
                                  if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0++
                                  {
                                      diff_index = 3.0;
                                  }
                                  if( pseg == nseg && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) ) //N0--
                                  {
                                      diff_index = 4.0;
                                  }
                                  if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0+-
                                  {
                                      diff_index = 5.0;
                                  }
                              }
            
                              h_realsim_coracc[ivar]->Fill(diff_index);
                          }
                          else // Accidentals IBDs
                          {
                              if( vars[ivar] == "X" || vars[ivar] == "Y" )
                              {
                                  if( pseg == nseg - m ) //N+
                                  {
                                      diff_index = 1.0;
                                  }
                                  if( pseg == nseg + m ) //N-
                                  {
                                      diff_index = 2.0;
                                  }
                                  if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0++
                                  {
                                      diff_index = 3.0;
                                  }
                                  if( pseg == nseg && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) && !LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) ) //N0--
                                  {
                                      diff_index = 4.0;
                                  }
                                  if( pseg == nseg && LiveSegments(pseg, direction_pos, Form("Data_period%s", period.c_str())) && LiveSegments(pseg, direction_neg, Form("Data_period%s", period.c_str())) ) //N0+-
                                  {
                                      diff_index = 5.0;
                                  }
                              }
            
                              h_realsim_acc[ivar]->Fill(diff_index, xRx/100.0);
                          }
            
                      }//end of if Esmear
            
                 }//end loop over variables

            }//end loop over the root file

            f->Close();

       }//end of loop over files

       cout << "\n";
    }
    //----------------------------------------------------

    //----------------------------------------------------
    // Fill the perfect simulation: histograms
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
           cout << "Looking at Perfect Simulation file " << countlines << " / 20" << endl;

           string line;
           getline(file, line);
           TString st = Form("/project/prospect/tmp/mmo58/prospect_bundle/MyWork/NeutrinoDirectionality/Simulations/H5Files/%s/Jun_Perfect_IBD_2020.root", line.data()); // Yale Meitner server
           //cout << "st = " << st << endl;
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
                     int pseg = Th->GetLeaf("maxseg")->GetValue(0); //maximum-energy segment number
                     int nseg = Th->GetLeaf("n_seg")->GetValue(0); //neutron capture segment number
                     double delayed = Th->GetLeaf("n_xyz")->GetValue(ivar); //n_xyz is the neutron position in mm. Array of dim 3 (x, y, z)
                     double prompt = Th->GetLeaf("xyz")->GetValue(ivar); //xyz is the prompt position in mm. Array of dim 3 (x, y, z)
                     double diff_index = delayed - prompt; //separation of prompt and delayed events

                     // Exclude segment 32. The population (intensity) of this segment is 50% lower than the neighbor segments. An issue with the calibration might cause this
                     if( pseg == 32 || nseg == 32 )
                         continue;

                     int m = 0;
                     string direction_pos, direction_neg; // right/left/up/down direction
                     if( vars[ivar] == "X" )
                     {
                         m = 1;
                         direction_pos = "right";
                         direction_neg = "left";
                     }
                     if( vars[ivar] == "Y" )
                     {
                         m = 14;
                         direction_pos = "up";
                         direction_neg = "down";
                     }

                     if( 0.8 < Esmear && Esmear < 7.2 )
                     {
                         // Correlated IBDs are IBDs where the neutron capture occurs between 1 and 120 us after the prompt signal
                         if( 1.0 < ncapt_dt/1000.0 && ncapt_dt/1000.0 < 120.0 ) // IBD candidates
                         {
                             if( vars[ivar] == "X" || vars[ivar] == "Y" )
                             {
                                 if( pseg == nseg - m ) //N+
                                 {
                                     diff_index = 1.0;
                                 }
                                 if( pseg == nseg + m ) //N-
                                 {
                                     diff_index = 2.0;
                                 }
                                 if( pseg == nseg && LiveSegments(pseg, direction_pos, "PerfSim") && !LiveSegments(pseg, direction_neg, "PerfSim") ) //N0++
                                 {
                                     diff_index = 3.0;
                                 }
                                 if( pseg == nseg && LiveSegments(pseg, direction_neg, "PerfSim") && !LiveSegments(pseg, direction_pos, "PerfSim") ) //N0--
                                 {
                                     diff_index = 4.0;
                                 }
                                 if( pseg == nseg && LiveSegments(pseg, direction_pos, "PerfSim") && LiveSegments(pseg, direction_neg, "PerfSim") ) //N0+-
                                 {
                                     diff_index = 5.0;
                                 }
                             }

                             h_perfsim_coracc[ivar]->Fill(diff_index);
                         }
                         else // Accidentals IBDs
                         {
                             if( vars[ivar] == "X" || vars[ivar] == "Y" )
                             {
                                 if( pseg == nseg - m ) //N+
                                 {
                                     diff_index = 1.0;
                                 }
                                 if( pseg == nseg + m ) //N-
                                 {
                                     diff_index = 2.0;
                                 }
                                 if( pseg == nseg && LiveSegments(pseg, direction_pos, "PerfSim") && !LiveSegments(pseg, direction_neg, "PerfSim") ) //N0++
                                 {
                                     diff_index = 3.0;
                                 }
                                 if( pseg == nseg && LiveSegments(pseg, direction_neg, "PerfSim") && !LiveSegments(pseg, direction_pos, "PerfSim") ) //N0--
                                 {
                                     diff_index = 4.0;
                                 }
                                 if( pseg == nseg && LiveSegments(pseg, direction_pos, "PerfSim") && LiveSegments(pseg, direction_neg, "PerfSim") ) //N0+-
                                 {
                                     diff_index = 5.0;
                                 }
                             }

                             h_perfsim_acc[ivar]->Fill(diff_index, xRx/100.0);
                         }

                     }//end of if Esmear

                }//end loop over variables

           }//end loop over the root file

           // Returns the next character in the input sequence, without extracting it: The character is left as the next character to be extracted from the stream
           file.peek();
           f->Close();

       }//end of while

       cout << "\n";
    }
    //----------------------------------------------------

    cout << " I am going to calculate IBD signal events for data and simulation" << endl;
    double total_IBD, total_IBD_err;

    for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
    {
         // Data IBD events
         // IBD events = (Correlated - Accidental)_{reactor on} - AtmScale*livetimeOn/livetimeOff*(Correlated - Accidental)_{reactor off}
         h_data_IBD[ivar] = (TH1D*)h_data_coracc_on[ivar]->Clone(Form("Data_IBD_%s", vars[ivar].c_str()));
         h_data_IBD[ivar]->Add(h_data_acc_on[ivar], -1.0);
         h_data_IBD[ivar]->Add(h_data_coracc_off[ivar], -atm_scaling*livetime_on/livetime_off);
         h_data_IBD[ivar]->Add(h_data_acc_off[ivar], atm_scaling*livetime_on/livetime_off);
         total_IBD = h_data_IBD[ivar]->IntegralAndError(0, h_data_IBD[ivar]->GetNbinsX() + 1, total_IBD_err);
         //cout << Form("  Data IBD signal events in %s: %0.3f ± %0.3f", vars[ivar].c_str(), total_IBD, total_IBD_err) << endl;
         cout << Form("  Data IBD signal events in %s: %0.0f ± %0.0f", vars[ivar].c_str(), total_IBD, total_IBD_err) << endl;

         // Simulation IBD events
         // IBD events = (Correlated - Accidental)
         h_realsim_IBD[ivar] = (TH1D*)h_realsim_coracc[ivar]->Clone(Form("RealSim_IBD_%s", vars[ivar].c_str()));
         h_realsim_IBD[ivar]->Add(h_realsim_acc[ivar], -1.0);
         total_IBD = h_realsim_IBD[ivar]->IntegralAndError(0, h_realsim_IBD[ivar]->GetNbinsX() + 1, total_IBD_err);
         cout << Form("  Realistic Simulation IBD signal events in %s: %0.3f ± %0.3f", vars[ivar].c_str(), total_IBD, total_IBD_err) << endl;
         //cout << Form("  Realistic Simulation IBD signal events in %s: %0.0f ± %0.0f", vars[ivar].c_str(), total_IBD, total_IBD_err) << endl;

         // Simulation IBD events
         // IBD events = (Correlated - Accidental)
         h_perfsim_IBD[ivar] = (TH1D*)h_perfsim_coracc[ivar]->Clone(Form("PerfSim_IBD_%s", vars[ivar].c_str()));
         h_perfsim_IBD[ivar]->Add(h_perfsim_acc[ivar], -1.0);
         total_IBD = h_perfsim_IBD[ivar]->IntegralAndError(0, h_perfsim_IBD[ivar]->GetNbinsX() + 1, total_IBD_err);
         cout << Form("  Perfect Simulation IBD signal events in %s: %0.3f ± %0.3f", vars[ivar].c_str(), total_IBD, total_IBD_err) << endl;
         //cout << Form("  Perfect Simulation IBD signal events in %s: %0.0f ± %0.0f", vars[ivar].c_str(), total_IBD, total_IBD_err) << endl;
    }

    cout << " I calculated IBD signal events" << endl;

    // Save stuff into the root file
    f_output->cd();
    for( unsigned int ivar = 0; ivar != vars.size(); ++ivar )
    {
         h_data_IBD[ivar]->Write();
         h_realsim_IBD[ivar]->Write();
         h_perfsim_IBD[ivar]->Write();
    }
    f_output->Close();

    printf("\nProgram Complete!\n\n");

    return 0;

} // End Program


int main( int argc, char** argv )
{
    if( argc == 2 )
    {
        cout << Form("Running neutrino directionality hists (Unbiased segment correction) for data %s", argv[1]) << endl;
        return UnbiasedSegCorrPRLHists( argv[1] );
    }

    else
    {
        std::cerr << "Incorrect input parameters. Use:"     << std::endl;
        std::cerr << "  ./UnbiasedSegCorrPRLHists [period]" << std::endl;
        std::cerr << "  ./UnbiasedSegCorrPRLHists period1"  << std::endl;
        return -1;
    }
}
