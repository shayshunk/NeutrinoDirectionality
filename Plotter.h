
bool ApplyAxisStyle( TH1 *h, bool centerXTitle = true, bool centerYTitle = true, bool centerZTitle = true );
double Chi2DataMC( const TH1 *dataHist, const TH1 *mcHist, int & ndf, bool useOnlyShapeErrors = false );

void AddPlotLabel( const char* label,
                   const double x,
                   const double y,
                   const double size = 0.05,
                   const int color = 1,
                   const int font = 62,
                   const int align = 22,
                   const double angle = 0 );

void DecodePosition( const string& opts,
                     double size,
                     int &align,
                     double &xLabel,
                     double &yLabel );

void AddChi2Label( const TH1* dataHist,
                   const TH1* mcHist,
                   const string& opts,
                   double size = 0.04,
                   double yOffset = 0.0,
                   bool useOnlyShapeErrors = false,
                   bool useProb = false );


//-- marker settings
int data_marker = 20;
int ratio_marker = 20;
double data_marker_size = 1.3;
double ratio_marker_size = 1.3;
//double data_marker_size = 1.0;
//double ratio_marker_size = 1.0;

//-- line settings
//int data_line_width = 1;
int data_line_width = 3;
int data_line_style = 1;
int mc_line_width = 3;
int mc_line_style = 1;
int ratio_line_width = 3;

//-- color settings
int data_color = 1;
int mc_color   = 2;
int ratio_color = 1;

//-- title settings
int title_font = 62;
double title_size = 0.06;

//-- axis options
bool hist_min_zero    = true;
bool axis_draw_grid_x = false;
bool axis_draw_grid_y = false;
int axis_max_digits   = 3;
int axis_title_font_x = 62;
int axis_title_font_y = 62;
int axis_title_font_z = 62;
double axis_title_offset_x = 1.15;
double axis_title_offset_y = 1.2;
double axis_title_offset_z = .75;
double axis_title_size_x = 0.06;
double axis_title_size_y = 0.06;
double axis_title_size_z = 0.06;

//-- axis label options
double axis_label_font = 42;
double axis_label_size = 0.05;

//-- margins
double extra_top_margin = -0.02; //negative means go closer to edge

//-- correlation
bool draw_corr_max1 = false;
bool draw_corr_red_blue = true;

// Do you want to consider correlations to under/overflow bins in chi2 calculation?
bool chi2_use_overflow_err = false;


// Utility to set up basic root environment
void SetRootEnv()
{
     //gStyle->SetPalette(palette_style);

     // Canvas Styles
     gStyle->SetCanvasDefW(900);
     gStyle->SetCanvasDefH(750);
     gStyle->SetOptStat(0000);
     gStyle->SetOptFit(0000);
     gStyle->SetOptTitle(0);
     gStyle->SetCanvasColor(0);
     gStyle->SetPadBorderMode(0);
     gStyle->SetFrameBorderMode(0);
     gStyle->SetCanvasBorderMode(0);
     gStyle->SetPadTopMargin(0.09);
     gStyle->SetPadBottomMargin(0.15);
     gStyle->SetPadLeftMargin(0.15);
     gStyle->SetPadRightMargin(0.15);
     gStyle->SetFrameLineWidth(2);
     gStyle->SetHistLineWidth(2);

     // Axis Styles
     gStyle->SetHistMinimumZero( hist_min_zero );
     gStyle->SetTitleOffset( axis_title_offset_x, "X" );
     gStyle->SetTitleSize( axis_title_size_x, "X" );
     gStyle->SetTitleFont( axis_title_font_x, "X" );
     gStyle->SetTitleOffset( axis_title_offset_y, "Y" );
     gStyle->SetTitleSize( axis_title_size_y, "Y" );
     gStyle->SetTitleFont( axis_title_font_y, "Y" );
     gStyle->SetTitleOffset( axis_title_offset_z, "Z" );
     gStyle->SetTitleSize( axis_title_size_z, "Z" );
     gStyle->SetTitleFont( axis_title_font_z, "Z" );
     gStyle->SetLabelFont( axis_label_font, "XYZ" );
     gStyle->SetLabelSize( axis_label_size, "XYZ" );
     TGaxis::SetMaxDigits(axis_max_digits);
     gStyle->SetPadGridX( axis_draw_grid_x );
     gStyle->SetPadGridY( axis_draw_grid_y );

     // Marker Styles
     gStyle->SetMarkerStyle(data_marker);
     gStyle->SetMarkerSize(data_marker_size);
     gStyle->SetMarkerColor(data_color);
     gStyle->SetEndErrorSize(2);
     gStyle->SetErrorX(0.5);
}

bool ApplyAxisStyle( TH1 *h,
                     bool centerXTitle /*= true*/,
                     bool centerYTitle /*= true*/,
                     bool centerZTitle /*= true*/ )
{
     //!Set the X axis
     h->GetXaxis()->CenterTitle( centerXTitle );
     h->GetXaxis()->SetTitleOffset( axis_title_offset_x );
     h->GetXaxis()->SetTitleSize( axis_title_size_x );
     h->GetXaxis()->SetTitleFont( axis_title_font_x );
     h->GetXaxis()->SetLabelFont( axis_label_font );
     h->GetXaxis()->SetLabelSize( axis_label_size );

     //!Set the Y axis
     h->GetYaxis()->CenterTitle( centerYTitle );
     h->GetYaxis()->SetTitleOffset( axis_title_offset_y );
     h->GetYaxis()->SetTitleSize( axis_title_size_y );
     h->GetYaxis()->SetTitleFont( axis_title_font_y );
     h->GetYaxis()->SetLabelFont( axis_label_font );
     h->GetYaxis()->SetLabelSize( axis_label_size );

     //! Set the Z axis
     if( h->GetZaxis() != NULL )
     {
         h->GetZaxis()->CenterTitle( centerZTitle );
         h->GetZaxis()->SetTitleOffset( axis_title_offset_z );
         h->GetZaxis()->SetTitleSize( axis_title_size_z );
         h->GetZaxis()->SetTitleFont( axis_title_font_z );
         h->GetZaxis()->SetLabelFont( axis_label_font );
         h->GetZaxis()->SetLabelSize( axis_label_size );
     }

     return true;
}

// Easily add Latex labels on plots
void AddPlotLabel( const char* label,
                   const double x,
                   const double y,
                   const double size /*= 0.05*/,
                   const int color /*= 1*/,
                   const int font /*= 62*/,
                   const int align /*= 22*/,
                   const double angle /*= 0*/ )
{
     TLatex *latex = new TLatex( x, y, label );
     latex->SetNDC();
     latex->SetTextSize(size);
     latex->SetTextColor(color);
     latex->SetTextFont(font);
     latex->SetTextAlign(align);
     latex->SetTextAngle(angle);
     latex->Draw();
}

void AddHistoTitle( const char* title,
                    double titleSize,
                    int titleFont )
{
     AddPlotLabel(title, 0.5, 1 - gStyle->GetPadTopMargin() - extra_top_margin, titleSize, 1, titleFont, 21, 0.0);
}
/*
string PlotDir( string dir )
{
     string PLOTS_ROOT = getenv("PWD");
     string plotdir( PLOTS_ROOT + "/" + dir );

     if( 0 == system( Form("test -d %s", plotdir.c_str()) ) )
         system( Form("rm -r %s", plotdir.c_str()) ); 

     system( Form("mkdir -m 755 -p %s", plotdir.c_str()) );

     return plotdir;
}

string HistDir( string dir )
{
     string HISTS_ROOT = getenv("PWD");
     string histdir( HISTS_ROOT + "/" + dir );

     if( 0 == system( Form("test -d %s", histdir.c_str()) ) )
         system( Form("rm -r %s", histdir.c_str()) ); 

     system( Form("mkdir -m 755 -p %s", histdir.c_str()) );

     return histdir;
}
*/
string PlotDir( string dir )
{
     string PLOTS_ROOT = getenv("PWD");
     string plotdir( PLOTS_ROOT + "/" + dir );

     if( 0 != system( Form("test -d %s", plotdir.c_str()) ) )
         system( Form("mkdir -m 755 -p %s", plotdir.c_str()) );

     return plotdir;
}

string HistDir( string dir )
{
     string HISTS_ROOT = getenv("PWD");
     string histdir( HISTS_ROOT + "/" + dir );

     if( 0 != system( Form("test -d %s", histdir.c_str()) ) )
         system( Form("mkdir -m 755 -p %s", histdir.c_str()) );

     return histdir;
}

// Supply a comma-separated list of formats you want to print
void MultiPrint( TCanvas *c,
                 const string print_topdir,
                 const string& typeStr )
{
     string name = string(c->GetName());

     vector<string> types;
     size_t i = 0;
     size_t j = typeStr.find(',');
     while( j != string::npos )
     {
            types.push_back( typeStr.substr(i, j - i) );
            i = ++j;
            j = typeStr.find(',', i);
     }
     if( j == string::npos )
         types.push_back( typeStr.substr(i, typeStr.size()) );

     //we don't need an info statement here...
     const int oldVerbosity = gErrorIgnoreLevel;
     gErrorIgnoreLevel = kWarning;
     for( vector<string>::const_iterator itType = types.begin(); itType != types.end(); ++itType )
     {
          if( print_topdir.empty() )
              c->SaveAs( Form("%s.%s", name.c_str(), itType->c_str()), itType->c_str() );
          else
              c->SaveAs( Form("%s/%s.%s", print_topdir.c_str(), name.c_str(), itType->c_str()), itType->c_str() );
     }

     gErrorIgnoreLevel = oldVerbosity;
}

// Calculate the chi2 between two histograms
// The chi2 could also be calculated as tmp_h_data_IBD->Chi2Test(tmp_h_realsim_IBD, "WW P CHI2"). "WW" when both histograms are weighted. See ROOT documentation
// https://root.cern.ch/doc/master/classTH1.html#ab7d63c7c177ccbf879b5dc31f2311b27
// https://root-forum.cern.ch/t/compare-histograms-with-chi2test-which-option/37861
double Chi2DataMC( const TH1 * dataHist, const TH1 * mcHist, int & ndf, bool useOnlyShapeErrors )
{
    TH1D *tmpData = (TH1D*)dataHist->Clone("tmp_data_chi2");
    TH1D *tmpMC = (TH1D*)mcHist->Clone("tmp_mc_chi2");

    if( tmpData->GetSumw2N() == 0 )
        tmpData->Sumw2();
    if( tmpMC->GetSumw2N() == 0 )
        tmpMC->Sumw2();

    double chi2 = 0; 
    ndf = 0; 

    const int lowBin  = tmpMC->GetXaxis()->GetFirst();
    const int highBin = tmpMC->GetXaxis()->GetLast();

    for( int i = lowBin; i <= highBin; ++i )
    {
         if( tmpData->GetBinError(i)*tmpData->GetBinError(i) + tmpMC->GetBinError(i)*tmpMC->GetBinError(i) > 0 )
         {
             chi2 += (tmpData->GetBinContent(i) - tmpMC->GetBinContent(i))
                    *(tmpData->GetBinContent(i) - tmpMC->GetBinContent(i))
                    /(tmpData->GetBinError(i)*tmpData->GetBinError(i) + tmpMC->GetBinError(i)*tmpMC->GetBinError(i));

             ndf += 1;
         }
    }

    // If this is a shape comparison, subtract one degree of freedom
    if( useOnlyShapeErrors )
        ndf -= 1;

    return chi2;
}

// Get statistical error matrix
TMatrixD GetStatErrorMatrix( const TH1D *h )
{
   const int highBin = h->GetNbinsX() + 1;
   const int lowBin = 0;
   TMatrixD covmx(highBin + 1, highBin + 1);

   // Statistical error
   for( int iBin = lowBin; iBin <= highBin; ++iBin )
        covmx[iBin][iBin] = h->GetBinError(iBin);

   return covmx*covmx;
}

double Chi2DataMC_Test( const TH1D *dataHist, const TH1D *mcHist, int & ndf )
{
    TH1D *tmpData = (TH1D*)dataHist->Clone("tmp_data_chi2");
    TH1D *tmpMC = (TH1D*)mcHist->Clone("tmp_mc_chi2");

    if( tmpData->GetSumw2N() == 0 )
        tmpData->Sumw2();
    if( tmpMC->GetSumw2N() == 0 )
        tmpMC->Sumw2();

    const int lowBin  = tmpMC->GetXaxis()->GetFirst();
    const int highBin = tmpMC->GetXaxis()->GetLast();

    // Defining Error Matrix dimensions
    int Nbins = highBin - lowBin + 1;

    // Get the covariance matrix
    TMatrixD covMatrix(Nbins, Nbins);
    {
       const int NbinsTotal = tmpMC->GetNbinsX() + 2; //+2 for under/overflow
       TMatrixD covMatrixTmp(NbinsTotal, NbinsTotal);

       covMatrixTmp = GetStatErrorMatrix(tmpData);
       //covMatrixTmp = GetStatErrorMatrix(tmpMC);

       // Select only covariance elements in the histogram range
       // Unless using the contributions to covariance from overflow 
       if( chi2_use_overflow_err )
       {
           covMatrix = covMatrixTmp;
       }
       else
       {
           for( int i = 0; i != Nbins; ++i )
           {
                for( int j = 0; j != Nbins; ++j )
                {
                     covMatrix[i][j] = covMatrixTmp[i + lowBin][j + lowBin];
                }
           }
       }

    }// end of block to get covariance

    // Now, we invert the covariance error Matrix and store the result in "errorMatrix"
    // Note: TDecompSVD can handle singular matrices
    // ROOT can't seem to handle small entries with lots of zeros? Suggested scaling the histogram and then rescaling the inverted matrix
    //covMatrix *= 1e80;

    TDecompSVD error(covMatrix);
    TMatrixD errorMatrix(covMatrix);
    if( !error.Invert(errorMatrix) )
    {
         Warning("Chi2DataMC", "Cannot invert total covariance matrix. Using statistical errors only for Chi2 calculation.");

         return 0;
    }

    //errorMatrix *= 1e80;

    // Calculating chi2
    // Under/overflow bins not taken into account in the chi2 calculation
    ndf = 0;
    double chi2 = 0.;
    for( int i = lowBin; i <= highBin; ++i )
    {
         const int iErrBin = i - lowBin;
         const double x_data_i = tmpData->GetBinContent(i);
         const double x_mc_i   = tmpMC->GetBinContent(i);

         for( int j = lowBin; j <= highBin; ++j )
         {
              const int jErrBin = j - lowBin;
              const double x_data_j = tmpData->GetBinContent(j);
              const double x_mc_j   = tmpMC->GetBinContent(j);
              const double chi2_ij = (x_data_i - x_mc_i) * errorMatrix[iErrBin][jErrBin] * (x_data_j - x_mc_j);
              chi2 += chi2_ij;
         }

         ++ndf;
    }

    return chi2;
}

void DrawErrorMatrix( const TMatrixD & matrix, const TAxis * axis, const double maximum ) 
{
    // Create a 2D histogram with the matrix elements
    TH2D *h2D = new TH2D("h_matrix", Form("matrix;%s;%s", axis->GetTitle(), axis->GetTitle()),
                          axis->GetNbins(), axis->GetXbins()->GetArray(),
                          axis->GetNbins(), axis->GetXbins()->GetArray() );

    for( int i = 0; i < matrix.GetNrows(); i++ )
    {
         for( int j = 0; j < matrix.GetNcols(); j++ )
         {
              h2D->SetBinContent(h2D->GetBin(i, j), matrix[i][j]);
         }
    }

    h2D->GetXaxis()->SetRange(axis->GetFirst(), axis->GetLast());
    h2D->GetYaxis()->SetRange(axis->GetFirst(), axis->GetLast());
    h2D->SetTitleSize(axis->GetTitleSize(), "XY");
    h2D->GetXaxis()->SetNdivisions(axis->GetNdivisions());
    h2D->GetYaxis()->SetNdivisions(axis->GetNdivisions());

    // Set minimum and maximum by type
    double min = h2D->GetBinContent( h2D->GetMinimumBin() );
    double max = h2D->GetBinContent( h2D->GetMaximumBin() );

    // if maximum is supplied, use it, and assume symmetric
    if( maximum > 0.0 )
    {
        max = maximum;
        min = -maximum;
    }

    // For correlation matrix fix maximum to 1
    // This test usually means this is a correlation matrix
    if( draw_corr_max1 && 0 <= min && min < 1. && 0 < max && max < 1. )
    {
        h2D->SetMinimum( floor( min*10. )/10. );
        h2D->SetMaximum( 1.0 );
    }
    else
    {
        // Separate minimum and maximum if they are very close
        if( (max - min) < 0.1 )
        {
            int middle = floor( 0.5 + (min + max)/2.0 );
            h2D->SetMinimum( middle - 0.5 );
            h2D->SetMaximum( middle + 0.5 );
        }

        // If minimum is negative and maximum is positive
        // enter z axis in zero and make limits symetric
        else if( min*max < 0 )
        {
            const double absmax = std::max( fabs(min), fabs(max) );
            h2D->SetMinimum( -absmax );
            h2D->SetMaximum(  absmax );
        }

        // I don't understand why one would want to always separate by at least 1,
        // but I don't want to force the change.
        // I think the symmetric min/max is what we normally would want
        if( draw_corr_red_blue )
        {
            const double absmax = std::max( fabs(min), fabs(max) );
            h2D->SetMinimum( -absmax );
            h2D->SetMaximum(  absmax );
        }
    }

    // Draw the 2d histogram
    //if( draw_corr_red_blue )
        //SetCorrelationPalette();

    // Make the number of sigfigs on the z axis be the same, why this is not the default I do not know
    h2D->GetZaxis()->SetDecimals(true);

    h2D->DrawCopy("colz");

    delete h2D;
}

// Decodes a string to determine location, alignment of plot label
// If opts is a mixture of two strings like TR-TC, then use the average of those two positions and align of the first
void DecodePosition( const string & opts,
                     double size,
                     int & align,
                     double & xLabel,
                     double & yLabel )
{
     size_t dashLoc = opts.find("-");
     if( dashLoc != string::npos )
     {
         const string opts1 = opts.substr(0, dashLoc);
         int align1;
         double x1, y1;
         DecodePosition( opts1, size, align1, x1, y1 );

         const string opts2 = opts.substr(dashLoc+1);
         int align2;
         double x2, y2;
         DecodePosition( opts2, size, align2, x2, y2 );

         align = align1;
         xLabel = (x1 + x2 ) / 2.;
         yLabel = (y1 + y2) / 2.;
         return;
     }

     const double xLeft  = gStyle->GetPadLeftMargin() + 0.03;
     const double xCenter = .5;
     const double xRight = 1 - gStyle->GetPadRightMargin() - 0.025;

     const double yBottom = gStyle->GetPadBottomMargin() + size/2.;
     const double yCenter = .5;
     const double yTop    = 1 - gStyle->GetPadTopMargin() - size/2.;

     // Default is TC (top center)
     if( opts == "TC" || opts == "")
     {
         align = 23;
         xLabel = xCenter;
         yLabel = yTop;
     }
     else if( opts == "TR" )
     {
         align = 33;
         xLabel = xRight;
         yLabel = yTop;
     }
     else if( opts == "TL" )
     {
         align = 13;
         xLabel = xLeft;
         yLabel = yTop;
     }
     else if( opts == "BR" )
     {
         align = 31;
         xLabel = xRight;
         yLabel = yBottom;
     }
     else if( opts == "BL" )
     {
         align = 11;
         xLabel = xLeft;
         yLabel = yBottom;
     }
     else if( opts == "BC" )
     {
         align = 21;
         xLabel = xCenter;
         yLabel = yBottom;
     }
     else if( opts == "L" )
     {
         align = 12;
         xLabel = xLeft;
         yLabel = yCenter;
     }
     else if( opts == "C" )
     {
         align = 22;
         xLabel = xCenter;
         yLabel = yCenter;
     }
     else if( opts == "R" )
     {
         align = 32;
         xLabel = xRight;
         yLabel = yCenter;
     }
     else
     {
         Warning("DecodePosition", Form("Position option '%s' is not valid. No values have been set.", opts.c_str()));
     }
}

// Writes the chi2/ndf between two histograms on the plot
void AddChi2Label( const TH1 * dataHist,
                   const TH1 * mcHist,
                   const string & opts,
                   double size,
                   double yOffset,
                   bool useOnlyShapeErrors,
                   bool useProb )
{
     int align;
     double xLabel, yLabel;
     DecodePosition( opts, size, align, xLabel, yLabel );
     yLabel += yOffset;

     int ndf;
     double chi2 = Chi2DataMC( dataHist, mcHist, ndf, useOnlyShapeErrors );

     char *words;
     if( useProb )
     {
         double prob = dataHist->Chi2Test(mcHist, "WW P");
         if( prob < 0.01 )
             words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f, Prob = %3.2e", chi2, ndf, chi2/(Double_t)ndf, prob);
         else
             words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f, Prob = %3.2f", chi2, ndf, chi2/(Double_t)ndf, prob);
     }
     else
         words = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chi2, ndf, chi2/(Double_t)ndf);

     AddPlotLabel( words, xLabel, yLabel, size, 1, 62, align );
}
