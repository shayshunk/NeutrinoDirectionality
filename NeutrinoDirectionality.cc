
#include "NeutrinoDirectionality.h"

using std::cout;
using std::string;
using std::array;
using std::ifstream;
using std::getline;

int main()
{
	// Ignore Warnings
	gErrorIgnoreLevel = kError;

	// Defining the 3 directions for calculations
	array<char, 3> axes{'x', 'y', 'z'};

	// Set up what we're measuring. Check enums in header for what the ints are
	array<float, 5> phi, phiError;
	array<float, 5> theta, thetaError;

	// Need histograms for counting each variable. Check enums in header for what the ints are
	// Don't need an array for the true reactor direction
	array<array<TH1D*, 4>, 4> histogram;

	// Set up histograms for all 3 directions
	for (const char &axis : axes)
	{
		for (int i = 0; i < 4; i++) // Dataset
		{
			for (int j = 0; j < 5; j++) // Signal set
			{
				string dataset = DatasetToString(i);
				string signalSet = SignalToString(j);
				string direction = string(1, axis);
				string histogramName = dataset + "_" + signalSet + "_" + axis;
				histogram[i][j] = new TH1D(histogramName.c_str(), dataset.c_str(), bins, -histogramMax, histogramMax);
			}
		}
	}

	// Set up our output file
	TFile *outputFile = new TFile("Directionality.root", "recreate");
	outputFile->Close();

	return 0;
}