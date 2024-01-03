#include <iostream>
#include <string>
#include <array>
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
	array<string, 3> axes{"x", "y", "z"};

	// Set up our output file
	//TFile *outputFile = new TFile("Directionality.root", "recreate");

	// Set up what we're measuring. Check enums in header for what the ints are
	array<float, 5> phi, phiError;
	array<float, 5> theta, thetaError;

	// Need histograms for counting each variable. Check enums in header for what the ints are
	// Don't need an array for the true reactor direction
	array<array<TH1D*, 4>, 4> histogram;

	cout << "Checking the enum: " << Data << "\n";

	return 0;
}