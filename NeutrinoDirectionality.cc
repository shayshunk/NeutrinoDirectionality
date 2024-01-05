
#include "NeutrinoDirectionality.h"

using std::cout, std::string, std::ifstream, std::array, std::getline;

bool FillHistogram(array<array<array<std::shared_ptr<TH1D>, 3>, 5>, 4> &histogram, TreeValues &currentEntry)
{
	// Applying energy cut
	if (currentEntry.Esmear < 0.8 || currentEntry.Esmear > 7.4)
	{
		return false;
	}

	// Calculate neutron displacement
	double diffIndex = currentEntry.delayedPosition - currentEntry.promptPosition; 
	int signalSet;

	if (currentEntry.nCaptTime > pow(10, 3) && currentEntry.nCaptTime < 120 * pow(10, 3)) // Correlated Dataset
	{
		// Figure out whether the reactor is on and assign signal index
		if (currentEntry.reactorOn)
			{ signalSet = CorrelatedReactorOn; }
		else
			{ signalSet = CorrelatedReactorOff; }
		
		histogram[currentEntry.dataSet][signalSet][currentEntry.direction]->Fill(diffIndex);
		if (currentEntry.promptSegment == currentEntry.neutronSegment)
		{
			bool posDirection = false, negDirection = false;

			if (currentEntry.direction == x)
			{
				posDirection = checkNeighbor(detectorConfig, currentEntry.period, currentEntry.promptSegment, 'r');
				negDirection = checkNeighbor(detectorConfig, currentEntry.period, currentEntry.promptSegment, 'l');
			}
			else if (currentEntry.direction == y)
			{
				posDirection = checkNeighbor(detectorConfig, currentEntry.period, currentEntry.promptSegment, 'u');
				negDirection = checkNeighbor(detectorConfig, currentEntry.period, currentEntry.promptSegment, 'd');
			}

			if (posDirection && !negDirection)
				histogram[currentEntry.dataSet + 1][signalSet][currentEntry.direction]->Fill(segmentWidth);
			else if (!posDirection && negDirection)
				histogram[currentEntry.dataSet + 1][signalSet][currentEntry.direction]->Fill(-segmentWidth);
			else if (posDirection && negDirection)
				histogram[currentEntry.dataSet + 1[signalSet][currentEntry.direction]]->Fill(0.0);
		}
	}
	return true;
}

bool SetUpHistograms(array<array<array<std::shared_ptr<TH1D>, 3>, 5>, 4> &histogram, int dataSet, int period = 0)
{
	// Declaring some variables for use later
	int lineCounter = 0, totalLines = 0;
	bool reactorOn = true;
	double livetimeOff = 0, livetimeOn = 0;

	// Combining names into file list name
	string fileList = Form(dataPath, std::to_string(period), std::to_string(period));

	// Opening and checking file list
	ifstream file;
	file.open(fileList, ifstream::in);

	if ( !(file.is_open() && file.good()))
	{
		cout << "File list not found! Exiting.\n";
		return false;
	}

	while (file.good() && !file.eof())
	{
		lineCounter++;

		if (dataSet == Data || dataSet == DataUnbiased)
		{
			totalLines = 4000;

			if (lineCounter % 200 == 0)
				cout << "Looking at file: " << lineCounter << "/" << totalLines << '\n';	
		}
		else if (dataSet == Sim || dataSet == SimUnbiased)
		{
			totalLines = 20;

			cout << "Looking at file: " << lineCounter << "/" << totalLines << '\n';
		}

		// Reading file list
		string line;
		getline(file, line);

		// Combining names into root file name
		TString rootFilename = Form(fileName, std::to_string(period), line.data());	

		if (rootFilename.Contains(" 0"))
		{
			rootFilename.ReplaceAll(" 0", "");
			reactorOn = false;
		}
		else if (rootFilename.Contains(" 1"))
		{
			rootFilename.ReplaceAll(" 1", "");
			reactorOn = true;
		}

		// Open the root file
		auto rootFile = std::make_unique<TFile>(rootFilename);

		// Declare deadtime correction coefficient (veto deadtime correction)
		double xRx;

		// Going into empty scope to let the pointers die out for safety
		{
			TVectorD *runtime = (TVectorD*) rootFile->Get("runtime");
			TVectorD *promptVeto = (TVectorD*) rootFile->Get("accumulated/P2kIBDPlugin.tveto_prompt"); // prompt veto deadtime
			TVectorD *delayedVeto = (TVectorD*) rootFile->Get("accumulated/P2kIBDPlugin.tveto_delayed"); // delayed veto deadtime
			xRx = runtime->Max() / (runtime->Max() - promptVeto->Max()) * runtime->Max() / (runtime->Max() - delayedVeto->Max());

			if (reactorOn) 
				{ livetimeOn += runtime->Max() / xRx; }
			else 
				{ livetimeOff += runtime->Max() / xRx; }
		}

		// Grab tree and cast to unique pointer
		auto tree = std::unique_ptr<TTree> (static_cast<TTree*>(rootFile->Get("P2kIBDPlugin/Tibd")));

		long nEntries = tree->GetEntries();

		for (long i = 0; i < nEntries; i++)
		{
			tree->GetEntry(i);

			for (int direction = 0; direction < 3; direction++)
			{
				// Intializing struct of relevant values
				TreeValues currentEntry;

				// Grabbing relevant values from the tree entry
				currentEntry.Esmear = tree->GetLeaf("Esmear")->GetValue(0);
				currentEntry.nCaptTime = tree->GetLeaf("ncapt_dt")->GetValue(0);
				currentEntry.promptSegment = tree->GetLeaf("maxseg")->GetValue(0);
				currentEntry.neutronSegment = tree->GetLeaf("n_seg")->GetValue(0);
				currentEntry.promptPosition = tree->GetLeaf("xyz")->GetValue(direction);
				currentEntry.delayedPosition = tree->GetLeaf("n_xyz")->GetValue(direction);
				
				// Copying some loop values into current entry
				currentEntry.dataSet = dataSet;
				currentEntry.period = period;
				currentEntry.direction = direction;
				currentEntry.reactorOn = reactorOn;

				FillHistogram(histogram, currentEntry);
			}
		}

	}

	return true;
}

int main()
{
	// Ignore Warnings
	gErrorIgnoreLevel = kError;

	// Set up what we're measuring. Check enums in header for what the ints are
	array<float, 5> phi, phiError;
	array<float, 5> theta, thetaError;

	// Need histograms for counting each variable. Check enums in header for what the ints are
	// Don't need an array for the true reactor direction
	array<array<array<std::shared_ptr<TH1D>, 3>, 5>, 4> histogram;

	// Set up histograms for all 3 directions
	for (int i = 0; i < 4; i++) // Dataset
	{
		for (int j = 0; j < 5; j++) // Signal set
		{
			// No reactor off for simulations
			if ((i == Sim || i == SimUnbiased) && (j == CorrelatedReactorOff || j == AccidentalReactorOff)) continue;

			for (int c = 0; c < 3; c++)
			{
				string dataset = DatasetToString(i);
				string signalSet = SignalToString(j);
				string axis = AxisToString(c);
				string histogramName = dataset + "_" + signalSet + "_" + axis;
				histogram[c][i][j] = std::make_shared<TH1D>(histogramName.c_str(), dataset.c_str(), bins, -histogramMax, histogramMax);
			}
		}
	}

	// Set up our output file
	auto outputFile = std::make_unique<TFile>("Directionality.root", "recreate");
	outputFile->Close();

	return 0;
}