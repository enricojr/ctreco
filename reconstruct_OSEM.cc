#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <thread>

#include "libraries/CTProjectorClass.hh"
#include "libraries/FourierClass.hh"
#include "libraries/frameClass.hh"

using namespace std;

/**************************/
/**************************/
/* functions declarations */
/**************************/
/**************************/

typedef struct  {
	string sinogramFileName;
	int projectionSize;
	string reconstructionFileName;
	int OSEM_subsetSize;
	const char * chi2FileName;
	double chi2Cut_subsets;
	double chi2Cut_fullSet;
	string rotationAxisChi2FileName;
	int rotationAxisRange;
	string fixedSinogramFileName;
	bool searchAxis;
	int sinogramShift;
} threadArgs;

void getParameters(const char *configurationFileName, string &sinogramFileName, int &projectionSize, string &reconstructionFileName, int &OSEM_subsetSize, string &chi2FileNameFormat, double &chi2Cut_subsets, double &chi2Cut_fullSet, string &rotationAxisChi2FileName, int &rotationAxisRange, string &fixedSinogramFileName, bool &searchAxis, int &manualSinogramShift, const char *configurationFileName2, int &startSlice, int &endSlice, string &offsetFileName);
void printParameters(const char *configurationFileName, string sinogramFileName, int projectionSize, string reconstructionFileName, int OSEM_subsetSize, string chi2FileNameFormat, double chi2Cut_subsets, double chi2Cut_fullSet, string rotationAxisChi2FileName, int rotationAxisRange, string fixedSinogramFileName, bool searchAxis, int manualSinogramShift, int startSlice, int endSlice);
//void reconstruct(int threadId, string sinogramFileName, int projectionSize, string reconstructionFileName, int OSEM_subsetSize, const char *chi2FileName, double chi2Cut_subsets, double chi2Cut_fullSet, string rotationAxisChi2FileName, int rotationAxisRange, string fixedSinogramFileName, bool searchAxis, int manualSinogramShift);
void reconstruct(int threadId, threadArgs targs);
int findAxisShift (string rotationAxisChi2FileName, int rotationAxisRange, int manualSinogramShift, bool searchAxis, int projectionSize, string sinogramFileName);
vector<int> loadAxisShiftsOffsets(string offsetsFilename);
vector<int> Split(string toparse);
void DumpVector(vector<int>, int);
static const int GetConcurrentThreads();

/********/
/********/
/* main */
/********/
/********/
int main(int argc, char *argv[]){

	if(argc < 3){
		cout << " - " << argv[0] << " usage: " << endl;
		cout << "   " << argv[0] << " <configuration file name> <makeSinogram configuration file name>" << endl;
		return 0;
	}
	// Number of threads to be used
	static const int num_threads = GetConcurrentThreads();
	cout << "\t-   concurentThreads    = " << num_threads << endl;

	////////////////////////
	// parameters setting //
	////////////////////////
	int default_int = -777;
	int default_double = -888.666;
	string default_string = "EMPTY";
	string sinogramFileName = default_string;
	int projectionSize = default_int;
	string reconstructionFileNameFormat = default_string;
	int OSEM_subsetSize = default_int;
	string chi2FileNameFormat = default_string;
	double chi2Cut_subsets = default_double;
	double chi2Cut_fullSet = default_double;
	string rotationAxisChi2FileName = default_string;
	int rotationAxisRange = default_int;
	bool searchAxis = true;
	int manualSinogramShift = default_int;
	string sinogramFileNameFormat = default_string;
	string rotationAxisChi2FileNameFormat = default_string;
	string fixedSinogramFileNameFormat = default_string;
	string fixedSinogramFileName = default_string;
	int startSlice = default_int;
	int endSlice = default_int;
	int shift = default_int;
	string offsetFileName = default_string;

	////////////////////////////////////////////////
	// getting parameters from configuration file //
	////////////////////////////////////////////////
	getParameters(argv[1], sinogramFileNameFormat, projectionSize, reconstructionFileNameFormat, OSEM_subsetSize, chi2FileNameFormat, chi2Cut_subsets, chi2Cut_fullSet, rotationAxisChi2FileNameFormat, rotationAxisRange, fixedSinogramFileNameFormat, searchAxis, manualSinogramShift, argv[2], startSlice, endSlice, offsetFileName);
	printParameters(argv[1], sinogramFileNameFormat, projectionSize, reconstructionFileNameFormat, OSEM_subsetSize, chi2FileNameFormat, chi2Cut_subsets, chi2Cut_fullSet, rotationAxisChi2FileNameFormat, rotationAxisRange, fixedSinogramFileNameFormat, searchAxis, manualSinogramShift, startSlice, endSlice);
	int axisShift = 0;
	int count = 0;

	////////////////////////////////////////////////////////////////////////
	// Previous to reconstruction load the offsets from "offsetsFilename" //
	// If this is a OSEM reco job.  If this is run to find the axis then  //
	//  it is not needed												  //
	////////////////////////////////////////////////////////////////////////
	vector<int> offsets;
	if( !searchAxis ) {
		offsets = loadAxisShiftsOffsets(offsetFileName);
		cout << "\t- Offsets loaded. Number of slices in offsets file = " << offsets.size() << endl;
		if( offsets.size() == 0 ) {
			cout << " - ERROR!!! - couldn't read the offsets. Check on file " << offsetFileName << endl;
			return 1;
		}
		DumpVector(offsets, 10);
	}

	if(searchAxis) {

		ofstream offsetFile(offsetFileName.c_str());
		if(!offsetFile){
			cout << " - ERROR!!! - cannot open file " << offsetFileName << endl;
			return 1;
		}

		for(int selectedSlice = startSlice; selectedSlice <= endSlice; selectedSlice++){
			int size = 1000;
			char sinogramFileName[size];
			sprintf(sinogramFileName, sinogramFileNameFormat.c_str(), selectedSlice);
			char rotationAxisChi2FileName[size];
			sprintf(rotationAxisChi2FileName, rotationAxisChi2FileNameFormat.c_str(), selectedSlice);
			int currentShift = findAxisShift(rotationAxisChi2FileName, rotationAxisRange, manualSinogramShift, searchAxis, projectionSize, sinogramFileName);
			axisShift += currentShift;
			count++;
			offsetFile << selectedSlice << "\t" << currentShift << endl;
		}
		offsetFile.close();
		shift = axisShift / count;
		if( ( abs((double)axisShift/count) - floor(abs((double)axisShift/count)) ) > 0.5 ){ // to properly round the shift
			if(shift<0) shift--;
			else shift++;
		}
		cout << " - average shift = " << shift << endl;
		return 0;
	} else {
		//shift = manualSinogramShift;
		shift = 0; // The shift will be taken from the offset file on a slice per slice basis
	}
	cout << "------------------------------------------------------------------" << endl;
	cout << "----------------- \t axis shift = " << shift << "\t -----------------" << endl;
	cout << "------------------------------------------------------------------" << endl;

	////////////////////
	// reconstruction //
	////////////////////

	int threadCntr = 0;
	thread t[num_threads];
	threadArgs targs; // the arguments

	for(int selectedSlice = startSlice; selectedSlice <= endSlice; selectedSlice++) {

		// This is the exact shift for this particular slice
		shift = offsets[selectedSlice];

		int size = 1000;
		char sinogramFileName[size];
		sprintf(sinogramFileName, sinogramFileNameFormat.c_str(), selectedSlice);
		char reconstructionFileName[size];
		sprintf(reconstructionFileName, reconstructionFileNameFormat.c_str(), selectedSlice);
		char fixedSinogramFileName[size];
		sprintf(fixedSinogramFileName, fixedSinogramFileNameFormat.c_str(), selectedSlice);
		char chi2FileName[size];
		sprintf(chi2FileName, chi2FileNameFormat.c_str(), selectedSlice);

		// Launch the reconstruction for a slice in a thread
		targs.sinogramFileName = sinogramFileName;
		targs.projectionSize = projectionSize;
		targs.reconstructionFileName = reconstructionFileName;
		targs.OSEM_subsetSize = OSEM_subsetSize;
		targs.chi2FileName = chi2FileName;
		targs.chi2Cut_subsets = chi2Cut_subsets;
		targs.chi2Cut_fullSet = chi2Cut_fullSet;
		targs.rotationAxisChi2FileName = rotationAxisChi2FileName;
		targs.rotationAxisRange = rotationAxisRange;
		targs.fixedSinogramFileName = fixedSinogramFileName;
		targs.searchAxis = searchAxis;
		targs.sinogramShift = shift;

		t[threadCntr] = thread(reconstruct, threadCntr, targs);
		//t[threadCntr] = thread(reconstruct, threadCntr, sinogramFileName, projectionSize,
		//		reconstructionFileName, OSEM_subsetSize, chi2FileName, chi2Cut_subsets, chi2Cut_fullSet,
		//		rotationAxisChi2FileName, rotationAxisRange, fixedSinogramFileName, searchAxis, shift);
		threadCntr++;

		// In between launching the threads and joining them
		// extra work can be done

		// Pick up the threads when the maximum thread occ has been reached
		//  or if the selectedSlice has reached maximum
		if(threadCntr == num_threads || selectedSlice == endSlice) {
			// All launched.  Join them to the main thread.
			for (int i = 0; i < threadCntr; ++i) {
				t[i].join();
			}
			threadCntr = 0;
		}


	}

	cout << "--------------------------------------------------------------------------------------------------" << endl;
	cout << "----------------- \t new projection size after finding rotation axis = " << projectionSize + abs(shift) << "\t -----------------" <<  endl;
	cout << "--------------------------------------------------------------------------------------------------" << endl;

	return 0;
}

/*************************/
/*************************/
/* functions definitions */
/*************************/
/*************************/
void getParameters(const char *configurationFileName, string &sinogramFileNameFormat, int &projectionSize, string &reconstructionFileNameFormat, int &OSEM_subsetSize, string &chi2FileNameFormat, double &chi2Cut_subsets, double &chi2Cut_fullSet, string &rotationAxisChi2FileName, int &rotationAxisRange, string &fixedSinogramFileNameFormat, bool &searchAxis, int &manualSinogramShift, const char *configurationFileName2, int &startSlice, int &endSlice, string &offsetFileName){

	// reading makeSinogram configuration file
	ifstream configurationFile2;
	configurationFile2.open(configurationFileName2);
	if(!configurationFile2){
		cout << " - ERROR!!! - getParameters(): cannot open file " << configurationFileName2 << endl;
		return ;
	}

	string line;
	getline(configurationFile2, line);
	getline(configurationFile2, line);
	getline(configurationFile2, line);
	getline(configurationFile2, line);
	getline(configurationFile2, line);
	getline(configurationFile2, line);

	getline(configurationFile2, line);
	startSlice = atoi(line.c_str());

	getline(configurationFile2, line);
	endSlice = atoi(line.c_str());

	configurationFile2.close();

	// reading reconstruct_OSEM configuration file
	ifstream configurationFile;
	configurationFile.open(configurationFileName);
	if(!configurationFile){
		cout << " - ERROR!!! - getParameters(): cannot open file " << configurationFileName << endl;
		return ;
	}

	// format of the input sinogram file name
	getline(configurationFile, sinogramFileNameFormat);

	// projection size
	getline(configurationFile, line);
	projectionSize = atoi(line.c_str());

	// format of the reconstruction file name
	getline(configurationFile, reconstructionFileNameFormat);

	// OSEM subset size
	getline(configurationFile, line);
	OSEM_subsetSize = atoi(line.c_str());

	// format of the chi2 file name
	getline(configurationFile, chi2FileNameFormat);

	// chi2 cut on subset iterations
	getline(configurationFile, line);
	chi2Cut_subsets = atof(line.c_str());

	// chi2 cut on full set iterations
	getline(configurationFile, line);
	chi2Cut_fullSet = atof(line.c_str());

	// format of the file names for the chi2 of the rotation axis finding routine
	getline(configurationFile, rotationAxisChi2FileName);

	// range around the center where to look for the rotation axis
	getline(configurationFile, line);
	rotationAxisRange = atoi(line.c_str());

	// format of the sinogram after shift of the rotation axis
	getline(configurationFile, fixedSinogramFileNameFormat);

	// boolean: search axis
	getline(configurationFile, line);
	searchAxis = atoi(line.c_str());

	// manual axis shift
	getline(configurationFile, line);
	manualSinogramShift = atoi(line.c_str());

	// name of the offset file
	getline(configurationFile, line);
	offsetFileName = line;

	configurationFile.close();

	return ;
}

void printParameters(const char *configurationFileName, string sinogramFileName, int projectionSize, string reconstructionFileName, int OSEM_subsetSize, string chi2FileNameFormat, double chi2Cut_subsets, double chi2Cut_fullSet, string rotationAxisChi2FileName, int rotationAxisRange, string fixedSinogramFileName, bool searchAxis, int manualSinogramShift, int /*startSlice*/, int /*endSlice*/){

	cout << " - printParameters(): parameter values from " << configurationFileName << endl;
	cout << "\t-                            file name sinogram = " << sinogramFileName << endl;
	cout << "\t-                               projection size = " << projectionSize << endl;
	cout << "\t-                      file name reconstruction = " << reconstructionFileName << endl;
	cout << "\t-                              OSEM subset size = " << OSEM_subsetSize << endl;
	cout << "\t-                           OSEM chi2 file name = " << chi2FileNameFormat << endl;
	cout << "\t-            OSEM chi2 cut on subset iterations = " << chi2Cut_subsets << endl;
	cout << "\t-          OSEM chi2 cut on full set iterations = " << chi2Cut_fullSet << endl;
	cout << "\t-        file name chi2 of rotation axis finder = " << rotationAxisChi2FileName << endl;
	cout << "\t-                           rotation axis range = " << rotationAxisRange << endl;
	cout << "\t- file name sinogram after fixing rotation axis = " << fixedSinogramFileName << endl;
	if(searchAxis) cout << "\t-             automatic search of rotation axis = YES" << endl;
	else{
		cout << "\t-             automatic search of rotation axis = NO" << endl;
		cout << "\t-                         manual sinogram shift = " << manualSinogramShift << endl;
	}

	return ;
}

void reconstruct(int threadId, threadArgs targs) {

	//void reconstruct(int threadId, string sinogramFileName, int projectionSize, string reconstructionFileName, int OSEM_subsetSize, const char *chi2FileName, double chi2Cut_subsets, double chi2Cut_fullSet, string /*rotationAxisChi2FileName*/, int /*rotationAxisRange*/, string fixedSinogramFileName, bool /*searchAxis*/, int sinogramShift){

	cout << "*** Launched in thread " << threadId << " ***" << endl;
	cout << "producing " << targs.reconstructionFileName << endl;
	cout << "producing " << targs.fixedSinogramFileName << endl;

	frameClass *frame = new frameClass(targs.projectionSize, targs.projectionSize);
	frame -> set_to(0.);
	CTProjectorClass *projector = new CTProjectorClass(frame);
	delete frame;
	projector -> importSinogram(targs.sinogramFileName.c_str(), targs.projectionSize);
	projector -> shiftSinogram_zeroPadding(targs.sinogramShift);
	targs.projectionSize = projector -> get_projectionSize();
	cout << " ---------- SHIFT = " << targs.sinogramShift << endl;
	cout << " - reconstruct(): new projection size after finding rotation axis = " << targs.projectionSize << endl;
	projector -> writeSinogram(targs.fixedSinogramFileName.c_str());

	//  return ;

	//////////////////
	// running OSEM //
	//////////////////
	cout << " - reconstruct(): running OSEM" << endl;
	frame = new frameClass(targs.projectionSize, targs.projectionSize);
	frame -> set_to(0.);
	projector = new CTProjectorClass(frame);
	delete frame;
	projector -> importSinogram(targs.fixedSinogramFileName.c_str(), targs.projectionSize);
	projector -> OSEM(targs.OSEM_subsetSize, targs.chi2FileName, targs.chi2Cut_subsets, targs.chi2Cut_fullSet);
	projector -> writeFrameToFile(targs.reconstructionFileName.c_str());
	delete projector;

	return;
}

int findAxisShift (string rotationAxisChi2FileName, int rotationAxisRange, int /*manualSinogramShift*/, bool /*searchAxis*/, int projectionSize, string sinogramFileName){

	//////////////////////////
	// fixing rotation axis //
	//////////////////////////
	cout << " - reconstruct(): finding rotation axis" << endl;
	frameClass *frame = new frameClass(projectionSize, projectionSize);
	frame -> set_to(0.);
	CTProjectorClass *projector = new CTProjectorClass(frame);
	delete frame;
	projector -> importSinogram(sinogramFileName.c_str(), projectionSize);
	cout << " - findAxisShift(): writing " << rotationAxisChi2FileName << endl;
	int shift = projector -> fixRotationAxis_zeroPadding(rotationAxisChi2FileName.c_str(), rotationAxisRange);
	//  shift = projector -> fixRotationAxis(rotationAxisChi2FileName.c_str(), rotationAxisRange);
	delete projector;

	return shift;
}

vector<int> loadAxisShiftsOffsets(string offsetsFilename) {

	vector<int> offsets;
	ifstream ifs (offsetsFilename.c_str(), std::ifstream::in);

	char temp[1024];
	string cmpString;
	//int lineCntr = 0;
	vector<int> oneLine;

	while ( ifs.good() ) {

		ifs.getline(temp, 1024);

		// avoid parsing in the EOF
		if ( ifs.eof() ) break;

		cmpString = string(temp);
		oneLine = Split(cmpString);

		// The second column is the offset
		offsets.push_back( oneLine[1] );

	}

	return offsets;
}

vector<int> Split(string toparse) {

	vector<int> split;
	istringstream iss (toparse);
	do {
		string sub;
		iss >> sub;
		split.push_back( atof(sub.c_str()) );
	} while (iss);

	return split;
}

// Dump a few values
void DumpVector(vector<int> v, int n_head_tail) {

	int N = v.size();
	bool dots = false;
	cout << "\t- Dump vector < ";
	for (int i = 0 ; i < N ; i++ ) {
		if( i < n_head_tail || i > N - n_head_tail - 1) {
			cout << v[i] << " ";
		}
		if ( !dots ) {
			if( i > n_head_tail && i < N - n_head_tail - 1 ) {
				cout << " ..... ";
				dots = true; // done
			}
		}

	}
	cout << "> " << endl;

}
static const int GetConcurrentThreads() {

	// Determine the number of concurrent threads
	cout << "\t-   Detecting number of cores using std::thread::hardware_concurrency ... ";
	unsigned concurentThreadsSupported = thread::hardware_concurrency();
	// If it is not detectable thread::hardware_concurrency() will return 0
	if ( concurentThreadsSupported == 0 ) {
		cout << "couldn't detect it.  Will try an OS specific call." << endl;
	} else {
		cout << "Successful." << endl;
		return concurentThreadsSupported * 2;
	}

	cout << "\t-   Trying sysconf( _SC_NPROCESSORS_ONLN ) [linux specific] ... ";
	int concurentThreadsSupported_linux = sysconf( _SC_NPROCESSORS_ONLN ); // this value can be negative if call fails

	if ( concurentThreadsSupported_linux <= 0 ) {
		concurentThreadsSupported_linux = 2;
		cout << "couldn't detect it. Manually setting " << concurentThreadsSupported_linux * 2 << " threads as value by default." << endl;
	} else {
		cout << "Successful." << endl;
		return concurentThreadsSupported_linux * 2;
	}

	return concurentThreadsSupported_linux * 2;
}
