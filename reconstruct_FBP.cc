#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>

#include "libraries/CTProjectorClass.hh"
#include "libraries/FourierClass.hh"
#include "libraries/frameClass.hh"
using namespace std;

/**************************/
/**************************/
/* functions declarations */
/**************************/
/**************************/
void getParameters(const char *configurationFileName, string &sinogramFileNameFormat, int &projectionSize, string &reconstructionFileNameFormat, string &offsetsFilename, const char *configurationFileName2, int &startSlice, int &endSlice);
void printParameters(const char *configurationFileName, string sinogramFileNameFormat, int projectionSize, string reconstructionFileName, string offsetsFilename, int startSlice, int endSlice);
void reconstruct(string sinogramFileName, int projectionSize, string reconstructionFileName, int offset);
vector<int> loadAxisShiftsOffsets(string offsetsFilename);
vector<int> Split(string toparse);
void DumpVector(vector<int>, int);

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

	//int time0 = time(0);

	////////////////////////
	// parameters setting //
	////////////////////////
	int default_int = -777;
	string default_string = "EMPTY";
	string sinogramFileName = default_string;
	int projectionSize = default_int;
	string reconstructionFileNameFormat = default_string;
	string rotationAxisChi2FileName = default_string;
	//int rotationAxisRange = default_int;
	string sinogramFileNameFormat = default_string;
	string rotationAxisChi2FileNameFormat = default_string;
	string fixedSinogramFileNameFormat = default_string;
	string fixedSinogramFileName = default_string;
	string offsetsFilename = default_string;
	//int manualSinogramShift = default_int;
	int startSlice = default_int;
	int endSlice = default_int;


	////////////////////////////////////////////////
	// getting parameters from configuration file //
	////////////////////////////////////////////////
	getParameters(argv[1], sinogramFileNameFormat, projectionSize, reconstructionFileNameFormat, offsetsFilename, argv[2], startSlice, endSlice);
	printParameters(argv[1], sinogramFileNameFormat, projectionSize, reconstructionFileNameFormat, offsetsFilename, startSlice, endSlice);

	////////////////////////////////////////////////////////////////////////
	// Previous to reconstruction load the offsets from "offsetsFilename" //
	////////////////////////////////////////////////////////////////////////
	vector<int> offsets = loadAxisShiftsOffsets(offsetsFilename);
	cout << "\t- Offsets loaded. Number of slices in offsets file = " << offsets.size() << endl;
	DumpVector(offsets, 10);

	////////////////////
	// reconstruction //
	////////////////////
	// Check that I have enough offsets for all the slices
	if ( endSlice > (int)offsets.size()-1 ) {
		cout << "\t- [ERR ] the indexes don't quite match.  Check the startSlice, endSlice parameters.  Giving up." << endl;
		return 0;
	}
	// Now reconstruct
	for (int selectedSlice = startSlice; selectedSlice <= endSlice; selectedSlice++) {

		int size = 1000;
		char sinogramFileName[size];
		sprintf(sinogramFileName, sinogramFileNameFormat.c_str(), selectedSlice);
		char reconstructionFileName[size];
		sprintf(reconstructionFileName, reconstructionFileNameFormat.c_str(), selectedSlice);
		// Fetch the offset
		int offset = offsets[selectedSlice];
		// check the index is right
		reconstruct(sinogramFileName, projectionSize, reconstructionFileName, offset);
	}

	return 0;
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
/*************************/
/*************************/
/* functions definitions */
/*************************/
/*************************/
void getParameters(const char *configurationFileName, string &sinogramFileNameFormat, int &projectionSize, string &reconstructionFileNameFormat, string &offsetsFilename,  const char *configurationFileName2, int &startSlice, int &endSlice){

	ifstream configurationFile2; // this is the makeSinogram configuration file
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


	ifstream configurationFile; // this is the FBP reconstruction configuration file
	configurationFile.open(configurationFileName);
	if(!configurationFile){
		cout << " - ERROR!!! - getParameters(): cannot open file " << configurationFileName << endl;
		return ;
	}

	// file name format of the (centered axis) sinogram
	getline(configurationFile, sinogramFileNameFormat);

	// projection size
	getline(configurationFile, line);
	projectionSize = atoi(line.c_str());

	// file name format of the reconstruction
	getline(configurationFile, reconstructionFileNameFormat);

	// file containing the offsets after rotation axis finding
	getline(configurationFile, offsetsFilename);


	configurationFile.close();

	return ;
}

void printParameters(const char *configurationFileName, string sinogramFileNameFormat, int projectionSize, string reconstructionFileName, string offsetsFilename, int startSlice, int endSlice){

	cout << " - printParameters(): parameter values from " << configurationFileName << endl;
	cout << "\t-       file name format sinogram = " << sinogramFileNameFormat << endl;
	cout << "\t-                 projection size = " << projectionSize << endl;
	cout << "\t- file name format reconstruction = " << reconstructionFileName << endl;
	cout << "\t-                offsets filename = " << offsetsFilename << endl;
	cout << "\t-                     start slice = " << startSlice << endl;
	cout << "\t-                       end slice = " << endSlice << endl;

	return ;
}

void reconstruct(string sinogramFileName, int projectionSize, string reconstructionFileName, int offset){

	/////////////////
	// running FBP //
	/////////////////
	projectionSize = 512;
	cout << " - reconstruct(): reconstructing sinogram " << sinogramFileName << endl;
	frameClass *frame = new frameClass(projectionSize, projectionSize);
	frame -> set_to(0.);
	CTProjectorClass *projector = new CTProjectorClass(frame);
	delete frame;
	projector -> importSinogram(sinogramFileName.c_str(), projectionSize);

	// Shift first
	projector -> shiftSinogram_zeroPadding(offset);
	projectionSize = projector -> get_projectionSize();
	cout << " - reconstruct(): new projection size after finding rotation axis = " << projectionSize << endl;

	projector -> FBP();
	projector -> writeFrameToFile(reconstructionFileName.c_str());
	delete projector;
	cout << " - reconstruct(): written reconstruction " << reconstructionFileName << endl << endl;

	return ;
}
