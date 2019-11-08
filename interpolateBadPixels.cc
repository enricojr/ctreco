#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include "libraries/frameClass.hh"
using namespace std;

/************************/
/************************/
/* function definitions */
/************************/
/************************/
void readConfigurationFile(char *configurationFileName, string &fileNameFormat, string &OBfileName, int &CTposMin, int &CTposMax, int &CTposStep, int &nx, int &ny);
bool OBFileOK(string OBfileName, int nx, int ny);

/********/
/********/
/* main */
/********/
/********/
int main(int argc, char *argv[]){

	if(argc < 2){
		cout << " - " << argv[0] << " usage: " << endl;
		cout << "   " << argv[0] << " <configuration file name>" << endl;
		return 0;
	}

	char *configurationFileName = argv[1];
	string fileNameFormat = "";
	string OBfileName = "";
	int CTposMin = -1;
	int CTposMax = -1;
	int CTposStep = -1;
	int nx = -1;
	int ny = -1;

	readConfigurationFile(configurationFileName, fileNameFormat, OBfileName, CTposMin, CTposMax, CTposStep, nx, ny);

	cout << " - fileNameFormat = " << fileNameFormat << endl;
	cout << " - OBfileName = " << OBfileName << endl;
	cout << " - CTposMin = " << CTposMin << endl;
	cout << " - CTposMax = " << CTposMax << endl;
	cout << " - CTposStep = " << CTposStep << endl;
	cout << " - nx = " << nx << endl;
	cout << " - ny = " << ny << endl;

	// Having the OBfile here is very important.
	// If it is not present the data will be overwritten with zeroes
	// Checking that it is available first and that it is not a zero matrix.
	if ( ! OBFileOK(OBfileName, nx, ny) ) {
		cout << " - Something went wrong with the OB file ... giving up." << endl;
		return 1;
	}


	frameClass *OBframe = new frameClass(nx, ny);
	OBframe -> readFromFile_ASCIIMatrix(OBfileName.c_str());
	OBframe -> binarize(1);

	int size = 1000;
	for(int CTpos=CTposMin; CTpos<=CTposMax; CTpos+=CTposStep){
		cout << " - processing CTpos = " << CTpos << endl;

		char fileName[size];
		sprintf(fileName, fileNameFormat.c_str(), CTpos);
		frameClass *frame = new frameClass(nx, ny);
		frame -> readFromFile_ASCIIMatrix(fileName);

		frameClass *newFrame = new frameClass(nx, ny);
		newFrame -> multiplyFrames(frame, OBframe);
		delete frame;

		newFrame -> interpolateBadPixels(0.);

		newFrame -> writeToFile_ASCIIMatrix(fileName);
		delete newFrame;
	}

	OBframe -> readFromFile_ASCIIMatrix(OBfileName.c_str());
	OBframe -> interpolateBadPixels(0.);
	OBframe -> writeToFile_ASCIIMatrix(OBfileName.c_str());

	delete OBframe;

	return 0;
}

/*************************/
/*************************/
/* function declarations */
/*************************/
/*************************/

void readConfigurationFile(char *configurationFileName, string &fileNameFormat, string &OBfileName, int &CTposMin, int &CTposMax, int &CTposStep, int &nx, int &ny){

	ifstream configurationFile;
	configurationFile.open(configurationFileName);
	if(configurationFile == NULL){
		cout << " - ERROR!!! - readConfigurationFile(): cannot open configuration file " << configurationFileName << endl;
		return ;
	}

	string line;

	// format of the input frames files
	getline(configurationFile, line);
	fileNameFormat = line;

	// format of the input open beam files
	getline(configurationFile, line);
	OBfileName = line;

	// minimum angle position
	getline(configurationFile, line);
	CTposMin = atoi(line.c_str());

	// maximum angle position
	getline(configurationFile, line);
	CTposMax = atoi(line.c_str());

	// angle step
	getline(configurationFile, line);
	CTposStep = atoi(line.c_str());

	// number of pixels nx
	getline(configurationFile, line);
	nx = atoi(line.c_str());

	// number of pixels ny
	getline(configurationFile, line);
	ny = atoi(line.c_str());

	configurationFile.close();

	return ;
}

bool OBFileOK(string fileName, int nx, int ny){

	ifstream file;
	file.open(fileName.c_str());
	if(file == NULL)  {
		cout << " - ERROR!!! - cannot open file " << fileName << endl;
		return false;
	}

	//  Now check that this is not a file full with zeroes
	long int addEntries = 0;
	long int entry = 0;
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			file >> entry;
			addEntries += entry;
		}
	}
	file.close();

	if ( addEntries == 0 ) {
		cout << " - ERROR!!! - " << fileName << " seems to be a matrix full of zeroes." << endl;
		return false;
	}

	return true;
}
