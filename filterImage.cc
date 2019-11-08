#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include "libraries/frameClass.hh"
using namespace std;

void readConfigurationFile(char *configurationFileName, string &fileNameFormat, int &nx, int &ny, int &sigma, string &outputFileNameFormat, const char *configurationFileName2, int &startSlice, int &endSlice);

int main(int argc, char *argv[]){

  if(argc < 3){
    cout << " - " << argv[0] << " usage: " << endl;
    cout << "   " << argv[0] << " <configuration file name> <makeSinogram configuration file name>" << endl;
    return 0;
  }

  string fileNameFormat = "";
  int nx = -1;
  int ny = -1;
  int sigma = -1;
  string outputFileNameFormat = "";
  int startSlice = -1;
  int endSlice = -1;

  readConfigurationFile(argv[1], fileNameFormat, nx, ny, sigma, outputFileNameFormat, argv[2], startSlice, endSlice);

  cout << " - fileNameFormat = " << fileNameFormat << endl;
  cout << " - nx = " << nx << endl;
  cout << " - ny = " << ny << endl;
  cout << " - sigma = " << sigma << endl;
  cout << " - outputFileNameFormat = " << outputFileNameFormat << endl;
  cout << " - startSlice = " << startSlice << endl;
  cout << " - endSlice = " << endSlice << endl;

  for(int i=startSlice; i<=endSlice; i++){
    char fileNameIn[1000];
    sprintf(fileNameIn, fileNameFormat.c_str(), i);
    cout << " - reading image " << fileNameIn << endl;
    frameClass *frame = new frameClass(nx, ny);
    frame -> readFromFile_ASCIIMatrix(fileNameIn);
    frameClass *filtered = frame -> applyGaussFilter(sigma);
    delete frame;
    char fileNameOut[1000];
    sprintf(fileNameOut, outputFileNameFormat.c_str(), i);
    cout << " - writing image " << fileNameOut << endl;
    filtered -> writeToFile_ASCIIMatrix(fileNameOut);
    delete filtered;
  }

  return 0;
}

void readConfigurationFile(char *configurationFileName, string &fileNameFormat, int &nx, int &ny, int &sigma, string &outputFileNameFormat, const char *configurationFileName2, int &startSlice, int &endSlice){

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

  // number of pixels nx
  getline(configurationFile, line);
  nx = atoi(line.c_str());

  // number of pixels ny
  getline(configurationFile, line);
  ny = atoi(line.c_str());

  // Gaussian filter sigma
  getline(configurationFile, line);
  sigma = atoi(line.c_str());

  // format of the output frames files
  getline(configurationFile, line);
  outputFileNameFormat = line;

  configurationFile.close();

  ifstream configurationFile2; // this is the makeSinogram configuration file
  configurationFile2.open(configurationFileName2);
  if(configurationFile2 == NULL){
    cout << " - ERROR!!! - getParameters(): cannot open file " << configurationFileName2 << endl;
    return ;
  }

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

  return ;
}
