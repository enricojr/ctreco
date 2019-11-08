#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include "libraries/frameClass.hh"
using namespace std;

void readConfigurationFile(char *configurationFileName, string &fileNameFormat, string &OBfileName, int &CTposMin, int &CTposMax, int &CTposStep, int &nx, int &ny, string &outputFileNameFormat);

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
  string outputFileNameFormat = "";

  readConfigurationFile(configurationFileName, fileNameFormat, OBfileName, CTposMin, CTposMax, CTposStep, nx, ny, outputFileNameFormat);

  cout << " - fileNameFormat = " << fileNameFormat << endl;
  cout << " - OBfileName = " << OBfileName << endl;
  cout << " - CTposMin = " << CTposMin << endl;
  cout << " - CTposMax = " << CTposMax << endl;
  cout << " - CTposStep = " << CTposStep << endl;
  cout << " - nx = " << nx << endl;
  cout << " - ny = " << ny << endl;
  cout << " - outputFileNameFormat = " << outputFileNameFormat << endl;

  frameClass *OBframe = new frameClass(nx, ny);
  OBframe -> readFromFile_ASCIIMatrix(OBfileName.c_str());

  int size = 1000;
  for(int CTpos=CTposMin; CTpos<=CTposMax; CTpos+=CTposStep){
    cout << " - processing CTpos = " << CTpos << endl;

    char fileName[size];
    sprintf(fileName, fileNameFormat.c_str(), CTpos);
    frameClass *frame = new frameClass(nx, ny);
    frame -> readFromFile_ASCIIMatrix(fileName);

    frameClass *normalizedFrame = new frameClass(nx, ny);
    normalizedFrame -> divideFrames(frame, OBframe);
    delete frame;

    ////////////////////////////////
    ////////////////////////////////
    // THIS IS VERY IMPORTANT !!! //
    ////////////////////////////////
    ////////////////////////////////
    normalizedFrame -> minusLog();
    ////////////////////////////////
    ////////////////////////////////
    ////////////////////////////////
    ////////////////////////////////

    char outputFileName[size];
    sprintf(outputFileName, outputFileNameFormat.c_str(), CTpos);
    normalizedFrame -> writeToFile_ASCIIMatrix(outputFileName);
    delete normalizedFrame;
  }

  delete OBframe;

  return 0;
}

void readConfigurationFile(char *configurationFileName, string &fileNameFormat, string &OBfileName, int &CTposMin, int &CTposMax, int &CTposStep, int &nx, int &ny, string &outputFileNameFormat){

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

  // format of the output frames files
  getline(configurationFile, line);
  outputFileNameFormat = line;

  configurationFile.close();

  return ;
}
