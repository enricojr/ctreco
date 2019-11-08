#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include "libraries/frameClass.hh"
using namespace std;

void readConfigurationFile(char *configurationFileName, string &fileNameFormat, int &CTposMin, int &CTposMax, int &CTposStep, int &frameMin, int &frameMax, int &frameStep, int &nx, int &ny, string &outputFileNameFormat);

int main(int argc, char *argv[]){

  if(argc < 2){
    cout << " - " << argv[0] << " usage: " << endl;
    cout << "   " << argv[0] << " <configuration file name>" << endl;
    return 0;
  }

  char *configurationFileName = argv[1];
  string fileNameFormat = "";
  int CTposMin = -1;
  int CTposMax = -1;
  int CTposStep = -1;
  int frameMin = -1;
  int frameMax = -1;
  int frameStep = -1;
  int nx = -1;
  int ny = -1;
  string outputFileNameFormat = "";

  readConfigurationFile(configurationFileName, fileNameFormat, CTposMin, CTposMax, CTposStep, frameMin, frameMax, frameStep, nx, ny, outputFileNameFormat);

  cout << " - fileNameFormat = " << fileNameFormat << endl;
  cout << " - CTposMin = " << CTposMin << endl;
  cout << " - CTposMax = " << CTposMax << endl;
  cout << " - CTposStep = " << CTposStep << endl;
  cout << " - frameMin = " << frameMin << endl;
  cout << " - frameMax = " << frameMax << endl;
  cout << " - frameMStep = " << frameStep << endl;
  cout << " - nx = " << nx << endl;
  cout << " - ny = " << ny << endl;
  cout << " - outputFileNameFormat = " << outputFileNameFormat << endl;

  int size = 1000;
  for(int CTpos=CTposMin; CTpos<=CTposMax; CTpos+=CTposStep){
    cout << " - processing CTpos = " << CTpos << endl;
    frameClass *summedFrame = new frameClass(nx, ny);
    summedFrame -> set_to(0.);
    for(int frameNumber=frameMin; frameNumber<=frameMax; frameNumber+=frameStep){
      char fileName[size];
      sprintf(fileName, fileNameFormat.c_str(), CTpos, frameNumber);
      frameClass *frame = new frameClass(nx, ny);
      frame -> readFromFile_ASCIIMatrix(fileName);
      summedFrame -> addFrame(frame);
      delete frame;
    }
    char outputFileName[size];
    sprintf(outputFileName, outputFileNameFormat.c_str(), CTpos);
    summedFrame -> writeToFile_ASCIIMatrix(outputFileName);
    delete summedFrame;
  }

  return 0;
}

void readConfigurationFile(char *configurationFileName, string &fileNameFormat, int &CTposMin, int &CTposMax, int &CTposStep, int &frameMin, int &frameMax, int &frameStep, int &nx, int &ny, string &outputFileNameFormat){

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

  // minimum angle position
  getline(configurationFile, line);
  CTposMin = atoi(line.c_str());

  // maximum angle position
  getline(configurationFile, line);
  CTposMax = atoi(line.c_str());

  // angle step
  getline(configurationFile, line);
  CTposStep = atoi(line.c_str());

  // minimum frame number
  getline(configurationFile, line);
  frameMin = atoi(line.c_str());

  // maximum frame number
  getline(configurationFile, line);
  frameMax = atoi(line.c_str());

  // frame number step
  getline(configurationFile, line);
  frameStep = atoi(line.c_str());

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
