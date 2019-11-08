#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include "libraries/CTProjectorClass.hh"
#include "libraries/FourierClass.hh"
#include "libraries/frameClass.hh"
using namespace std;

/**************************/
/**************************/
/* functions declarations */
/**************************/
/**************************/
void getParameters(const char *configurationFileName, string &fileNameFormat, int &nx, int &ny, double &angleMin, double &angleMax, double &angleStep, int &startSlice, int &endSlice, int &sliceDirection, string &sinogramFileNameFormat, int &lowerRange, int &upperRange, double &rescaleAngle, int &correctQuadCross);
void printParameters(const char *configurationFileName, string fileNameFormat, int nx, int ny, double angleMin, double angleMax, double angleStep, int startSlice, int endSlice, int sliceDirection, const char *sinogramFileNameFormat, int lowerRange, int upperRange, double rescaleAngle, int correctQuadCross);
void makeSinogram(const char *configurationFileName, vector<frameClass *> frames, int nx, int ny, double angleMin, double angleMax, double angleStep, int selectedSlice, int sliceDirection, char *sinogramFileName, int lowerRange, int upperRange, double rescaleAngle, int correctQuadCross);
vector<frameClass *> loadFrames(string fileNameFormat, int nx, int ny, double angleMin, double angleMax, double angleStep);

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

  //int time0 = time(0);

  ////////////////////////
  // parameters setting //
  ////////////////////////
  int default_int = -777;
  int default_double = -777.888;
  string fileNameFormat = "EMPTY";
  int nx = default_int;
  int ny = default_int;
  double angleMin = default_double;
  double angleMax = default_double;
  double angleStep = default_double;
  int startSlice = default_int;
  int endSlice = default_int;
  int sliceDirection = default_int; // 0 x, 1 y
  string sinogramFileNameFormat = "EMPTY";
  int lowerRange = default_int;
  int upperRange = default_int;
  double rescaleAngle = default_double;
  int correctQuadCross = default_int;

  ////////////////////////////////////////////////
  // getting parameters from configuration file //
  ////////////////////////////////////////////////
  getParameters(argv[1], fileNameFormat, nx, ny, angleMin, angleMax, angleStep, startSlice, endSlice, sliceDirection, sinogramFileNameFormat, lowerRange, upperRange, rescaleAngle, correctQuadCross);
  printParameters(argv[1], fileNameFormat, nx, ny, angleMin, angleMax, angleStep, startSlice, endSlice, sliceDirection, sinogramFileNameFormat.c_str(), lowerRange, upperRange, rescaleAngle, correctQuadCross);

  ////////////////////
  // loading frames //
  ////////////////////
  vector<frameClass *> frames = loadFrames(fileNameFormat, nx, ny, angleMin, angleMax, angleStep);

  //////////////////////////////////////
  // getting data and making sinogram //
  //////////////////////////////////////
  for(int selectedSlice = startSlice; selectedSlice <= endSlice; selectedSlice++){
    int size = 1000;
    char sinogramFileName[size];
    sprintf(sinogramFileName, sinogramFileNameFormat.c_str(), selectedSlice);
    cout << " - processing slice " << selectedSlice << ", file name " << sinogramFileName << endl;
    makeSinogram(argv[1], frames, nx, ny, angleMin, angleMax, angleStep, selectedSlice, sliceDirection, sinogramFileName, lowerRange, upperRange, rescaleAngle, correctQuadCross);
  }

  //////////////
  // deleting //
  //////////////
  for(unsigned int i=0; i<frames.size(); i++){
    delete frames[i];
  }
}

/*************************/
/*************************/
/* functions definitions */
/*************************/
/*************************/
void getParameters(const char *configurationFileName, string &fileNameFormat, int &nx, int &ny, double &angleMin, double &angleMax, double &angleStep, int &startSlice, int &endSlice, int &sliceDirection, string &sinogramFileNameFormat, int &lowerRange, int &upperRange, double &rescaleAngle, int &correctQuadCross){

  ifstream configurationFile;
  configurationFile.open(configurationFileName);
  if(configurationFile == NULL){
    cout << " - ERROR!!! - getParameters(): cannot open file " << configurationFileName << endl;
    return ;
  }

  string line;

  // file name format of input raw frames
  getline(configurationFile, line);
  fileNameFormat = line.c_str();

  // number of pixels nx
  getline(configurationFile, line);
  nx = atoi(line.c_str());

  // number of pixels ny
  getline(configurationFile, line);
  ny = atoi(line.c_str());

  // minimum angle position
  getline(configurationFile, line);
  angleMin = atof(line.c_str());

  // maximum angle position
  getline(configurationFile, line);
  angleMax = atof(line.c_str());

  // angle step
  getline(configurationFile, line);
  angleStep = atof(line.c_str());

  // starting slice
  getline(configurationFile, line);
  startSlice = atoi(line.c_str());

  // ending slice
  getline(configurationFile, line);
  endSlice = atoi(line.c_str());

  // slice direction: 0=x, 1=y
  getline(configurationFile, line);
  sliceDirection = atoi(line.c_str());

  // format of the output sinogram file names
  getline(configurationFile, sinogramFileNameFormat);

  // data range: minimum pixel row (column)
  getline(configurationFile, line);
  lowerRange = atoi(line.c_str());

  // data range: maximum pixel row (column)
  getline(configurationFile, line);
  upperRange = atoi(line.c_str());

  // angle rescale factor
  getline(configurationFile, line);
  rescaleAngle = atof(line.c_str());

  // quad cross artifact correction
  getline(configurationFile, line);
  correctQuadCross = atoi(line.c_str());

  configurationFile.close();

  return ;
}

void printParameters(const char *configurationFileName, string fileNameFormat, int nx, int ny, double angleMin, double angleMax, double angleStep, int startSlice, int endSlice, int sliceDirection, const char *sinogramFileNameFormat, int lowerRange, int upperRange, double rescaleAngle, int correctQuadCross){

  cout << " - printParameters(): parameter values from " << configurationFileName << endl;
  cout << "\t-  format of data files names = " << fileNameFormat << endl;
  cout << "\t-                         nx = " << nx << endl;
  cout << "\t-                          ny = " << ny << endl;
  cout << "\t-                   angle min = " << angleMin << endl;
  cout << "\t-                   angle max = " << angleMax << endl;
  cout << "\t-                  angle step = " << angleStep << endl;
  cout << "\t-                 start slice = " << startSlice << endl;
  cout << "\t-                   end slice = " << endSlice << endl;
  cout << "\t-             slice direction =";
  if(sliceDirection == 0) cout << " x " << endl;
  else if(sliceDirection == 1) cout << " y " << endl;
  else{
    cout << " undefined!" << endl;
    cout << " - ERROR!!! - printParameters(): slice direction is undefined, must be either 0 (for x) or 1 (for y)" << endl;
    return ;
  }
  cout << "\t-          sinogram file name = " << sinogramFileNameFormat << endl;
  cout << "\t-           lower slice range = " << lowerRange << endl;
  cout << "\t-           upper slice range = " << upperRange << endl;
  cout << "\t-               rescale angle = " << rescaleAngle << endl;
  cout << "\t- correct quad cross artifact =";
  if(correctQuadCross == 0) cout << " NO" << endl;
  else cout << " YES, at position " << correctQuadCross << endl;

  return ;
}

void makeSinogram(const char * /*configurationFileName*/, vector<frameClass *> frames, int nx, int ny, double angleMin, double angleMax, double angleStep, int selectedSlice, int sliceDirection, char *sinogramFileName, int lowerRange, int upperRange, double rescaleAngle, int correctQuadCross){

  frameClass *frame = new frameClass(nx, ny);
  frame -> set_to(0.);
  CTProjectorClass *projector = new CTProjectorClass(frame);
  delete frame;
  projector -> selectSinogramFromFrames(frames, nx, ny, angleMin, angleMax, angleStep, selectedSlice, sliceDirection, lowerRange, upperRange, rescaleAngle);
  if(correctQuadCross != 0){
    projector -> correctQuadCrossArtifact(correctQuadCross);
    cout << " - makeSinogram(): new size after correction for quad cross artifact = " << projector -> get_projectionSize() << endl;
  }
  projector -> writeSinogram(sinogramFileName);
  cout << " - makeSinogram(): written sinogram " << sinogramFileName << endl; 
  delete projector;

  return ;
}

vector<frameClass *> loadFrames(string fileNameFormat, int nx, int ny, double angleMin, double angleMax, double angleStep){

  vector<frameClass *> frames;

  for(double angle=angleMin; angle<=angleMax; angle += angleStep){

    char fileName[1000];
    sprintf(fileName, fileNameFormat.c_str(), (int)angle);
    cout << " - loadFrames(): loading file " << fileName << endl;

    frameClass *frame = new frameClass(nx, ny);
    frame -> readFromFile_ASCIIMatrix(fileName);

    frames.push_back(frame);

  }

  return frames;
}
