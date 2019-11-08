#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
using namespace std;

void applyStyle(){
  gStyle -> SetPaperSize(20, 20);
  gStyle -> SetPadTopMargin(0.2);
  gStyle -> SetPadRightMargin(0.2);
  gStyle -> SetPadBottomMargin(0.16);
  gStyle -> SetPadLeftMargin(0.16);
  gStyle -> SetPalette(51);
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1111);
  return ;
}

void plotSinogram(const char *fileName, int size, double max=-1){

  ////////////////////
  // getting angles //
  ////////////////////
  ifstream file;
  file.open(fileName);
  if(file == NULL){
    cout << " - plotSinogram(): cannot open file " << fileName << endl;
    return ;
  }
  string angleLine;
  getline(file, angleLine);
  istringstream iss(angleLine);
  vector<double> angles;
  while(iss){
    string sub;
    iss >> sub;
    double ang = atof(sub.c_str());
    angles.push_back(ang);
  }
  file.close();
  int nAngles = angles.size() - 2;
  double angle[nAngles];
  for(int i=0; i<nAngles; i++){
    angle[i] = angles[i+1];
  }

  double projection[nAngles][size];

  /////////////////////////
  // getting projections //
  /////////////////////////
  file.open(fileName);
  if(file == NULL){
    cout << " - plotSinogram(): cannot open file " << fileName << endl;
    return ;
  }
  string dummy;
  for(int i=0; i<nAngles+1; i++){
    file >> dummy;
  }

  for(int i=0; i<size; i++){
    file >> dummy;
    for(int j=0; j<nAngles; j++){
      file >> dummy;
      double ang = atof(dummy.c_str());    
      projection[j][i] = ang;
    }
  }
  file.close();

  //////////////
  // plotting //
  //////////////
  applyStyle();

  TH2F *h_sinogram = new TH2F("h_sinogram", "", nAngles, angle[0], angle[nAngles-1], size, 0, size);
  for(int i=0; i<nAngles; i++){
    for(int j=0; j<size; j++){
      h_sinogram -> SetBinContent(i+1, j+1, projection[i][j]);
    }
  }

  TCanvas *c_sinogram = new TCanvas("c_sinogram", "c_sinogram", 0, 0, 700, 700);
  //  c_sinogram -> SetLogz();
  if(max != -1.) h_sinogram -> SetMaximum(max);
  h_sinogram -> GetXaxis() -> SetTitle("Angle [deg]");
  h_sinogram -> GetXaxis() -> SetTitleSize(0.05);
  h_sinogram -> GetYaxis() -> SetTitle("Projection coordinate [pixel]");
  h_sinogram -> GetYaxis() -> SetTitleSize(0.05);
  h_sinogram -> GetZaxis() -> SetTitle("-log#left(#frac{N}{N_{0}}#right)");
  h_sinogram -> GetZaxis() -> SetTitleSize(0.05);
  h_sinogram -> Draw("colz");

  return ;
}

void plotProjectionFromSinogram(const char *fileName, int size, int projectionIndex){

  ////////////////////
  // getting angles //
  ////////////////////
  ifstream file;
  file.open(fileName);
  if(file == NULL){
    cout << " - plotProjectionFromSinogram(): cannot open file " << fileName << endl;
    return ;
  }
  string angleLine;
  getline(file, angleLine);
  istringstream iss(angleLine);
  vector<double> angles;
  while(iss){
    string sub;
    iss >> sub;
    double ang = atof(sub.c_str());
    angles.push_back(ang);
  }
  file.close();
  int nAngles = angles.size() - 2;
  double angle[nAngles];
  for(int i=0; i<nAngles; i++){
    angle[i] = angles[i+1];
  }

  if(projectionIndex >= nAngles || projectionIndex <0){
    cout << " - ERROR!!! - plotProjectionFromSinogram(): projectionIndex (" << projectionIndex << ") out of boundaries " << nAngles << endl;
    return ;
  }

  double projection[nAngles][size];

  /////////////////////////
  // getting projections //
  /////////////////////////
  file.open(fileName);
  if(file == NULL){
    cout << " - plotProjectionFromSinogram(): cannot open file " << fileName << endl;
    return ;
  }
  string dummy;
  for(int i=0; i<nAngles+1; i++){
    file >> dummy;
  }

  for(int i=0; i<size; i++){
    file >> dummy;
    for(int j=0; j<nAngles; j++){
      file >> dummy;
      double ang = atof(dummy.c_str());    
      projection[j][i] = ang;
    }
  }
  file.close();

  //////////////
  // plotting //
  //////////////
  applyStyle();

  TH1F *h_projection = new TH1F("h_projection", "", size, 0, size);
  for(int i=0; i<size; i++){
    h_projection -> SetBinContent(i+1, projection[projectionIndex][i]);
  }

  TCanvas *c_projection = new TCanvas("c_projection", "c_projection", 0, 0, 700, 700);
  c_projection -> SetLogz();
  h_projection -> GetXaxis() -> SetTitle("Projection coordinate [pixel]");
  h_projection -> GetXaxis() -> SetTitleSize(0.05);
  h_projection -> GetYaxis() -> SetTitle("-log#left(#frac{I}{I_{0}}#right)");
  h_projection -> GetYaxis() -> SetTitleSize(0.05);
  h_projection -> SetLineWidth(2);
  h_projection -> Draw("l");

  return ;
}

void plotProjectionFromFile(const char *fileName){

  /////////////////////
  // collecting data //
  /////////////////////
  ifstream file;
  file.open(fileName);
  if(file == NULL){
    cout << " - ERROR!!! - plotProjectionFromFile(): cannot open file " << fileName << endl;
    return ;
  }
  vector<double> values;
  double val;
  while(file >> val){
    values.push_back(val);
  }
  file.close();

  //////////////
  // plotting //
  //////////////
  applyStyle();

  TH1F *h_projection = new TH1F("h_projection", "", values.size(), 0, values.size());
  for(unsigned int i=0; i<values.size(); i++){
    h_projection -> SetBinContent(i+1, values[i]);
  }

  TCanvas *c_projection = new TCanvas("c_projection", "c_projection", 0, 0, 700, 700);
  //  c_projection -> SetLogz();
  h_projection -> GetXaxis() -> SetTitle("Projection coordinate [pixel]");
  h_projection -> GetXaxis() -> SetTitleSize(0.05);
  h_projection -> GetYaxis() -> SetTitle("-log#left(#frac{I}{I_{0}}#right)");
  h_projection -> GetYaxis() -> SetTitleSize(0.05);
  h_projection -> SetLineWidth(2);
  h_projection -> Draw("l");

  return ;
}
