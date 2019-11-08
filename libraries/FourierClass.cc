#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>
#include <fftw3.h>
#include "FourierClass.hh"
using namespace std;

//////////////////////////
//////////////////////////
// 1D Fourier transform //
//////////////////////////
//////////////////////////
FourierClass_1D::FourierClass_1D(int nData_){ 
  nData = nData_;
  dataIn = NULL;
  dataOut = NULL;
  dataBack = NULL;
  dataModule = NULL;
  dataPhase = NULL;
}

FourierClass_1D::~FourierClass_1D(){
  if(dataIn != NULL) fftw_free(dataIn);
  if(dataOut != NULL) fftw_free(dataOut);
  if(dataBack != NULL) fftw_free(dataBack);
  if(dataModule != NULL) delete dataModule;
  if(dataPhase != NULL) delete dataPhase;
};

void FourierClass_1D::init(){
  if(dataIn == NULL) dataIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nData);
  return ;
}

void FourierClass_1D::getRealData(int nRealData, double *realData){
  if(nRealData != nData){
    cout << " - ERROR!!! - FourierClass_1D::getRealData(): incompatible data size" << endl;
    return ;
  }
  init();
  for(int i=0; i<nData; i++){
    dataIn[i][0] = realData[i];
    dataIn[i][1] = 0.;
  }
  return ;
}

void FourierClass_1D::logarithm(){
  for(int i=0; i<nData; i++){
    dataIn[i][0] = log(dataIn[i][0]);
    dataIn[i][1] = log(dataIn[i][1]);
  }
  return ;
}

void FourierClass_1D::fft(){
  if(dataOut == NULL) dataOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nData);
  fftw_plan plan = fftw_plan_dft_1d(nData, dataIn, dataOut, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan); 
  fftw_destroy_plan(plan);
  for(int i=0; i<nData; i++){
    dataOut[i][0] /= nData;
    dataOut[i][1] /= nData;
  }
  return ;
}

void FourierClass_1D::rfft(){
  if(dataBack == NULL) dataBack = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nData);
  fftw_plan plan = fftw_plan_dft_1d(nData, dataOut, dataBack, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan); 
  fftw_destroy_plan(plan);
  return ;
}

void FourierClass_1D::module(){
  if(dataModule == NULL) dataModule = new double[nData];
  for(int i=0; i<nData; i++){
    dataModule[i] = sqrt(dataOut[i][0] * dataOut[i][0] + dataOut[i][1] * dataOut[i][1]);
  }
  return ;
}

void FourierClass_1D::phase(){
  if(dataPhase == NULL) dataPhase = new double[nData];
  for(int i=0; i<nData; i++){
    dataPhase[i] = atan(dataOut[i][1] / dataOut[i][0]);
  }
  return ;
}

void FourierClass_1D::writeIn(const char *fileName){
  ofstream file;
  file.open(fileName);
  if(file == NULL){
    cout << " - ERROR!!! - FourierClass_1D::writeIn(): cannot open file " << fileName << endl;
    return ;
  }
  for(int i=0; i<nData; i++){
    file << i << "\t" << dataIn[i][0] << "\t" << dataIn[i][1] << endl;
  }
  file.close();
  return ;
}

void FourierClass_1D::writeOut(const char *fileName){
  ofstream file;
  file.open(fileName);
  if(file == NULL){
    cout << " - ERROR!!! - FourierClass_1D::writeOut(): cannot open file " << fileName << endl;
    return ;
  }
  for(int i=0; i<nData; i++){
    file << i << "\t" << dataOut[i][0] << "\t" << dataOut[i][1] << endl;
  }
  file.close();
  return ;
}

void FourierClass_1D::writeBack(const char *fileName){
  ofstream file;
  file.open(fileName);
  if(file == NULL){
    cout << " - ERROR!!! - FourierClass_1D::writeBack(): cannot open file " << fileName << endl;
    return ;
  }
  for(int i=0; i<nData; i++){
    file << i << "\t" << dataBack[i][0] << "\t" << dataBack[i][1] << endl;
  }
  file.close();
  return ;
}


void FourierClass_1D::writeModule(const char *fileName){
  ofstream file;
  file.open(fileName);
  if(file == NULL){
    cout << " - ERROR!!! - FourierClass_1D::writeModule(): cannot open file " << fileName << endl;
    return ;
  }
  for(int i=0; i<nData; i++){
    file << i << "\t" << dataModule[i] << endl;
  }
  file.close();
  return ;
}

void FourierClass_1D::writePhase(const char *fileName){
  ofstream file;
  file.open(fileName);
  if(file == NULL){
    cout << " - ERROR!!! - FourierClass_1D::writePhase(): cannot open file " << fileName << endl;
    return ;
  }
  for(int i=0; i<nData; i++){
    file << i << "\t" << dataPhase[i] << endl;
  }
  file.close();
  return ;
}

void FourierClass_1D::filterDataOut(){
  if(dataOut == NULL){
    cout << " - ERROR!!! - FourierClass_1D::filterDataOut(): dataOut is not initialized" << endl;
    return ;
  }
  for(int i=1; i<nData; i++){
    double frequency = -1.;
    if(i < nData / 2.) frequency = (double) i;
    else frequency = (double) (nData - i);
    dataOut[i][0] *= fabs(frequency);
    dataOut[i][1] *= fabs(frequency);
  }
  return ;
}

//////////////////////////
//////////////////////////
// 2D Fourier transform //
//////////////////////////
//////////////////////////
/*
FourierClass_2D::FourierClass_2D(int nDataX_, int nDataY_){ 
  nDataX = nDataX_;
  nDataY = nDataY_;
  dataIn = NULL;
  dataOut = NULL;
  dataBack = NULL;
}

FourierClass_2D::~FourierClass_2D(){
  if(dataIn != NULL) fftw_free(dataIn);
  if(dataOut != NULL) fftw_free(dataOut);
  if(dataBack != NULL) fftw_free(dataBack);
};

void FourierClass_2D::init(){
  if(dataIn == NULL) dataIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nDataX * nDataY);
}

void FourierClass_2D::fft(){
  if(dataOut == NULL) dataOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nDataX * nDataY);
  fftw_plan plan = fftw_plan_dft_1d(nDataX, nDataY, dataIn, dataOut, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan); 
  fftw_destroy_plan(plan);
}

void FourierClass_2D::rfft(){
  if(dataBack == NULL) dataBack = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nDataX * nDataY);
  fftw_plan plan = fftw_plan_dft_1d(nDataX, nDataY, dataOut, dataBack, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan); 
  fftw_destroy_plan(plan);
}

void FourierClass_2D::writeIn(const char *fileName){
  ofstream file;
  file.open(fileName);
  for(int i=0; i<nData; i++){
    file << i << "\t" << dataIn[i][0] << "\t" << dataIn[i][1] << endl;
  }
  file.close();
}

void FourierClass_2D::writeOut(const char *fileName){
  ofstream file;
  file.open(fileName);
  for(int i=0; i<nData; i++){
    file << i << "\t" << dataOut[i][0] << "\t" << dataOut[i][1] << endl;
  }
  file.close();
}

void FourierClass_2D::writeBack(const char *fileName){
  ofstream file;
  file.open(fileName);
  for(int i=0; i<nData; i++){
    file << i << "\t" << dataBack[i][0] << "\t" << dataBack[i][1] << endl;
  }
  file.close();
}
*/
