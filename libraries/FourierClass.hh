#ifndef FOURIERCLASS_HH
#define FOURIERCLASS_HH

#include <fftw3.h>

class FourierClass_1D{
public:

  FourierClass_1D(int nData_);
  ~FourierClass_1D();
  void init();
  void getRealData(int nRealData, double *realData);
  void logarithm();
  void fft();
  void rfft();
  void module();
  void phase();
  void writeIn(const char *fileName);
  void writeOut(const char *fileName);
  void writeBack(const char *fileName);
  void writeModule(const char *fileName);
  void writePhase(const char *fileName);
  void filterDataOut();

  int nData;
  fftw_complex *dataIn;
  fftw_complex *dataOut;
  fftw_complex *dataBack;
  double *dataModule;
  double *dataPhase;
};

/*
class FourierClass_2D{
public:

  FourierClass_2D(int nDataX_, int nDataY_);
  ~FourierClass_2D();
  void init();
  void fft();
  void rfft();
  void writeIn(const char *fileName);
  void writeOut(const char *fileName);
  void writeBack(const char *fileName);

  int nDataX;
  int nDataY;
  fftw_complex *dataIn;
  fftw_complex *dataOut;
  fftw_complex *dataBack;
};
*/

#endif
