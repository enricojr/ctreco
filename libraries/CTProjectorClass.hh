#ifndef CTPROJECTORCLASS_CC
#define CTPROJECTORCLASS_CC

#include "frameClass.hh"
//#include "crossSection.hh"
#include <vector>
using namespace std;

////////////////////
// spectrum class //
////////////////////
class spectrumClass{
public:
  spectrumClass(); // empty constructor
  ~spectrumClass(); // empty destructor
  int readFromFile(const char *fileName); // reads spectrum from file
  void readFromFile(const char *fileName, double Emax); // reads spectrum from file up to Emax
  void readFromFile(const char *fileName, double Emin, double Emax); // reads spectrum from file from Emin up to Emax
  int writeToFile(const char *fileName); // writes spectrum to file
  void print(); // prints spectrum on terminal
  void createTriangularSpectrum(double Emax, double Estep, double N); // creates an ideal bremsstrahlung spectrum up to Emax, with step Estep, with area normalized to N
  void rescaleSpectrum(double factor); // rescales spectrum
  double integrate(double Emin, double Emax, double dE); // calculates number of photons between Emin and Emax, with step dE
  double integrate(double Emin, double Emax); // calculates number of photons between Emin and Emax
  void attenuate(int Z, double thickness); // attenuates input spectrum via photoelectric absorption
  void attenuateWax(double thickness); // attenuates input spectrum via photoelectric absorption through wax
  spectrumClass *clone(); // clones spectrum to new spectrumClass object
  void createSignalToThicknessCalibration(const char *fileName, int Z, double thickness_max, double thickness_step, double THL, double Emax, double dE); // creates file containing signal-to-thickness calibration data, for element Z
  double getValueAtEnergy(double E); // returns intensity value of bin just after energy E
  int addSpectrum(spectrumClass *spectrum, double weight); // adds spectrum to current one
  spectrumClass *getSpectumWithGivenDE(double dE); // returns a spectrum corresponding to the present one, rebinned in order to have dE as specified by the argument

  vector<double> energy; // energy array
  vector<double> intensity; // intensity array
private:
};

//////////////////////
// projection class //
//////////////////////
class projectionClass{
public:

  //////////////////////////////
  // constructor / destructor //
  //////////////////////////////
  projectionClass(int size_); // initializes size to size_
  ~projectionClass(); // deletes projection array

  ///////////////
  // functions //
  ///////////////
  void init(); // allocates projection array
  void set_to(double val); // sets all values of projection array to val
  void minusLogarithm(); // turns projection into its logarithm; negative numbers are set to zero
  void getSliceFromFrame(frameClass *frame, int nx, int ny, int selectedSlice, int sliceDirection, int lowerRange); // gets slice from a specific row/column of the given frame

  ///////////////////////
  // getters / setters //
  ///////////////////////
  void set_point(int i, double val); // set value of point i of projection array to val
  double get_point(int i) const; // returns value of point i of projection array
  double *get_projection() const; // returns pointer to projection array
  int get_size() const; // returns value of size

  /////////
  // I/O //
  /////////
  void writeToFile(const char *fileName); // writes projection to file

private:
  double default_int; // default value for int variables
  double default_double; // default value for double variables
  int size; // size of the projection array
  double *projection; // projection array
};

///////////////////////
// OSEM_subset class //
///////////////////////
class OSEMsubsetClass{

public:

  OSEMsubsetClass(int nProjections_, int size_);
  ~OSEMsubsetClass();
  void assignAngleIndex(int index, int angle);
  void assignProjection(int index, projectionClass *projection, int angle);
  void exportAsSinogram(const char *fileName);

  int nProjections;
  int size;
  int *angleIndex;
  projectionClass **OSEMprojection;

private:
};

////////////////////////
// CT projector class //
////////////////////////
class CTProjectorClass{

public:

  //////////////////////////////
  // constructor / destructor //
  //////////////////////////////
  CTProjectorClass(frameClass *frame_); // creates projector object assigning frame_ to frame
  ~CTProjectorClass(); // destroys objects, deletes frame

  ///////////////
  // functions //
  ///////////////
  void init(); // allocates frame
  void importFrame(frameClass *frame_); // deletes frame and copies frame_ to frame
  void importAngles(int nAngles_, double *angle_); // reinitializes angles and copies angle_ to angle
  void project_sum(); // projects frame by adding values along y direction
  void project_attenuate(double pixelSize, double Emin, double Emax, double dE); // projects frame by calculating attenuation through object. Object values must be in units of Z (atomic number), pixel size in units of um. Number of detected photons is calculated between Emin (threshold) and Emax (spectrum end), with energy step dE. Projections are normalized ar I/I_0 -> the logarithm has to be taken separately: use logarithmOfSinogram() function
  void project_attenuate_withNoise(double pixelSize, double Emin, double Emax, double dE); // same as project_attenuate(), but including random Poisson noise
  void importSinogram(const char *fileName, int projectionSize_); // imports angles and projections from sinogram file
  void selectSinogramFromFrames(vector<frameClass *> frames, int nx, int ny, double angleMin, double angleMax, double angleStep, int selectedSlice, int sliceDirection, int lowerRange, int upperRange, double rescaleAngle); // selects sinogram from a given slice out of a set of CT scan frames
  int fixRotationAxis(const char *chi2FileName, int rotationAxisRange); // looks for the rotation axis in the sinogram in a range comprised between rotationAxisRange from the center of the frame, and translates data accordingly; IMPORTANT: assumes that first and last projections are related via a 180 deg rotation!!!
  int fixRotationAxis_zeroPadding(const char *chi2FileName, int rotationAxisRange); // looks for the rotation axis in the sinogram in a range comprised between rotationAxisRange from the center of the frame, and translates data accordingly using zero padding; IMPORTANT: assumes that first and last projections are related via a 180 deg rotation!!!
  void shiftSinogram(int shift); // shifts the sinogram and reduces its size accordingly
  void shiftSinogram_zeroPadding(int shift); // shifts the sinogram and increases its size accordingly
  void signalToThicknessCalibration(const char *fileName_STT); // applies signal-to-thickness calibration to sinogram
  void BP(); // runs geometrical backprojection calculation
  void FBP(); // runs FBP calculation. WARNING: assumes that rotation axis is at (nx/2., ny/2.)
  void eliminateNegatives(); // sets to zero all negative numbers contained in the current sinogram
  void OSEM(int subsetSize, const char *chi2File_name, double chi2CUt_subsets, double chi2Cut_fullSet); // runs OSEM reconstruction, writes chi2 into file
  void OSEM_initializeFrame(); // OSEM initialization of frame
  void OSEM_projectSubset(OSEMsubsetClass *OSEMsubset); // calls OSEM_projectAtAngle for all projections of the given subset
  void OSEM_attenuateSubset(OSEMsubsetClass *OSEMsubset, double pixelSize, double Emin, double Emax, double dE); // calls OSEM_projectAtAngle for all projections of the given subset
  void OSEM_projectAtAngle(double angleValue, projectionClass *OSEM_projection); // calculates projection at angle using simple summation algorithm
  void OSEM_attenuateAtAngle(double angleValue, projectionClass *OSEM_projection, double pixelSize, double Emin, double Emax, double dE); // calculates projection at angle using simple summation algorithm
  void OSEM_calculateCorrections(OSEMsubsetClass *OSEMsubset); // calls OSEM_calculateCorrection for all projections of the given subset
  void OSEM_calculateCorrection(projectionClass *OSEM_projection, int projectionIndex); // calculates correction factors by taking the ration between the calculated OSEM_projection and the corresponding original projection
  double OSEM_calculateChi2(OSEMsubsetClass *OSEMsubset);
  frameClass *OSEM_backprojectSubset(OSEMsubsetClass *OSEMsubset); // calculates frame from backprojection of OSEMsubset
  void OSEM_applyCorrections(frameClass *OSEM_BPFrame); // apply corrections to current frame
  void renormalizeFrame(); // renormalizes frame using frameClass::normalizeToMaximum()
  void logarithmOfSinogram(); // turns the current sinogram into its logarithm
  void computeAverageFromCircle(double radius, double &mean, double &RMS); // calls omonimous function from frameClass
  double get_minimumCTnumber(); // calls get_minimum() from frameClass
  double get_maximumCTnumber(); // calls get_maximum() from frameClass
  void copyProjectionsFromOtherProjector(CTProjectorClass *otherProjector); // copies projections from another projector
  double getAverageOfSinogram(); // returns the average value of the sinogram
  double getAverageOfSinogramInRange(double min, double max); // returns the average value of the sinogram, counting only those values in range
  void correctQuadCrossArtifact(int position); // splits pixels at position and at position+1 into 3 pixels

  ///////////
  // I / O //
  ///////////
  void writeFrameToFile(const char *fileName); // writes frame to file in ASCII matrix format
  void writeSinogram(const char *fileName); // writes sinogram to file
  void readSpectrumFromFile(const char *fileName); // initializes spectrum object and reads spectrum from file

  ///////////////////////
  // getters / setters //
  ///////////////////////
  frameClass *getFrame() const; // returns pointer to frame
  int get_nAngles() const; // returns nAngles
  double get_projectionValue(int projectionIndex, int positionIndex) const; // returns value of projectionIndex projection at position positionIndex
  void set_projectionValue(int projectionIndex, int positionIndex, double value); // sets value of projectionIndex projection at position positionIndex
  int get_projectionSize() const; // returns projectionSize
  int get_nx() const; // returns nx
  int get_ny() const; // returns ny
  double *get_angle() const; // returns pointer to angle array
  spectrumClass *get_spectrum() const; // returns pointer to spectrum

private:
  int default_int; // default value for integer variables
  double default_double; // default value for double variables
  int nx; // number of pixels x
  int ny; // number of pixels y
  frameClass *frame; // frameClass object
  int nAngles; // number of projection angles
  double *angle; // array of projection angles, in degrees
  int projectionSize; // size of projections
  projectionClass **projection; // array of projections
  spectrumClass *spectrum; // spectrum used to compute projection through attenuation
};

#endif
