#ifndef FRAMECLASS_CC
#define FRAMECLASS_CC

#include <iostream>
#include <vector>
using namespace std;

class frameClass{

public:
  //////////////////////////////
  // constructors/destructors //
  //////////////////////////////
  frameClass(); // empty constructor
  frameClass(int nx_, int ny_); // constructor with call to initialize() function
  ~frameClass(); // destructor

  /////////////////////
  // getters/setters //
  /////////////////////
  void set_nx(int nx_); // sets nx
  void set_ny(int ny_); // sets ny
  int get_nx() const; // gets nx
  int get_ny() const; // gets ny
  void set_to(double val); // sets full frame to value val
  void set_point(int x, int y, double val); // sets point (x,y) to value val
  double get_point(int x, int y) const; // gets value of point (x,y)
  //double ** get_frame(); // returns pointer to frame2D
  double *  get_frame(); // returns pointer to frame
  bool is_initialized(); // returns true if frame is correctly initialized
  void get_row(int index, double *row, int arraySize); // writes row on array
  void get_column(int index, double *column, int arraySize); // writes column on array
  int get_default_double(){return default_double;};

  ///////////////
  // functions //
  ///////////////
  void initialize(); // allocates memory for the frame
  void multiplyBy(double val); // multiplyes all frame by factor val
  void addValue(double val); // adds value val to all pixels
  void addFrame(frameClass *other); // adds other frame to current frame
  void normalizeToMaximum(); // renormalizes frame in order for maximum value to be set to 1
  frameClass *rebin(int nBinX, int nBinY); // rebins frame by grouping nBinX pixels in x and nBinY in y
  frameClass *addRow(int row); // adds row at position row
  frameClass *addColumn(int column); // adds column at position column
  void interpolateBadPixels(); // interpolates pixels with value default_double
  void interpolateBadPixels(double bad_pixel_value); // interpolates pixels with value bad_pixel_value
  void interpolateBadPixelsSolid(); // interpolates pixels with value default_double by most present value in the neighborhood
  frameClass *rotate(double angle); // rotates frame by angle, interpolation done by averaging
  frameClass *rotateSolid(double angle); // rotates frame by angle, interpolation done by most present value in the neighborhood
  void rotateCurrent(double angle); // rotates current frame by angle
  frameClass *squarize(); // creates a square frame out of a rectangular frame
  frameClass *makeCircle_in(double valueOutsideRadius); // makes a circular frame within original square frame (data loss outside radius)
  frameClass *makeCircle_out(double valueOutsideRadius); // makes a circular frame containing original square frame (new data added ouside original frame)
  frameClass *selection(int bottomLeftX, int bottomLeftY, int topRightX, int topRightY); // creates a new frame from a selection of the current frame
  frameClass *clone(); // clones current frame
  frameClass *flipHorizontally(); // flips frame around the vertical axis
  frameClass *flipVertically(); // flips frame around the horizontal axis
  void copyFromFrame(frameClass *other); // copies this frame from other frame
  void computeAverageFromCircle(double radius, double &mean, double &RMS); // returns the mean and RMS values from a circle of radius radius in the center of the frame
  double get_minimum(); // returns minimum value in frame
  double get_maximum(); // returns maximum value in frame
  void setToZeroBelow(double threshold); // sets element to zero if below threshold
  void binarize(double threshold); // sets to zero elements below threshold, to one elements above
  void subtractFrames(frameClass *frame1, frameClass *frame2); // current frame is set as frame1-frame2
  void sumFrames(frameClass *frame1, frameClass *frame2); // current frame is set as frame1+frame2
  void divideFrames(frameClass *frame1, frameClass *frame2); // current frame is set as frame1/frame2; if frame2 is zero, value is set to zero
  void multiplyFrames(frameClass *frame1, frameClass *frame2); // current frame is set as frame1*frame2
  void minusLog(); // each value is substituted by its -log, with particular attention to 0 and negative values
  void multiplyByFrame(frameClass *frame); // multiplies current frame by given frame
  void getStats(double &mean, double &rms, int &entries); // computes mean and rms
  void getStatsInRange(double &mean, double &rms, int &entries, double min, double max); // computes mean and rms in specified range
  void setTo_inCircle(double radius, double val); // sets to val all pixels in circle
  frameClass *applyGaussFilter(int sigma); // applies Gaussian filter, sigma is given in units of pixels; does not work on images with negative values

  /////////
  // I/O //
  /////////
  bool readFromFile_ASCIIMatrix(const char* fileName); // reads frame from ASCII matrix file
  void writeToFile_ASCIIMatrix(const char* fileName) const; // writes frame to ASCII matrix file
  void readFromFile_ASCIITable(const char* fileName); // reads frame from XYC ASCII table
  void writeToFile_ASCIITable(const char* fileName) const; // writes frame to XYC ASCII table
  void readFromFile_binaryTable(const char* fileName); // reads frame from binary table file

private:
  int default_int; // default value for integer variables
  double default_double; // default value for double variables
  int _nx; // number of pixels x
  int _ny; // number of pixels y
  double ** _frame2D; // for compatibility
  double *  _frame;	   // frame pointer
  bool initialized; // is set to true if frame is correctly initialized

  double getMostOccurentValue(vector<double> array);
};


#endif
