#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <math.h>
#include <sstream>
#include <vector>
#include <string>
#include "frameClass.hh"
#include "FourierClass.hh"
//#include "crossSection.hh"
#include "CTProjectorClass.hh"
//#include "mucal.h"
//#include "mucal.c"
#include <random>
#include <time.h>
using namespace std;

///////////////////////////////
// spectrum class definition //
///////////////////////////////
spectrumClass::spectrumClass(){
}

spectrumClass::~spectrumClass(){
}

int spectrumClass::readFromFile(const char *fileName){
	/****************************/
	/* reads spectrum from file */
	/****************************/
	/* format:       */
	/*               */
	/* #Energy dN/dE */
	/* .       .     */
	/*****************/
	ifstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - spectrumClass::readFromFile(): cannot open file " << fileName << endl;
		return -1;
	}
	string dummy;
	file >> dummy;
	file >> dummy;
	while(file >> dummy){
		energy.push_back(atof(dummy.c_str()));
		file >> dummy;
		intensity.push_back(atof(dummy.c_str()));
	}
	file.close();

	return 0;
}

void spectrumClass::readFromFile(const char *fileName, double Emax){
	/****************************/
	/* reads spectrum from file */
	/****************************/
	/* format:       */
	/*               */
	/* #Energy dN/dE */
	/* .       .     */
	/*****************/
	ifstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - spectrumClass::readFromFile(): cannot open file " << fileName << endl;
		return ;
	}
	string dummy;
	file >> dummy;
	file >> dummy;
	while(file >> dummy){
		if(atof(dummy.c_str()) > Emax) return ;
		energy.push_back(atof(dummy.c_str()));
		file >> dummy;
		intensity.push_back(atof(dummy.c_str()));
	}
	file.close();

	return ;
}

void spectrumClass::readFromFile(const char *fileName, double Emin, double Emax){
	/****************************/
	/* reads spectrum from file */
	/****************************/
	/* format:       */
	/*               */
	/* #Energy dN/dE */
	/* .       .     */
	/*****************/
	ifstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - spectrumClass::readFromFile(): cannot open file " << fileName << endl;
		return ;
	}
	string dummy;
	file >> dummy;
	file >> dummy;
	while(file >> dummy){
		if(atof(dummy.c_str()) < Emin) continue ;
		if(atof(dummy.c_str()) > Emax) return ;
		energy.push_back(atof(dummy.c_str()));
		file >> dummy;
		intensity.push_back(atof(dummy.c_str()));
	}
	file.close();

	return ;
}

double spectrumClass::getValueAtEnergy(double E){

	for(unsigned int i=0; i<energy.size(); i++){
		if(energy[i] > E) return intensity[i];
	}

	return -777.;
}

int spectrumClass::writeToFile(const char *fileName){
	/*****************************/
	/* writes spectrum from file */
	/*****************************/
	/* format:       */
	/*               */
	/* #Energy dN/dE */
	/* .       .     */
	/*****************/
	ofstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - spectrumClass::writeToFile(): cannot open file " << fileName << endl;
		return -1;
	}
	file << "#Energy\tdN/dE" << endl;
	for(int i=0; i<(int)energy.size(); i++){
		file << energy[i] << "\t" << intensity[i] << endl;
	}
	file.close();

	return 0;
}

void spectrumClass::print(){
	/******************/
	/* couts spectrum */
	/******************/
	for(int i=0; i<(int)energy.size(); i++){
		cout << i << "\t" << energy[i] << "\t" << intensity[i] << endl;
	}
	return ;
}

void spectrumClass::createTriangularSpectrum(double Emax, double Estep, double N){

	if(energy.size() != 0){
		cout << " - ERROR - spectrumClass::createTriangularSpectrum(): spectrum already exists" << endl;
		return ;
	}

	for(double E=0.; E<Emax; E+=Estep){
		energy.push_back(E);
		intensity.push_back((2 * N / Emax) * (1 - E / Emax));
	}
	energy.push_back(Emax);
	intensity.push_back(0.);

	return ;
}

void spectrumClass::rescaleSpectrum(double factor){

	for(int i=0; i<(int)intensity.size(); i++){
		intensity[i] *= factor;
	}

	return ;
}

double spectrumClass::integrate(double Emin, double Emax, double dE){
	/*******************************************/
	/* integrates spectrum (number of photons) */
	/*******************************************/
	double integral = 0.;

	for(int i=1; i<(int)energy.size()-1; i++){

		if(energy[i-1] < Emin && energy[i] > Emin){
			integral += (intensity[i-1] * (energy[i] - Emin) / (energy[i] - energy[i-1]));
		}
		else if(energy[i] < Emax && energy[i+1] > Emax){
			integral += (intensity[i] * (Emax - energy[i]) / (energy[i+1] - energy[i]));
		}
		else if(energy[i] > Emin && energy[i] < Emax){
			integral += intensity[i];
		}

	}

	return integral * dE;
}

double spectrumClass::integrate(double Emin, double Emax){
	/*******************************************/
	/* integrates spectrum (number of photons) */
	/*******************************************/
	double integral = 0.;

	for(int i=1; i<(int)energy.size()-1; i++){

		double dE = energy[i+1] - energy[i];

		if(energy[i-1] < Emin && energy[i] > Emin){
			integral += dE * (intensity[i-1] * (energy[i] - Emin) / (energy[i] - energy[i-1]));
		}
		else if(energy[i] < Emax && energy[i+1] > Emax){
			integral += dE * (intensity[i] * (Emax - energy[i]) / (energy[i+1] - energy[i]));
		}
		else if(energy[i] > Emin && energy[i] < Emax){
			integral += dE * intensity[i];
		}

	}

	return integral;
}

void spectrumClass::attenuate(int /*Z*/, double /*thickness*/){
	//   /**********************************************************/
	//   /* attenuates input spectrum via photoelectric absorption */
	//   /**********************************************************/
	//   if(Z == 0){
	//     cout << " - ERROR!!! - spectrumClass::attenuate: Z cannot be null" << endl;
	//     return ;
	//   }

	//   elementData *element = new elementData();
	//   for(int i=0; i<energy.size(); i++){
	//     if(Z == 1) element -> collectWax(Z, energy[i]);
	//     else element -> collectElementData(Z, energy[i]);
	//     double exponent = - (element -> density) * (element -> photoelectricCrossSection) * thickness / 10000.; // factor 10000 translates um in cm
	//     if(exponent > 0.){
	//       cout << " - ERROR!!! - spectrumClass::attenuate(): positive exponent" << endl;
	//       return ;
	//     }
	//     intensity[i] = intensity[i] * exp(exponent); // factor 10000 translates um in cm
	//   }
	//   delete element;

	return ;
}

void spectrumClass::attenuateWax(double thickness){

	double density = 0.9;

	////////////////////////////////////////////////////////////////////////
	// surrogate function that fits NIST XCOM data for compound C32H64O2: //
	// sigma(E) = exp(B) * pow(E,A)                                       //
	////////////////////////////////////////////////////////////////////////
	double B = 7.8;
	double A = -3.105;

	for(int i=0; i<(int)energy.size(); i++){
		double photoelectricCrossSection = exp(B) * pow(energy[i], A);
		double exponent = - density * photoelectricCrossSection * thickness / 10000.; // factor 10000 translates um in cm
		if(exponent > 0.){
			cout << " - ERROR!!! - spectrumClass::attenuateWax(): positive exponent" << endl;
			return ;
		}
		intensity[i] = intensity[i] * exp(exponent);
	}

	return ;
}

spectrumClass *spectrumClass::clone(){

	spectrumClass *newSpectrum = new spectrumClass();
	for(int i=0; i<(int)energy.size(); i++){
		newSpectrum -> energy.push_back(energy[i]);
		newSpectrum -> intensity.push_back(intensity[i]);
	}
	return newSpectrum;

}

void spectrumClass::createSignalToThicknessCalibration(const char *fileName, int Z, double thickness_max, double thickness_step, double THL, double Emax, double dE){

	//////////////////////////////////////////////////
	// thickness_max and thickness_step in units of mm
	//////////////////////////////////////////////////

	ofstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - spectrumClass::createSignalToThicknessCalibration(): cannot open file " << fileName << endl;
		return ;
	}

	for(double thickness = 0.; thickness <= thickness_max; thickness += thickness_step){
		spectrumClass *tmpSpectrum = clone();
		tmpSpectrum -> attenuate(Z, thickness);
		file << thickness << "\t" << tmpSpectrum -> integrate(THL, Emax, dE) << endl;
		delete tmpSpectrum;
	}

	file.close();

	return ;
}

int spectrumClass::addSpectrum(spectrumClass *spectrum, double weight){

	for(unsigned int i=0; i<spectrum -> energy.size(); i++){
		for(unsigned int j=0; j<energy.size()-1; j++){
			if(spectrum -> energy[i] >= energy[j] && spectrum -> energy[i] < energy[j+1]){
				intensity[j] += weight * spectrum -> intensity[i];
				break;
			}
		}
		if(spectrum -> energy[i] >= energy[energy.size() - 1]){
			intensity[energy.size() - 1] += weight * spectrum -> intensity[i];
			break;
		}
	}

	return 0;
}

spectrumClass *spectrumClass::getSpectumWithGivenDE(double dE){

	spectrumClass *spectrum = new spectrumClass();

	for(double E=energy[0]; E<=energy[energy.size() - 1]; E+=dE){
		double Enext = E+dE;
		double totalIntensity = 0.;
		bool lastBinReached = false;
		for(unsigned int i=0; i<energy.size(); i++){
			if(energy[i] >= E && energy[i] < Enext){
				totalIntensity += intensity[i];
				if(i == energy.size() - 1) lastBinReached = true;
			}
		}
		spectrum -> energy.push_back(E);
		spectrum -> intensity.push_back(totalIntensity);
		if(lastBinReached) return spectrum;
	}

	return spectrum;
}

//////////////////////
// projection class //
//////////////////////

projectionClass::projectionClass(int size_){
	default_int = -777;
	default_double = -777.;
	projection = NULL;
	size = size_;
	init();
}

projectionClass::~projectionClass(){
	if(projection != NULL) delete [] projection;
}

void projectionClass::init(){
	projection = new double[size];
	set_to(default_double);
	return ;
}

void projectionClass::set_to(double val){
	for(int i=0; i<size; i++){
		projection[i] = val;
	}
	return ;
}

void projectionClass::minusLogarithm(){
	for(int i=0; i<size; i++){
		projection[i] = -log(projection[i]);
		if(projection[i] < 0) projection[i] = 0.;
	}
	return ;
}

void projectionClass::getSliceFromFrame(frameClass *frame, int /*nx*/, int /*ny*/, int selectedSlice, int sliceDirection, int lowerRange){

	if(sliceDirection == 0){
		for(int i=0; i<size; i++){
			set_point(i, frame -> get_point(lowerRange + i, selectedSlice));
		}
	}
	else if(sliceDirection == 1){
		for(int i=0; i<size; i++){
			set_point(i, frame -> get_point(selectedSlice, lowerRange + i));
		}
	}
	else{
		cout << " - ERROR!!! - projectionClass::getSliceFromFrames(): slice direction is undefined, must be either 0 (for x) or 1 (for y)" << endl;
		return ;
	}

	return ;
}

void projectionClass::set_point(int i, double val){
	if(i < 0 || i >= size){
		cout << " - ERROR!!! - projectionClass::set_point(): index (" << i << ") out of range (" << size << ")" << endl;
		return ;
	}
	projection[i] = val;
	return ;
}

double projectionClass::get_point(int i) const{
	if(i < 0 || i >= size){
		cout << " - ERROR!!! - projectionClass::get_point(): index (" << i << ") out of range (" << size << ")" << endl;
		return default_double;
	}
	return projection[i];
}

double *projectionClass::get_projection() const{
	return projection;
}

int projectionClass::get_size() const{
	return size;
}

void projectionClass::writeToFile(const char *fileName){

	ofstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " -ERROR!!! - projectionClass::writeToFile(): cannot open file " << fileName << endl;
		return ;
	}
	for(int i=0; i<size; i++){
		file << projection[i] << endl;
	}
	file.close();

	return ;
}

///////////////////////
// OSEM subset class //
///////////////////////

OSEMsubsetClass::OSEMsubsetClass(int nProjections_, int size_){
	nProjections = nProjections_;
	size = size_;
	angleIndex = new int[nProjections];
	OSEMprojection = new projectionClass*[nProjections];
	for(int i=0; i<nProjections; i++){
		OSEMprojection[i] =  new projectionClass(size);
	}
}

OSEMsubsetClass::~OSEMsubsetClass(){
	delete [] angleIndex;
	for(int i=0; i<nProjections; i++){
		delete OSEMprojection[i];
	}
	delete [] OSEMprojection;
}

void OSEMsubsetClass::assignAngleIndex(int index, int angle){
	if(index >= nProjections){
		cout << " - ERROR!!! - OSEMsubsetClass::assignAngleIndex(): projection index out of bounds" << endl;
		return ;
	}
	angleIndex[index] = angle;
	return ;
}

void OSEMsubsetClass::assignProjection(int index, projectionClass *projection, int angle){
	if(projection -> get_size() != size){
		cout << " - ERROR!!! - OSEMsubsetClass::assignProjection(): incompatible sizes" << endl;
		return ;
	}
	if(index >= nProjections){
		cout << " - ERROR!!! - OSEMsubsetClass::assignProjection(): projection index out of bounds" << endl;
		return ;
	}
	angleIndex[index] = angle;
	for(int i=0; i<size; i++){
		OSEMprojection[index] -> set_point(i, projection -> get_point(i));
	}
	return ;
}

void OSEMsubsetClass::exportAsSinogram(const char *fileName){

	ofstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - OSEMsubsetClass::exportAsSinogram(): cannot open file " << fileName << endl;
		return ;
	}
	file << "#\t";
	for(int j=0; j<nProjections; j++){
		file << angleIndex[j] << "\t";
	}
	file << endl;
	for(int i=0; i<size; i++){
		file << i << "\t";
		for(int j=0; j<nProjections; j++){
			file << OSEMprojection[j] -> get_point(i) << "\t";
		}
		file << endl;
	}
	file.close();

	return ;
}

////////////////////////
// CT projector class //
////////////////////////

CTProjectorClass::CTProjectorClass(frameClass *frame_){
	default_int = -777;
	default_double = -777.;
	nx = frame_ -> get_nx();
	ny = frame_ -> get_ny();
	nAngles = default_int;
	if(nx < ny) projectionSize = nx;
	else projectionSize = ny;
	angle = NULL;
	spectrum = NULL;
	projection = NULL;
	init();
	importFrame(frame_);
}

CTProjectorClass::~CTProjectorClass(){
	if(frame != NULL) delete frame;
	if(spectrum != NULL) delete spectrum;
	if(projection != NULL){
		for(int i=0; i<nAngles; i++){
			delete projection[i];
		}
	}
}

void CTProjectorClass::init(){
	frame = new frameClass(nx, ny);
	frame -> set_to(default_double);
	return ;
}

void CTProjectorClass::importFrame(frameClass *frame_){
	if(nx != frame_ -> get_nx() || ny != frame_ -> get_ny()){
		cout << " - ERROR!!! - CTProjectorClass::importFrame(): incompatible frame sizes" << endl;
		return ;
	}
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			frame -> set_point(i, j, frame_ -> get_point(i, j));
		}
	}
	return ;
}

void CTProjectorClass::importAngles(int nAngles_, double *angle_){
	if(angle != NULL) delete angle;
	if(projection != NULL){
		for(int i=0; i<nAngles; i++){
			delete projection[i];
		}
	}
	nAngles = nAngles_;
	angle = new double[nAngles];
	projection = new projectionClass*[nAngles];
	for(int i=0; i<nAngles; i++){
		projection[i] = new projectionClass(projectionSize);
	}
	for(int i=0; i<nAngles; i++){
		angle[i] = angle_[i];
	}
	return ;
}

void CTProjectorClass::project_sum(){
	for(int i=0; i<nAngles; i++){
		frameClass *rotatedFrame = frame -> rotate(angle[i]);
		projection[i] -> set_to(0.);
		for(int x=0; x<nx; x++){
			double sum = 0.;
			for(int y=0; y<ny; y++){
				sum += rotatedFrame -> get_point(x, y);
			}
			projection[i] -> set_point(x, sum);
		}
		delete rotatedFrame;
	}
	return ;
}

void CTProjectorClass::project_attenuate(double pixelSize, double Emin, double Emax, double dE){

	bool debug = false;

	double openBeamNormalization = spectrum -> integrate(Emin, Emax, dE);
	if(debug) cout << "****** openBeamNormalization = " << openBeamNormalization << endl;

	for(int i=0; i<nAngles; i++){
		if(debug) cout << "Angle = " << i << endl;
		frameClass *rotatedFrame = frame -> rotate(angle[i]);
		if(debug) cout << "Rotated frame" << endl;
		projection[i] -> set_to(0.);
		if(debug) cout << "Projection " << i << " set to zero" << endl;
		if(debug) cout << "Looping on x:" << endl;
		for(int x=0; x<nx; x++){
			if(debug) cout << "\tx = " << x << endl;
			int nZ = 100; // number of chemical elements
			int Zarray[nZ];
			if(debug) cout << "\tZarray set to zero" << endl;
			for(int z=0; z<nZ; z++){
				Zarray[z] = 0;
			}
			for(int y=0; y<ny; y++){
				int Z = rotatedFrame -> get_point(x, y);
				if(Z != 0) Zarray[Z]++;
			}

			spectrumClass *tmpSpectrum = spectrum -> clone();
			if(debug) cout << "\tContent of Zarray:" << endl;
			for(int z=0; z<nZ; z++){
				if(Zarray[z] != 0){
					if(debug) cout << "\t\t" << Zarray[z] << " pixels at Z = " << z << endl;
					tmpSpectrum -> attenuate(z, Zarray[z] * pixelSize);
				}
			}

			if(debug) cout << "\tIntegral from " << Emin << " to " << Emax << " = " << tmpSpectrum -> integrate(Emin, Emax, dE) << endl;
			if(debug) cout << "\tRenormalized to open beam = " << tmpSpectrum -> integrate(Emin, Emax, dE) / openBeamNormalization << endl;

			/////////////////////////////// CRUCIAL!!!
			// here include per pixel calibration: Emin becomes actual pixel threshold
			projection[i] -> set_point(x, tmpSpectrum -> integrate(Emin, Emax, dE) / openBeamNormalization);
			//////////////////////////////////////////

			delete tmpSpectrum;
		}
		delete rotatedFrame;
	}

	return ;
}

void CTProjectorClass::project_attenuate_withNoise(double /*pixelSize*/, double /*Emin*/, double /*Emax*/, double /*dE*/){

	// bool debug = false;

	// /////////////////////////////////////////////////////////////////////////////////////////////
	// // calculating open beam normalization factor = total number of photons in source spectrum //
	// /////////////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////// CRUCIAL!!!
	// // here include per pixel calibration: Emin becomes actual pixel threshold
	// double openBeamNormalization = spectrum -> integrate(Emin, Emax, dE);
	// //////////////////////////////////////////
	// if(debug) cout << "****** openBeamNormalization = " << openBeamNormalization << endl;

	// ////////////////////////////
	// // looping on view angles //
	// ////////////////////////////
	// for(int i=0; i<nAngles; i++){

	//   ////////////////////
	//   // rotating frame //
	//   ////////////////////
	//   if(debug) cout << "Angle = " << i << endl;
	//   frameClass *rotatedFrame = frame -> rotate(angle[i]);
	//   if(debug) cout << "Rotated frame" << endl;

	//   //////////////////////////
	//   // resetting projection //
	//   //////////////////////////
	//   projection[i] -> set_to(0.);
	//   if(debug) cout << "Projection " << i << " set to zero" << endl;
	//   if(debug) cout << "Looping on x:" << endl;

	//   //////////////////////////
	//   // looping on ray lines //
	//   //////////////////////////
	//   for(int x=0; x<nx; x++){

	//     /////////////////////////////////////////////////////
	//     // counting occurrences of Z values along ray line //
	//     /////////////////////////////////////////////////////
	//     if(debug) cout << "\tx = " << x << endl;
	//     int nZ = 100; // number of chemical elements
	//     int Zarray[nZ];
	//     if(debug) cout << "\tZarray set to zero" << endl;
	//     for(int z=0; z<nZ; z++){
	// 	Zarray[z] = 0;
	//     }
	//     for(int y=0; y<ny; y++){
	// 	int Z = rotatedFrame -> get_point(x, y);
	// 	if(Z != 0) Zarray[Z]++;
	//     }

	//     /////////////////////////////////////////////////////
	//     // attenuating spectrum according to Z occurrences //
	//     /////////////////////////////////////////////////////
	//     spectrumClass *tmpSpectrum = spectrum -> clone();
	//     if(debug) cout << "\tContent of Zarray:" << endl;
	//     for(int z=0; z<nZ; z++){
	// 	if(Zarray[z] != 0){
	// 	  if(debug) cout << "\t\t" << Zarray[z] << " pixels at Z = " << z << endl;
	// 	  tmpSpectrum -> attenuate(z, Zarray[z] * pixelSize);
	// 	}
	//     }

	//     if(debug) cout << "\tIntegral from " << Emin << " to " << Emax << " = " << tmpSpectrum -> integrate(Emin, Emax, dE) << endl;

	//     //////////////////////////////////////////////////////////////////////
	//     // calculating number of photons corresponding to detected spectrum //
	//     //////////////////////////////////////////////////////////////////////
	//     /////////////////////////////// CRUCIAL!!!
	//     // here include per pixel calibration: Emin becomes actual pixel threshold
	//     double nPhotons = tmpSpectrum -> integrate(Emin, Emax, dE);
	//     ///////////////////////////////////////////

	//     ///////////////////////////////
	//     // randomizing photon counts //
	//     ///////////////////////////////
	//     srand(x+i);
	//     int rand1 = rand() % 1000 + 1;
	//     srand(x-i);
	//     int rand2 = rand() % 1000 + 1;
	//     srand(rand1 + rand2);
	//     default_random_engine generator;
	//     poisson_distribution<int> distribution(nPhotons);
	//     int iSecret = rand() % 1000 + 1;
	//     if(debug) cout << "\tRandom number = " << iSecret << endl;
	//     int nPhotonsPoisson = distribution(generator);
	//     for(int iRandom=0; iRandom<iSecret; iRandom++){
	// 	nPhotonsPoisson = distribution(generator);
	//     }
	//     if(nPhotonsPoisson == 0) cout << " - WARNING!!! - CTProjectorClass::project_attenuate_withNoise(): 0 PHOTONS!!!" << endl;
	//     if(nPhotonsPoisson < 0) cout << " - WARNING!!! - CTProjectorClass::project_attenuate_withNoise(): NEGATIVE COUNT!!!" << endl;
	//     if(debug) cout << "\tWith Poisson noise: " << nPhotonsPoisson << endl;

	//     //////////////////////////////
	//     // setting projection value //
	//     //////////////////////////////
	//     projection[i] -> set_point(x, nPhotonsPoisson / openBeamNormalization);

	//     ////////////////////////////////////////
	//     // deleting temporary spectrum object //
	//     ////////////////////////////////////////
	//     delete tmpSpectrum;
	//   }

	//   /////////////////////////////////////
	//   // deleting temporary frame object //
	//   /////////////////////////////////////
	//   delete rotatedFrame;
	// }

	return ;
}

frameClass *CTProjectorClass::getFrame() const{
	return frame;
}

int CTProjectorClass::get_nAngles() const{
	return nAngles;
}

double *CTProjectorClass::get_angle() const{
	return angle;
}

spectrumClass *CTProjectorClass::get_spectrum() const{
	return spectrum;
}

double CTProjectorClass::get_projectionValue(int projectionIndex, int positionIndex) const{
	if(projectionIndex < 0 || projectionIndex > nAngles - 1){
		cout << " - ERROR!!! - CTProjectorClass::get_projectionValue(): index out of bound"<< endl;
		return default_double;
	}
	if(positionIndex < 0 || positionIndex > projectionSize - 1){
		cout << " - ERROR!!! - CTProjectorClass::get_projectionValue(): index out of bound"<< endl;
		return default_double;
	}
	return projection[projectionIndex] -> get_point(positionIndex);
}

void CTProjectorClass::set_projectionValue(int projectionIndex, int positionIndex, double value){
	if(projectionIndex < 0 || projectionIndex > nAngles - 1){
		cout << " - ERROR!!! - CTProjectorClass::get_projectionValue(): index out of bound"<< endl;
		return ;
	}
	if(positionIndex < 0 || positionIndex > projectionSize - 1){
		cout << " - ERROR!!! - CTProjectorClass::get_projectionValue(): index out of bound"<< endl;
		return ;
	}
	projection[projectionIndex] -> set_point(positionIndex, value);
	return ;
}

int CTProjectorClass::get_projectionSize() const{
	return projectionSize;
}

int CTProjectorClass::get_nx() const{
	return nx;
}

int CTProjectorClass::get_ny() const{
	return ny;
}

void CTProjectorClass::writeFrameToFile(const char *fileName){
	frame -> writeToFile_ASCIIMatrix(fileName);
	return ;
}

void CTProjectorClass::writeSinogram(const char *fileName){
	ofstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - CTProjectorClass::writeSinogram(): cannot open file " << fileName << endl;
		return ;
	}
	file << "#\t";
	for(int j=0; j<nAngles; j++){
		file << angle[j] << "\t";
	}
	file << endl;
	for(int i=0; i<projectionSize; i++){
		file << i << "\t";
		for(int j=0; j<nAngles; j++){
			file << projection[j] -> get_point(i) << "\t";
		}
		file << endl;
	}
	file.close();
	return ;
}

void CTProjectorClass::readSpectrumFromFile(const char *fileName){

	spectrum = new spectrumClass();
	spectrum -> readFromFile(fileName);

	return ;
}

void CTProjectorClass::importSinogram(const char *fileName, int projectionSize_){

	//////////////////////////////////////////
	// freeing angle and projector pointers //
	//////////////////////////////////////////
	projectionSize = projectionSize_;
	if(angle != NULL) delete angle;
	if(projection != NULL){
		for(int i=0; i<nAngles; i++){
			delete projection[i];
		}
	}

	////////////////////
	// getting angles //
	////////////////////
	ifstream file;
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - CTProjectorClass::importSinogram(): cannot open file " << fileName << endl;
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
	nAngles = angles.size() - 2;
	angle = new double[nAngles];
	projection = new projectionClass*[nAngles];
	for(int i=0; i<nAngles; i++){
		projection[i] = new projectionClass(projectionSize);
		angle[i] = angles[i+1];
	}

	/////////////////////////
	// getting projections //
	/////////////////////////
	file.open(fileName);
	if(file == NULL){
		cout << " - ERROR!!! - CTProjectorClass::importSinogram(): cannot open file " << fileName << endl;
		return ;
	}
	string dummy;
	for(int i=0; i<nAngles+1; i++){
		file >> dummy;
	}

	for(int i=0; i<projectionSize; i++){
		file >> dummy;
		for(int j=0; j<nAngles; j++){
			file >> dummy;
			double ang = atof(dummy.c_str());
			projection[j] -> set_point(i, ang);
		}
	}

	file.close();

	return ;
}

void CTProjectorClass::selectSinogramFromFrames(vector<frameClass *> frames, int nx, int ny, double angleMin, double angleMax, double angleStep, int selectedSlice, int sliceDirection, int lowerRange, int upperRange, double rescaleAngle){

	///////////////////////////////
	// assigning projection size //
	///////////////////////////////
	projectionSize = upperRange - lowerRange;
	if(projectionSize <= 0){
		cout << " - ERROR!!! - CTProjectorClass::selectSinogramFromFrames(): projection size (" << projectionSize << ") is null or negative: check lower and upper range" << endl;
		return ;
	}
	if(sliceDirection == 0){
		if(projectionSize > nx){
			cout << " - ERROR!!! - CTProjectorClass::selectSinogramFromFrames(): projection size (" << projectionSize << ") more than frame size (" << nx << ") along x" << endl;
			return ;
		}
	}
	else if(sliceDirection == 1){
		if(projectionSize > ny){
			cout << " - ERROR!!! - CTProjectorClass::selectSinogramFromFrames(): projection size (" << projectionSize << ") more than frame size (" << ny << ") along y" << endl;
			return ;
		}
	}
	else {
		cout << " - ERROR!!! - CTProjectorClass::selectSinogramFromFrames(): slice direction is undefined, must be either 0 (for x) or 1 (for y)" << endl;
		return ;
	}

	//////////////////////////////////////////
	// freeing angle and projector pointers //
	//////////////////////////////////////////
	if(angle != NULL) delete angle;
	if(projection != NULL){
		for(int i=0; i<nAngles; i++){
			delete projection[i];
		}
	}

	//////////////////////////////////////////////
	// assigning angles and getting projections //
	//////////////////////////////////////////////
	nAngles = 0;
	for(double ang = angleMin; ang <= angleMax; ang += angleStep){
		nAngles ++;
	}
	angle = new double[nAngles];
	projection = new projectionClass*[nAngles];
	for(int i=0; i<nAngles; i++){
		angle[i] = angleMin + i * angleStep;
		projection[i] = new projectionClass(projectionSize);
		projection[i] -> getSliceFromFrame(frames[i], nx, ny, selectedSlice, sliceDirection, lowerRange);
		angle[i] /= rescaleAngle;
	}

	return ;
}

int CTProjectorClass::fixRotationAxis(const char *chi2FileName, int rotationAxisRange){

	projectionClass *reverted180 = new projectionClass(projectionSize);
	for(int i=0; i<projectionSize; i++){
		reverted180 -> set_point(projectionSize - 1 - i, projection[nAngles - 1] -> get_point(i));
	}

	ofstream file;
	file.open(chi2FileName);
	if(file == NULL){
		cout << " - ERROR!!! - CTProjectorClass::fixRotationAxis(): cannot open file " << chi2FileName << endl;
		return default_int;
	}

	double chi2Min = 10000000;
	int offset = default_int;

	// shift left
	for(int i=0; i<projectionSize-1; i++){

		double chi2 = 0.;

		for(int j=0; j<i+1; j++){
			double diff = projection[0] -> get_point(j) - reverted180 -> get_point(projectionSize -1 - i + j);
			chi2 += (diff * diff);
		}

		chi2 /= ((double) (i+1));
		chi2 = sqrt(chi2);

		file << i - projectionSize + 1 << "\t" << chi2 << endl;
		if(fabs(i - projectionSize + 1) < rotationAxisRange && chi2 < chi2Min){
			chi2Min = chi2;
			offset = i - projectionSize + 1;
		}

	}

	// shift right
	for(int i=0; i<projectionSize; i++){

		double chi2 = 0.;

		for(int j=0; j<projectionSize - i; j++){
			double diff = projection[0] -> get_point(i + j) - reverted180 -> get_point(j);
			chi2 += (diff * diff);
		}

		chi2 /= ((double) (projectionSize - i));
		chi2 = sqrt(chi2);

		file << i << "\t" << chi2 << endl;
		if(i < rotationAxisRange && chi2 < chi2Min){
			chi2Min = chi2;
			offset = i;
		}

	}

	file.close();

	cout << " - CTProjectorClass::fixRotationAxis():         chi2 minimum = " << chi2Min << endl;
	cout << " - CTProjectorClass::fixRotationAxis(): rotation axis offset = " << offset << endl;

	////////////////////////
	// adjusting sinogram //
	////////////////////////
	shiftSinogram(offset);

	return offset;
}

int CTProjectorClass::fixRotationAxis_zeroPadding(const char *chi2FileName, int rotationAxisRange){

	projectionClass *reverted180 = new projectionClass(projectionSize);
	for(int i=0; i<projectionSize; i++){
		reverted180 -> set_point(projectionSize - 1 - i, projection[nAngles - 1] -> get_point(i));
	}

	ofstream file;
	file.open(chi2FileName);
	if(file == NULL){
		cout << " - ERROR!!! - CTProjectorClass::fitRotationAxis_zeroPadding(): cannot open file " << chi2FileName << endl;
		return default_int;
	}

	double chi2Min = 10000000;
	int offset = default_int;

	// shift left
	for(int i=0; i<projectionSize-1; i++){

		double chi2 = 0.;

		for(int j=0; j<i+1; j++){
			double diff = projection[0] -> get_point(j) - reverted180 -> get_point(projectionSize -1 - i + j);
			chi2 += (diff * diff);
		}

		chi2 /= ((double) (i+1));
		chi2 = sqrt(chi2);

		file << i - projectionSize + 1 << "\t" << chi2 << endl;
		if(fabs(i - projectionSize + 1) < rotationAxisRange && chi2 < chi2Min){
			chi2Min = chi2;
			offset = i - projectionSize + 1;
		}

	}

	// shift right
	for(int i=0; i<projectionSize; i++){

		double chi2 = 0.;

		for(int j=0; j<projectionSize - i; j++){
			double diff = projection[0] -> get_point(i + j) - reverted180 -> get_point(j);
			chi2 += (diff * diff);
		}

		chi2 /= ((double) (projectionSize - i));
		chi2 = sqrt(chi2);

		file << i << "\t" << chi2 << endl;
		if(i < rotationAxisRange && chi2 < chi2Min){
			chi2Min = chi2;
			offset = i;
		}

	}

	file.close();

	cout << " - CTProjectorClass::fitRotationAxis_zeroPadding():         chi2 minimum = " << chi2Min << endl;
	cout << " - CTProjectorClass::fitRotationAxis_zeroPadding(): rotation axis offset = " << offset << endl;


	////////////////////////
	// adjusting sinogram //
	////////////////////////
	shiftSinogram_zeroPadding(offset);

	return offset;
}

void CTProjectorClass::shiftSinogram(int shift){

	if(shift >= projectionSize){
		cout << " - ERROR!!! - CTProjectorClass::shiftSinogram(): shift size is more than projection size" << endl;
		return ;
	}
	int newSize = projectionSize - fabs(shift);

	projectionClass **newProjection = new projectionClass*[nAngles];
	for(int i=0; i<nAngles; i++){
		newProjection[i] = new projectionClass(newSize);
		for(int j=0; j<newSize; j++){
			if(shift > 0) newProjection[i] -> set_point(j, projection[i] -> get_point(shift + j));
			else newProjection[i] -> set_point(j, projection[i] -> get_point(j));
		}
		projectionClass *tmp = projection[i];
		projection[i] = newProjection[i];
		delete tmp;
	}
	projectionSize = newSize;

	return ;
}

void CTProjectorClass::shiftSinogram_zeroPadding(int shift){

	if(shift >= projectionSize){
		cout << " - ERROR!!! - CTProjectorClass::shiftSinogram_zeroPadding(): shift size is more than projection size" << endl;
		return ;
	}
	int newSize = projectionSize + fabs(shift);

	projectionClass **newProjection = new projectionClass*[nAngles];
	for(int i=0; i<nAngles; i++){
		newProjection[i] = new projectionClass(newSize);
		if(shift<0){
			for(int j=0; j<fabs(shift); j++){
				newProjection[i] -> set_point(j, 0.);
			}
			for(int j=fabs(shift); j<newSize; j++){
				newProjection[i] -> set_point(j, projection[i] -> get_point(j-fabs(shift)));
			}
		}
		else{
			for(int j=0; j<projectionSize; j++){
				newProjection[i] -> set_point(j, projection[i] -> get_point(j));
			}
			for(int j=projectionSize; j<newSize; j++){
				newProjection[i] -> set_point(j, 0.);
			}
		}

		projectionClass *tmp = projection[i];
		projection[i] = newProjection[i];
		delete tmp;
	}
	projectionSize = newSize;

	return ;
}

void CTProjectorClass::signalToThicknessCalibration(const char *fileName_STT){

	//////////////////////////////////////////////////
	// getting signal-to-thickness calibration data //
	//////////////////////////////////////////////////
	vector<double> thickness;
	vector<double> photonCounts;
	ifstream file;
	file.open(fileName_STT);
	if(file == NULL){
		cout << " - ERROR!!! - CTProjectorClass::signalToThicknessCalibration(): cannot open file " << fileName_STT << endl;
		return ;
	}
	double val;
	while(file >> val){
		thickness.push_back(val);
		file >> val;
		photonCounts.push_back(val);
	}
	file.close();

	//////////////////////////
	// calibrating sinogram //
	//////////////////////////
	for(int i=0; i<nAngles; i++){
		for(int j=0; j<projectionSize; j++){

			double nPhotons = projection[i] -> get_point(j);
			if(nPhotons > photonCounts[0]){
				cout << " - WARNING!!! - CTProjectorClass::signalToThicknessCalibration(): strange situation!!!: nPhotons = " << nPhotons << ", photonCounts = " << photonCounts[0] << endl;
				projection[i] -> set_point(j, thickness[0]);
			}
			else if(nPhotons < photonCounts[photonCounts.size()-1]){
				cout << " - WARNING!!! - CTProjectorClass::signalToThicknessCalibration(): STT calibration has to be extended to further thicknesses" << endl;
				projection[i] -> set_point(j, thickness[photonCounts.size()-1]);
			}
			else{
				for(int k=0; k<(int)photonCounts.size(); k++){
					if(nPhotons < photonCounts[k] && nPhotons >= photonCounts[k+1]){
						projection[i] -> set_point(j, thickness[k]);
						break;
					}
				}
			}

		}
	}

	return ;
}

void CTProjectorClass::BP(){

	frame -> set_to(0.);

	for(int k=0; k<nAngles; k++) {

		frameClass *frameTmp = new frameClass(nx, ny);
		frameTmp -> set_to(0.);
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				frameTmp -> set_point(i, j, projection[k] -> get_point(i));
			}
		}
		frameClass *frameTmpRotated = frameTmp -> rotate(-angle[k]);
		frameClass *frameTmpRotatedCircle = frameTmpRotated -> makeCircle_in(0.);
		frame -> addFrame(frameTmpRotated);

		delete frameTmp;
		delete frameTmpRotated;
		delete frameTmpRotatedCircle;
	}
	frame -> multiplyBy(1./((double)get_nAngles() * ny));
	frame -> setToZeroBelow(0.); // added later

	return ;
}

void CTProjectorClass::FBP(){

	//////////////////////////////////////
	// setting to zero negative numbers //
	//////////////////////////////////////
	cout << " - CTProjectorClass::FBP(): setting to zero negative numbers" << endl;
	eliminateNegatives();

	//////////////////////
	// FT and filtering //
	//////////////////////
	FourierClass_1D **fft = new FourierClass_1D*[nAngles];
	for(int i=0; i<nAngles; i++){
		//////////////////////
		// getting sinogram //
		//////////////////////
		fft[i] = new FourierClass_1D(projectionSize);
		fft[i] -> getRealData(projectionSize, projection[i] -> get_projection());
		//////////////////////////////////////
		// Fourier transform of projections //
		//////////////////////////////////////
		fft[i] -> fft();
		///////////////
		// filtering //
		///////////////
		fft[i] -> filterDataOut();
		fft[i] -> rfft();
	}

	////////////////////
	// backprojection //
	////////////////////
	double origin = default_double;
	int size = projectionSize;
	if(size%2 == 0) origin = size/2 - 1;
	else origin = size/2;
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			double df = 0.;
			for(int k=0; k<nAngles; k++){
				double cosTheta = cos(M_PI * angle[k] / 180.);
				double sinTheta = sin(M_PI * angle[k] / 180.);
				double t = ((i-origin) * cosTheta - (j-origin) * sinTheta) + projectionSize / 2.;
				int point = t;
				double decimal = t - point;
				if(t < 0 || t > projectionSize - 1) continue;
				else if(t == 0){
					df += fft[k] -> dataBack[0][0];
				}
				else if(t == projectionSize - 1){
					df += fft[k] -> dataBack[projectionSize - 1][0];
				}
				else{
					df += ((fft[k] -> dataBack[point][0] + decimal * (fft[k] -> dataBack[point + 1][0] - fft[k] -> dataBack[point][0])));
				}
			}
			//      frame -> set_point(i, j, M_PI * df / (nAngles * projectionSize)); // my guess -> reproduces good CT numbers, does not depend on nAngles, depends weakly on projectionSize
			frame -> set_point(i, j, M_PI * df / nAngles); // from the book
		}
	}

	for(int i=0; i<nAngles; i++){
		delete fft[i];
	}

	frameClass *tmpFrame = frame -> makeCircle_in(0.);
	frame -> copyFromFrame(tmpFrame);
	delete tmpFrame;

	return ;
}

void CTProjectorClass::eliminateNegatives(){

	for(int i=0; i<nAngles; i++){
		for(int j=0; j<projectionSize; j++){
			if(projection[i] -> get_point(j) < 0.) projection[i] -> set_point(j, 0.);
		}
	}

	return ;
}

void CTProjectorClass::OSEM(int subsetSize, const char *chi2File_name, double chi2Cut_subsets, double chi2Cut_fullSet){

	//////////////
	// controls //
	//////////////
	if(subsetSize > nAngles){
		cout << " - ERROR!!! - CTProjectorClass::OSEM() subset size (" << subsetSize << ") is greater than number of projections (" << nAngles << ")" << endl;
		return ;
	}
	if(subsetSize <=0){
		cout << " - ERROR!!! - CTProjectorClass::OSEM() subset size (" << subsetSize << ") is not acceptable" << endl;
		return ;
	}

	/////////////
	// options //
	/////////////
	bool debug = false;
	bool enableTemporaryOutput = false;

	//////////////////////////////////////
	// setting to zero negative numbers //
	//////////////////////////////////////
	cout << " - CTProjectorClass::OSEM(): setting to zero negative numbers" << endl;
	eliminateNegatives();

	///////////////////////
	// opening chi2 file //
	///////////////////////
	ofstream chi2File;
	chi2File.open(chi2File_name);
	if(chi2File == NULL){
		cout << " - ERROR!!! - CTProjectorClass::OSEM(): cannot open file " << chi2File_name<< endl;
		return ;
	}

	//////////////////////
	// creating subsets //
	//////////////////////
	if(debug) cout << " ******* Creating subsets" << endl;
	int nSubsets = nAngles / subsetSize;
	if(debug) cout << " - nSubsets = " << nSubsets<< endl;
	int nExtraProjections = nAngles - nSubsets * subsetSize;
	if(debug) cout << " - nExtraProjections = " << nExtraProjections << endl;
	vector<int> subset[nSubsets];
	for(int i=0; i<nSubsets; i++){
		for(int j=0; j<subsetSize; j++){
			subset[i].push_back(i + nSubsets * j);
		}
	}
	for(int i=0; i<nExtraProjections; i++){
		subset[nSubsets - 1].push_back(nAngles - 1 - i);
	}
	if(debug){
		cout << " - Subsets:" << endl;
		for(int i=0; i<nSubsets; i++){
			cout << "\tsubset[" << i << "] = { ";
			for(int j=0; j<(int)subset[i].size(); j++){
				cout << subset[i][j];
				if(j != (int)subset[i].size() - 1) cout << " , ";
			}
			cout << " }" << endl;
		}
	}

	///////////////////////
	// iterating subsets //
	///////////////////////
	if(debug) cout << " ******* Iterating subsets" << endl;
	cout << "\tInitializing frame" << endl;
	OSEM_initializeFrame();
	if(enableTemporaryOutput) frame -> writeToFile_ASCIIMatrix("idealDetector_data/tmp/OSEM_initialization.txt");
	if(debug) cout << " - Starting iteration" << endl;

	int countSubIteration = 0;
	double chi2Old = 0.;
	double chi2Diff = 1000.;
	int iteration = 0;
	bool stopIterationSubsets = false;

	while(stopIterationSubsets == false) {

		if(debug) cout << " - Iteration " << iteration << endl;
		double chi2Average = 0.;

		for(int i=0; i<nSubsets; i++) {
			if(debug) cout << "\tsubset " << i << endl;
			if(debug) cout << "\tsubiteration " << countSubIteration << endl;
			countSubIteration ++;

			////////////////////////////////
			// creating OSEMsubset object //
			////////////////////////////////
			if(debug) cout << "\t\tcreating OSEMsubset object" << endl;
			OSEMsubsetClass *OSEMsubset = new OSEMsubsetClass(subset[i].size(), nx);
			for(int j=0; j<(int)subset[i].size(); j++){
				OSEMsubset -> assignAngleIndex(j, subset[i][j]);
			}

			///////////////////////
			// projecting subset //
			///////////////////////
			if(debug) cout << "\t\tprojecting frame onto subset" << endl;
			OSEM_projectSubset(OSEMsubset);

			if(enableTemporaryOutput){
				for(int j=0; j<(int)subset[i].size(); j++){
					stringstream ii_iteration;
					ii_iteration << iteration;
					stringstream ii_subset;
					ii_subset << i;
					stringstream ii_angle;
					ii_angle << subset[i][j];
					string tmp_fileName_projection = "idealDetector_data/tmp/OSEM_projection_iteration_" + ii_iteration.str() + "_subset_" + ii_subset.str() + "_angle_" + ii_angle.str() +".txt";
					OSEMsubset -> OSEMprojection[j] -> writeToFile(tmp_fileName_projection.c_str());
				}
			}

			/////////////////////////////
			// calculating corrections //
			/////////////////////////////
			if(debug) cout << "\t\tcalculating corrections" << endl;
			OSEM_calculateCorrections(OSEMsubset);

			if(enableTemporaryOutput){
				for(int j=0; j<(int)subset[i].size(); j++){
					stringstream ii_iteration;
					ii_iteration << iteration;
					stringstream ii_subset;
					ii_subset << i;
					stringstream ii_angle;
					ii_angle << subset[i][j];
					string tmp_fileName_correction = "idealDetector_data/tmp/OSEM_correction_iteration_" + ii_iteration.str() + "_subset_" + ii_subset.str() + "_angle_" + ii_angle.str() +".txt";
					OSEMsubset -> OSEMprojection[j] -> writeToFile(tmp_fileName_correction.c_str());
				}
			}

			//////////////////////
			// calculating chi2 //
			//////////////////////
			if(debug) cout << "\t\tcalculating chi2" << endl;
			double chi2 = OSEM_calculateChi2(OSEMsubset);
			if(chi2 < 0){
				cout << " - ERROR!!! - CTProjectorClass::OSEM(): irregular chi2" << endl;
				return ;
			}
			cout << "\titeration " << iteration << " subset " << i << "/" << nSubsets << endl;
			chi2Average += chi2;

			////////////////////////////////
			// backprojecting corrections //
			////////////////////////////////
			if(debug) cout << "\t\tbackprojecting corrections" << endl;
			frameClass *OSEM_BPFrame = OSEM_backprojectSubset(OSEMsubset);

			if(enableTemporaryOutput){
				stringstream ii_iteration;
				ii_iteration << iteration;
				stringstream ii_subset;
				ii_subset << i;
				string tmp_fileName = "idealDetector_data/tmp/OSEM_frame_correction_iteration_" + ii_iteration.str() + "_subset_" + ii_subset.str() + ".txt";
				OSEM_BPFrame -> writeToFile_ASCIIMatrix(tmp_fileName.c_str());
			}

			//////////////////////////
			// applying corrections //
			//////////////////////////
			if(debug) cout << "\t\tapplying corrections" << endl;
			OSEM_applyCorrections(OSEM_BPFrame);

			if(enableTemporaryOutput){
				stringstream ii_iteration;
				ii_iteration << iteration;
				stringstream ii_subset;
				ii_subset << i;
				string tmpFrame_fileName = "idealDetector_data/tmp/OSEM_frame_iteration_" + ii_iteration.str() + "_subset_" + ii_subset.str() + ".txt";
				frame -> writeToFile_ASCIIMatrix(tmpFrame_fileName.c_str());
			}

			////////////////////////////////
			// deleting temporary objects //
			////////////////////////////////
			delete OSEM_BPFrame;
			delete OSEMsubset;
		}
		chi2Average /= nSubsets;
		cout << "\tAverage chi2 data at iteration " << iteration << ":" << endl;
		cout << "\t\told chi2 = " << chi2Old << endl;
		cout << "\t\tnew chi2 = " << chi2Average << endl;
		chi2Diff = 2 * fabs(chi2Old - chi2Average) / (chi2Old + chi2Average);
		cout << "\t\t  % diff = " << chi2Diff << endl;
		cout << "\t\t   limit = " << chi2Cut_subsets << endl;
		chi2File << iteration << "\t" << chi2Average << endl;
		if(iteration > 2 && chi2Average > chi2Old && stopIterationSubsets == false){
			cout << "\t\tCHI2 DIVERGING" << endl;
			stopIterationSubsets = true;
		}
		if(chi2Diff < chi2Cut_subsets && stopIterationSubsets == false){
			cout << "\t\tCHI2 LIMIT REACHED" << endl;
			stopIterationSubsets = true;
		}
		if(stopIterationSubsets == true) cout << "\t\tLAST ITERATION" << endl;
		chi2Old = chi2Average;
		iteration ++;
	}

	////////////////////////
	// iterating full set //
	////////////////////////
	if(debug) cout << "******* iterating FULL SET" << endl;
	stopIterationSubsets = false;

	while(stopIterationSubsets == false){

		countSubIteration ++;

		/////////////////////////////////////////////
		// creating OSEMsubset object for full set //
		/////////////////////////////////////////////
		if(debug) cout << "\t\tcreating OSEMsubset object for full set" << endl;
		OSEMsubsetClass *OSEMsubset = new OSEMsubsetClass(nAngles, nx);
		for(int j=0; j<nAngles; j++){
			OSEMsubset -> assignAngleIndex(j, j);
		}

		//////////////////////////
		// projecting  full set //
		//////////////////////////
		if(debug) cout << "\t\tprojecting frame onto full set" << endl;
		OSEM_projectSubset(OSEMsubset);

		if(enableTemporaryOutput){
			for(int j=0; j<nAngles; j++){
				stringstream ii_iteration;
				ii_iteration << iteration;
				stringstream ii_angle;
				ii_angle << j;
				string tmp_fileName_projection = "idealDetector_data/tmp/OSEM_fullSet_projection_iteration_" + ii_iteration.str() + "_angle_" + ii_angle.str() +".txt";
				OSEMsubset -> OSEMprojection[j] -> writeToFile(tmp_fileName_projection.c_str());
			}
		}

		/////////////////////////////
		// calculating corrections //
		/////////////////////////////
		if(debug) cout << "\t\tcalculating corrections" << endl;
		OSEM_calculateCorrections(OSEMsubset);

		//////////////////////
		// calculating chi2 //
		//////////////////////
		if(debug) cout << "\t\tcalculating chi2" << endl;
		double chi2 = OSEM_calculateChi2(OSEMsubset);
		if(chi2 < 0){
			cout << " - ERROR!!! - CTProjectorClass::OSEM(): irregular chi2" << endl;
			return ;
		}
		chi2File << iteration << "\t" << chi2 << endl;
		cout << "\tfull set iteration " << iteration << endl;
		cout << "\t\told chi2 = " << chi2Old << endl;
		cout << "\t\tnew chi2 = " << chi2 << endl;
		chi2Diff = 2 * fabs(chi2Old - chi2) / (chi2Old + chi2);
		cout << "\t\t  % diff = " << chi2Diff << endl;
		cout << "\t\t   limit = " << chi2Cut_fullSet << endl;
		if(chi2Old - chi2 < 0. && stopIterationSubsets == false){
			cout << "\t\tCHI2 DIVERGING" << endl;
			stopIterationSubsets = true;
		}
		if(chi2Diff < chi2Cut_fullSet && stopIterationSubsets == false){
			cout << "\t\tCHI2 LIMIT REACHED" << endl;
			stopIterationSubsets = true;
		}
		if(stopIterationSubsets == true) cout << "\t\tLAST ITERATION" << endl;
		chi2Old = chi2;

		////////////////////////////////
		// backprojecting corrections //
		////////////////////////////////
		if(debug) cout << "\t\tbackprojecting corrections" << endl;
		frameClass *OSEM_BPFrame = OSEM_backprojectSubset(OSEMsubset);

		//////////////////////////
		// applying corrections //
		//////////////////////////
		if(debug) cout << "\t\tapplying corrections" << endl;
		OSEM_applyCorrections(OSEM_BPFrame);

		if(enableTemporaryOutput){
			stringstream ii_iteration;
			ii_iteration << iteration;
			string tmpFrame_fileName = "idealDetector_data/tmp/OSEM_frame_fullSet_iteration_" + ii_iteration.str() + ".txt";
			frame -> writeToFile_ASCIIMatrix(tmpFrame_fileName.c_str());
		}

		////////////////////////////////
		// deleting temporary objects //
		////////////////////////////////
		delete OSEM_BPFrame;
		delete OSEMsubset;

		iteration ++;
	}
	chi2File.close();

	return ;
}

void CTProjectorClass::OSEM_initializeFrame(){

	BP();

	return ;
}

void CTProjectorClass::OSEM_projectAtAngle(double angleValue, projectionClass *projection){

	frameClass *rotatedFrame = frame -> rotate(angleValue);
	projection -> set_to(0.);
	for(int x=0; x<nx; x++){
		double sum = 0.;
		for(int y=0; y<ny; y++){
			double point = rotatedFrame -> get_point(x, y);
			if(std::isinf(point) || std::isnan(point)){
				cout << " - ERROR!!! - CTProjectorClass::OSEM_projectAtAngle(): non regular data point" << endl;
				return ;
			}
			sum += point;
		}
		projection -> set_point(x, sum);
	}
	delete rotatedFrame;

	return ;
}

void CTProjectorClass::OSEM_attenuateAtAngle(double angleValue, projectionClass *projection, double pixelSize, double Emin, double Emax, double dE){

	///////////////////////////////////////
	// computing open beam normalization //
	///////////////////////////////////////
	double openBeamNormalization = spectrum -> integrate(Emin, Emax, dE);

	////////////////////
	// rotating frame //
	////////////////////
	frameClass *rotatedFrame = frame -> rotate(angleValue);

	/////////////////////////////
	// initializing projection //
	/////////////////////////////
	projection -> set_to(0.);

	////////////////////////////////
	// starting loop on ray lines //
	////////////////////////////////
	for(int x=0; x<nx; x++){

		//////////////////////////////////////
		// counting occurrences of Z values //
		//////////////////////////////////////
		int nZ = 100; // number of chemical elements
		int Zarray[nZ];
		for(int z=0; z<nZ; z++){
			Zarray[z] = 0;
		}
		for(int y=0; y<ny; y++){
			int Z = rotatedFrame -> get_point(x, y);
			if(Z != 0) Zarray[Z]++;
		}

		/////////////////
		// attenuating //
		/////////////////
		spectrumClass *tmpSpectrum = spectrum -> clone();
		for(int z=0; z<nZ; z++){
			if(Zarray[z] != 0){
				tmpSpectrum -> attenuate(z, Zarray[z] * pixelSize);
			}
		}

		////////////////////////////////////////////////////////////////
		// calculating number of photons and normalizing to open beam //
		////////////////////////////////////////////////////////////////
		double value = -log(tmpSpectrum -> integrate(Emin, Emax, dE) / openBeamNormalization);
		projection -> set_point(x, value);

		/////////////////////////////////
		// deleting temporary spectrum //
		/////////////////////////////////
		delete tmpSpectrum;
	}

	//////////////////////////////
	// deleting temporary frame //
	//////////////////////////////
	delete rotatedFrame;

	return ;
}

void CTProjectorClass::OSEM_projectSubset(OSEMsubsetClass *OSEMsubset){

	for(int i=0; i<OSEMsubset -> nProjections; i++){
		OSEM_projectAtAngle(angle[OSEMsubset -> angleIndex[i]], OSEMsubset -> OSEMprojection[i]);
	}

	return ;
}

void CTProjectorClass::OSEM_attenuateSubset(OSEMsubsetClass *OSEMsubset, double pixelSize, double Emin, double Emax, double dE){

	for(int i=0; i<OSEMsubset -> nProjections; i++){
		OSEM_attenuateAtAngle(angle[OSEMsubset -> angleIndex[i]], OSEMsubset -> OSEMprojection[i], pixelSize, Emin, Emax, dE);
		// temporary output
		stringstream ii_angle;
		ii_angle << OSEMsubset -> angleIndex[i];
		string fileName = "idealDetector_data/tmp/specOSEM_projection_angle_" + ii_angle.str() +".txt";
		OSEMsubset -> OSEMprojection[i] -> writeToFile(fileName.c_str());
		// end temporary output
	}

	return ;
}

void CTProjectorClass::OSEM_calculateCorrections(OSEMsubsetClass *OSEMsubset){

	for(int i=0; i<OSEMsubset -> nProjections; i++){
		OSEM_calculateCorrection(OSEMsubset -> OSEMprojection[i], OSEMsubset -> angleIndex[i]);
		// // temporary output
		// stringstream ii_angle;
		// ii_angle << OSEMsubset -> angleIndex[i];
		// string fileName = "idealDetector_data/tmp/OSEM_correction_angle_" + ii_angle.str() +".txt";
		// OSEMsubset -> OSEMprojection[i] -> writeToFile(fileName.c_str());
		// // end temporary output
	}

	return ;
}

void CTProjectorClass::OSEM_calculateCorrection(projectionClass *OSEM_projection, int projectionIndex){

	for(int i=0; i<OSEM_projection -> get_size(); i++){
		double ratio = projection[projectionIndex] -> get_point(i) / OSEM_projection -> get_point(i); // expectation / data (???)
		if(std::isinf(ratio) || std::isnan(ratio) || ratio > 10){
			//      cout << " - WARNING!!! - CTProjectorClass::OSEM_calculateCorrection(): correction is not a number: temporary solution by setting ratio = 0" << endl;
			ratio = 0.; // temporary solution
		}
		OSEM_projection -> set_point(i, ratio);
	}

	return ;
}

double CTProjectorClass::OSEM_calculateChi2(OSEMsubsetClass *OSEMsubset){

	double chi2_mean = 0.;

	for(int i=0; i<OSEMsubset -> nProjections; i++){
		double chi2= 0.;
		for(int j=0; j<OSEMsubset -> OSEMprojection[i] -> get_size(); j++){
			double dist = OSEMsubset -> OSEMprojection[i] -> get_point(j) - 1.;
			if(std::isinf(dist) || std::isnan(dist)){
				cout << " - ERROR!!! - CTProjectorClass::OSEM_calculateChi2(): distance from 1. is not a regular number: " << dist << endl;
				return -777.888;
			}
			chi2 += (dist * dist);
		}
		chi2 /= OSEMsubset -> OSEMprojection[i] -> get_size();
		chi2_mean += chi2;
	}

	chi2_mean /= OSEMsubset -> nProjections;
	if(std::isinf(chi2_mean) || std::isnan(chi2_mean)){
		cout << " - ERROR!!! - CTProjectorClass::OSEM_calculateChi2(): chi2 is not a regular number" << endl;
		return -777.888;
	}

	return chi2_mean;
}

frameClass *CTProjectorClass::OSEM_backprojectSubset(OSEMsubsetClass *OSEMsubset){

	frameClass *OSEM_BPFrame = new frameClass(nx, ny);

	OSEM_BPFrame -> set_to(0.);
	for(int k=0; k<OSEMsubset -> nProjections; k++){

		frameClass *frameTmp = new frameClass(nx, ny);
		frameTmp -> set_to(0.);
		for(int i=0; i<nx; i++){
			for(int j=0; j<ny; j++){
				frameTmp -> set_point(i, j, OSEMsubset -> OSEMprojection[k] -> get_point(i));
			}
		}

		frameClass *frameTmpRotated = frameTmp -> rotate(-angle[OSEMsubset -> angleIndex[k]]);
		frameClass *frameTmpRotatedCircle = frameTmpRotated -> makeCircle_in(0.);
		OSEM_BPFrame -> addFrame(frameTmpRotated);

		delete frameTmp;
		delete frameTmpRotated;
		delete frameTmpRotatedCircle;
	}
	OSEM_BPFrame -> multiplyBy(1./OSEMsubset -> nProjections);
	OSEM_BPFrame -> setToZeroBelow(0.);

	return OSEM_BPFrame;
}

void CTProjectorClass::OSEM_applyCorrections(frameClass *OSEM_BPFrame){

	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			double newValue = frame -> get_point(i, j) * OSEM_BPFrame -> get_point(i, j);
			frame -> set_point(i, j, newValue);
		}
	}

	return ;
}

void CTProjectorClass::renormalizeFrame(){
	frame -> normalizeToMaximum();
	return ;
}

void CTProjectorClass::logarithmOfSinogram(){
	for(int i=0; i<nAngles; i++){
		projection[i] -> minusLogarithm();
	}
	return ;
}

void CTProjectorClass::computeAverageFromCircle(double radius, double &mean, double &RMS){
	frame -> computeAverageFromCircle(radius, mean, RMS);
	return ;
}

double CTProjectorClass::get_minimumCTnumber(){
	return frame -> get_minimum();
}

double CTProjectorClass::get_maximumCTnumber(){
	return frame -> get_maximum();
}

void CTProjectorClass::copyProjectionsFromOtherProjector(CTProjectorClass *otherProjector){

	/////////////////////////
	// check compatibility //
	/////////////////////////
	if(otherProjector -> get_nAngles() != get_nAngles()){
		cout << " - ERROR!!! - CTProjectorClass::copyProjectionsFromOtherProjector(): incompatible number of angles" << endl;
		return ;
	}
	if(otherProjector -> get_projectionSize() != get_projectionSize()){
		cout << " - ERROR!!! - CTProjectorClass::copyProjectionsFromOtherProjector(): incompatible projection size" << endl;
		return ;
	}

	/////////////
	// copying //
	/////////////
	for(int i=0; i<nAngles; i++){
		for(int j=0; j<projectionSize; j++){
			projection[i] -> set_point(j, otherProjector -> projection[i] -> get_point(j));
		}
	}

	return ;
}

double CTProjectorClass::getAverageOfSinogram(){

	double average = 0.;

	for(int i=0; i<nAngles; i++){
		for(int j=0; j<projectionSize; j++){
			average += projection[i] -> get_point(j);
		}
	}

	return average / (nAngles * projectionSize);
}

double CTProjectorClass::getAverageOfSinogramInRange(double min, double max){

	double average = 0.;
	double count = 0.;

	for(int i=0; i<nAngles; i++){
		for(int j=0; j<projectionSize; j++){
			double val = projection[i] -> get_point(j);
			if(val>=min && val<=max){
				average += val;
				count ++;
			}
		}
	}

	return average / (count);
}

void CTProjectorClass::correctQuadCrossArtifact(int position){

	projectionSize += 4;

	//////////////////////////////
	// creating new projections //
	//////////////////////////////
	projectionClass **newProjection = new projectionClass*[nAngles];
	for(int i=0; i<nAngles; i++){
		newProjection[i] = new projectionClass(projectionSize);
	}

	/////////////////////////////
	// filling new projections //
	/////////////////////////////
	for(int i=0; i<projectionSize; i++){
		for(int j=0; j<nAngles; j++){
			if(i < position) newProjection[j] -> set_point(i, projection[j] -> get_point(i));
			else if(i == position || i == position+1 || i == position+2) newProjection[j] -> set_point(i, projection[j] -> get_point(position));
			else if(i == position+3 || i == position+4 || i == position+5) newProjection[j] -> set_point(i, projection[j] -> get_point(position+1));
			else newProjection[j] -> set_point(i, projection[j] -> get_point(i-4));
		}
	}

	//////////////////////////////////////////
	// freeing angle and projector pointers //
	//////////////////////////////////////////
	for(int i=0; i<nAngles; i++){
		delete projection[i];
	}
	delete projection;

	projection = newProjection;

	return ;
}
