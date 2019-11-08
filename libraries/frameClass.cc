#include "frameClass.hh"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

frameClass::frameClass(){
	default_int = -777;
	default_double = -777.;
	_nx = default_int;
	_ny = default_int;
	_frame = NULL;
	_frame2D = NULL;
	initialized = false;
}

frameClass::frameClass(int nx_, int ny_){
	default_int = -777;
	default_double = -777.;
	_nx = nx_;
	_ny = ny_;
	initialize();
}

frameClass::~frameClass(){

	delete [] _frame;

	for(int i=0; i<_nx; i++){
		delete [] _frame2D[i];
	}
	delete [] _frame2D;
}

void frameClass::set_nx(int nx_){
	_nx = nx_;
	return ;
}

void frameClass::set_ny(int ny_){
	_ny = ny_;
	return ;
}

int frameClass::get_nx() const{
	return _nx;
}

int frameClass::get_ny() const{
	return _ny;
}

void frameClass::set_to(double val){
	if(!initialized) cout << " - ERROR!!! - frameClass::set_to(): frame is not yet allocated" << endl;
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			_frame[(j*_nx)+i] = val;
		}
	}
	return ;
}

void frameClass::set_point(int x, int y, double val){
	//if(!initialized) cout << " - ERROR!!! - frameClass::set_point(): frame is not yet allocated" << endl;
	//if(x >= _nx) cout << " - ERROR!!! - frameClass::set_point(): x index (" << x << ") out of bound (" << _nx << ")" << endl;
	//if(y >= _ny) cout << " - ERROR!!! - frameClass::set_point(): y index (" << y << ") out of bound (" << _ny << ")" << endl;
	_frame[(y*_nx)+x] = val;
	return ;
}

double frameClass::get_point(int x, int y) const{
	//if(!initialized) cout << " - ERROR!!! - frameClass::get_point(): frame is not yet allocated" << endl;
	//if(x >= _nx) cout << " - ERROR!!! - frameClass::get_point(): x index (" << x << ") out of bound (" << _nx << ")" << endl;
	//if(y >= _ny) cout << " - ERROR!!! - frameClass::get_point(): y index (" << y << ") out of bound (" << _ny << ")" << endl;
	return _frame[(y*_nx)+x];
	//return _frame[1];
}

/*
double ** frameClass::get_frame(){

	// Only if requested copy the contents in the 2D array and return it
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			_frame2D[i][j] = _frame[(i*_nx)+j];
		}
	}
	return _frame2D;
}
*/

double * frameClass::get_frame(){
	return _frame;
}

bool frameClass::is_initialized(){
	if(initialized) return true;
	else return false;
}

void frameClass::setToZeroBelow(double threshold){
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			if(_frame[(j*_nx)+i] < threshold) _frame[(j*_nx)+i] = 0.;
		}
	}
	return ;
}

void frameClass::binarize(double threshold){
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			if(_frame[(j*_nx)+i] < threshold) _frame[(j*_nx)+i] = 0.;
			else _frame[(j*_nx)+i] = 1.;
		}
	}
	return ;
}

void frameClass::subtractFrames(frameClass *frame1, frameClass *frame2){

	if(_nx != frame1 -> _nx || _ny != frame1 -> _ny){
		cout << " - ERROR!!! - frameClass::subtractFrames(): incompatible size with first frame" << endl;
		return ;
	}
	if(_nx != frame2 -> _nx || _ny != frame2 -> _ny){
		cout << " - ERROR!!! - frameClass::subtractFrames(): incompatible size with second frame" << endl;
		return ;
	}

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double diff = frame1 -> get_point(i, j) - frame2 -> get_point(i, j);
			set_point(i, j, diff);
		}
	}

	return ;
}

void frameClass::sumFrames(frameClass *frame1, frameClass *frame2){

	if(_nx != frame1 -> _nx || _ny != frame1 -> _ny){
		cout << " - ERROR!!! - frameClass::sumFrames(): incompatible size with first frame" << endl;
		return ;
	}
	if(_nx != frame2 -> _nx || _ny != frame2 -> _ny){
		cout << " - ERROR!!! - frameClass::sumFrames(): incompatible size with second frame" << endl;
		return ;
	}

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double sum = frame1 -> get_point(i, j) + frame2 -> get_point(i, j);
			set_point(i, j, sum);
		}
	}

	return ;
}

void frameClass::divideFrames(frameClass *frame1, frameClass *frame2){

	if(_nx != frame1 -> _nx || _ny != frame1 -> _ny){
		cout << " - ERROR!!! - frameClass::divideFrames(): incompatible size with first frame" << endl;
		return ;
	}
	if(_nx != frame2 -> _nx || _ny != frame2 -> _ny){
		cout << " - ERROR!!! - frameClass::divideFrames(): incompatible size with second frame" << endl;
		return ;
	}

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double div = frame1 -> get_point(i, j) / frame2 -> get_point(i, j);
			if(frame2 -> get_point(i, j) == 0.) div = 0.;
			set_point(i, j, div);
		}
	}

	return ;
}

void frameClass::multiplyFrames(frameClass *frame1, frameClass *frame2){

	if(_nx != frame1 -> _nx || _ny != frame1 -> _ny){
		cout << " - ERROR!!! - frameClass::multiplyFrames(): incompatible size with first frame" << endl;
		return ;
	}
	if(_nx != frame2 -> _nx || _ny != frame2 -> _ny){
		cout << " - ERROR!!! - frameClass::multiplyFrames(): incompatible size with second frame" << endl;
		return ;
	}

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double product = frame1 -> get_point(i, j) * frame2 -> get_point(i, j);
			set_point(i, j, product);
		}
	}

	return ;
}

void frameClass::multiplyByFrame(frameClass *frame_){

	if(get_nx() != frame_ -> get_nx() || get_ny() != frame_ -> get_ny()){
		cout << " - ERROR!!! - frameClass::multiplyByFrame(): incompatible size with first frame" << endl;
		return ;
	}

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double product = frame_ -> get_point(i, j) * get_point(i, j);
			set_point(i, j, product);
		}
	}

	return ;
}

void frameClass::minusLog(){
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			if(get_point(i, j) > 0.) set_point(i, j, -log(get_point(i, j)));
			else set_point(i, j, 0.);
		}
	}
	return ;
}

void frameClass::initialize(){


	_frame2D = new double*[_nx];
	for(int i=0; i<_nx; i++){
		_frame2D[i] = new double[_ny];
		for(int j=0; j<_ny; j++){
			_frame2D[i][j] = default_double;
		}
	}

	_frame = new double[ _nx*_ny ];
	for(int i = 0 ; i < _nx*_ny ; i++) {
		_frame[i] = default_double;
	}

	initialized = true;

	return;
}


void frameClass::multiplyBy(double val){
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			_frame[(j*_nx)+i] *= val;
		}
	}
	return ;
}

void frameClass::addValue(double val){
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			_frame[(j*_nx)+i] += val;
		}
	}
	return ;
}

void frameClass::addFrame(frameClass *other){
	if(other -> get_nx() != _nx) cout << " - ERROR!!! - frameClass::addFrame(): incompatible x size: " << _nx << " vs " << other -> get_nx() << endl;
	if(other -> get_ny() != _ny) cout << " - ERROR!!! - frameClass::addFrame(): incompatible y size: " << _ny << " vs " << other -> get_ny() << endl;
	if(!initialized) cout << " - ERROR!!! - frameClass::addFrame(): current frame is not yet allocated" << endl;
	if(!other -> is_initialized()) cout << " - ERROR!!! - frameClass::addFrame(): other frame is not yet allocated" << endl;
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			_frame[(j*_nx)+i] += other -> get_point(i, j);
		}
	}
	return ;
}

void frameClass::normalizeToMaximum(){
	double maxVal = -1000000.;
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			if(get_point(i, j) > maxVal) maxVal = get_point(i, j);
		}
	}

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			set_point(i, j, get_point(i, j) / maxVal);
		}
	}

	return ;
}

frameClass *frameClass::rebin(int nBinX, int nBinY){
	if(nBinX < 1 || nBinY < 1){
		cout << " - ERROR!!! - frameClass::rebin(): invalid rebinning value" << endl;
		return NULL;
	}
	int nxNew = _nx / nBinX;
	if((nxNew+1) * nBinX < _nx) nxNew ++;
	int nyNew = _ny / nBinY;
	if((nyNew+1) * nBinY < _ny) nyNew ++;
	frameClass *newFrame = new frameClass(nxNew, nyNew);
	for(int i=0; i<nxNew; i++){
		for(int j=0; j<nyNew; j++){
			double newVal = 0.;
			double countFilled = 0;
			for(int ii=0; ii<nBinX; ii++){
				for(int jj=0; jj<nBinY; jj++){
					int coordX = i*nBinX+ii;
					int coordY = j*nBinY+jj;
					if(coordX<_nx && coordY<_ny){
						newVal += get_point(coordX, coordY);
						countFilled ++;
					}
				}
			}
			newFrame -> set_point(i, j, 1. * newVal * nBinX * nBinY / countFilled);
		}
	}
	return newFrame;
}

frameClass *frameClass::addRow(int row){
	frameClass *newFrame = new frameClass(_nx+1, _ny);
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			if(i < row) newFrame -> set_point(i, j, get_point(i, j));
			else if(i >= row) newFrame -> set_point(i+1, j, get_point(i, j));
		}
	}
	newFrame -> interpolateBadPixels();
	return newFrame;
}

frameClass *frameClass::addColumn(int column){
	frameClass *newFrame = new frameClass(_nx, _ny+1);
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			if(j < column) newFrame -> set_point(i, j, get_point(i, j));
			else if(j >= column) newFrame -> set_point(i, j+1, get_point(i, j));
		}
	}
	newFrame -> interpolateBadPixels();
	return newFrame;
}

void frameClass::rotateCurrent(double angle){
	if(_nx != _ny){
		cout << " - ERROR!!! - frameClass::rotateCurrent(): cannot rotate rectangular matrix: squarize frame before rotation" << endl;
		return ;
	}

	frameClass *newFrame = new frameClass(_nx, _ny);

	double origin = default_double;
	int size = _nx;
	if(size%2 == 0) origin = size/2 - 1;
	else origin = size/2;
	angle = angle * M_PI / 180.;
	for(int j=0; j<size; j++){
		for(int k=0; k<size; k++){
			double Z = get_point(j, k);
			double J = (double)(j - origin);
			double K = (double)(k - origin);
			double R = sqrt(J * J + K * K);
			//////////////
			// cylinder //
			//////////////
			if(R > size / sqrt(2.)) Z = default_double;
			//////////////
			double Phi = default_double;
			bool center = false;
			if(J==0){
				if(K==0){
					center = true;
				}
				else if(K>0) Phi = M_PI/2.;
				else Phi = 1.5*M_PI;
			}
			else if(J>0){
				if(K==0) Phi = 0.;
				else if(K>0) Phi = atan(K/J);
				else Phi = 2*M_PI + atan(K/J);
			}
			else{
				if(K==0) Phi = M_PI;
				else Phi = M_PI + atan(K/J);
			}
			double newPhi = Phi + angle;
			if(newPhi == 2*M_PI) newPhi = 0;
			if(newPhi > 2*M_PI) newPhi-=2*M_PI;
			double newJ = R * cos(newPhi);
			double newK = R * sin(newPhi);
			int newJ_int = newJ;
			int newK_int = newK;
			double decimalJ = fabs(newJ - newJ_int);
			double decimalK = fabs(newK - newK_int);
			if(decimalJ>=0.5){
				if(newJ > 0) newJ_int+=1;
				else if(newJ < 0) newJ_int-=1;
			}
			if(decimalK>=0.5){
				if(newK > 0) newK_int+=1;
				else if(newK < 0) newK_int-=1;
			}
			int newj = newJ_int + origin;
			int newk = newK_int + origin;
			if(center){
				newj = origin;
				newk = origin;
				continue;
			}
			if(newj<0) continue;
			if(newk<0) continue;
			if(newj>=size) continue;
			if(newk>=size) continue;
			newFrame -> set_point(newj, newk, Z);
		}
	}
	newFrame -> interpolateBadPixels();

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			set_point(i, j, newFrame -> get_point(i, j));
		}
	}
	delete newFrame;

	return ;
}

frameClass *frameClass::rotate(double angle){
	if(_nx != _ny){
		cout << " - ERROR!!! - frameClass::rotate(): cannot rotate rectangular matrix: squarize frame before rotation" << endl;
		return NULL;
	}

	frameClass *newFrame = new frameClass(_nx, _ny);

	double origin = default_double;
	int size = _nx;
	if(size%2 == 0) origin = size/2 - 1;
	else origin = size/2;
	angle = angle * M_PI / 180.;
	for(int j=0; j<size; j++){
		for(int k=0; k<size; k++){
			double Z = get_point(j, k);
			double J = (double)(j - origin);
			double K = (double)(k - origin);
			double R = sqrt(J * J + K * K);
			//////////////
			// cylinder //
			//////////////
			if(R > size / sqrt(2.)) Z = default_double;
			//////////////
			double Phi = default_double;
			bool center = false;
			if(J==0){
				if(K==0){
					center = true;
				}
				else if(K>0) Phi = M_PI/2.;
				else Phi = 1.5*M_PI;
			}
			else if(J>0){
				if(K==0) Phi = 0.;
				else if(K>0) Phi = atan(K/J);
				else Phi = 2*M_PI + atan(K/J);
			}
			else{
				if(K==0) Phi = M_PI;
				else Phi = M_PI + atan(K/J);
			}
			double newPhi = Phi + angle;
			if(newPhi == 2*M_PI) newPhi = 0;
			if(newPhi > 2*M_PI) newPhi-=2*M_PI;
			double newJ = R * cos(newPhi);
			double newK = R * sin(newPhi);
			int newJ_int = newJ;
			int newK_int = newK;
			double decimalJ = fabs(newJ - newJ_int);
			double decimalK = fabs(newK - newK_int);
			if(decimalJ>=0.5){
				if(newJ > 0) newJ_int+=1;
				else if(newJ < 0) newJ_int-=1;
			}
			if(decimalK>=0.5){
				if(newK > 0) newK_int+=1;
				else if(newK < 0) newK_int-=1;
			}
			int newj = newJ_int + origin;
			int newk = newK_int + origin;
			if(center){
				newj = origin;
				newk = origin;
				continue;
			}
			if(newj<0) continue;
			if(newk<0) continue;
			if(newj>=size) continue;
			if(newk>=size) continue;
			newFrame -> set_point(newj, newk, Z);
		}
	}
	newFrame -> interpolateBadPixels();

	return newFrame;
}

frameClass *frameClass::rotateSolid(double angle){
	if(_nx != _ny){
		cout << " - ERROR!!! - frameClass::rotateSolid(): cannot rotate rectangular matrix: squarize frame before rotation" << endl;
		return NULL;
	}

	frameClass *newFrame = new frameClass(_nx, _ny);

	double origin = default_double;
	int size = _nx;
	if(size%2 == 0) origin = size/2 - 1;
	else origin = size/2;
	angle = angle * M_PI / 180.;
	for(int j=0; j<size; j++){
		for(int k=0; k<size; k++){
			double Z = get_point(j, k);
			double J = (double)(j - origin);
			double K = (double)(k - origin);
			double R = sqrt(J * J + K * K);
			//////////////
			// cylinder //
			//////////////
			if(R > size / sqrt(2.)) Z = default_double;
			//////////////
			double Phi = default_double;
			bool center = false;
			if(J==0){
				if(K==0){
					center = true;
				}
				else if(K>0) Phi = M_PI/2.;
				else Phi = 1.5*M_PI;
			}
			else if(J>0){
				if(K==0) Phi = 0.;
				else if(K>0) Phi = atan(K/J);
				else Phi = 2*M_PI + atan(K/J);
			}
			else{
				if(K==0) Phi = M_PI;
				else Phi = M_PI + atan(K/J);
			}
			double newPhi = Phi + angle;
			if(newPhi == 2*M_PI) newPhi = 0;
			if(newPhi > 2*M_PI) newPhi-=2*M_PI;
			double newJ = R * cos(newPhi);
			double newK = R * sin(newPhi);
			int newJ_int = newJ;
			int newK_int = newK;
			double decimalJ = fabs(newJ - newJ_int);
			double decimalK = fabs(newK - newK_int);
			if(decimalJ>=0.5){
				if(newJ > 0) newJ_int+=1;
				else if(newJ < 0) newJ_int-=1;
			}
			if(decimalK>=0.5){
				if(newK > 0) newK_int+=1;
				else if(newK < 0) newK_int-=1;
			}
			int newj = newJ_int + origin;
			int newk = newK_int + origin;
			if(center){
				newj = origin;
				newk = origin;
				continue;
			}
			if(newj<0) continue;
			if(newk<0) continue;
			if(newj>=size) continue;
			if(newk>=size) continue;
			newFrame -> set_point(newj, newk, Z);
		}
	}

	newFrame -> interpolateBadPixelsSolid();

	return newFrame;
}

frameClass *frameClass::squarize(){
	if(_nx == _ny) return NULL;
	int size = default_int;
	if(_nx>_ny) size = _nx;
	else size = _ny;
	double stepX = 1. * _nx / size;
	double stepY = 1. * _ny / size;
	frameClass *newFrame = new frameClass(size, size);
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			double ii = 1. * _nx * i / size;
			double jj = 1. * _ny * j / size;
			double iii = ii + stepX;
			double jjj = jj + stepY;
			int II = ii;
			int JJ = jj;
			int III = iii;
			int JJJ = jjj;
			double iiFrac = (III - ii) / stepX;
			double jjFrac = (JJJ - jj) / stepY;
			double iiiFrac = (iii - III) / stepX;
			double jjjFrac = (jjj - JJJ) / stepY;
			if(stepX == 1.){
				if(JJ == JJJ) newFrame -> set_point(i, j, get_point(II, JJ));
				else{
					if(JJJ == _ny) newFrame -> set_point(i, j, get_point(II, JJ));
					else newFrame -> set_point(i, j, get_point(II, JJ) * jjFrac + get_point(II, JJJ) * jjjFrac);
				}
			}
			else if(stepY == 1.){
				if(II == III) newFrame -> set_point(i, j, get_point(II, JJ));
				else{
					if(III == _nx) newFrame -> set_point(i, j, get_point(II, JJ));
					else newFrame -> set_point(i, j, get_point(II, JJ) * iiFrac + get_point(III, JJ) * iiiFrac);
				}
			}
		}
	}
	return newFrame;
}

void frameClass::interpolateBadPixels(){

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double val = get_point(i, j);
			if(isnan(val) || isinf(val)) set_point(i, j, default_double);
		}
	}
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double val = get_point(i, j);
			if(val == default_double){
				int countNeighbours = 0;
				double average = 0.;
				double neighbour_val = default_double;
				// left block
				if(i!=0){
					if(j!=0){
						neighbour_val = get_point(i-1, j-1);
						if(neighbour_val!=default_double){
							++countNeighbours;
							average += neighbour_val;
						}
					}
					if(j!=_ny-1){
						neighbour_val = get_point(i-1, j+1);
						if(neighbour_val!=default_double){
							++countNeighbours;
							average += neighbour_val;
						}
					}
					neighbour_val = get_point(i-1, j);
					if(neighbour_val!=default_double){
						++countNeighbours;
						average += neighbour_val;
					}
				}
				// right block
				if(i!=_nx-1){
					if(j!=0){
						neighbour_val = get_point(i+1, j-1);
						if(neighbour_val!=default_double){
							++countNeighbours;
							average += neighbour_val;
						}
					}
					if(j!=_ny-1){
						neighbour_val = get_point(i+1, j+1);
						if(neighbour_val!=default_double){
							++countNeighbours;
							average += neighbour_val;
						}
					}
					neighbour_val = get_point(i+1, j);
					if(neighbour_val!=default_double){
						++countNeighbours;
						average += neighbour_val;
					}
				}
				// central block
				if(j!=0){
					neighbour_val = get_point(i, j-1);
					if(neighbour_val!=default_double){
						++countNeighbours;
						average += neighbour_val;
					}
				}
				if(j!=_ny-1){
					neighbour_val = get_point(i, j+1);
					if(neighbour_val!=default_double){
						++countNeighbours;
						average += neighbour_val;
					}
				}
				// averaging
				if(countNeighbours != 0) average = average / countNeighbours;
				set_point(i, j, average);
			}
		}
	}
	return ;
}

void frameClass::interpolateBadPixelsSolid(){
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double val = get_point(i, j);
			if(isnan(val) || isinf(val)) set_point(i, j, default_double);
		}
	}
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			double val = get_point(i, j);
			if(val == default_double){
				vector<double> neighbour_val;
				// left block
				if(i!=0){
					if(j!=0){
						if(get_point(i-1, j-1) != default_double) neighbour_val.push_back(get_point(i-1, j-1));
					}
					if(j!=_ny-1){
						if(get_point(i-1, j+1) != default_double) neighbour_val.push_back(get_point(i-1, j+1));
					}
					if(get_point(i-1, j) != default_double) neighbour_val.push_back(get_point(i-1, j));
				}
				// right block
				if(i!=_nx-1){
					if(j!=0){
						if(get_point(i+1, j-1) != default_double) neighbour_val.push_back(get_point(i+1, j-1));
					}
					if(j!=_ny-1){
						if(get_point(i+1, j+1) != default_double) neighbour_val.push_back(get_point(i+1, j+1));
					}
					if(get_point(i+1, j) != default_double) neighbour_val.push_back(get_point(i+1, j));
				}
				// central block
				if(j!=0){
					if(get_point(i, j-1) != default_double) neighbour_val.push_back(get_point(i, j-1));
				}
				if(j!=_ny-1){
					if(get_point(i, j+1) != default_double) neighbour_val.push_back(get_point(i, j+1));
				}
				// determining most present value
				// if not, putting random value
				set_point(i, j, getMostOccurentValue(neighbour_val));
				neighbour_val.clear();
			}
		}
	}
	return ;
}

double frameClass::getMostOccurentValue(vector<double> array){

	vector<double> val;
	vector<int> count;

	for(unsigned int i=0; i<array.size(); i++){

		bool found = false;
		for(unsigned int j=0; j<val.size(); j++){
			if(array[i] == val[j]){
				found = true;
				count[j] ++;
				break;
			}
		}
		if(found == false){
			val.push_back(array[i]);
			count.push_back(1);
		}

	}

	int indexMax = default_int;
	int countMax = -1;
	for(unsigned int i=0; i<count.size(); i++){
		if(count[i] > countMax){
			countMax = count[i];
			indexMax = i;
		}
	}

	val.clear();
	count.clear();

	if(indexMax == default_int) return 0;
	else return val[indexMax];
}

void frameClass::interpolateBadPixels(double bad_pixel_value){

	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			if(get_point(i, j) == bad_pixel_value) set_point(i, j, default_double);
		}
	}

	interpolateBadPixels();

	return ;
}

void frameClass::setTo_inCircle(double radius, double val){

	double originX = default_double;
	if(_nx%2 == 0) originX = _nx/2 - 1;
	else originX = _nx/2;
	double originY = default_double;
	if(_ny%2 == 0) originY = _ny/2 - 1;
	else originY = _ny/2;

	for(int j=0; j<_nx; j++){
		for(int k=0; k<_ny; k++){
			double J = (double)(j - originX);
			double K = (double)(k - originY);
			double R = sqrt(J * J + K * K);
			if(R < radius){
				set_point(j, k, val);
			}
		}
	}

	return ;
}

frameClass *frameClass::makeCircle_in(double valueOutsideRadius){

	double originX = default_double;
	if(_nx%2 == 0) originX = _nx/2 - 1;
	else originX = _nx/2;
	double originY = default_double;
	if(_ny%2 == 0) originY = _ny/2 - 1;
	else originY = _ny/2;

	int size = default_int;
	if(_nx<_ny) size = _nx;
	else size = _ny;
	double origin = default_double;
	if(size%2 == 0) origin = size/2 - 1;
	else origin = size/2;
	frameClass *newFrame = new frameClass(size, size);
	newFrame -> set_to(valueOutsideRadius);
	for(int j=0; j<_nx; j++){
		for(int k=0; k<_ny; k++){
			double J = (double)(j - originX);
			double K = (double)(k - originY);
			double R = sqrt(J * J + K * K);
			if(R < size / 2){
				int newJ = J + origin;
				int newK = K + origin;
				newFrame -> set_point(newJ, newK, get_point(j, k));
			}
		}
	}
	newFrame -> interpolateBadPixels();
	return newFrame;
}

frameClass *frameClass::makeCircle_out(double valueOutsideRadius){
	double originX = default_double;
	if(_nx%2 == 0) originX = _nx/2 - 1;
	else originX = _nx/2;
	double originY = default_double;
	if(_ny%2 == 0) originY = _ny/2 - 1;
	else originY = _ny/2;
	int size = sqrt(_nx * _nx + _ny * _ny);
	double origin = default_double;
	if(size%2 == 0) origin = size/2 - 1;
	else origin = size/2;
	frameClass *newFrame = new frameClass(size, size);
	newFrame -> set_to(valueOutsideRadius);
	for(int j=0; j<_nx; j++){
		for(int k=0; k<_ny; k++){
			double J = (double)(j - originX);
			double K = (double)(k - originY);
			int newJ = J + origin;
			int newK = K + origin;
			newFrame -> set_point(newJ, newK, get_point(j, k));
		}
	}
	return newFrame;
}

frameClass *frameClass::selection(int bottomLeftX, int bottomLeftY, int topRightX, int topRightY){
	if(bottomLeftX < 0 || bottomLeftY < 0){
		cout << " - ERROR!!! - frameClass::selection(): invalid selection coordinates" << endl;
		return NULL;
	}
	if(topRightX > _nx || topRightY > _ny){
		cout << " - ERROR!!! - frameClass::selection(): invalid selection coordinates" << endl;
		return NULL;
	}
	int nxNew = topRightX - bottomLeftX;
	int nyNew = topRightY - bottomLeftY;
	if(nxNew <= 0 || nyNew <= 0){
		cout << " - ERROR!!! - frameClass::selection(): invalid selection coordinates" << endl;
		return NULL;
	}
	frameClass *newFrame = new frameClass(nxNew, nyNew);
	for(int i=0; i<nxNew; i++){
		for(int j=0; j<nyNew; j++){
			newFrame -> set_point(i, j, get_point(bottomLeftX + i, bottomLeftY + j));
		}
	}
	return newFrame;
}

frameClass *frameClass::clone(){
	frameClass *newFrame = new frameClass(_nx, _ny);
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			newFrame -> set_point(i, j, get_point(i, j));
		}
	}
	return newFrame;
}

frameClass *frameClass::flipHorizontally(){
	frameClass *newFrame = new frameClass(_nx, _ny);
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			newFrame -> set_point(i, j, get_point(_nx-i-1, j));
		}
	}
	return newFrame;
}

frameClass *frameClass::flipVertically(){
	frameClass *newFrame = new frameClass(_nx, _ny);
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			newFrame -> set_point(i, j, get_point(i, _ny-j-1));
		}
	}
	return newFrame;
}

void frameClass::copyFromFrame(frameClass *other){
	if(other -> _nx != _nx || other -> _ny != _ny){
		cout << " - ERROR!!! - frameClass::copyFromFrame(): incompatible frame size" << endl;
		return ;
	}
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			set_point(i, j, other -> get_point(i, j));
		}
	}
	return ;
}

void frameClass::computeAverageFromCircle(double radius, double &mean, double &RMS){

	mean = 0.;
	RMS = 0.;

	int N = 0;
	double xC = _nx / 2.;
	double yC = _ny / 2.;
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			double dX = i - xC;
			double dY = j - yC;
			double dist = sqrt(dX * dX + dY * dY);
			if(dist < radius){
				N ++;
				mean += get_point(i, j);
				RMS += (get_point(i, j) * get_point(i, j));
			}
		}
	}

	mean /= N; // E[x]
	RMS /= N; // E[x^2]

	RMS -= (mean * mean); // E[x^2] - (E[x])^2

	RMS = sqrt(RMS);

	return ;
}

double frameClass::get_minimum(){

	double min = 100000000000.;
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			if(get_point(i, j) < min) min = get_point(i, j);
		}
	}

	return min;
}

double frameClass::get_maximum(){

	double max = -100000000000.;
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			if(get_point(i, j) > max) max = get_point(i, j);
		}
	}

	return max;
}

bool frameClass::readFromFile_ASCIIMatrix(const char* fileName){
	if(!initialized) {
		cout << " - ERROR!!! - frameClass::readFromFile_ASCIIMatrix(): frame is not yet allocated" << endl;
		return false;
	}
	ifstream file;
	file.open(fileName);
	if(file == NULL) {
		cout << " - ERROR!!! - frameClass::readFromFile_ASCIIMatrix(): cannot open file " << fileName << endl;
		return false;
	}
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			file >> _frame[(j*_nx)+i];
		}
	}
	file.close();
	return true;
}

void frameClass::writeToFile_ASCIIMatrix(const char* fileName) const{
	if(!initialized) cout << " - ERROR!!! - frameClass::writeToFile_ASCIIMatrix(): frame is not yet allocated" << endl;
	ofstream file;
	file.open(fileName);
	if(file == NULL) cout << " - ERROR!!! - frameClass::writeToFile_ASCIIMatrix(): cannot open file " << fileName << endl;
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			file << _frame[(j*_nx)+i] << " ";
		}
		file << endl;
	}
	file.close();
	return ;
}

void frameClass::readFromFile_ASCIITable(const char* fileName){
	if(!initialized) cout << " - ERROR!!! - frameClass::readFromFile_ASCIITable(): frame is not yet allocated" << endl;
	ifstream file;
	file.open(fileName);
	if(file == NULL) cout << " - ERROR!!! - frameClass::readFromFile_ASCIITable(): cannot open file " << fileName << endl;
	set_to(0.);
	int x = default_int;
	int y = default_int;
	double val = default_double;
	while(file >> x){
		file >> y;
		file >> val;
		_frame[(y*_nx)+x] = val;
	}
	file.close();
	return ;
}

void frameClass::writeToFile_ASCIITable(const char* fileName) const{
	if(!initialized) cout << " - ERROR!!! - frameClass::writeToFile_ASCIITable(): frame is not yet allocated" << endl;
	ofstream file;
	file.open(fileName);
	if(file == NULL) cout << " - ERROR!!! - frameClass::writeToFile_ASCIITable(): cannot open file " << fileName << endl;
	for(int i=0; i<_nx; i++){
		for(int j=0; j<_ny; j++){
			if(_frame[(j*_nx)+i] == 0.) continue;
			file << i << "\t" << j << "\t" << _frame[(j*_nx)+i] << endl;
		}
	}
	file.close();
	return ;
}

void frameClass::readFromFile_binaryTable(const char* fileName){
	if(!initialized) cout << " - ERROR!!! - frameClass::readFromFile_binaryTable(): frame is not yet allocated" << endl;
	ifstream file;
	file.open(fileName, ios::in | ios::binary);
	if(file == NULL) cout << " - ERROR!!! - frameClass::readFromFile_binaryTable(): cannot open file " << fileName << endl;
	while(!file.eof()){
		int dumb = -1;
		int x1 = file.get();
		int x2 = file.get();
		dumb = file.get();
		dumb = file.get();
		int y1 = file.get();
		int y2 = file.get();
		dumb = file.get();
		dumb = file.get();
		int c1 = file.get();
		int c2 = file.get();
		int x = x1 + 256 * x2;
		int y = y1 + 256 * y2;
		int c = c1 + 256 * c2;
		if(file.eof()) break;
		if(x>get_nx()) cout << " - ERROR!!! - frameClass::readFrameFile_binaryTable(): x index value (" << x << ") out of bounds (0" << " - " << get_nx() << ")" << endl;
		if(y>get_ny()) cout << " - ERROR!!! - frameClass::readFrameFile_binaryTable(): y index value (" << y << ") out of bounds (0" << " - " << get_ny() << ")" << endl;
		_frame[(y*_nx)+x] = c;
		c = dumb;
	}
	file.close();

	return ;
}

void frameClass::get_row(int index, double *row, int arraySize){

	if(index >=get_nx()){
		cout << " - ERROR!!! - frameClass::get_row(): index (" << index << ") out of bound (" << get_nx() << ")" << endl;
		return ;
	}
	if(arraySize != get_nx()){
		cout << " - ERROR!!! - frameClass::get_row(): array size (" << arraySize << ") is incompatible (" << get_nx() << ")" << endl;
		return ;
	}

	for(int i=0; i<get_nx(); i++){
		//row[i] = frame[i][index];
		row[i] = _frame[(index*_nx) + i];
	}

	return ;
}

void frameClass::get_column(int index, double *column, int arraySize){

	if(index >=get_ny()){
		cout << " - ERROR!!! - frameClass::get_column(): index (" << index << ") out of bound (" << get_ny() << ")" << endl;
		return ;
	}
	if(arraySize != get_ny()){
		cout << " - ERROR!!! - frameClass::get_column(): array size (" << arraySize << ") is incompatible (" << get_ny() << ")" << endl;
		return ;
	}

	for(int i=0; i<get_ny(); i++){
		//column[i] = frame[index][i];
		column[i] = _frame[(i*_nx) + index];
	}

	return ;
}

void frameClass::getStats(double &mean, double &rms, int &entries){

	entries = 0;

	mean = 0.;
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			mean += get_point(i, j);
			entries ++;
		}
	}
	mean /= (double)entries;

	rms = 0.;
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			double res = (get_point(i, j) - mean);
			rms += (res * res);
		}
	}

	rms /= (double)entries;
	rms = sqrt(rms);

	return ;
}

void frameClass::getStatsInRange(double &mean, double &rms, int &entries, double min, double max){

	entries = 0;

	mean = 0.;
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			if(get_point(i, j) < min || get_point(i, j) > max) continue;
			mean += get_point(i, j);
			entries ++;
		}
	}
	mean /= (double)entries;

	rms = 0.;
	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){
			if(get_point(i, j) < min || get_point(i, j) > max) continue;
			double res = (get_point(i, j) - mean);
			rms += (res * res);
		}
	}

	rms /= (double)entries;
	rms = sqrt(rms);

	return ;
}

frameClass *frameClass::applyGaussFilter(int sigma){

	bool debug = false;

	int size = 2 * 2 * sigma + 1;
	double mu = 2 * sigma;
	double kernel[size][size];
	double norm = 0.;
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			kernel[i][j] = exp(-((i - mu) * (i - mu) + (j - mu) * (j - mu)) / (2 * sigma * sigma));
			norm += kernel[i][j];
		}
	}
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			kernel[i][j] /= norm;
		}
	}

	if(debug){
		cout << " - frameClass::applyGaussFilter(): Kernel:" << endl;
		for(int i=0; i<size; i++){
			cout << "\t";
			for(int j=0; j<size; j++){
				cout << kernel[i][j] << " ";
			}
			cout << endl;
		}
	}

	frameClass *filtered = new frameClass(get_nx(), get_ny());
	filtered -> set_to(0.);

	for(int i=0; i<get_nx(); i++){
		for(int j=0; j<get_ny(); j++){

			double val = 0.;

			norm = 0.;
			for(int kx=0; kx<size; kx++){
				int kPosX = kx - (int)mu;
				if(i + kPosX < 0 || i + kPosX > get_nx() - 1) continue;
				for(int ky=0; ky<size; ky++){
					int kPosY = ky - (int)mu;
					if(j + kPosY < 0 || j + kPosY > get_ny() - 1) continue;
					val += kernel[kx][ky] * get_point(i + kPosX, j + kPosY);
					norm += kernel[kx][ky];
				}
			}

			filtered -> set_point(i, j, val / norm);

		}
	}

	return filtered;
}

