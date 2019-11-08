#include <iostream>
#include <stdlib.h>

#include "libraries/frameClass.cc"
#include "libraries/frameClass.hh"

#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iterator>
#include <vector>

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream

#define __line_width_guess 12000

using namespace std;

void applyStyle(){
	gStyle -> SetPaperSize(20, 20);
	gStyle -> SetPadTopMargin(0);
	gStyle -> SetPadRightMargin(0);
	gStyle -> SetPadBottomMargin(0);
	gStyle -> SetPadLeftMargin(0);
	//  gStyle -> SetPalette(51);
	gStyle -> SetOptStat(0);
	//  gStyle -> SetOptTitle(0);
	gStyle -> SetOptFit(0);
	gStyle->SetTitleW(0.6);

}

int plotFrame(const char* fileName, const char* imageName, int sliceNumber, int size, double min = 0., double max = 0., bool logscale = false) {
	return 0;
}

int plotFrame(const char* fileName, const char* imageName, int sliceNumber, double min = 0., double max = 0., bool logscale = false) {

	applyStyle();

	//std::ifstream 	file("/home/wesley/storage/RJ45_CT/data1_reco/reconstruction/OSEM/frame_288.txt");
	//std::cout << std::distance(std::istream_iterator<double>(file), std::istream_iterator<double>());

	ifstream ifs(fileName, ifstream::in);
	char c;
	int nCR = 0;
	while ( ifs.good() ) {

		c = ifs.get();
		//cout << temp << endl;
		if ( c == 0xa ) nCR++;

	}
	ifs.close();

	// Rewind
	ifs.open(fileName, ifstream::in);
	char tempLine[__line_width_guess];
	ifs.getline(tempLine, __line_width_guess);
	//cout << tempLine << endl;
	int itr = 0;
	int rowCntr = 0;
	while( tempLine[itr] != 0xa ) {

		if( tempLine[itr] == ' ' ) rowCntr++;

		itr++;
	}

	// number of lines in y
	int sizey = nCR;
	// number of lines in x
	int sizex = rowCntr;
	// Report
	cout << "sizex = " << sizex << ", sizey = " << sizey << endl;

	int size = sizey;


	frameClass *frame = new frameClass(size, size);
	frame -> readFromFile_ASCIIMatrix(fileName);

	char title[1000];
	sprintf(title, "Frame %03d", sliceNumber);

	TH2F *hFrame = new TH2F("frame", title, size, 0, size, size, 0, size);
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			hFrame -> SetBinContent(i, j, frame -> get_point(i, j));
		}
	}

	TCanvas *cFrame = new TCanvas("cFrame", "cFrame", 0, 0, 700, 700);
	if(logscale) cFrame -> SetLogz();
	if(min==0. && max==0.){
		hFrame -> SetMinimum(frame -> get_minimum());
		hFrame -> SetMaximum(frame -> get_maximum());
	}
	else{
		hFrame -> SetMinimum(min);
		hFrame -> SetMaximum(max);
	}
	hFrame -> GetXaxis() -> SetTitle("x");
	hFrame -> GetYaxis() -> SetTitle("y");
	hFrame -> GetXaxis() -> SetTitleSize(0.05);
	hFrame -> GetYaxis() -> SetTitleSize(0.05);
	hFrame -> GetYaxis() -> SetTitleOffset(1.44);
	hFrame -> Draw("colz");

	cFrame -> SaveAs(imageName);

	//delete frame;
	//delete hFrame;
	//delete cFrame;

	return 0;
}

int drawImages(int sliceMin, int sliceMax, int sliceStep, int size, double min, double max, const char *fileNameFormat, const char*imageNameFormat, bool logscale = false){

	for(int sliceNumber = sliceMin; sliceNumber<=sliceMax; sliceNumber+=sliceStep){
		char fileName[1000];
		sprintf(fileName, fileNameFormat, sliceNumber);
		char imageName[1000];
		sprintf(imageName, imageNameFormat, sliceNumber);
		plotFrame(fileName, imageName, sliceNumber, size, min, max, logscale);
	}

	return 0;
}

int drawImagesCT(int sliceMin, int sliceMax, int sliceStep, int size, double min, double max, const char *fileNameFormat, const char*imageNameFormat, bool logscale = false){

	int count = 0;
	for(int sliceNumber = sliceMin; sliceNumber<=sliceMax; sliceNumber+=sliceStep){
		char fileName[1000];
		sprintf(fileName, fileNameFormat, sliceNumber);
		char imageName[1000];
		sprintf(imageName, imageNameFormat, count);
		plotFrame(fileName, imageName, sliceNumber, size, min, max, logscale);
		count ++;
	}

	return 0;
}

int demo(){

	int sliceMin = 0;
	int sliceMax = 63;
	int sliceStep = 1;
	int size = 65;
	int min = 0;
	int max = 800.;
	const char *fileNameFormatFBP = "demo/reconstruction/FBP/frame_%d.txt";
	const char *imageNameFormatFBP = "demo/images/FBP/frame_%03d.png";
	bool logscale = false;
	drawImages(sliceMin, sliceMax, sliceStep, size, min, max, fileNameFormatFBP, imageNameFormatFBP, logscale);

	max = 8.;
	const char *fileNameFormatOSEM = "demo/reconstruction/OSEM/frame_%d.txt";
	const char *imageNameFormatOSEM = "demo/images/OSEM/frame_%03d.png";
	drawImages(sliceMin, sliceMax, sliceStep, size, min, max, fileNameFormatOSEM, imageNameFormatOSEM, logscale);

	sliceMin = 0;
	sliceMax = 288000;
	sliceStep = 3200;
	size = 64;
	min = -3;
	max = 260;
	const char *fileNameFormatCT = "demo/data/raw/frame_%d.txt";
	const char *imageNameFormatCT = "demo/images/CT/frame_%03d.png";
	drawImagesCT(sliceMin, sliceMax, sliceStep, size, min, max, fileNameFormatCT, imageNameFormatCT, logscale);

	return 0;
}
