int plotSlice(const char *fileName,
	      const unsigned int size){

  ifstream file(fileName);
  if(!file){
    cout << __PRETTY_FUNCTION__ << ": ERROR!!! - cannot open file " << fileName << endl;
    return 1;
  }

  TH2F *h2 = new TH2F("h2", "h2", size, 0, size, size, 0, size);
  double val = -1.;
  for(unsigned int i=0; i<size; i++){
    for(unsigned int j=0; j<size; j++){
      file >> val;
      h2 -> SetBinContent(i+1, j+1, val);
    }
  }

  file.close();

  gStyle -> SetPaperSize(20, 20);
  gStyle -> SetPadTopMargin(0.15);
  gStyle -> SetPadRightMargin(0.2);
  gStyle -> SetPadBottomMargin(0.18);
  gStyle -> SetPadLeftMargin(0.13);
  gStyle -> SetOptTitle(0);
  gStyle -> SetOptFit(1111);
  gStyle -> SetTitleW(0.6);
  gStyle -> SetOptStat(0);
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetNumberContours(255);
  gStyle -> SetPalette(55);

  TCanvas *cc = new TCanvas("cc", "cc", 0, 0, 1000, 1000);
  h2 -> Draw("colz");
  
  return 0;
}
