#include "TH1D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TF1.h"
#include "TFile.h"
#include <TStopwatch.h>

using namespace std;
double DeltaPhi(double phi1, double phi2); // relative angle

int main(int argc, char **argv)
{

	TROOT root("flow","run mc");

	if ( argc<3 ) {
		cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		cout<<"+  "<<argv[0]<<" <outputFile> <Nevt> <random seed> "<<endl;
		cout << endl << endl;
		exit(1);
	}

	// CONSTANT
	char *outFile = argv[1];
	Int_t numberEvents= atoi(argv[2]);
	Int_t random_seed = atoi(argv[3]);
	TFile *fout = new TFile(outFile,"recreate");

	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout << "Outfile : "<< outFile << endl;
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	const double pi = TMath::Pi();
	// Parameter
	TRandom *myRandom = new TRandom(random_seed);
	const int NS=2;
	TH1D *hphi[NS];
	for (int i=0;i<NS;i++) hphi[i] = new TH1D(Form("h_phi%02d",i),Form("h_phi%02d",i), 300, 0, 2*pi); 
        TH1D *hDphiDiff = new TH1D("hDphiDiff","hDphi", 300, 0, 3*pi);
        TH1D *hDphiRS = new TH1D("hDphiRS","hDphiRS", 300, 0, 3*pi);
        TH1D *hatan2 = new TH1D("hatan2","hatan2", 300, -2*pi, 2*pi);
	cout<<"Events="<<numberEvents <<" RNDM seed "<< random_seed << endl;

	TStopwatch timer;
	timer.Start();

	double phi[NS]={-999};
	int ieout = numberEvents/10;
	if (ieout<1) ieout=1;
	for(int ievt =0; ievt<numberEvents; ievt++){ 
		if(ievt%ieout == 0) cout << " Event: " << ievt << " " << ievt*100.0/numberEvents <<"% "<< endl;
		for (int i=0;i<NS;i++) {
			phi[i] = myRandom->Uniform(0, 2*pi);
			hphi[i]->Fill(phi[i]);
		}
		hDphiDiff->Fill(phi[0]-phi[1]);
		hDphiRS->Fill(DeltaPhi(phi[0],phi[1]));
		hatan2->Fill(atan2(sin(phi[0]-phi[1]), cos(phi[0]-phi[1])));

	} // end of the events
	cout <<"Total number of event processed = "<< numberEvents << endl;
	cout <<"Writing histograms into a root file"<<endl;
	fout->cd();

	for (int i=0;i<NS;i++) hphi[i]->Write();
	hDphiDiff->Write();
	hDphiRS->Write();
	hatan2->Write();

	fout->Close();
	timer.Print();
	cout <<"Sucessfully finished."<<endl;
}
double DeltaPhi(double phi1, double phi2) {
  // dphi
  double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
  return res>0 ? res : 2.*TMath::Pi()+res ;
}
