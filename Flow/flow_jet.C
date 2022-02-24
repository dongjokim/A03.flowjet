#include <Pythia8/Pythia.h>
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
#include <Pythia8/Pythia.h>
#include <TStopwatch.h>



using namespace std;
using namespace Pythia8; 


int main(int argc, char **argv)
{

  TROOT root("flow_jet","run mc");

  if ( argc<7 ) {
        cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
        cout<<"+  "<<argv[0]<<" <outputFile> <N_ch in midrapidity><eta_min><eta_max><pTHardMin> <random seed> <Nevt>"<<endl;
        cout << endl << endl;
        exit(1);
  }

  // CONSTANT
  char *outFile = argv[1];
  Int_t nParMidRap  = atoi(argv[2]); // vary with centrality replaced by raw multiplicity
  Float_t eta_min   = atof(argv[3]); 
  Float_t eta_max   = atof(argv[4]);
  Float_t pTHatMin =  atof(argv[5]);
  Int_t random_seed = atoi(argv[6]);
  Int_t numberEvents= atoi(argv[7]);
  TFile *fout = new TFile(outFile,"recreate");

  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout << "Outfile : "<< outFile << endl;
  cout << "Inputs : "<<endl;
  cout <<"\t nParMidRap = "<< nParMidRap << endl;
  cout <<"\t pTHardMin = "<< pTHatMin << endl;
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  const double pi = TMath::Pi();
  // Parameter
  TRandom *myRandom = new TRandom(random_seed);

  const int NH = 4; // 1: first 2:2nd
  double input_vn[NH] = {100, 0.0,0.3,0.2}; //[0] is underlying

  const int NBINS=200;
  double LogBinsX[201], LimL=0.1, LimH=110;
  double logBW = (log(LimH)-log(LimL))/NBINS;
  for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);
  TH1D *hChargedPtMid = new TH1D("hChargedPtMid"," ", NBINS, LogBinsX ); hChargedPtMid->Sumw2();
  TH1D *hChargedPt= new TH1D("hChargedPt"," ", NBINS, LogBinsX ); hChargedPt->Sumw2();
  TH1D *hChargedEta = new TH1D("hChargedEta"," ", 200, -10,10 ); hChargedEta->Sumw2();

  LimL=0.1;  LimH=300;
  logBW = (log(LimH)-log(LimL))/NBINS;
  for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);
  TH1D *hJetPt    = new TH1D("hJetPt"," ", NBINS, LogBinsX ); hJetPt->Sumw2();
  /// PYTHIA 
  //---------------------
  //Pythia initialization 
  //---------------------
  Pythia pythia;   // Generator.
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Read in commands from external file.
  pythia.readFile("config.cmnd");   

  // Extract settings to be used in the main program.
  //int    numberEvents  = pythia.mode("Main:numberOfEvents");
  //int    nList   = pythia.mode("Main:numberToList");
  int    nShow   = pythia.mode("Main:timesToShow");
  bool   showCS  = pythia.flag("Main:showChangedSettings");
  bool   showCPD = pythia.flag("Main:showChangedParticleData");

  cout<<"Events="<<numberEvents <<" RNDM seed "<< random_seed << endl;
  cout<<"PTHARDMin = "<< pTHatMin << endl;

  pythia.readString(Form("PhaseSpace:pTHatMin = %f", pTHatMin ));
  pythia.readString(Form("Random:seed=%02d",random_seed));
  // Initialize. Beam parameters set in .cmnd file.
  pythia.init();

  // List changed data.
  if (showCS)  pythia.settings.listChanged();
  if (showCPD) pdt.listChanged();

  TStopwatch timer;
  timer.Start();

  double psi_in[NH]={-999};
  int ieout = numberEvents/10;
  if (ieout<1) ieout=1;
  for(int ievt =0; ievt<numberEvents; ievt++){ 
    if(nShow > 0 && ievt%ieout == 0) cout << " Event: " << ievt << " " << ievt*100.0/numberEvents <<"% "<< endl;

    // PYTHIA
    if(!pythia.next()){
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    // change ..
    psi_in[0] = 0.;
    psi_in[1] = 0.;
    psi_in[2] = myRandom->Uniform(-pi/2, pi/2); // generate random event plane psi from -pi/2 ~ pi/2
    psi_in[3] = myRandom->Uniform(-pi/3,  pi/3);  // random psi from -pi/4 ~ pi/4

    TF1 *dNdphi_bulk = new TF1("dNdphi_bulk", "[0]*(1+2*[1]*(x- [4]) +2*[2]*cos(2*(x-[5]))+2*[3]*cos(3*(x-[6])))",0 ,2* pi);

    dNdphi_bulk->SetParameter(0, input_vn[0]); // underyn const(?) - 1 over 2Pi dN..... and so on.
    dNdphi_bulk->SetParameter(1, input_vn[1]);  //v1  - ignore directed flow
    dNdphi_bulk->SetParameter(2, input_vn[2]);  // v2
    dNdphi_bulk->SetParameter(3, input_vn[3]);  // v3
    dNdphi_bulk->SetParameter(4, psi_in[1]); // v1 event plane
    dNdphi_bulk->SetParameter(5, psi_in[2]); // v2 event plane
    dNdphi_bulk->SetParameter(6, psi_in[3]); // v3 event palne

    //generatation particle with random angle phi(which followed dNdphi folmuar)
    vector <double> phi_bulk;
    vector<double> phi_jet_mid;
    for( int ip=0;ip<nParMidRap;ip++ ){
      phi_bulk.push_back(dNdphi_bulk->GetRandom());
      phi_jet_mid.push_back(dNdphi_bulk->GetRandom() );
    }
    
    // Adding jets
    for ( int ip=0;ip<event.size();ip++ ){
      Particle & p = event[ip];
      if( !p.isFinal() || !p.isCharged() ) continue;
      if( p.eta() > eta_min && p.eta() < eta_max ){ phi_jet_mid.push_back( p.phi() ); hChargedPtMid->Fill(p.pT());}
      hChargedPt->Fill(p.pT());
      hChargedEta->Fill(p.eta());
    }
   
    delete dNdphi_bulk;

  } // end of the events
  cout <<"Total number of event processed = "<< numberEvents << endl;
  cout <<"Writing histograms into a root file"<<endl;

  fout->cd();

  hChargedPt->Write();hChargedPtMid->Write();
  hChargedEta->Write();
  fout->Close();
  timer.Print();
  cout <<"Sucessfully finished."<<endl;
}
