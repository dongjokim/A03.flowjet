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
#include <Pythia8/Pythia.h>
#include <TStopwatch.h>

#include "fastjet/ClusterSequence.hh"
#include <iostream>

#include <TClonesArray.h>

using namespace std;
using namespace Pythia8; 
using namespace fastjet;

class MyUserInfo : public PseudoJet::UserInfoBase{
  public:
    // default ctor
    MyUserInfo(const int & pdg_id_in,const vector<bool> &pType) :
      _pdg_id(pdg_id_in){ ispid = pType;}

    /// access to the PDG id
    int pdg_id() const { return _pdg_id;}
    void PrintIsPid() const { 
      for(unsigned int i=0;i<ispid.size();i++) {
       cout << ispid.at(i)<<":";
      }
      cout << endl;
    }
    bool IsType(int i) const { return ispid[i];}

  protected:
    int _pdg_id;         // the associated pdg id
    vector<bool> ispid;
};

int main(int argc, char **argv) {

  if(argc<4){
    cout<<"usage: pythia pythia.config card.input <output.root> [random_seed]"<<endl;exit(1);//usage of the code, requieres 4 inputs
  }
  TStopwatch timer; 
  timer.Start();   

  char* pythiaconfig  = argv[1];
  double pTHatMin     = atof(argv[2]);
  double pTHatMax     = atof(argv[3]);
  TString outputs = argv[4];
  Int_t random_seed = argc>5 ? atoi(argv[5]) : 0;//placing the inputs into variables


  TFile *fout = new TFile(outputs.Data(),"RECREATE");
  fout->cd();//opening of the output file

  //---------------------
  //Pythia initialization 
  //---------------------
  Pythia pythia;   // Generator.
  //Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Read in commands from external file.
  pythia.readFile(pythiaconfig);   

  // Extract settings to be used in the main program.
  int    nEvent  = pythia.mode("Main:numberOfEvents");
  bool   showCS  = pythia.flag("Main:showChangedSettings");
  bool   showCPD = pythia.flag("Main:showChangedParticleData");
  double energy  = pythia.mode("Beams:eCM");

  pythia.readString(Form("PhaseSpace:pTHatMin ==%f",pTHatMin));
  pythia.readString(Form("PhaseSpace:pTHatMax ==%f",pTHatMax));
  cout<<"Events="<<nEvent <<" RNDM seed "<< random_seed << endl;

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed=%02d",random_seed));

  // Initialize. Beam parameters set in .cmnd file.
  pythia.init();

  // List changed data. 
  if (showCS)  pythia.settings.listChanged();
  if (showCPD) pdt.listChanged();


  TH1D *hCrossSectionInfo = new TH1D("hCrossSection","CrossSectionInfo",6,0,6);

  //------------------------------------------------------------------
  // Define jet reconstruction
  //------------------------------------------------------------------
  vector<PseudoJet> finalparticles;
  vector<PseudoJet> chparticles;

  double PartMinPtCutForJet = 0.5;// atlas paper
  double etaMaxCutForJet = 1.2;///
  double coneR = 0.6;
  double etaMaxCutForPart = etaMaxCutForJet+coneR;///
  double MinJetPt = 20.; // Min Jet Pt cut to disregard low pt jets

  JetDefinition jet_def(antikt_algorithm, coneR); 

  TClonesArray *inputList = new TClonesArray("AliJBaseTrack",1500);
  //--------------------------------------------------------
  //         B e g i n    e v e n t    l o o p.
  //--------------------------------------------------------
  cout<<"Let's start" <<endl; 
  int ieout = nEvent/20;
  if (ieout<1) ieout=1;
  int EventCounter = 0;

  for(int iEvent = 0; iEvent < nEvent; ++iEvent) {//begin event loop

    if (!pythia.next()) continue;
    finalparticles.clear();
    chparticles.clear();
    if(iEvent % ieout == 0) cout << iEvent << "\t" << int(float(iEvent)/nEvent*100) << "%" << endl ;
    int icent = 0;
    
    for (int i = 0; i < pythia.event.size(); ++i) {//loop over all the particles in the event
      // Building input particle list for the jet reconstruction
      // Make sure "applying same cut done in ATLAS paper ?  flavor of the parciles ,eta , min pt cut 
      if (pythia.event[i].isFinal() && TMath::Abs(pythia.event[i].eta()) < etaMaxCutForPart && pythia.event[i].pT()>PartMinPtCutForJet){
        // Building the charged particle list to be used for manual constituent construction instead of using fastjet consistuents.
        if( pythia.event[i].isCharged() && pythia.event[i].isHadron() ) { // Only check if it is charged and hadron since the acceptance was checked in the previous if
          chparticles.push_back( PseudoJet(pythia.event[i].px(), pythia.event[i].py() , pythia.event[i].pz(), pythia.event[i].e() ) );
        }

        finalparticles.push_back( particle );

        TLorentzVector lvparticle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        AliJBaseTrack track( lvparticle );
        track.SetID(pythia.event[i].id());
        track.SetParticleType(kJHadron);
        track.SetTrackEff(1.);
        new ((*inputList)[inputList->GetEntriesFast()]) AliJBaseTrack(track);

      } // end of finalparticles
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Run the clustering, Reconstruct jets
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ClusterSequence cs(finalparticles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Loop over jets and fill various histos 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for (unsigned i = 0; i < jets.size(); i++) {
      // jet eta cut
      if(fabs(jets[i].eta())>etaMaxCutForJet) continue; 

      vector<PseudoJet> constituents = jets[i].constituents();
      fhistos->fhJetPseudorapidity->Fill(jets[i].eta());  
      fhistos->fhJetRapidity->Fill(jets[i].rap());  
      fhistos->fhJetPt->Fill(jets[i].pt());

      int ptbin = fcard->GetBin(kTriggType, jets[i].pt());
      // Only for the jet pt bin registered in the card
      if(ptbin<0) continue;

      TLorentzVector lvjet = TLorentzVector(jets[i].px(),jets[i].py(),jets[i].pz(),jets[i].e());//creating the jet lorentz vector for the deltar calculation

      //Charged particle loop
      int counter = 0;
      for(unsigned k = 0; k<chparticles.size() ; k++){
        TLorentzVector lvpart = TLorentzVector(chparticles[k].px(),chparticles[k].py(),chparticles[k].pz(),chparticles[k].e());//particle lorentz vector
        double dist=lvjet.DeltaR(lvpart);//delta R distance=sqrt(...)
        if (dist<coneR){
          counter++;
          double dotproductch=dotproduct(lvjet,lvpart); 
          double zch=dotproductch/(jets[i].modp2());
          double xich=TMath::Log(1./zch);
          fhistos->fhZch[ptbin]->Fill(zch);
          fhistos->fhZchLogBin[ptbin]->Fill(zch);
          fhistos->fhXich[ptbin]->Fill(xich);
        }
      }
      fhistos->fhJetNch[ptbin]->Fill(counter); 
      fhistos->fhTriggPtBin[ptbin]->Fill(jets[i].pt());
      fhistos->fhJetPseudorap[ptbin]->Fill(jets[i].eta());
      fhistos->fhJetrap[ptbin]->Fill(jets[i].rap());

      //Consistuents loop
      int Np=0,Nn=0;
      for (unsigned j = 0; j < constituents.size(); j++) {
        TLorentzVector lvfinalpart = TLorentzVector(constituents[j].px(),constituents[j].py(),constituents[j].pz(),constituents[j].e()); //
        double dotproductall=dotproduct(lvjet,lvfinalpart); 
        double z=dotproductall/(jets[i].modp2());
        double xi=TMath::Log(1./z);
        fhistos->fhZ[ptbin]->Fill(z);
        fhistos->fhZLogBin[ptbin]->Fill(z);
        fhistos->fhXi[ptbin]->Fill(xi);
      } // end of constituents loop
      fhistos->fhJetN[ptbin]->Fill(constituents.size());// number of all particle per jet histo
      double nr = double(Nn)/double(Np);
      fhistos->fhChargeRatioInJet[ptbin]->Fill(nr);
      if(DEBUG) printf("Charge info in a jet ::: pT=%.1f Np=%d Nn=%d %.1f\n",jets[i].pt(),Np,Nn,nr);

    }//end of the jet loop
    EventCounter++;
  }//event loop

  long nTried = pythia.info.nTried();
  long nAccepted = pythia.info.nAccepted();
  double sigmaGen = pythia.info.sigmaGen();
  //double sigmaErr = pythia.info.sigmaErr();
  hCrossSectionInfo->Fill(0.5,nTried);
  hCrossSectionInfo->Fill(1.5,nAccepted);
  hCrossSectionInfo->Fill(2.5,sigmaGen);
  hCrossSectionInfo->Fill(3.5,EventCounter);
  hCrossSectionInfo->Fill(4.5,energy);
  hCrossSectionInfo->Fill(5.5,1);

  fout->Write();
  fcard->WriteCard(fout); // Write card into the file
  fout->Close();
  cout << EventCounter << " events are analyzed successfully."<< endl;
  timer.Print(); 
  return 0;
}


// CHECK FUNCTION
double dotproduct(TLorentzVector jets,TLorentzVector particles){//takes in two Lorentz vectors and transform it into TVector3 and makes the dot product
  double jetpx=jets.Px();
  double jetpy=jets.Py();
  double jetpz=jets.Pz();

  double partpx=particles.Px();
  double partpy=particles.Py();
  double partpz=particles.Pz();

  TVector3 jet=TVector3(jetpx,jetpy,jetpz);
  TVector3 particle=TVector3(partpx,partpy,partpz);

  double dotproduct= particle.Dot(jet);

  //delete jet;
  //delete particle;

  return dotproduct;
}
