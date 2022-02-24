//===========================================================
// JHistos.h
//===========================================================

#ifndef JHISTOS_H
#define JHISTOS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TFile.h>

#include  "AliJCard.h"
#include <TClonesArray.h>
#include <TFormula.h>
#define kPtDim 15
#define kPTypeDim 15

class JHistos {

	public:
		JHistos(AliJCard* cardP); 
		virtual ~JHistos(){;}	  //destructor

		// ALICE methods =====================================================
		void CreateQAHistos();
		void CreateFFHistos();
		void CreateDiJetHistos();

	public:
		AliJCard  *fcard;
		TH1D *fhiCentr; // Event Counting
		TH1D *fhJetPt;
		TH1D *fhJetPtATLAS;
		TH1D *fhJetPtATLASEta;
		TH1D *fhJetPtATLASInvar;
		TH1D *fhJetPtATLASEtaInvar;
		TH1D *fhJetPseudorapidity;
		TH1D *fhJetRapidity;
		TH1D *fhchPartPt;
		TH1D *fhchPartEta;
		TH1D *fhchPartPhi;
		TH1D *fhJetPartEta;
		TH1D *fhJetPartPt;
		TH1D *fhJetPartPhi;

		TH1D *fhJetPseudorap[kPtDim];
		TH1D *fhJetrap[kPtDim];
		TH1D *fhTriggPtBin[kPtDim]; // counting the number triggers in bins
		TH1D *fhJetNch[kPtDim];
		TH1D *fhJetN[kPtDim];
		TH1D *fhChargeRatioInJet[kPtDim];
		TH1D *fhZch[kPtDim];
		TH1D *fhZchCon[kPtDim];
		TH1D *fhXich[kPtDim];
		TH1D *fhZ[kPtDim];
		TH1D *fhZPType[kPTypeDim][kPtDim]; // particle type
		TH1D *fhZLogBin[kPtDim];
		TH1D *fhZchLogBin[kPtDim];
		TH1D *fhXi[kPtDim];

		TH1D *fhDiJetInvM;
		TH1D *fhDiJetInvMBins[kPtDim];
		TH1D *fhDiJetPtPair[kPtDim];
		TH1D *fhDiJetPtPairRaw[kPtDim];
		TH1D *fhDiJetDeltaPhi[kPtDim];

	protected:

};

#endif

