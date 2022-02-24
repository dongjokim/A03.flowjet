#include <TVector3.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TFormula.h>


#include  "JHistos.h"
#include  "inputbinning.h"


//______________________________________________________________________________
JHistos::JHistos(AliJCard* cardP)   // constructor
{   // constructor
    fcard=cardP;
}

//__________________________________________________________________________________________________
void JHistos::CreateQAHistos() {
	// Log pT bins
        int NBINS=500;
        double LogBinsX[NBINS+1], LimL=0.1, LimH=500;
        double logBW = (log(LimH)-log(LimL))/NBINS;
        for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);

	int Nbins = 100;
	double rapMin=-2; double rapMax=2;

	fhJetPt=new TH1D("fhJetPt","Jet p_{T} spectra",NBINS, LogBinsX);
	fhJetPt->Sumw2();

	fhJetPtATLAS=new TH1D("fhJetPtATLAS","Jet p_{T} spectra 0<y<0.3",numAtlasptbin,atlasptbins );
	fhJetPtATLAS->Sumw2();
	fhJetPtATLASEta=new TH1D("fhJetPtATLASEta","Jet p_{T} spectra 0<#eta<0.3",numAtlasptbin,atlasptbins );
	fhJetPtATLASEta->Sumw2();
	fhJetPtATLASInvar=new TH1D("fhJetPtATLASInvar","Jet p_{T} 1/pt spectra 0<y<0.3",numAtlasptbin,atlasptbins );
	fhJetPtATLASInvar->Sumw2();
	fhJetPtATLASEtaInvar=new TH1D("fhJetPtATLASEtaInvar","Jet p_{T} 1/pt spectra 0<#eta<0.3",numAtlasptbin,atlasptbins );
	fhJetPtATLASEtaInvar->Sumw2();

        fhJetRapidity=new TH1D("fhJetRapidity","Jet Rapidity",Nbins,rapMin,rapMax);
	fhJetRapidity->Sumw2();

        fhJetPseudorapidity=new TH1D("fhJetPseudorapidity","Jet Pseudorapidity",Nbins,rapMin,rapMax);
	fhJetPseudorapidity->Sumw2();

	fhchPartPt=new TH1D("fhPartPt","Charged particles p_{T}",NBINS, LogBinsX);
	fhchPartPt->Sumw2();

	fhchPartEta=new TH1D("fhPartEta","Charged particles rapidity",Nbins,rapMin,rapMax);
	fhchPartEta->Sumw2();

	fhJetPartPt=new TH1D("fhJetPartPt","Jet particles p_{T}",NBINS, LogBinsX);
	fhJetPartPt->Sumw2();

	fhJetPartEta=new TH1D("fhJetPartEta","Jet particles rapidity",Nbins,rapMin,rapMax);
	fhJetPartEta->Sumw2();

	fhJetPartPhi=new TH1D("fhJetPartPhi","Jet particles #phi",Nbins,rapMin,rapMax);
	fhJetPartPhi->Sumw2();
	fhchPartPhi=new TH1D("fhchPartPhi","Jet Charged particles #phi",Nbins,rapMin,rapMax);
	fhchPartPhi->Sumw2();
	
	fhiCentr = new TH1D("hiCentr","ith centrality", 10, -0.5, 9.5);
	fhiCentr->Sumw2();
}

//__________________________________________________________________________________________________
void JHistos::CreateFFHistos() {
  const int bins = NBins-1;
	//double lr = 0.;
	//double ur = 1.0;
	
	int nrapbins = 500;
	int nxibins = 100;
	double lxi = 0.;
	double hxi = 25.;
	
	double ptbw=10/100.0;  //see hPt histo below, let's make 10 bins per 1GeV/c

	int nZ = 100;
	double zLow = 0.01, zHigh = 1.1;

	double logBinsZ[101];
	double logZ = (log(zHigh)-log(zLow))/nZ;
	for(int ij=0;ij<=nZ;ij++) logBinsZ[ij]=zLow*exp(ij*logZ);

	for(int hit=0; hit < fcard->GetNoOfBins(kTriggType); hit++){
		float pTt1 = fcard->GetBinBorder(kTriggType, hit);
		float pTt2 = fcard->GetBinBorder(kTriggType, hit + 1);
		char htit[100], hname[100];

		fhTriggPtBin[hit] = new TH1D(Form("hTriggPtBin%02d",hit),"", (int)TMath::Ceil((pTt2-pTt1)/ptbw),pTt1, pTt2);

		sprintf(htit, "pTt: %3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hZ%02d", hit);
		fhZ[hit] = new TH1D(hname, htit, bins, xbins);
		fhZ[hit]->Sumw2();

		sprintf(hname, "hZLogBin%02d", hit);
		fhZLogBin[hit] = new TH1D(hname, htit, nZ, logBinsZ);
		fhZLogBin[hit]->Sumw2();

		for(int ity=0;ity<kPTypeDim;ity++) {
			sprintf(hname, "hZPType%02d%02d", ity, hit);
			fhZPType[ity][hit] = new TH1D(hname, htit, nZ, logBinsZ);
			fhZPType[ity][hit]->Sumw2();
		}

		sprintf(hname, "hZchLogBin%02d", hit);
		fhZchLogBin[hit] = new TH1D(hname, htit, nZ, logBinsZ);
		fhZchLogBin[hit]->Sumw2();

		sprintf(htit, "PseudoRap: %3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hPseudoRap%02d", hit);
		fhJetPseudorap[hit]=new TH1D(hname, htit, nrapbins, -5, 5);
		fhJetPseudorap[hit]->Sumw2();

		sprintf(htit, "Rap: %3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hRap%02d", hit);
		fhJetrap[hit]=new TH1D(hname, htit, nrapbins, -5, 5);
		fhJetrap[hit]->Sumw2();

		sprintf(htit, "#xi: %3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hXi%02d", hit);
		fhXi[hit]=new TH1D(hname, htit, nxibins, lxi, hxi);
		fhXi[hit]->Sumw2();

		sprintf(htit, "pTt: %3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hZch%02d", hit);
		fhZch[hit] = new TH1D(hname, htit, bins, xbins);
		fhZch[hit]->Sumw2();
		sprintf(hname, "hZchCon%02d", hit);
		fhZchCon[hit] = new TH1D(hname, htit, bins, xbins);
		fhZchCon[hit]->Sumw2();

		sprintf(htit, "#xi: %3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hXich%02d", hit);
		fhXich[hit]=new TH1D(hname, htit, nxibins, lxi, hxi);
		fhXich[hit]->Sumw2();

		sprintf(htit, "Ncht: %3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hNch%02d", hit);
		fhJetNch[hit]=new TH1D(hname, htit, 101, -0.5 , 100.5);
		fhJetNch[hit]->Sumw2();

		sprintf(htit, "Nt: %3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hN%02d", hit);
		fhJetN[hit]=new TH1D(hname, htit, 101, -0.5 , 100.5);
		fhJetN[hit]->Sumw2();

		sprintf(htit, "%3.1f-%3.1f", pTt1, pTt2);
		sprintf(hname, "hChargeRatioInJet%02d", hit);
		fhChargeRatioInJet[hit]=new TH1D(hname, htit, 88 , -0.1 , 2.1);
		fhChargeRatioInJet[hit]->Sumw2();
	}
}
//__________________________________________________________________________________________________
void JHistos::CreateDiJetHistos() {
	int nBINS2=100;
	double logBinsX2[nBINS2+2], limL2=0.1, LimH2=2000;
	double logBW2 = (log(LimH2)-log(limL2))/nBINS2;
	for(int ij=0;ij<=nBINS2;ij++) logBinsX2[ij+1]=limL2*exp(ij*logBW2);
	logBinsX2[0]=0;

	int nDPhi=1000;double xDPhi0=-1; double xDPhi1=TMath::TwoPi()+1;
	double mjjbw=10/100.0;  //see hPt histo below, let's make 10 bins per 1GeV/c

	fhDiJetInvM = new TH1D("hJetInvM","",nBINS2, logBinsX2 );

	for(int him=0; him < fcard->GetNoOfBins(kDiJetType); him++){
		float m1 = fcard->GetBinBorder(kDiJetType, him);
		float m2 = fcard->GetBinBorder(kDiJetType, him + 1);
		char htit[100], hname[100];
		//== DiJets
		sprintf(htit, "%3.1f-%3.1f", m1, m2);
		sprintf(hname, "hDiJetInvMBins%02d", him);
		fhDiJetInvMBins[him] = new TH1D(hname,htit,(int)TMath::Ceil((m2-m1)/mjjbw), m1, m2 );
		fhDiJetInvMBins[him]->Sumw2();
		sprintf(htit, "%3.1f-%3.1f", m1, m2);
		sprintf(hname, "hDiJetPtPair%02d", him);
		fhDiJetPtPair[him] = new TH1D(hname,htit,nBINS2, logBinsX2 );
		fhDiJetPtPair[him]->Sumw2();
		sprintf(hname, "hDiJetPtPairRaw%02d", him);
		fhDiJetPtPairRaw[him] = new TH1D(hname,htit,nBINS2, logBinsX2 );
		fhDiJetPtPairRaw[him]->Sumw2();
		sprintf(hname, "hDiJetDeltaPhi%02d", him);
		fhDiJetDeltaPhi[him] = new TH1D(hname,htit,nDPhi, xDPhi0, xDPhi1);
		fhDiJetDeltaPhi[him]->Sumw2();
	}
}
