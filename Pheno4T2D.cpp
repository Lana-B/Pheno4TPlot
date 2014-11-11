//============================================================================
// Name        : HelloROOT2.cpp
// Author      : Lana
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"

#include "TStyle.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TLegend.h"
#include "TFunction.h"
#include "TF1.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//#include "../macros/Style.C"
using namespace std;

void Pheno4T2D(){

	const int nPoints = 12;
	const int nPointscoup = 5;

	double dSGMass[nPoints] = {350, 375, 400, 425, 450, 475, 500, 600, 700, 800, 900, 1000}; //sgluon masses
	double dSRefficiency[nPoints] = {0.003758, 0.0040968, 0.004889, 0.0058141, 0.00689055, 0.0074328, 0.00839955, 0.0100639, 0.0127462, 0.0152723, 0.0175777, 0.0173885}; //efficiencies in SR28 for nominal 0.03 luL3x3 lUR3x3
	double dSGKfactor[nPoints] = {1.7, 1.7, 1.7, 1.7, 1.75, 1.8, 1.8, 1.8, 1.9, 2.0, 2.1, 2.2}; //K factor for NLO
	double dSG_BR_tt[nPoints] = {0.0298, 0.2121, 0.3326, 0.4022, 0.4426, 0.4658, 0.4784, 0.4783, 0.4472, 0.4083, 0.3689, 0.3319}; //slguon branching ratios to ttbar
	double dSGxsec2[nPoints] = {5.411, 177.18, 305.9, 308.23, 254.8, 194.9, 147.99, 39.05, 10.77, 3.09, 0.922, 0.284}; // sgluon cross section in fb

	double coupling1[nPoints] = {0.0130262, 0.0187498, 0.0248832, 0.0299049. 0.0335501, 0.0359585, 0.0373754, 0.0373562, 0.0340081, 0.0304148, 0.0273303, 0.0248388}; //weights for 0.01 luL3x3 lUR3x3
	double coupling2[nPoints] = {0.204233, 0.253832, 0.297214, 0.327543, 0.347328, 0.359539, 0.366436, 0.366344, 0.3497, 0.330413. 0.312486, 0.296728}; //weights for 0.02 luL3x3 lUR3x3
	double coupling4[nPoints] = {3.01901, 2.3287, 1.99494, 1.83384, 1.74911, 1.70298, 1.67871,  1.67903, 1.73982, 1.82069, 1.90835, 1.99669}; //weights for 0.04 luL3x3 lUR3x3
	double coupling5[nPoints] = {6.95959, 4.06873, 3.04732, 2.62383, 2.41687, 2.30882, 2.2532, 2.25391, 2.39485, 2.59075, 2.81438, 3.05215}; //weights in SR28 for 0.05 luL3x3 lUR3x3

	double couplingPoints[nPointscoup] = {0.001, 0.002, 0.003, 0.004, 0.005}; //luL3x3 lUR3x3

	double dLumi = 19.5; //fb^-1
	double SM4T_SReff = 0.0078; //Search region efficiency for Standard Model 4 top
	double SM4T_xsec = 1.3; //fb
	double SM4T_dilepBR_SS = 0.2962*0.5; //dilepton branching ratio * same sign

	double SM_nEvents = SM4T_SReff*SM4T_xsec*dLumi;  //predicted number of SM4T events 

	cout<<"SM number of events : "<<SM_nEvents<<endl;

	double dSG_nEvents[nPoints * nPointscoup];
	double dmassfull[nPoints * nPointscoup];
	double dcoupfull[nPoints * nPointscoup];

	for(int i=0; i<nPoints; i++){ //Fill sgluon mass array
		dmassfull[i] = dSGMass[i];
		dmassfull[i + nPoints] = dSGMass[i];
		dmassfull[i + 2*nPoints] = dSGMass[i];
		dmassfull[i + 3*nPoints] = dSGMass[i];
		dmassfull[i + 4*nPoints] = dSGMass[i];
	}

	for(int i=0; i<nPoints; i++){ //Fill coupling points array
		dcoupfull[i] = couplingPoints[0];
		dcoupfull[i + nPoints] = couplingPoints[1];
		dcoupfull[i + 2*nPoints] = couplingPoints[2];
		dcoupfull[i + 3*nPoints] = couplingPoints[3];
		dcoupfull[i + 4*nPoints] = couplingPoints[4];

	}

	for (int i=0; i<nPoints; i++){  //Fill nEvents for every combination of coupling vs mass
		double nEvents1 = dSRefficiency[i]*dSGxsec2[i]*dSGKfactor[i]*dLumi*coupling1[i];   //weighted by coupling 1
		dSG_nEvents[i] = nEvents1;
		double nEvents2 = dSRefficiency[i]*dSGxsec2[i]*dSGKfactor[i]*dLumi*coupling2[i];   //weighted by coupling 2
		dSG_nEvents[i + nPoints] = nEvents2;
		double nEvents3 = dSRefficiency[i]*dSGxsec2[i]*dSGKfactor[i]*dLumi;				   //nominal
		dSG_nEvents[i + 2*nPoints] = nEvents3;
		double nEvents4 = dSRefficiency[i]*dSGxsec2[i]*dSGKfactor[i]*dLumi*coupling4[i];   //weighted by coupling 4
		dSG_nEvents[i + 3*nPoints] = nEvents4;
		double nEvents5 = dSRefficiency[i]*dSGxsec2[i]*dSGKfactor[i]*dLumi*coupling5[i];   //weighted by coupling 5
		dSG_nEvents[i + 4*nPoints] = nEvents5;
	}

	TFile *f1 = new TFile("Pheno4T2D.root", "RECREATE");
	TCanvas *c2 = new TCanvas("c1","nevents", 600, 600);
	c2->SetLogz();

	TGraph2D *gnEvents = new TGraph2D();
	for(i =0; i<(nPoints*nPointscoup); i++){
		gnEvents->SetPoint(i, dmassfull[i], dcoupfull[i], dSG_nEvents[i]);
	}

	gnEvents->Draw("COLZ");

	TH2D* h_theory = gnEvents->GetHistogram();
    c2->SetRightMargin(0.15);
    c2->SetLeftMargin(0.11);
    h_theory->SetLineWidth(4);
    h_theory->SetLineStyle(5);
    h_theory->SetLineColor(kBlack);
    h_theory->GetYaxis()->SetTitleOffset(1.48);
    h_theory->SetTitle("");
    h_theory->SetXTitle("M_{sgluon} (GeV)");
    h_theory->SetYTitle("Effective sgluon-top quark coupling");
    h_theory->SetZTitle("events");
    h_theory->DrawClone("colz");
    double contours[2];
    contours[0] = 5.6.;   //contour at 5.6 events
    contours[1] = 0.0013;
    h_theory->SetContour(1,contours);
    h_theory->DrawClone("cont1same");
	//gnEvents->SetTitle("");
	//gnEvents->GetYaxis()->SetTitle("coupling");
	//gnEvents->GetXaxis()->SetTitle("sgluon mass");

	//TLine *line = new TLine(300,0.2,1065,0.2);
	//line->SetLineColor(kRed);
	//line->Draw();

    f1->cd();
	h_theory->Write("h_theory",TObject::kOverwrite);

	c2->SaveAs("SG_nEvents2D.pdf");
	delete c2;

}
