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

void Pheno4T(){

	const int nPoints = 12;

	double dSGMass[nPoints] = {350, 375, 400, 425, 450, 475, 500, 600, 700, 800, 900, 1000};
	double dSRefficiency[nPoints] = {0.003758, 0.0040968, 0.004889, 0.0058141, 0.00689055, 0.0074328, 0.00839955, 0.0100639, 0.0127462, 0.0152723, 0.0175777, 0.0173885};
	double dSGxsec[nPoints] = {0, 0, 1625, 0, 0, 0, 358.1, 94.9, 28.4, 9.26, 3.22, 1.17}; //cross section in fb
	double dSGKfactor[nPoints] = {1.7, 1.7, 1.7, 1.7, 1.75, 1.8, 1.8, 1.8, 1.9, 2.0, 2.1, 2.2};
	double dSG_BR_tt[nPoints] = {0.0298, 0.2121, 0.3326, 0.4022, 0.4426, 0.4658, 0.4784, 0.4783, 0.4472, 0.4083, 0.3689, 0.3319};
	double dSGxsec2[nPoints] = {5.411, 177.18, 305.9, 308.23, 254.8, 194.9, 147.99, 39.05, 10.77, 3.09, 0.922, 0.284}; //cross section in fb

	double dLumi = 19.5; //fb^-1
	double SM4T_SReff = 0.0078;
	double SM4T_xsec = 1.3; //fb
	double SM4T_dilepBR_SS = 0.2962*0.5; //dilepton branching ration * same sign

	double SM_nEvents = SM4T_SReff*SM4T_xsec*dLumi;

	cout<<"SM number of events : "<<SM_nEvents<<endl;

	double dSG_nEvents[nPoints];

	for (int i=0; i<nPoints; i++){
		double nEvents = dSRefficiency[i]*dSGxsec2[i]*dSGKfactor[i]*dLumi;
		dSG_nEvents[i] = nEvents;
		cout<<"ev: "<<nEvents<<endl;
	}

	TCanvas *c1 = new TCanvas("c1","test", 600, 600);

	TGraph *gEfficiency = new TGraph(nPoints, dSGMass, dSRefficiency);
	gEfficiency->Draw("a3l");
	gEfficiency->SetTitle("");

		TLine *line2 = new TLine(300,0.0078,1065,0.0078);
	line2->SetLineColor(kRed);
	line2->Draw();

	c1->SaveAs("SG_eff.pdf");
	delete c1;

	TCanvas *c2 = new TCanvas("c1","nevents", 600, 600);
	TGraph *gnEvents = new TGraph(nPoints, dSGMass, dSG_nEvents);
	gnEvents->Draw("a3l");
	gnEvents->SetTitle("");
	gnEvents->GetYaxis()->SetTitle("number of events");
	gnEvents->GetXaxis()->SetTitle("sgluon mass");

	TLine *line = new TLine(300,0.2,1065,0.2);
	line->SetLineColor(kRed);
	line->Draw();

	c2->SaveAs("SG_nEvents.pdf");
	delete c2;

	//return 0;
}
