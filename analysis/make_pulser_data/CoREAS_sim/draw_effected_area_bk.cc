#include "TTree.h"
#include <iostream>
#include "TGraph.h" 
#include "TVirtualFFT.h" 
#include "TFile.h"
#include <fstream>
#include <math.h>
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include <iostream>
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TChain.h"
#include <TBits.h>
#include <ctime>

#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>

using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8

#define HALF_WINDOW_SPACING 75



//TRandom *r3;

//.x make_coreas_data.cc(1, "", 20210408, 18969, 18970)
void draw_effected_area(){


	TFile *f_ntuple3 = new TFile("ntuple6.root");
	TNtuple *ntuple3 = (TNtuple *) f_ntuple3->Get("ntuple");
	ntuple3->SetMarkerStyle(8);
	//ntuple3->SetMarkerColor(3);
	//ntuple3->SetLineColor(3);
	//*0.26195951349
	//ntuple3->Draw("EffectedArea*1.0e-6:ShowerTheta", "", "P");

	

	
	
	TFile *f_ntuple2 = new TFile("ntuple5.root");
	TNtuple *ntuple2 = (TNtuple *) f_ntuple2->Get("ntuple");
	ntuple2->SetMarkerStyle(8);
	ntuple2->SetMarkerColor(2);
	ntuple2->SetLineColor(2);
	//*0.26195951349
	//ntuple2->Draw("EffectedArea*1.0e-6:ShowerTheta", "", "P same");


	TFile *f_ntuple1 = new TFile("ntuple4.root");
	TNtuple *ntuple1 = (TNtuple *) f_ntuple1->Get("ntuple");
	ntuple1->SetMarkerStyle(8);
	ntuple1->SetMarkerColor(3);
	ntuple1->SetLineColor(3);
	//*0.26195951349
	//"ShowerTheta>14.8&&ShowerTheta<15.2"
	//ntuple1->Draw("EffectedArea*1.0e-6:ShowerTheta", "", "P same");
		
	
	
	float _Energy, _ShowerTheta, _ShowerPhi, _EffectedArea;
	double Energy[3][9][19], ShowerTheta[3][9][19], ShowerPhi[3][9][19], EffectedArea_tp[3][9][19], EffectedArea_pt[3][19][9], theta_list[10], phi_list[19];;
	double Sum_phi[3][10] = {0.};
	double Sum_theta[3][19] = {0.};
	double Sum_phitheta[3] = {0.};
	
	ntuple1->SetBranchAddress("Energy", &_Energy);
	ntuple1->SetBranchAddress("ShowerTheta", &_ShowerTheta);
	ntuple1->SetBranchAddress("ShowerPhi", &_ShowerPhi);
	ntuple1->SetBranchAddress("EffectedArea", &_EffectedArea);
	
	ntuple2->SetBranchAddress("Energy", &_Energy);
	ntuple2->SetBranchAddress("ShowerTheta", &_ShowerTheta);
	ntuple2->SetBranchAddress("ShowerPhi", &_ShowerPhi);
	ntuple2->SetBranchAddress("EffectedArea", &_EffectedArea);
	
	ntuple3->SetBranchAddress("Energy", &_Energy);
	ntuple3->SetBranchAddress("ShowerTheta", &_ShowerTheta);
	ntuple3->SetBranchAddress("ShowerPhi", &_ShowerPhi);
	ntuple3->SetBranchAddress("EffectedArea", &_EffectedArea);
	
	for(int i =0; i<ntuple1->GetEntries(); i++){
		int i_phi = i%19;
		int i_theta = i/19;
		
		
	
		ntuple1->GetEntry(i);
		Energy[0][i_theta][i_phi] = _Energy;
		//cout<<_Energy<<endl;
		ShowerTheta[0][i_theta][i_phi] = 90 - _ShowerTheta;
		ShowerPhi[0][i_theta][i_phi] = _ShowerPhi;
		EffectedArea_tp[0][i_theta][i_phi] = _EffectedArea;
		EffectedArea_pt[0][i_phi][i_theta] = _EffectedArea;
		//Sum_phi[0][i_theta] += _EffectedArea;
		//Sum_theta[0][i_phi] += _EffectedArea*TMath::Cos(_ShowerTheta);
		//Sum_phitheta[0] += _EffectedArea;
		//if(i_theta==2){cout<<"  "<<i_phi<<"  "<<EffectedArea_tp[0][i_theta][i_phi]<<endl;}
		
		theta_list[i_theta] = 90 - _ShowerTheta;
		phi_list[i_phi] = _ShowerPhi;
		
		ntuple2->GetEntry(i);
		Energy[1][i_theta][i_phi] = _Energy;
		ShowerTheta[1][i_theta][i_phi] = 90 - _ShowerTheta;
		ShowerPhi[1][i_theta][i_phi] = _ShowerPhi;
		EffectedArea_tp[1][i_theta][i_phi] = _EffectedArea;
		EffectedArea_pt[1][i_phi][i_theta] = _EffectedArea;
		//Sum_phi[1][i_theta] += _EffectedArea;
		//Sum_theta[1][i_phi] += _EffectedArea*TMath::Cos(_ShowerTheta);
		//Sum_phitheta[1] += _EffectedArea;
		
		ntuple3->GetEntry(i);
		Energy[2][i_theta][i_phi] = _Energy;
		ShowerTheta[2][i_theta][i_phi] = 90 - _ShowerTheta;
		ShowerPhi[2][i_theta][i_phi] = _ShowerPhi;
		EffectedArea_tp[2][i_theta][i_phi] = _EffectedArea;
		EffectedArea_pt[2][i_phi][i_theta] = _EffectedArea;
		//Sum_phi[2][i_theta] += _EffectedArea;
		//Sum_theta[2][i_phi] += _EffectedArea*TMath::Cos(_ShowerTheta);
		//Sum_phitheta[2] += _EffectedArea;
		
		//cout<<i_phi<<"  "<<i_theta<<"   "<<_ShowerPhi<<"   "<<_ShowerTheta<<endl;
	
	}
	
	//TGraph *wyn = new TGraph(9, theta_list, EffectedArea_pt[0][0]);
	//wyn->Draw();
	
	for(int e=0; e<3; e++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			TGraph *wyn = new TGraph(9, theta_list, EffectedArea_pt[e][i_phi]);
			for(double theta = 87.9; theta>=40.1; theta-=0.2){
				Sum_theta[e][i_phi] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*wyn->Eval(theta)/(TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((88)*TMath::DegToRad()));
			}
			delete wyn;
			wyn = NULL;
		}
		

		for(int i_theta = 0; i_theta<9; i_theta++){
			for(int i_phi = 0; i_phi<19; i_phi++){
				Sum_phi[e][i_theta] += EffectedArea_tp[e][i_theta][i_phi]/19.;
				//if(i_theta==2){cout<<e<<"  "<<i_phi<<"  "<<EffectedArea_tp[e][i_theta][i_phi]/19.<<endl;}
			}
			//cout<<Sum_phi[e][2]<<endl;
		}
		
	
	}
	
	
	/*
	double a[2] = {-70, 30};
	double b[2] = {1.0e5, 3.0e7};
	TGraph *wyn = new TGraph(2, a, b);
	wyn->Draw("AP");
	wyn->SetMarkerColor(0);
	
	TGraph *cyc = new TGraph(19, phi_list, Sum_theta[2]);
	cyc->SetMarkerStyle(8);
	cyc->Draw("same");
	//cyc->SetMarkerColor(2);
	//cyc->SetLineColor(2);
	
	TGraph *cyc1 = new TGraph(19, phi_list, Sum_theta[1]);
	
	cyc1->Draw("same");
	cyc1->SetMarkerStyle(8);
	cyc1->SetMarkerColor(2);
	cyc1->SetLineColor(2);
	
	TGraph *cyc2 = new TGraph(19, phi_list, Sum_theta[0]);
	cyc2->SetMarkerStyle(8);
	cyc2->Draw("same");
	cyc2->SetMarkerColor(3);
	cyc2->SetLineColor(3);
	
	*/
	
	
	
	double a[2] = {-70, 30};
	double b[2] = {1.0e5, 3.0e7};
	TGraph *wyn = new TGraph(2, a, b);
	//wyn->Draw("AP");
	wyn->SetMarkerColor(0);
	
	
	TGraph *cyc[3];
	
	theta_list[9] = 89;
	Sum_phi[0][9] = 0;
	Sum_phi[1][9] = 0;
	Sum_phi[2][9] = 0;
	
	cyc[2] = new TGraph(10, theta_list, Sum_phi[2]);
	//cyc->SetMarkerStyle(8);
	cyc[2]->Draw("");
	//cyc->SetMarkerColor(2);
	//cyc->SetLineColor(2);
	
	cyc[1] = new TGraph(10, theta_list, Sum_phi[1]);
	
	
	cyc[1]->SetMarkerStyle(8);
	cyc[1]->SetMarkerColor(2);
	cyc[1]->SetLineColor(2);
	cyc[1]->Draw("same");
	
	cyc[0] = new TGraph(10, theta_list, Sum_phi[0]);
	cyc[0]->SetMarkerStyle(8);
	cyc[0]->SetMarkerColor(3);
	cyc[0]->SetLineColor(3);
	cyc[0]->Draw("same");
	
	
	int e = 2;
	
	double sh_theta[5] = {40.9, 54.1, 65.25, 75.45, 85.2};
	double sh_acceptance[3][5] = {0.};
	/*
	for(int e=0; e<3; e++){
		for(double theta = 48.2; theta>=33.6; theta-=0.2){
			sh_acceptance[e][0] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta)*(3.14159/2.);
			}
		for(double theta = 60.0; theta>=48.2; theta-=0.2){
			sh_acceptance[e][1] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta)*(3.14159/2.);
			}
		for(double theta = 70.5; theta>=60.0; theta-=0.2){
			sh_acceptance[e][2] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta)*(3.14159/2.);
			}
			for(double theta = 80.4; theta>=70.5; theta-=0.2){
			sh_acceptance[e][3] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta)*(3.14159/2.);
			}
			for(double theta = 90.0; theta>=80.4; theta-=0.2){
			sh_acceptance[e][4] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta)*(3.14159/2.);
			}
			
		}
	*/		
	
	
	sh_acceptance[e][0] = Sum_phi[e][0]*3.14159/2.*(TMath::Cos(33.6*TMath::DegToRad())-TMath::Cos(48.2*TMath::DegToRad()));
	
	sh_acceptance[e][1] = (Sum_phi[e][1]*TMath::Sin(50.*TMath::DegToRad()) + Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()))/(TMath::Sin(50.*TMath::DegToRad()) + TMath::Sin(60.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())));

	//cout<<(Sum_phi[e][1]+Sum_phi[e][2])/2.*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())))<<endl;
	//cout<<sh_acceptance[1]<<endl;
	
	sh_acceptance[e][2] = (Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()) + Sum_phi[e][3]*TMath::Sin(65.*TMath::DegToRad()) + Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()))/(TMath::Sin(60.*TMath::DegToRad()) + TMath::Sin(65.*TMath::DegToRad()) + TMath::Sin(70.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(60.0*TMath::DegToRad())-TMath::Cos(70.5*TMath::DegToRad())));

	sh_acceptance[e][3] = (Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()) + Sum_phi[e][5]*TMath::Sin(75.*TMath::DegToRad()) + Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()))/(TMath::Sin(70.*TMath::DegToRad()) + TMath::Sin(75.*TMath::DegToRad()) + TMath::Sin(80.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(70.7*TMath::DegToRad())-TMath::Cos(80.4*TMath::DegToRad())));
	
	

	sh_acceptance[e][4] = (Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()) + Sum_phi[e][7]*TMath::Sin(85.*TMath::DegToRad()) + Sum_phi[e][8]*TMath::Sin(88.*TMath::DegToRad()))/(TMath::Sin(80.*TMath::DegToRad()) + TMath::Sin(85.*TMath::DegToRad()) + TMath::Sin(88.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(80.4*TMath::DegToRad())-TMath::Cos(90.0*TMath::DegToRad())));
	

	TGraph *hxj = new TGraph(5, sh_theta, sh_acceptance[e]);
	hxj->SetMarkerStyle(8);
	//hxj->Draw("APL");
	hxj->SetMarkerColor(3);
	hxj->SetLineColor(3);	


			
}

