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


	TNtuple *ntuple[7] = {NULL};
	
	TFile *f_ntuple_1750 = new TFile("ntuple_with_pps_1750.root");
	ntuple[0] = (TNtuple *) f_ntuple_1750->Get("ntuple");

	TFile *f_ntuple_1775 = new TFile("ntuple_with_pps_1775.root");
	ntuple[1] = (TNtuple *) f_ntuple_1775->Get("ntuple");
	
	TFile *f_ntuple_1800 = new TFile("ntuple_with_pps_1800.root");
	ntuple[2] = (TNtuple *) f_ntuple_1800->Get("ntuple");
	
	TFile *f_ntuple_1825 = new TFile("ntuple_with_pps_1825.root");
	ntuple[3] = (TNtuple *) f_ntuple_1825->Get("ntuple");
	
	TFile *f_ntuple_1850 = new TFile("ntuple_with_pps_1850.root");
	ntuple[4] = (TNtuple *) f_ntuple_1850->Get("ntuple");
	
	TFile *f_ntuple_1875 = new TFile("ntuple_with_pps_1875.root");
	ntuple[5] = (TNtuple *) f_ntuple_1875->Get("ntuple");
	
	TFile *f_ntuple_real = new TFile("ntuple_with_pps_real_1850.root");
	ntuple[6] = (TNtuple *) f_ntuple_real->Get("ntuple");
	

	double DSR_radius[9] = {750, 750, 710, 782, 870, 1020, 1720, 2550, 3200};
	


	
	
	float _Energy, _ShowerTheta, _ShowerPhi, _EffectedArea;
	double Energy[7][9][19], ShowerTheta[7][9][19], ShowerPhi[7][9][19], EffectedArea_tp[7][9][19], EffectedArea_pt[7][19][9], theta_list[10], phi_list[19];
	double Sum_phi[7][10] = {0.};
	double Sum_theta[7][19] = {0.};
	double Sum_phitheta[7] = {0.};
	
	double Error_EffectedArea_tp[7][9][19], Error_EffectedArea_pt[7][19][9];
	double Error_Sum_phi[7][10] = {0.};
	double Error_Sum_theta[7][19] = {0.};
	double Error_Sum_phitheta[7] = {0.}; 
	
	double Error_x[20] = {0.};
	
	for(int k=0; k<7; k++){
		ntuple[k]->SetBranchAddress("Energy", &_Energy);
		ntuple[k]->SetBranchAddress("ShowerTheta", &_ShowerTheta);
		ntuple[k]->SetBranchAddress("ShowerPhi", &_ShowerPhi);
		ntuple[k]->SetBranchAddress("EffectedArea", &_EffectedArea);
	}
	

	
	for(int i =0; i<ntuple[0]->GetEntries(); i++){
		int i_phi = i%19;
		int i_theta = i/19;
		
		
		for(int k=0; k<6; k++){
			ntuple[k]->GetEntry(i);
			Energy[k][i_theta][i_phi] = _Energy;
			ShowerTheta[k][i_theta][i_phi] = 90 - _ShowerTheta;
			ShowerPhi[k][i_theta][i_phi] = _ShowerPhi;
			
			double _Error_EffectedArea = sqrt((_EffectedArea/(3.14159*DSR_radius[i_theta]*DSR_radius[i_theta]))*600)/600.*(3.14159*DSR_radius[i_theta]*DSR_radius[i_theta]);
			
			EffectedArea_tp[k][i_theta][i_phi] = _EffectedArea*1.0e-6;
			EffectedArea_pt[k][i_phi][i_theta] = _EffectedArea*1.0e-6;
			theta_list[i_theta] = 90 - _ShowerTheta;
			phi_list[i_phi] = _ShowerPhi;
			
			Error_EffectedArea_tp[k][i_theta][i_phi] = _Error_EffectedArea*1.0e-6;
			Error_EffectedArea_pt[k][i_phi][i_theta] = _Error_EffectedArea*1.0e-6;
		}
		
		bool real = true;
		if(real){
			ntuple[6]->GetEntry(i+0*ntuple[0]->GetEntries());
			cout<<_Energy<<endl;
			if(_Energy!=18.25) continue;
			Energy[3][i_theta][i_phi] = _Energy;
			ShowerTheta[3][i_theta][i_phi] = 90 - _ShowerTheta;
			ShowerPhi[3][i_theta][i_phi] = _ShowerPhi;
			
			double _Error_EffectedArea = sqrt((_EffectedArea/(3.14159*DSR_radius[i_theta]*DSR_radius[i_theta]))*600)/600.*(3.14159*DSR_radius[i_theta]*DSR_radius[i_theta]);
			
			EffectedArea_tp[3][i_theta][i_phi] = _EffectedArea*1.0e-6;
			EffectedArea_pt[3][i_phi][i_theta] = _EffectedArea*1.0e-6;
			theta_list[i_theta] = 90 - _ShowerTheta;
			phi_list[i_phi] = _ShowerPhi;
			
			Error_EffectedArea_tp[3][i_theta][i_phi] = _Error_EffectedArea*1.0e-6;
			Error_EffectedArea_pt[3][i_phi][i_theta] = _Error_EffectedArea*1.0e-6;

		
		}
		
		
	
	}
	
	//TGraph *wyn = new TGraph(9, theta_list, EffectedArea_pt[3][0]);
	//wyn->Draw();
	
	for(int k=0; k<7; k++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			TGraph *wyn = new TGraph(9, theta_list, EffectedArea_pt[k][i_phi]);
			for(double theta = 87.9; theta>=40.1; theta-=0.2){
				Sum_theta[k][i_phi] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*wyn->Eval(theta)/(TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((88)*TMath::DegToRad()));
			}
			delete wyn;
			wyn = NULL;
			
			
			
			TGraph *hxj = new TGraph(9, theta_list, Error_EffectedArea_pt[k][i_phi]);
			for(double theta = 87.9; theta>=40.1; theta-=0.2){
				Error_Sum_theta[k][i_phi] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*hxj->Eval(theta)/(TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((88)*TMath::DegToRad()));
			}
			delete hxj;
			hxj = NULL;
		}
		

		for(int i_theta = 0; i_theta<9; i_theta++){
			for(int i_phi = 0; i_phi<19; i_phi++){
				Sum_phi[k][i_theta] += EffectedArea_tp[k][i_theta][i_phi]/19.;
				//if(i_theta==2){cout<<e<<"  "<<i_phi<<"  "<<EffectedArea_tp[e][i_theta][i_phi]/19.<<endl;}
				Error_Sum_phi[k][i_theta] += Error_EffectedArea_tp[k][i_theta][i_phi]/19.;
			}
			//cout<<Sum_phi[e][2]<<endl;
		}
		
	
	}
	
	
	
	TGraphErrors *cyc[7] = {NULL};
	
	/*
	double a[2] = {-70, 30};
	double b[2] = {1.0e-5, 4.0e0};
	TGraph *wyn = new TGraph(2, a, b);
	wyn->Draw("AP");
	wyn->SetMarkerColor(0);
	
	for(int k=0; k<6; k++){
		//cyc[k] = new TGraph(9, theta_list, Sum_phi[k]);
		//cyc[k]->SetMarkerStyle(8);
		//cyc[k]->SetMarkerColor(k+1);
		//cyc[k]->SetLineColor(k+1);
		//cyc[k]->Draw("same");
		
		cyc[k] = new TGraphErrors(19, phi_list, Sum_theta[k], Error_x, Error_Sum_theta[k]);
		cyc[k]->SetMarkerStyle(8);
		cyc[k]->SetMarkerColor(k+1);
		cyc[k]->SetLineColor(k+1);
		cyc[k]->SetFillColor(k+1);
		cyc[k]->SetFillStyle(3005);
		cyc[k]->Draw("PL same 3");
		
	}
	*/
	
	double a[2] = {0, 90};
	double b[2] = {1.0e-5, 20.0e0};
	TGraph *wyn = new TGraph(2, a, b);
	wyn->Draw("AP");
	wyn->SetMarkerColor(0);
	for(int k=0; k<7; k++){
		//Sum_phi[k][6] = 0.;
 		//Sum_phi[k][7] = 0.;
		//Sum_phi[k][8] = 0.;
		//cyc[k] = new TGraph(9, theta_list, Sum_phi[k]);
		//cyc[k]->SetMarkerStyle(8);
		//cyc[k]->SetMarkerColor(k+1);
		//cyc[k]->SetLineColor(k+1);
		//cyc[k]->Draw("same");
		
		cyc[k] = new TGraphErrors(9, theta_list, Sum_phi[k], Error_x, Error_Sum_phi[k]);
		cyc[k]->SetMarkerStyle(8);
		cyc[k]->SetMarkerColor(k+1);
		cyc[k]->SetLineColor(k+1);
		cyc[k]->SetFillColor(k+1);
		cyc[k]->SetFillStyle(3005);
		//cyc[k]->Draw("PL same 3");
		
	}

	
	
	

	
	int e = 5;
	
	double sh_theta[5] = {40.9, 54.1, 65.25, 75.45, 85.2};
	double sh_acceptance[6][5] = {0.};
	double Error_sh_acceptance[6][5] = {0.};
	
	for(int e=0; e<6; e++){
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
		
		TGraph *hxj = new TGraph(5, sh_theta, sh_acceptance[e]);
		hxj->SetMarkerStyle(9);
		//if(e==3)
		//hxj->Draw("PL same");
		hxj->SetMarkerColor(e+1);
		hxj->SetLineColor(e+1);
	}
	
	
	for(int e=0; e<6; e++){
/*	
		sh_acceptance[e][0] = Sum_phi[e][0]*3.14159/2.*(TMath::Cos(33.6*TMath::DegToRad())-TMath::Cos(48.2*TMath::DegToRad()));
	
		sh_acceptance[e][1] = (Sum_phi[e][1]*TMath::Sin(50.*TMath::DegToRad()) + Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()))/(TMath::Sin(50.*TMath::DegToRad()) + TMath::Sin(60.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())));

		//cout<<(Sum_phi[e][1]+Sum_phi[e][2])/2.*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())))<<endl;
		//cout<<sh_acceptance[1]<<endl;
	
		sh_acceptance[e][2] = (Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()) + Sum_phi[e][3]*TMath::Sin(65.*TMath::DegToRad()) + Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()))/(TMath::Sin(60.*TMath::DegToRad()) + TMath::Sin(65.*TMath::DegToRad()) + TMath::Sin(70.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(60.0*TMath::DegToRad())-TMath::Cos(70.5*TMath::DegToRad())));

		sh_acceptance[e][3] = (Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()) + Sum_phi[e][5]*TMath::Sin(75.*TMath::DegToRad()) + Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()))/(TMath::Sin(70.*TMath::DegToRad()) + TMath::Sin(75.*TMath::DegToRad()) + TMath::Sin(80.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(70.7*TMath::DegToRad())-TMath::Cos(80.4*TMath::DegToRad())));


		sh_acceptance[e][4] = (Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()) + Sum_phi[e][7]*TMath::Sin(85.*TMath::DegToRad()) + Sum_phi[e][8]*TMath::Sin(88.*TMath::DegToRad()))/(TMath::Sin(80.*TMath::DegToRad()) + TMath::Sin(85.*TMath::DegToRad()) + TMath::Sin(88.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(80.4*TMath::DegToRad())-TMath::Cos(90.0*TMath::DegToRad())));
	
*/





		Error_sh_acceptance[e][0] = Error_Sum_phi[e][0]*3.14159/2.*(TMath::Cos(33.6*TMath::DegToRad())-TMath::Cos(48.2*TMath::DegToRad()));
	
		Error_sh_acceptance[e][1] = (Error_Sum_phi[e][1]*TMath::Sin(50.*TMath::DegToRad()) + Error_Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()))/(TMath::Sin(50.*TMath::DegToRad()) + TMath::Sin(60.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())));

		//cout<<(Sum_phi[e][1]+Sum_phi[e][2])/2.*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())))<<endl;
		//cout<<sh_acceptance[1]<<endl;
	
		Error_sh_acceptance[e][2] = (Error_Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()) + Error_Sum_phi[e][3]*TMath::Sin(65.*TMath::DegToRad()) + Error_Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()))/(TMath::Sin(60.*TMath::DegToRad()) + TMath::Sin(65.*TMath::DegToRad()) + TMath::Sin(70.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(60.0*TMath::DegToRad())-TMath::Cos(70.5*TMath::DegToRad())));

		Error_sh_acceptance[e][3] = (Error_Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()) + Error_Sum_phi[e][5]*TMath::Sin(75.*TMath::DegToRad()) + Error_Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()))/(TMath::Sin(70.*TMath::DegToRad()) + TMath::Sin(75.*TMath::DegToRad()) + TMath::Sin(80.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(70.7*TMath::DegToRad())-TMath::Cos(80.4*TMath::DegToRad())));


		Error_sh_acceptance[e][4] = (Error_Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()) + Error_Sum_phi[e][7]*TMath::Sin(85.*TMath::DegToRad()) + Error_Sum_phi[e][8]*TMath::Sin(88.*TMath::DegToRad()))/(TMath::Sin(80.*TMath::DegToRad()) + TMath::Sin(85.*TMath::DegToRad()) + TMath::Sin(88.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(80.4*TMath::DegToRad())-TMath::Cos(90.0*TMath::DegToRad())));
		
		
		TGraph *hxj = new TGraphErrors(5, sh_theta, sh_acceptance[e], Error_x, Error_sh_acceptance[e]);
		hxj->SetMarkerStyle(47);
		if(e==2) 
		hxj->Draw("PL same 3");
		hxj->SetMarkerColor(e+1);
		hxj->SetLineColor(e+1);
		hxj->SetFillColor(e+1);
		hxj->SetFillStyle(3005);
		
		hxj->SetLineWidth(3);
		
		

	}
	
	double sh_tm_theta[6] = {23.5, 40.9, 54.1, 65.25, 75.45, 85.2};
	double sh_tm_acceptance_low[9][6] = {{1.04e-5, 0.000655, 7.68e-5, 9.12e-7, 1.0e-12, 1.0e-12}, {0.000674, 0.002119, 0.001387, 0.001111, 6.71e-6, 1.0e-12}, {0.003701, 0.010341, 0.017542, 0.017035, 0.001286, 1.0e-12, }, {0.005421, 0.014283, 0.028896, 0.055125, 0.044884, 0.000904}, {0.00794, 0.021544, 0.055125, 0.088177, 0.136966, 0.00349}, {0.012333, 0.028896, 0.067702, 0.133003, 0.269074, 0.083148}, {0.020921, 0.038756, 0.083148, 0.163348, 0.470036, 0.649211}, {0.030644, 0.047599, 0.09351, 0.246389, 0.668554, 2.16349}, {0.083148, 0.163348, 0.360892, 0.751863, 1.56639, 2.901746}};
	
	
	double sh_tm_acceptance[9][6] = {{0.000134, 0.001725, 0.000617, 8.38e-5, 1.0e-12, 1.0e-12}, {0.001286, 0.003811, 0.004286, 0.004820, 9.27e-6, 1.0e-12}, {0.004545, 0.014283, 0.020921, 0.031557, 0.00771, 1.0e-12, }, {0.006096, 0.015598, 0.034462, 0.071796, 0.133003, 0.004286}, {0.009469, 0.024229, 0.061994, 0.108295, 0.206594, 0.007487}, {0.014283, 0.032497, 0.073935, 0.145249,0.340311,0.145249}, {0.023528,0.042325, 0.088177, 0.173227, 0.528607, 1.008424}, {0.034462, 0.050477, 0.102119, 0.246389, 0.708986, 2.580225}, {0.09351, 0.183703, 0.39412, 0.821088, 1.86811, 3.669974}};
	
	
	
	double sh_tm_acceptance_high[9][6] = {{0.000599, 0.002926, 0.002383, 0.000904, 5.0e-6, 1.0e-12}, {0.002119, 0.005749, 0.008671, 0.012333, 0.000271, 2.548e-6}, {0.005421, 0.017035, 0.024229, 0.049017, 0.027248, 0.000191}, {0.007487, 0.018065, 0.039911, 0.088177, 0.133003,  0.009469}, {0.010649, 0.02646, 0.06772, 0.12179, 0.269074, 0.018603}, {0.017035, 0.034462, 0.078407, 0.154033, 0.39412, 0.25373}, {0.024951, 0.044884, 0.090804, 0.189176, 0.577277, 1.3134}, {0.037635, 0.055125, 0.108295, 0.285347, 0.774264, 2.901746}, {0.10219, 0.200617, 0.417955, 0.870744, 2.040112, 4.507297}};

	double sh_tm_acceptance_error_high[9][6];
	double sh_tm_acceptance_error_low[9][6];
	for(int k = 0; k<9; k++){
		for(int i = 0; i<6; i++){
			sh_tm_acceptance_error_high[k][i] = sh_tm_acceptance_high[k][i] - sh_tm_acceptance[k][i];
			sh_tm_acceptance_error_low[k][i] =  sh_tm_acceptance[k][i] - sh_tm_acceptance_low[k][i];
		}
	}

	TMultiGraph *mg = new TMultiGraph();
	for(int k = 0; k<9; k++){
		TGraph *hxj = new TGraphAsymmErrors(6, sh_tm_theta, sh_tm_acceptance[k], Error_x, Error_x, sh_tm_acceptance_error_low[k], sh_tm_acceptance_error_high[k]);
		hxj->SetMarkerStyle(10);
		//if(k==3 || k==4) 
		//hxj->Draw("PL same 3");
		hxj->SetMarkerColor(k+1);
		hxj->SetLineColor(k+1);
		hxj->SetFillColor(k+1);
		hxj->SetFillStyle(3005);
		hxj->SetLineWidth(3);
		if(k==3 || k==4)
		mg->Add(hxj);
		
	}
	mg->Draw("PL");
	
	



}



















