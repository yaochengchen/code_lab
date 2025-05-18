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

double PowLawProductIntegral(double logE1, double logE2, double logF1, double logF2, double logA1, double logA2);

//TRandom *r3;

//.x make_coreas_data.cc(1, "", 20210408, 18969, 18970)
void draw_effected_area(){


	TFile * f_ntuple_analysis = new TFile("/home/cyc/software/TAROGE-4_analysis/draw/ntuple_analysis_efficiency.root");
	float a_e, a_t, a_p, a_ef;
	double analysis_efficiency[6][9][19] = {0.};
	TNtuple *analysis_ntuple = (TNtuple *) f_ntuple_analysis->Get("ntuple");
	analysis_ntuple->SetBranchAddress("e", &a_e);
	analysis_ntuple->SetBranchAddress("i_theta", &a_t);
	analysis_ntuple->SetBranchAddress("i_phi", &a_p);
	analysis_ntuple->SetBranchAddress("efficiency", &a_ef);
	for(int i=0; i<analysis_ntuple->GetEntries(); i++){
		analysis_ntuple->GetEntry(i);
		analysis_efficiency[(int)a_e][(int)a_t][(int)a_p] = a_ef;
	}




	TNtuple *ntuple[6] = {NULL};
	
	
	TFile *f_ntuple_real = new TFile("ntuple_position_factor_1750.root");
	ntuple[0] = (TNtuple *) f_ntuple_real->Get("ntuple");
	
	f_ntuple_real = new TFile("ntuple_position_factor_1850.root");
	ntuple[1] = (TNtuple *) f_ntuple_real->Get("ntuple");
	

	double DSR_radius[9] = {750, 750, 710, 782, 870, 1020, 1720, 2550, 3200};


	
	float _Energy, _ShowerTheta, _ShowerPhi, _EffectedArea;
	double Energy[6][9][19], ShowerTheta[6][9][19], ShowerPhi[6][9][19], EffectedArea_tp[6][9][19], EffectedArea_pt[6][19][9], theta_list[10], phi_list[19];
	double Sum_phi[6][10] = {0.};
	double Sum_theta[6][19] = {0.};
	double Sum_phitheta[6] = {0.};
	
	double EffectedArea[9][19][6];
	
	double Error_EffectedArea_tp[6][9][19], Error_EffectedArea_pt[6][19][9];
	double Error_Sum_phi[6][10] = {0.};
	double Error_Sum_theta[6][19] = {0.};
	double Error_Sum_phitheta[6] = {0.}; 
	
	double Error_x[20] = {0.};
	
	for(int k=0; k<2; k++){
		ntuple[k]->SetBranchAddress("Energy", &_Energy);
		ntuple[k]->SetBranchAddress("ShowerTheta", &_ShowerTheta);
		ntuple[k]->SetBranchAddress("ShowerPhi", &_ShowerPhi);
		ntuple[k]->SetBranchAddress("EffectedArea", &_EffectedArea);
	}
	

	
	for(int i =0; i<ntuple[0]->GetEntries(); i++){
		int i_phi = i%19;
		int i_theta = (i/19)%9;
		
		
		for(int n=0; n<2; n++){
			ntuple[n]->GetEntry(i);
			int k = 0;
			switch((int) round(_Energy*100)){
				case 1750:
					k = 0;
					break;
				case 1775:
					k = 1;
					break;
				case 1800:
					k = 2;
					break;
				case 1825:
					k = 3;
					break;
				case 1850:
					k = 4;
					break;
				case 1875:
					k = 5;
					break;
			}
			Energy[k][i_theta][i_phi] = _Energy;
			ShowerTheta[k][i_theta][i_phi] = 90 - _ShowerTheta;
			ShowerPhi[k][i_theta][i_phi] = _ShowerPhi;
			
			//_EffectedArea *= analysis_efficiency[k][i_theta][i_phi];
			//_EffectedArea *= pow(pow(10,(5-k)*0.25),2.3);
			
			double _Error_EffectedArea = sqrt((_EffectedArea/(3.14159*DSR_radius[i_theta]*DSR_radius[i_theta]))*1000)/1000.*(3.14159*DSR_radius[i_theta]*DSR_radius[i_theta]);
			
			EffectedArea[i_theta][i_phi][k] = _EffectedArea*1.0e-6;
			EffectedArea_tp[k][i_theta][i_phi] = _EffectedArea*1.0e-6;
			EffectedArea_pt[k][i_phi][i_theta] = _EffectedArea*1.0e-6;
			theta_list[i_theta] = 90 - _ShowerTheta;
			phi_list[i_phi] = _ShowerPhi;
			
			Error_EffectedArea_tp[k][i_theta][i_phi] = _Error_EffectedArea*1.0e-6;
			Error_EffectedArea_pt[k][i_phi][i_theta] = _Error_EffectedArea*1.0e-6;
		}
	
	
	}
	
	for(int i_theta = 0; i_theta<9; i_theta++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			for(int k=0; k<5; k++){
				double a = PowLawProductIntegral(17.5+0.25*k,17.75+0.25*k,-3.3*(17.5+0.25*k), -3.3*(17.75+0.25*k),log10(EffectedArea[i_theta][i_phi][k]*1.0e6+1.0e-12),log10(EffectedArea[i_theta][i_phi][k+1]*1.0e6+1.0e-12));
				double b = PowLawProductIntegral(17.5+0.25*k,17.75+0.25*k,-3.3*(17.5+0.25*k), -3.3*(17.75+0.25*k),log10(1),log10(1));
				//cout<<a<<"  "<<b<<endl;
				double EAA = a/b/1.0e6;
				//EffectedArea_tp[k][i_theta][i_phi] = EAA;
				//EffectedArea_pt[k][i_phi][i_theta] = EAA;
			
			}
				
		}
	}
	
	//TGraph *wyn = new TGraph(9, theta_list, EffectedArea_pt[3][0]);
	//wyn->Draw();
	
	for(int k=0; k<6; k++){
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
	
	
	
	TGraphErrors *cyc[6] = {NULL};
	
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
	
	TCanvas* c_A = new TCanvas("c_2", "area acceptance", 1200, 800);
	c_A->SetLogy();
	double a[2] = {0, 90};
	double b[2] = {1.0e-5, 10.0e0};
	TGraph *wyn = new TGraph(2, a, b);
	wyn->Draw("AP");
	wyn->SetMarkerColor(0);
	
	for(int k=0; k<6; k++){
		//Sum_phi[k][6] = 0.;
 		//Sum_phi[k][6] = 0.;
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
		cyc[k]->Draw("PL same 3");
		
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
	
	double fuck[5] = {0.};
	for(int e=0; e<6; e++){
		for(int f=0; f<5; f++){
			fuck[f] += sh_acceptance[e][f];
		}
	}
	
	TGraph *hxj = new TGraph(5, sh_theta, fuck);
		hxj->SetMarkerStyle(9);
		//hxj->Draw("PL same");
	
	
	for(int e=0; e<6; e++){
/*	
		sh_acceptance[e][0] = Sum_phi[e][0]*3.14159/2.*(TMath::Cos(33.6*TMath::DegToRad())-TMath::Cos(48.2*TMath::DegToRad()));
	
		sh_acceptance[e][1] = (Sum_phi[e][1]*TMath::Sin(50.*TMath::DegToRad()) + Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()))/(TMath::Sin(50.*TMath::DegToRad()) + TMath::Sin(60.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())));

		//cout<<(Sum_phi[e][1]+Sum_phi[e][2])/2.*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())))<<endl;
		//cout<<sh_acceptance[1]<<endl;
	
		sh_acceptance[e][2] = (Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()) + Sum_phi[e][3]*TMath::Sin(65.*TMath::DegToRad()) + Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()))/(TMath::Sin(60.*TMath::DegToRad()) + TMath::Sin(65.*TMath::DegToRad()) + TMath::Sin(70.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(60.0*TMath::DegToRad())-TMath::Cos(70.5*TMath::DegToRad())));

		sh_acceptance[e][3] = (Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()) + Sum_phi[e][5]*TMath::Sin(75.*TMath::DegToRad()) + Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()))/(TMath::Sin(70.*TMath::DegToRad()) + TMath::Sin(75.*TMath::DegToRad()) + TMath::Sin(80.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(70.7*TMath::DegToRad())-TMath::Cos(80.4*TMath::DegToRad())));


		sh_acceptance[e][4] = (Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()) + Sum_phi[e][6]*TMath::Sin(85.*TMath::DegToRad()) + Sum_phi[e][8]*TMath::Sin(88.*TMath::DegToRad()))/(TMath::Sin(80.*TMath::DegToRad()) + TMath::Sin(85.*TMath::DegToRad()) + TMath::Sin(88.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(80.4*TMath::DegToRad())-TMath::Cos(90.0*TMath::DegToRad())));
	
*/





		Error_sh_acceptance[e][0] = Error_Sum_phi[e][0]*3.14159/2.*(TMath::Cos(33.6*TMath::DegToRad())-TMath::Cos(48.2*TMath::DegToRad()));
	
		Error_sh_acceptance[e][1] = (Error_Sum_phi[e][1]*TMath::Sin(50.*TMath::DegToRad()) + Error_Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()))/(TMath::Sin(50.*TMath::DegToRad()) + TMath::Sin(60.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())));

		//cout<<(Sum_phi[e][1]+Sum_phi[e][2])/2.*(3.14159/2.*(TMath::Cos(48.2*TMath::DegToRad())-TMath::Cos(60.0*TMath::DegToRad())))<<endl;
		//cout<<sh_acceptance[1]<<endl;
	
		Error_sh_acceptance[e][2] = (Error_Sum_phi[e][2]*TMath::Sin(60.*TMath::DegToRad()) + Error_Sum_phi[e][3]*TMath::Sin(65.*TMath::DegToRad()) + Error_Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()))/(TMath::Sin(60.*TMath::DegToRad()) + TMath::Sin(65.*TMath::DegToRad()) + TMath::Sin(70.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(60.0*TMath::DegToRad())-TMath::Cos(70.5*TMath::DegToRad())));

		Error_sh_acceptance[e][3] = (Error_Sum_phi[e][4]*TMath::Sin(70.*TMath::DegToRad()) + Error_Sum_phi[e][5]*TMath::Sin(75.*TMath::DegToRad()) + Error_Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()))/(TMath::Sin(70.*TMath::DegToRad()) + TMath::Sin(75.*TMath::DegToRad()) + TMath::Sin(80.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(70.7*TMath::DegToRad())-TMath::Cos(80.4*TMath::DegToRad())));


		Error_sh_acceptance[e][4] = (Error_Sum_phi[e][6]*TMath::Sin(80.*TMath::DegToRad()) + Error_Sum_phi[e][6]*TMath::Sin(85.*TMath::DegToRad()) + Error_Sum_phi[e][8]*TMath::Sin(88.*TMath::DegToRad()))/(TMath::Sin(80.*TMath::DegToRad()) + TMath::Sin(85.*TMath::DegToRad()) + TMath::Sin(88.*TMath::DegToRad()))*(3.14159/2.*(TMath::Cos(80.4*TMath::DegToRad())-TMath::Cos(90.0*TMath::DegToRad())));
		
		
		TGraph *hxj = new TGraphErrors(5, sh_theta, sh_acceptance[e], Error_x, Error_sh_acceptance[e]);
		hxj->SetMarkerStyle(47);
		//if(e==5) 
		//hxj->Draw("L same 3");
		hxj->SetMarkerColor(2);
		hxj->SetLineColor(2);
		hxj->SetFillColor(2);
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
			sh_tm_acceptance_error_high[k][i] = (sh_tm_acceptance_high[k][i] - sh_tm_acceptance[k][i]);
			sh_tm_acceptance_error_low[k][i] =  (sh_tm_acceptance[k][i] - sh_tm_acceptance_low[k][i]);
			sh_tm_acceptance[k][i] *= 1.0;//1.6;
		}
	}

	TMultiGraph *mg = new TMultiGraph();
	for(int k = 0; k<9; k++){
		TGraph *hxj = new TGraphAsymmErrors(6, sh_tm_theta, sh_tm_acceptance[k], Error_x, Error_x, sh_tm_acceptance_error_low[k], sh_tm_acceptance_error_high[k]);
		hxj->SetMarkerStyle(10);
		//if(k==7) 
		//hxj->Draw("L same 3");
		hxj->SetMarkerColor(1);
		hxj->SetLineColor(1);
		hxj->SetFillColor(1);
		hxj->SetFillStyle(3005);
		hxj->SetLineWidth(3);
		//if(k==6 || k==6)
		//mg->Add(hxj);
		
	}
	//mg->Draw("PL same");
	
	



}












/*
	Ref: Peter Gorham,  ANITA note
	integration  area acceptance * flux = event rate
	piecewise power law

	E: energy
	F: flux
	A: area acceptance
*/
double PowLawProductIntegral(double logE1, double logE2, double logF1, double logF2, double logA1, double logA2)
{

	if(logE1 == logE2)  return 0.; //zero width
	double DLogE = logE2 - logE1; //for variable logE bin
	double DLogF = logF2 - logF1;

	double alpha = ( logA2 - logA1 )/ DLogE; //power index of area acceptance
	double beta =  DLogF / DLogE; //power index of flux
	double ab1 = alpha + beta + 1.; //+1 for integral
	double rate = 0.;

	if(ab1==0.) rate = DLogE * pow( 10., -(alpha+beta)*logE1 ) ;
	else{
		// [(E2)^(a+b+1) - (E1)^(a+b+1)]/(a+b+1)   from energy integral
		//steep cutoff in Ae(E) may cause large ab1  and 10^ab1 explode
	 	//rate =  ( pow( 10.,  ab1*logEs[e+1] ) - pow( 10., ab1*logEs[e] ) ) / ab1; // eV^ab1

	 	//better, smaller index; extract E1
	 	rate = pow( 10., logE1 ) * ( pow( 10., ab1 * DLogE ) - 1. ) / ab1; // eV^ab1
	 }

	 rate *=  pow(10., logA1 + logF1);// pow(10., logA1) * pow(10., logF1 ) ;

	 if(rate != rate)  //NaN
	 	printf("\tWarning: rate = NaN logE %.3f-%.3f logF= %.3f ~ %.3f  logA= %.3f ~ %.3f  %.3f + %.3f +1 = %.3f  %.3g\n", logE1, logE2, logF1, logF2, logA1, logA2, alpha, beta, ab1, rate);

	 return rate;
}












