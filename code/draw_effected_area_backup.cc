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
	double analysis_efficiency[8][10][19] = {0.};
	TNtuple *analysis_ntuple = (TNtuple *) f_ntuple_analysis->Get("ntuple");
	analysis_ntuple->SetBranchAddress("e", &a_e);
	analysis_ntuple->SetBranchAddress("i_theta", &a_t);
	analysis_ntuple->SetBranchAddress("i_phi", &a_p);
	analysis_ntuple->SetBranchAddress("efficiency", &a_ef);
	for(int i=0; i<analysis_ntuple->GetEntries(); i++){
		analysis_ntuple->GetEntry(i);
		analysis_efficiency[(int)a_e][(int)a_t][(int)a_p] = a_ef;
	}



	TFile *f_ntuple_real[16] = {NULL};
	TNtuple *ntuple[16] = {NULL};
	
	
	f_ntuple_real[0] = new TFile("ntuple_position_factor_with_pps_1750.root");
	ntuple[0] = (TNtuple *) f_ntuple_real[0]->Get("ntuple");
	
	f_ntuple_real[1] = new TFile("ntuple_position_factor_with_pps_1775.root");
	ntuple[1] = (TNtuple *) f_ntuple_real[1]->Get("ntuple");
	
	f_ntuple_real[2] = new TFile("ntuple_position_factor_with_pps_1800.root");
	ntuple[2] = (TNtuple *) f_ntuple_real[2]->Get("ntuple");
	
	f_ntuple_real[3] = new TFile("ntuple_position_factor_with_pps_1825.root");
	ntuple[3] = (TNtuple *) f_ntuple_real[3]->Get("ntuple");
	
	f_ntuple_real[4] = new TFile("ntuple_position_factor_with_pps_1850.root");
	ntuple[4] = (TNtuple *) f_ntuple_real[4]->Get("ntuple");
	
	f_ntuple_real[5] = new TFile("ntuple_position_factor_with_pps_1875.root");
	ntuple[5] = (TNtuple *) f_ntuple_real[5]->Get("ntuple");
	
	f_ntuple_real[6] = new TFile("ntuple_position_factor_with_pps_1900.root");
	ntuple[6] = (TNtuple *) f_ntuple_real[6]->Get("ntuple");
	
	f_ntuple_real[7] = new TFile("ntuple_position_factor_with_pps_1950.root");
	ntuple[7] = (TNtuple *) f_ntuple_real[7]->Get("ntuple");
	
	f_ntuple_real[8] = new TFile("ntuple_position_factor_with_pps_17500.root");
	ntuple[8] = (TNtuple *) f_ntuple_real[8]->Get("ntuple");
	
	f_ntuple_real[9] = new TFile("ntuple_position_factor_with_pps_17750.root");
	ntuple[9] = (TNtuple *) f_ntuple_real[9]->Get("ntuple");
	
	f_ntuple_real[10] = new TFile("ntuple_position_factor_with_pps_18000.root");
	ntuple[10] = (TNtuple *) f_ntuple_real[10]->Get("ntuple");
	
	f_ntuple_real[11] = new TFile("ntuple_position_factor_with_pps_18250.root");
	ntuple[11] = (TNtuple *) f_ntuple_real[11]->Get("ntuple");
	
	f_ntuple_real[12] = new TFile("ntuple_position_factor_with_pps_18500.root");
	ntuple[12] = (TNtuple *) f_ntuple_real[12]->Get("ntuple");
	
	f_ntuple_real[13] = new TFile("ntuple_position_factor_with_pps_18750.root");
	ntuple[13] = (TNtuple *) f_ntuple_real[13]->Get("ntuple");
	
	f_ntuple_real[14] = new TFile("ntuple_position_factor_with_pps_19000.root");
	ntuple[14] = (TNtuple *) f_ntuple_real[14]->Get("ntuple");
	
	f_ntuple_real[15] = new TFile("ntuple_position_factor_with_pps_19500.root");
	ntuple[15] = (TNtuple *) f_ntuple_real[15]->Get("ntuple");
	
	
	
	

	


	
	float _Energy, _ShowerTheta, _ShowerPhi, _EffectedArea, _DSR_radius, _Triggered_count, _Total_count;
	double Energy[8][29][19], ShowerTheta[8][29][19], ShowerPhi[8][29][19], EffectedArea_tp[8][29][19], EffectedArea_pt[8][19][29], theta_list[29], phi_list[19], DSR_radius[8][29][19];
	double Ave_phi[8][29] = {0.};
	double Ave_theta[8][19] = {0.};
	double Ave_phitheta[8] = {0.};
	
	double EffectedArea[29][19][8];
	
	double Error_EffectedArea_tp[8][29][19], Error_EffectedArea_pt[8][19][29];
	double Error_Ave_phi[8][29] = {0.};
	double Error_Ave_theta[8][19] = {0.};
	double Error_Ave_phitheta[8] = {0.}; 
	
	double Error_x[29] = {0.};
	
	for(int k=0; k<16; k++){
		ntuple[k]->SetBranchAddress("Energy", &_Energy);
		ntuple[k]->SetBranchAddress("ShowerTheta", &_ShowerTheta);
		ntuple[k]->SetBranchAddress("ShowerPhi", &_ShowerPhi);
		ntuple[k]->SetBranchAddress("EffectedArea", &_EffectedArea);
		ntuple[k]->SetBranchAddress("DSR_radius", &_DSR_radius);
		ntuple[k]->SetBranchAddress("Triggered_count", &_Triggered_count);
		ntuple[k]->SetBranchAddress("Total_count", &_Total_count);
	}
	

	
	double fitting_APT[19][9] = {0.};
	double fitting_EAPT[19][9] = {0.};
	double fitting_theta[9] = {50., 40., 30., 25., 20., 15., 10., 5., 2.};
	for(int i =0; i<ntuple[7]->GetEntries(); i++){
			int i_phi = i%19;
			int	i_theta = (i/19)%9;
			ntuple[7]->GetEntry(i);
			fitting_APT[i_phi][i_theta] = _EffectedArea*1.0e-6;
			double _Error_EffectedArea = 0.;
			if(_Triggered_count>0){_Error_EffectedArea = _EffectedArea/sqrt(_Triggered_count);}
			//cout<<_Error_EffectedArea<<endl;
			fitting_EAPT[i_phi][i_theta] = _Error_EffectedArea*1.0e-6;
	}
	
	
		
	for(int n=0; n<16; n++){
		if(n==7){continue;}
		for(int i =0; i<ntuple[n]->GetEntries(); i++){
			int i_phi = i%19;
			int i_theta = (i/19)%28;
			
			if((n==5)&&(i_theta>10)){i_theta+=1;}//since 1875 62 degree forgot to run
			if((n==6)&&(i_theta>19)){i_theta+=1;}//since 1900 76 degree forgot to run
			
			
			if(n>7){i_theta = 28;}
		
			ntuple[n]->GetEntry(i);
			int k = 0;
			switch((int) round(_Energy)){
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
				case 1900:
					k = 6;
					break;
				case 1950:
					k = 7;
					break;
				case 17500:
					k = 0;
					break;
				case 17750:
					k = 1;
					break;
				case 18000:
					k = 2;
					break;
				case 18250:
					k = 3;
					break;
				case 18500:
					k = 4;
					break;
				case 18750:
					k = 5;
					break;
				case 19000:
					k = 6;
					break;
				case 19500:
					k = 7;
					break;
			}
			Energy[k][i_theta][i_phi] = _Energy/100.;
			ShowerTheta[k][i_theta][i_phi] = 90 - _ShowerTheta;
			ShowerPhi[k][i_theta][i_phi] = _ShowerPhi;
			
			//_EffectedArea *= analysis_efficiency[k][i_theta][i_phi];
			//_EffectedArea *= pow(pow(10,(5-k)*0.25),2.3);
			
			double _Error_EffectedArea = 0.;
			if(_Triggered_count>0){_Error_EffectedArea = _EffectedArea/sqrt(_Triggered_count);}
			//cout<<_Error_EffectedArea<<endl;
			
			EffectedArea[i_theta][i_phi][k] = _EffectedArea*1.0e-6;
			EffectedArea_tp[k][i_theta][i_phi] = _EffectedArea*1.0e-6;
			EffectedArea_pt[k][i_phi][i_theta] = _EffectedArea*1.0e-6;
			theta_list[i_theta] = 90 - _ShowerTheta;
			phi_list[i_phi] = _ShowerPhi;
			
			Error_EffectedArea_tp[k][i_theta][i_phi] = _Error_EffectedArea*1.0e-6;
			Error_EffectedArea_pt[k][i_phi][i_theta] = _Error_EffectedArea*1.0e-6;
			//double fuck=_Error_EffectedArea*1.0e-6;
			//if(fuck!=fuck){cout<<fuck<<endl;}
		}
	
	
	}
	
	for(int i=0; i<19; i++){
		for(int j=0; j<29; j++){
			if(j<4||j>21){
		
				EffectedArea[j][i][0] = 0;
				EffectedArea_tp[0][j][i] = 0;
				EffectedArea_pt[0][i][j] = 0;
			
				Error_EffectedArea_tp[0][j][i] = 0;
				Error_EffectedArea_pt[0][i][j] = 0;
			}
		}
	}
	
	for(int i=0; i<19; i++){
		for(int j=12; j<13; j++){
			Energy[2][j][i] = 18.0;
			ShowerTheta[2][j][i] = ShowerTheta[0][j][i];
			ShowerPhi[2][j][i] = ShowerPhi[0][j][i];
			//cout<<ShowerTheta[0][j][i]<<endl;
			double A = (EffectedArea[11][i][2] + EffectedArea[13][i][2])/2.;
			double EA = (Error_EffectedArea_tp[2][11][i] + Error_EffectedArea_tp[2][13][i])/2.;
			
			EffectedArea[j][i][2] = A;
			EffectedArea_tp[2][j][i] = A;
			EffectedArea_pt[2][i][j] = A;
			
			Error_EffectedArea_tp[2][j][i] = EA;
			Error_EffectedArea_pt[2][i][j] = EA;
		}
		EffectedArea[28][i][2] = 0;
		EffectedArea_tp[2][28][i] = 0;
		EffectedArea_pt[2][i][28] = 0;
			
		Error_EffectedArea_tp[2][28][i] = 0;
		Error_EffectedArea_pt[2][i][28] = 0;
	}
		
	for(int i=0; i<19; i++){
		for(int j=11; j<12; j++){
			Energy[5][j][i] = 18.75;
			ShowerTheta[5][j][i] = ShowerTheta[0][j][i];
			ShowerPhi[5][j][i] = ShowerPhi[0][j][i];
			
			double A = (EffectedArea[10][i][5] + EffectedArea[12][i][5])/2.;
			double EA = (Error_EffectedArea_tp[5][10][i] + Error_EffectedArea_tp[5][12][i])/2.;
			
			EffectedArea[j][i][5] = A;
			EffectedArea_tp[5][j][i] = A;
			EffectedArea_pt[5][i][j] = A;
			
			Error_EffectedArea_tp[5][j][i] = EA;
			Error_EffectedArea_pt[5][i][j] = EA;
		}
	}
	
	for(int i=0; i<19; i++){
		for(int j=20; j<21; j++){
			Energy[6][j][i] = 19.0;
			ShowerTheta[6][j][i] = ShowerTheta[0][j][i];
			ShowerPhi[6][j][i] = ShowerPhi[0][j][i];
			
			double A = (EffectedArea[19][i][6] + EffectedArea[21][i][6])/2.;
			double EA = (Error_EffectedArea_tp[6][19][i] + Error_EffectedArea_tp[6][21][i])/2.;
			
			EffectedArea[j][i][6] = A;
			EffectedArea_tp[6][j][i] = A;
			EffectedArea_pt[6][i][j] = A;
			
			Error_EffectedArea_tp[6][j][i] = EA;
			Error_EffectedArea_pt[6][i][j] = EA;
		}
	}
	
	for(int i=0; i<19; i++){
		TGraph *hxj = new TGraph(9, fitting_theta, fitting_APT[i]);
		TGraph *syt = new TGraph(9, fitting_theta, fitting_EAPT[i]);
		for(int j=0; j<28; j++){
			Energy[7][j][i] = 19.5;
			ShowerTheta[7][j][i] = ShowerTheta[0][j][i];
			ShowerPhi[7][j][i] = ShowerPhi[0][j][i];
			
			double A = hxj->Eval(90-ShowerTheta[7][j][i]);
			double EA = syt->Eval(90-ShowerTheta[7][j][i]);
			
			EffectedArea[j][i][7] = A;
			EffectedArea_tp[7][j][i] = A;
			EffectedArea_pt[7][i][j] = A;
			
			Error_EffectedArea_tp[7][j][i] = EA;
			Error_EffectedArea_pt[7][i][j] = EA;
		}
	}
	
	
	TFile *f = new TFile("/media/cyc/For_Linux/CoREAS_Sim/proton/Auger5045-combined-spectrum-data-2019.root", "read"); //Dir_Sim

	TGraphAsymmErrors *gr_EJ = (TGraphAsymmErrors*) f->Get("gr_EJ");
	TGraphAsymmErrors *gr_logJ = (TGraphAsymmErrors*) f->Get("gr_logJ");
	
	


	
	for(int k=0; k<8; k++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			TGraph *wyn = new TGraph(29, theta_list, EffectedArea_pt[k][i_phi]);
			for(double theta = 89.9; theta>=40.1; theta-=0.2){
				Ave_theta[k][i_phi] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*wyn->Eval(theta)/(TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((90)*TMath::DegToRad()));
				//Ave_theta[k][i_phi] *= (TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((90)*TMath::DegToRad()));
			}
			delete wyn;
			wyn = NULL;
			
			
			
			TGraph *hxj = new TGraph(29, theta_list, Error_EffectedArea_pt[k][i_phi]);
			for(double theta = 89.99; theta>=40.01; theta-=0.02){
				Error_Ave_theta[k][i_phi] += (TMath::Cos((theta-0.01)*TMath::DegToRad())-TMath::Cos((theta+0.01)*TMath::DegToRad()))*hxj->Eval(theta)/(TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((90)*TMath::DegToRad()));
				//Error_Ave_theta[k][i_phi] *= (TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((90)*TMath::DegToRad()));
			}
			delete hxj;
			hxj = NULL;
		}
		

		for(int i_theta = 0; i_theta<29; i_theta++){
			for(int i_phi = 0; i_phi<19; i_phi++){
				Ave_phi[k][i_theta] += EffectedArea_tp[k][i_theta][i_phi]/19.;
				//if(i_theta==2){cout<<e<<"  "<<i_phi<<"  "<<EffectedArea_tp[e][i_theta][i_phi]/19.<<endl;}
				Error_Ave_phi[k][i_theta] += Error_EffectedArea_tp[k][i_theta][i_phi]/19.;
				//Ave_phi[k][i_theta] *= (TMath::Pi()/2.);
				//Error_Ave_phi[k][i_theta] *= (TMath::Pi()/2.);
			}
			//cout<<Ave_phi[e][2]<<endl;
		}
		
	
	}
	

// sum phi first
///////////////////////////////////////////////////////////////////////////////////////////////////////////	
		
	
	TGraphErrors *cyc[8] = {NULL};
	TGraphAsymmErrors *cycA[8] = {NULL};
	
	
	TCanvas* c_A = new TCanvas("c_1", "area acceptance", 1200, 800);
	c_A->SetLogy();
	double a[2] = {0, 90};
	double b[2] = {1.0e-5, 10.0e10};
	TGraph *wyn = new TGraph(2, a, b);
	wyn->Draw("AP");
	wyn->SetMarkerColor(0);
	
	for(int k=0; k<8; k++){
		cyc[k] = new TGraphErrors(29, theta_list, Ave_phi[k], Error_x, Error_Ave_phi[k]);
		cyc[k]->SetMarkerStyle(8);
		cyc[k]->SetMarkerColor(k+1);
		cyc[k]->SetLineColor(k+1);
		cyc[k]->SetFillColor(k+1);
		cyc[k]->SetFillStyle(3005);
		//cyc[k]->Draw("PL same 3");
		
	}
	
	for(int i=0; i<8; i++){
		//cout<<theta_list[i]<<endl;
		//cout<<Energy[i][0][0]<<endl;
	}

	
	double Bin_E_list[8] = {17.375, 17.625, 17.875, 18.125, 18.375, 18.625, 18.875, 19.25};
	double Ave_phi_Bin_E[3][8][29] = {0.};//lower bound, normal, upper bound
	double Bin_E_Ave_phi[3][29][8] = {0.};
	
	for(int i_theta = 0; i_theta<29; i_theta++){
			for(int k=0; k<8; k++){
				double a_energy = 17.25;
				if(k>0){a_energy = Energy[k-1][0][0];}
				
				double b_energy = Energy[k][0][0];
				
				for(int l=0; l<3; l++){
					double a_EA = 1.0e-12;
					if(k>0){
						a_EA = Ave_phi[k-1][i_theta]*1.0e6 + (l-1.)*Error_Ave_phi[k-1][i_theta]*1.0e6 + 1.0e-12;
					}
				
					double b_EA = Ave_phi[k][i_theta]*1.0e6 + (l-1.)*Error_Ave_phi[k][i_theta] + 1.0e-12;
					
					if(a_EA<0){a_EA=1.0e-12;}
					if(b_EA<0){b_EA=1.0e-12;}
				
					double a = PowLawProductIntegral(a_energy, b_energy, gr_logJ->Eval(a_energy), gr_logJ->Eval(b_energy), log10(a_EA), log10(b_EA));
					double b = PowLawProductIntegral(a_energy, b_energy, gr_logJ->Eval(a_energy), gr_logJ->Eval(b_energy), log10(1), log10(1));
					//cout<<a<<"  "<<b<<endl;
					//double EA = a/b/1.0e6;
					double EA = a*180*24*60*60*4;
					//EA = 1;//b_EA;
					Ave_phi_Bin_E[l][k][i_theta] = EA;
					Bin_E_Ave_phi[l][i_theta][k] = EA;
				}
			}		
	}
	
	
	double BinEdge_theta_list[7] = {40, 50.33, 59.29, 67.48, 75.20, 82.66, 90};
	double BinCenter_theta_list[6] = {0.};
	
	double Bin_theta_Bin_E[3][6][8] = {0.};
	double Bin_E_Bin_theta[3][8][6] = {0.};
	
	for(int k=0; k<8; k++){
		for(int l=0; l<3; l++){
			TGraph *hxj = new TGraph(29, theta_list, Ave_phi_Bin_E[l][k]);
			for(int i=0; i<6; i++){
				BinCenter_theta_list[i] = (BinEdge_theta_list[i] + BinEdge_theta_list[i+1])/2.;
				double acceptance = 0;
				for(double theta = BinEdge_theta_list[i]+0.01; theta <= BinEdge_theta_list[i+1]-0.01; theta += 0.02){
					acceptance += (TMath::Cos((theta-0.01)*TMath::DegToRad())-TMath::Cos((theta+0.01)*TMath::DegToRad()))*hxj->Eval(theta);
				}
				acceptance *= (TMath::Pi()/2.);
				Bin_theta_Bin_E[l][i][k] = acceptance;
				Bin_E_Bin_theta[l][k][i] = acceptance;
			}
			delete hxj;
			hxj = NULL;
		}
		
	}
	
	double Errors_Bin_theta_Bin_E[2][6][8] = {0.};
	double Errors_Bin_E_Bin_theta[2][8][6] = {0.};
	for(int k=0; k<8; k++){
		for(int i=0; i<6; i++){
			Errors_Bin_theta_Bin_E[0][i][k] = Bin_theta_Bin_E[1][i][k] - Bin_theta_Bin_E[0][i][k];
			Errors_Bin_theta_Bin_E[1][i][k] = Bin_theta_Bin_E[2][i][k] - Bin_theta_Bin_E[1][i][k];
			Errors_Bin_E_Bin_theta[0][k][i] = Bin_theta_Bin_E[1][i][k] - Bin_theta_Bin_E[0][i][k];
			Errors_Bin_E_Bin_theta[1][k][i] = Bin_theta_Bin_E[2][i][k] - Bin_theta_Bin_E[1][i][k];
		}
	}
	
	for(int k=0; k<8; k++){
		cycA[k] = new TGraphAsymmErrors(6, BinCenter_theta_list, Bin_E_Bin_theta[1][k], Error_x, Error_x, Errors_Bin_E_Bin_theta[0][k], Errors_Bin_E_Bin_theta[1][k]);
		cycA[k]->SetMarkerStyle(8);
		cycA[k]->SetMarkerColor(k+1);
		cycA[k]->SetLineColor(k+1);
		cycA[k]->SetFillColor(k+1);
		cycA[k]->SetFillStyle(3005);
		cycA[k]->Draw("PL same 3");
		
	}
	
	
	for(int k=0; k<6; k++){
		cycA[k] = new TGraphAsymmErrors(8, Bin_E_list, Bin_theta_Bin_E[1][k], Error_x, Error_x, Errors_Bin_theta_Bin_E[0][k], Errors_Bin_theta_Bin_E[1][k]);
		cycA[k]->SetMarkerStyle(8);
		cycA[k]->SetMarkerColor(k+1);
		cycA[k]->SetLineColor(k+1);
		cycA[k]->SetFillColor(k+1);
		cycA[k]->SetFillStyle(3005);
		//cycA[k]->Draw("PL same 3");
		
	}
	
	double sum_1 = 0;
	for(int k=0; k<6; k++){
		for(int i=0; i<8; i++){
			sum_1 += Bin_theta_Bin_E[1][k][i];
		}
	}
	cout<<sum_1<<endl;
	
////////////////////////////////////////////////////////////////////////////////////////////////



// sum theta first
///////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	
	
	TGraphErrors *syt[8] = {NULL};
	TGraphAsymmErrors *sytA[8] = {NULL};
	
	
	TCanvas* c_B = new TCanvas("c_2", "area acceptance", 1200, 800);
	c_B->SetLogy();
	double c[2] = {-70, 20};
	double d[2] = {1.0e-5, 10.0e10};
	TGraph *ysm = new TGraph(2, c, d);
	ysm->Draw("AP");
	ysm->SetMarkerColor(0);
	
	for(int k=0; k<8; k++){
		syt[k] = new TGraphErrors(19, phi_list, Ave_theta[k], Error_x, Error_Ave_theta[k]);
		syt[k]->SetMarkerStyle(8);
		syt[k]->SetMarkerColor(k+1);
		syt[k]->SetLineColor(k+1);
		syt[k]->SetFillColor(k+1);
		syt[k]->SetFillStyle(3005);
		//syt[k]->Draw("PL same 3");
		
	}
	
	

	
	//double Bin_E_list[8] = {17.375, 17.625, 17.875, 18.125, 18.375, 18.625, 18.875, 19.25};
	double Ave_theta_Bin_E[3][8][19] = {0.};//lower bound, normal, upper bound
	double Bin_E_Ave_theta[3][19][8] = {0.};
	
	for(int i_phi = 0; i_phi<19; i_phi++){
			for(int k=0; k<8; k++){
				double a_energy = 17.25;
				if(k>0){a_energy = Energy[k-1][0][0];}
				
				double b_energy = Energy[k][0][0];
				
				for(int l=0; l<3; l++){
					double a_EA = 1.0e-12;
					if(k>0){
						a_EA = Ave_theta[k-1][i_phi]*1.0e6 + (l-1)*Error_Ave_theta[k-1][i_phi]*1.0e6 + 1.0e-12;
					}
				
					double b_EA = Ave_theta[k][i_phi]*1.0e6 + (l-1)*Error_Ave_theta[k][i_phi]*1.0e6 + 1.0e-12;
					
					if(a_EA<0){a_EA=1.0e-12;}
					if(b_EA<0){b_EA=1.0e-12;}
					
					double a = PowLawProductIntegral(a_energy, b_energy, gr_logJ->Eval(a_energy), gr_logJ->Eval(b_energy), log10(a_EA), log10(b_EA));
					double b = PowLawProductIntegral(a_energy, b_energy, gr_logJ->Eval(a_energy), gr_logJ->Eval(b_energy), log10(1), log10(1));
					//cout<<a<<"  "<<b<<endl;
					double EA = a/b/1.0e6;
					//EA = 1;//b_EA;
					Ave_theta_Bin_E[l][k][i_phi] = EA;
					Bin_E_Ave_theta[l][i_phi][k] = EA;
				}
			}		
	}
	
	
	//double BinEdge_phi_list[7] = {0.};
	double BinCenter_phi_list[6] = {0.};
	
	double Bin_phi_Bin_E[3][6][8] = {0.};
	double Bin_E_Bin_phi[3][8][6] = {0.};
	
	for(int k=0; k<8; k++){
		for(int l=0; l<3; l++){
			for(int i=0; i<6; i++){
				BinCenter_phi_list[i] = (phi_list[1+i*3] + phi_list[2+i*3])/2.;
				double acceptance = (0.5*Ave_theta_Bin_E[l][k][i*3] + Ave_theta_Bin_E[l][k][i*3+1] + Ave_theta_Bin_E[l][k][i*3+2] + 0.5*Ave_theta_Bin_E[l][k][i*3+3])*(5*TMath::DegToRad());
			
				acceptance *= (TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((90)*TMath::DegToRad()));
			
				Bin_phi_Bin_E[l][i][k] = acceptance;
				Bin_E_Bin_phi[l][k][i] = acceptance;
				//cout<<BinCenter_phi_list[i]<<"   "<<acceptance<<endl;
			}
		}
	}
	
	
	double Errors_Bin_phi_Bin_E[2][6][8] = {0.};
	double Errors_Bin_E_Bin_phi[2][8][6] = {0.};
	for(int k=0; k<8; k++){
		for(int i=0; i<6; i++){
			Errors_Bin_phi_Bin_E[0][i][k] = Bin_phi_Bin_E[1][i][k] - Bin_phi_Bin_E[0][i][k];
			Errors_Bin_phi_Bin_E[1][i][k] = Bin_phi_Bin_E[2][i][k] - Bin_phi_Bin_E[1][i][k];
			Errors_Bin_E_Bin_phi[0][k][i] = Bin_phi_Bin_E[1][i][k] - Bin_phi_Bin_E[0][i][k];
			Errors_Bin_E_Bin_phi[1][k][i] = Bin_phi_Bin_E[2][i][k] - Bin_phi_Bin_E[1][i][k];
		}
	}
	
	
	for(int k=0; k<8; k++){
		sytA[k] = new TGraphAsymmErrors(6, BinCenter_phi_list, Bin_E_Bin_phi[1][k], Error_x, Error_x, Errors_Bin_E_Bin_phi[0][k], Errors_Bin_E_Bin_phi[1][k]);
		sytA[k]->SetMarkerStyle(8);
		sytA[k]->SetMarkerColor(k+1);
		sytA[k]->SetLineColor(k+1);
		sytA[k]->SetFillColor(k+1);
		sytA[k]->SetFillStyle(3005);
		//sytA[k]->Draw("PL same 3");
		
	}
	
	
	for(int k=0; k<6; k++){
		sytA[k] = new TGraphAsymmErrors(8, Bin_E_list, Bin_phi_Bin_E[1][k], Error_x, Error_x, Errors_Bin_phi_Bin_E[0][k], Errors_Bin_phi_Bin_E[1][k]);
		sytA[k]->SetMarkerStyle(8);
		sytA[k]->SetMarkerColor(k+1);
		sytA[k]->SetLineColor(k+1);
		sytA[k]->SetFillColor(k+1);
		sytA[k]->SetFillStyle(3005);
		sytA[k]->Draw("PL same 3");
		
	}
	
	double sum_2 = 0;
	for(int k=0; k<6; k++){
		for(int i=0; i<8; i++){
			sum_2 += Bin_phi_Bin_E[1][k][i];
		}
	}
	cout<<sum_2<<endl;
	
	
	
////////////////////////////////////////////////////////////////////////////////////////////////

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












