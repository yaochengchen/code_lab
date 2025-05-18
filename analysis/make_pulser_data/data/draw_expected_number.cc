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

//double candidates_angle[11][2]={{6.7, 46.8},{16.6,-26.9},{12.4,24.7},{6.1,18.4},{26.4,5.5},{7.4,-5.6},{22.3,0.7},{13.3,20.2},{24.6,7.},{-0.75,27.2},{24.6,4.3}};

double candidates_angle[57][2]={{6.7, 20.8}, {16.6, -52.9}, {12.4, -1.3}, {6.1, -7.6}, {7.4, -31.6}, {13.3, -5.8}, {24.6, -19}, {-0.75, 1.2}, {2.61, -22.56}, {25.9, -31.3}, {14, -11.5}, {29, -13.9}, {8.7, -11.8}, {25.7, -31}, {7.4, -31.6}, {4.6, -16.69}, {25.3, -4.1}, {28.3, -7.6}, {9.7 -30.7}, {28.7, -12.1}, {9.3, -59.5}, {27.7, -44.5}, {34.4, -42.7}, {22, -62.5}, {10, -15.1}, {21, -16.3}, {11.7, -12.4}, {27.7, -43}, {0.5, -56.65}, {31.6, -23.8}, {15.4, 0.2}, {9.6, 1.1}, {2.1, 18.9}, {18.3, -54.7}, {28.3, -50.2}, {27.3, -1.9}, {2.9, -24.86}, {51, 7.4}, {19.3, -45.4}, {30.6, -12.1}, {12, -1.3}, {27, -7.9}, {17.7, -29.2}, {11, -37.6}, {9.3, -10.9}, {9.7, -24.7}, {-0.5, 24.2}, {15.7, -49.9}, {41.6, -32.2}, {30.7, -11.2}, {3.1, -31.1}, {31.7, 13.1}, {41.6, 21.8}, {31.4, 5.6}, {18.3, -43.9}, {21.6, 12.8}, {18.7, 11.9}};

//.x make_coreas_data.cc(1, "", 20210408, 18969, 18970)
void draw_expected_number(){


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
		//cout<<a_t<<"   "<<(int)a_t<<endl;
		analysis_efficiency[(int)a_e][(int)a_t][(int)a_p] = a_ef;
	}
	for(int i_theta = 0; i_theta<9; i_theta++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			for(int k=0; k<6; k++){
				//analysis_efficiency[k][i_theta][i_phi] = 1;
				//if(k==1&&i_theta>4)analysis_efficiency[k][i_theta][i_phi] = 0.15;
			}
		}
	}




	TNtuple *ntuple[6] = {NULL};
	
	
	TFile *f_ntuple_real = new TFile("ntuple_position_factor_1750.root");
	ntuple[0] = (TNtuple *) f_ntuple_real->Get("ntuple");
	
	f_ntuple_real = new TFile("ntuple_position_factor_1850.root");
	ntuple[1] = (TNtuple *) f_ntuple_real->Get("ntuple");
	

	double DSR_radius[9] = {750, 750, 710, 782, 870, 1020, 1720, 2550, 3200};
	double Energy_list[6] = {pow(10,17.5), pow(10,17.75), pow(10,18.00), pow(10,18.25), pow(10,18.50), pow(10,18.75)};


	
	float _Energy, _ShowerTheta, _ShowerPhi, _EffectedArea;
	double Energy[6][9][19], ShowerTheta[6][9][19], ShowerPhi[6][9][19], EffectedArea_tp[6][9][19], EffectedArea_pt[6][19][9], theta_list[10], phi_list[19];
	double Sum_phi[6][10] = {0.};
	double Sum_theta[6][19] = {0.};
	double Sum_phitheta[6] = {0.};
	
	double Error_EffectedArea_tp[6][9][19], Error_EffectedArea_pt[6][19][9];
	double Error_Sum_phi[6][10] = {0.};
	double Error_Sum_theta[6][19] = {0.};
	double Error_Sum_phitheta[6] = {0.}; 
	
	double EffectedArea[9][19][6];
	double Error_EffectedArea[9][19][6];
	
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
			
		
			
			double _Error_EffectedArea = sqrt((_EffectedArea/(3.14159*DSR_radius[i_theta]*DSR_radius[i_theta]))*1000)/1000.*(3.14159*DSR_radius[i_theta]*DSR_radius[i_theta]);
			
			EffectedArea[i_theta][i_phi][k] = _EffectedArea*1.0e-6;
			
			theta_list[i_theta] = 90 - _ShowerTheta;
			phi_list[i_phi] = _ShowerPhi;
			
			Error_EffectedArea[i_theta][i_phi][k] = _Error_EffectedArea*1.0e-6;
			
			
			EffectedArea[i_theta][i_phi][k] = (EffectedArea[i_theta][i_phi][k]+0*Error_EffectedArea[i_theta][i_phi][k])*analysis_efficiency[k][i_theta][i_phi];
			
		}
	
	
	}
	
	
	TFile *f = new TFile("/media/cyc/For_Linux/CoREAS_Sim/proton/Auger5045-combined-spectrum-data-2019.root", "read"); //Dir_Sim

	TGraphAsymmErrors *gr_EJ = (TGraphAsymmErrors*) f->Get("gr_EJ");
	TGraphAsymmErrors *gr_logJ = (TGraphAsymmErrors*) f->Get("gr_logJ");
	
	
	//TCanvas *c = new TCanvas("c_cr", "cr spec", 0, 0, 1500, 2*800);
	//c->Divide(1,2);
	//c->cd(1);	gPad->SetGrid(1,1);
	//gPad->SetLogy(1);	
	//gr_EJ->Draw("ap");

	//c->cd(2);	gPad->SetGrid(1,1);
	//gr_logJ->Draw("ap");
	
	
	
	//PowLawProductIntegral(double logE1, double logE2, double logF1, double logF2, double logA1, double logA2)
	double c = PowLawProductIntegral(17.5, 17.75, gr_logJ->Eval(17.5), gr_logJ->Eval(17.75), log10(1000), log10(2000))*180*24*60*60;
	cout<<c<<endl;
	c=0;
	for(double E = 17.5; E<17.75; E+=0.001){
		c += pow(10,gr_logJ->Eval(E+0.0005)) * pow(10, (log10(2000./1000.)/0.25)*(E+0.0005-17.5)+log10(1000.)) *(pow(10,E+0.001)-pow(10,E)) * 180*24*60*60;
	}
	cout<<c<<endl;
	
	
	
	double EffectedArea_Flux[9][19][6];
	double EffectedArea_Flux_tp[6][9][19];
	double EffectedArea_Flux_pt[6][19][9];
	 
	for(int i_theta = 0; i_theta<9; i_theta++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			TGraph *wyn = new TGraph(6, Energy_list, EffectedArea[i_theta][i_phi]);
			for(int k=0; k<5; k++){
			
				EffectedArea_Flux[i_theta][i_phi][k] = PowLawProductIntegral(17.5+0.25*k, 17.75+0.25*k, gr_logJ->Eval(17.5+0.25*k), gr_logJ->Eval(17.75+0.25*k), log10(EffectedArea[i_theta][i_phi][k]*1.0e6+1.0e-12), log10(EffectedArea[i_theta][i_phi][k+1]*1.0e6+1.0e-12)) * 180*24*60*60*4.8;
				//for(double E = 17.5+0.25*k; E<17.75+0.25*k; E+=0.001){
				//	EffectedArea_Flux[i_theta][i_phi][k] += (pow(10,gr_logJ->Eval(E+0.0005)) * wyn->Eval(pow(10,E+0.0005))*1.0e6)*(pow(10,E+0.001)-pow(10,E)) * 180*24*60*60;
					//EffectedArea_Flux[i_theta][i_phi][k] = (pow(pow(10,0.25*k),-3.3)*pow(10,gr_logJ->Eval(17.5+0.25*k*0)) * wyn->Eval(pow(10,17.5+0.25*k))*1.0e6)*(pow(10,17.625+0.25*k)-pow(10,17.375+0.25*k)) * 180*24*60*60;
					//EffectedArea_Flux[i_theta][i_phi][k] = wyn->Eval(pow(10,17.5+0.25*k));
					EffectedArea_Flux_tp[k][i_theta][i_phi] = EffectedArea_Flux[i_theta][i_phi][k];
					EffectedArea_Flux_pt[k][i_phi][i_theta] = EffectedArea_Flux[i_theta][i_phi][k];
				//}
			}
			delete wyn;
			wyn = NULL;
		}
	}
			
	//TGraph *wyn = new TGraph(6, Energy_list, EffectedArea[3][10]);
	//wyn->Draw("APL");
			

	
	for(int k=0; k<5; k++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			TGraph *wyn = new TGraph(9, theta_list, EffectedArea_Flux_pt[k][i_phi]);
			for(double theta = 87.9; theta>=40.1; theta-=0.2){
				Sum_theta[k][i_phi] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*wyn->Eval(theta)/(TMath::Cos((40)*TMath::DegToRad())-TMath::Cos((88)*TMath::DegToRad())) *5*TMath::DegToRad();
			}
			delete wyn;
			wyn = NULL;
		}
	}
	
	
	
	TGraph *cyc[5] = {NULL};
	
	
	double a[2] = {40, 90};//{40, 90};//{-68, 30};
	double b[2] = {1.0e-25, 6e1};
	TGraph *wyn = new TGraph(2, a, b);
	wyn->Draw("AP");
	wyn->SetMarkerColor(0);
	
	for(int k=0; k<5; k++){
		cyc[k] = new TGraph(19, phi_list, Sum_theta[k]);
		cyc[k]->SetMarkerStyle(8);
		cyc[k]->SetMarkerColor(k+1);
		cyc[k]->SetLineColor(k+1);
		//cyc[k]->Draw("PL same");

	}

	
	
	
	double Sum_e_theta[19] = {0.};
	for(int k=0; k<5; k++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			Sum_e_theta[i_phi] += Sum_theta[k][i_phi];
			
		}
	}
	
	TGraph *abc3 = new TGraph(19, phi_list, Sum_e_theta);
	//abc3->Draw("PL same");
	
	
	double total_num = 0;
	
	for(int k=0; k<5; k++){
		for(int i_phi = 0; i_phi<19; i_phi++){
			total_num += Sum_theta[k][i_phi]*5*TMath::DegToRad();
			
		}
	}
	
	cout<<"total number:  "<<total_num<<endl;
	




	for(int k=0; k<5; k++){
		for(int i_theta = 0; i_theta<9; i_theta++){
			for(int i_phi = 0; i_phi<19; i_phi++){
				Sum_phi[k][i_theta] += EffectedArea_Flux_tp[k][i_theta][i_phi]*5*TMath::DegToRad();
			}
		}
	}
	
	for(int k=0; k<5; k++){
		Sum_phi[k][0] *= (10*TMath::DegToRad());//35~45
		Sum_phi[k][1] *= (10*TMath::DegToRad());//45~55
		Sum_phi[k][2] *= (7.5*TMath::DegToRad());//55~57.5
		Sum_phi[k][3] *= (5*TMath::DegToRad());//57.5~62.5
		Sum_phi[k][4] *= (5*TMath::DegToRad());//62.5~67.5
		Sum_phi[k][5] *= (5*TMath::DegToRad());//67.5~72.5
		Sum_phi[k][6] *= (5*TMath::DegToRad());//72.5~77.5
		Sum_phi[k][7] *= (5*TMath::DegToRad());//77.5~82.5
		Sum_phi[k][8] *= (4*TMath::DegToRad());//82.5~86.5
		Sum_phi[k][9] *= (3.5*TMath::DegToRad());//86.5~90
	}
	
	
	for(int k=0; k<5; k++){
		cyc[k] = new TGraph(9, theta_list, Sum_phi[k]);
		cyc[k]->SetMarkerStyle(8);
		cyc[k]->SetMarkerColor(k+1);
		cyc[k]->SetLineColor(k+1);
		cyc[k]->Draw("PL same");

	}
	
	TH1F *energy_distribution = new TH1F("h1", "h1", 12, 17, 20);
	
	double Sum_e_phi[9] = {0.};
	for(int k=0; k<5; k++){
		for(int i_theta = 0; i_theta<9; i_theta++){
			Sum_e_phi[i_theta] += Sum_phi[k][i_theta];
			energy_distribution->Fill(17.625 + k*0.25, Sum_phi[k][i_theta]);
		}
	}
	
	TGraph *abc = new TGraph(9, theta_list, Sum_e_phi);
	abc->Draw("PL same");
	
	double Sum_all = 0;
	for(int i=0; i<9; i++){
		Sum_all += Sum_e_phi[i];
	}
	cout<<"all:  "<<Sum_all<<endl;
	
	
	TCanvas *canvas1 = new TCanvas("c1", "c1");
	canvas1->cd();
	energy_distribution->Draw();
	for(int i=1; i<13; i++){
		cout<<energy_distribution->GetBinCenter(i)<<"   "<<energy_distribution->GetBinContent(i)<<endl;
	}
	

	TCanvas *canvas2 = new TCanvas("c2", "c2");
	canvas2->cd();
	



	double sh_theta[5] = {40.9, 54.1, 65.25, 75.45, 85.2};
	double sh_acceptance[6][5] = {0.};
	for(int e=4; e>=0; e--){
		for(double theta = 48.2; theta>=33.6; theta-=0.2){
			sh_acceptance[e][0] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta);
			}
		for(double theta = 60.0; theta>=48.2; theta-=0.2){
			sh_acceptance[e][1] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta);
			}
		for(double theta = 70.5; theta>=60.0; theta-=0.2){
			sh_acceptance[e][2] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta);
			}
			for(double theta = 80.4; theta>=70.5; theta-=0.2){
			sh_acceptance[e][3] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta);
			}
			for(double theta = 88.0; theta>=80.4; theta-=0.2){
			sh_acceptance[e][4] += (TMath::Cos((theta-0.1)*TMath::DegToRad())-TMath::Cos((theta+0.1)*TMath::DegToRad()))*cyc[e]->Eval(theta);
			}
		
		TGraph *hxj = new TGraph(5, sh_theta, sh_acceptance[e]);
		hxj->SetMarkerStyle(9);
		//if(e==4) hxj->Draw("APL");
		hxj->GetYaxis()->SetRangeUser(0,5);
		//hxj->Draw("PL same");
		hxj->SetMarkerColor(e+1);
		hxj->SetLineColor(e+1);
		hxj->SetLineWidth(3);
	}
	

	double Sum_e[9] = {0.};
	for(int k=0; k<5; k++){
		for(int i_theta = 0; i_theta<5; i_theta++){
			Sum_e[i_theta] += sh_acceptance[k][i_theta];
		}
	}
	
	TGraph *abc2 = new TGraph(5, sh_theta, Sum_e);
	//abc2->Draw("PL same");
	abc2->SetLineWidth(3);	
	
	
	
	double candidates_theta[57]={0};//{0,0,4,3,3};
	double candidates_phi[57]={0.};
	//double candidates_theta_error[5]={0,0,2,1.732,1.732};
	double all_one[57];
	for(int i=0; i<57; i++){
		candidates_phi[i] = candidates_angle[i][1];//-26;
		candidates_theta[i] = 90-candidates_angle[i][0];
		all_one[i] = 1;
	}
	
	double phi_number[19] = {0.};
	double phi_number_error[19] = {0.};
	for(int k=0; k<19; k++){
		for(int i=0; i<57; i++){
			if((candidates_phi[i]>(phi_list[k]-2.5))&&(candidates_phi[i]<=(phi_list[k]+2.5))){
				phi_number[k] ++;
			}
		}
		phi_number_error[k] = sqrt(phi_number[k]);
	}
	
	double theta_number[5] = {0.};
	double theta_number_error[5] = {0.};
	
	for(int i=0; i<57; i++){
		if((candidates_theta[i]>=33.6)&&(candidates_theta[i]<48.2)){
			theta_number[0] ++;
			//cout<<i<<"  "<<candidates_theta[i]<<endl;
		}
		
		if((candidates_theta[i]>=48.2)&&(candidates_theta[i]<60.0)){
			theta_number[1] ++;
		}
		
		if((candidates_theta[i]>=60.0)&&(candidates_theta[i]<70.5)){
			theta_number[2] ++;
		}
		
		if((candidates_theta[i]>=70.5)&&(candidates_theta[i]<80.4)){
			theta_number[3] ++;
		}
		
		if((candidates_theta[i]>=80.4)&&(candidates_theta[i]<88.0)){
			theta_number[4] ++;
		}
		
	}
		
	for(int k=0; k<5; k++){
		theta_number_error[k] = sqrt(theta_number[k]);
		cout<<theta_number[k]<<endl;
	}
	
	
	TGraphErrors *syt = new TGraphErrors(5, sh_theta, theta_number,0,theta_number_error);
	//syt->Draw("P same");
	syt->SetMarkerStyle(8);
	syt->SetMarkerSize(1);



	/*
	double theta_number[9] = {0.};
	double theta_number_error[9] = {0.};
	
	for(int i=0; i<57; i++){
		if((candidates_theta[i]>=35)&&(candidates_theta[i]<45)){
			theta_number[0] ++;
			//cout<<i<<"  "<<candidates_theta[i]<<endl;
		}
		
		if((candidates_theta[i]>=45)&&(candidates_theta[i]<55)){
			theta_number[1] ++;
		}
		
		if((candidates_theta[i]>=55)&&(candidates_theta[i]<62.5)){
			theta_number[2] ++;
		}
		
		if((candidates_theta[i]>=62.5)&&(candidates_theta[i]<67.5)){
			theta_number[3] ++;
		}
		
		if((candidates_theta[i]>=67.5)&&(candidates_theta[i]<72.5)){
			theta_number[4] ++;
		}
		
		if((candidates_theta[i]>=72.5)&&(candidates_theta[i]<77.5)){
			theta_number[5] ++;
		}
		
		if((candidates_theta[i]>=77.5)&&(candidates_theta[i]<82.5)){
			theta_number[6] ++;
		}
		
		if((candidates_theta[i]>=82.5)&&(candidates_theta[i]<86.5)){
			theta_number[7] ++;
		}
		
		if((candidates_theta[i]>=86.5)&&(candidates_theta[i]<90)){
			theta_number[8] ++;
		}
		
	}
		
	for(int k=0; k<9; k++){
		theta_number_error[k] = sqrt(theta_number[k]);
		cout<<theta_number[k]<<endl;
	}
	
	
	TGraphErrors *syt = new TGraphErrors(9, theta_list, theta_number,0,theta_number_error);
	syt->Draw("P same");
	syt->SetMarkerStyle(8);
	syt->SetMarkerSize(1);
	*/
	
	
	
	TGraphErrors *ysm = new TGraphErrors(19, phi_list, phi_number,0,phi_number_error);
	//ysm->Draw("P same");
	ysm->SetMarkerStyle(8);
	ysm->SetMarkerSize(1);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	double sh_tm_acceptance[9][6] = {{0.000134, 0.001725, 0.000617, 8.38e-5, 1.0e-12, 1.0e-12}, {0.001286, 0.003811, 0.004286, 0.004820, 9.27e-6, 1.0e-12}, {0.004545, 0.014283, 0.020921, 0.031557, 0.00771, 1.0e-12, }, {0.006096, 0.015598, 0.034462, 0.071796, 0.133003, 0.004286}, {0.009469, 0.024229, 0.061994, 0.108295, 0.206594, 0.007487}, {0.014283, 0.032497, 0.073935, 0.145249,0.340311,0.145249}, {0.023528,0.042325, 0.088177, 0.173227, 0.528607, 1.008424}, {0.034462, 0.050477, 0.102119, 0.246389, 0.708986, 2.580225}, {0.09351, 0.183703, 0.39412, 0.821088, 1.86811, 3.669974}};
	
	double sh_tm_energy[8] = {17.125, 17.375, 17.625, 17.875, 18.125, 18.375, 18.75, 19.25};
	double sh_tm_event[9][6] = {0.};
	double sh_tm_acceptance_te[6][9] = {0.};
	double sh_tm_theta[6] = {23.5, 40.9, 54.1, 65.25, 75.45, 85.2};
	
	for(int i=0; i<6; i++){
		for(int e=0; e<8; e++){
			sh_tm_acceptance_te[i][e] = sh_tm_acceptance[e][i];
		}
	}
	
	for(int i=0; i<6; i++){
		TGraph *wyn = new TGraph(9, sh_tm_energy, sh_tm_acceptance_te[i]);
		for(int e=0; e<6; e++){
			for(double E = 17.0+0.25*e; E<17.25+0.25*e; E+=0.001){
				
				sh_tm_event[e][i] += (pow(10,gr_logJ->Eval(E+0.0005)) * wyn->Eval(E+0.0005)*1.0e6)*(pow(10,E+0.001)-pow(10,E)) * 180*24*60*60;	
				//gr_EJ
				//sh_tm_event[e][i] += (gr_EJ->Eval(E+0.0005) * wyn->Eval(E+0.0005)*1.0e6)*(pow(10,E+0.001)-pow(10,E)) * 180*24*60*60;
			}
			
			sh_tm_event[e][i] = PowLawProductIntegral(17.0+0.25*e, 17.25+0.25*e, gr_logJ->Eval(17.0+0.25*e), gr_logJ->Eval(17.25+0.25*e), log10(sh_tm_acceptance_te[i][e]), log10(sh_tm_acceptance_te[i][e])) * 180*24*60*60;
		}
	}
	
	
	
	TMultiGraph *mg = new TMultiGraph();
	for(int k = 0; k<6; k++){
		TGraph *hxj = new TGraph(6, sh_tm_theta, sh_tm_event[k]);
		hxj->SetMarkerStyle(10);
		//if(k==3 || k==4) 
		//hxj->Draw("PL same 3");
		hxj->SetMarkerColor(k+1);
		hxj->SetLineColor(k+1);
		hxj->SetFillColor(k+1);
		hxj->SetFillStyle(3005);
		hxj->SetLineWidth(3);
		//if(k==6 || k==6)
		mg->Add(hxj);
		//cout<<"f"<<endl;
		
	}
	
	
	
	double sh_Sum_e[6] = {0.};
	for(int k=0; k<6; k++){
		for(int i_theta = 0; i_theta<6; i_theta++){
			sh_Sum_e[i_theta] += sh_tm_event[k][i_theta];
		}
	}
	
	TGraph *abc1 = new TGraph(6, sh_tm_theta, sh_Sum_e);
	mg->Add(abc1);
	//mg->Draw("APL");
	
	
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




