#include <cstring>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST
#include <iostream>
#include <vector>
#include "../Reconstruction_tchain_new.h"
#include "TComplex.h"

#define Speed_Of_Light 2.99792458e8
#define WINDOW_SPACING 150
#define HALF_WINDOW_SPACING 75
#define HALF_WINDOW_WIDTH 76
#define PI TMath::Pi()

using namespace std;

template <typename T>
T** new_Array2D(int row, int col);

template <typename T> 
void delete_Array2D(T **arr, int row, int col);

void processing();
void initial();
void draw_candidate_event_coherent_sum_deconvoluted_CLEAN(int date, int run_num, int event_num, double in_theta, double in_phi, int polarity, int filter);
void get_candidate_file(string path, vector<string> &filenames);

void get_FFT(int channel, bool label[2049]);
void find_CW(int channel, bool label[2049]);
//HALF_WINDOW_WIDTH
void SW_fft(double SW_x[1500], int Max_i, double SW_re[751], double SW_im[751]);

void Planck_taper_windowing(double epsilon, int N, double x[], double w_x[]);

double coor_ant[8][3] = {{328914.252, 2697625.953, 709.389},{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949},{328936.041, 2697639.685, 704.949}};
double delay_ant[8] = {0,-92.4,53.2,-30.3,84.6,59.5,-1523.1,-743.1};


struct Ground_Respond{
  vector<vector<vector<vector<double>>>> resp_mag;
  vector<vector<vector<vector<double>>>> resp_pha;
};
// 0 for low elevation, 1 for high elevation
Ground_Respond gd_responds[2] = {{vector<vector<vector<vector<double>>>> (141, vector<vector<vector<double>>>(143, vector<vector<double>>(4, vector<double>(2049, 0)))), vector<vector<vector<vector<double>>>> (141, vector<vector<vector<double>>>(143, vector<vector<double>>(4, vector<double>(2049, 0))))}, {vector<vector<vector<vector<double>>>> (43, vector<vector<vector<double>>>(15, vector<vector<double>>(4, vector<double>(2049, 0)))), vector<vector<vector<vector<double>>>> (43, vector<vector<vector<double>>>(15, vector<vector<double>>(4, vector<double>(2049, 0))))}};

void read_responses_high();
void deconvolute_response(double x[4096], int ch, int is_high, int pos_theta, int pos_phi);

Double_t noise_mag[8][2049] = {0};
void get_noise_level(double noise_mag[8][2049]);

Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j);
Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j);
void coherent_add(int pol, double in_theta, double in_phi, double f_x[8][1500], double added_x[15000]);

Double_t t[1500]={0};
Double_t f_t[15000]={0};

	Double_t x[8][4096]={0};//8 antennas
	double Re_ant[8][2049] = {0.};
	double Im_ant[8][2049] = {0.};



double FFTCrossCorrelate(double ws1_re[2049], double ws1_im[2049], double ws2_re[2049], double ws2_im[2049], double ws_xcor[409600], double &tdoa_best, double &n_peak, int polarity);
void CLEAN_deconvolution(double Re_ant[8][2049], double Im_ant[8][2049], double ANT_V_theta_re[2049], double ANT_V_theta_im[2049], double E_theta[4][4096], double E_phi[4][4096], double ws_dec[409600], int polarity, int filter);
void CLEAN_FFT(double x[409600]);
void CLEAN_simulation(double x[4096], double E_field_x[4096], double ANT_V_theta_re[2049], double ANT_V_theta_im[2049], int channel);


void get_response();
void Quotient(double a_re, double a_im, double b_re, double b_im, double &o_re, double &o_im);
void Invert_complex(TComplex a, TComplex b, TComplex c, TComplex d, TComplex M, TComplex N, TComplex &x, TComplex &y);
void Deconvolute_antenna(double Re_ant[8][2049], double Im_ant[8][2049], double ANT_H_theta_re[2049], double ANT_H_theta_im[2049], double ANT_H_phi_re[2049], double ANT_H_phi_im[2049], double ANT_V_theta_re[2049], double ANT_V_theta_im[2049], double ANT_V_phi_re[2049], double ANT_V_phi_im[2049], double E_theta[4][4096], double E_phi[4][4096]);
void max_X_correlation_delay(double delay_theta[4], double delay_phi[4], double E_theta[4][4096], double E_phi[4][4096]);
double get_delay();
void Adjust_delay(double Re_ant[8][2049], double Im_ant[8][2049], double H_V_delay, double in_theta, double in_phi);

  Double_t ** FEE_re = new_Array2D<double>(8, 2049);
  Double_t ** FEE_im = new_Array2D<double>(8, 2049);
  
  double ANT_V_theta_re[37][73][2049] = {0.};//theta: -90, -85, ..., 90; phi:0, 5, ..., 360
  double ANT_V_theta_im[37][73][2049] = {0.};
  double ANT_V_phi_re[37][73][2049] = {0.};//theta: -90, -85, ..., 90; phi:0, 5, ..., 360
  double ANT_V_phi_im[37][73][2049] = {0.};
  double ANT_H_theta_re[37][73][2049] = {0.};//theta: -90, -85, ..., 90; phi:0, 5, ..., 360
  double ANT_H_theta_im[37][73][2049] = {0.};
  double ANT_H_phi_re[37][73][2049] = {0.};//theta: -90, -85, ..., 90; phi:0, 5, ..., 360
  double ANT_H_phi_im[37][73][2049] = {0.};
  
  double H_ant4_re[2049] = {0.};
  double H_ant4_im[2049] = {0.};
  double V_ant4_re[2049] = {0.};
  double V_ant4_im[2049] = {0.};
  
  Double_t * Correlator_time = new Double_t [409600];

  int N = 4096;
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"); 
  
  
int main(int argc, char *argv[]){
	initial();
	processing();
return 0;
}


void initial(){
double t4point_phi = 90 - 117.09;
double t4point_theta = -0.08;
double base_thick = 0.01;

//transfer t4 back side to front side
coor_ant[6][0] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Cos(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[6][1] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Sin(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[6][2] += TMath::Sin(t4point_theta*TMath::DegToRad()) * base_thick;

coor_ant[7][0] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Cos(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[7][1] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Sin(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[7][2] += TMath::Sin(t4point_theta*TMath::DegToRad()) * base_thick;



//set tower-2 to be origin
double set_t2_origin_photo[3] = {328919.115, 2697635.05, 707.86};//same as tower2
for(int i=0; i<8; i++){
for(int j=0; j<3; j++){
coor_ant[i][j] = coor_ant[i][j] - set_t2_origin_photo[j];
cout<<coor_ant[i][j]<<"  ";
}
cout<<endl;
}


get_response();


}

void processing(){
    cout<<"start 1"<<endl;
    string path = "/home/cyc/software/TAROGE-4_analysis/draw/candidate_events";
    vector<string> filenames;

	get_candidate_file(path, filenames);
	cout<<"start 2"<<endl;

	int date;
	int run_num;
	int event_num;
	double in_theta;
	double in_phi;
  
  	for(auto filename : filenames){
 
		cout<<"filename:  "<<filename<<endl;	
	
  		const char * c = filename.c_str();
  
  		FILE *fpin = fopen(c,"r");
  		while(1){
    	int p = fscanf(fpin,"%d %d %d %lf %lf", &date, &run_num, &event_num, &in_theta, &in_phi);
    	if (p!=5) break;
    	int polarity = 0;
    	int filter = 0;
    	draw_candidate_event_coherent_sum_deconvoluted_CLEAN(date, run_num, event_num, in_theta, in_phi, polarity, filter);
    	}
	
    }
}

void draw_candidate_event_coherent_sum_deconvoluted_CLEAN(int date, int run_num, int event_num, double in_theta, double in_phi, int polarity, int filter){






	const char * fname_evt;

	
	bool bGood;
	Int_t timeStamp_FPGA;
	double eventTime = 0;
	Short_t T1H[1500];
	Short_t T2H[1500];
	Short_t T3H[1500];
	Short_t T4H[1500];
	Short_t T1V[1500];
	Short_t T2V[1500];
	Short_t T3V[1500];
	Short_t T4V[1500];
	
	
	
 	
	Double_t E_field_x[8][4096]={0};//8 antennas
	


    fname_evt =  Form("/media/cyc/1p9TB/TAROGE4_DATA/data/%d/run%08d.root", date, run_num);
    //fname_evt = "/home/cyc/software/TAROGE-4_analysis/make_pulser_data/run20240101.root";

	TFile *file = new TFile(fname_evt);
    TTree *Tree_Muon = (TTree*) file->Get("t");
          Tree_Muon->SetBranchAddress("bGood",&bGood);
		  Tree_Muon->SetBranchAddress("timeStamp_FPGA",&timeStamp_FPGA);
		  Tree_Muon->SetBranchAddress("eventTime", &eventTime);
		  Tree_Muon->SetBranchAddress("T1H",T1H);
		  Tree_Muon->SetBranchAddress("T2H",T2H);
		  Tree_Muon->SetBranchAddress("T3H",T3H);
		  Tree_Muon->SetBranchAddress("T4H",T4H);
		  Tree_Muon->SetBranchAddress("T1V",T1V);
		  Tree_Muon->SetBranchAddress("T2V",T2V);
		  Tree_Muon->SetBranchAddress("T3V",T3V);
		  Tree_Muon->SetBranchAddress("T4V",T4V);
		  
	Tree_Muon->GetEntry(event_num);
	
	for (int i=0;i<15000;i++)
	{
	  f_t[i]=i*0.08;
	  for(int ch=0; ch<8; ch++){
	  	x[ch][i] = 0.;
	  }
	  
	}
  for (int i=0;i<1500;i++)
	{
	  if(T4H[i]>32512) cout<<T4H[i]<<endl;
	  t[i]=i*0.8;
	  x[0][i]= (float)T1H[i]*500/32512;
	  x[2][i]= (float)T2H[i]*500/32512;
	  x[4][i]= (float)T3H[i]*500/32512;
	  x[6][i]= (float)T4H[i]*500/32512;
	  //int k = i+3;
	  //if(k>=1500){k-=1500;}
	  x[1][i]= (float)T1V[i]*500/32512;  //caution cyc!!
	  x[3][i]= (float)T2V[i]*500/32512;
	  x[5][i]= (float)T3V[i]*500/32512;
	  x[7][i]= (float)T4V[i]*500/32512;
	  //x[1][i]= 0;//(float)T1V[i]*500/32512;  //caution cyc!!
	  //x[3][i]= 0;//(float)T2V[i]*500/32512;
	  //x[5][i]= 0;//(float)T3V[i]*500/32512;
	  //x[7][i]= 0;//(float)T4V[i]*500/32512;
	  //if(x[6][i]>500) cout<<x[6][i]<<endl;	  
	}
	/*
	x[1][497] = x[3][497] = x[5][497] = x[7][497] = -0.2;
	x[1][498] = x[3][498] = x[5][498] = x[7][498] = 0.2;
	x[1][499] = x[3][499] = x[5][499] = x[7][499] = 0.7;
	x[1][500] = x[3][500] = x[5][500] = x[7][500] = 1.;
	x[1][501] = x[3][501] = x[5][501] = x[7][501] = 0.7;
	x[1][502] = x[3][502] = x[5][502] = x[7][502] = 0.2;
	x[1][503] = x[3][503] = x[5][503] = x[7][503] = -0.25;
	*/
  	
	int i_theta = round(in_theta/5.) + 18;
	double tmpt_phi = in_phi-17;
	if(tmpt_phi<0){tmpt_phi+=360;}
	
	double tempt_phi = in_phi+26;
	if(tempt_phi<0){tempt_phi+=360;}
	
	int i_phi = round(tempt_phi/5.);
	//i_theta = 18;
	//i_phi = 0;	
	
	//CLEAN_simulation(x[1], E_field_x[1], ANT_V_theta_re[i_theta][i_phi], ANT_V_theta_im[i_theta][i_phi], 1);
	//CLEAN_simulation(x[3], E_field_x[3], ANT_V_theta_re[i_theta][i_phi], ANT_V_theta_im[i_theta][i_phi], 3);
	//CLEAN_simulation(x[5], E_field_x[5], ANT_V_theta_re[i_theta][i_phi], ANT_V_theta_im[i_theta][i_phi], 5);
	//CLEAN_simulation(x[7], E_field_x[7], ANT_V_theta_re[i_theta][i_phi], ANT_V_theta_im[i_theta][i_phi], 7);
	
	TCanvas *c4 = new TCanvas("c4","Canvas Example",1500,1000);
	c4->cd();
	double *t = new double[4096];
	for(int i=0; i<4096; i++){t[i] = i*0.8;}
	TGraph *cyc = new TGraph(4096, t, x[1]);
	cyc->Draw();
	string c4_fig_name = "./photos/" + std::to_string(1) + "_CLEAN_waveform_" + std::to_string(1) + ".png";
    c4->SaveAs(c4_fig_name.c_str());
    
    
	

	bool label[2049] = {false};
	for(int channel = 0; channel<8; channel++){
		//bool label[2049] = {false};
		find_CW(channel, label);
   	}
   	
	for(int channel = 0; channel<8; channel++){
		//bool label[2049] = {false};
   		get_FFT(channel, label);
   	}
   	
   	TCanvas *c5 = new TCanvas("c5","Canvas Example",1500,1000);
	c5->cd();
	TGraph *wyn = new TGraph(4096, t, x[1]);
	wyn->Draw();
	string c5_fig_name = "./photos/" + std::to_string(1) + "_CLEAN_waveform_" + std::to_string(2) + ".png";
    c5->SaveAs(c5_fig_name.c_str());



	in_phi = in_phi - 26.;
	
	// V_pol signal need to shift toward left side.
	double H_V_delay = get_delay();//ns
	cout<<"HV delay "<<H_V_delay<<endl;
	Adjust_delay(Re_ant, Im_ant, H_V_delay, in_theta, in_phi);
	
	//cout<<ANT_H_theta_re[18][0][700]<<endl;
	double E_theta[4][4096] = {0.};
	double E_phi[4][4096] = {0.};

	//Deconvolute_antenna(Re_ant, Im_ant, ANT_H_theta_re[i_theta][i_phi], ANT_H_theta_im[i_theta][i_phi], ANT_H_phi_re[i_theta][i_phi], ANT_H_phi_im[i_theta][i_phi], ANT_V_theta_re[i_theta][i_phi], ANT_V_theta_im[i_theta][i_phi], ANT_V_phi_re[i_theta][i_phi], ANT_V_phi_im[i_theta][i_phi], E_theta, E_phi);
	
	double *ws_dec = new double[409600];
	for(int i=0; i<409600; i++){ws_dec[i] = 0;}
	CLEAN_deconvolution(Re_ant, Im_ant, ANT_V_theta_re[i_theta][i_phi], ANT_V_theta_im[i_theta][i_phi],  E_theta, E_phi, ws_dec, polarity, filter);
	
	double delay_theta[4]={0.};
	double delay_phi[4]={0.};
	//max_X_correlation_delay(delay_theta, delay_phi, E_theta, E_phi);

	TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
	c1->Divide(2, 4);
	
	double *clean_time = new double[4096];
	for(int i=0; i<4096; i++){
		clean_time[i] = 0.8*i;
	
	}
   
   for(int ant = 0; ant < 4; ant++){
   		TGraph *cyc = new TGraph(4096, clean_time, E_theta[ant]);
   		TGraph *wyn = new TGraph(4096, clean_time, E_phi[ant]);
   		c1->cd(2*ant+1);
   		cyc->Draw("APL");
   		cyc->SetLineColor(ant+1);
   		c1->cd(2*ant+2);
   		wyn->Draw("APL");
   		wyn->SetLineColor(ant+1);
   		
   		}
   		
   	string c1_fig_name = "./" + std::to_string(date) + "_individual_waveform_" + std::to_string(run_num) + "_" + std::to_string(event_num) + "_cyc_label_filter" + ".png";
    c1->SaveAs(c1_fig_name.c_str());


    
    double added_E_theta[15000] = {0.};
    double added_E_phi[15000] = {0.};
    
    for(int ant=0; ant<4; ant++){
    	for(int i=0; i<4096; i++){
    		added_E_theta[i] += E_theta[ant][i]/4.;
    		added_E_phi[i] += E_phi[ant][i]/4.;
    	}
    }
    
    double max_E_theta = 0.;
    int max_E_theta_i = -1;
    double max_E_phi = 0.;
    int max_E_phi_i = -1;
    	
    for(int i=0; i<15000; i++){
    	if(fabs(added_E_theta[i]) > max_E_theta){
    		max_E_theta = fabs(added_E_theta[i]);
    		max_E_theta_i = i;
    	}
    	if(fabs(added_E_phi[i]) > max_E_phi){
    		max_E_phi = fabs(added_E_phi[i]);
    		max_E_phi_i = i;
    	}
    }    


    double RMS_theta = 0;
    double RMS_phi = 0;
    for(int i=0; i<15000; i++){
    	if((i<(max_E_theta_i-400))||(i>(max_E_theta_i+400))){
    		RMS_theta += pow(added_E_theta[i],2);
    	}
    	if((i<(max_E_phi_i-400))||(i>(max_E_phi_i+400))){
    		RMS_phi += pow(added_E_phi[i],2);
    	}
    }
    
    RMS_theta = sqrt(RMS_theta/(15000.-801.));
    RMS_phi = sqrt(RMS_phi/(15000.-801.));
    
    double theta_phi_ratio = max_E_phi/max_E_theta;
    double ratio_error = theta_phi_ratio*(RMS_theta/max_E_theta + RMS_phi/max_E_phi);
    

    TCanvas *c3 = new TCanvas("c3","Canvas Example",1500,1000);
    c3->cd();
    //c3->Divide(2, 1);
    
    TGraph *FW_waveform1 = new TGraph(4096, clean_time, added_E_theta);
    
    FW_waveform1->GetXaxis()->SetRangeUser(150,450);
    //FW_waveform1->GetYaxis()->SetRangeUser(-1.6,1.6);
    FW_waveform1->GetXaxis()->SetTitle("ns");
    FW_waveform1->GetYaxis()->SetTitle("mV/m");
    string title_name = std::to_string(theta_phi_ratio) + "+-" + std::to_string(ratio_error);
    //FW_waveform1->SetTitle(title_name.c_str());
    FW_waveform1->SetLineColor(2);

    //c3->cd(2);
    TGraph *FW_waveform2 = new TGraph(4096, clean_time, added_E_phi);
    
    FW_waveform2->GetXaxis()->SetRangeUser(150,450);
    //FW_waveform2->GetYaxis()->SetRangeUser(-1.6,1.6);
    FW_waveform2->GetXaxis()->SetTitle("ns");
    FW_waveform2->GetYaxis()->SetTitle("mV/m");
    FW_waveform2->SetTitle(title_name.c_str());
    FW_waveform2->SetLineColor(4);
    
    FW_waveform2->Draw("APL");
    FW_waveform1->Draw("PL same");

	TLegend * leg = new TLegend(0.75,0.75,0.9,0.9);
	leg->AddEntry(FW_waveform1,"E_{#theta}","LP");
	leg->AddEntry(FW_waveform2,"E_{#phi}","LP");
	leg->Draw();

        
    string fig_name = "./" + std::to_string(date) + "_waveform_" + std::to_string(run_num) + "_" + std::to_string(event_num) + "_cyc_label_filter" + ".png";
    c3->SaveAs(fig_name.c_str());
    
    //delete FW_waveform1;
    //delete FW_waveform2;
    //delete c3;



//delete Tree_Muon;
//file->Close();
//delete file;


}




int filter_N = 4096;
Int_t filter_n_size = filter_N;
TVirtualFFT *filter_fft_forward = TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
TVirtualFFT *filter_fft_back = TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");
double *time_signal = new double [4096];

void find_CW(int channel, bool label[2049]){

   filter_fft_forward->SetPoints(x[channel]);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(Re_ant[channel],Im_ant[channel]);
   

   
   // H-pol antenna 4  for HV_delay
   if(channel==6){
 	 for (Int_t k=327; k<330; k++){//0 means DC component 
        H_ant4_re[k] = Re_ant[channel][k]; 
        H_ant4_im[k] = Im_ant[channel][k];
                             }
   }

   // V-pol antenna 4  for HV_delay
   if(channel==7){
 	 for (Int_t k=327; k<330; k++){//0 means DC component 
        V_ant4_re[k] = Re_ant[channel][k]; 
        V_ant4_im[k] = Im_ant[channel][k];
                             }
   }
   
   
 for (Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
       if((k<606)||(k>1147)) {Re_ant[channel][k] = 0.; Im_ant[channel][k] = 0.;}
                             }
                               	     	   
  	   
  	   double power[2049] = {0.};
  	   for(Int_t k=0; k<2049; k++){
  	   		power[k] = Re_ant[channel][k] * Re_ant[channel][k] + Im_ant[channel][k] * Im_ant[channel][k];
  	   		}
  	   		
  	   	
  	   	// caution!! sort function will change value in power array!
  	   	double median_1 = TMath::Median(199, &power[655]);
  	   	
  	   	cout<<median_1<<endl;

  	   	double median_2 = TMath::Median(167, &power[915]);
  	   	
  	   	cout<<median_2<<endl;
  	   	
  	
	
	for(Int_t k=0; k<885; k++){
		if(power[k]/median_1>10){
			//cout<<k<<endl;
			label[k-3] = true;
			label[k-3] = true;
			label[k-2] = true;
			label[k-2] = true;
			label[k-1] = true;
			label[k-1] = true;
			label[k] = true;
			label[k] = true;
			label[k+1] = true;
			label[k+1] = true;
			label[k+2] = true;
			label[k+2] = true;
			label[k+3] = true;
			label[k+3] = true;
			
			if(power[k-4]/median_1>5){
				label[k-4] = true;
				label[k-4] = true;
				label[k-5] = true;
				label[k-5] = true;
				label[k-6] = true;
				label[k-6] = true;
				label[k-7] = true;
				label[k-7] = true;
			}
			if(power[k+4]/median_1>5){
				label[k+4] = true;
				label[k+4] = true;
				label[k+5] = true;
				label[k+5] = true;
				label[k+6] = true;
				label[k+6] = true;
				label[k+7] = true;
				label[k+7] = true;
			}
			if(power[k-8]/median_1>5){
				label[k-8] = true;
				label[k-8] = true;
				label[k-9] = true;
				label[k-9] = true;
				label[k-10] = true;
				label[k-10] = true;
				label[k-11] = true;
				label[k-11] = true;
			}
			if(power[k+8]/median_1>5){
				label[k+8] = true;
				label[k+8] = true;
				label[k+9] = true;
				label[k+9] = true;
				label[k+10] = true;
				label[k+10] = true;
				label[k+11] = true;
				label[k+11] = true;
			}
		}
	
	}

	for(Int_t k=885; k<2049; k++){
		if(power[k]/median_2>10){
			//cout<<k<<endl;
			label[k-3] = true;
			label[k-3] = true;
			label[k-2] = true;
			label[k-2] = true;
			label[k-1] = true;
			label[k-1] = true;
			label[k] = true;
			label[k] = true;
			label[k+1] = true;
			label[k+1] = true;
			label[k+2] = true;
			label[k+2] = true;
			label[k+3] = true;
			label[k+3] = true;
			
			if(power[k-4]/median_2>5){
				label[k-4] = true;
				label[k-4] = true;
				label[k-5] = true;
				label[k-5] = true;
				label[k-6] = true;
				label[k-6] = true;
				label[k-7] = true;
				label[k-7] = true;
			}
			if(power[k+4]/median_2>5){
				label[k+4] = true;
				label[k+4] = true;
				label[k+5] = true;
				label[k+5] = true;
				label[k+6] = true;
				label[k+6] = true;
				label[k+7] = true;
				label[k+7] = true;
			}
			if(power[k-8]/median_2>5){
				label[k-8] = true;
				label[k-8] = true;
				label[k-9] = true;
				label[k-9] = true;
				label[k-10] = true;
				label[k-10] = true;
				label[k-11] = true;
				label[k-11] = true;
			}
			if(power[k+8]/median_2>5){
				label[k+8] = true;
				label[k+8] = true;
				label[k+9] = true;
				label[k+9] = true;
				label[k+10] = true;
				label[k+10] = true;
				label[k+11] = true;
				label[k+11] = true;
			}
		}
	
	}



}

void get_FFT(int channel, bool label[2049]){

  	   
	double new_x[4096] = {0.};
	for(Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
       if((k<606)||(k>1147)) {Re_ant[channel][k] = 0.; Im_ant[channel][k] = 0.;}
       //if((k<327)||(k>329)) {Re_ant[k] = 0.; Im_ant[k] = 0.;}
       else{
       //Re_ant[channel][k]*= (!label[k]);
       //Im_ant[channel][k]*= (!label[k]);
       //Re_ant[channel][k]/= sqrt(4096.);//normalization
       //Im_ant[channel][k]/= sqrt(4096.);
       //cout<<Re_ant[k]<<"   "<<Im_ant[k]<<"   "<<endl;
       //FEE_re[0][k] = 1;
       //FEE_im[0][k] = 0;
       //cout<<log10(FEE_re[0][k]*FEE_re[0][k] + FEE_im[0][k]*FEE_im[0][k])<<endl;
       //Quotient(Re_ant[channel][k], Im_ant[channel][k], FEE_re[channel][k], FEE_im[channel][k], Re_ant[channel][k], Im_ant[channel][k]);
 
       }
    }
    
    fft_back->SetPointsComplex(Re_ant[channel], Im_ant[channel]);
    fft_back->Transform();
    Correlator_time = fft_back->GetPointsReal();
    
    double max = 0;
    double max_i = 0;
    for(Int_t k=0; k<1500; k++){
        new_x[k] = Correlator_time[k]/4096.;
        if(fabs(new_x[k])>max){max = fabs(new_x[k]); max_i = k;}
    }
    
    for(Int_t k=0; k<1500; k++){
        if((k<max_i-40)||(k>max_i+40)){
        	new_x[k] = 0;
        }
        x[channel][k] = new_x[k]; 
    }
    
   filter_fft_forward->SetPoints(new_x);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(Re_ant[channel],Im_ant[channel]);

   
   
}





template <typename T>
T** new_Array2D(int row, int col)  
{  
    int size = sizeof(T);  
    int point_size = sizeof(T*);  
    //先申请内存，其中sizeof(T*) * row表示存放row个行指针   
    T **arr = (T **) malloc(point_size * row + size * row * col);  
    if (arr != NULL)  
    {     
        T *head = (T*)((long)arr + point_size * row);  
        for (int i = 0; i < row; ++i)  
        {  
            arr[i] =  (T*)((long)head + i * col * size);  
            for (int j = 0; j < col; ++j)  
                new (&arr[i][j]) T{}; // {} means initial to 0
        }  
    }  
    return (T**)arr;  
}  
//释放二维数组   
template <typename T>  
void delete_Array2D(T **arr, int row, int col)  
{  
    for (int i = 0; i < row; ++i)  
        for (int j = 0; j < col; ++j)  
            arr[i][j].~T();  
    if (arr != NULL)  
        free((void**)arr);  
}


void Planck_taper_windowing(double epsilon, int N, double x[], double w_x[]){
	
	double _w_x[1500] = {0.};
	
	double weighting[1500] = {0.};
	weighting[0] = 0.;
	double epsilon_N = epsilon*N;
	for(int i=1; i<epsilon_N; i++){
		double factor = 1./(1 + exp(epsilon_N/i - epsilon_N/(epsilon_N-i)));
		weighting[i] = factor;
		}
	for(int i=epsilon_N; i<=N/2.; i++){
		weighting[i] = 1.;
		}
	
	for(int i=0; i<=N/2.; i++){
		weighting[N-i] = weighting[i];
		}
	for(int i = 0; i<N; i++){
		//_w_x[i] = weighting[i]*x[i];
		double new_value = weighting[i]*x[i];
		w_x[i] = new_value;
		//cout<<weighting[i]*x[i]<<endl;
		//cout<<w_x[i]<<endl;
		//w_x[i] = weighting[i];// * x[i];
		//cout<<i<<"   "<<w_x[i]<<endl;
		}
	}



void read_responses_high(){

        Int_t azimuth = 0.;
        Int_t elevation = 0.;
        Double_t responses_magnitude_1H[2049] = {0.};
        Double_t responses_phase_1H[2049] = {0.};
        Double_t responses_magnitude_2H[2049] = {0.};
        Double_t responses_phase_2H[2049] = {0.};
        Double_t responses_magnitude_3H[2049] = {0.};
        Double_t responses_phase_3H[2049] = {0.};
        Double_t responses_magnitude_4H[2049] = {0.};
        Double_t responses_phase_4H[2049] = {0.};
    
    const char * res_fname =  "/home/cyc/software/TAROGE-4_pulser/pulser/H_pol/ground_responses_high.root";
    TFile * res_file = new TFile(res_fname);
    TTree * responsesTree = (TTree*) res_file->Get("ground_responses");
    
        responsesTree->SetBranchAddress("azimuth", &azimuth);
        responsesTree->SetBranchAddress("elevation", &elevation);
        
        responsesTree->SetBranchAddress("responses_magnitude_1H", &responses_magnitude_1H);
        responsesTree->SetBranchAddress("responses_phase_1H", &responses_phase_1H);
        
        responsesTree->SetBranchAddress("responses_magnitude_2H", &responses_magnitude_2H);
        responsesTree->SetBranchAddress("responses_phase_2H", &responses_phase_2H);
        
        responsesTree->SetBranchAddress("responses_magnitude_3H", &responses_magnitude_3H);
        responsesTree->SetBranchAddress("responses_phase_3H", &responses_phase_3H);
        
        responsesTree->SetBranchAddress("responses_magnitude_4H", &responses_magnitude_4H);
        responsesTree->SetBranchAddress("responses_phase_4H", &responses_phase_4H);
        
   for(int entry = 0; entry < responsesTree->GetEntries(); entry++){
   		responsesTree->GetEntry(entry);
   		for(int i = 0; i<2049; i++){
   			if(elevation<0 || elevation>=43 || azimuth<0 || azimuth>=15){
   			cout<<elevation<<"   "<<azimuth<<endl;
   			}
   			//cout<<responses_magnitude_1H[i]<<endl;
   			
   			gd_responds[1].resp_mag[elevation][azimuth][0][i] = responses_magnitude_1H[i];
   			gd_responds[1].resp_pha[elevation][azimuth][0][i] = responses_phase_1H[i];
   			
   			gd_responds[1].resp_mag[elevation][azimuth][1][i] = responses_magnitude_2H[i];
   			gd_responds[1].resp_pha[elevation][azimuth][1][i] = responses_phase_2H[i];
   			
   			gd_responds[1].resp_mag[elevation][azimuth][2][i] = responses_magnitude_3H[i];
   			gd_responds[1].resp_pha[elevation][azimuth][2][i] = responses_phase_3H[i];
   			
   			gd_responds[1].resp_mag[elevation][azimuth][3][i] = responses_magnitude_4H[i];
   			gd_responds[1].resp_pha[elevation][azimuth][3][i] = responses_phase_4H[i];
   			
   			//if(responses_magnitude_1H[i]!=0) cout<<resp_mag[elevation][azimuth][0][i]<<"   "<<responses_magnitude_1H[i]<<endl;
   			
   			}
   		}
	
}


int decon_N = 4096;
Int_t decon_n_size = decon_N;
TVirtualFFT *decon_fft_forward = TVirtualFFT::FFT(1, &decon_n_size,"R2C ES K");
TVirtualFFT *decon_fft_backward = TVirtualFFT::FFT(1, &decon_N, "C2RBACKWARD M K");
Double_t * fft_back_decon_x = new Double_t[4096];

void deconvolute_response(double x[4096], int ch, int is_high, int pos_theta, int pos_phi){

	  
   int ant = ch/2;
   double re[2049] = {0};
   double im[2049] = {0};


   decon_fft_forward->SetPoints(x);
   decon_fft_forward->Transform();
   decon_fft_forward->GetPointsComplex(re,im);
   
   for(int i = 606; i<1148; i++){
        //if((i<606)||(i>1147)) {ant_re[ant*2][i] = 0; ant_im[ant*2][i] = 0; continue;}
   		double magnitude = sqrt(re[i] * re[i] + im[i] * im[i]);
   		
   		//h1->SetBinContent(i-606, log10(magnitude));

     	double S = magnitude - noise_mag[ch][i];
     	double SNR = S/noise_mag[ch][i];
     	double factor = 1./(1. + 1./pow(SNR,2));
     	//cout<<noise_mag[ch][i]<<endl;
     	//cout<<pos_theta<<"  "<<pos_phi<<"   "<<gd_responds[is_high].resp_mag[pos_theta][pos_phi][ant][i]<<endl;
     	   		
   		double phase = TMath::ATan(im[i]/(re[i]+1.0e-12));
     	if(re[i]<0){phase += TMath::Pi();}
     	if(phase>TMath::Pi()) {phase -= 2*TMath::Pi();}
     	double new_mag = factor * magnitude * gd_responds[is_high].resp_mag[pos_theta][pos_phi][ant][i];

     	double new_pha = phase + gd_responds[is_high].resp_pha[pos_theta][pos_phi][ant][i];
     	new_pha = fmod(new_pha, 2*PI);
   		if(new_pha>PI) new_pha = new_pha - 2*PI;
   		if(new_pha<(-1*PI)) new_pha = new_pha + 2*PI;
   		

     	
     	// TO DO
     	re[i] = new_mag*TMath::Cos(new_pha);
     	im[i] = new_mag*TMath::Sin(new_pha);

     	     				 }  
     	     				 			
   decon_fft_backward->SetPointsComplex(re, im);
   decon_fft_backward->Transform();
   fft_back_decon_x = decon_fft_backward->GetPointsReal();

	for(int i=0;i<1500;i++){
		x[i] = fft_back_decon_x[i]/4096.;
  						}
   
        
}

void get_noise_level(double noise_mag[8][2049]){

	float h1, h2, h3, h4, v1, v2, v3, v4;
   
   FILE *fpin = fopen("/home/cyc/software/TAROGE-4_pulser/pulser/H_pol/noise_level.txt","r");// 600m
   cout<<fpin<<endl;
  for(int i=0; i<2049;i++){

   int p = fscanf(fpin,"%*f %f %f %f %f %f %f %f %f", &h1, &v1, &h2, &v2, &h3, &v3, &h4, &v4);
   noise_mag[0][i] = h1;
   noise_mag[1][i] = v1;
   noise_mag[2][i] = h2;
   noise_mag[3][i] = v2;
   noise_mag[4][i] = h3;
   noise_mag[5][i] = v3;
   noise_mag[6][i] = h4;
   noise_mag[7][i] = v4;
   
   }
}





//result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
void coherent_add(int pol, double in_theta, double in_phi, double f_x[8][1500], double added_x[15000]){	

   double R = 10000;
   
   // interpolation
   double _x[4][15000] = {0.};
   double shifted_x[4][15000] = {0.};
   //double added_x[15000] = {0.};
   
   for(int ant = 0; ant < 4; ant++){
   		int ch = 2*ant + pol;
   		
   		TSpline3 *s = new TSpline3("grs", t, f_x[ch], 1500);
   		for(int i=0; i<15000; i++){
   			_x[ant][i] = s->Eval(0.08*i);
   			}
   		delete s;
   		s = nullptr;
   		}
   	
   	// spliter for oscillator	
   	for(int i=0; i<15000; i++){
   		_x[2*3+pol][i] = sqrt(2.);
   		}   
     
   for(int ant = 0; ant < 4; ant++){
   		int ch = 2*ant + pol;    	

    	double delta_t = expectedTimeDiff_angle_with_R(in_theta, in_phi, R, coor_ant[ch], coor_ant[2*3+pol]) + (-delay_ant[2*3+pol] + delay_ant[ch])/1000.;
    	int shift_points = round(delta_t/0.08);
    	cout<<"shift_points:   "<<shift_points<<endl;

    		for(int k=0; k<15000; k++){
    			int l = k + shift_points;
    			if(l<0){l+=15000;}
    			if(l>=15000){l-=15000;}
   				shifted_x[ant][l] = _x[ant][k];
   			}
			
   }



   
   for(int ant = 0; ant < 4; ant++){
   		int ch = 2*ant + pol; 
   		for(int k=0; k<15000; k++){
   			added_x[k] += shifted_x[ant][k]/4.; 
   		}
   }

	TGraph *FW_waveform1 = new TGraph(15000, f_t, shifted_x[0]);
    //FW_waveform1->Draw("APL");
    
	TGraph *FW_waveform2 = new TGraph(15000, f_t, shifted_x[2]);
    //FW_waveform2->Draw("PL same");
    FW_waveform2->SetLineColor(2);
   
}




Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j)
         {
//cout<<"output:   "<<R<<endl;
   theta = theta * TMath::DegToRad();
   phi = phi * TMath::DegToRad();
   double drone_xyz[3] = {0};
   //R=10000000;
   drone_xyz[0] = TMath::Cos(theta)*TMath::Cos(phi)*R;
   drone_xyz[1] = TMath::Cos(theta)*TMath::Sin(phi)*R;
   drone_xyz[2] = TMath::Sin(theta)*R;
   return expectedTimeDiff_coord(drone_xyz, antenna_i, antenna_j);
         }
         
                  
         
Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j)
         {
         double sum_i = 0;
         double sum_j = 0;
    for(int k=0; k<3; k++){
        sum_i += pow((drone_xyz[k] - antenna_i[k]), 2);
        sum_j += pow((drone_xyz[k] - antenna_j[k]), 2);
        }

   double dt = (sqrt(sum_j) - sqrt(sum_i))/(Speed_Of_Light);
   //cout<<dt/1e-9<<"wtf?   "<<endl;


   // s to ns
   return dt/1.0e-9;
         }
         
         
         
         
void get_response(){
  
  // in unit of MHz
  double FreqBin = 1./((4096)*0.8e-9)/1.0e06;

  TFile* f0 = new TFile("/home/cyc/software/shih-hao/TAM/response/Taroge4/Taroge4-FEEresponse.root");
  
  TGraph* gr_response_mag[8] = {NULL};
  TGraph* gr_response_phase[8] = {NULL};
  
  //need to check the mapping order
  // caution!! mapping order T1H T2H T3H T4H T1V T2V T3 T4V
  gr_response_mag[0] = (TGraph*)f0->Get("gr_mag_0");
  gr_response_phase[0] = (TGraph*)f0->Get("gr_phase_0");

  gr_response_mag[1] = (TGraph*)f0->Get("gr_mag_4");
  gr_response_phase[1] = (TGraph*)f0->Get("gr_phase_4");

  gr_response_mag[2] = (TGraph*)f0->Get("gr_mag_1");
  gr_response_phase[2] = (TGraph*)f0->Get("gr_phase_1");
  
  gr_response_mag[3] = (TGraph*)f0->Get("gr_mag_5");
  gr_response_phase[3] = (TGraph*)f0->Get("gr_phase_5");

  gr_response_mag[4] = (TGraph*)f0->Get("gr_mag_2");
  gr_response_phase[4] = (TGraph*)f0->Get("gr_phase_2");
  
  gr_response_mag[5] = (TGraph*)f0->Get("gr_mag_6");
  gr_response_phase[5] = (TGraph*)f0->Get("gr_phase_6");

  gr_response_mag[6] = (TGraph*)f0->Get("gr_mag_3");
  gr_response_phase[6] = (TGraph*)f0->Get("gr_phase_3");

  gr_response_mag[7] = (TGraph*)f0->Get("gr_mag_7");
  gr_response_phase[7] = (TGraph*)f0->Get("gr_phase_7");
  
  //Double_t * FEE_mag = new Double_t [1000];
  //Double_t * FEE_phase = new Double_t [1000];
  //double FEE_mag[8][2049] = {0.};
  //double FEE_phase[8][2049] = {0.};
  
  for(int j=0; j<8; j++){
      for(int i=0; i<2049; i++){
      	  //From MHz to Hz
  	      double fre = FreqBin*i*1.0e6;
  	      
  	  	  double FEE_mag = gr_response_mag[j]->Eval(fre, 0, "S");
  	  	  double FEE_phase = gr_response_phase[j]->Eval(fre, 0, "S");
  	  	  //cout<<FEE_mag<<"  "<<FEE_phase<<endl;
  	  	  
  	  	  // change from dB to magnitude, sqrt is because from gain in power to voltage
          FEE_mag = sqrt(pow(10, FEE_mag/10.));
          FEE_phase = FEE_phase*TMath::DegToRad();
          
          //cout<<FEE_mag<<"   "<<FEE_phase<<endl<<endl;
          FEE_re[j][i] = FEE_mag*cos(FEE_phase);
          FEE_im[j][i] = FEE_mag*sin(FEE_phase);
          
  	  	  }
  	  
  	  }

  
  TFile* f1 = new TFile("/home/cyc/software/shih-hao/TAM/response/Taroge4/Taroge3-LPDAresponse-Vpol.root");
  TTree * LPDA_response_V = (TTree*) f1->Get("AntTree");

  TFile* f2 = new TFile("/home/cyc/software/shih-hao/TAM/response/Taroge4/Taroge3-LPDAresponse-Hpol.root");
  TTree * LPDA_response_H = (TTree*) f2->Get("AntTree");
  
  double thetas = 0;
  double phis = 0;
  // 50 MHz to 650 MHz, in unit of MHz
  double frequencies[121] = {0.};
  double Re_Hs_theta[111] = {0.};
  double Im_Hs_theta[111] = {0.};
  double Re_Hs_phi[111] = {0.};
  double Im_Hs_phi[111] = {0.};
 
  //double ANT_V_re[37][73][2049] = {0.};//theta: -90, -85, ..., 90; phi:0, 5, ..., 360
  //double ANT_V_im[37][73][2049] = {0.};
  //double ANT_H_re[37][73][2049] = {0.};//theta: -90, -85, ..., 90; phi:0, 5, ..., 360
  //double ANT_H_im[37][73][2049] = {0.};



// For V-pol
  
  LPDA_response_V->SetBranchAddress("thetas", &thetas);
  LPDA_response_V->SetBranchAddress("phis", &phis);
  LPDA_response_V->SetBranchAddress("frequencies", &frequencies);
  LPDA_response_V->SetBranchAddress("Re_Hs_theta", &Re_Hs_theta);
  LPDA_response_V->SetBranchAddress("Im_Hs_theta", &Im_Hs_theta);
  LPDA_response_V->SetBranchAddress("Re_Hs_phi", &Re_Hs_phi);
  LPDA_response_V->SetBranchAddress("Im_Hs_phi", &Im_Hs_phi);


  int event_number = LPDA_response_V->GetEntries();
  cout<<LPDA_response_V<<"  "<<event_number<<endl;

  for (int entry = 0; entry<event_number; entry++){

      LPDA_response_V->GetEntry(entry);
      cout<<entry<<endl;
      //continue;
      //if(thetas==0 && phis==0) break;
      // frequencies in unit of MHz
      TGraph * lpda_theta_re = new TGraph(111, frequencies, Re_Hs_theta);
      TGraph * lpda_theta_im = new TGraph(111, frequencies, Im_Hs_theta);
      TGraph * lpda_phi_re = new TGraph(111, frequencies, Re_Hs_phi);
      TGraph * lpda_phi_im = new TGraph(111, frequencies, Im_Hs_phi);
      // For faster evaluation
      lpda_theta_re->SetBit(true);
      lpda_theta_im->SetBit(true);
      lpda_phi_re->SetBit(true);
      lpda_phi_im->SetBit(true);
      
      int i_theta = thetas/5 + 18;//-90/5=-18
      int i_phi = phis/5;
      //if(i_theta!=15||i_phi!=5) continue;

      
      for(int i=0; i<2049; i++){
  	      double fre = FreqBin*i;
  	      
  	  	  ANT_V_theta_re[i_theta][i_phi][i] = lpda_theta_re->Eval(fre, 0, "S");
  	  	  ANT_V_theta_im[i_theta][i_phi][i] = lpda_theta_im->Eval(fre, 0, "S");
  	  	  ANT_V_phi_re[i_theta][i_phi][i] = lpda_phi_re->Eval(fre, 0, "S");
  	  	  ANT_V_phi_im[i_theta][i_phi][i] = lpda_phi_im->Eval(fre, 0, "S");
  	  	   	  
  	  	  }

    delete lpda_theta_re;
	lpda_theta_re = nullptr;
    delete lpda_theta_im;
	lpda_theta_im = nullptr;
    delete lpda_phi_re;
	lpda_phi_re = nullptr;
    delete lpda_phi_im;
	lpda_phi_im = nullptr;
    }





// For H-pol
  
  LPDA_response_H->SetBranchAddress("thetas", &thetas);
  LPDA_response_H->SetBranchAddress("phis", &phis);
  LPDA_response_H->SetBranchAddress("frequencies", &frequencies);
  LPDA_response_H->SetBranchAddress("Re_Hs_theta", &Re_Hs_theta);
  LPDA_response_H->SetBranchAddress("Im_Hs_theta", &Im_Hs_theta);
  LPDA_response_H->SetBranchAddress("Re_Hs_phi", &Re_Hs_phi);
  LPDA_response_H->SetBranchAddress("Im_Hs_phi", &Im_Hs_phi);


  event_number = LPDA_response_V->GetEntries();
  cout<<LPDA_response_H<<"  "<<event_number<<endl;

  for (int entry = 0; entry<event_number; entry++){

      LPDA_response_H->GetEntry(entry);
      cout<<entry<<endl;
      //continue;
      //if(thetas==0 && phis==0) break;
      // frequencies in unit of MHz
      TGraph *lpda_theta_re = new TGraph(111, frequencies, Re_Hs_theta);
      TGraph *lpda_theta_im = new TGraph(111, frequencies, Im_Hs_theta);
      TGraph *lpda_phi_re = new TGraph(111, frequencies, Re_Hs_phi);
      TGraph *lpda_phi_im = new TGraph(111, frequencies, Im_Hs_phi);
      // For faster evaluation
      lpda_theta_re->SetBit(true);
      lpda_theta_im->SetBit(true);
      lpda_phi_re->SetBit(true);
      lpda_phi_im->SetBit(true);
      
      int i_theta = thetas/5 + 18;//-90/5=-18
      int i_phi = phis/5;
      //if(i_theta!=15||i_phi!=5) continue;

      
      for(int i=0; i<2049; i++){
  	      double fre = FreqBin*i;
  	  	  ANT_H_theta_re[i_theta][i_phi][i] = lpda_theta_re->Eval(fre, 0, "S");
  	  	  ANT_H_theta_im[i_theta][i_phi][i] = lpda_theta_im->Eval(fre, 0, "S");
  	  	  ANT_H_phi_re[i_theta][i_phi][i] = lpda_phi_re->Eval(fre, 0, "S");
  	  	  ANT_H_phi_im[i_theta][i_phi][i] = lpda_phi_im->Eval(fre, 0, "S");
  	  	  //cout<<ANT_H_theta_re[i_theta][i_phi][i]<<endl;
  	  	  }
    delete lpda_theta_re;
	lpda_theta_re = nullptr;
    delete lpda_theta_im;
	lpda_theta_im = nullptr;
    delete lpda_phi_re;
	lpda_phi_re = nullptr;
    delete lpda_phi_im;
	lpda_phi_im = nullptr;
    }
	
}

void Quotient(double a_re, double a_im, double b_re, double b_im, double &o_re, double &o_im){
	double tmp = b_re*b_re + b_im*b_im;

		if(tmp ==0. ){	//inf
			o_re = 0.;
			o_im = 0.;
		}else{
			o_re = ( a_re*b_re + a_im*b_im ) / tmp;
			o_im = ( a_im*b_re - a_re*b_im ) / tmp;
		}
		//cout<<a_re<<"  "<<a_im<<"   "<<b_re<<"  "<<b_im<<"   "<<o_re<<"   "<<o_im<<endl;

}

//(a b)(x)=(M)-->y=(dM-bN)/(ad-bc)
//(c d)(y)=(N)-->x=(cM-aN)/(bc-ad)
void Invert_complex(TComplex a, TComplex b, TComplex c, TComplex d, TComplex M, TComplex N, TComplex &x, TComplex &y){
    //cout<<a.Re()<<"  "<<b.Re()<<"   "<<c.Im()<<"   "<<M.Im()<<endl;
	x = (d*M-b*N)/(a*d-b*c);
	y = (c*M-a*N)/(b*c-a*d);
	//cout<<"Invert_complex  "<<x.Re()<<"   "<<x.Im()<<endl;
}

void Deconvolute_antenna(double Re_ant[8][2049], double Im_ant[8][2049], double ANT_H_theta_re[2049], double ANT_H_theta_im[2049], double ANT_H_phi_re[2049], double ANT_H_phi_im[2049], double ANT_V_theta_re[2049], double ANT_V_theta_im[2049], double ANT_V_phi_re[2049], double ANT_V_phi_im[2049], double E_theta[4][4096], double E_phi[4][4096]){


   	double E_theta_re[4][2049] = {0.};
   	double E_theta_im[4][2049] = {0.};
   	double E_phi_re[4][2049] = {0.};
   	double E_phi_im[4][2049] = {0.};
   	
   	//TComplex c_R_s(-0.996671,0.00266341);
   	//TComplex c_R_p(-0.132986,-0.344515);
   		
   for(int ant=0; ant<4; ant++){

   		for(Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
   		    //1147
       		if((k<606)||(k>1148)) {continue;}
       		//if((k<327)||(k>329)) {continue;}
       		else{
       		
       		//ANT_H_theta_re[k] = 0.;
       		//ANT_H_theta_im[k] = 0.;
       		//ANT_V_phi_re[k] = 0.;
       		//ANT_V_phi_im[k] = 0.;
       		
       		E_theta_re[ant][k] = Re_ant[2*ant+1][k];
       		E_theta_im[ant][k] = Im_ant[2*ant+1][k];
       		E_phi_re[ant][k] = Re_ant[2*ant][k];
       		E_phi_im[ant][k] = Im_ant[2*ant][k];
       		
       			//cout<<ANT_H_theta_re[k]<<"    "<<ANT_H_theta_im[k]<<endl;
       			TComplex x(0, 0);// E_theta
       			TComplex y(0, 0);// E_phi
       			TComplex M(Re_ant[2*ant][k], Im_ant[2*ant][k]);// H-pol
       			TComplex N(Re_ant[2*ant+1][k], Im_ant[2*ant+1][k]);// V-pol
       			TComplex a(ANT_H_theta_re[k], ANT_H_theta_im[k]);
       			TComplex b(ANT_H_phi_re[k], ANT_H_phi_im[k]);
       			TComplex c(ANT_V_theta_re[k], ANT_V_theta_im[k]);
       			TComplex d(ANT_V_phi_re[k], ANT_V_phi_im[k]);
       			Invert_complex(a, b, c, d, M, N, x, y);
       			
       			//x = x*c_R_p;
       			//y = y*c_R_s;
       			
       			E_theta_re[ant][k] = x.Re();
       			E_theta_im[ant][k] = x.Im();
       			E_phi_re[ant][k] = y.Re();
       			E_phi_im[ant][k] = y.Im();
       			
       		}
                             		}
                             
  		filter_fft_back->SetPointsComplex(E_theta_re[ant], E_theta_im[ant]);
  		filter_fft_back->Transform();
  		time_signal = filter_fft_back->GetPointsReal();


		for(int i=0;i<4096;i++){
			E_theta[ant][i] = (time_signal[i]/4096.);
  								}
  						

  		filter_fft_back->SetPointsComplex(E_phi_re[ant], E_phi_im[ant]);
  		filter_fft_back->Transform();
  		time_signal = filter_fft_back->GetPointsReal();


		for(int i=0;i<4096;i++){
			E_phi[ant][i] = (time_signal[i]/4096.);
  								}
  		}

}




void CLEAN_deconvolution(double Re_ant[8][2049], double Im_ant[8][2049], double ANT_V_theta_re[2049], double ANT_V_theta_im[2049], double E_theta[4][4096], double E_phi[4][4096], double ws_dec[409600], int polarity, int filter){

	
	//double ws_dec[409600];
	
	double response_re[4][2049] = {0.};
	double response_im[4][2049] = {0.};

	for(int ant=0; ant<4; ant++){
		int pol = 1;
		int channel = ant*2+pol;
		
   		for(Int_t k=0; k<2049; k++){// 0 means DC component // from 185-350 MHz
   		    //1147
       		if((k<606)||(k>1148)) {continue;}
       			
       		TComplex ANT_V(ANT_V_theta_re[k], ANT_V_theta_im[k]);
       		TComplex FEE(FEE_re[channel][k], FEE_im[channel][k]);
       		TComplex response = ANT_V*FEE;//caution!!! cyc
       		
       		response_re[ant][k] = response.Re();
       		response_im[ant][k] = response.Im();
       	}
       	
       	double tdoa_best = 0.;
       	//double ws_xcor[409600];
       	double *ws_xcor = new double[4096];
		for(int i=0; i<4096; i++){ws_xcor[i] = 0; ws_dec[i]=0;}
       	double tmp_n_peak;
       	
       	double xcorPeak = FFTCrossCorrelate(Re_ant[channel], Im_ant[channel], response_re[ant], response_im[ant], ws_xcor, tdoa_best, tmp_n_peak, 0);
       	
		const double Threshold = 0.05 * fabs(xcorPeak); //#### of peak xcor coeff
		
		//double ws_hh_xcor[409600];
		double *ws_hh_xcor = new double[4096];
		for(int i=0; i<4096; i++){ws_hh_xcor[i] = 0;}
		
		double tmp_tdoa;
		double n_peak_hh;
		const double Norm_hh = FFTCrossCorrelate(response_re[ant], response_im[ant], response_re[ant], response_im[ant], ws_hh_xcor, tmp_tdoa, n_peak_hh, 1);
		
		cout<<xcorPeak<<"   "<<Norm_hh<<"    "<<tdoa_best<<"    "<<tmp_tdoa<<"    "<<tmp_n_peak<<"  "<<n_peak_hh<<endl;
		
		int step = 0, n_shift=0;
		//if(0)
		while( fabs(xcorPeak) > Threshold)
		{
			int n_peak_hy = 0;
			double max_xcorPeak = TMath::MaxElement(3000, &ws_xcor[1000]);
			double min_xcorPeak = TMath::MinElement(3000, &ws_xcor[1000]);
			if(step==0){
				if(polarity==1){xcorPeak = max_xcorPeak; n_peak_hy = TMath::LocMax(3000, &ws_xcor[1000]);}
				else if(polarity==-1){xcorPeak = min_xcorPeak; n_peak_hy = TMath::LocMin(3000, &ws_xcor[1000]);}
				else if(polarity==0){
					if(max_xcorPeak>fabs(min_xcorPeak)){
						xcorPeak = max_xcorPeak; n_peak_hy = TMath::LocMax(3000, &ws_xcor[1000]);
					}
					else{
						xcorPeak = min_xcorPeak; n_peak_hy = TMath::LocMin(3000, &ws_xcor[1000]);
					}
				}
			}
			
			else{
				if(max_xcorPeak>fabs(min_xcorPeak)){
					xcorPeak = max_xcorPeak; n_peak_hy = TMath::LocMax(3000, &ws_xcor[1000]);
				}
				else{
					xcorPeak = min_xcorPeak; n_peak_hy = TMath::LocMin(3000, &ws_xcor[1000]);
				}
			}
			
			double scale = xcorPeak / Norm_hh;
			
			n_peak_hy+=1000;
			int n_shift = n_peak_hy - n_peak_hh;
			
			const int N_min = 0;//max( 0, n_shift); //skip (n-n_shift) <0
			const int N_max = 4096;//min( 409600, 409600 + n_shift ); //the shorter one
			for(int n=N_min; n< N_max; n++)
			{
				int i = n - n_shift;
				if(i<0) i += 4096;
				if(i>=4096) i -= 4096;
				ws_xcor[n] += ((-1)*scale) * ws_hh_xcor[i];
			}
	
			if(n_shift >= 0 ){ws_dec[n_shift] += scale;} //+-
		
			else{ws_dec[ 4096 +n_shift] += scale;}
		
			//if(scale>1) break;
			
			cout<<step<<"  "<<n_peak_hy<<"   "<<n_peak_hh<<"    "<<n_shift<<"   "<<xcorPeak<<"   "<<Norm_hh<<"   "<<scale<<"   "<<ws_xcor[n_peak_hy]<<endl;
			/*
			TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
			c1->cd();
			double *t = new double[4096];
			for(int i=0; i<4096; i++){t[i] = i*0.8;}
			TGraph *cyc = new TGraph(4096, t, ws_dec);
			cyc->Draw();
			string c1_fig_name = "./photos/" + std::to_string(ant) + "_CLEAN_waveform_" + std::to_string(step) + ".png";
    		c1->SaveAs(c1_fig_name.c_str());
    		*/
    		step++;
	
		
     	}
     	
     	TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
			c1->cd();
			double *t = new double[4096];
			for(int i=0; i<4096; i++){t[i] = i*0.8;}
			TGraph *cyc = new TGraph(4096, t, ws_dec);
			cyc->Draw();
			string c1_fig_name = "./photos/" + std::to_string(ant) + "_CLEAN_waveform_" + std::to_string(step) + ".png";
    		c1->SaveAs(c1_fig_name.c_str());
     	
     	if(filter) {CLEAN_FFT(ws_dec);}
     	
     	for(int i=0; i<4096; i++){
     		E_theta[ant][i] = ws_dec[i*1];
     		//E_phi[ant][i] = ws_xcor[i*1];
     		
     	}
     	
     }
     
     


}



void max_X_correlation_delay(double delay_theta[4], double delay_phi[4], double E_theta[4][4096], double E_phi[4][4096]){
	// added on 24.0424
	double E_theta_re[4][2049] = {0.};
   	double E_theta_im[4][2049] = {0.};
   	double E_phi_re[4][2049] = {0.};
   	double E_phi_im[4][2049] = {0.};
  	
  	double *time_signal = new double[4096]{0.};
  	int sum_N = 4096;
  	TVirtualFFT *fft_sum = TVirtualFFT::FFT(1, &sum_N, "C2RBACKWARD M K");	
  	double Sum_x_theta[4] = {0};
  	double Sum_x_phi[4] = {0};
  	for(int ant = 0; ant<4; ant++){
  	
  		double max_theta = 0;
  		double max_phi = 0;
  		int max_i_theta, max_i_phi;
  		for(int i=200;i<1300;i++){
  			if(fabs(E_theta[ant][i])>max_theta){max_theta=fabs(E_theta[ant][i]); max_i_theta=i;}
  			if(fabs(E_phi[ant][i])>max_phi){max_phi=fabs(E_phi[ant][i]); max_i_phi=i;}
  		}
  		double _E_theta[4096] = {0.};
  		double _E_phi[4096] = {0.};
  		for(int i=0;i<4096;i++){
  			_E_theta[i] = ((i<max_i_theta-100)||(i>max_i_theta+100))?0.:E_theta[ant][i];
  			_E_phi[i] = ((i<max_i_phi-100)||(i>max_i_phi+100))?0.:E_phi[ant][i];
  		}
  		
  		double temp_re[2049]={0.};
  		double temp_im[2049]={0.};
  		filter_fft_forward->SetPoints(_E_theta);
   		filter_fft_forward->Transform();
   		filter_fft_forward->GetPointsComplex(temp_re,temp_im);
   		for (Int_t k=0; k<2049; k++){
   			E_theta_re[ant][k] = temp_re[k];
   			E_theta_im[ant][k] = temp_re[k];
   		}
   		
   		filter_fft_forward->SetPoints(_E_phi);
   		filter_fft_forward->Transform();
   		filter_fft_forward->GetPointsComplex(temp_re,temp_im);
  		for (Int_t k=0; k<2049; k++){
   			E_phi_re[ant][k] = temp_re[k];
   			E_phi_im[ant][k] = temp_re[k];
   		}
  		
          fft_sum->SetPointsComplex(E_theta_re[ant],E_theta_re[ant]);
          fft_sum->Transform();
          time_signal = fft_sum->GetPointsReal();
          
          for(int i = 0; i <4096; i++){
          Sum_x_theta[ant] += (time_signal[i]/(4096.))*(time_signal[i]/(4096.));
          
          }
          
          fft_sum->SetPointsComplex(E_phi_re[ant],E_phi_re[ant]);
          fft_sum->Transform();
          time_signal = fft_sum->GetPointsReal();
          
          for(int i = 0; i <4096; i++){
          Sum_x_phi[ant] += (time_signal[i]/(4096.))*(time_signal[i]/(4096.));
          
          }
       }
    //delete []time_signal;
    //delete fft_sum;
    
    double ** correlation_pair_theta = new_Array2D<double>(6, 409600);
    double ** correlation_pair_phi = new_Array2D<double>(6, 409600);
    double *Correlator_time = new double[409600]{0.};
    double Correlator_re[204801] = {0};            
    double Correlator_im[204801] = {0};
    int N = 409600;
  	TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");	
   	int i = 0;
    int j = 1;

   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 0; j = 1;
        break;
      case 1:
        i = 0; j = 2;
        break;
      case 2:
        i = 0; j = 3;
        break;
      case 3:
        i = 1; j = 2;
        break;
      case 4:
        i = 1; j = 3;
        break;
      case 5:
        i = 2; j = 3;
        break;

    }


        for (Int_t k=1; k<2049; k++){//0 means DC component
            Correlator_re[k] = (E_theta_re[i][k] * E_theta_re[j][k] + E_theta_im[i][k] * E_theta_im[j][k]);               
            Correlator_im[k] = (E_theta_re[i][k] * E_theta_im[j][k] - E_theta_im[i][k] * E_theta_re[j][k]);
        }//k  

          fft_back->SetPointsComplex(Correlator_re,Correlator_im);
          fft_back->Transform();
          Correlator_time = fft_back->GetPointsReal();

         double factor = sqrt(Sum_x_theta[i])*sqrt(Sum_x_theta[j])*4096.;
         double max = 0;
         int mm;
                       for (Int_t k=0; k<409600; k++){
                        //cout<<Correlator_time[k]<<endl;
                correlation_pair_theta[pair][k] = Correlator_time[k]/factor;
                if (correlation_pair_theta[pair][k] > max){max = correlation_pair_theta[pair][k]; mm = k;}
                                                    }//k 
		//cout<<"mm:  "<<mm<<"   ";
		if(mm>204799) {mm -= 409600;}
		cout<<"mm:  "<<mm<<endl;
     }// pair loop




   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 0; j = 1;
        break;
      case 1:
        i = 0; j = 2;
        break;
      case 2:
        i = 0; j = 3;
        break;
      case 3:
        i = 1; j = 2;
        break;
      case 4:
        i = 1; j = 3;
        break;
      case 5:
        i = 2; j = 3;
        break;

    }


        

        for (Int_t k=1; k<2049; k++){//0 means DC component
            Correlator_re[k] = (E_phi_re[i][k] * E_phi_re[j][k] + E_phi_im[i][k] * E_phi_im[j][k]);               
            Correlator_im[k] = (E_phi_re[i][k] * E_phi_im[j][k] - E_phi_im[i][k] * E_phi_re[j][k]);
        }//k  

          fft_back->SetPointsComplex(Correlator_re,Correlator_im);
          fft_back->Transform();
          Correlator_time = fft_back->GetPointsReal();

         double factor = sqrt(Sum_x_phi[i])*sqrt(Sum_x_phi[j])*4096.;
         double max = 0;
         int mm;
                       for (Int_t k=0; k<409600; k++){
                        //cout<<Correlator_time[k]<<endl;
                correlation_pair_phi[pair][k] = Correlator_time[k]/factor;
                if (correlation_pair_phi[pair][k] > max){max = correlation_pair_phi[pair][k]; mm = k;}
                                                    }//k 
		//cout<<"mm:  "<<mm<<"   ";
		if(mm>204799) {mm -= 409600;}
		cout<<"mm:  "<<mm<<endl;

     }// pair loop
     
     //delete[] Correlator_time;
     //delete fft_back; 
     //delete[] correlation_pair_theta;
     //delete[] correlation_pair_phi;
     
     //result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
     //double delay_theta[4] = {0.};
     //double delay_phi[4] = {0.};
     double pair_delay[6] = {0.};
     double max_theta = 0;
     double max_phi = 0;
     for(int i=-1000; i<1000; i+=10){
     	for(int j=-1000; j<1000; j+=10){
     		for(int k=-1000; k<1000; k+=10){
     			pair_delay[0] = i;
     			pair_delay[1] = j;
     			pair_delay[2] = k;
     			pair_delay[3] = j-i;
     			pair_delay[4] = k-i;
     			pair_delay[5] = k-j;
     			double sum_theta = 0;
     			double sum_phi = 0;
     			for(int pair=0; pair<6; pair++){
     				int m = (-1)*pair_delay[pair];
     				if(m<0) m+= 409600;
     				sum_theta += correlation_pair_theta[pair][m];
     				sum_phi += correlation_pair_phi[pair][m];
     			}
     			if(sum_theta>max_theta){
     			max_theta=sum_theta;delay_theta[1]=i;delay_theta[2]=j;delay_theta[3]=k;
     			}
     			if(sum_phi>max_phi){
     			max_phi=sum_phi;delay_phi[1]=i;delay_phi[2]=j;delay_phi[3]=k;
     			}
     		}
     	}
     }
     
     for(int ant=0; ant<4; ant++){
     	cout<<delay_theta[ant]<<endl;
     	cout<<delay_phi[ant]<<endl;
     }
     
     

	// added on 24.0424
}


        
         
//result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
// i: Hpol; j: Vpol
double get_delay(){	






        double Correlator_re[204801] = {0};            
        double Correlator_im[204801] = {0};

                  	   for (Int_t k=1; k<2049; k++){//0 means DC component
            Correlator_re[k] = (H_ant4_re[k] * V_ant4_re[k] + H_ant4_im[k] * V_ant4_im[k]); 
                   
            Correlator_im[k] = (H_ant4_re[k] * V_ant4_im[k] - H_ant4_im[k] * V_ant4_re[k]);
            //cout<<H_ant4_re[k]<<endl;
                
                                                    }//k  



          fft_back->SetPointsComplex(Correlator_re,Correlator_im);
          fft_back->Transform();
          Correlator_time = fft_back->GetPointsReal();
          

         double max = 0;
         int mm;
                       for (Int_t k=0; k<409600; k++){
                if (Correlator_time[k] > max){max = Correlator_time[k]; mm = k;}
                //cout<<Correlator_time<<endl<<endl;

                                                    }//k 

if(mm>204799) {mm -= 409600;}


double delay_time = mm/409600.*0.8*4096.;
double correlation_value = max;

cout<<"HV delay: "<<delay_time<<"   "<<correlation_value<<endl;
return delay_time;
}




//result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
// i: Hpol; j: Vpol
double FFTCrossCorrelate(double ws1_re[2049], double ws1_im[2049], double ws2_re[2049], double ws2_im[2049], double ws_xcor[409600], double &tdoa_best, double &n_peak, int polarity){	




        double Correlator_re[2049] = {0};            
        double Correlator_im[2049] = {0};
        

		
		/*
        for (Int_t k=1; k<2049; k++){//0 means DC component
            Correlator_re[k] = ws1_re[k]; 
                   
            Correlator_im[k] = ws1_im[k];
            //cout<<H_ant4_re[k]<<endl;
                                                    }//k  

        fft_back->SetPointsComplex(Correlator_re,Correlator_im);
        fft_back->Transform();
        Correlator_time = fft_back->GetPointsReal();
        
        double sum_1 = 0;
        for(Int_t k=0; k<4096; k++){
        	sum_1 += pow(Correlator_time[k]/4096.,2);
        	//ws_xcor[k] = Correlator_time[k]/4096.;
        }
                
		

        for (Int_t k=1; k<2049; k++){//0 means DC component
            Correlator_re[k] = ws2_re[k]; 
                   
            Correlator_im[k] = ws2_im[k];
            //cout<<H_ant4_re[k]<<endl;
                                                    }//k  

        fft_back->SetPointsComplex(Correlator_re,Correlator_im);
        fft_back->Transform();
        Correlator_time = fft_back->GetPointsReal();
        
        double sum_2 = 0;
        for(Int_t k=0; k<4096; k++){
        	sum_2 += pow(Correlator_time[k]/4096.,2);
        }
        */
        
        double sum_1=1, sum_2=1;
        
                
        

        for (Int_t k=1; k<2049; k++){//0 means DC component
            Correlator_re[k] = (ws1_re[k] * ws2_re[k] + ws1_im[k] * ws2_im[k]); 
                   
            Correlator_im[k] = (ws1_re[k] * ws2_im[k] - ws1_im[k] * ws2_re[k]);
            //cout<<H_ant4_re[k]<<endl;
                
                                                    }//k  



        fft_back->SetPointsComplex(Correlator_re,Correlator_im);
        fft_back->Transform();
        Correlator_time = fft_back->GetPointsReal();
          

         double peak = 0;
         int mm;
        
        for(Int_t k=0; k<4096; k++){
        	int i = 2047+k;
        	if(i>=4096){i-=4096;}
        	ws_xcor[i] = Correlator_time[k]/(sqrt(sum_1*sum_2)*4096.);
        }
        
        
        for(Int_t k=0; k<4096; k++){	
        	if(polarity==1){
            	if (ws_xcor[k] > peak){peak = ws_xcor[k]; mm = k;}
            }
            if(polarity==0){
            	if (fabs(ws_xcor[k]) > fabs(peak)){peak = ws_xcor[k]; mm = k;}
            }
            if(polarity==-1){
            	if (ws_xcor[k] < peak){peak = ws_xcor[k]; mm = k;}
            }
        }//k 
		
n_peak=mm;
if(mm>2047) {mm -= 4096;}


tdoa_best = mm/4096.*0.8*4096.;
double correlation_value = peak;

return correlation_value;

}

void CLEAN_FFT(double x[409600]){
	int N = 4096;
	TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C ES K");
	double *re = new double[4096];
	double *im = new double[4096];
	for(int i=0; i<4096; i++){re[i] = 0; im[i] = 0;}
	fft_forward->SetPoints(x);
   	fft_forward->Transform();
   	fft_forward->GetPointsComplex(re,im);
   	
   	for (Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
       if((k<606)||(k>1147)) {re[k] = 0.; im[k] = 0.;}
    }
    
    fft_back->SetPointsComplex(re,im);
    fft_back->Transform();
    Correlator_time = fft_back->GetPointsReal();
    
    for(Int_t k=0; k<4096; k++){
        x[k] = Correlator_time[k]/4096.;
    }
	
}


void CLEAN_simulation(double x[4096], double E_field_x[4096], double ANT_V_theta_re[2049], double ANT_V_theta_im[2049], int channel){
	int N = 4096;
	TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C ES K");
	double *re = new double[4096];
	double *im = new double[4096];
	for(int i=0; i<4096; i++){re[i] = 0; im[i] = 0;}
	fft_forward->SetPoints(x);
   	fft_forward->Transform();
   	fft_forward->GetPointsComplex(re,im);

   	for (Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
       if((k<606)||(k>1147)) {re[k] = 0.; im[k] = 0.;}
    }
    
    fft_back->SetPointsComplex(re,im);
    fft_back->Transform();
    Correlator_time = fft_back->GetPointsReal();
    
    for(Int_t k=0; k<4096; k++){
        E_field_x[k] = Correlator_time[k]/4096.;
    }
  
  
  double response_re[2049]={0};
  double response_im[2049]={0};
    
    
    
   	
   	for (Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
   	
   		TComplex ANT_V(ANT_V_theta_re[k], ANT_V_theta_im[k]);
  		TComplex FEE(FEE_re[channel][k], FEE_im[channel][k]);
  		TComplex response = FEE*ANT_V;
  		response_re[k] = response.Re();
  		response_im[k] = response.Im();
   	
       if((k<606)||(k>1147)) {re[k] = 0.; im[k] = 0.;}
       else{
       		TComplex x(re[k], im[k]);
       		cout<<TComplex::Abs(x)<<endl;
       		TComplex y(response_re[k], response_im[k]);
       		TComplex z = x*y;
       		re[k] = z.Re();
       		im[k] = z.Im();
       }
 
    }
    
    fft_back->SetPointsComplex(re,im);
    fft_back->Transform();
    Correlator_time = fft_back->GetPointsReal();
    
    TRandom3 *r3;
	r3 = new TRandom3();
    for(Int_t k=0; k<1500; k++){
        x[k] = Correlator_time[k]/4096. + r3->Gaus(0,10);
    }
	
}

void Adjust_delay(double Re_ant[8][2049], double Im_ant[8][2049], double H_V_delay, double in_theta, double in_phi){

	// MHz
	double FreqBin = 1./((4096)*0.8e-9)/1.0e06;
	


	// H-pol antenna 1 T1H as the standar
	for(int channel=1; channel<8; channel++){
		
		// ns
		double delta_t = expectedTimeDiff_angle_with_R(in_theta, in_phi, 100000, coor_ant[0], coor_ant[channel]);
		cout<<delta_t<<"   "<<H_V_delay<<endl;
		//V-pol
		//caution: should be + 
		if(channel%2==1){delta_t += H_V_delay;}
		
		for(int k=0; k<2049; k++){
			double mag = sqrt(Re_ant[channel][k]*Re_ant[channel][k] + Im_ant[channel][k]*Im_ant[channel][k]);
			double phase = TMath::ATan2(Im_ant[channel][k], Re_ant[channel][k]);		

			
			// MHz to Hz
			double delay_phase = 2*TMath::Pi()*k*FreqBin*1.0e06*delta_t*1.0e-9;
			phase = phase + delay_phase;
			Re_ant[channel][k] = mag*TMath::Cos(phase);
			Im_ant[channel][k] = mag*TMath::Sin(phase);
		}
	}
	
	
	/*
	for(int ant=0; ant<4; ant++){
		for(int k=0; k<2049; k++){
			//V-pol
			double mag = sqrt(Re_ant[2*ant+1][k]*Re_ant[2*ant+1][k] + Im_ant[2*ant+1][k]*Im_ant[2*ant+1][k]);
			double phase = TMath::ATan2(Im_ant[2*ant+1][k], Re_ant[2*ant+1][k]);
			// MHz to Hz
			double delay_phase = 2*TMath::Pi()*k*FreqBin*1.0e06*H_V_delay*1.0e-9;
			phase = phase + delay_phase;
			Re_ant[2*ant+1][k] = mag*TMath::Cos(phase);
			Im_ant[2*ant+1][k] = mag*TMath::Sin(phase);
		}
	}
	*/

}
         
         
void get_candidate_file(string path, vector<string> &filenames){
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(path.c_str()))){
        cout<<"Folder doesn't Exist!"<<endl;
        return;
    }
    while((ptr = readdir(pDir))!=0) {
    //cout<<"1  "<<std::strstr(ptr->d_name, ".")<<endl;//cannot cout, will cause system crazy because cannot cout<<NULL<<endl;
       if (std::strstr(ptr->d_name, ".txt")&&!std::strstr(ptr->d_name, ".txt~")){
            filenames.push_back(path + "/" + ptr->d_name);
    }
    }
    closedir(pDir);
}
