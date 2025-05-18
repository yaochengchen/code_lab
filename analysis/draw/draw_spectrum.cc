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


int get_FFT(int ch, double x[1500], double FW_x[1500], double filtered_x[1500], double FW_re[751], double FW_im[751], double deep_filtered_x[1500]);
int filter_time(double filtered_x[1500], double SW_x[1500], int Max_i);
//HALF_WINDOW_WIDTH
void SW_fft(double SW_x[1500], int Max_i, double SW_re[751], double SW_im[751]);

void Planck_taper_windowing(double epsilon, int N, double x[], double w_x[]);



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



void draw_spectrum(int date, int run_num, int event_num){


 //get_noise_level(noise_mag);

 //read_responses_high();



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
	
	Double_t t[4096]={0};
 	Double_t x[8][4096]={0};//8 antennas
	Double_t new_x[8][4096]={0};//8 antennas
 	//Double_t frequency[751]={0};//8 antennas
	//Double_t powspm[8][751]={0};//8 antennas
	


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

// CYC!! ONLY SEE THE FIRST PULSE
  for (int i=615;i<1500;i++)
	{
	  if(T4H[i]>32512) cout<<T4H[i]<<endl;
	  t[i]=i*0.8;
	  x[0][i]= (float)T1H[i]*500/32512.;
	  x[2][i]= (float)T2H[i]*500/32512.;
	  x[4][i]= (float)T3H[i]*500/32512.;
	  x[6][i]= sqrt(2)*((float)T4H[i]*500/32512.);
	  //int k = i+3;
	  //if(k>=1500){k-=1500;}
	  x[1][i]= (float)T1V[i]*500/32512.;
	  x[3][i]= (float)T2V[i]*500/32512.;
	  x[5][i]= (float)T3V[i]*500/32512.;
	  x[7][i]= sqrt(2)*((float)T4V[i]*500/32512.);
	  //if(x[6][i]>500) cout<<x[6][i]<<endl;	  
	}
	
	TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
	c1->Divide(2, 4);
	TCanvas *c3 = new TCanvas("c3","Canvas Example",1500,1000);
	c3->Divide(2, 4);
	TCanvas *c2 = new TCanvas("c2","Canvas Example",1500,1000);
	c2->Divide(2, 4);
	for(int channel = 0; channel<=7; channel++){
	
	
	double FW_x[8][4096] = {0.};
	double filtered_x[8][4096] = {0.};
	double deep_filtered_x[8][4096] = {0.};
	double FW_re[8][2049] = {0.};
	double FW_im[8][2049] = {0.};
	double SW_x[8][4096] = {0.};
	
	//HALF_WINDOW_WIDTH
	double SW_re[8][2049] = {0.};
	double SW_im[8][2049] = {0.};
	
	int Max_i = get_FFT(channel, x[channel], FW_x[channel], filtered_x[channel], FW_re[channel], FW_im[channel], deep_filtered_x[channel]);
	
	double FW_power_spectrum[8][2049] = {0.};
    double FW_frequency[2049] = {0.};
    double T = 4096*0.8*1.0e-9;
    double n_size = 4096.;
    for(int k=0; k<2049; k++){
    	FW_frequency[k] = k*625./2048.;
    	double power = 2.*(T*1.0e06)*(FW_re[channel][k]*FW_re[channel][k]+FW_im[channel][k]*FW_im[channel][k])/(n_size*n_size*50.*1000.);
    	if(power<=0){power = 1.0e-12;}
    	FW_power_spectrum[channel][k] = 10*log10(power);
    	}
    	
	c2->cd(channel+1);
    TGraph *FW_spectrum = new TGraph(2049, FW_frequency, FW_power_spectrum[channel]);
    FW_spectrum->Draw("APL");
    FW_spectrum->GetYaxis()->SetRangeUser(-80,-10);
    FW_spectrum->GetXaxis()->SetRangeUser(150,400);

    FW_spectrum->GetXaxis()->SetTitle("MHz");
    FW_spectrum->GetYaxis()->SetTitle("dBm/MHz");	
    	
    	
	
	}
	
	
delete Tree_Muon;
file->Close();
delete file;


}




int filter_N = 4096;
Int_t filter_n_size = filter_N;
TVirtualFFT *filter_fft_forward = TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
TVirtualFFT *filter_fft_back = TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");
Double_t * fft_back_x = new Double_t[4096];

int get_FFT(int ch, double x[4096], double FW_x[4096], double filtered_x[4096], double FW_re[2049], double FW_im[2049], double deep_filtered_x[4096]){


   filter_fft_forward->SetPoints(x);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(FW_re, FW_im);
   
   FW_re[0] = 0.;
   FW_im[0] = 0.;
   filter_fft_back->SetPointsComplex(FW_re, FW_im);
   filter_fft_back->Transform();
   fft_back_x = filter_fft_back->GetPointsReal();

	for(int i=0;i<4096;i++){
		x[i] = fft_back_x[i]/4096.;
		//cout<<x[i]<<endl;
  						}
  						
  						
  						   
	double data[4096] = {0.};
	//for(int i=0; i<1500; i++){
	//	FW_x[i] = (0.42-0.5*cos(2*PI*(i)/1500.)+0.08*cos(4*PI*(i)/1500.))*x[i];
	//}
	//Planck_taper_windowing(0.2, 1500, x, data);
	for(int i=0; i<4096; i++){data[i]=x[i];}
	

	
	for(int i=0; i<4096; i++){
		FW_x[i] = data[i];
		//cout<<data[i]<<endl;
	}	
	

   filter_fft_forward->SetPoints(FW_x);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(FW_re, FW_im);

	return 1;
   
   
}




int filter_time(double filtered_x[1500], double SW_x[1500], int Max_i){
    //cout<<Max_i<<endl;
	/*
    
           for (int i=0;i<1500;i++)
	{

   	    double factor = (0.42-0.5*cos(2*PI*(i)/1500.)+0.08*cos(4*PI*(i)/1500.));
		if(factor<0.1){factor=0.1;}
		factor = 1;
		
   	       SW_x[i] = ((i<Max_i-HALF_WINDOW_SPACING)||(i>Max_i+HALF_WINDOW_SPACING))?0:(0.42-0.5*cos(2*PI*(i-Max_i+HALF_WINDOW_SPACING)/WINDOW_SPACING)+0.08*cos(2*PI*(i-Max_i+HALF_WINDOW_SPACING)/HALF_WINDOW_SPACING))*filtered_x[i]/factor;
   	       //SW_x[i] = ((i<Max_i-HALF_WINDOW_SPACING)||(i>Max_i+HALF_WINDOW_SPACING))?0:filtered_x[i];
 	                               
	}
	*/
	Planck_taper_windowing(0.1, 150, &filtered_x[Max_i-HALF_WINDOW_SPACING], &SW_x[Max_i-HALF_WINDOW_SPACING]);

return Max_i;

}



   //
   Double_t * sXcor_time = new Double_t[WINDOW_SPACING];

   int sXcor_N = WINDOW_SPACING;
   Int_t sXcor_size = sXcor_N;		
   TVirtualFFT *sXcor_fft_forward = TVirtualFFT::FFT(1, &sXcor_size,"R2C ES K");
   
   TVirtualFFT *sXcor_fft_backward = TVirtualFFT::FFT(1, &sXcor_N, "C2RBACKWARD M K");

//HALF_WINDOW_WIDTH
void SW_fft(double SW_x[1500], int Max_i, double SW_re[751], double SW_im[751]){
   sXcor_fft_forward->SetPoints(&SW_x[Max_i-HALF_WINDOW_SPACING]);//
   sXcor_fft_forward->Transform();
   sXcor_fft_forward->GetPointsComplex(SW_re, SW_im); 
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
