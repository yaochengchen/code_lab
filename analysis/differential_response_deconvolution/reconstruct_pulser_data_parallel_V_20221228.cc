//20240417 change x[4096] new_x[4096] to 4097; since filter_n_size & decon_n_size is 4097, if not change, some times wierd number after FFT will appears. It is betters to change TAROGE-4_pulser folder: reconstruct_pulser_data_parallel.cc, store_response_high.cc, store_response.cc, pulser_data_first_reconstruction.cc, get_noise_spectrum.cc, change n_size to be n instead of n+1.;

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
#include <TBits.h>
//#include <thread>
#include <unistd.h>
using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8
#define N_Threads 10


template <typename T>
T** new_Array2D(int row, int col);

template <typename T> 
void delete_Array2D(T **arr, int row, int col);

template <typename T>
T*** new_Array3D(int height, int row, int col);

template <typename T> 
void delete_Array3D(T ***p, int height, int row, int col);
 

void thread_function(int n_thread, int is_high, int nt_j, double R, int lower_thetapos, int upper_thetapos, int lower_phipos, int upper_phipos);

void deconvolute_response(int thread, int is_high, int pos_theta, int pos_phi);
 
 
 double get_delay_from_phase(int thread, int ant_i, int ant_j, double delay, TH1F *h1);
 
 
 
void get_cross_correlation(int thread, int nt_j);
int filter_time(double f_x[4096], double new_x[4097], int Max_i);
int filter_fre(int thread, double x[4097], double f_x[4096]);
Double_t expectedTimeDiff_angle(double theta, double phi, double *xyz);
Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j);
Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j);

void get_noise_level(double noise_mag[8][2049]);

void read_responses_low();
void read_responses_high();
bool SetDataTreeBranches(TTree* tr);

void low_elevation_reconstruction(int n_thread, int nt_j, double R, int lower_thetapos, int upper_thetapos, int lower_phipos, int upper_phipos, double max_theta[3], double max_phi[3], double max[3]);
void high_elevation_reconstruction(int n_thread, int nt_j, double R, int lower_thetapos, int upper_thetapos, int lower_phipos, int upper_phipos, double max_theta[3], double max_phi[3], double max[3]);




//std::thread threads[N_Threads];

double ** Correlator_time = new_Array2D<double>(N_Threads, 409600);
double ** time_signal = new_Array2D<double>(N_Threads, 4096);
int N = 409600;
TVirtualFFT *fft_back[N_Threads] = {NULL};// TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
int sum_N = 4096;
TVirtualFFT *fft_sum[N_Threads] = {NULL};// TVirtualFFT::FFT(1, &sum_N, "C2RBACKWARD M K");

int filter_N = 4096;
Int_t filter_n_size = filter_N+1;
TVirtualFFT *filter_fft_forward[N_Threads] = {NULL};// TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
TVirtualFFT *filter_fft_back[N_Threads] = {NULL};// TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");
Double_t ** filtered_x = new_Array2D<double>(N_Threads, 4096);

int decon_N = 4096;
Int_t decon_n_size = decon_N+1;		
TVirtualFFT *decon_fft_forward[N_Threads] = {NULL};// TVirtualFFT::FFT(1, &decon_n_size,"R2C ES K");
   
//********************************************
//base center   2021 photogrametry  Tower-4 is back side
double coor_ant[8][3] = {{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949},{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949}};

  double delay_ant[8] = {0,53.2,84.6,-1523.1,-92.4,-30.3,59.5,-743.1};
  
  TBits * gTrigBits = NULL;
  TBits * gOverVoltBits = NULL;
  //TBits gTrigBits(64);
  //TBits gOverVoltBits(8);
  
  int gEventNum;
  
  double ave_cro_H;
  double ave_y_peak;
  bool bGood;
  double eventTime = 0;
  Int_t timeStamp_FPGA;

  double Attenuation = 0;
  double drn_azi;
  double drn_ele;
  double rel_x;
  double rel_y;
  double rel_z;
  
  double rough_X_cor;
  double rough_theta;
  double rough_phi;

double constructed_theta = 0;
double constructed_phi = 0;
double constructed_X_cor = 0;

double tmpt_constructed_theta[50000]={0};
double tmpt_constructed_phi[50000]={0};
double tmpt_constructed_X_cor[50000]={0};

double cal_timestamp[50000]={0};
double cal_event[50000]={0};
int cal_Event_N = 0;

int flag_check[50000] = {0};

  Double_t t[1500];
  Short_t T1H[1500];
  Short_t T2H[1500];
  Short_t T3H[1500];
  Short_t T4H[1500];
  Short_t T1V[1500];
  Short_t T2V[1500];
  Short_t T3V[1500];
  Short_t T4V[1500];
 //******************************************************
 
 double drone_azi[50000] = {0};
 double drone_ele[50000] = {0}; 
 
 Double_t noise_mag[8][2049] = {0};


double ***x = new_Array3D<double>(N_Threads, 8, 4097);
 Double_t ***new_x = new_Array3D<Double_t>(N_Threads, 8, 4097);
 Double_t ***f_x = new_Array3D<Double_t>(N_Threads, 8, 4096);

double *** ant_re = new_Array3D<double>(N_Threads, 8, 2049);
double *** ant_im = new_Array3D<double>(N_Threads, 8, 2049);
  
  //double resp_mag[141][143][4][2049] = {0};
  //double resp_pha[141][143][4][2049] = {0};
  //double resp_mag_high[43][15][4][2049] = {0};
  //double resp_pha_high[43][15][4][2049] = {0};

struct Ground_Respond{
  vector<vector<vector<vector<double>>>> resp_mag;
  vector<vector<vector<vector<double>>>> resp_pha;
};
// 0 for low elevation, 1 for high elevation
Ground_Respond gd_responds[2] = {{vector<vector<vector<vector<double>>>> (141, vector<vector<vector<double>>>(181, vector<vector<double>>(4, vector<double>(2049, 0)))), vector<vector<vector<vector<double>>>> (141, vector<vector<vector<double>>>(181, vector<vector<double>>(4, vector<double>(2049, 0))))}, {vector<vector<vector<vector<double>>>> (36, vector<vector<vector<double>>>(18, vector<vector<double>>(4, vector<double>(2049, 0)))), vector<vector<vector<vector<double>>>> (36, vector<vector<vector<double>>>(18, vector<vector<double>>(4, vector<double>(2049, 0))))}};
    
double *** cyc_mag = new_Array3D<double>(N_Threads, 4, 2049);
double *** cyc_pha = new_Array3D<double>(N_Threads, 4, 2049);



  


 double *** correlation_pair = new_Array3D<double>(N_Threads, 6, 409600);	
 //double correlation_pair[6][409600] = {0};
 
 

//double ** Xcor_tdoa = new_Array2D<double>(6, 50000);
//double ** Xcor_max = new_Array2D<double>(50000, 6);

//double Expect_Xcor_tdoa[6][10000]={0};



void initial(){

//coordinates
//x,y,z: east, north, height; unit: m
//feed center  photographmetric
//  double coor_ant[8][3] = {{-4.979,-9.003,7.810},{0,0,6.218},{4.547,8.170,4.381},{16.970,4.388,3.282},{-4.979,-9.003,7.810},{0,0,6.218},{4.547,8.170,4.381},{16.970,4.388,3.282}};

//base center  photographmetric
//  double coor_ant[8][3] = {{-4.988,-9.007,1.61},{0,0,0},{4.560,8.171,-1.855},{16.994,4.372,-2.967},{-4.988,-9.007,1.61},{0,0,0},{4.560,8.171,-1.855},{16.994,4.372,-2.967}};


//base center  total station
//double coor_ant[8][3] = {{-18.365,-17.390,5.783},{-13.428,-8.345,4.262},{-8.927,-0.146,2.477},{3.486,-3.835,1.365},{-18.365,-17.390,5.783},{-13.428,-8.345,4.262},{-8.927,-0.146,2.477},{3.486,-3.835,1.365}};


double t4point_phi = 90 - 117.09;
double t4point_theta = -0.08;
double base_thick = 0.01;

//transfer t4 back side to front side
coor_ant[3][0] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Cos(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[3][1] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Sin(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[3][2] += TMath::Sin(t4point_theta*TMath::DegToRad()) * base_thick;

coor_ant[7][0] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Cos(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[7][1] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Sin(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[7][2] += TMath::Sin(t4point_theta*TMath::DegToRad()) * base_thick;



//set tower-2 to be origin
double set_t2_origin_photo[3] = {328919.115, 2697635.05, 707.86};//same as tower2
double set_t2_origin[3] = {-13.428,-8.345,4.262};//same as tower2
for(int i=0; i<8; i++){
for(int j=0; j<3; j++){
coor_ant[i][j] = coor_ant[i][j] - set_t2_origin_photo[j];// + set_t2_origin[j];
cout<<coor_ant[i][j]<<"  ";
}
cout<<endl;
}


// Initialization of fft 
for(int thread=0; thread<N_Threads; thread++){
fft_back[thread] = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
fft_sum[thread] = TVirtualFFT::FFT(1, &sum_N, "C2RBACKWARD M K");

filter_fft_forward[thread] = TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
filter_fft_back[thread] = TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");
	
decon_fft_forward[thread] = TVirtualFFT::FFT(1, &decon_n_size,"R2C ES K");
}



 get_noise_level(noise_mag);

 //cout<<"here1"<<endl;

 read_responses_low();
 read_responses_high();


}


 
//int main(int argc, char *argv[]) {
 void reconstruct_pulser_data(string fname_evt, int fnum_evt, double distance_X_max, double rough_theta, double rough_phi){
 cout<<fnum_evt<<endl;



double event_time[50000] = {0};



//double relative_coord[50000][3] = {0};
//double ** relative_coord = new_Array2D<double>(50000, 3);


 
 
 
//double Sum_template = 0;
//double Re_template[2049] = {0};
//double Im_template[2049] = {0};





  Int_t event_number;




double timestamp[50000]={0};

double event[50000]={0};



 //cout<<"here2"<<endl;
 
    TFile *file = new TFile(fname_evt.c_str(), "read");
    TTree *tr_event = (TTree*) file->Get("t");
//cout<<44444<<endl;	
SetDataTreeBranches(tr_event);

//cout<<55555<<endl;
		  
 event_number = tr_event->GetEntries();
 //cout<<"input:   "<<1<<endl;
 
 int n_thread = 0;
 
for(int nt_j=fnum_evt; nt_j<fnum_evt+1; nt_j++){
//cout<<"input:   "<<nt_j<<endl;
      
      tr_event->GetEntry(nt_j);
      
      //if(!bGood){continue;}
      //if(rough_X_cor<4.7){continue;}
      





     event_time[nt_j] = eventTime;
     

     
	 double R = distance_X_max;
     //cout<<"R:  "<<R<<endl;


     int upper_phipos = 0;
     int lower_phipos = 0;

     int upper_thetapos = 0;
     int lower_thetapos = 0;
     
     rough_phi -= 26; //change back to T4 coordinates (phi=0 point to east)

for(int high_low=0; high_low<2; high_low++){          
     bool good_event = true;
     int is_high = 0;
     // low elevation events
     //if(rough_theta<6.6){
     if(high_low==0){
     	cout<<"low decovolution: "<<endl;
     	int m_phipos = (int) round((rough_phi - (-25.))/0.6) + 90;
     	int m_thetapos = (int) round(rough_theta*10) + 70;
     	
     	upper_phipos = min(181, m_phipos + 9);//[,)
     	lower_phipos = max(0, m_phipos - 8);
     	
     	upper_thetapos = min(141, m_thetapos + 51);
     	lower_thetapos = max(0, m_thetapos - 50);
     	     
     	if((R<700.)||(m_phipos<0)||(m_phipos>=181)||(m_thetapos<0)||(m_thetapos>=141)||(!bGood)){
     		good_event = false;
     	}
     }
     
     //high elevation events
     else{
     	cout<<"high decovolution: "<<endl;
     	is_high = 1;
     	int m_phipos = (int) round((rough_phi - (-25.) - 3.)/6.) + 9;
     	int m_thetapos = (int) round(rough_theta) - 7;
     	
     	upper_phipos = 17;//min(18, m_phipos + 3);
     	lower_phipos = 0;//max(0, m_phipos - 2);
     	
     	upper_thetapos = 35;//min(36, m_thetapos + 6);
     	lower_thetapos = 0;//max(0, m_thetapos - 5);    
     	
     	if((R<700.)||(m_phipos<0)||(m_phipos>=18)||(m_thetapos<0)||(m_thetapos>=36)||(!bGood)){
     		good_event = false;
     	}
     }
     

     

     
  for (int i=0;i<1500;i++)
	{
	  //t[i]=i*0.8;
	  x[n_thread][0][i]= (float)T1H[i]*500/32512;
	  x[n_thread][2][i]= (float)T2H[i]*500/32512;
	  x[n_thread][4][i]= (float)T3H[i]*500/32512;
	  x[n_thread][6][i]= (float)T4H[i]*500/32512;
	  x[n_thread][1][i]= (float)T1V[i]*500/32512;
	  x[n_thread][3][i]= (float)T2V[i]*500/32512;
	  x[n_thread][5][i]= (float)T3V[i]*500/32512;
	  x[n_thread][7][i]= (float)T4V[i]*500/32512;	  	  
	}      		
//threads[n_thread] = std::thread(thread_function, n_thread, is_high, nt_j, R, lower_thetapos, upper_thetapos, lower_phipos, upper_phipos);

thread_function(n_thread, is_high, nt_j, R, lower_thetapos, upper_thetapos, lower_phipos, upper_phipos);

n_thread ++;

if(n_thread>=N_Threads){
	for(int i = 0; i<N_Threads; i++){
		//threads[i].join();
		}
	n_thread = 0;
}

}// high low loop

}// number loop


file->Close();   
    
              }// end of main function
              

     

void thread_function(int n_thread, int is_high, int nt_j, double R, int lower_thetapos, int upper_thetapos, int lower_phipos, int upper_phipos){              
              
 	 double max[3] = {0};
     double max_theta[3] = {0};
     double max_phi[3] = {0};   		

if(!is_high){
low_elevation_reconstruction(n_thread, nt_j, R, lower_thetapos, upper_thetapos, lower_phipos, upper_phipos, max_theta, max_phi, max);
}
else{
high_elevation_reconstruction(n_thread, nt_j, R, lower_thetapos, upper_thetapos, lower_phipos, upper_phipos, max_theta, max_phi, max);
}

tmpt_constructed_theta[nt_j] = max_theta[0];
tmpt_constructed_phi[nt_j] = max_phi[0];
tmpt_constructed_X_cor[nt_j] = max[0];
cout<<max_theta[0]<<"   "<<max_phi[0]<<"    "<<max[0]<<endl;

}             



	
	

	
//result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
void get_cross_correlation(int thread, int nt_j){	






for(int pair = 0; pair <8; pair++){


//606 1147
 for (Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
       if((k<606)||(k>1147)) {ant_re[thread][pair][k] = 0.; ant_im[thread][pair][k] = 0.;}
                             }
  
          
double FreqBin = 1250./4096.;
   int WeakNotchHalfWidth = round(3./FreqBin);  //MHz +-
   int StrongNotchHalfWidth = round(6./FreqBin);  //MHz +-


  const int NKnownCW = 14;
  double   CWPeak[ NKnownCW ] = {
    244.2/FreqBin, 249.2/FreqBin, 200./FreqBin, 260.8/FreqBin, 300/FreqBin, 255.8/FreqBin, 244.2/FreqBin, 260.85/FreqBin, 200/FreqBin, 232.5/FreqBin, 244.2/FreqBin, 249.2/FreqBin, 262.5/FreqBin, 300/FreqBin
  };  

  int CWPeakBin[ NKnownCW ]={0};
  for(int j=0;j<NKnownCW;j++) CWPeakBin[j] = TMath::Nint( round(CWPeak[j]) );
  
  for(int i=0; i<NKnownCW; i++){
//cout<<"(Filter) notch: "<< CWPeak[i]<<endl;
   for(int n= -WeakNotchHalfWidth; n <= WeakNotchHalfWidth; n++){
    //cout<<"Notch: "<< CWPeak[i]+n <<endl;
    //ant_re[thread][pair][ CWPeakBin[i]+n ] *= 0.;
    //ant_im[thread][pair][ CWPeakBin[i]+n ] *= 0.;
   }
  }


}










double Sum_x[8] = {0};


for(int channel = 0; channel<8; channel++){

          fft_sum[thread]->SetPointsComplex(ant_re[thread][channel],ant_im[thread][channel]);
          fft_sum[thread]->Transform();
          time_signal[thread] = fft_sum[thread]->GetPointsReal();
          
          for(int i = 0; i <4096; i++){
          Sum_x[channel] += (time_signal[thread][i]/(4096.))*(time_signal[thread][i]/(4096.));
          
          
          }
       }
	//cout<<Sum_x[0]<<"  "<<Sum_x[1]<<"  "<<Sum_x[4]<<"  "<<Sum_x[5]<<"  "<<Sum_x[7]<<"  "<<endl;







   
   	int i = 0;
    int j = 2;

   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 1; j = 3;
        break;
      case 1:
        i = 1; j = 5;
        break;
      case 2:
        i = 1; j = 7;
        break;
      case 3:
        i = 3; j = 5;
        break;
      case 4:
        i = 3; j = 7;
        break;
      case 5:
        i = 5; j = 7;
        break;

    }







        double Correlator_re[204801] = {0};            
        double Correlator_im[204801] = {0};

                  	   for (Int_t k=1; k<2049; k++){//0 means DC component
            Correlator_re[k] = (ant_re[thread][i][k] * ant_re[thread][j][k] + ant_im[thread][i][k] * ant_im[thread][j][k]); 
                   
            Correlator_im[k] = (ant_re[thread][i][k] * ant_im[thread][j][k] - ant_im[thread][i][k] * ant_re[thread][j][k]);
            //cout<<ant_re[i][k]<<endl;
                
                                                    }//k  


        //Double_t * Correlator_time = new Double_t [409600];

          fft_back[thread]->SetPointsComplex(Correlator_re,Correlator_im);
          fft_back[thread]->Transform();
          Correlator_time[thread] = fft_back[thread]->GetPointsReal();
          




         double factor = sqrt(Sum_x[i])*sqrt(Sum_x[j])*4096.;
         double max = 0;
         int mm;
                       for (Int_t k=0; k<409600; k++){
                        //cout<<Correlator_time[k]<<endl;
                correlation_pair[thread][pair][k] = Correlator_time[thread][k]/factor;
                if (correlation_pair[thread][pair][k] > max){max = correlation_pair[thread][pair][k]; mm = k;}
                //cout<<correlation_pair[pair][k]<<endl;
                                                    }//k 

//cout<<max<<"  "<<mm<<endl;

if(mm>204799) {mm -= 409600;}

//cout<<"max crosscorrelation:   pair max mm:   "<<pair<<"   "<<max<<"  "<<mm<<endl;

     }// pair loop



}
	
	
	


int filter_time(double f_x[4096], double new_x[4097], int Max_i){

//int Max_i = 0;
//double Max = 0.;

           //for (int i=0;i<1500;i++)
	//{

	  //if(fabs(x[i])>Max) {Max = fabs(x[i]);  Max_i = i;} 	  

	//}
	//cout<<Max_i<<endl;
	
//hamming window
double a0 = 0.53836;
double a1 = 1. - a0;

 	for(int i=0;i<1500;i++)
		{
		new_x[i] = ((i<Max_i-75)||(i>Max_i+50))?0.:f_x[i];//(a0 - a1*cos(2.*PI*(i-Max_i+200)/(400.-1.)))*f_x[i];
		}

/*
double max=0.;
int max_i = 0;
 	for(int i=300;i<1200;i++)
	{
		if((i>Max_i+50)){
			if(fabs(f_x[i])>max){
				max = fabs(f_x[i]);
				max_i = i;
			}
		}
	}

 	for(int i=0;i<1500;i++)
	{
		new_x[i] = ((i<max_i-50)||(i>max_i+50))?0.:f_x[i];//(a0 - a1*cos(2.*PI*(i-Max_i+200)/(400.-1.)))*f_x[i];
	}
*/   	
double Re_ant[2049];
double Im_ant[2049];

           for (int i=1500;i<4096;i++)
	{
   	       new_x[i] = 0;	                               
	}
	
return Max_i;

}





int filter_fre(int thread, double x[4097], double f_x[4096]){


double window_x[4096] = {0};


//blackman window
double a0 = 7938./18608.;
double a1 = 9240./18608.;
double a2 = 1430./18608.;
 	for(int i=0;i<1500;i++)
		{
		window_x[i] = (a0 - a1*cos(2.*PI*i/1500.) + a2*cos(4.*PI*i/1500.))*x[i];
		}
   	
double Re_ant[2049];
double Im_ant[2049];

/*
//hamming window
double a0 = 0.53836;
double a1 = 1. - a0;
 	for(int i=0;i<1500;i++)
		{
		window_x[i] = (a0 - a1*cos(2.*PI*i/(1500.-1.)))*f_x[i];
		}
*/

   filter_fft_forward[thread]->SetPoints(x);
   filter_fft_forward[thread]->Transform();
   filter_fft_forward[thread]->GetPointsComplex(Re_ant,Im_ant);

 for (Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
       if((k<606)||(k>1147)) {Re_ant[k] = 0.; Im_ant[k] = 0.;}
                             }
  
          
          
double FreqBin = 1250./4096.;
   int WeakNotchHalfWidth = round(3./FreqBin);  //MHz +-
   int StrongNotchHalfWidth = round(6./FreqBin);  //MHz +-


  const int NKnownCW = 14;
  double   CWPeak[ NKnownCW ] = {
    244.2/FreqBin, 249.2/FreqBin, 200./FreqBin, 260.8/FreqBin, 300/FreqBin, 255.8/FreqBin, 244.2/FreqBin, 260.85/FreqBin, 200/FreqBin, 232.5/FreqBin, 244.2/FreqBin, 249.2/FreqBin, 262.5/FreqBin, 300/FreqBin
  };


  int CWPeakBin[ NKnownCW ]={0};
  for(int j=0;j<NKnownCW;j++) CWPeakBin[j] = TMath::Nint( round(CWPeak[j]) );
  
  for(int i=0; i<NKnownCW; i++){
//cout<<"(Filter) notch: "<< CWPeak[i]<<endl;
   for(int n= -WeakNotchHalfWidth; n <= WeakNotchHalfWidth; n++){
    //cout<<"Notch: "<< CWPeak[i]+n <<endl;
    Re_ant[ CWPeakBin[i]+n ] *= 0.;
    Im_ant[ CWPeakBin[i]+n ] *= 0.;
   }
  }

//Double_t * filtered_x = new Double_t [1500];

  filter_fft_back[thread]->SetPointsComplex(Re_ant,Im_ant);
  filter_fft_back[thread]->Transform();
  filtered_x[thread] = filter_fft_back[thread]->GetPointsReal();



int Max_i = 0;
double Max = 0.;

//inversed hamming windows
/*
 	for(int i=0;i<1500;i++)
		{
		f_x[i] = filtered_x[i]/4096./(a0 - a1*cos(2.*PI*i/(1500.-1.)));
		}
*/
for(int i=0;i<4096;i++){
f_x[i] = filtered_x[thread][i]/4096.;
  						}

for(int i=200;i<1000;i++){
if(fabs(filtered_x[thread][i])>Max) {Max = fabs(filtered_x[thread][i]);  Max_i = i;}
  						}
//delete fft_back;
//fft_back = 0;
//delete fft_forward;					
//fft_forward = 0;
cout<<Max_i<<endl;						
return Max_i;

}

//theta:-pi/2.~pi/2;  phi:0~2*pi, x-axis direction is phi = 0;
//xyz: antenna2 - antenna1
//return: relative time between signal arrive antenna2 and antenna1, if it is minus, means signal arrive antenna2 first.
Double_t expectedTimeDiff_angle(double theta, double phi, double *xyz)
         {
   double n[3] = {0};

   phi = phi*TMath::DegToRad();
   theta = theta*TMath::DegToRad();
   
   n[0] = (-1)*TMath::Cos(theta)*TMath::Cos(phi);
   n[1] = (-1)*TMath::Cos(theta)*TMath::Sin(phi);
   n[2] = (-1)*TMath::Sin(theta);

   double dt = (n[0]*xyz[0] + n[1]*xyz[1] + n[2]*xyz[2])/(Speed_Of_Light);


   // s to ns
   return dt/1e-9;
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
         
         
         


     	     						


       	     						
void deconvolute_response(int thread, int is_high, int pos_theta, int pos_phi){

  
	for(int ant = 0; ant < 4; ant ++){
	
    
 		int max_i = filter_fre(thread, x[thread][ant*2+1], f_x[thread][ant*2+1]);
	    filter_time(f_x[thread][ant*2+1], new_x[thread][ant*2+1], max_i);  

   double re[2049] = {0};
   double im[2049] = {0};


   decon_fft_forward[thread]->SetPoints(new_x[thread][ant*2+1]);
   decon_fft_forward[thread]->Transform();
   decon_fft_forward[thread]->GetPointsComplex(re,im);
   
   for(int i = 606; i<1148; i++){
        //if((i<606)||(i>1147)) {ant_re[ant*2][i] = 0; ant_im[ant*2][i] = 0; continue;}
   		double magnitude = sqrt(re[i] * re[i] + im[i] * im[i]);
   		
   		//h1->SetBinContent(i-606, log10(magnitude));

     	double S = magnitude - noise_mag[ant*2+1][i];
     	double SNR = S/noise_mag[ant*2+1][i];
     	double factor = 1./(1. + 1./pow(SNR,2));
     	
     	   		
   		double phase = TMath::ATan(im[i]/(re[i]+1.0e-12));
     	if(re[i]<0){phase += TMath::Pi();}
     	if(phase>TMath::Pi()) {phase -= 2*TMath::Pi();}
     	double new_mag = factor * magnitude * gd_responds[is_high].resp_mag[pos_theta][pos_phi][ant][i];

     	double new_pha = phase + gd_responds[is_high].resp_pha[pos_theta][pos_phi][ant][i];
     	new_pha = fmod(new_pha, 2*PI);
   		if(new_pha>PI) new_pha = new_pha - 2*PI;
   		if(new_pha<(-1*PI)) new_pha = new_pha + 2*PI;
   		
     	cyc_mag[thread][ant][i] = new_mag;
     	cyc_pha[thread][ant][i] = new_pha;
     	
     	// TO DO //Do not using deconvolute
     	//re[i] = new_mag*TMath::Cos(new_pha);
     	//im[i] = new_mag*TMath::Sin(new_pha);
     	ant_re[thread][ant*2+1][i] = re[i];
     	ant_im[thread][ant*2+1][i] = im[i];
     	     				 }  
     	     				 			}
     
}




double delta_t_phase[2048] = {0};

double get_delay_from_phase(int thread, int ant_i, int ant_j, double delay, TH1F * h1){

   	double mean = 0;
   	double m_count = 0;
   	for(int i = 606; i<1148; i++){
   		double extress = 2*PI*(delay*1.0e-9)*(i*(1250e6)/4096.);
   		extress = fmod(extress, 2*PI); 
   		//if(extress>PI) {extress = extress - 2 * PI;}
   		//if(extress<(-1*PI)) {extress = extress + 2*PI;}
   		double delta_phase = cyc_pha[thread][ant_j][i] - cyc_pha[thread][ant_i][i];
   		//if(delta_phase>PI) delta_phase = delta_phase - 2*PI;
   		//if(delta_phase<(-1*PI)) delta_phase = delta_phase + 2*PI;
   		delta_phase  = delta_phase + extress;
   		delta_phase = fmod(delta_phase, 2*PI);
   		if(delta_phase>PI) delta_phase = delta_phase - 2*PI;
   		if(delta_phase<(-1*PI)) delta_phase = delta_phase + 2*PI;
   		delta_t_phase[i] = (delta_phase)/(2*PI*(i*(1250e6)/4096.));
   		delta_t_phase[i] *= 1.0e9;//s to ns
   		if(fabs(delta_t_phase[i])<2){
   		double factor = 10*log10(sqrt(cyc_mag[thread][ant_i][i]*cyc_mag[thread][ant_j][i]));
   		//cout<<factor<<"  "<<cyc_pha[ant_i][i]<<"   "<<cyc_pha[ant_j][i]<<endl;
   		mean += factor*delta_t_phase[i];
   		m_count += factor;
   		}
   		//h1->Fill(delta_t_phase[i]);
   	}
   	if(m_count>0) {mean /= m_count;}
   	else {mean = 0;}
   	
   	//cout<<"mean is: "<<mean<<endl;
   	double sigma = 0;
   	double fine_mean = 0;
   	double count = 0;
   	   	for(int i = 606; i<1148; i++){
   	   	
   	   		if((delta_t_phase[i]>mean-0.5)&&(delta_t_phase[i]<mean+0.5)){
   	   			double factor = 10*log10(sqrt(cyc_mag[thread][ant_i][i]*cyc_mag[thread][ant_j][i]));
   	   			//double factor = sqrt(cyc_mag[ant_i][i]*cyc_mag[ant_j][i]);
   	   			count += factor;
   	   			//sigma += factor*pow((delta_t_phase[i] - mean), 2);
   	   			fine_mean += factor*delta_t_phase[i];
   	   		}
   		
   		}
		
		if(count>0) {fine_mean /= count;}
		else {fine_mean = 0;}

   	   	for(int i = 606; i<1148; i++){
   	   	
   	   		if((delta_t_phase[i]<fine_mean-0.5)||(delta_t_phase[i]>fine_mean+0.5)){
   	   			for(int ant = 0; ant<4; ant++){
   	   			    ant_re[thread][ant*2][i] = 0.;
   	   			    ant_im[thread][ant*2][i] = 0.;
   	   			}
   	   		}
   		
   		}		
		

	return fabs(fine_mean);
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



template <typename T>
T*** new_Array3D(int height, int row, int col)  
{  
   int i,j,k;   // p[height][row][col]

   T ***p;
p = new T **[height];
for(i=0; i<height; i++)
{
   p[i]=new T *[row];
   for(j=0; j<row; j++)
   p[i][j]=new T[col];
}

//finish creating use p[i][j][k] to access the data  
for(i=0; i<height; i++)  
{
   for(j=0; j<row; j++)  
         {
    for(k=0;k<col;k++)
    {
     p[i][j][k]=0;
     //cout<<p[i][j][k]<<" ";
    }
        //cout<<endl;
   }
    //cout<<endl;
}
return p;
}


template <typename T>  
void delete_Array3D(T ***p, int height, int row, int col)
{
     //free the memory
     int i,j,k; 
for(i=0; i<height; i++)
{
       for(j=0; j<col; j++)
       {  
           delete [] p[i][j];  
       }  
   }      
   for(i=0; i<height; i++)  
{      
      delete [] p[i];  
    }  
    delete [] p;  

}

void get_noise_level(double noise_mag[8][2049]){

	float h1, h2, h3, h4, v1, v2, v3, v4;
   
   FILE *fpin = fopen("/home/cyc/software/TAROGE-4_analysis/differential_response_deconvolution/noise_level.txt","r");// 600m
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


void read_responses_low(){

        Int_t azimuth = 0.;
        Int_t elevation = 0.;
        Double_t responses_magnitude_1V[2049] = {0.};
        Double_t responses_phase_1V[2049] = {0.};
        Double_t responses_magnitude_2V[2049] = {0.};
        Double_t responses_phase_2V[2049] = {0.};
        Double_t responses_magnitude_3V[2049] = {0.};
        Double_t responses_phase_3V[2049] = {0.};
        Double_t responses_magnitude_4V[2049] = {0.};
        Double_t responses_phase_4V[2049] = {0.};
    
    const char * res_fname =  "/home/cyc/software/TAROGE-4_pulser/pulser/V_pol/ground_responses.root";
    TFile * res_file = new TFile(res_fname);
    TTree * responsesTree = (TTree*) res_file->Get("ground_responses");
    
        responsesTree->SetBranchAddress("azimuth", &azimuth);
        responsesTree->SetBranchAddress("elevation", &elevation);
        
        responsesTree->SetBranchAddress("responses_magnitude_1V", &responses_magnitude_1V);
        responsesTree->SetBranchAddress("responses_phase_1V", &responses_phase_1V);
        
        responsesTree->SetBranchAddress("responses_magnitude_2V", &responses_magnitude_2V);
        responsesTree->SetBranchAddress("responses_phase_2V", &responses_phase_2V);
        
        responsesTree->SetBranchAddress("responses_magnitude_3V", &responses_magnitude_3V);
        responsesTree->SetBranchAddress("responses_phase_3V", &responses_phase_3V);
        
        responsesTree->SetBranchAddress("responses_magnitude_4V", &responses_magnitude_4V);
        responsesTree->SetBranchAddress("responses_phase_4V", &responses_phase_4V);
        
   for(int entry = 0; entry < responsesTree->GetEntries(); entry++){
   		responsesTree->GetEntry(entry);
   		for(int i = 0; i<2049; i++){
   			if(elevation<0 || elevation>140 || azimuth<0 || azimuth>180){
   			//cout<<elevation<<"   "<<azimuth<<endl;
   			}
   			
   			gd_responds[0].resp_mag[elevation][azimuth][0][i] = responses_magnitude_1V[i];
   			gd_responds[0].resp_pha[elevation][azimuth][0][i] = responses_phase_1V[i];
   			
   			gd_responds[0].resp_mag[elevation][azimuth][1][i] = responses_magnitude_2V[i];
   			gd_responds[0].resp_pha[elevation][azimuth][1][i] = responses_phase_2V[i];
   			
   			gd_responds[0].resp_mag[elevation][azimuth][2][i] = responses_magnitude_3V[i];
   			gd_responds[0].resp_pha[elevation][azimuth][2][i] = responses_phase_3V[i];
   			
   			gd_responds[0].resp_mag[elevation][azimuth][3][i] = responses_magnitude_4V[i];
   			gd_responds[0].resp_pha[elevation][azimuth][3][i] = responses_phase_4V[i];
   			
   			//if(responses_magnitude_1V[i]!=0) cout<<resp_mag[elevation][azimuth][0][i]<<"   "<<responses_magnitude_1V[i]<<endl;
   			
   			}
   		}
	
}


void read_responses_high(){

        Int_t azimuth = 0.;
        Int_t elevation = 0.;
        Double_t responses_magnitude_1V[2049] = {0.};
        Double_t responses_phase_1V[2049] = {0.};
        Double_t responses_magnitude_2V[2049] = {0.};
        Double_t responses_phase_2V[2049] = {0.};
        Double_t responses_magnitude_3V[2049] = {0.};
        Double_t responses_phase_3V[2049] = {0.};
        Double_t responses_magnitude_4V[2049] = {0.};
        Double_t responses_phase_4V[2049] = {0.};
    
    const char * res_fname =  "/home/cyc/software/TAROGE-4_pulser/pulser/V_pol/ground_responses_high.root";
    TFile * res_file = new TFile(res_fname);
    TTree * responsesTree = (TTree*) res_file->Get("ground_responses");
    
        responsesTree->SetBranchAddress("azimuth", &azimuth);
        responsesTree->SetBranchAddress("elevation", &elevation);
        
        responsesTree->SetBranchAddress("responses_magnitude_1V", &responses_magnitude_1V);
        responsesTree->SetBranchAddress("responses_phase_1V", &responses_phase_1V);
        
        responsesTree->SetBranchAddress("responses_magnitude_2V", &responses_magnitude_2V);
        responsesTree->SetBranchAddress("responses_phase_2V", &responses_phase_2V);
        
        responsesTree->SetBranchAddress("responses_magnitude_3V", &responses_magnitude_3V);
        responsesTree->SetBranchAddress("responses_phase_3V", &responses_phase_3V);
        
        responsesTree->SetBranchAddress("responses_magnitude_4V", &responses_magnitude_4V);
        responsesTree->SetBranchAddress("responses_phase_4V", &responses_phase_4V);
        
   for(int entry = 0; entry < responsesTree->GetEntries(); entry++){
   		responsesTree->GetEntry(entry);
   		for(int i = 0; i<2049; i++){
   			if(elevation<0 || elevation>=36 || azimuth<0 || azimuth>=18){
   			//cout<<elevation<<"   "<<azimuth<<endl;
   			}
   			
   			gd_responds[1].resp_mag[elevation][azimuth][0][i] = responses_magnitude_1V[i];
   			gd_responds[1].resp_pha[elevation][azimuth][0][i] = responses_phase_1V[i];
   			
   			gd_responds[1].resp_mag[elevation][azimuth][1][i] = responses_magnitude_2V[i];
   			gd_responds[1].resp_pha[elevation][azimuth][1][i] = responses_phase_2V[i];
   			
   			gd_responds[1].resp_mag[elevation][azimuth][2][i] = responses_magnitude_3V[i];
   			gd_responds[1].resp_pha[elevation][azimuth][2][i] = responses_phase_3V[i];
   			
   			gd_responds[1].resp_mag[elevation][azimuth][3][i] = responses_magnitude_4V[i];
   			gd_responds[1].resp_pha[elevation][azimuth][3][i] = responses_phase_4V[i];
   			
   			//if(responses_magnitude_1V[i]!=0) cout<<resp_mag[elevation][azimuth][0][i]<<"   "<<responses_magnitude_1V[i]<<endl;
   			
   			}
   		}
	
}


bool SetDataTreeBranches(TTree* tr){

  printf("SetDataTree SetBranchAddresses(): %s\t%s\n", tr->GetName(), tr->GetTitle() );
  tr-> SetBranchAddress("event",&gEventNum);
  tr-> SetBranchAddress("bGood",&bGood);     
  tr-> SetBranchAddress("timeStamp_FPGA", &timeStamp_FPGA); //TAROGE-4, FPGA 24-bit counter    
  tr-> SetBranchAddress("triggerBits",&gTrigBits ); //TAROGE-4  TBits, align with ChIDs    Hpol and Vpol
  tr-> SetBranchAddress("overVoltBits",&gOverVoltBits ); //TAROGE-4  TBits(kNCh)   overvoltage or not
  tr-> SetBranchAddress("eventTime",&eventTime); //PC timestamp

  tr-> SetBranchAddress("T1H", &T1H);
  tr-> SetBranchAddress("T2H", &T2H);
  tr-> SetBranchAddress("T3H", &T3H);
  tr-> SetBranchAddress("T4H", &T4H);
  tr-> SetBranchAddress("T1V", &T1V);
  tr-> SetBranchAddress("T2V", &T2V);
  tr-> SetBranchAddress("T3V", &T3V);
  tr-> SetBranchAddress("T4V", &T4V);
 


 
  return true;
}


void low_elevation_reconstruction(int n_thread, int nt_j, double R, int lower_thetapos, int upper_thetapos, int lower_phipos, int upper_phipos, double max_theta[3], double max_phi[3], double max[3]){

for(int thetapos = lower_thetapos; thetapos<upper_thetapos; thetapos++){
	for(int phipos = lower_phipos; phipos<upper_phipos; phipos++){



deconvolute_response(n_thread, 0, thetapos, phipos);


    double theta = (thetapos-70)/10.;
    double phi = (phipos - 90)*0.6 - 25.;


   for(int pair = 0; pair <6; pair++){

     int i = 0;
     int j = 1;

    switch(pair){
      case 0:
        i = 4; j = 5;
        break;
      case 1:
        i = 4; j = 6;
        break;
      case 2:
        i = 4; j = 7;
        break;
      case 3:
        i = 5; j = 6;
        break;
      case 4:
        i = 5; j = 7;
        break;
      case 5:
        i = 6; j = 7;
        break;

    }

    double xyz[3] = {0};
    for(int k=0; k<3; k++){xyz[k] = coor_ant[j][k] - coor_ant[i][k];}

    double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;//ns

    //get_delay_from_phase(thread, i, j, delta_t, h1[pair]);
	
    }//pair loop



get_cross_correlation(n_thread, nt_j);



double inter_map[5][13] = {0};
double sum_X_coe = 0;
   for(int pair = 0; pair <6; pair++){

     int i = 0;
     int j = 1;

    switch(pair){
      case 0:
        i = 4; j = 5;
        break;
      case 1:
        i = 4; j = 6;
        break;
      case 2:
        i = 4; j = 7;
        break;
      case 3:
        i = 5; j = 6;
        break;
      case 4:
        i = 5; j = 7;
        break;
      case 5:
        i = 6; j = 7;
        break;

    }

    double xyz[3] = {0};
    for(int k=0; k<3; k++){xyz[k] = coor_ant[j][k] - coor_ant[i][k];}
    

          
          
          

    for(int i_theta=0; i_theta<5; i_theta++){
      for(int i_phi=0; i_phi<13; i_phi++){
          double theta = -0.3 + (0.15)*i_theta + (thetapos-70)/10.;
          double phi = -0.6 + (0.1)*i_phi + (phipos - 90)*0.6 - 25.;
          double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;

          delta_t = delta_t/(0.8*4096./409600.);
          if(delta_t<0) delta_t += 409600;
          if(delta_t<0||delta_t>409600){continue;cout<<"error"<<endl;}

          inter_map[i_theta][i_phi] += correlation_pair[n_thread][pair][(int) delta_t];
        }
      }
      
    if(pair == 5) {
     
     for (int ci=0; ci<5; ci++){
      for (int cj=0; cj<13; cj++){
      
      		
          if(inter_map[ci][cj] > max[0])
          {
           max[2] = max[1]; max[1] = max[0]; max[0] = inter_map[ci][cj];
           max_theta[2] = max_theta[1]; max_theta[1] = max_theta[0]; max_theta[0] = -0.3 + (0.15)*ci + (thetapos-70)/10.;
           max_phi[2] = max_phi[1]; max_phi[1] = max_phi[0]; max_phi[0] = -0.6 + (0.1)*cj + (phipos - 90)*0.6 -25;
          }
     }
   }

   
   
     }


    }//pair loop

//cout<<thetapos<<"   "<<phipos<<"   "<<sum_X_coe<<endl;
//cout<<endl;

	}//posphi loop
}//postheta loop

}
















void high_elevation_reconstruction(int n_thread, int nt_j, double R, int lower_thetapos, int upper_thetapos, int lower_phipos, int upper_phipos, double max_theta[3], double max_phi[3], double max[3]){

for(int thetapos = lower_thetapos; thetapos<upper_thetapos; thetapos++){
	for(int phipos = lower_phipos; phipos<upper_phipos; phipos++){



deconvolute_response(n_thread, 1, thetapos, phipos);


    double theta = thetapos + 7.;
    double phi = (phipos - 9.)*6. + 3 - 25.;


   for(int pair = 0; pair <6; pair++){

     int i = 0;
     int j = 1;

    switch(pair){
      case 0:
        i = 4; j = 5;
        break;
      case 1:
        i = 4; j = 6;
        break;
      case 2:
        i = 4; j = 7;
        break;
      case 3:
        i = 5; j = 6;
        break;
      case 4:
        i = 5; j = 7;
        break;
      case 5:
        i = 6; j = 7;
        break;

    }

    double xyz[3] = {0};
    for(int k=0; k<3; k++){xyz[k] = coor_ant[j][k] - coor_ant[i][k];}

    double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;//ns

    //get_delay_from_phase(thread, i, j, delta_t, h1[pair]);
    	
    }//pair loop



get_cross_correlation(n_thread, nt_j);



double inter_map[5][21] = {0};
double sum_X_coe = 0;
   for(int pair = 0; pair <6; pair++){

     int i = 0;
     int j = 1;

    switch(pair){
      case 0:
        i = 4; j = 5;
        break;
      case 1:
        i = 4; j = 6;
        break;
      case 2:
        i = 4; j = 7;
        break;
      case 3:
        i = 5; j = 6;
        break;
      case 4:
        i = 5; j = 7;
        break;
      case 5:
        i = 6; j = 7;
        break;

    }

    double xyz[3] = {0};
    for(int k=0; k<3; k++){xyz[k] = coor_ant[j][k] - coor_ant[i][k];}
    

          
          
          

    for(int i_theta=0; i_theta<5; i_theta++){
      for(int i_phi=0; i_phi<21; i_phi++){
          double theta = -0.6 + 0.3*i_theta + thetapos + 7.;
          double phi = -3. + (0.3)*i_phi + (phipos - 9.)*6. + 3 - 25.;
          double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;

          delta_t = delta_t/(0.8*4096./409600.);
          if(delta_t<0) delta_t += 409600;
          if(delta_t<0||delta_t>409600){continue;cout<<"error"<<endl;}

          inter_map[i_theta][i_phi] += correlation_pair[n_thread][pair][(int) delta_t];
        }
      }
      
    if(pair == 5) {
     
     for (int ci=0; ci<5; ci++){
      for (int cj=0; cj<21; cj++){
      		
          if(inter_map[ci][cj] > max[0])
          {
           max[2] = max[1]; max[1] = max[0]; max[0] = inter_map[ci][cj];
           max_theta[2] = max_theta[1]; max_theta[1] = max_theta[0]; max_theta[0] = -0.6 + 0.3*ci + thetapos + 7.;
           max_phi[2] = max_phi[1]; max_phi[1] = max_phi[0]; max_phi[0] = -3. + 0.3*cj + (phipos - 9.)*6. + 3 - 25.;
          }
     }
   }

   
   
     }


    }//pair loop

//cout<<thetapos<<"   "<<phipos<<"   "<<sum_X_coe<<endl;
//cout<<endl;

	}//posphi loop
}//postheta loop

}
