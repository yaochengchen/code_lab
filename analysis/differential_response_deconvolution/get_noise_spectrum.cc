//20240417 change RECONSTRUCT_PULSER_DATA_PARALLEL.CC x[4096] new_x[4096] to 4097; since filter_n_size & decon_n_size is 4097, if not change, some times wierd number after FFT will appears. It is betters to change get_noise_spectrum.cc, change n_size to be n instead of n+1.;


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
using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8


int filter_fre(double x[4096], double f_x[4096]);
int filter_time(double x[4096], double new_x[4096], int Max_i);
void get_spectrum();


 Double_t x[8][4096]={0};//8 antennas
 Double_t new_x[8][4096]={0};//8 antennas
 Double_t f_x[8][4096]={0};//8 antennas
 
 double Re_template[2049] = {0};
double Im_template[2049] = {0};

	double tmp_mag[8][2049] = {0};
	double tmp_pha[8][2049] = {0};


void get_noise_spectrum(){


  Double_t t[1500];
  Short_t T1H[1500];
  Short_t T2H[1500];
  Short_t T3H[1500];
  Short_t T4H[1500];
  Short_t T1V[1500];
  Short_t T2V[1500];
  Short_t T3V[1500];
  Short_t T4V[1500];
  bool bGood;
  Int_t timeStamp_FPGA; 

  Int_t event_number;
  const char * fname_evt;
  TFile *file;
  TTree *Tree_Muon;
  double eventTime = 0;
    TBits* triggerBitsccc = 0;
    
    int count = 0;
    
    for(int number = 27131; number< 27177; number++){
  
    fname_evt =  Form("/media/cyc/For_Linux/TAROGE4_DATA/data/20211212/run%08d.root",number);//1212
    file = new TFile(fname_evt);
    Tree_Muon = (TTree*) file->Get("t");
          Tree_Muon->SetBranchAddress("bGood",&bGood);
		  Tree_Muon->SetBranchAddress("timeStamp_FPGA",&timeStamp_FPGA);
		  Tree_Muon->SetBranchAddress("triggerBits",&triggerBitsccc);
		  Tree_Muon->SetBranchAddress("eventTime", &eventTime);
		  Tree_Muon->SetBranchAddress("T1H",T1H);
		  Tree_Muon->SetBranchAddress("T2H",T2H);
		  Tree_Muon->SetBranchAddress("T3H",T3H);
		  Tree_Muon->SetBranchAddress("T4H",T4H);
		  Tree_Muon->SetBranchAddress("T1V",T1V);
		  Tree_Muon->SetBranchAddress("T2V",T2V);
		  Tree_Muon->SetBranchAddress("T3V",T3V);
		  Tree_Muon->SetBranchAddress("T4V",T4V);



	
	

	for(int j = 0; j < Tree_Muon->GetEntries(); j++){

      Tree_Muon->GetEntry(j);
      
      if(!bGood){continue;}
      if(triggerBitsccc->CountBits() != 0) {continue;}
      



  for (int i=0;i<1500;i++)
	{
	  t[i]=i*0.8;
	  x[0][i]= (float)T1H[i]*500/32512;
	  x[2][i]= (float)T2H[i]*500/32512;
	  x[4][i]= (float)T3H[i]*500/32512;
	  x[6][i]= (float)T4H[i]*500/32512;
	  x[1][i]= (float)T1V[i]*500/32512;
	  x[3][i]= (float)T2V[i]*500/32512;
	  x[5][i]= (float)T3V[i]*500/32512;
	  x[7][i]= (float)T4V[i]*500/32512;	  	  
	}
	
	get_spectrum();
	count ++;


}//j loop
}//run loop

cout<<count<<endl;
   for(int ant = 0; ant < 8; ant ++){
   		for(int i = 1; i<2049; i++){
     		tmp_mag[ant][i] /= (double) count;
     	     						}
     	     				}

  
  const char * out_filename = "./noise_level.txt";
  fstream out_put;
  out_put.open(out_filename); 
   
   for(int i = 0; i<2049; i++){
   		out_put<<i<<"  ";
   		for(int ant = 0; ant < 8; ant ++){	
     		out_put<<tmp_mag[ant][i]<<"  ";
     	     						}
     	out_put<<endl;
     	     				}
//cout<<out_put<<endl;
out_put.close();
}// main function end


   int N = 4096;
   Int_t n_size = N+1;
   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &n_size,"R2C ES K");


void get_spectrum(){
	for(int ant = 0; ant < 8; ant ++){
	
   
   	filter_fre(x[ant], f_x[ant]);
   	int max_i = 500;
	filter_time(f_x[ant], new_x[ant], max_i);
   //filter_time(x[ant], new_x[ant], max_i);  


   fft_forward->SetPoints(new_x[ant]);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(Re_template,Im_template);

   for(int i = 1; i<2049; i++){
   		double magnitude = sqrt(Re_template[i] * Re_template[i] + Im_template[i] *  Im_template[i]);
   		double phase = TMath::ATan(Im_template[i]/(Re_template[i]+1.0e-12));
     	if(Re_template[i]<0){phase += TMath::Pi();}
     	if(phase>TMath::Pi()) {phase -= 2*TMath::Pi();}
     	tmp_mag[ant][i] += magnitude;
     	tmp_pha[ant][i] = phase;// - 2*PI*(10e-9)*(i*(1250e6)/4096.);
     	     						}

	}
}




int filter_time(double f_x[4096], double new_x[4096], int Max_i){

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
		new_x[i] = ((i<Max_i-100)||(i>Max_i+100))?0:(a0 - a1*cos(2.*PI*(i-Max_i+100)/(200.-1.)))*f_x[i];
		}
   	
double Re_ant[2049];
double Im_ant[2049];

           for (int i=1500;i<4096;i++)
	{
   	       new_x[i] = 0;	                               
	}
	
return Max_i;

}



int filter_N = 4096;
Int_t filter_n_size = filter_N+1;
TVirtualFFT *filter_fft_forward = TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
TVirtualFFT *filter_fft_back = TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");
Double_t * filtered_x = new Double_t [4096];

int filter_fre(double x[4096], double f_x[4096]){

//int N = 1500;
//Int_t n_size = N+1;
   //TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &n_size,"R2C ES K");
      
   //TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");

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



   filter_fft_forward->SetPoints(x);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(Re_ant,Im_ant);

 for (Int_t k=0; k<2049; k++){//0 means DC component // from 185-350 MHz
       if((k<606)||(k>1147)) {Re_ant[k] = 0.; Im_ant[k] = 0.;}
                             }
  
          
/*          
double FreqBin = 1250./4096.;
   int WeakNotchHalfWidth = round(3./FreqBin);  //MHz +-
   int StrongNotchHalfWidth = round(6./FreqBin);  //MHz +-


  const int NKnownCW = 8;
  double   CWPeak[ NKnownCW ] = {
    165./FreqBin, 187.5/FreqBin, 244.167/FreqBin, 249.167/FreqBin, 272.5/FreqBin, 
    299.165/FreqBin, 200./FreqBin, 262./FreqBin  //in cal-pulser, not found in forced trig
      //187.5, 272.5 are the strongest
  };  //5.e6, 105.e6, 398.333e6, 483.333e6
  //{ 180.83e6/FreqBin, 217.5e6/FreqBin, 272.5e6/FreqBin, 276.65e6/FreqBin, 313.6e6/FreqBin, 500.e6/FreqBin};
  // { 104.9e6/FreqBin, 107.7e6/FreqBin, 115.05e6/FreqBin, 125.83e6/FreqBin, 144.5e6/FreqBin, 146.68e6/FreqBin, 155.84e6/FreqBin, 157.55e6/FreqBin, 161.6e6/FreqBin, 165.e6/FreqBin, 120.e6/FreqBin, 124.e6/FreqBin, 128.e6/FreqBin, 141.e6/FreqBin, 245.e6/FreqBin, 254.e6/FreqBin, 283.e6/FreqBin}; ################# 500 is osc
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
*/
//Double_t * filtered_x = new Double_t [1500];

  filter_fft_back->SetPointsComplex(Re_ant,Im_ant);
  filter_fft_back->Transform();
  filtered_x = filter_fft_back->GetPointsReal();



int Max_i = 0;
double Max = 0.;

for(int i=0;i<4096;i++){
f_x[i] = filtered_x[i]/4096.;
  						}

for(int i=200;i<1300;i++){
if(fabs(filtered_x[i])>Max) {Max = fabs(filtered_x[i]);  Max_i = i;}
  						}
//delete fft_back;
//fft_back = 0;
//delete fft_forward;					
//fft_forward = 0;
  						
return Max_i;

}







