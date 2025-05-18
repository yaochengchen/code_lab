#include <string>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST

#define Speed_Of_Light 2.99792458e8
#define WINDOW_SPACING 150
#define HALF_WINDOW_SPACING 75
#define HALF_WINDOW_WIDTH 76
#define PI TMath::Pi()

template <typename T>
T** new_Array2D(int row, int col);

template <typename T> 
void delete_Array2D(T **arr, int row, int col);



//base center   2021 photogrametry  Tower-4 is back side
double coor_ant[8][3] = {{328914.252, 2697625.953, 709.389},{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949},{328936.041, 2697639.685, 704.949}};
double delay_ant[8] = {0,-92.4,53.2,-30.3,84.6,59.5,-1523.1,-743.1};
int peak_position[8] = {0};


int get_FFT(int ch, double x[1500], double FW_x[1500], double filtered_x[1500], double FW_re[751], double FW_im[751], double deep_filtered_x[1500]);
int filter_time(double filtered_x[1500], double SW_x[1500], int Max_i);
//HALF_WINDOW_WIDTH
void SW_fft(double SW_x[1500], int Max_i, double SW_re[751], double SW_im[751]);

void Planck_taper_windowing(double epsilon, int N, double x[], double w_x[]);



double Correlator_re[2][6][HALF_WINDOW_WIDTH] = {0};            
double Correlator_im[2][6][HALF_WINDOW_WIDTH] = {0};
double correlation_pair[2][6][32768] = {0};
 
	double FW_x[8][1500] = {0.};
	double filtered_x[8][1500] = {0.};
	double deep_filtered_x[8][1500] = {0.};
	double FW_re[8][751] = {0.};
	double FW_im[8][751] = {0.};
	double SW_x[8][1500] = {0.};
	double SUM_x[8] = {0};
	
	//HALF_WINDOW_WIDTH
	double SW_re[8][751] = {0.};
	double SW_im[8][751] = {0.};
	
void get_cross_correlation(int pol, double tdoa[2][6], double Xcor[2][6]);
double real_cross_correlation(int pol, int pair, int i, int j);
int reconstruction(int pol, double event_theta[2], double event_phi[2], double x_value[2]);
Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j);                  
Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j); 




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



TH2F *h2 = new TH2F("h2","h2", 3600, -180, 180, 1800, -90, 90);

void draw_interferometric_map(int date, int run_num, int event_num){

 //get_noise_level(noise_mag);

 //read_responses_high();

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
	
	Double_t t[1500]={0};
 	Double_t x[8][1500]={0};//8 antennas
	Double_t new_x[8][1500]={0};//8 antennas
 	//Double_t frequency[751]={0};//8 antennas
	//Double_t powspm[8][751]={0};//8 antennas
	


    fname_evt =  Form("/media/cyc/For_Linux/TAROGE4_DATA/data/%d/run%08d.root", date, run_num);
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

  for (int i=0;i<1500;i++)
	{
	  t[i]=i*0.8;
	  x[0][i]= (float)T1H[i]*500/32512;
	  x[2][i]= (float)T2H[i]*500/32512;
	  x[4][i]= (float)T3H[i]*500/32512;
	  x[6][i]= (float)T4H[i]*500/32512;
	  //int k = i+3;
	  //if(k>=1500){k-=1500;}
	  x[1][i]= (float)T1V[i]*500/32512;
	  x[3][i]= (float)T2V[i]*500/32512;
	  x[5][i]= (float)T3V[i]*500/32512;
	  x[7][i]= (float)T4V[i]*500/32512;	  	  
	}
	
	


	
	for(int ch = 0; ch<8; ch++){
		int Max_i = get_FFT(ch, x[ch], FW_x[ch], filtered_x[ch], FW_re[ch], FW_im[ch], deep_filtered_x[ch]);
		peak_position[ch] = Max_i;
		// full or short waveform
		//filter_time(deep_filtered_x[ch], SW_x[ch], Max_i);
		
		for(int i=0; i<1500; i++){
			//cout<<i<<"  "<<SW_x[ch][i]<<endl;
			//SUM_x[ch] += pow(SW_x[ch][i], 2);
			// full or short waveform
			SUM_x[ch] += pow(deep_filtered_x[ch][i], 2);
			
		}
		
		// full or short waveform
		//SW_fft(SW_x[ch], Max_i, SW_re[ch], SW_im[ch]);
	
	}
	

double tdoa[2][6] = {0.};
double Xcor[2][6] = {0.};
double event_theta[2] = {0.};
double event_phi[2] = {0.};
double x_value[2] = {0.};

get_cross_correlation(1, tdoa, Xcor);
reconstruction(1, event_theta, event_phi, x_value);
    
gStyle->SetPalette(55);
h2->Draw("colorz");
h2->GetZaxis()->SetRangeUser(-3, 5.5);



delete Tree_Muon;
file->Close();
delete file;


}




int filter_N = 1500;
Int_t filter_n_size = filter_N;
TVirtualFFT *filter_fft_forward = TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
TVirtualFFT *filter_fft_back = TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");
Double_t * fft_back_x = new Double_t[1500];

int get_FFT(int ch, double x[1500], double FW_x[1500], double filtered_x[1500], double FW_re[751], double FW_im[751], double deep_filtered_x[1500]){


   filter_fft_forward->SetPoints(x);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(FW_re, FW_im);
   
   FW_re[0] = 0.;
   FW_im[0] = 0.;
   filter_fft_back->SetPointsComplex(FW_re, FW_im);
   filter_fft_back->Transform();
   fft_back_x = filter_fft_back->GetPointsReal();

	for(int i=0;i<1500;i++){
		x[i] = fft_back_x[i]/1500.;
		//cout<<x[i]<<endl;
  						}
  						
  						
  						   
	double data[4096] = {0.};
	//for(int i=0; i<1500; i++){
	//	FW_x[i] = (0.42-0.5*cos(2*PI*(i)/1500.)+0.08*cos(4*PI*(i)/1500.))*x[i];
	//}
	Planck_taper_windowing(0.2, 1500, x, data);
	
	int pos_theta = (int) round(17.2) - 7;
	int pos_phi = (int) round((-39.7 - (-25.))/6.) + 7;
	//deconvolute_response(data, ch, 1, pos_theta, pos_phi);
	

	for(int i=0; i<1500; i++){
		FW_x[i] = data[i];
		//cout<<data[i]<<endl;
	}	
	

   filter_fft_forward->SetPoints(FW_x);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(FW_re, FW_im);



 for (Int_t k=0; k<751; k++){//0 means DC component // from 185-350 MHz
       if((k<222)||(k>420)) {FW_re[k] = 0.; FW_im[k] = 0;}
                             }



	
  	//filter_fft_back->SetPointsComplex(FW_re, FW_im);
  	//filter_fft_back->Transform();
  	//fft_back_x = filter_fft_back->GetPointsReal();

	//for(int i=0;i<1500;i++){
	//	filtered_x[i] = fft_back_x[i]/1500.;
  	//					}		  


   //filter_fft_forward->SetPoints(filtered_x);
   //filter_fft_forward->Transform();
   //filter_fft_forward->GetPointsComplex(FW_re, FW_im);
   

  	   
  	   double power[751] = {0.};
  	   for(Int_t k=0; k<751; k++){
  	   		power[k] = FW_re[k] * FW_re[k] + FW_im[k] * FW_im[k];
  	   		}



  	   	// caution!! sort function will change value in power array!
  	   	double median_1 = TMath::Median(73, &power[240]);
  	   	
  	   	cout<<median_1<<endl;

  	   	double median_2 = TMath::Median(61, &power[335]);
  	   	
  	   	cout<<median_2<<endl;
  	   	
  	   
  	double deep_re[751] = {0.};  
  	double deep_im[751] = {0.}; 
 	for(Int_t k=0; k<751; k++){//0 means DC component // from 185-350 MHz
    	deep_re[k] = FW_re[k];
    	deep_im[k] = FW_im[k];
                             }

	
	for(Int_t k=0; k<324; k++){
		if(power[k]/median_1>10){
			//cout<<k<<endl;
			deep_re[k-1] = 0.;
			deep_im[k-1] = 0.;
			deep_re[k] = 0.;
			deep_im[k] = 0.;
			deep_re[k+1] = 0.;
			deep_im[k+1] = 0.;
			
			if(power[k-2]/median_1>5){
				deep_re[k-2] = 0.;
				deep_im[k-2] = 0.;
			}
			if(power[k+2]/median_1>5){
				deep_re[k+2] = 0.;
				deep_im[k+2] = 0.;
			}
			if(power[k-3]/median_1>5){
				deep_re[k-3] = 0.;
				deep_im[k-3] = 0.;
			}
			if(power[k+3]/median_1>5){
				deep_re[k+3] = 0.;
				deep_im[k+3] = 0.;
			}
		}
	
	}

	for(Int_t k=324; k<751; k++){
		if(power[k]/median_2>10){
			//cout<<k<<endl;
			deep_re[k-1] = 0.;
			deep_im[k-1] = 0.;
			deep_re[k] = 0.;
			deep_im[k] = 0.;
			deep_re[k+1] = 0.;
			deep_im[k+1] = 0.;
			
			if(power[k-2]/median_2>5){
				deep_re[k-2] = 0.;
				deep_im[k-2] = 0.;
			}
			if(power[k+2]/median_2>5){
				deep_re[k+2] = 0.;
				deep_im[k+2] = 0.;
			}
			if(power[k-3]/median_2>5){
				deep_re[k-3] = 0.;
				deep_im[k-3] = 0.;
			}
			if(power[k+3]/median_2>5){
				deep_re[k+3] = 0.;
				deep_im[k+3] = 0.;
			}
		}
	
	}
	


	//cout<<"here:  "<<deep_re[328]<<endl;
	//deep_re[328] = 0.;
	for(Int_t k=0; k<751; k++){
    	FW_re[k] = deep_re[k];
    	FW_im[k] = deep_im[k];
    }	   
  	   
  	filter_fft_back->SetPointsComplex(deep_re, deep_im);
  	filter_fft_back->Transform();
  	fft_back_x = filter_fft_back->GetPointsReal();

	for(int i=0;i<1500;i++){
		deep_filtered_x[i] = fft_back_x[i]/1500.;
  						}  	   
  	   
  	   
   int Max_i = 0;
   double Max = 0.;

   for(int i=200;i<1300;i++){
	   if(fabs(deep_filtered_x[i])>Max) {
		   Max = fabs(deep_filtered_x[i]);  Max_i = i; //cout<<"i: "<<i<<"  max:  "<<Max<<endl;
		   }
  	   }  	   
 	   
  	   
   return Max_i;
   
   
   
}




int filter_time(double filtered_x[1500], double SW_x[1500], int Max_i){

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






	
Double_t ** Correlator_time = new_Array2D<double>(12, 32768);

int N = 32768;
//TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
TVirtualFFT *fft_back[12] = {TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), 
TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), 
TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K")};

	
//result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
void get_cross_correlation(int pol, double tdoa[2][6], double Xcor[2][6]){	



   	int i = 0;
    int j = 2;
    
   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 0+pol; j = 2+pol;
        break;
      case 1:
        i = 0+pol; j = 4+pol;
        break;
      case 2:
        i = 0+pol; j = 6+pol;
        break;
      case 3:
        i = 2+pol; j = 4+pol;
        break;
      case 4:
        i = 2+pol; j = 6+pol;
        break;
      case 5:
        i = 4+pol; j = 6+pol;
        break;

    }

       // full or short waveform
       //for (Int_t k=1; k<HALF_WINDOW_WIDTH; k++){//0 means DC component
       //     Correlator_re[pol][pair][k] = (SW_re[i][k] * SW_re[j][k] + SW_im[i][k] * SW_im[j][k]);            
       //     Correlator_im[pol][pair][k] = (SW_re[i][k] * SW_im[j][k] - SW_im[i][k] * SW_re[j][k]);
       //}//k 

       for (Int_t k=1; k<751; k++){//0 means DC component
            Correlator_re[pol][pair][k] = (FW_re[i][k] * FW_re[j][k] + FW_im[i][k] * FW_im[j][k]);            
            Correlator_im[pol][pair][k] = (FW_re[i][k] * FW_im[j][k] - FW_im[i][k] * FW_re[j][k]);
       }//k 
                                                           
    tdoa[pol][pair] = real_cross_correlation(pol, pair, i, j);

     }// pair loop


}
   


double real_cross_correlation(int pol, int pair, int i, int j){

        double Fine_re[16385] = {0};            
        double Fine_im[16385] = {0};

        // full or short waveform
        //for (Int_t k=1; k<HALF_WINDOW_WIDTH; k++){//0 means DC component
        for(Int_t k=1; k<751; k++){//0 means DC component
            Fine_re[k] = Correlator_re[pol][pair][k]; 
                   
            Fine_im[k] = Correlator_im[pol][pair][k];
            //cout<<ant_re[i][k]<<endl;
            
            //if((pol==0)&&(pair==5)&&(k==100)) cout<<Fine_re[k]<<endl;
                
        }//k 


          int which_fft = pol*6+pair;

          fft_back[which_fft]->SetPointsComplex(Fine_re, Fine_im);
          fft_back[which_fft]->Transform();
          Correlator_time[which_fft] = fft_back[which_fft]->GetPointsReal();


//if((pol==0)&&(pair==5)) cout<<Correlator_time[0]<<endl;


         // full or short waveform
         //double factor = sqrt(SUM_x[i])*sqrt(SUM_x[j])*WINDOW_SPACING;
         double factor = sqrt(SUM_x[i])*sqrt(SUM_x[j])*1500.;
         double max = 0;
         int mm;

                       for (Int_t k=0; k<32768; k++){//32768
                        //cout<<Correlator_time[k]<<endl;
                correlation_pair[pol][pair][k] = Correlator_time[which_fft][k]/factor;
                if (correlation_pair[pol][pair][k] > max){max = correlation_pair[pol][pair][k]; mm = k;}

                                                    }//k 

//cout<<max<<"  "<<mm<<endl;

if(mm>16383) {mm -= 32768;}
//if(mm>149) {mm -= 300;}
cout<<max<<endl;
return mm*(0.8*WINDOW_SPACING/32768.);

}	
	


int reconstruction(int pol, double event_theta[2], double event_phi[2], double x_value[2]){

double R = 10000;

double inter_max[2] = {0};
   //0 for H-pol 1 for V-pol
   //for(int pol = 0; pol<2; pol++){
      
   	int i = 0;
    int j = 2;
    
//double inter_map[56][100] = {0};//theta: -6~50; phi: -76~24
Double_t ** inter_map = new_Array2D<double>(1800, 3600);

     double max[3] = {0};
     int max_theta[3] = {0};
     int max_phi[3] = {0};
  
   for(int pair = 0; pair <6; pair++){

  //threads[pol*6+pair].join();

    switch(pair){
      case 0:
        i = 0+pol; j = 2+pol;
        break;
      case 1:
        i = 0+pol; j = 4+pol;
        break;
      case 2:
        i = 0+pol; j = 6+pol;
        break;
      case 3:
        i = 2+pol; j = 4+pol;
        break;
      case 4:
        i = 2+pol; j = 6+pol;
        break;
      case 5:
        i = 4+pol; j = 6+pol;
        break;

    }



    for(int i_theta=0; i_theta<1800; i_theta++){
      for(int i_phi=0; i_phi<3600; i_phi++){
          double theta = (i_theta - 900.)/10.;//cyc
          double phi = i_phi/10.;
          double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          
          // full or short waveform
          //delta_t = delta_t - (peak_position[j] - peak_position[i])*0.8;
          //delta_t = delta_t/(0.8*WINDOW_SPACING/32768.);
          delta_t = delta_t/(0.8*1500./32768.);
          if(delta_t<0) delta_t += 32768;
          
          if(delta_t<0||delta_t>32768) continue; //cout<<"delta_t error!  "<<delta_t<<endl;

          //cout<<delta_t<<endl;
          //cout<<correlation_pair[pol][pair][(int) delta_t]<<endl;
          inter_map[i_theta][i_phi] += correlation_pair[pol][pair][(int) round(delta_t)];
          //inter_map[theta+90][phi+90] += theta+phi;
          int correct_phi = i_phi + 260;
    	  if(correct_phi>1800) {correct_phi -= 3600;}
    	  correct_phi += 1800;
          h2->SetBinContent(correct_phi, i_theta, inter_map[i_theta][i_phi]);
        }
      }
      
    if(pair == 5) {


     for (int ci=0; ci<1800; ci++){
      for (int cj=0; cj<3600; cj++){
          if(inter_map[ci][cj] > max[0]){
          		if((abs(ci-max_theta[0])<=2)&&(abs(cj-max_phi[0])<=1)){
          			max[0] = inter_map[ci][cj]; max_theta[0] = ci; max_phi[0] = cj;
          		}
          		else if((abs(ci-max_theta[1])<=2)&&(abs(cj-max_phi[1])<=1)){
          			max[1] = max[0]; max_theta[1] = max_theta[0]; max_phi[1] = max_phi[0];
          			max[0] = inter_map[ci][cj]; max_theta[0] = ci; max_phi[0] = cj;
          		}
          		else{
          			max[2] = max[1]; max_theta[2] = max_theta[1]; max_phi[2] = max_phi[1];
          			max[1] = max[0]; max_theta[1] = max_theta[0]; max_phi[1] = max_phi[0];
          			max[0] = inter_map[ci][cj]; max_theta[0] = ci; max_phi[0] = cj;
          		}
          		//cout<<a<<"   "<<max_theta[0]<<"   "<<max_theta[1]<<"    "<<max_phi[0]<<"   "<<max_phi[1]<<endl;
          }
     }
   }

     }
     
     	

	}//pair loop


if(max[0]==0){ 
	event_theta[pol] = -1000;
	event_phi[pol] = -1000;
	x_value[pol] = -6000;	
	return false;
	}




//////////////////////////////////////////////////////////////////////
//fine search
/*
double fine_inter_map[3][40][20] = {0};
int fine_max_theta = 0;
int fine_max_phi = 0;
double fine_max = 0;
int max_w = 0;

   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 0+pol; j = 2+pol;
        break;
      case 1:
        i = 0+pol; j = 4+pol;
        break;
      case 2:
        i = 0+pol; j = 6+pol;
        break;
      case 3:
        i = 2+pol; j = 4+pol;
        break;
      case 4:
        i = 2+pol; j = 6+pol;
        break;
      case 5:
        i = 4+pol; j = 6+pol;
        break;

    }



	for(int w = 0; w < 3; w++){
    	for(int i_theta=0; i_theta<40; i_theta++){
      	for(int i_phi=0; i_phi<20; i_phi++){
          double theta = -2 + i_theta/10. + max_theta[w]-6;//cyc
          double phi = -1 + i_phi/10. + max_phi[w];
          double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          
          delta_t = delta_t - (peak_position[j] - peak_position[i])*0.8;
          delta_t = delta_t/(0.8*WINDOW_SPACING/32768.);
          if(delta_t<0) delta_t += 32768;
          if(delta_t<0||delta_t>32768) continue; 
          fine_inter_map[w][i_theta][i_phi] += correlation_pair[pol][pair][(int) round(delta_t)];
        }
      }
      
    if(pair == 5) {


     for (int ci=0; ci<40; ci++){
      for (int cj=0; cj<20; cj++){
          if(fine_inter_map[w][ci][cj] > fine_max){
          			fine_max = fine_inter_map[w][ci][cj]; fine_max_theta = ci; fine_max_phi = cj; max_w = w;
          }
       }
     }
     }
     }//w loop  
     	

	}//pair loop
	
     inter_max[pol] = fine_max;
     
     event_theta[pol] = -2 + fine_max_theta/10. + max_theta[max_w]-6;//cyc
     event_phi[pol] = -1 + fine_max_phi/10. + max_phi[max_w];
     x_value[pol] = inter_max[pol];
*/
	cout<<(max_theta[0]-900)/10.<<"  "<<(max_phi[0])/10.<<"   "<<max[0]/6.<<endl;
	//}//pol loop
	return true;
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
