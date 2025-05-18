#include "TTree.h"
#include <iostream>
#include "TGraph.h" 
#include "TVirtualFFT.h" 
#include "TFile.h"
#include <fstream>
#include <math.h>
using namespace std;
 
 
 
bool Is_pulser();
void get_template(double &Sum_template, double Re_template[751], double Im_template[751]); 
double get_cross_correlation(double x[1500], double Sum_template, double Re_template[751], double Im_template[751], int &x_peak);
int filter_time(double x[1500]);
void filter_fre(double x[1500]);
 
 Double_t x[8][4096]={0};//8 antennas
 
 
//int main(int argc, char *argv[]) {
 void get_trigger(){
 
 
 TH1F * abc = new TH1F("abc","abc",7,35,70);
 
 
   double Drone_time[10000] = {0};
 double Total_attenuation[10000] = {0};
 double Drone_azimuth[10000] = {0};
 double Drone_elevation[10000] = {0};
 double relative_x[10000] = {0};
 double relative_y[10000] = {0};
 double relative_z[10000] = {0};
 
 float _drn_time = 0;
 float _tot_att = 0;
 float _drn_azimuth = 0;
 float _drn_elevation = 0;
 float _relative_x = 0;
 float _relative_y = 0;
 float _relative_z = 0;
 
   FILE *fpin = fopen("/home/cyc/programming/TAROGE_4_pulser_good_timing/2nd_theodolite_origin.csv","r"); 
  for(int i=0; i<1989;i++){

   int p = fscanf(fpin,"%*f %f %*f %f %f %*f %*f %*f %f %f %f %*f %f", &_drn_time, &_drn_azimuth, &_drn_elevation, &_relative_x, &_relative_y, &_relative_z, &_tot_att);
   
 //18:18:13.450006000//run 2113 event 0 time
 Drone_time[i] = (18*60*60+18*60+14.+_drn_time)-(18*60*60+18*60+13.45) +1.24;// 
 
 Total_attenuation[i] = _tot_att/1.;
 Drone_azimuth[i] = _drn_azimuth;
 Drone_elevation[i] = _drn_elevation;
 relative_x[i] = _relative_x;
 relative_y[i] = _relative_y;
 relative_z[i] = _relative_z;
 

 if((i%2==1)&&(Drone_time[i]>90)&&(Drone_time[i]<710)) {abc->Fill(Total_attenuation[i]);}
 
 if(i%2) {Drone_time[i] += 0.499;}//for old DGPS data is + 0.999, for new drone data is + 0.499

                         }

  fclose(fpin);
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 TFile *ntuple_file = new TFile("ntuple.root","RECREATE");

TNtuple *ntuple = new TNtuple("ntuple","ntuple","run:event:cross_correlation:V_p:time:time_stamp:attenuation:flag:drn_azimuth:drn_elevation:rel_x:rel_y:rel_z",320000);
 
 
 int number1 = 2112;
 int number2 = 2114;
 
 
 
 double Sum_template = 0;
double Re_template[751] = {0};
double Im_template[751] = {0};

get_template(Sum_template, Re_template, Im_template);
 
 
 
 
 
 

 
 
 
//int number1 =  atoi (argv[1]);
//int number2 =  atoi (argv[2]);
 
  bool bGood;
  Int_t timeStamp_FPGA; 

  Int_t event_number;
  const char * fname_evt;
  TFile *file;
  TTree *Tree_Muon;
int Event_N = 0;
double eventTime = 0;
double eventTime_start = 0;

double time_prev = 0;
double time_accum = 0;

double timestamp[100000]={0};

double event[100000]={0};


int pulser_event_count=0;


double cal_timestamp[100000]={0};
double cal_event[100000]={0};
int cal_Event_N = 0;

int flag_check[10000] = {0};

  Double_t t[1500];
  Short_t T1H[1500];
  Short_t T2H[1500];
  Short_t T3H[1500];
  Short_t T4H[1500];
  Short_t T1V[1500];
  Short_t T2V[1500];
  Short_t T3V[1500];
  Short_t T4V[1500];
  



     TCanvas * multi_waveform = new TCanvas("multi_waveform","multi_waveform",0,0,1800,1100);
     multi_waveform -> Divide(2,4);
     multi_waveform -> Draw();
     
     TGraph * draw_waveform;


  for(int number = number1; number<number2; number++){
  

  
    fname_evt =  Form( "/home/cyc/programming/TAROGE_4_pulser_good_timing/run%08d.root",number);
    file = new TFile(fname_evt);
    Tree_Muon = (TTree*) file->Get("t");
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
		  
 event_number = Tree_Muon->GetEntries();

     for (int j=0;j<event_number;j++){

      Tree_Muon->GetEntry(j);
      
       if(number==number1&&j==0) {eventTime_start = eventTime;}
      
//cout<<"good?: "<<bGood<<endl;
if(bGood!=1) continue;
//cout<<timeStamp_FPGA<<endl;
     double time = timeStamp_FPGA * 3.2e-7;// 16*20ns,see repeaterR in verilog code;
 

      
      if(time<time_prev) time_accum += 1;
      
      time_prev = time;
     
     timestamp[Event_N] = fmod(time, 0.1);//(((int)(time*100000))%((int)(10000)))/((double) (100000));
     

     
//     cout<<timestamp[Event_N]<<endl;
     
     // cannot use this since tune_threshold will comsume more than 1 second
     //event[Event_N] = time_accum + time;//eventTime - eventTime_start;
     event[Event_N] = eventTime - eventTime_start;
//    cout<<"fuck:  "<<time_accum + time - (eventTime - eventTime_start)<<endl;

     double timestamp_fit = (0.05775-0.053)/(90.-619.3)*(event[Event_N]-90.) + 0.05775;
     double timestamp_low = timestamp_fit - 0.005;//0.005
     double timestamp_high = timestamp_fit + 0.001;//0.0002;
   
   //&&(time>0.3)&&(time<0.4)  
     if((89<event[Event_N])&&(event[Event_N]<720)&&(timestamp_low<timestamp[Event_N])&&(timestamp[Event_N]<timestamp_high)){
     //if(true){
    //89 720 
     //if((89<event[Event_N])&&(event[Event_N]<720)){
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


double original_x[8][1500] = {0};
           for (int k=0;k<8;k++){
           for (int i=0;i<1500;i++)
	{
  	  original_x[k][i] = x[k][i];
	}
                                 }


bool _Is_pulser = Is_pulser();

int x_peak[4] = {0};

double cro_1H = get_cross_correlation(x[0], Sum_template, Re_template, Im_template, x_peak[0]);
double cro_2H = get_cross_correlation(x[2], Sum_template, Re_template, Im_template, x_peak[1]);
double cro_3H = get_cross_correlation(x[4], Sum_template, Re_template, Im_template, x_peak[2]);
double cro_4H = get_cross_correlation(x[6], Sum_template, Re_template, Im_template, x_peak[3]);


int delta_10 = fabs(x_peak[1] - x_peak[0]);
int delta_20 = fabs(x_peak[2] - x_peak[0]);
int delta_30 = fabs(x_peak[3] - x_peak[0]);

double ave_cro_H = (cro_1H + cro_2H + cro_3H + cro_4H)/4.;

double ave_y_peak = (fabs(original_x[0][x_peak[0]]) + fabs(original_x[2][x_peak[1]]) + fabs(original_x[4][x_peak[2]]))/3.;

double _sum = 0;
for(int k=0;k<1500;k++){if((k<(x_peak[2]-50))||(k>(x_peak[2]+50))) _sum += pow(original_x[4][k],2);}

cout<<sqrt(_sum/1500.)<<endl;

//cout<<x[0][x_peak[0]]<<"  "<<x[2][x_peak[1]]<<"  "<<x[4][x_peak[2]]<<"  "<<x[6][x_peak[3]]<<endl;

_Is_pulser = false;
if((ave_cro_H>0.5)&&(delta_10<200)&&(delta_20<200)&&(delta_30<200))
//if(ave_cro_H>0.5)
//if(true)
				  {_Is_pulser = true; 
				   //cout<<cro_1H<<"  "<<cro_2H<<"  "<<cro_3H<<"  "<<endl;
                   //cout<<timestamp[Event_N]<<"  "<<time<<endl;
                   //cout<<delta_10<<" "<<delta_20<<" "<<delta_30<<endl;
                   cal_timestamp[cal_Event_N] = timestamp[Event_N];
                   cal_event[cal_Event_N] = event[Event_N];
                   cal_Event_N++;
                  }

   if(_Is_pulser){
   		/*
    	if((timestamp[Event_N]>0.0545)&&(timestamp[Event_N]<0.0554)&&(event[Event_N]<400)){
			for(int k=0; k<8;k++){
				multi_waveform->cd(k+1);
				draw_waveform = new TGraph(1500,t,x[k]);
				draw_waveform->Draw("APL");
				draw_waveform->GetYaxis()->SetRangeUser(-200,200);
				multi_waveform->Update();
	                     }
	                     
			usleep(1000000);	
		}
		*/             
	                     
	                     
//     cout<<number<<"  "<<j<<endl;

	//sleep(1);

double Attenuation = 0;
int flag = 1;
float drn_azi;
float drn_ele;
float rel_x;
float rel_y;
float rel_z;

  for(int i=0; i<1989;){

      if((Drone_time[i]<event[Event_N])&&(event[Event_N]<Drone_time[i+1]))
          {
            Attenuation = Total_attenuation[i];
            drn_azi = Drone_azimuth[i];
            drn_ele = Drone_elevation[i];
            flag_check[pulser_event_count] = i;
            rel_x = relative_x[i];
            rel_y = relative_y[i];
            rel_z = relative_z[i];
            //event[Event_N] += 18.;//request from CY, LEAP second between GPS and UTC
            break;
          }
 
  i+=2;
                      }

if(flag_check[pulser_event_count]==flag_check[pulser_event_count-1]) {flag=0;};

//"run:event:cross_correlation:V_p:time:time_stamp:attenuation:flag:drn_azimuth:drn_elevation"
cout<<number<<"  "<<j<<endl;
     ntuple->Fill(number,j,ave_cro_H,ave_y_peak,event[Event_N],time,Attenuation,flag,drn_azi,drn_ele,rel_x,rel_y,rel_z);
			
     pulser_event_count++;													
          	    }
     		
     																					}


     
     
     Event_N++;
                                  }// j loop
                                                      
                                                      }// number loop 
                                                      
     TGraph * cyc = new TGraph(cal_Event_N,cal_event,cal_timestamp);
     cyc->SetMarkerStyle(7);
     cyc->SetMarkerColor(2);
     cyc->GetYaxis()->SetRangeUser(0,0.1);
     cyc->Draw("AP");

     
     //60 to 660 for pulser run
     double fit_time[2] = {90,619.3};
     double fit_stamp[2] = {0.05775,0.053};
     TGraph * wyn = new TGraph(2,fit_time,fit_stamp);
     wyn->Draw("PL same");
     


//     multi_waveform->Draw();
     cout<<pulser_event_count<<endl;
 
 
 
 
 
      TGraph * hxj = new TGraph(1989,Drone_time,Total_attenuation);
     hxj->SetMarkerStyle(6);
     hxj->SetMarkerColor(4);
     hxj->SetLineColor(3);
     //hxj->GetYaxis()->SetRangeUser(0,0.1);
     hxj->Draw("PL same");
     
//abc->Draw();

     
   ntuple_file->cd();
   ntuple->Write();
   ntuple_file->Close();
   delete ntuple_file;
 
              }// end of main function
              
              
              
              
/*              
              void get_sub_band_waveform(int Q)  
{ 



  fft_own->SetPoints(x[Q]);
  fft_own->Transform();


  fft_own->GetPointsComplex(re_full,im_full);



    for(int i=1;i<2049;i++)
      {
	re_band[j][i]= 2*re_full[i]*bandfactor[j][i];

	im_band[j][i]= 2*im_full[i]*bandfactor[j][i];
      }


  
  fft_back1->SetPointsComplex(re_band[0],im_band[0]);
  fft_back1->Transform();
  fft_back1->GetPointsComplex(v_tmp_re[0],v_tmp_im[0]);


}

*/


bool Is_pulser()
{

bool Is_pulser = false;

double Max[8]={0};
int Max_i[8]={0};
double Sum[8]={0};
double _Rms[8]={0};
double _Snr[8]={0};

           for (int i=0;i<1500;i++)
	{
            for(int j=0; j<8; j++){
	  if(fabs(x[j][i])>Max[j]) {Max[j] = fabs(x[j][i]);  Max_i[j] = i;} 	  
 	                               }
	}
	
	
           for (int i=0;i<1500;i++)
	{
            for(int j=0; j<8; j++){
   	        if((i<Max_i[j]-50)||(i>Max_i[j]+50))
 	  Sum[j] += pow(x[j][i],2);
 	                               }
	}	
	
	
	 for(int j=0; j<8; j++){
	_Rms[j] = sqrt(Sum[j]/1400.);
	
	_Snr[j] = Max[j]/_Rms[j];
						  }
						  
	double Ave_H_Max = (Max[0]+Max[2]+Max[4])/3.;
	double Ave_V_Max = (Max[1]+Max[3]+Max[5])/3.;
	
	double Ave_H_Snr = (_Snr[0]+_Snr[2]+_Snr[4])/3.;
	
//if((Ave_H_Max>1*Ave_V_Max)&&(Ave_H_Snr>1)) Is_pulser=true;
if((Ave_H_Snr>1.5)) Is_pulser=true;
Is_pulser=true;
return Is_pulser;

}
	
	
	
	

double get_cross_correlation(double x[1500], double Sum_template, double Re_template[751], double Im_template[751], int &x_peak){	

filter_fre(x);
x_peak = filter_time(x);

int N = 1500;

   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C");
      
   TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
   
   	
double Re_ant[751];
double Im_ant[751];



   fft_forward->SetPoints(x);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(Re_ant,Im_ant);

double Correlator_re[751] = {0};            
double Correlator_im[751] = {0};

          	   for (Int_t k=1; k<751; k++){//0 means DC component
    Correlator_re[k] = (Re_template[k] * Re_ant[k] + Im_template[k] * Im_ant[k])/1500.; 
           
    Correlator_im[k] = (Re_template[k] * Im_ant[k] - Im_template[k] * Re_ant[k])/1500.;
        
                                            }//k  
  
          

  
double cro_max = 0;

Double_t * Correlator_time = new Double_t [1500];

  fft_back->SetPointsComplex(Correlator_re,Correlator_im);
  fft_back->Transform();
  Correlator_time = fft_back->GetPointsReal();
  

     
     cro_max = Correlator_time[0];
     
    for (int k=1; k<1500; k++){
  
        if(Correlator_time[k] > cro_max) {cro_max = Correlator_time[k];}
   		                         }

 double Sum_x = 0;
  for(int i=0;i<1500;i++){Sum_x += x[i]*x[i];}

 double factor = sqrt(Sum_template)*sqrt(Sum_x);
 
// cout<<cro_max/factor<<endl;
 
     cro_max /= factor;


   delete fft_forward;
   delete fft_back;
   
//   delete[] Correlator_time; // when delete fft_back, it already delete this
   	
	return cro_max;
}
	
	
	


void get_template(double &Sum_template, double Re_template[751], double Im_template[751]){	

 int N = 1500;
 double x[1500] = {0};
 Short_t T2V[1500];
   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C");



   TFile * file_template = new TFile("/home/cyc/programming/TAROGE_4_pulser_good_timing/run00002101.root");
   TTree * Tree_template = (TTree*) file_template->Get("t");

		  Tree_template->SetBranchAddress("T2V",T2V);

		  


      Tree_template->GetEntry(1);

Sum_template=0;

           for (int i=0;i<1500;i++) {x[i]= (float)T2V[i]*500/32512; }
           
	filter_fre(x);	  
	int x_peak = filter_time(x);  

           for (int i=0;i<1500;i++) {Sum_template += x[i]*x[i]; }

   fft_forward->SetPoints(x);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(Re_template,Im_template);
 
   
   delete fft_forward;
   
   //Tree_Muon->Close();
   delete Tree_template;
   file_template->Close();
   delete file_template;
}	



int filter_time(double x[1500]){

int Max_i = 0;
double Max = 0.;

           for (int i=0;i<1500;i++)
	{

	  if(fabs(x[i])>Max) {Max = fabs(x[i]);  Max_i = i;} 	  

	}
	//cout<<Max_i<<endl;
	
           for (int i=0;i<1500;i++)
	{

   	       x[i] = ((i<Max_i-100)||(i>Max_i+100))?0:x[i];
 	                               
	}

return Max_i;

}




void filter_fre(double x[1500]){

int N = 1500;

   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C");
      
   TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
   
   	
double Re_ant[751];
double Im_ant[751];



   fft_forward->SetPoints(x);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(Re_ant,Im_ant);


 for (Int_t k=0; k<751; k++){//0 means DC component // from 185-350 MHz
       if((k<222)||(k>420)) {Re_ant[k] = 0.; Im_ant[k] = 0;}
                             }
  
          

  


Double_t * filtered_x = new Double_t [1500];

  fft_back->SetPointsComplex(Re_ant,Im_ant);
  fft_back->Transform();
  filtered_x = fft_back->GetPointsReal();

for(int i=0;i<1500;i++){
x[i] = filtered_x[i]/1500.;
  						}


}


