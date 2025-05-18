//20210526, change the antenna coordinates and delay order, from T1H T2H....T3V T4V -->> T1H T1V....T4H T4V


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
using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8
 

 
void get_cross_correlation(int nt_j, double tdoa[2][6], double Xcor[2][6]);
int filter_time(double x[4096], double new_x[4096], int Max_i);
int filter_fre(double x[4096], double &power_sum);
Double_t expectedTimeDiff_angle(double theta, double phi, double *xyz);
Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j);
Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j);

void get_simple_Xcor(double tdoa[2][6], double Xcor[2][6], double saSNR[2]);
Double_t * filtered_x = new Double_t [1500];
 
 Double_t x[8][4096]={0};//8 antennas
 Double_t new_x[8][4096]={0};//8 antennas
  double ant_re[8][2049] = {0};
  double ant_im[8][2049] = {0};
  


 double correlation_pair[2][6][32768] = {0};
 
 
double Xcor_tdoa[6][10000]={0};
double Xcor_max[10000][6]={0};
double Sum_x[8] = {0};
double Correlator_re[2][6][151] = {0};            
double Correlator_im[2][6][151] = {0};
int peak_position[8] = {0};

double Expect_Xcor_tdoa[6][20000]={0};

int event_record[20000][2] = {0};// run, entry

 
int main(int argc, char *argv[]) {

string date_name = argv[1];
string run_list = argv[2];

cout<<endl<<date_name<<endl;
cout<<run_list<<endl;

//void Real_CRselect_tchain(string date_name, string run_list){

int which_date = std::stoi(date_name);

TChain *tr_event = new TChain("t");
ifstream inFile(run_list);
	string str;
	int fNFile = 0;
	vector<Long64_t> fvNEntries;
	vector<int> run_number;
	
	int which_run = 0;
	int which_event = 0;

	while(getline(inFile, str)){

		cout<<"reading: "<< str <<endl;
		if( str.find(".root") != std::string::npos)
		{
		    int nPos_1 = str.find_last_of("n");
            int nPos_2 = str.find_last_of(".");
            string run_num="";
            run_num = str.substr(nPos_1+1,nPos_2-nPos_1-1);
            run_number.push_back(std::stoi(run_num));
            
			fNFile +=  tr_event->AddFile(str.c_str());

			fvNEntries.push_back(tr_event->GetEntries());
		}


	}

	inFile.close();
	
	//cout<<fNFile<<endl;
	//for(auto context : fvNEntries){
		//cout<<context<<endl;
		//}
	//for(auto context : run_number){
		//cout<<context<<endl;
		//}




  const   char * out_filename;

  out_filename =  Form("./candidate_events/%08d.txt", which_date);
  fstream out_put;
  out_put.open(out_filename);// without ios::out   clean the file before write




//base center   2021 photogrametry  Tower-4 is back side
double coor_ant[8][3] = {{328914.252, 2697625.953, 709.389},{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949},{328936.041, 2697639.685, 704.949}};

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

  double delay_ant[8] = {0,-92.4,53.2,-30.3,84.6,59.5,-1523.1,-743.1};
 

 
  


TH2F *h2 = new TH2F("h2","h2", 25, -26.2, -23.7, 110, -6, 5);
//TCanvas *c4 = new TCanvas("c4","Canvas Example",1500,1000);

TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
c1->Divide(3,2);

TCanvas *c2 = new TCanvas("c2","Canvas Example",1500,1000);
c2->Divide(2,1);

//TCanvas *c3 = new TCanvas("c3","Canvas Example",1500,1000);

//c3->Divide(1,3);

float nt_run;
float nt_event;
float nt_time;
float drn_azimuth;
float drn_elevation;




double event_time[20000] = {0};
double event_theta[2][20000] = {0};
double event_phi[2][20000] = {0};
double delta_HV_theta[20000] = {0};
double delta_HV_phi[20000] = {0};



int nt_j = 0;

 
 double Sum_template = 0;
double Re_template[2049] = {0};
double Im_template[2049] = {0};

 

  bool bGood;
  Int_t timeStamp_FPGA; 

  Int_t event_number;

int Event_N = 0;
double eventTime = 0;
double eventTime_start = 0;

double time_prev = 0;
double time_accum = 0;

double timestamp[20000]={0};

double event[20000]={0};


int pulser_event_count=0;




double cal_timestamp[20000]={0};
double cal_event[20000]={0};
int cal_Event_N = 0;

int flag_check[20000] = {0};

  Double_t t[1500];
  Short_t T1H[1500];
  Short_t T2H[1500];
  Short_t T3H[1500];
  Short_t T4H[1500];
  Short_t T1V[1500];
  Short_t T2V[1500];
  Short_t T3V[1500];
  Short_t T4V[1500];
  



int count_first_cut = 0;

for(int rot_deg = 0; rot_deg<1; rot_deg++){



 int number1 = 19948;
 int number2 = 20007;

 //int number1 = 16955;
 //int number2 = 16957;

double rotate_degree = 0;
nt_j=0;

          tr_event->SetBranchAddress("bGood",&bGood);
		  tr_event->SetBranchAddress("timeStamp_FPGA",&timeStamp_FPGA);
		  tr_event->SetBranchAddress("eventTime", &eventTime);
		  tr_event->SetBranchAddress("T1H",T1H);
		  tr_event->SetBranchAddress("T2H",T2H);
		  tr_event->SetBranchAddress("T3H",T3H);
		  tr_event->SetBranchAddress("T4H",T4H);
		  tr_event->SetBranchAddress("T1V",T1V);
		  tr_event->SetBranchAddress("T2V",T2V);
		  tr_event->SetBranchAddress("T3V",T3V);
		  tr_event->SetBranchAddress("T4V",T4V);
		  
 event_number = tr_event->GetEntries();
 
 cout<<event_number<<endl;
 
for (int entry = 0; entry<event_number; entry++){



      tr_event->GetEntry(entry);

if(entry==0) eventTime_start = eventTime;

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

     event_time[nt_j] = eventTime - eventTime_start;
     

double R = 10000;
     
double tdoa[2][6] = {0};
double Xcor[2][6] = {0};
double saSNR[2] = {0};

get_simple_Xcor(tdoa, Xcor, saSNR);

//cout<<tdoa[0][1]<<"   "<<Xcor[0][1]<<endl;

double H_Xcor = 0;
double V_Xcor = 0;
	
for(int pair = 0; pair <6; pair++){
    H_Xcor += Xcor[0][pair]/6.;
    V_Xcor += Xcor[1][pair]/6.;
    }
    
    int flag = 0;
    if((H_Xcor>0.7)&&(V_Xcor<0.7)&&(saSNR[0]>3.)) {flag = 1;}
    if((H_Xcor<0.7)&&(V_Xcor>0.7)&&(saSNR[1]>3.)) {flag = 2;}
    if((H_Xcor>0.7)&&(V_Xcor>0.7)&&(saSNR[0]>3.)&&(saSNR[1]>3.)) {flag = 3;}
    
    if(flag>0) {count_first_cut++; //cout<<H_Xcor<<endl;
    }
    else {continue;}
    
	
//continue;

get_cross_correlation(nt_j, tdoa, Xcor);







double inter_max[2] = {0};
   //0 for H-pol 1 for V-pol
   for(int pol = 0; pol<2; pol++){
      
   	int i = 0;
    int j = 2;
    
double inter_map[56][100] = {0};//theta: -6~50; phi: -76~24

     double max[3] = {0};
     int max_theta[3] = {0};
     int max_phi[3] = {0};
  
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



    for(int i_theta=0; i_theta<56; i_theta++){
      for(int i_phi=0; i_phi<100; i_phi++){
          double theta = i_theta - 6.;
          double phi = -76 + i_phi;
          double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          
          delta_t = delta_t - (peak_position[j] - peak_position[i])*0.8;
          delta_t = delta_t/(0.8*300./32768.);
          if(delta_t<0) delta_t += 32768;
          if(delta_t<0||delta_t>32768) continue; //cout<<"delta_t error!  "<<delta_t<<endl;

          //cout<<delta_t<<endl;
          //cout<<correlation_pair[pol][pair][(int) delta_t]<<endl;
          inter_map[i_theta][i_phi] += correlation_pair[pol][pair][(int) round(delta_t)];
          //inter_map[theta+90][phi+90] += theta+phi;
          //h2->SetBinContent(theta+90, phi+90, correlation_pair[pol][pair][(int) delta_t]);
        }
      }
      
    if(pair == 5) {


     for (int ci=0; ci<56; ci++){
      for (int cj=0; cj<100; cj++){
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





//////////////////////////////////////////////////////////////////////
//fine search

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
          double theta = -2 + i_theta/10. + max_theta[w]-6;
          double phi = -1 + i_phi/10. + max_phi[w]-76;
          double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          
          delta_t = delta_t - (peak_position[j] - peak_position[i])*0.8;
          delta_t = delta_t/(0.8*300./32768.);
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
     
     event_theta[pol][nt_j] = -2 + fine_max_theta/10. + max_theta[max_w]-6;
     event_phi[pol][nt_j] = -1 + fine_max_phi/10. + max_phi[max_w]-76;

	
	}//pol loop

delta_HV_theta[nt_j] = event_theta[0][nt_j] - event_theta[1][nt_j];
delta_HV_phi[nt_j] = event_phi[0][nt_j] - event_phi[1][nt_j];



//cout<<inter_max<<endl<<endl;
//&&(fabs(delta_HV_theta[nt_j])<10)&&(fabs(delta_HV_phi[nt_j])<10)
if((inter_max[0]>4.2)||(inter_max[1]>4.2)) {
//if((inter_max[0]>4.2)&&(inter_max[1]>4.2)&&(fabs(delta_HV_theta[nt_j])<6)&&(fabs(delta_HV_phi[nt_j])<6)) {
cout<<delta_HV_theta[nt_j]<<"    "<<delta_HV_phi[nt_j]<<endl;


	int run_N = fvNEntries.size();
	int run_pos = 0;
	for( ; run_pos<run_N; run_pos++){
		if(entry < fvNEntries.at(run_pos)) {break;}
	}
	 //cout<<run_pos<<endl;
	 which_run = run_number.at(run_pos);
	 
	 if(run_pos==0){
	 	which_event = entry;
	 	}
	 else{
	 	which_event = entry-fvNEntries.at(run_pos-1);
	 	}

event_record[nt_j][0] = which_run;
event_record[nt_j][1] = which_event;

nt_j++; 
//cout<<number<<"   "<<entry<<endl;
}
//nt_j++;

}// event loop

cout<<nt_j<<endl;

double time_difference[20000] = {0};
for(int i=0; i<nt_j-1; i++){
time_difference[i] = event_time[i+1] - event_time[i];
}
/*
//first event;
if(((event_time[1] - event_time[0]) > 10)||(nt_j==1)){//if only one event, store
	cout<<"candidate:  "<<event_record[0][0]<<"   "<<event_record[0][1]<<endl;
	cout<<event_time[0]<<"    "<<event_time[1]<<endl<<endl;
	//date run event thetaH phiH thetaV phiV
	out_put<<which_date<<"    "<<event_record[0][0]<<"    "<<event_record[0][1]<<"    "<<event_theta[0][0]<<"    "<<event_phi[0][0]<<"    "<<event_theta[1][0]<<"    "<<event_phi[1][0]<<endl;
}

for(int i=1; i<nt_j-1; i++){
	if((time_difference[i-1]>10)&&(time_difference[i]>10)&&(nt_j>2)){
		cout<<"candidate:  "<<event_record[i][0]<<"   "<<event_record[i][1]<<endl;
		cout<<event_time[i-1]<<"    "<<event_time[i]<<"    "<<event_time[i+1]<<endl<<endl;
		out_put<<which_date<<"    "<<event_record[i][0]<<"    "<<event_record[i][1]<<"    "<<event_theta[0][i]<<"    "<<event_phi[0][i]<<"    "<<event_theta[1][i]<<"    "<<event_phi[1][i]<<endl;
	}
}

//last event;
if(((event_time[nt_j-1] - event_time[nt_j-2]) > 10)&&(nt_j>1)){//(nt_j>1) exclude only one event.
	cout<<"candidate:  "<<event_record[nt_j-1][0]<<"   "<<event_record[nt_j-1][1]<<endl;
	cout<<event_time[nt_j-2]<<"    "<<event_time[nt_j-1]<<endl<<endl;
	out_put<<which_date<<"    "<<event_record[nt_j-1][0]<<"    "<<event_record[nt_j-1][1]<<"    "<<event_theta[0][nt_j-1]<<"    "<<event_phi[0][nt_j-1]<<"    "<<event_theta[1][nt_j-1]<<"    "<<event_phi[1][nt_j-1]<<endl;
}
*/

for(int i=0; i<nt_j; i++){
		out_put<<which_date<<"    "<<event_record[i][0]<<"    "<<event_record[i][1]<<"    "<<event_theta[0][i]<<"    "<<event_phi[0][i]<<"    "<<event_theta[1][i]<<"    "<<event_phi[1][i]<<endl;
}
   
/*
cout<<count_first_cut<<endl;

TGraph * theta_phi_map[2];

     for(int pol = 0; pol<2; pol++){
theta_phi_map[pol] = new TGraph(nt_j, event_phi[pol], event_theta[pol]);
theta_phi_map[pol]->SetMarkerStyle(kCircle);
}
theta_phi_map[0]->Draw("AP");
theta_phi_map[1]->Draw("P same");
theta_phi_map[0]->SetMarkerColor(kRed);
theta_phi_map[1]->SetMarkerColor(kBlue);
*/


     TGraph * cyc[2];
     TGraph * wyn[2];
     TGraph * hxj[2];
     
     for(int pol = 0; pol<2; pol++){
     
     c2->cd(pol+1);
     
     if(pol==0) {
     cyc[pol]= new TGraph(nt_j, event_theta[0], delta_HV_theta);

     
     
     cyc[pol]->Draw("AP");
     cyc[pol]->GetYaxis()->SetTitle("elevation");

     }
     else {
     cyc[pol]= new TGraph(nt_j, event_theta[0], delta_HV_phi);

     cyc[pol]->Draw("AP");
     cyc[pol]->GetYaxis()->SetTitle("azimuth");
     }

     cyc[pol]->SetMarkerStyle(kCircle);
     cyc[pol]->SetMarkerColor(kRed);
     cyc[pol]->GetXaxis()->SetTitle("time:s");
      
    
     }


}//rot_deg loop

  out_put.close(); 
    
              }// end of main function
              
              


Double_t * sXcor_time = new Double_t [300];

   int sXcor_N = 300;
   Int_t sXcor_size = sXcor_N+1;		
   TVirtualFFT *sXcor_fft_forward = TVirtualFFT::FFT(1, &sXcor_size,"R2C ES K");
   
   TVirtualFFT *sXcor_fft_backward = TVirtualFFT::FFT(1, &sXcor_N, "C2RBACKWARD M K");
   
void get_simple_Xcor(double tdoa[2][6], double Xcor[2][6], double saSNR[2]){

   double re[8][151] = {0};
   double im[8][151] = {0};

//double Sum_x[8] = {0};
double sSNR[8] = {0};
  
	for(int ant = 0; ant < 8; ant ++){
		Sum_x[ant] = 0;
	
    	int max_i = filter_fre(x[ant], sSNR[ant]);
    	peak_position[ant] = max_i;
    	filter_time(filtered_x, new_x[ant], max_i);  
    	

   sXcor_fft_forward->SetPoints(&new_x[ant][max_i-150]);
   sXcor_fft_forward->Transform();
   sXcor_fft_forward->GetPointsComplex(re[ant],im[ant]);

    	  for(int i = 0; i <300; i++){
          Sum_x[ant] += pow(new_x[ant][max_i-150+i], 2);
          
          }
   
   }
   
   saSNR[0] = (sSNR[0] + sSNR[2] +sSNR[4] +sSNR[6])/4.;
   saSNR[0] = (sSNR[1] + sSNR[3] +sSNR[5] +sSNR[7])/4.;
              


   //0 for H-pol 1 for V-pol
   for(int pol = 0; pol<2; pol++){
      
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
    








        //double Correlator_re[pol][pair][101] = {0};            
        //double Correlator_im[pol][pair][101] = {0};

                  	   for (Int_t k=1; k<151; k++){//0 means DC component
            Correlator_re[pol][pair][k] = (re[i][k] * re[j][k] + im[i][k] * im[j][k]); 
                   
            Correlator_im[pol][pair][k] = (re[i][k] * im[j][k] - im[i][k] * re[j][k]);
            //cout<<ant_re[i][k]<<endl;
            
            //if((pol==0)&&(pair==5)&&(k==100)) cout<<Correlator_re[pol][pair][k]<<endl;
                
                                                    }//k  


        //Double_t * Correlator_time = new Double_t [409600];

          sXcor_fft_backward->SetPointsComplex(Correlator_re[pol][pair],Correlator_im[pol][pair]);
          sXcor_fft_backward->Transform();
          sXcor_time = sXcor_fft_backward->GetPointsReal();
          

//if((pol==0)&&(pair==5)) cout<<sXcor_time[0]<<endl;


         double factor = sqrt(Sum_x[i])*sqrt(Sum_x[j])*300.;
         double max = 0;
         int mm;
                       for (Int_t k=0; k<300; k++){
                sXcor_time[k] = sXcor_time[k]/factor;
                if (sXcor_time[k] > max){max = sXcor_time[k]; mm = k;}
                                                    }//k 

//cout<<max<<"  "<<mm<<endl;

if(mm>149) {mm -= 300;}

//cout<<mm<<endl;

tdoa[pol][pair] = mm*0.8;
Xcor[pol][pair] = max;
//cout<<pol<<"    "<<pair<<"   "<<factor<<"    "<<Xcor[pol][pair]<<"  "<<tdoa[pol][pair]<<endl;

     }// pair loop

    }// pol loop

}




	
	
Double_t * Correlator_time = new Double_t [32768];

int N = 32768;
TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");


	
//result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
void get_cross_correlation(int nt_j, double tdoa[2][6], double Xcor[2][6]){	



   //0 for H-pol 1 for V-pol
   for(int pol = 0; pol<2; pol++){
      
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




        double Fine_re[16385] = {0};            
        double Fine_im[16385] = {0};

                  	   for (Int_t k=1; k<151; k++){//0 means DC component
            Fine_re[k] = Correlator_re[pol][pair][k]; 
                   
            Fine_im[k] = Correlator_im[pol][pair][k];
            //cout<<ant_re[i][k]<<endl;
            
            //if((pol==0)&&(pair==5)&&(k==100)) cout<<Fine_re[k]<<endl;
                
                                                    }//k 




          fft_back->SetPointsComplex(Fine_re, Fine_im);
          fft_back->Transform();
          Correlator_time = fft_back->GetPointsReal();


//if((pol==0)&&(pair==5)) cout<<Correlator_time[0]<<endl;


         double factor = sqrt(Sum_x[i])*sqrt(Sum_x[j])*300.;
         double max = 0;
         int mm;

                       for (Int_t k=0; k<32768; k++){//32768
                        //cout<<Correlator_time[k]<<endl;
                correlation_pair[pol][pair][k] = Correlator_time[k]/factor;
                if (correlation_pair[pol][pair][k] > max){max = correlation_pair[pol][pair][k]; mm = k;}

                                                    }//k 

//cout<<max<<"  "<<mm<<endl;

if(mm>16383) {mm -= 32768;}
//if(mm>149) {mm -= 300;}

//cout<<mm<<endl;

tdoa[pol][pair] = mm/32768.*0.8*300.;//32768.
Xcor[pol][pair] = max;
//cout<<pol<<"    "<<pair<<"   "<<factor<<"    "<<Xcor[pol][pair]<<"  "<<tdoa[pol][pair]<<endl;

     }// pair loop

    }// pol loop

}
   


	
	
	


int filter_time(double x[4096], double new_x[4096], int Max_i){

//int Max_i = 0;
//double Max = 0.;

           //for (int i=0;i<1500;i++)
	//{

	  //if(fabs(x[i])>Max) {Max = fabs(x[i]);  Max_i = i;} 	  

	//}
	//cout<<Max_i<<endl;
	
           for (int i=0;i<1500;i++)
	{

   	       new_x[i] = ((i<Max_i-100)||(i>Max_i+100))?0:(0.42-0.5*cos(2*PI*(i-Max_i+100)/200.)+0.08*cos(2*PI*(i-Max_i+100)/200.))*x[i];
   	       
   	       //x[i] = ((i<Max_i-100)||(i>Max_i+100))?0:x[i];
 	                               
	}

           for (int i=1500;i<4096;i++)
	{
   	       new_x[i] = 0;	                               
	}
	
return Max_i;

}



int filter_N = 1500;
Int_t filter_n_size = filter_N+1;
TVirtualFFT *filter_fft_forward = TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
TVirtualFFT *filter_fft_back = TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");
//Double_t * filtered_x = new Double_t [1500];

int filter_fre(double x[4096], double &SNR){

//int N = 1500;
//Int_t n_size = N+1;
   //TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &n_size,"R2C ES K");
      
   //TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
   
   	
double Re_ant[751];
double Im_ant[751];



   filter_fft_forward->SetPoints(x);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(Re_ant,Im_ant);


 for (Int_t k=0; k<751; k++){//0 means DC component // from 185-350 MHz
       if((k<222)||(k>420)) {Re_ant[k] = 0.; Im_ant[k] = 0;}
                             }
  
          
double FreqBin = 5/6.;
   int WeakNotchHalfWidth = 3./FreqBin;  //MHz +-
   int StrongNotchHalfWidth = 6./FreqBin;  //MHz +-


  const int NKnownCW = 8;
  double   CWPeak[ NKnownCW ] = {
    165./FreqBin, 187.5/FreqBin, 244.167/FreqBin, 249.167/FreqBin, 272.5/FreqBin, 
    299.165/FreqBin, 200./FreqBin, 262./FreqBin  //in cal-pulser, not found in forced trig
      //187.5, 272.5 are the strongest
  };
  int CWPeakBin[ NKnownCW ]={0};
  for(int j=0;j<NKnownCW;j++) CWPeakBin[j] = TMath::Nint( CWPeak[j] );
  
  for(int i=0; i<NKnownCW; i++){
//cout<<"(Filter) notch: "<< CWPeak[i]<<endl;
   for(int n= -WeakNotchHalfWidth; n <= WeakNotchHalfWidth; n++){
    //cout<<"Notch: "<< CWPeak[i]+n <<endl;
    Re_ant[ CWPeakBin[i]+n ] *= 0.;
    Im_ant[ CWPeakBin[i]+n ] *= 0.;
   }
  }

//Double_t * filtered_x = new Double_t [1500];

  filter_fft_back->SetPointsComplex(Re_ant,Im_ant);
  filter_fft_back->Transform();
  filtered_x = filter_fft_back->GetPointsReal();



int Max_i = 0;
double Max = 0.;

for(int i=0;i<1500;i++){
	filtered_x[i] /= 1500.;
  						}

for(int i=200;i<1300;i++){
if(fabs(filtered_x[i])>Max) {Max = fabs(filtered_x[i]);  Max_i = i;}
  						}

double power_sum = 0.;  						
for(int i=200;i<1300;i++){
	if((i<Max_i-150)||(i<Max_i+150)){
		power_sum += pow(filtered_x[i],2);
		}
						}
power_sum = sqrt(power_sum/(1500.-300.-1.));
SNR = Max/power_sum;
//delete fft_back;
//fft_back = 0;
//delete fft_forward;					
//fft_forward = 0;
  						
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
         
         
         


   						
     	     						
