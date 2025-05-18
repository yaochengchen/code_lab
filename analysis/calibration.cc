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
using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8
 
 void read_pulser_template();
 void get_response(double delay[4], int pos_theta, int pos_phi);
double cyc_mag[4][2049] = {0};
double cyc_pha[4][2049] = {0};
 void deconvolute_response(int pos_theta, int pos_phi);
 
 
 double get_delay_from_phase(int ant_i, int ant_j, double delay, TH1F *h1);
 
 
 
void get_cross_correlation(int nt_j, double tdoa[2][6], double Xcor[2][6]);
int filter_time(double x[4096], double new_x[4096], int Max_i);
int filter_fre(double x[4096]);
Double_t expectedTimeDiff_angle(double theta, double phi, double *xyz);
Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j);
Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j);

void get_simple_Xcor(double tdoa[2][6], double Xcor[2][6]);
Double_t * filtered_x = new Double_t [1500];
 
 Double_t x[8][4096]={0};//8 antennas
 Double_t new_x[8][4096]={0};//8 antennas
  double ant_re[8][2049] = {0};
  double ant_im[8][2049] = {0};
  
  double tmp_mag[2049] = {0};
  double tmp_pha[2049] = {0};
  double resp_mag[110][5][4][2049] = {0};//  theta:-6 -5.9 .... -5   phi:-24 -24.5 ...-26
  double resp_pha[110][5][4][2049] = {0};

  double resp_re[110][5][4][2049] = {0};//  theta:-6 -5.9 .... -5   phi:-24 -24.5 ...-26
  double resp_im[110][5][4][2049] = {0}; 

 double correlation_pair[2][6][32768] = {0};
 
 
double Xcor_tdoa[6][10000]={0};
double Xcor_max[10000][6]={0};
double Sum_x[8] = {0};
double Correlator_re[2][6][151] = {0};            
double Correlator_im[2][6][151] = {0};
int peak_position[8] = {0};

double Expect_Xcor_tdoa[6][20000]={0};

 
//int main(int argc, char *argv[]) {
 void calibration(){
 



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
 
  //delay_ant[2] += ((-500 + 31*900./100.) + (-50 + 54*100./100.));
  //delay_ant[4] += ((-500 + 14*900./100.) + (-50 + 62*100./100.));
  //delay_ant[6] += ((-500 + 37*900./100.) + (-50 + 70*100./100.));
 
  
TFile *cyc_ntuple_file = new TFile("cyc_ntuple_dense.root","RECREATE");

TNtuple *cyc_ntuple = new TNtuple("cyc_ntuple","cyc_ntuple","time:theta:phi:R:delta_theta:delta_phi:X_cor",320000);
 
 
TFile *first_ntuple_file = new TFile("/home/cyc/20210121/first_flight/ntuple.root");
//TFile *first_ntuple_file = new TFile("../second_third_flight/ntuple.root");
TNtuple * first_ntuple = (TNtuple*) first_ntuple_file->Get("ntuple");

TFile *third_ntuple_file = new TFile("/home/cyc/20210121/second_third_flight/ntuple.root");
TNtuple * third_ntuple = (TNtuple*) third_ntuple_file->Get("ntuple");

TH2F *h2 = new TH2F("h2","h2", 25, -26.2, -23.7, 110, -6, 5);
TCanvas *c4 = new TCanvas("c4","Canvas Example",1500,1000);

TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
c1->Divide(3,2);

TCanvas *c2 = new TCanvas("c2","Canvas Example",1500,1000);
c2->Divide(2,1);

//TCanvas *c3 = new TCanvas("c3","Canvas Example",1500,1000);

//c3->Divide(1,3);

TRandom3 *r = new TRandom3();

read_pulser_template();

float nt_run;
float nt_event;
float nt_time;
float drn_azimuth;
float drn_elevation;
// ntuple->SetBranchAddress("run",&nt_run);
// ntuple->SetBranchAddress("event",&nt_event);
// ntuple->SetBranchAddress("time",&nt_time);
// ntuple->SetBranchAddress("drn_azimuth",&drn_azimuth);
// ntuple->SetBranchAddress("drn_elevation",&drn_elevation);

float rel_x;
float rel_y;
float rel_z;
// ntuple->SetBranchAddress("rel_x",&rel_x);
// ntuple->SetBranchAddress("rel_y",&rel_y);
// ntuple->SetBranchAddress("rel_z",&rel_z);



first_ntuple->SetBranchAddress("run",&nt_run);
first_ntuple->SetBranchAddress("event",&nt_event);
first_ntuple->SetBranchAddress("time",&nt_time);
first_ntuple->SetBranchAddress("drn_azimuth",&drn_azimuth);
first_ntuple->SetBranchAddress("drn_elevation",&drn_elevation);
first_ntuple->SetBranchAddress("rel_x",&rel_x);
first_ntuple->SetBranchAddress("rel_y",&rel_y);
first_ntuple->SetBranchAddress("rel_z",&rel_z);


third_ntuple->SetBranchAddress("run",&nt_run);
third_ntuple->SetBranchAddress("event",&nt_event);
third_ntuple->SetBranchAddress("time",&nt_time);
third_ntuple->SetBranchAddress("drn_azimuth",&drn_azimuth);
third_ntuple->SetBranchAddress("drn_elevation",&drn_elevation);
third_ntuple->SetBranchAddress("rel_x",&rel_x);
third_ntuple->SetBranchAddress("rel_y",&rel_y);
third_ntuple->SetBranchAddress("rel_z",&rel_z);



const Int_t first_nt_nentries = (Int_t)first_ntuple->GetEntries() - 1;
const Int_t third_nt_nentries = (Int_t)third_ntuple->GetEntries() - 1;





double event_time[20000] = {0};
double event_theta[20000] = {0};
double dec_event_theta[20000] = {0};
double event_phi[20000] = {0};
double drone_azi[20000] = {0};
double drone_ele[20000] = {0};

double relative_coord[20000][3] = {0};

int nt_j = 0;
//cout<<nt_nentries<<endl;
//ntuple->GetEntry(nt_j);
//cout<<nt_run<<"  "<<nt_event<<endl;
//
   
// int number1 = 16946;
// int number2 = 16948;
 
 
 
 double Sum_template = 0;
double Re_template[2049] = {0};
double Im_template[2049] = {0};

 
TH1F *h1[6];
for(int i=0; i<6; i++){
 h1[i] = new TH1F("h1","h1", 100, -2, 2);
}
 //double rotate_degree = 1;

 
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
  



double r_deg[40] = {0};
double r_mean[40] = {0};
double r_sigma[40] = {0};

bool first_time = true;
bool second_time = true;

int count_first_cut = 0;

for(int rot_deg = 0; rot_deg<1; rot_deg++){

 int number1 = 16955;
 int number2 = 16957;
int hahaha = 16956;

double rotate_degree = 0;
nt_j=0;
int entry = 0;

  for(int number = number1; number<number2; number++){
  //continue;

  
    fname_evt =  Form("/home/cyc/20210121/run%08d.root",number);
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
 
 
 
for (; entry<third_nt_nentries; entry++){

      third_ntuple->GetEntry(entry);
      drn_azimuth = 90 - drn_azimuth;
      
      
      if((drn_azimuth<-50) || (drn_azimuth>50) || (drn_elevation<-6) || (drn_elevation>50)) continue;

      //if((nt_time<1362.6)&&(nt_time>1362.4)) {cout<<nt_event<<"   "<<number<<"   "<<nt_j<<endl;}
      
      //if((nt_time<1550)||(nt_time>1551)) {continue;}
      //if((nt_time<1650)||(nt_time>1850)) {continue;}
      
//cout<<"print here 1"<<endl;
        					 
      
      if(((int) nt_run)==hahaha){hahaha++; break;}
      //continue;
//cout<<nt_run<<"  "<<nt_event<<"   "<<number<<"   "<<nt_j<<endl;   


   
      int j = (int) nt_event;
      
//for (int j=0;j<event_number;j++){ 
  
      Tree_Muon->GetEntry(j);



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

     event_time[nt_j] = nt_time;
     

     

     drone_ele[nt_j] = drn_elevation;
     drone_azi[nt_j] = drn_azimuth - rotate_degree;
     
     relative_coord[nt_j][0] = rel_x;
     relative_coord[nt_j][1] = rel_y;
     relative_coord[nt_j][2] = rel_z;
     
     double R = sqrt(rel_x*rel_x + rel_y*rel_y + rel_z*rel_z);
     //cout<<"input:   "<<R<<endl;
     
double tdoa[2][6] = {0};
double Xcor[2][6] = {0};

get_simple_Xcor(tdoa, Xcor);

//cout<<tdoa[0][1]<<"   "<<Xcor[0][1]<<endl;

double H_Xcor = 0;
double V_Xcor = 0;
	
for(int pair = 0; pair <6; pair++){
    H_Xcor += Xcor[0][pair]/6.;
    V_Xcor += Xcor[1][pair]/6.;
    }
    
    int flag = 0;
    if((H_Xcor>0.7)&&(V_Xcor<0.7)) {flag = 1;}
    //if((H_Xcor<0.7)&&(V_Xcor>0.7)) {flag = 2;}
    if((H_Xcor>0.7)&&(V_Xcor>0.7)) {flag = 3;}
    
    if(flag>0) count_first_cut++;
    
	
//continue;

get_cross_correlation(nt_j, tdoa, Xcor);

//cout<<tdoa[0][1]<<"   "<<Xcor[0][1]<<endl<<endl;

//get_cross_correlation(nt_j);

for(int pair = 0; pair <6; pair++){

int pol = 0;
int i = 0;
int j = 2;
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
    
   double delta_t = expectedTimeDiff_angle_with_R(drn_elevation, drn_azimuth, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
   
   delta_t = delta_t - (peak_position[j] - peak_position[i])*0.8;
   
Expect_Xcor_tdoa[pair][nt_j] = tdoa[0][pair] - delta_t;

//cout<<pair<<"  "<<tdoa[0][pair]<<"    "<<delta_t<<endl;

}
//cout<<endl;
/*
double test_Xcor = 0;
for(int pair = 0; pair <6; pair++){

int pol = 0;
int i = 0;
int j = 2;
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

          double delta_t = expectedTimeDiff_angle_with_R(11.3, -25.8, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          
          delta_t = delta_t - (peak_position[j] - peak_position[i])*0.8;
          delta_t = delta_t/(0.8*300./32768.);
          if(delta_t<0) delta_t += 32768;
          if(delta_t<0||delta_t>32768) continue; //cout<<"delta_t error!  "<<delta_t<<endl;

          //cout<<delta_t<<endl;
          //cout<<correlation_pair[pol][pair][(int) delta_t]<<endl;
          test_Xcor += correlation_pair[pol][pair][(int) delta_t];

}

cout<<test_Xcor<<endl;
*/



/*
   //0 for H-pol 1 for V-pol
   for(int pol = 0; pol<1; pol++){
      
   	int i = 0;
    int j = 2;
    
double inter_map[56][100] = {0};//theta: -6~50; phi: -50~50
double inter_max = 0;
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
          double phi = -50 + i_phi;
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
          double phi = -1 + i_phi/10. + max_phi[w]-50;
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
	
     inter_max = fine_max;
     
     event_theta[nt_j] = -2 + fine_max_theta/10. + max_theta[max_w]-6;
     event_phi[nt_j] = -1 + fine_max_phi/10. + max_phi[max_w]-50;
	
	//cout<<max_theta[0]-6<<"   "<<max_phi[0]-50<<endl;
	cout<<pol<<"  "<<drn_elevation<<"   "<<drn_azimuth<<"   "<<event_theta[nt_j]<<"   "<<event_phi[nt_j]<<"   "<<inter_max<<"   "<<max_w<<endl;
	
	}//pol loop

*/     
	
nt_j++;

}// number loop

}// run loop
     

cout<<count_first_cut<<endl;










     TGraph * cyc[2];
     TGraph * wyn[2];
     TGraph * hxj[2];
     
     for(int pol = 0; pol<2; pol++){
     
     c2->cd(pol+1);
     
     if(pol==0) {
     cyc[pol]= new TGraph(nt_j, event_time, event_theta);
     wyn[pol] = new TGraph(nt_j, event_time, drone_ele);
     hxj[pol]= new TGraph(nt_j, event_time, dec_event_theta);
     cyc[pol]->Draw("AP");
     cyc[pol]->GetYaxis()->SetTitle("elevation");
     hxj[pol]->Draw("P same");
     hxj[pol]->SetMarkerStyle(kCircle);
     hxj[pol]->SetMarkerColor(kBlue);
     }
     else {
     cyc[pol]= new TGraph(nt_j, event_time, event_phi);
     wyn[pol] = new TGraph(nt_j, event_time, drone_azi);
     cyc[pol]->Draw("AP");
     cyc[pol]->GetYaxis()->SetTitle("azimuth");
     }

     cyc[pol]->SetMarkerStyle(kCircle);
     cyc[pol]->SetMarkerColor(kRed);
     cyc[pol]->GetXaxis()->SetTitle("time:s");
     

     

     wyn[pol]->Draw("PL same");
     
     }










c4->cd();
TGraph * abc[6];
     for(int pair = 0; pair<6; pair++){
     abc[pair] = new TGraph(nt_j, event_time, Expect_Xcor_tdoa[pair]);//event_time drone_azi

     if(pair==0){
     abc[pair]->Draw("AP");
     abc[pair]->GetXaxis()->SetTitle("time: s");//time:s
     abc[pair]->GetYaxis()->SetTitle("ns");
     abc[pair]->GetYaxis()->SetRangeUser(-2,2);
     abc[pair]->GetXaxis()->SetRangeUser(1000,1550);
     }
     else{abc[pair]->Draw("P");}
     
     
     abc[pair]->SetMarkerStyle(kCircle);
     abc[pair]->SetMarkerColor(pair+1);
     }






}//rot_deg loop



     







   cyc_ntuple_file->cd();
   cyc_ntuple->Write();
   cyc_ntuple_file->Close();
   delete cyc_ntuple_file;
   
    
              }// end of main function
              
              


Double_t * sXcor_time = new Double_t [300];

   int sXcor_N = 300;
   Int_t sXcor_size = sXcor_N+1;		
   TVirtualFFT *sXcor_fft_forward = TVirtualFFT::FFT(1, &sXcor_size,"R2C ES K");
   
   TVirtualFFT *sXcor_fft_backward = TVirtualFFT::FFT(1, &sXcor_N, "C2RBACKWARD M K");
   
void get_simple_Xcor(double tdoa[2][6], double Xcor[2][6]){

   double re[8][151] = {0};
   double im[8][151] = {0};

//double Sum_x[8] = {0};
  
	for(int ant = 0; ant < 8; ant ++){
		Sum_x[ant] = 0;
	
    	int max_i = filter_fre(x[ant]);
    	peak_position[ant] = max_i;
    	filter_time(filtered_x, new_x[ant], max_i);  
    	

   sXcor_fft_forward->SetPoints(&new_x[ant][max_i-150]);
   sXcor_fft_forward->Transform();
   sXcor_fft_forward->GetPointsComplex(re[ant],im[ant]);

    	  for(int i = 0; i <300; i++){
          Sum_x[ant] += pow(new_x[ant][max_i-150+i], 2);
          
          
          }
   
   }         
              


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

int filter_fre(double x[4096]){

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
         
         
         

void read_pulser_template(){

 int N = 4096;
 Int_t n_size = N+1;
 double x[4096] = {0};
 double new_x[4096] = {0};
 double t[4096] = {0};
 double new_t[4096] = {0};

 Short_t T1TH[1500];
 double Re_template[2049] = {0};
 double Im_template[2049] = {0};
   TVirtualFFT *fft_forward_tmp = TVirtualFFT::FFT(1, &n_size,"R2C ES K");



   TFile * file_template = new TFile("/home/cyc/20210121/run-rbm00019207.root");
   TTree * Tree_template = (TTree*) file_template->Get("t");

		  Tree_template->SetBranchAddress("T1TH",T1TH);

		  


      Tree_template->GetEntry(1);



           for (int i=0;i<1500;i++) {x[i]= (float)T1TH[i]*500/32512;}
           for (int i=0;i<4096;i++) {t[i] = i*0.8;}
           
	int max_i = filter_fre(x);
	filter_time(x, new_x, max_i);  


   fft_forward_tmp->SetPoints(new_x);
   fft_forward_tmp->Transform();
   fft_forward_tmp->GetPointsComplex(Re_template,Im_template);
   
   for(int i = 1; i<2049; i++){
   		double magnitude = sqrt(Re_template[i] * Re_template[i] + Im_template[i] *  Im_template[i]);
   		double phase = TMath::ATan(Im_template[i]/(Re_template[i]+1.0e-12));
     	if(Re_template[i]<0){phase += TMath::Pi();}
     	if(phase>TMath::Pi()) {phase -= 2*TMath::Pi();}
     	tmp_mag[i] = magnitude;
     	tmp_pha[i] = phase;// - 2*PI*(10e-9)*(i*(1250e6)/4096.);
     	//Re_template[i] = magnitude*TMath::Cos(tmp_pha[i]);
     	//Im_template[i] = magnitude*TMath::Sin(tmp_pha[i]);
     	     						}
     	     						
   /*  	     						
   TVirtualFFT *fft_back_tmp = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");

    double * time_signal = new double [4096];

          fft_back_tmp->SetPointsComplex(Re_template, Im_template);
          fft_back_tmp->Transform();
          time_signal = fft_back_tmp->GetPointsReal(); 
          
   for (int i=0;i<4096;i++) {time_signal[i] = time_signal[i]/1500.; new_t[i] = 0.8*i*1500./4096.;}      
           
   TGraph * draw_template = new TGraph(4096,t,x);
   draw_template->Draw();

   TGraph * draw_delay = new TGraph(4096,t,time_signal);
   draw_delay->Draw();
   draw_delay->SetLineColor(2);
   delete fft_back_tmp;
   */   
   delete fft_forward_tmp;
   fft_forward_tmp = 0;
   
   delete Tree_template;
   file_template->Close();
   delete file_template;
}         

//TH1F *h1 = new TH1F("h1","h1", 542, 606*625/2048., 1148*625/2048.);
//TCanvas *c5 = new TCanvas("c5","Canvas Example",1500,1000);
        
void get_response(double delay[4], int pos_theta, int pos_phi){

   int max_i = filter_fre(x[3*2]);
   filter_time(x[3*2], new_x[3*2], max_i);  

 double Re_template[2049] = {0};
 double Im_template[2049] = {0};
 
   int N = 4096;
   Int_t n_size = N+1;
   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &n_size,"R2C ES K");
   fft_forward->SetPoints(new_x[3*2]);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(Re_template,Im_template);

   for(int i = 1; i<2049; i++){
   		double magnitude = sqrt(Re_template[i] * Re_template[i] + Im_template[i] *  Im_template[i]);
   		double phase = TMath::ATan(Im_template[i]/(Re_template[i]+1.0e-12));
     	if(Re_template[i]<0){phase += TMath::Pi();}
     	if(phase>TMath::Pi()) {phase -= 2*TMath::Pi();}
     	tmp_mag[i] = magnitude;
     	tmp_pha[i] = phase;// - 2*PI*(10e-9)*(i*(1250e6)/4096.);
     	//Re_template[i] = magnitude*TMath::Cos(tmp_pha[i]);
     	//Im_template[i] = magnitude*TMath::Sin(tmp_pha[i]);
     	     						}



	for(int ant = 0; ant < 4; ant ++){
	
 		int max_i = filter_fre(x[ant*2]);
	    filter_time(x[ant*2], new_x[ant*2], max_i);  

   double re[2049] = {0};
   double im[2049] = {0};
   int N = 4096;
   Int_t n_size = N+1;
   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &n_size,"R2C ES K");
   fft_forward->SetPoints(new_x[ant*2]);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(re,im);
   
   for(int i = 606; i<1148; i++){
   		double magnitude = sqrt(re[i] * re[i] + im[i] * im[i]);
   		double phase = TMath::ATan(im[i]/(re[i]+1.0e-12));
   		
   		//h1->SetBinContent(i-606, log10(magnitude));
   		
     	if(re[i]<0){phase += TMath::Pi();}
     	if(phase>TMath::Pi()) {phase -= 2*TMath::Pi();}
     	double tmp_delay_phase = tmp_pha[i] - 2*PI*(delay[ant]*1e-9)*(i*(1250e6)/4096.);
     	double response_mag = tmp_mag[i]/magnitude;
     	double response_pha = tmp_delay_phase - phase;
     	
     	//resp_mag[pos_theta][pos_phi][ant][i] = response_mag;
     	//resp_pha[pos_theta][pos_phi][ant][i] = response_pha;
     	
     	resp_re[pos_theta][pos_phi][ant][i] += response_mag*TMath::Cos(response_pha);
     	resp_im[pos_theta][pos_phi][ant][i] += response_mag*TMath::Sin(response_pha);
     	//resp_re[pos_theta][pos_phi][ant][i] = response_mag*TMath::Cos(response_pha);
     	//resp_im[pos_theta][pos_phi][ant][i] = response_mag*TMath::Sin(response_pha);
     	
     	resp_mag[pos_theta][pos_phi][ant][i] = sqrt(pow(resp_re[pos_theta][pos_phi][ant][i],2) + pow(resp_im[pos_theta][pos_phi][ant][i],2));
     	resp_pha[pos_theta][pos_phi][ant][i] = TMath::ATan(resp_im[pos_theta][pos_phi][ant][i]/(resp_re[pos_theta][pos_phi][ant][i]+1.0e-12));
     	if(resp_re[pos_theta][pos_phi][ant][i]<0){resp_pha[pos_theta][pos_phi][ant][i] += TMath::Pi();}
     	if(resp_pha[pos_theta][pos_phi][ant][i]>TMath::Pi()) {resp_pha[pos_theta][pos_phi][ant][i] -= 2*TMath::Pi();}

     	     				 }
/*
  c5->cd();
  h1->Draw();
  //h1->GetXaxis()->SetRangeUser(3, 6);
  h1->SetStats(kFALSE);
  h1->GetXaxis()->SetTitle("MHz");
  h1->GetYaxis()->SetTitle("dB");
  char * name = Form("tmpt_theta_%d__phi_%d__ant_%d__.png", pos_theta, pos_phi, ant);
  c5->SaveAs(name);
*/
     	     				 			}
}
     	     						

   int decon_N = 4096;
   Int_t decon_n_size = decon_N+1;		
   TVirtualFFT *decon_fft_forward = TVirtualFFT::FFT(1, &decon_n_size,"R2C ES K");
       	     						
void deconvolute_response(int pos_theta, int pos_phi){

   //int N = 4096;
   //Int_t n_size = N+1;		
   //TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &n_size,"R2C ES K");
//double cyc_mag[4][2049] = {0};
//double cyc_pha[4][2049] = {0};   
	for(int ant = 0; ant < 4; ant ++){
	
    
 		int max_i = filter_fre(x[ant*2]);
	    filter_time(x[ant*2], new_x[ant*2], max_i);  

   double re[2049] = {0};
   double im[2049] = {0};


   decon_fft_forward->SetPoints(new_x[ant*2]);
   decon_fft_forward->Transform();
   decon_fft_forward->GetPointsComplex(re,im);
   
   for(int i = 606; i<1148; i++){
        //if((i<606)||(i>1147)) {ant_re[ant*2][i] = 0; ant_im[ant*2][i] = 0; continue;}
   		double magnitude = sqrt(re[i] * re[i] + im[i] * im[i]);
   		
   		//h1->SetBinContent(i-606, log10(magnitude));
   		
   		double phase = TMath::ATan(im[i]/(re[i]+1.0e-12));
     	if(re[i]<0){phase += TMath::Pi();}
     	if(phase>TMath::Pi()) {phase -= 2*TMath::Pi();}
     	//double new_mag = magnitude * resp_mag[pos_theta][pos_phi][ant][i];
     	double new_mag = magnitude;
     	//double new_mag = magnitude / resp_mag[pos_theta][pos_phi][ant][i];
     	//double new_mag = tmp_mag[i];
     	//double new_mag = magnitude * (resp_mag[pos_theta][pos_phi][ant][i]<3?resp_mag[pos_theta][pos_phi][ant][i]:3);
     	//double new_mag = (magnitude>pow(10,1.7)?magnitude:0) * (resp_mag[pos_theta][pos_phi][ant][i]<3?resp_mag[pos_theta][pos_phi][ant][i]:3);
     	double new_pha = phase + resp_pha[pos_theta][pos_phi][ant][i];
     	new_pha = fmod(new_pha, 2*PI);
   		if(new_pha>PI) new_pha = new_pha - 2*PI;
   		if(new_pha<(-1*PI)) new_pha = new_pha + 2*PI;
   		
     	cyc_mag[ant][i] = new_mag;
     	cyc_pha[ant][i] = new_pha;
     	
     	re[i] = new_mag*TMath::Cos(new_pha);
     	im[i] = new_mag*TMath::Sin(new_pha);
     	ant_re[ant*2][i] = re[i];
     	ant_im[ant*2][i] = im[i];
     	     				 }
/*
  c5->cd();
  h1->Draw();
  //h1->GetXaxis()->SetRangeUser(3, 6);
  h1->SetStats(kFALSE);
  h1->GetXaxis()->SetTitle("MHz");
  h1->GetYaxis()->SetTitle("dB");
  char * name = Form("event_theta_%d__phi_%d__ant_%d__.png", pos_theta, pos_phi, ant);
  c5->SaveAs(name);
*/  
     	     				 			}
	//delete fft_forward;
	//fft_forward = 0;
/*	
	for(int ant = 0; ant < 4; ant ++){
   for(int i = 606; i<1148; i++){
        if((i<606)||(i>1147)) {ant_re[ant*2][i] = 0; ant_im[ant*2][i] = 0; continue;}
        double new_mag = pow(cyc_mag[0][i] * cyc_mag[1][i] * cyc_mag[2][i] * cyc_mag[3][i], 0.25);

     	ant_re[ant*2][i] = new_mag*TMath::Cos(cyc_pha[ant][i]);
     	ant_im[ant*2][i] = new_mag*TMath::Sin(cyc_pha[ant][i]);
     	}
     }
*/     
}

double delta_t_phase[2048] = {0};

double get_delay_from_phase(int ant_i, int ant_j, double delay, TH1F * h1){

   	double mean = 0;
   	double m_count = 0;
   	for(int i = 606; i<1148; i++){
   		double extress = 2*PI*(delay*1.0e-9)*(i*(1250e6)/4096.);
   		extress = fmod(extress, 2*PI); 
   		//if(extress>PI) {extress = extress - 2 * PI;}
   		//if(extress<(-1*PI)) {extress = extress + 2*PI;}
   		double delta_phase = cyc_pha[ant_j][i] - cyc_pha[ant_i][i];
   		//if(delta_phase>PI) delta_phase = delta_phase - 2*PI;
   		//if(delta_phase<(-1*PI)) delta_phase = delta_phase + 2*PI;
   		delta_phase  = delta_phase + extress;
   		delta_phase = fmod(delta_phase, 2*PI);
   		if(delta_phase>PI) delta_phase = delta_phase - 2*PI;
   		if(delta_phase<(-1*PI)) delta_phase = delta_phase + 2*PI;
   		delta_t_phase[i] = (delta_phase)/(2*PI*(i*(1250e6)/4096.));
   		delta_t_phase[i] *= 1.0e9;//s to ns
   		if(fabs(delta_t_phase[i])<2){
   		double factor = 10*log10(sqrt(cyc_mag[ant_i][i]*cyc_mag[ant_j][i]));
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
   	   			double factor = 10*log10(sqrt(cyc_mag[ant_i][i]*cyc_mag[ant_j][i]));
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
   	   			    ant_re[ant*2][i] = 0.;
   	   			    ant_im[ant*2][i] = 0.;
   	   			}
   	   		}
   		
   		}		
		

	return fabs(fine_mean);
}
   						
     	     						
