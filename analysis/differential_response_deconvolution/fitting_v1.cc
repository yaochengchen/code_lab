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
#define Earth_Radius 6371000.
#define Drone_Altitude 704.7683
//From google earth, T4 location is about 686 m

Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j);
Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j);



void reflection_find_source(double direct_angle, double reflect_angle, double T4_Altitude, double &distance_X_max, double &delay, double &distance_t4_reflection_point){
	
	double antenna_altitude = T4_Altitude;
	//double direct_angle = 2.35;//degrees
	//double reflect_angle = -3.1;//degrees
	double c = Earth_Radius; // Earth radius
	double a = c + antenna_altitude;// T4
	double C = (90 +reflect_angle)*TMath::DegToRad(); // CR coming angle
	
	double b = a*TMath::Cos(C) - sqrt(pow(a*TMath::Cos(C),2) + (c*c-a*a));
	
	distance_t4_reflection_point = b;
	//cout<<b<<"  "<<a*TMath::Cos(C)<<"   "<<(pow(a*TMath::Cos(C),2) + (c*c-a*a))<<"   "<<90 +reflect_angle<<endl;
	
	double A = TMath::ACos((b*b+c*c-a*a)/(2.*b*c))*TMath::RadToDeg();
	//cout<<A-90<<endl;
	TVector2 T4_position(0., antenna_altitude);
	TVector2 reflection_point(b*TMath::Cos(reflect_angle*TMath::DegToRad()), b*TMath::Sin(reflect_angle*TMath::DegToRad())+T4_position.Y());
	//cout<<reflection_point.X()<<"  "<<reflection_point.Y()<<endl;
	
	
	double differ_angle = fabs(reflect_angle) - (A-90);
	double reflection_point_source_angle = (A-90) - differ_angle;
	//cout<<"reflection_point_source_angle:   "<<reflection_point_source_angle<<endl;
	
	double k_reflect = TMath::Tan(reflection_point_source_angle*TMath::DegToRad());
	double b_reflect = reflection_point.Y() - k_reflect*reflection_point.X();//y-kx
	
	double k_direct = TMath::Tan(direct_angle*TMath::DegToRad());
	double b_direct = T4_position.Y() - k_direct*T4_position.X();//y-kx
	
	double source_x = (b_direct-b_reflect)/(k_reflect-k_direct);
	double source_y = k_direct*source_x + b_direct;
	TVector2 source_position(source_x, source_y);
	//cout<<source_x<<"   "<<source_y<<endl;
	
	double S_d = (source_position-T4_position).Mod();
	double S_r_1 = (reflection_point-T4_position).Mod();
	double S_r_2 = (source_position-reflection_point).Mod();
	//cout<<(S_d-S_r_1-S_r_2)/Speed_Of_Light*1.0e9<<endl;
	//cout<<S_d<<endl;
	
	distance_X_max = S_d;
	delay = (S_d-S_r_1-S_r_2)/Speed_Of_Light*1.0e9;
	
	if(source_x<0){delay=-100000000;}
	
	//cout<<distance_X_max<<"   "<<delay<<endl;


}






void source_find_reflection(double direct_angle, double Distance_Xmax, double T4_Altitude, double &delay, double &reflected_angle, double &distance_t4_reflection_point){

	//double direct_angle = 2.35;//degrees
	//double Distance_Xmax = 100e3;
	
	
	double trial_delay;
	double trial_distance;
	double trial_reflect_angle;
	
	double distance_error;
	
	double trial_reflect_angle_s = (-1)*direct_angle - 10.;
	double trial_reflect_angle_b = (-1)*TMath::ACos(Earth_Radius/(Earth_Radius+T4_Altitude))*TMath::RadToDeg()-0.000001;
	double epsilon = Distance_Xmax/10000.;
	//cout<<"trial_reflect_angle_b : "<<trial_reflect_angle_b<<endl;
		
	int iteration = 100;
	while(iteration){
		
		//cout<<90-TMath::ACos(Earth_Radius/(Earth_Radius+T4_Altitude))*TMath::RadToDeg()<<endl;
		
		trial_reflect_angle = (trial_reflect_angle_s + trial_reflect_angle_b)/2.;
		if(trial_reflect_angle>0){cout<<trial_reflect_angle<<endl;}
		
		reflection_find_source(direct_angle, trial_reflect_angle, T4_Altitude, trial_distance, trial_delay, distance_t4_reflection_point);
		distance_error = trial_distance - Distance_Xmax;
		//cout<<direct_angle<<"   "<<Distance_Xmax<<"   "<<trial_reflect_angle<<"    "<<distance_error<<"    "<<trial_delay<<endl;
		if(fabs(distance_error)<epsilon){break;}
		if(trial_delay == -100000000){trial_reflect_angle_b = trial_reflect_angle;}
		else if(distance_error<0){trial_reflect_angle_s = trial_reflect_angle;}
		else if(distance_error>0){trial_reflect_angle_b = trial_reflect_angle;}
		
		//cout<<direct_angle<<"   "<<Distance_Xmax<<"   "<<trial_reflect_angle<<"    "<<distance_error<<"    "<<trial_delay<<endl;
		iteration --;
		
	}
	//cout<<distance_error<<"    "<<trial_delay<<endl;
	
	delay = trial_delay;
	reflected_angle = trial_reflect_angle;
	if(!iteration) delay = -100000000;

}






//base center   2021 photogrametry  Tower-4 is back side
double coor_ant[8][3] = {{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949},{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949}};

  double delay_ant[8] = {0,53.2,84.6,-1523.1,-92.4,-30.3,59.5,-743.1};


void initial(){




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


 
}
 
//int main(int argc, char *argv[]) {
 void fitting(){
 	initial();
 	double measured_delays[4] = {-225.952, -225.184, -224.640, -224.368};
 	double measured_antenna_delays[6] = {-3.1521681, -6.0424853, -47.342842, -2.8903172, -44.190674, -41.300357};
 	double measured_antenna_second_delays[6] = {-3.6396090, -7.0978297, -48.734993, -3.4582206, -45.095384, -41.637163};
 	
 	double min_chi = 1.0e12;
 	double min_R, min_theta, min_phi, min_z;
 	for(double theta=2.0; theta<=3; theta+=0.1){
 		for(double phi=-23; phi<=-22.2; phi+=0.1){
 			for(double R=10000; R<=500000; R+=10){
 				for(double z=665; z<=730; z+=1){
 					double delay = 0.;
 					double reflected_angle[4] = {0.};
 					double distance_t4_reflection_point[4] = {0.};
 					double chi_square = 0.;
 					for(int ant=0; ant<4; ant++){
 						double antenna_altitude = z + coor_ant[ant][2];
 						source_find_reflection(theta, R, antenna_altitude, delay, reflected_angle[ant], distance_t4_reflection_point[ant]);		
 						chi_square += (pow((delay-measured_delays[ant])/0.04,2));
 						//cout<<delay<<endl;
 					}
 					
 					for(int pair = 0; pair <6; pair++){

     					int i = 0;
     					int j = 1;

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
							
							// direct signal
          					double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          					chi_square += (pow((delta_t-measured_antenna_delays[pair])/0.04,2));
          					
          					
          					
          					// reflected signal
          					double reflected_theta_angle = reflected_angle[1];
          					double distance_reflection = distance_t4_reflection_point[1];
          					
          					double delta_t_second = expectedTimeDiff_angle_with_R(reflected_theta_angle, phi, distance_reflection, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          					chi_square += (pow((delta_t_second-measured_antenna_second_delays[pair])/0.08,2));
          					
          			}// pair loop
          			
          			//cout<<chi_square<<endl;
          			
 					if(chi_square < min_chi){
 						min_chi = chi_square;
 						min_R = R;
 						min_theta = theta;
 						min_phi = phi;
 						min_z = z;
 						cout<<min_chi<<"  "<<min_R<<"  "<<min_theta<<"  "<<min_phi<<"  "<<min_z<<endl;
 						cout<<"   "<<reflected_angle[0]<<"   "<<reflected_angle[1]<<"   "<<reflected_angle[2]<<"   "<<reflected_angle[3]<<endl;
 						cout<<"   "<<distance_t4_reflection_point[0]<<"   "<<distance_t4_reflection_point[1]<<"   "<<distance_t4_reflection_point[2]<<"   "<<distance_t4_reflection_point[3]<<endl;
 						
 					}
 				}
 			}
 		}
 	}
 	
 	



    
              }// end of main function
              

     

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
         
         
         


     	     						

