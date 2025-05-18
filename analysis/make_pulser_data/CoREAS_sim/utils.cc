// add (-1) to E_theta_station to match definition of HFSS, which is downward


#include "TComplex.h"
#include "TMath.h"
#include <dirent.h>
#include <vector>
#include "constant_value.h"
#include <TVector3.h>
using namespace std;
  
  
template <typename T>
T** new_Array2D(int row, int col);

template <typename T> 
void delete_Array2D(T **arr, int row, int col);

template <typename T>
T*** new_Array3D(int height, int row, int col);

template <typename T> 
void delete_Array3D(T ***p, int height, int row, int col);
  
  
#define Speed_Of_Light 2.99792458e8  
Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j);  
  
  
//base center   2021 photogrametry after Tower-4 back side correction
double coor_ant[8][3] = {{-4.863, -9.097, 1.529}, {-4.863, -9.097, 1.529}, {0, 0, 0}, {0, 0, 0}, {4.449, 8.261, -1.78}, {4.449, 8.261, -1.78}, {16.9349, 4.63045, -2.91101}, {16.9349, 4.63045, -2.91101}};  
  
  
  
  
  
  
  
  
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
      //cout<<entry<<endl;
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


	for(int i_theta=0; i_theta<37; i_theta++){
	  for(int i=0; i<2049; i++){
  	      double fre = FreqBin*i;
  	      
  	  	  ANT_V_theta_re[i_theta][72][i] = ANT_V_theta_re[i_theta][0][i];
  	  	  ANT_V_theta_im[i_theta][72][i] = ANT_V_theta_im[i_theta][0][i];
  	  	  ANT_V_phi_re[i_theta][72][i] = ANT_V_phi_re[i_theta][0][i];
  	  	  ANT_V_phi_im[i_theta][72][i] = ANT_V_phi_im[i_theta][0][i];
  	  	   	  
  	  }
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
      //cout<<entry<<endl;
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
    
	for(int i_theta=0; i_theta<37; i_theta++){
	  for(int i=0; i<2049; i++){
  	      double fre = FreqBin*i;
  	      
  	  	  ANT_H_theta_re[i_theta][72][i] = ANT_H_theta_re[i_theta][0][i];
  	  	  ANT_H_theta_im[i_theta][72][i] = ANT_H_theta_im[i_theta][0][i];
  	  	  ANT_H_phi_re[i_theta][72][i] = ANT_H_phi_re[i_theta][0][i];
  	  	  ANT_H_phi_im[i_theta][72][i] = ANT_H_phi_im[i_theta][0][i];
  	  	   	  
  	  }
  	}
	
}




void get_candidate_file(string path, std::vector<string> &filenames){
	DIR *pDir;
	struct dirent* ptr;
	if(!(pDir = opendir(path.c_str()))){
		cout<<"Folder doesn't Exist!"<<endl;
		return;
	}
	//cout<<"open folder"<<endl;
	//cout<<readdir(pDir)<<endl;
	while((ptr = readdir(pDir))!=0){
		if(strstr(ptr->d_name, ".dat")&&!strstr(ptr->d_name, ".dat~")){
			filenames.push_back(path+"/"+ptr->d_name);
			//cout<<ptr->d_name<<endl;
		}
	}
	closedir(pDir);
}


   	int FFT_N_coreas = 4096*8;
   	Int_t iFFT_N_coreas = FFT_N_coreas;		
   	TVirtualFFT *FFT_forward_coreas_theta = TVirtualFFT::FFT(1, &FFT_N_coreas,"R2C ES K");
   	TVirtualFFT *FFT_forward_coreas_phi = TVirtualFFT::FFT(1, &FFT_N_coreas,"R2C ES K");
    
    TVirtualFFT *FFT_backward_coreas = TVirtualFFT::FFT(1, &iFFT_N_coreas, "C2RBACKWARD M K");
   	Double_t * FFT_time_coreas = new Double_t [4096*8];
   	double theta_re_coreas[4096*8] = {0.};
   	double theta_im_coreas[4096*8] = {0.};
   	double phi_re_coreas[4096*8] = {0.};
   	double phi_im_coreas[4096*8] = {0.};
   	double FFT_re_coreas[4096*8] = {0.};
   	double FFT_im_coreas[4096*8] = {0.};
   	double binning_MHz = 1./(4096*8*0.1*1.0e-9)*1.0e-6;

void coreas_covolution(double t_E_theta[], double t_E_phi[], double x[], int pos_theta, int pos_phi, int channel, double delay){

	FFT_forward_coreas_theta->SetPoints(t_E_theta);
	FFT_forward_coreas_theta->Transform();
	FFT_forward_coreas_theta->GetPointsComplex(theta_re_coreas, theta_im_coreas);
	
	FFT_forward_coreas_phi->SetPoints(t_E_phi);
	FFT_forward_coreas_phi->Transform();
	FFT_forward_coreas_phi->GetPointsComplex(phi_re_coreas, phi_im_coreas);
		
	for(int k=0; k<2048*8+1; k++){
		if((k<500)||(k>2000)) {
			theta_re_coreas[k] = 0;
			theta_im_coreas[k] = 0;
			phi_re_coreas[k] = 0;
			phi_im_coreas[k] = 0;
		}
		
		else{
			TComplex f_E_theta(theta_re_coreas[k], theta_im_coreas[k]);
			TComplex f_E_phi(phi_re_coreas[k], phi_im_coreas[k]);
			TComplex A_theta;
			TComplex A_phi;
			//H-pol
			if(channel%2==0){
				A_theta = TComplex(ANT_H_theta_re[pos_theta][pos_phi][k], ANT_H_theta_im[pos_theta][pos_phi][k]);
				A_phi = TComplex(ANT_H_phi_re[pos_theta][pos_phi][k], ANT_H_phi_im[pos_theta][pos_phi][k]);
			}
			//V-pol
			else{
				A_theta = TComplex(ANT_V_theta_re[pos_theta][pos_phi][k], ANT_V_theta_im[pos_theta][pos_phi][k]);
				A_phi = TComplex(ANT_V_phi_re[pos_theta][pos_phi][k], ANT_V_phi_im[pos_theta][pos_phi][k]);
			}
			TComplex A_waveform = f_E_theta*A_theta + f_E_phi*A_phi;
			
			TComplex FEE(FEE_re[channel][k], FEE_im[channel][k]);
			TComplex result = A_waveform*FEE;
			FFT_re_coreas[k] = result.Re();
			FFT_im_coreas[k] = result.Im();
			//cout<<k<<"  "<<sqrt(ANT_V_theta_re[18][0][k]*ANT_V_theta_re[18][0][k]+ANT_V_theta_im[18][0][k]*ANT_V_theta_im[18][0][k])<<"  "<<sqrt(FEE_re[0][k]*FEE_re[0][k]+FEE_im[0][k]*FEE_im[0][k])<<endl;
		}
		
	}

	FFT_backward_coreas->SetPointsComplex(FFT_re_coreas, FFT_im_coreas);
	FFT_backward_coreas->Transform();
	FFT_time_coreas = FFT_backward_coreas->GetPointsReal();
	int delay_int = round(delay/0.1);
	for(int i=0; i<1500; i++){
		int j = i*8 - delay_int;
		if(j<0){x[i] = 0.;}
		else{
			x[i] = 1000.*FFT_time_coreas[j]/(4096.*8.);// *1000: V to mV
		}
	}

}






void initial_antenna_direction(){
	// antenna pointing, east, x
	antenna_direction[0][0].SetMagThetaPhi(1, antenna_1_theta, antenna_1_phi);
	antenna_direction[1][0].SetMagThetaPhi(1, antenna_2_theta, antenna_2_phi);
	antenna_direction[2][0].SetMagThetaPhi(1, antenna_3_theta, antenna_3_phi);
	antenna_direction[3][0].SetMagThetaPhi(1, antenna_4_theta, antenna_4_phi);
	
	// antenna upward, z
	antenna_direction[0][2].SetMagThetaPhi(1, antenna_1_theta-90*TMath::DegToRad(), antenna_1_phi);
	antenna_direction[1][2].SetMagThetaPhi(1, antenna_2_theta-90*TMath::DegToRad(), antenna_2_phi);
	antenna_direction[2][2].SetMagThetaPhi(1, antenna_3_theta-90*TMath::DegToRad(), antenna_3_phi);
	antenna_direction[3][2].SetMagThetaPhi(1, antenna_4_theta-90*TMath::DegToRad(), antenna_4_phi);
	
	// antenna north, y; y=zXx
	antenna_direction[0][1] = antenna_direction[0][2].Cross(antenna_direction[0][0]);
	antenna_direction[1][1] = antenna_direction[1][2].Cross(antenna_direction[1][0]);
	antenna_direction[2][1] = antenna_direction[2][2].Cross(antenna_direction[2][0]);
	antenna_direction[3][1] = antenna_direction[3][2].Cross(antenna_direction[3][0]);
	
	//cout<<antenna_direction[0][1].Mag()<<"   "<<antenna_direction[1][1].Mag()<<"   "<<antenna_direction[2][1].Mag()<<"   "<<antenna_direction[3][1].Mag()<<"   "<<endl;
	
	
}

void initial_coreas_coordinate(){
	coreas_x.SetMagThetaPhi(1, 90*TMath::DegToRad(), (90-MagneticInclination)*TMath::DegToRad());
	coreas_y.SetMagThetaPhi(1, 90*TMath::DegToRad(), (180-MagneticInclination)*TMath::DegToRad());
	coreas_z.SetMagThetaPhi(1, 0, 0);
	}

void CoreasAngleToStation(double coreas_theta, double coreas_phi, double &station_theta, double &station_phi){
	station_theta = 90 - coreas_theta; // from zenith to elevation angle
	station_phi = (coreas_phi - 180) + 90 - MagneticInclination;
}


void StationAngleToAntenna(double station_theta, double station_phi, double antenna_theta[4], double antenna_phi[4], TVector3 E_theta_station[4], TVector3 E_phi_station[4]){
	station_theta = 90 -station_theta; // from elevation to zenith angle
	TVector3 station_angle;
	station_angle.SetMagThetaPhi(1, station_theta*TMath::DegToRad(), station_phi*TMath::DegToRad());
	//double antenna_theta[4] = {0.};
	//double antenna_phi[4] = {0.};
	
	for(int i=0; i<4; i++){
		double z_length = station_angle*antenna_direction[i][2];
		antenna_theta[i] = TMath::RadToDeg() * TMath::ACos(z_length);
		TVector3 station_phi_plane = station_angle - z_length*antenna_direction[i][2];
		station_phi_plane = station_phi_plane.Unit();
		antenna_phi[i] = TMath::RadToDeg() * TMath::ATan2(station_phi_plane*antenna_direction[i][1], station_phi_plane*antenna_direction[i][0]);
		
		antenna_theta[i] = 90 - antenna_theta[i]; // from zenith to elevation angle 
		if(antenna_phi[i]<0){antenna_phi[i]+=360;}
		
		//cout<<antenna_theta[i]<<"   "<<antenna_phi[i]<<endl;
		//cout<<station_angle.z()<<"   "<<z_length<<"   "<<antenna_upward[i].z()<<endl;
		//cout<<station_phi_plane.x()<<"  "<<station_phi_plane.y()<<"  "<<station_phi_plane.z()<<endl;
	}
	
	// phi-plane of each antenna is different due to antenna pointing elevation
	for(int i=0; i<4; i++){
		double a_p = (antenna_phi[i] + 90.)*TMath::DegToRad(); // Phi angle for E_phi vector in antenna coordinate.
		E_phi_station[i] = TMath::Cos(a_p)*antenna_direction[i][0] + TMath::Sin(a_p)*antenna_direction[i][1];
		E_phi_station[i].Unit();
		
		// see definition of E_theta in HFSS, it is downward
		E_theta_station[i] = (-1)*station_angle.Cross(E_phi_station[i]);
		//cout<<E_theta_station.Mag()<<"   "<<E_theta_station.x()<<"    "<<E_theta_station.y()<<"   "<<E_theta_station.z()<<endl;
	}
	
	
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
         
         
bool endsWith(const std::string& str, const std::string& suffix) {
    return str.length() >= suffix.length() &&
           str.substr(str.length() - suffix.length()) == suffix;
}



   int N = 1500;
   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C");
      
   TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
   
void filter_fre(double x[1500], bool save_100MHz){

double Re_ant[751];
double Im_ant[751];

   fft_forward->SetPoints(x);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(Re_ant,Im_ant);
   
   // remove DC part
   Re_ant[0] = 0.; Im_ant[0] = 0;

	if(!save_100MHz){
 		for(Int_t k=119; k<122; k++){//remove 100MHz
       		Re_ant[k] = 0.;
       		Im_ant[k] = 0;
                             }
                    }

Double_t * filtered_x = new Double_t [1500];

  fft_back->SetPointsComplex(Re_ant,Im_ant);
  fft_back->Transform();
  filtered_x = fft_back->GetPointsReal();

for(int i=0;i<1500;i++){
x[i] = filtered_x[i]/1500.;
  						}


}
