using namespace std;

void test(){


int a[8][20][500];
for (int i = 0; i < 8; ++i){
	for (int j = 0; j < 20; ++j){
		for (int k = 0; k < 100; ++k){
			a[i][j][k] = pow(-1, k) * (i+j);
		}
	}
}

for (int i = 0; i < 8; ++i){
	for (int j = 0; j < 20; ++j){
		sort(a[i][j], a[i][j]+100);
		cout<<a[i][j][50]<<endl;
	}
}		

sort(a[1][1], a[1][1]+100);
		for (int k = 0; k < 100; ++k){
			cout<<a[1][1][k]<<endl;
		}
		
}
/*
int filter_N = 1500;
Int_t filter_n_size = filter_N+1;
TVirtualFFT *filter_fft_forward = TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
TVirtualFFT *filter_fft_back = TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");

void test(){

TRandom3 *r3 = new TRandom3();
   	
double Re_ant[751];
double Im_ant[751];


double filtered_x[1500] = {0.};

Double_t * fft_back_x = new Double_t[1500];

double x[1500] = {0.};
double t[1500] = {0.};
double y[1500] = {0.};

for(Int_t k=0; k<1500; k++){
	t[k] = 0.8*k;
	//x[k] = r3->Gaus(0, 20);

}



		x[500] += -100;
		x[503] += 200;
		x[506] += -300;
		x[509] += 200;
		x[512] += -100;
		x[500] += -100;
		x[503] += 200;
		x[506] += -300;
		x[509] += 200;
		x[512] += -100;
		
   filter_fft_forward->SetPoints(x);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(Re_ant, Im_ant);


 for(Int_t k=0; k<751; k++){//0 means DC component // from 185-350 MHz
       //if((k<222)||(k>420)) {Re_ant[k] = 0.; Im_ant[k] = 0;}
                             }
  
          
double FreqBin = 5/6.;
   int WeakNotchHalfWidth = 3./FreqBin;  //MHz +-
   int StrongNotchHalfWidth = 6./FreqBin;  //MHz +-


  const int NKnownCW = 11;
  double   CWPeak[ NKnownCW ] = {
    165./FreqBin, 187.5/FreqBin, 200./FreqBin, 217.5/FreqBin, 244.167/FreqBin, 249.167/FreqBin, 262./FreqBin, 273.5/FreqBin, 276.65/FreqBin, 299.165/FreqBin, 313.6/FreqBin  //in cal-pulser, not found in forced trig
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


  filter_fft_back->SetPointsComplex(Re_ant,Im_ant);
  filter_fft_back->Transform();
  fft_back_x = filter_fft_back->GetPointsReal();



int Max_i = 0;
double Max = 0.;

for(int i=0;i<1500;i++){
	filtered_x[i] = fft_back_x[i]/1500.;
	y[i] = filtered_x[i] - x[i];
	//cout<<filtered_x[i]<<endl;
  						}

TGraph *cyc = new TGraph(1500, t, x);
cyc->Draw("ALP");

TGraph *wyn = new TGraph(1500, t, filtered_x);
wyn->Draw("LP same");
wyn->SetLineColor(2);
			
}
*/
