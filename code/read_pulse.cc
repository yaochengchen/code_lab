void cal(TGraph * &gr_m, TGraph * &gr_p, char * filename);
void get_band_factor();
void get_sub_band_strength(double x[4096], double sub_band_power[8]);

double template_sub_band_power[8] = {0.};
double template_v_pp = 828.131;
double signal_strength[4096]={0};
double bandfactor[8][4096][2]={0};//8 bands
  Double_t *re_full = new Double_t[4096];
  Double_t *im_full = new Double_t[4096];
  double **v_tmp_re;
  double **v_tmp_im;

  double re_band[8][4096] ={0};

  double im_band[8][4096] ={0};
  
  Int_t n_size = 4096;
  TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &n_size, "R2C ES K");

  TVirtualFFT *fft_back1 = TVirtualFFT::FFT(1, &n_size, "C2CBACKWARD M K");
  TVirtualFFT *fft_back2 = TVirtualFFT::FFT(1, &n_size, "C2CBACKWARD M K");
  TVirtualFFT *fft_back3 = TVirtualFFT::FFT(1, &n_size, "C2CBACKWARD M K");
  TVirtualFFT *fft_back4 = TVirtualFFT::FFT(1, &n_size, "C2CBACKWARD M K");
  TVirtualFFT *fft_back5 = TVirtualFFT::FFT(1, &n_size, "C2CBACKWARD M K");
  TVirtualFFT *fft_back6 = TVirtualFFT::FFT(1, &n_size, "C2CBACKWARD M K");
  TVirtualFFT *fft_back7 = TVirtualFFT::FFT(1, &n_size, "C2CBACKWARD M K");
  TVirtualFFT *fft_back8 = TVirtualFFT::FFT(1, &n_size, "C2CBACKWARD M K");
  
void read_pulse(){

TRandom *r3 = new TRandom3();

FILE *fpinput = fopen("/home/cyc/programming/ICRC/T4-DATA/ICRC_lab_test/0dB_pulse/template.txt","r");
float pulse_shape;
float pulse_time;

double pulse[2600]={0};
double time[4096]={0};

int maximum_point = 0;
float max_value = 0;

if (!fpinput) cout<<"Error while reading template file!!"<<endl;
int pulse_counter=0;
while(1){
 int p = fscanf(fpinput,"%f %*f %f %*f %*f", &pulse_time, &pulse_shape);//read channel B data, it is most stable
if (p!=2) break;

//cout<<pulse_shape<<endl;
pulse[pulse_counter] = pulse_shape;


      if (fabs(pulse[pulse_counter])>max_value) {
        maximum_point = pulse_counter; max_value = fabs(pulse[pulse_counter]);
               									}       
                                          
time[pulse_counter] = pulse_time;
pulse_counter+=1;

}
//cout<<pulse_counter<<" max: "<<maximum_point<<endl;
fclose(fpinput);

 TGraph * pulse_template = new TGraph(pulse_counter,time,pulse);
 //pulse_template->Draw();
 
 //double signal_strength[4096]={0};
 for(int i=0;i<105;i++){//400 points,
     
     signal_strength[i+500]=pulse_template->Eval(time[maximum_point]+0.8*(i-25.));// 0.8 ns 
 //    cout<<signal_strength[i]<<endl; 
//     cout<<maximum_point<<endl;
//     cout<<pulse[maximum_point]<<endl;
	//for obtain V_pp of template accurately.
	//cout<<pulse_template->Eval(2+i*0.005)<<"  "<<pulse_template->Eval(5.8+i*0.005)<<endl;
	//cout<<pulse_template->Eval(3.8+i*0.005)<<endl<<endl;
                       }
                       

	for(int i=0; i<4096; i++){
		time[i] = i*0.8;
	}
 TGraph * pulse_template1 = new TGraph(4096,time,signal_strength);
 pulse_template->Draw();
 
 

    v_tmp_re = new double*[8];
  for(int i=0;i<8;i++){
    v_tmp_re[i] = new double[4096];
  }

    v_tmp_im = new double*[8];
  for(int i=0;i<8;i++){
    v_tmp_im[i] = new double[4096];
  }            
            
	get_band_factor();
	
	get_sub_band_strength(signal_strength, template_sub_band_power);
	//cout<<template_sub_band_power[7]<<endl;
            
 }
 
 
 
 
 
 
 
void cal(TGraph * &gr_m, TGraph * &gr_p, char * filename){
  FILE *fpin = fopen(filename,"r"); 
  float frq;
  float mig;
  float pha;
  int p=1;
  double a[1001];
  double b[1001];
  double c[1001];

  int i=0;
  double offset = 0.;
  while(1){
    p = fscanf(fpin,"%f %*f %*f %f %f %*f %*f %*f %*f", &frq, &mig, &pha);
    pha += offset;
    if(i>0){
    	if(pha>c[i-1]){
    		offset -= 360.;
    	}
    }
    a[i]=frq/(1.0e6);
    b[i]=mig;
    c[i]=pha;
    
    i+=1;
    if (p!=3) break;
  }

  fclose(fpin);
  gr_m = new TGraph(i,a,b);
  gr_p = new TGraph(i,a,c);
  
}


// -----------------------------get band factors
void get_band_factor(){


  char * filename;


  TGraph * gr1_m;
  TGraph * gr1_p;
  filename = Form("/home/cyc/programming/Thesis/band_filter/sawfilters/20170516-sawfilter-16-avg/TA0574.S2P");;
  cal(gr1_m, gr1_p, filename);



  TGraph * gr2_m;
  TGraph * gr2_p;
  filename = Form("/home/cyc/programming/Thesis/band_filter/sawfilters/20170516-sawfilter-16-avg/TA0575.S2P");
  cal(gr2_m, gr2_p, filename);

  TGraph * gr3_m;
  TGraph * gr3_p;
  filename = Form("/home/cyc/programming/Thesis/band_filter/sawfilters/20170516-sawfilter-16-avg/TA0369.S2P");
  cal(gr3_m, gr3_p, filename);


  TGraph * gr4_m;
  TGraph * gr4_p;
  filename = Form("/home/cyc/programming/Thesis/band_filter/sawfilters/20170516-sawfilter-16-avg/TA1345.S2P");
  cal(gr4_m, gr4_p, filename);

 
  TGraph * gr5_m;
  TGraph * gr5_p;
  filename = Form("/home/cyc/programming/Thesis/band_filter/sawfilters/20170516-sawfilter-16-avg/TA0560.S2P");
  cal(gr5_m, gr5_p, filename);

  TGraph * gr6_m;
  TGraph * gr6_p;
  filename = Form("/home/cyc/programming/Thesis/band_filter/sawfilters/20170516-sawfilter-16-avg/TA0510.S2P");
  cal(gr6_m, gr6_p, filename);




  TGraph * gr7_m;
  TGraph * gr7_p;
  filename = Form("/home/cyc/programming/Thesis/band_filter/sawfilters/20170516-sawfilter-16-avg/TA0544.S2P");
  cal(gr7_m, gr7_p, filename);


  TGraph * gr8_m;
  TGraph * gr8_p;
  filename = Form("/home/cyc/programming/Thesis/band_filter/sawfilters/20170516-sawfilter-16-avg/TA0797.S2P");
  cal(gr8_m, gr8_p, filename);








  for(int i=1;i<2049;i++)
    {
      bandfactor[0][i][0] = sqrt(pow(10,((gr1_m->Eval((i)/3.2768))/10.)));
      bandfactor[1][i][0] = sqrt(pow(10,((gr2_m->Eval((i)/3.2768))/10.)));
      bandfactor[2][i][0] = sqrt(pow(10,((gr3_m->Eval((i)/3.2768))/10.)));
      bandfactor[3][i][0] = sqrt(pow(10,((gr4_m->Eval((i)/3.2768))/10.)));
      bandfactor[4][i][0] = sqrt(pow(10,((gr5_m->Eval((i)/3.2768))/10.)));
      bandfactor[5][i][0] = sqrt(pow(10,((gr6_m->Eval((i)/3.2768))/10.)));
      bandfactor[6][i][0] = sqrt(pow(10,((gr7_m->Eval((i)/3.2768))/10.)));
      bandfactor[7][i][0] = sqrt(pow(10,((gr8_m->Eval((i)/3.2768))/10.)));

      bandfactor[0][i][1] = gr1_p->Eval((i)/3.2768)*TMath::DegToRad();
      bandfactor[1][i][1] = gr2_p->Eval((i)/3.2768)*TMath::DegToRad();
      bandfactor[2][i][1] = gr3_p->Eval((i)/3.2768)*TMath::DegToRad();
      bandfactor[3][i][1] = gr4_p->Eval((i)/3.2768)*TMath::DegToRad();
      bandfactor[4][i][1] = gr5_p->Eval((i)/3.2768)*TMath::DegToRad();
      bandfactor[5][i][1] = gr6_p->Eval((i)/3.2768)*TMath::DegToRad();
      bandfactor[6][i][1] = gr7_p->Eval((i)/3.2768)*TMath::DegToRad();
      bandfactor[7][i][1] = gr8_p->Eval((i)/3.2768)*TMath::DegToRad();
      
    }
    
    for(int j=0;j<8;j++)
    {
    	for(int i=1;i<2049;i++)
    	{
    		// -5dB
    		if(bandfactor[j][i][0]<0.56234133){bandfactor[j][i][0]=0.;}
    	}
    }







}




void get_sub_band_strength(double x[4096], double sub_band_power[8])  
{ 
  fft_own->SetPoints(x);
  fft_own->Transform();


  fft_own->GetPointsComplex(re_full,im_full);


  for(int j=0;j<8;j++){
    for(int i=1;i<2049;i++)
      {
      
      	double magnitude = sqrt(re_full[i]*re_full[i] + im_full[i]*im_full[i]);
      	double phase = TMath::ATan2(im_full[i], re_full[i]);
      	
      	magnitude *= bandfactor[j][i][0];
      	phase += bandfactor[j][i][1];
      	
      	re_band[j][i] = magnitude * TMath::Cos(phase);
      	im_band[j][i] = magnitude * TMath::Sin(phase);
      	
	//re_band[j][i]= 2*re_full[i]*bandfactor[j][i][0];

	//im_band[j][i]= 2*im_full[i]*bandfactor[j][i][0];
      }
      re_band[j][0] = 0.;
      im_band[j][0] = 0.;
  }

  
  fft_back1->SetPointsComplex(re_band[0],im_band[0]);
  fft_back1->Transform();
  //cout<<"3"<<endl;
  fft_back1->GetPointsComplex(v_tmp_re[0],v_tmp_im[0]);
  //cout<<"4"<<endl;
  
  fft_back2->SetPointsComplex(re_band[1],im_band[1]);
  fft_back2->Transform();
  fft_back2->GetPointsComplex(v_tmp_re[1],v_tmp_im[1]);

  fft_back3->SetPointsComplex(re_band[2],im_band[2]);
  fft_back3->Transform();
  fft_back3->GetPointsComplex(v_tmp_re[2],v_tmp_im[2]);

  fft_back4->SetPointsComplex(re_band[3],im_band[3]);
  fft_back4->Transform();
  fft_back4->GetPointsComplex(v_tmp_re[3],v_tmp_im[3]);

  fft_back5->SetPointsComplex(re_band[4],im_band[4]);
  fft_back5->Transform();
  fft_back5->GetPointsComplex(v_tmp_re[4],v_tmp_im[4]);

  fft_back6->SetPointsComplex(re_band[5],im_band[5]);
  fft_back6->Transform();
  fft_back6->GetPointsComplex(v_tmp_re[5],v_tmp_im[5]);

  fft_back7->SetPointsComplex(re_band[6],im_band[6]);
  fft_back7->Transform();
  fft_back7->GetPointsComplex(v_tmp_re[6],v_tmp_im[6]);

  fft_back8->SetPointsComplex(re_band[7],im_band[7]);
  fft_back8->Transform();
  fft_back8->GetPointsComplex(v_tmp_re[7],v_tmp_im[7]);

  for(int j=0;j<8;j++){
  	/*
    for(int i=200;i<1300;i++)
      {
 double power_mw = (v_tmp_re[j][i]*v_tmp_re[j][i] + v_tmp_im[j][i]*v_tmp_im[j][i])/(4096.*4096.*50.*1000.); 
//4096:normalize, 50:resistance, 1000:change to mW
sub_band_power[j] += power_mw;
      }
     */
      
    double max = 0;
    int max_i = 0;
    for(int i=200;i<1300;i++)
      {
 		double power_mw = (v_tmp_re[j][i]*v_tmp_re[j][i] + v_tmp_im[j][i]*v_tmp_im[j][i])/(4096.*4096.*50.*1000.); 
//4096:normalize, 50:resistance, 1000:change to mW
		if(max<power_mw){max=power_mw; max_i=i;}
      }
      	//sub_band_power[j] = max;
    
    for(int i=0;i<1500;i++){
    	if((i<(max_i-31))||(i>(max_i+31))){continue;}
    	double power_mw = (v_tmp_re[j][i]*v_tmp_re[j][i] + v_tmp_im[j][i]*v_tmp_im[j][i])/(4096.*4096.*50.*1000.);
    	sub_band_power[j] += power_mw;
    }
     
                      }
}
