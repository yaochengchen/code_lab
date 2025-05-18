#include "dirent.h"
//function to calculate dot product of two vectors
double dot_product(double vector_a[3], double vector_b[3]) {
   double product = 0;
   for (int i = 0; i < 3; i++)
   product = product + vector_a[i] * vector_b[i];
   return product;
}
//function to calculate cross product of two vectors
void cross_product(double vector_a[3], double vector_b[3], double temp[3]) {
   temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
   temp[1] = (-1)*(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
   temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
}

struct sample_data {
    int radius = 0;
    int angle = 0;
    double strength = 0.;
};

bool mycompare(sample_data s1, sample_data s2){
    return s1.radius > s2.radius;
	}

void get_candidate_file(string path, vector<string> &filenames);



   	int FFT_N = 4096;
   	Int_t iFFT_N = FFT_N;		
   	TVirtualFFT *FFT_forward = TVirtualFFT::FFT(1, &FFT_N,"R2C ES K");
    
    TVirtualFFT *FFT_backward = TVirtualFFT::FFT(1, &iFFT_N, "C2RBACKWARD M K");
   	Double_t * FFT_time = new Double_t [4096];
   	double FFT_re[4096] = {0.};
   	double FFT_im[4096] = {0.};
   	double binning_MHz = 1./(4096*0.1*1.0e-9)*1.0e-6;

void filtering(double x[]);

void draw_cone(){



  double n_mag[3] = {36.4, 0, -25.8};//36.4, 0, -25.8
  double scale_mag = 0;
  for(int i=0; i<3; i++){
  	scale_mag += n_mag[i]*n_mag[i];
  }
  scale_mag = sqrt(scale_mag);
  for(int i=0; i<3; i++){
  	n_mag[i] /= scale_mag;
  }

  	double theta = (90-40)*TMath::DegToRad();
  	double phi = (90-180)*TMath::DegToRad();
  	double n_comein[3] = {0.};
	n_comein[0] = (-1)*TMath::Cos(theta)*TMath::Cos(phi);
	n_comein[1] = (-1)*TMath::Cos(theta)*TMath::Sin(phi);
	n_comein[2] = (-1)*TMath::Sin(theta);

	double P_geoemi[3] = {0.};
	cross_product(n_comein, n_mag, P_geoemi);
  double scale_P_geoemi = 0;
  for(int i=0; i<3; i++){
  	scale_P_geoemi += P_geoemi[i]*P_geoemi[i];
  }
  scale_P_geoemi = sqrt(scale_P_geoemi);
  for(int i=0; i<3; i++){
  	P_geoemi[i] /= scale_P_geoemi;
  	cout<<P_geoemi[i]<<endl;
  }

  double V_geoemi[3] = {0.};
  cross_product(n_comein, P_geoemi, V_geoemi);
  //cout<<sqrt(V_geoemi[0]*V_geoemi[0]+V_geoemi[1]*V_geoemi[1]+V_geoemi[2]*V_geoemi[2])<<endl;
  //cout<<V_geoemi[0]<<"  "<<V_geoemi[1]<<"  "<<V_geoemi[2]<<endl;
  double radiactive[3] = {(P_geoemi[0]+V_geoemi[0])/sqrt(2.), (P_geoemi[1]+V_geoemi[1])/sqrt(2.), (P_geoemi[2]+V_geoemi[2])/sqrt(2.)};




	double time[8000] = {0.};
	double E_x[8000] = {0.};
	double E_y[8000] = {0.};
	double E_z[8000] = {0.};

	double _time = 0.;
	double _E_x= 0.;
	double _E_y= 0.;
	double _E_z= 0.;
	TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1500);
	TCanvas *c2 = new TCanvas("c2","Canvas Example",1500,1500);
	//TCanvas *c3 = new TCanvas("c3","Canvas Example",1500,1500);

	TH2F* cone = new TH2F("cone","Polar coordinates",560,-5610,5590,560,-5610,5590);

	// check read_coreas.cc 353 lines; check until 56000
	string path = "/media/cyc/For_Linux/CoREAS_Sim/proton/89015/SIM000001_coreas";
	string listname = "/media/cyc/For_Linux/CoREAS_Sim/proton/10001/SIM000001_coreas/raw_pos_750_135.dat";
	
	vector<string> filenames;

	get_candidate_file(path, filenames);

	sample_data sample_strength[1600];
	int count_sample = 0;
	std::set<int> sample_radius;
	std::set<int> sample_angle;

	for(auto filename : filenames){
		//cout<<filename<<endl;
		const char * c = filename.c_str();
		FILE *fpin = fopen(c,"r");

  		int count=0;
  		while(1){
    		int p = fscanf(fpin,"%lf %lf %lf %lf", &_time, &_E_x, &_E_y, &_E_z);
    
			//double E_strength = sqrt(_E_x*_E_x);//sqrt(_E_x*_E_x +_E_y*_E_y +_E_z*_E_z);
			//if(E_strength>fabs(max_strength)){max_strength=_E_x;}//E_strength

			//double E_field[3] = {_E_x, _E_y, _E_z};
			//double E_strength = dot_product(E_field, P_geoemi);
			//if(fabs(E_strength)>fabs(max_strength)){max_strength=fabs(E_strength);}
    		E_x[count] = _E_x;
    		E_y[count] = _E_y;
    		E_z[count] = _E_z;
    		
    		if (p!=4) break;
			count++;
  		}


  		filtering(E_x);
  		filtering(E_y);
  		filtering(E_z);

		double max_strength = 0.;
  		for(int i=0; i<count; i++){

			double E_field[3] = {E_x[i], E_y[i], E_z[i]};
			//double E_strength = (dot_product(E_field, P_geoemi));
			double E_strength = sqrt(E_x[i]*E_x[i] + E_y[i]*E_y[i] + E_z[i]*E_z[i]);
			if(fabs(E_strength)>fabs(max_strength)){max_strength=E_strength;}

  		}   		


  		//cout<<max_strength<<endl;

  		int nPos_2 = filename.find_last_of("_");
  		int nPos_1 = filename.find_last_of("_", nPos_2-1);
  		int nPos_3 = filename.find(".dat");
  		//cout<<nPos_2<<"   "<<nPos_1<<endl;
  		double radius = (double) std::stoi(filename.substr(nPos_1+1, nPos_2-nPos_1-1));
  		double angle = (double) std::stoi(filename.substr(nPos_2+1, nPos_3-nPos_2-1));

  		sample_radius.insert((int) radius);
  		sample_angle.insert((int) angle);

  		sample_strength[count_sample].radius = radius;
  		sample_strength[count_sample].angle = angle;
  		sample_strength[count_sample].strength = max_strength;
  		count_sample++;

  		//cout<<angle<<"  ";
  		angle = angle*TMath::DegToRad();
  		//cout<<radius<<" "<<angle<<endl;
  		double x = radius*TMath::Cos(angle);
  		double y = radius*TMath::Sin(angle);
  		//cout<<x<<"  "<<y<<"   "<<angle<<"  "<<radius<<endl;
  		cone->SetBinContent(round(x/20.+560)+1, round(y/20.+560)+1, max_strength);
  		//cout<<round(x+1400)+1<<"  "<<round(y+1400)+1<<endl;

	}


	std::sort(sample_strength, sample_strength+1600, mycompare);
	

	sample_data axis_strength[1600];
	int count_axis = 0;
	std::set<int> axis_radius;
	std::set<int> axis_angle;

	for(int sample_angle=0; sample_angle<45; sample_angle+=45){

		double sample_raduis[20] = {0.};
		double sample_voltage[20] = {0.};
		int interpolate_count = 0;
		for(int i=0; i<count_sample; i++){
			int angle = sample_strength[i].angle;
			//cout<<angle<<endl;
			if(angle==sample_angle){
				sample_raduis[interpolate_count] = sample_strength[i].radius;
				sample_voltage[interpolate_count] = sample_strength[i].strength;
				interpolate_count++;
				//cout<<interpolate_count<<endl;
			}	

		}
		TGraph *along_radius = new TGraph(20, sample_raduis, sample_voltage);	

		double interpolate_radius[200] = {0.};
		double interpolate_voltage[200] = {0.};	

		for(int i=0; i<200; i++){
			double radius = 30+(1000.-30.)/200.*i;
			interpolate_radius[i] = radius;
			interpolate_voltage[i] = along_radius->Eval(radius,0,"S");

			axis_strength[count_axis].radius = interpolate_radius[i];
			axis_strength[count_axis].angle = sample_angle;
			axis_strength[count_axis].strength = interpolate_voltage[i];

			axis_radius.insert((int) radius);
  			axis_angle.insert((int) sample_angle);

			count_axis++;
		}
		//TGraph *along_radius_interpolate = new TGraph(200, interpolate_radius, interpolate_voltage);	
	

		//along_radius->Draw("APL");
		//along_radius_interpolate->SetLineColor(2);
		//along_radius_interpolate->Draw("same PL");
	}


	std::vector<int> list_radius;
	std::vector<int> list_angle;
	for(auto radius : sample_radius){
		list_radius.push_back(radius);
	}

	for(auto angle : sample_angle){
		list_angle.push_back(angle);
		//cout<<angle<<endl;
	}
	list_angle.push_back(360);

	for(int i=0; i<560; i++){
		for(int j=0; j<560; j++){
			double x = 20*i-5600;
			double y = 20*j-5600;
			double radius = std::sqrt(x*x + y*y);
      		double angle = std::atan2(y,x)*TMath::RadToDeg();
      		if(angle<0){angle+=360;}
      		//cout<<angle<<endl;

      		int l_radius = 0;
      		int l_angle = 0;
      		int h_radius = 0;
      		int h_angle = 0;
      		for(int k=0; k<list_radius.size(); k++){
      			int l=k+1;
      			//cout<<angle<<list_angle[k]<<list_angle[l]<<endl;
      			if((list_radius[k]<=radius)&&(list_radius[l]>=radius)){
      				l_radius = list_radius[k];
      				h_radius = list_radius[l];
      			}
      		}

      		for(int k=0; k<list_angle.size(); k++){
      			int l=k+1;
      			if((list_angle[k]<=angle)&&(list_angle[l]>=angle)){
      				l_angle = list_angle[k];
      				h_angle = list_angle[l];
      			}

      		}
      		

      		//if(l_radius==0){continue;}

      		double ll, lh, hl, hh = 0.;// radius angle
      		for(int k=0; k<count_sample; k++){
      			int radius = sample_strength[k].radius;
      			int angle = sample_strength[k].angle;
      			if((l_radius==radius)&&(l_angle==angle)){ll = sample_strength[k].strength;}
      			if((l_radius==radius)&&((h_angle%360)==angle)){lh = sample_strength[k].strength;}
      			if((h_radius==radius)&&(l_angle==angle)){hl = sample_strength[k].strength;}
      			if((h_radius==radius)&&((h_angle%360)==angle)){hh = sample_strength[k].strength;}
      		}
      		//cout<<ll<<" "<<lh<<" "<<hl<<" "<<hh<<endl;


      		double scale_radius_l = (h_radius-radius) / (h_radius - l_radius);
      		double scale_radius_h = (radius-l_radius) / (h_radius - l_radius);

      		double scale_angle_l = (h_angle-angle) / (h_angle - l_angle);
      		double scale_angle_h = (angle-l_angle) / (h_angle - l_angle);

      		double strength = ll*scale_radius_l*scale_angle_l + lh*scale_radius_l*scale_angle_h + hl*scale_radius_h*scale_angle_l + hh*scale_radius_h*scale_angle_h;
      		//cout<<ll<<" "<<lh<<" "<<hl<<" "<<hh<<endl;
      		//<<l_radius<<" "<<h_radius<<" "<<l_angle<<" "<<h_angle<<endl;
      		//cout<<strength<<endl;
			cone->SetBinContent(i+1, j+1, strength);
		}
	}
	c1->cd();
	cone->Draw("colz");
	cone->SetTitle("radiation pattern");
	cone->GetXaxis()->SetTitle("vXB");
	cone->GetYaxis()->SetTitle("vX(vXB)");
	cone->SetStats(0);
	//cone->GetXaxis()->SetRangeUser(-200, 200);
	//cone->GetYaxis()->SetRangeUser(-200, 200);




  
  const char * c = listname.c_str();
  double factor = 1;

/*
  //TH2F* delay_cone = new TH2F("delay cone","Polar coordinates",60,-150,150, 60,-150,150);
  TH2F* delay_cone = new TH2F("delay cone","Polar coordinates",3000,-1500.5,1499.5, 3000,-1500.5,1499.5);
  for(auto filename : filenames){
	//cout<<filename<<endl;
	const char * c = filename.c_str();
    int nPos_2 = filename.find_last_of("_");
  	int nPos_1 = filename.find_last_of("_", nPos_2-1);
  	int nPos_3 = filename.find(".dat");
  	//cout<<nPos_2<<"   "<<nPos_1<<endl;
  	double radius = (double) std::stoi(filename.substr(nPos_1+1, nPos_2-nPos_1-1));
  	double angle = (double) std::stoi(filename.substr(nPos_2+1, nPos_3-nPos_2-1));

  	if(angle==0||angle==180){continue;}
  	if(radius>100||radius<0){continue;}
  	double factor = 1;
  	if(angle>180){factor = -1;}
*/

  FILE *fpin = fopen(c,"r");

  int count=0;
  double cE_x[8000] = {0.};
  double cE_y[8000] = {0.};
  double cE_z[8000] = {0.};
  while(1){
    int p = fscanf(fpin,"%lf %lf %lf %lf", &_time, &_E_x, &_E_y, &_E_z);
    
	time[count] = _time;
	cE_x[count] = _E_x;
	cE_y[count] = _E_y;
	cE_z[count] = _E_z;
    
	
    if (p!=4) break;
	count++;
  }
  
  fclose(fpin);

  filtering(cE_x);
  filtering(cE_y);
  filtering(cE_z);

  TGraph * G_z = new TGraph(count, time, cE_z);
  G_z->SetLineColor(3);
  

  TGraph * G_y = new TGraph(count, time, cE_y);
  G_y->SetLineColor(2);
  

  TGraph * G_x = new TGraph(count, time, cE_x);
  //G_x->SetLineColor(3);
  
  c2->cd();
  G_z->Draw("apl");
  G_x->Draw("same pl");
  G_y->Draw("same pl");


  double mag_emi[8000] = {0.};
  double ele_emi[8000] = {0.};
  double E_v[8000] = {0.};
  for(int i=0; i<count; i++){
  	double E_field[3] = {E_x[i], E_y[i], E_z[i]};
  	mag_emi[i] = dot_product(E_field, P_geoemi);
  	ele_emi[i] = factor*dot_product(E_field, V_geoemi);
  	E_v[i] = dot_product(E_field, n_comein);//radiactive
  }

  //filtering(mag_emi);
	double frequency[72] = {0.};
	double G_phase[72] = {0.};
	double G_magnitude[72] = {0.};
	for(int i=73; i<145; i++){
		frequency[i-73] = i*binning_MHz;
		G_phase[i-73] = TMath::ATan2(FFT_im[i], FFT_re[i]);
		G_magnitude[i-73] = sqrt(pow(FFT_re[i],2) + pow(FFT_im[i],2));
	}

  //filtering(ele_emi);
	double A_phase[72] = {0.};
	double A_magnitude[72] = {0.};
	for(int i=73; i<145; i++){
		A_phase[i-73] = TMath::ATan2(FFT_im[i], FFT_re[i]);
		A_magnitude[i-73] = sqrt(pow(FFT_re[i],2) + pow(FFT_im[i],2));
	}

	double a = 0;
	TH1F * cyc = new TH1F("h1", "h1", 400, -2, 2);
	for(int i=0; i<72; i++){
		double delta_phase = G_phase[i] - A_phase[i];
		if(delta_phase>1*TMath::Pi()){delta_phase -= 2*TMath::Pi();}
		if(delta_phase<-1*TMath::Pi()){delta_phase += 2*TMath::Pi();}
		if(delta_phase>1*TMath::Pi()){delta_phase -= 2*TMath::Pi();}
		if(delta_phase<-1*TMath::Pi()){delta_phase += 2*TMath::Pi();}
		if(delta_phase>1*TMath::Pi()){delta_phase -= 2*TMath::Pi();}
		if(delta_phase<-1*TMath::Pi()){delta_phase += 2*TMath::Pi();}
		if(delta_phase>1*TMath::Pi()){delta_phase -= 2*TMath::Pi();}
		if(delta_phase<-1*TMath::Pi()){delta_phase += 2*TMath::Pi();}
		double delta_t = delta_phase/(2*TMath::Pi()*frequency[i]*1.0e6)*1.0e9;
		//cout<<delta_t<<"  "<<delta_phase<<endl;
		a+=delta_t;
		cyc->Fill(delta_t);
	}
	double delay = a/72.;
	cout<<delay<<endl;

	//c3->cd();
	//cyc->Draw();

  //filtering(E_v);

  //c2->cd();

  TGraph * G_mag = new TGraph(count, time, mag_emi);
  
  TGraph * G_ele = new TGraph(count, time, ele_emi);
  G_ele->SetLineColor(2);
  //G_ele->Draw("apl");
  G_mag->GetXaxis()->SetRangeUser(-5e-9, 20e-9);
  G_ele->GetYaxis()->SetRangeUser(-50e-8, 50e-8);
  //G_mag->Draw("same pl");

  TGraph * G_v = new TGraph(count, time, E_v);
  G_v->SetLineColor(3);
  //G_v->Draw("same pl");

	//G_ele->SetTitle("radiation pattern");
	//G_ele->GetXaxis()->SetTitle("vXB");
	//G_ele->GetYaxis()->SetTitle("(vXB)Xv");


/*
  		//cout<<angle<<"  ";
  		angle = angle*TMath::DegToRad();
  		//cout<<radius<<" "<<angle<<endl;
  		double x = radius*TMath::Cos(angle);
  		double y = radius*TMath::Sin(angle);
  		//cout<<x<<"  "<<y<<"   "<<angle<<"  "<<radius<<endl;
  		int i_x = round(x+1500.);
  		int j_y = round(y+1500.);
  		for(int i=i_x-2; i<i_x+3; i++){
  			for(int j=j_y-2; j<j_y+3; j++){
  				delay_cone->SetBinContent(i+1, j+1, delay);
  			}
  		}
  		
}
	c3->cd();
	delay_cone->Draw("colz");
	delay_cone->SetTitle("Delay b.t. Aska. & Geom.");
	delay_cone->GetXaxis()->SetTitle("vXB");
	delay_cone->GetYaxis()->SetTitle("vX(vXB)");
	delay_cone->SetStats(0);
	delay_cone->GetXaxis()->SetRangeUser(-200, 200);
	delay_cone->GetYaxis()->SetRangeUser(-200, 200);
*/
}





void get_candidate_file(string path, vector<string> &filenames){
	DIR *pDir;
	struct dirent* ptr;
	if(!(pDir = opendir(path.c_str()))){
		cout<<"Folder doesn't Exist!"<<endl;
		return;
	}
	//cout<<"open folder"<<endl;
	//cout<<readdir(pDir)<<endl;
	while((ptr = readdir(pDir))!=0){
		if(std::strstr(ptr->d_name, ".dat")&&!std::strstr(ptr->d_name, ".dat~")){
			filenames.push_back(path+"/"+ptr->d_name);
			//cout<<ptr->d_name<<endl;
		}
	}
	closedir(pDir);
}





void filtering(double x[]){

		FFT_forward->SetPoints(x);
		FFT_forward->Transform();
		FFT_forward->GetPointsComplex(FFT_re, FFT_im);
		

		for(int i=0; i<4096; i++){// 73 144
			if((i>145)||(i<73)){//(i>145)||(i<73);(i>33)
				FFT_re[i] = 0.;
				FFT_im[i] = 0.;
			}
		}

		FFT_backward->SetPointsComplex(FFT_re, FFT_im);
		FFT_backward->Transform();
		FFT_time = FFT_backward->GetPointsReal();
		for(int i=0; i<4096; i++){
			x[i] = FFT_time[i]/4096.;
		}

}
