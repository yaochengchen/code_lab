
#include "../Reconstruction_tchain_new.h"

int pol = 1;



	Double_t trigger_spectrum[8][20] = {0.};
	double pass_trigger_ratio[2] = {0.};
	double eventTime;
	int is_pps = 0;
	int V_bits = 0;
	int H_bits = 0;
	double SNR[8] = {0.};
	double pulse_power[8] = {0.};
    double reconstruct_cross_correlation[2] = {0.};
    double remain_pulse_cross_correlation[2] = {0.};
    double other_pol_cross_correlation[2] = {0.};
    double TDOA[6] = {0.};
    int multi_pulses[8] = {0};
    
    double reconstructed_theta[2] = {0.};
    double reconstructed_phi[2] = {0.};
	
	double *pps_ratio = new double[1000000];
	double *pulser_ratio = new double[1000000];
	double *RF_ratio = new double[1000000];
	double *pps_snr = new double[1000000];
	double *pulser_snr = new double[1000000];
	double *RF_snr = new double[1000000];
	int pps_count = 0;
	int RF_count = 0;
	int pulser_count = 0;
	int gRunNum=0;
	int gEventNum = 0;
	

Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j);
Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j);

void tree_cut(int which_date);


void read_tree(TTree * tree);	


//base center   2021 photogrametry  Tower-4 is back side
double coor_ant[8][3] = {{328914.252, 2697625.953, 709.389},{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949},{328936.041, 2697639.685, 704.949}};	
	
int main(int argc, char *argv[]){

	cout<<argv[1]<<endl;
	string date_name = argv[1];

	int which_date = std::stoi(date_name);
	tree_cut(which_date);
	return 1;
}

void tree_cut(int which_date){
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

TH1F * barson = new TH1F("RF_trigger", "RF_trigger", 400, -20, 20);
TH1F * barson_1 = new TH1F("RF_trigger_1", "RF_trigger_1", 100, 0, 1);

	
	
	
	double *s_time = new double[1000000];
	double *s_theta = new double[1000000];
	double *s_phi = new double[1000000];
	double *s_run = new double[1000000];
	double *s_event = new double[1000000];
	double *s_pol = new double[1000000];
	double *s_x_cor = new double[1000000];
	double *s_p_spec = new double[1000000];
	
	double *w_time = new double[1000000];
	double *w_theta = new double[1000000];
	double *w_phi = new double[1000000];
	double *w_run = new double[1000000];
	double *w_event = new double[1000000];
	
	double *m_time = new double[1000000];
	double *m_theta = new double[1000000];
	double *m_phi = new double[1000000];
	
    const  char * out_filename;
    
    out_filename =  Form("./candidate_events/%08d_H.sh", which_date);
    fstream out_put;
    out_put.open(out_filename, ios::out | ios::trunc);// without ios::out   clean the file before write

	// /media/cyc/For_Linux/TAROGE4_PROCESSED_DATA
	const  char * root_filename =  Form("/media/cyc/Data_disk/Linux/TAROGE4_PROCESSED_DATA/0823t_selection_%08d.root", which_date);
	TFile *RF_file = new TFile(root_filename);
    TTree *Tree_RF = (TTree*) RF_file->Get("selection");
    read_tree(Tree_RF);
    

	//root_filename =  Form("../0627t_selection_%08d.root", which_date);
	//TFile *RF_file_1 = new TFile(root_filename);
    //TTree *Tree_RF_1 = (TTree*) RF_file_1->Get("selection");
    //read_tree(Tree_RF_1);
    
      
    

    int Entries = Tree_RF->GetEntries();
    
    TH1F * RF_trigger = new TH1F("RF_trigger", "RF_trigger", 100, -50, 50);
    TH1F * pps_trigger = new TH1F("pps_trigger", "pps_trigger", 100, -50, 50);
    
    int s_count = 0;
    int w_count = 0;
    int m_count = 0;
    double start_time = 0;
    for(int entry=0; entry<Entries; entry++){
    	
    	//cout<<entry<<endl;
    	Tree_RF->GetEntry(entry);
    	if(entry==0) {start_time = eventTime; cout<<start_time<<endl;}
    	
    	if(gEventNum==6029&&gRunNum==16319){ 
    		reconstructed_theta[0] = 1.95;
    		reconstructed_theta[1] = 1.95;
    		reconstructed_phi[0] = 3.3-26;
    		reconstructed_phi[1] = 3.3-26;
    	}

       const time_t tm_start = (const time_t) eventTime;
       tm *tm_gmt = gmtime(&tm_start);
       int hour = tm_gmt->tm_hour+8;
       if(hour>=24){hour-=24;}
       //cout<<hour<<endl;
    	
    	if(H_bits>V_bits){pol = 0;}
    	else{pol = 1;}
    	//if(gRunNum==19350&&gEventNum==7223) cout<<pol<<"   "<<reconstruct_cross_correlation[pol]<<endl;



   	int i = 0;
    int j = 2;
    double E_TDOA[6] = {0.};
    
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
 	

    	double R = 10000;
    	double in_theta = reconstructed_theta[pol];
    	double in_phi = reconstructed_phi[pol];
    	
    	E_TDOA[pair] = expectedTimeDiff_angle_with_R(in_theta, in_phi, R, coor_ant[i], coor_ant[j]);
    	
	}


    	// boresight -26 degree
    	reconstructed_phi[pol] = reconstructed_phi[pol] + 26;
    	if(reconstructed_phi[pol]>180) {reconstructed_phi[pol] -= 360;}
    	
    	bool TDOA_correct = true;
    	if(pass_trigger_ratio[pol]>0.6&&reconstruct_cross_correlation[pol]>0.6){
			for(int k=0; k<6; k++){
    			//barson->Fill(E_TDOA[k] - TDOA[k]);
    			if(fabs(E_TDOA[k] - TDOA[k])>1.3){TDOA_correct = false;}
    		}
    		
    		
    		//Tree_RF_1->GetEntry(entry);
    		//barson_1->Fill(reconstruct_cross_correlation[pol]);
    		
    	}
    	
    	//selection->Draw("run:SNR[1]", "is_pps==0&&reconstruct_cross_correlation[1]>0.7&&remain_pulse_cross_correlation[1]<0.5&&other_pol_cross_correlation[1]>0.6&&pass_trigger_ratio[1]>0.8")

		double time = eventTime -start_time;
		
		double highest_power = 0;
		double lowest_power = 1.0e12;
		pulse_power[3] = 2*pulse_power[3];
		SNR[3] = sqrt(2)*SNR[3];
		double average_peak = 0.;
		double pulse_cycle = 0.;
		for(int ant=0; ant<4; ant++){
			if(pulse_power[ant]>highest_power){highest_power=pulse_power[ant];}
			if(pulse_power[ant]<lowest_power){lowest_power=pulse_power[ant];}
			average_peak += SNR[ant]/4.;
			
			pulse_cycle += multi_pulses[ant]/4.;
		}
		
		highest_power=pulse_power[3];
		
		
    	
    	//&&hour>=0&&hour<6
    	//&&remain_pulse_cross_correlation[pol]<0.5&&TDOA[1]>-50&&pulse_power[6+pol]/pulse_power[2+pol]<1
    	//if(gRunNum==19350&&gEventNum==7223) cout<<reconstructed_theta[pol]<<"  "<<reconstructed_phi[pol]<<endl;
    	if(pass_trigger_ratio[pol]>0.5&&reconstruct_cross_correlation[pol]>0.6){
    		//&&average_peak>40&&other_pol_cross_correlation[pol]>0.5&&(hour>=-22||hour<=60)&&time>-18000&&time<2160000&&pulse_cycle<300
    		
    	
    		//if(pass_trigger_ratio[pol]<0.7&&reconstruct_cross_correlation[pol]<0.7){
    			w_time[w_count] = eventTime;// - start_time;
    			w_theta[w_count] = reconstructed_theta[pol];
    			w_phi[w_count] = reconstructed_phi[pol];
    			w_run[w_count] = gRunNum;
				w_event[w_count] = gEventNum;
    			w_count++;
    			//}
    		//lowest_power/highest_power>0.5&&&&other_pol_cross_correlation[pol]>0.5&&TDOA_correct&&pulse_cycle<20&&remain_pulse_cross_correlation[pol]<0.5
    		if(lowest_power/highest_power>0.5&&fabs(TDOA[1])<60&&TDOA[5]<-17){
    		//if(true){
    			m_time[m_count] = eventTime;// - start_time;
    			m_theta[m_count] = reconstructed_theta[pol];
    			m_phi[m_count] = reconstructed_phi[pol];
    			m_count++; 
    			   			
    			//if(TDOA_correct){
    			if(pass_trigger_ratio[pol]>0.5&&reconstruct_cross_correlation[pol]>0.7){
    				s_time[s_count] = eventTime;// - start_time;
    				s_theta[s_count] = reconstructed_theta[pol];
    				s_phi[s_count] = reconstructed_phi[pol];
    				
					s_run[s_count] = gRunNum;
					s_event[s_count] = gEventNum;
					s_pol[s_count] = pol;
					s_x_cor[s_count] = reconstruct_cross_correlation[pol];
					s_p_spec[s_count] = pass_trigger_ratio[pol];
					
    				s_count++;
				}
			}
		}


    }


cout<<w_count<<"  "<<m_count<<"    "<<s_count<<endl;

	double f_time[10000] = {0.};
	double f_theta[10000] = {0.};
	double f_phi[10000] = {0.};
	int f_count = 0;
	
for(int k = 0; k<s_count; k++){
	double time = s_time[k];
	double theta = s_theta[k];
	double phi = s_phi[k];
	
	int reject_1 = 0;
	for(int j = 0; j<w_count; j++){
		//+-100s
		if(fabs(w_time[j] - time) < 600){
			if((fabs(w_phi[j] - phi)<5)&&(fabs(w_theta[j] - theta)<5)){
				reject_1 += 1;
			}
		}
	}


	int reject_2 = 0;
	for(int j = 0; j<s_count; j++){
		//+-100s
		if(fabs(s_time[j] - time) < 600){
			//if((fabs(m_phi[j] - phi)<5)&&(fabs(m_theta[j] - theta)<5)){
			//if(fabs(m_phi[j] - phi)<5){
				reject_2 += 1;
			//}
		}
	}

	//out_put<<"root -b -l <<EOF"<<endl;
	//out_put<<".L ../draw_candidate_event.cc+"<<endl;
	
	// itself
	if((reject_1==1)&&(reject_2==1)){
	//if(reject_2==1){
		//out_put<<"root -b -l <<EOF"<<endl;
		//out_put<<".x draw_candidate_event.cc("<<which_date<<", "<<s_run[k]<<", "<<s_event[k]<<")"<<endl;
		//<<", "<<theta<<", "<<phi<<" "<<"     "<<s_pol[k]<<"    "<<s_x_cor[k]<<"    "<<s_p_spec[k]<<endl;
		out_put<<which_date<<"    "<<s_run[k]<<"    "<<s_event[k]<<"    "<<theta<<"    "<<phi<<endl;
		//out_put<<".q"<<endl;
		//out_put<<"EOF"<<endl;
		f_time[f_count] = time;
		f_theta[f_count] = theta;
		f_phi[f_count] = phi;
		f_count++;
			
	}
	//out_put<<".x draw_candidate_event_coherent_sum.cc("<<which_date<<", "<<s_run[k]<<", "<<s_event[k]<<", "<<theta<<", "<<phi<<")"<<endl;

}


	
	for(int j = 0; j<s_count; j++){
		//out_put<<".x draw_candidate_event.cc("<<which_date<<", "<<s_run[j]<<", "<<s_event[j]<<")     "<<", "<<s_theta[j]<<", "<<s_phi[j]<<endl;
		//out_put<<"draw_candidate_event("<<which_date<<", "<<s_run[j]<<", "<<s_event[j]<<")"<<endl;
		//out_put<<"reconstruct_pulser_data(\"/media/cyc/1p9TB/TAROGE4_DATA/data/"<<which_date<<"/run000"<<s_run[j]<<".root\", "<<s_event[j]<<", 100000, "<<s_theta[j]<<", "<<s_phi[j]<<")"<<endl;
		//out_put<<".x draw_candidate_event.cc("<<which_date<<", "<<w_run[j]<<", "<<w_event[j]<<")"<<"  "<<w_time[j] - start_time<<"   "<<w_theta[j]<<"   "<<w_phi[j]<<endl;
	}
	
	//out_put<<".q"<<endl;
	//out_put<<"EOF"<<endl;
/*

TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
c1->cd();
double a[2] = {-180, 180};
double b[2] = {-10, 90};
TGraph *frame = new TGraph(2, a, b);
frame->SetMarkerColor(0);
frame->GetXaxis()->SetRangeUser(-180, 180);
frame->GetYaxis()->SetRangeUser(-30, 90);
frame->Draw("AP");
frame->GetXaxis()->SetTitle("Phi (degree)");
frame->GetYaxis()->SetTitle("Elevation (degree)");
frame->SetTitle("Angular distribution");

//TGraph *cyc = new TGraph(w_count, w_time, w_phi);
TGraph *cyc = new TGraph(w_count, w_phi, w_theta);
cyc->SetMarkerStyle(kCircle);
cyc->SetMarkerSize(0.5);
//cyc->Draw("P same");
//cyc->Draw("AP");


TGraph *ysm = new TGraph(m_count, m_phi, m_theta);
ysm->SetMarkerStyle(kCircle);
ysm->SetMarkerSize(0.5);
ysm->SetMarkerColor(4);
ysm->Draw("P same");

TGraph *wyn = new TGraph(s_count, s_phi, s_theta);
wyn->SetMarkerStyle(kCircle);
wyn->SetMarkerSize(0.5);
wyn->SetMarkerColor(3);
wyn->Draw("P same");


TGraph *hxj = new TGraph(f_count, f_phi, f_theta);
hxj->SetMarkerStyle(kCircle);
hxj->SetMarkerSize(0.5);
hxj->SetMarkerColor(2);
hxj->Draw("P same");

string fig_name = "./angular_distribution_" + std::to_string(which_date) + "_.png";
c1->SaveAs(fig_name.c_str());



TCanvas *c2 = new TCanvas("c2","Canvas Example",1500,1000);
c2->cd();

TGraph *cyc1 = new TGraph(w_count, w_time, w_phi);
cyc1->SetMarkerStyle(kCircle);
cyc1->SetMarkerSize(0.5);
cyc1->Draw("AP");
cyc1->GetXaxis()->SetTimeDisplay(1);
cyc1->GetXaxis()->SetNdivisions(503);
cyc1->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
cyc1->GetXaxis()->SetTimeOffset(60*60*8, "gmt");
cyc1->GetXaxis()->SetTitle("Time");
cyc1->GetYaxis()->SetTitle("Azimuth (degree)");

TGraph *ysm1 = new TGraph(m_count, m_time, m_phi);
ysm1->SetMarkerStyle(kCircle);
ysm1->SetMarkerSize(0.5);
ysm1->SetMarkerColor(4);
ysm1->Draw("P same");
ysm1->GetXaxis()->SetTimeDisplay(1);
ysm1->GetXaxis()->SetNdivisions(503);
ysm1->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
ysm1->GetXaxis()->SetTimeOffset(60*60*8, "gmt");

TGraph *wyn1 = new TGraph(s_count, s_time, s_phi);
wyn1->SetMarkerStyle(kCircle);
wyn1->SetMarkerSize(0.5);
wyn1->SetMarkerColor(3);
wyn1->Draw("P same");
wyn1->GetXaxis()->SetTimeDisplay(1);
wyn1->GetXaxis()->SetNdivisions(503);
wyn1->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
wyn1->GetXaxis()->SetTimeOffset(60*60*8, "gmt");

TGraph *hxj1 = new TGraph(f_count, f_time, f_phi);
hxj1->SetMarkerStyle(kCircle);
hxj1->SetMarkerSize(0.5);
hxj1->SetMarkerColor(2);
hxj1->Draw("P same");
wyn1->GetXaxis()->SetTimeDisplay(1);
wyn1->GetXaxis()->SetNdivisions(503);
wyn1->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
wyn1->GetXaxis()->SetTimeOffset(60*60*8, "gmt");

fig_name = "./phi_vs_time_" + std::to_string(which_date) + "_.png";
c2->SaveAs(fig_name.c_str());




TCanvas *c3 = new TCanvas("c3","Canvas Example",1500,1000);
c3->cd();

TGraph *cyc2 = new TGraph(w_count, w_time, w_theta);
cyc2->SetMarkerStyle(kCircle);
cyc2->SetMarkerSize(0.5);
cyc2->Draw("AP");
cyc2->GetXaxis()->SetTimeDisplay(1);
cyc2->GetXaxis()->SetNdivisions(503);
cyc2->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
cyc2->GetXaxis()->SetTimeOffset(60*60*8, "gmt");
cyc2->GetXaxis()->SetTitle("Time");
cyc2->GetYaxis()->SetTitle("Elevation (degree)");

TGraph *ysm2 = new TGraph(m_count, m_time, m_theta);
ysm2->SetMarkerStyle(kCircle);
ysm2->SetMarkerSize(0.5);
ysm2->SetMarkerColor(4);
ysm2->Draw("P same");
ysm2->GetXaxis()->SetTimeDisplay(1);
ysm2->GetXaxis()->SetNdivisions(503);
ysm2->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
ysm2->GetXaxis()->SetTimeOffset(60*60*8, "gmt");

TGraph *wyn2 = new TGraph(s_count, s_time, s_theta);
wyn2->SetMarkerStyle(kCircle);
wyn2->SetMarkerSize(0.5);
wyn2->SetMarkerColor(3);
wyn2->Draw("P same");
wyn2->GetXaxis()->SetTimeDisplay(1);
wyn2->GetXaxis()->SetNdivisions(503);
wyn2->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
wyn2->GetXaxis()->SetTimeOffset(60*60*8, "gmt");

TGraph *hxj2 = new TGraph(f_count, f_time, f_theta);
hxj2->SetMarkerStyle(kCircle);
hxj2->SetMarkerSize(0.5);
hxj2->SetMarkerColor(2);
hxj2->Draw("P same");
wyn2->GetXaxis()->SetTimeDisplay(1);
wyn2->GetXaxis()->SetNdivisions(503);
wyn2->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
wyn2->GetXaxis()->SetTimeOffset(60*60*8, "gmt");

fig_name = "./theta_vs_time_" + std::to_string(which_date) + "_.png";
c3->SaveAs(fig_name.c_str());
*/
//return 1;
//barson->Draw();
//barson_1->SetLineColor(2);
//barson_1->Draw("same");
cout<<"here"<<endl;
}





void read_tree(TTree * tree){

    tree->SetBranchAddress("eventTime", &eventTime); //PC timestamp
    tree->SetBranchAddress("is_pps", &is_pps);
    tree->SetBranchAddress("V_bits", &V_bits);
    tree->SetBranchAddress("H_bits", &H_bits);
    tree->SetBranchAddress("run", &gRunNum);
    tree->SetBranchAddress("event", &gEventNum);
    
	tree->SetBranchAddress("theta", reconstructed_theta);
	tree->SetBranchAddress("phi", reconstructed_phi);
	tree->SetBranchAddress("TDOA", TDOA);
	tree->SetBranchAddress("multi_pulses", multi_pulses);
	
	
    tree->SetBranchAddress("SNR", SNR);
    tree->SetBranchAddress("pulse_power", pulse_power);
	tree->SetBranchAddress("reconstruct_cross_correlation", reconstruct_cross_correlation);
	tree->SetBranchAddress("remain_pulse_cross_correlation", remain_pulse_cross_correlation);
	tree->SetBranchAddress("other_pol_cross_correlation", other_pol_cross_correlation);
    tree->SetBranchAddress("trigger_spectrum", trigger_spectrum);
    tree->SetBranchAddress("pass_trigger_ratio", pass_trigger_ratio);
    
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

