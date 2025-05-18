void count_trigger(){

TFile *ntuple_file = new TFile("./ntuple.root");

TNtuple * ntuple = (TNtuple*) ntuple_file->Get("ntuple");

float nt_run;
float nt_event;
float nt_time;
float drn_azimuth;
float drn_elevation;
ntuple->SetBranchAddress("run",&nt_run);
ntuple->SetBranchAddress("event",&nt_event);
ntuple->SetBranchAddress("time",&nt_time);
ntuple->SetBranchAddress("drn_azimuth",&drn_azimuth);
ntuple->SetBranchAddress("drn_elevation",&drn_elevation);

float rel_x;
float rel_y;
float rel_z;
ntuple->SetBranchAddress("rel_x",&rel_x);
ntuple->SetBranchAddress("rel_y",&rel_y);
ntuple->SetBranchAddress("rel_z",&rel_z);

float V_peak;
float flag;
ntuple->SetBranchAddress("V_p",&V_peak);
ntuple->SetBranchAddress("flag",&flag);

float distance;
float attenuation;
//ntuple->SetBranchAddress("distance",&distance);
ntuple->SetBranchAddress("attenuation",&attenuation);

float trig_eff;
ntuple->SetBranchAddress("trig_eff",&trig_eff);
float V_peak_4 = 0.;
ntuple->SetBranchAddress("V_p_4",&V_peak_4);
	
const Int_t nt_nentries = (Int_t)ntuple->GetEntries() - 1;

	double attenuation_list[1000] = {0.};
	double efficiency_list[1000] = {0.};
	double V_p_4_list[1000] = {0.};
	double time_list[1000] = {0.};
	int count = 0;
	
	for(int m=0; m<nt_nentries; m++){
      	ntuple->GetEntry(m);
      	if(((int) nt_run)!=2112){continue;}
      	double distance = sqrt(rel_x*rel_x+rel_y*rel_y+rel_z*rel_z);
      	if(nt_time<198||nt_time>405){continue;}
      	
      	attenuation_list[count] = (attenuation-73)*8+20;
      	V_p_4_list[count] = V_peak_4;
      	time_list[count] = nt_time;
      	count++;
    	  	
    }
    
    TGraph * gra = new TGraph(count, time_list, attenuation_list);
    TGraph * grb = new TGraph(count, time_list, V_p_4_list);
    grb->SetMarkerStyle(20);
    
    grb->Draw("AP");
    gra->Draw("PL");
      	
}
