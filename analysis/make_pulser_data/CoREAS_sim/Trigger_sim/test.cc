void test(){
	
TFile *ntuple_file = new TFile("./ntuple_20201218.root", "update");

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
auto trig_eff_branch = ntuple->Branch("trig_eff",&trig_eff,"trig_eff/F");
float V_peak_4 = 0.;
auto V_p_branch = ntuple->Branch("V_p_4",&V_peak_4,"V_p_4/F");
	
const Int_t nt_nentries = (Int_t)ntuple->GetEntries() - 1;




	for(int m=0; m<nt_nentries; m++){
      	ntuple->GetEntry(m);
      	
      	distance = sqrt(rel_x*rel_x+rel_y*rel_y+rel_z*rel_z);
      	if(nt_time<962||nt_time>1586){continue;}//198 405  
      	//if((attenuation<72)||(attenuation>74)){continue;}
      	//if(attenuation<76||attenuation>78){continue;}
      	
		V_peak_4 = 1000.;//max[6];
		V_p_branch->Fill();

	trig_eff = 777.;
	
	trig_eff_branch->Fill();
	//ntuple->Fill();
	
	} // m loop
	
	
	cout<<"fuck"<<endl;
	

	
	ntuple->Write("", TObject::kOverwrite); 
	ntuple_file->Close();

}
