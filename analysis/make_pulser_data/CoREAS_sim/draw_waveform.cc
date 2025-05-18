void draw_waveform(){
	int run = 22008;
	double scale_factor = 1.;
	
  double radius = 0;
  double pattern_phi = 0;
  double ShowerTheta;
  double ShowerPhi;
  Short_t coreas_T1H[1500];
  Short_t coreas_T2H[1500];
  Short_t coreas_T3H[1500];
  Short_t coreas_T4H[1500];
  Short_t coreas_T1V[1500];
  Short_t coreas_T2V[1500];
  Short_t coreas_T3V[1500];
  Short_t coreas_T4V[1500];
	
	const string coreas_root = Form("/media/cyc/For_Linux/CoREAS_Sim/proton/cyc/Coreas-t4-r%06d.root", run);
	//TFile *coreas_file = new TFile(coreas_root.c_str());
	TFile *coreas_file = new TFile("/home/cyc/software/TAROGE-4_analysis/make_pulser_data/data/make_coreas_data_with_pps_real1775.000000.root");
	//TFile *coreas_file = new TFile("/media/cyc/For_Linux/CoREAS_Sim/proton/cyc/Coreas-t4-r000001.root");
	TTree *EfieldTree = (TTree*) coreas_file->Get("t");
	
	//TTree *EfieldTree = (TTree*) coreas_file->Get("EfieldTree");
	EfieldTree->SetBranchAddress("radius", &radius);
	EfieldTree->SetBranchAddress("phi", &pattern_phi);
	EfieldTree->SetBranchAddress("ShowerTheta", &ShowerTheta);
	EfieldTree->SetBranchAddress("ShowerPhi", &ShowerPhi);
	
	
	EfieldTree->SetBranchAddress("T1V", coreas_T1V);
	EfieldTree->SetBranchAddress("T2V", coreas_T2V);
	EfieldTree->SetBranchAddress("T3V", coreas_T3V);
	EfieldTree->SetBranchAddress("T4V", coreas_T4V);
	EfieldTree->SetBranchAddress("T1H", coreas_T1H);
	EfieldTree->SetBranchAddress("T2H", coreas_T2H);
	EfieldTree->SetBranchAddress("T3H", coreas_T3H);
	EfieldTree->SetBranchAddress("T4H", coreas_T4H);
	
	EfieldTree->GetEntry(5200);
	
		double x_vpp[8][4096] = {0.};
		double t[1500] = {0.};
		for(int i=0; i<1500; i++){
			x_vpp[0][i] = coreas_T1H[i] * scale_factor;
			x_vpp[2][i] = coreas_T2H[i] * scale_factor;
			x_vpp[4][i] = coreas_T3H[i] * scale_factor;
			x_vpp[6][i] = coreas_T4H[i] * scale_factor;
			x_vpp[1][i] = coreas_T1V[i] * scale_factor;
			x_vpp[3][i] = coreas_T2V[i] * scale_factor;
			x_vpp[5][i] = coreas_T3V[i] * scale_factor;
			x_vpp[7][i] = coreas_T4V[i] * scale_factor;
		}

		
		for(int i=0; i<1500; i++){
			t[i] = i*0.8;
			x_vpp[0][i] = (x_vpp[0][i])*500./32512;
			x_vpp[2][i] = (x_vpp[2][i])*500./32512;
			x_vpp[4][i] = (x_vpp[4][i])*500./32512;
			x_vpp[6][i] = (x_vpp[6][i])*500./32512;
			x_vpp[1][i] = (x_vpp[1][i])*500./32512;
			x_vpp[3][i] = (x_vpp[3][i])*500./32512;
			x_vpp[5][i] = (x_vpp[5][i])*500./32512;
			x_vpp[7][i] = (x_vpp[7][i])*500./32512;
		}
		
		TGraph *cyc = new TGraph(1500, t, x_vpp[1]);
		cyc->Draw();
		
		TGraph *wyn = new TGraph(1500, t, x_vpp[7]);
		wyn->SetLineColor(2);
		wyn->Draw("same");
		cout<<ShowerTheta<<"  "<<ShowerPhi<<endl;
	
}
