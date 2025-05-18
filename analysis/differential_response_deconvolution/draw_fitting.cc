

 void draw_fitting(){
 	
 	TH2F *h2_theta_R = new TH2F("h2_theta_R","h2_theta_R", 201, 1.5, 3.5, 119, 10000, 600000);
 	
 	TH2F *h2_theta_phi = new TH2F("h2_theta_phi","h2_theta_phi", 201, 1.5, 3.5, 41, -22.8, -22.4);
 	
 	double _theta, _phi, _R, _altitude, _chi_square, _delays[4], _direct_delta_t[6], _reflected_delta_t[6];
 	
	TFile * res_file = new TFile("_fitting_result_delay_direct_HV.root");
    TTree * fittingTree = (TTree*) res_file->Get("fitting");
    
    fittingTree->SetBranchAddress("theta", &_theta);
    fittingTree->SetBranchAddress("phi", &_phi);
    fittingTree->SetBranchAddress("R", &_R);
    fittingTree->SetBranchAddress("altitude", &_altitude);
    fittingTree->SetBranchAddress("chi_square", &_chi_square);
    
    double chi_square_theta_R[201][119] = {0.};
    double chi_square_theta_phi[201][41] = {0.};
    for(int i=0; i<201; i++){for(int j=0; j<119; j++){chi_square_theta_R[i][j] = 1.0e100;}}
    for(int i=0; i<201; i++){for(int j=0; j<41; j++){chi_square_theta_phi[i][j] = 1.0e100;}}
    
    for(int entry = 0; entry < fittingTree->GetEntries(); entry++){
   		fittingTree->GetEntry(entry);
   		
   		int int_theta = round((_theta-1.5)/0.01);
   		int int_phi = round((_phi-(-22.8))/0.01);
   		int int_R = round((_R-10000)/5000.);
   		if(_chi_square < chi_square_theta_R[int_theta][int_R]){
   			chi_square_theta_R[int_theta][int_R] = _chi_square;
   			if(_chi_square>20.03833) continue;
   			if(_theta<2&&_chi_square<7) cout<<_theta<<"  "<<_phi<<"   "<<_R<<"    "<<_altitude<<"   "<<_chi_square<<endl;
   			//7.03833
   			//12.84863
   			//20.0625
   			//if(_altitude>686||_altitude<684) continue;
   			//if(_phi<-22.555||_phi>-22.545) continue;
   			h2_theta_R->SetBinContent(int_theta, int_R, _chi_square);
   		}
   		
   		if(_chi_square < chi_square_theta_phi[int_theta][int_phi]){
   			chi_square_theta_phi[int_theta][int_phi] = _chi_square;
   			if(_chi_square>100) continue;
   			//if(_R<50) continue;
   			//7.03833
   			//12.84863
   			//20.0625
   			//if(_altitude>686||_altitude<684) continue;
   			//if(_phi<-22.555||_phi>-22.545) continue;
   			h2_theta_phi->SetBinContent(int_theta, int_phi, _chi_square);
   		}
   		
   		
   		

   		
   	}
   	
   	for(int i=0; i<201; i++){
   		for(int j=0; j<119; j++){
   			//if(chi_square_list[i][j]<20)cout<<chi_square_list[i][j]<<endl;
   		}
   	}
	

	res_file->Close();
    
	gStyle->SetPalette(55);
	//h2_theta_R->Draw("colorz");
	//h2_theta_R->GetZaxis()->SetRangeUser(4.5, 7.1);
	
	
	h2_theta_phi->Draw("colorz");
	//h2_theta_phi->GetZaxis()->SetRangeUser(4.5, 7.1);
	//double contours[3] = {7.03833, 12.84863, 20.0625};
	
   	//h2_theta_phi->SetContour(1, contours);
   	//h2_theta_phi->Draw("cont3 same");
   	//h2_theta_phi->SetLineColor(kRed);

}// end of main function
              

     


     	     						

