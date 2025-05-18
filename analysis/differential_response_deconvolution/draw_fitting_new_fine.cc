
template <typename T>
T** new_Array2D(int row, int col);

template <typename T> 
void delete_Array2D(T **arr, int row, int col);

 void draw_fitting_new_fine(){
 	

 	TH2F *h2_theta_R = new TH2F("h2_theta_R","h2_theta_R", 31, 2.15, 2.75, 3501, 50000, 400000);

 	TH2F *h2_theta_phi = new TH2F("h2_theta_phi","h2_theta_phi", 31, 2.15, 2.75, 41, -22.75, -22.35);

 	TH2F *h2_altitude_R = new TH2F("h2_altitude_R","h2_altitude_R", 401, 678.2, 698.2, 3501, 50000, 400000);

 	double _theta, _phi, _R, _altitude, _chi_square, _delays[4], _direct_delta_t[6], _reflected_delta_t[6];

 	
 	
 	
	TFile * res_file = new TFile("/media/cyc/SSD_data_disk/_fitting_result_fine_delay_direct_HV_reflected.root");
    TTree * fittingTree = (TTree*) res_file->Get("fitting");

    
    
    fittingTree->SetBranchAddress("theta", &_theta);
    fittingTree->SetBranchAddress("phi", &_phi);
    fittingTree->SetBranchAddress("R", &_R);
    fittingTree->SetBranchAddress("altitude", &_altitude);
    fittingTree->SetBranchAddress("chi_square", &_chi_square);
    
    float chi_square_theta_R[31][3501] = {0.};
    float chi_square_theta_phi[31][41] = {0.};
    //double chi_square_altitude_R[401][3501] = {0.};
    float **chi_square_altitude_R = new_Array2D<float>(401, 3501);

    for(int i=0; i<31; i++){for(int j=0; j<3501; j++){chi_square_theta_R[i][j] = 1.0e100;}}
    for(int i=0; i<31; i++){for(int j=0; j<41; j++){chi_square_theta_phi[i][j] = 1.0e100;}}
    for(int i=0; i<401; i++){for(int j=0; j<3501; j++){chi_square_altitude_R[i][j] = 1.0e100;}}
    
    for(int entry = 0; entry < fittingTree->GetEntries(); entry++){
   		fittingTree->GetEntry(entry);
   		
   		int int_theta = round((_theta-2.15)/0.02);
   		int int_phi = round((_phi-(-22.75))/0.01);
   		int int_R = round((_R-50000)/100.);
   		int int_altitude = round((_altitude-678.2)/0.05);
   		if(_chi_square < chi_square_theta_R[int_theta][int_R]){
   			
   			//if(_chi_square>100) continue;
   			if(_R>350000) continue;
   			//if(_theta<2&&_chi_square<7) cout<<_theta<<"  "<<_phi<<"   "<<_R<<"    "<<_altitude<<"   "<<_chi_square<<endl;
   			//7.03833
   			//12.84863
   			//20.0625
   			//if(_altitude>689||_altitude<688) continue;
   			//if(_phi<-22.555||_phi>-22.545) continue;
   			//if(_theta>1.5) continue;
   			chi_square_theta_R[int_theta][int_R] = _chi_square;
   			h2_theta_R->SetBinContent(int_theta, int_R, _chi_square);
   		}
   		
   		if(_chi_square < chi_square_theta_phi[int_theta][int_phi]){
   			chi_square_theta_phi[int_theta][int_phi] = _chi_square;
   			//if(_chi_square>100) continue;
   			if(_R>350000) continue;
   			//7.03833
   			//12.84863
   			//20.0625
   			//if(_altitude>689||_altitude<688) continue;
   			//if(_phi<-22.555||_phi>-22.545) continue;
   			//if(_R<50000) continue;
   			chi_square_theta_phi[int_theta][int_phi] = _chi_square;
   			h2_theta_phi->SetBinContent(int_theta, int_phi, _chi_square);
   		}
   		
   		
   		if(_chi_square < chi_square_altitude_R[int_altitude][int_R]){
   			
   			//if(_chi_square>100) continue;
   			if(_R>350000) continue;
   			//7.03833
   			//12.84863
   			//20.0625
   			//if(_altitude>689||_altitude<688) continue;
   			//if(_phi<-22.555||_phi>-22.545) continue;
   			//if(_R<50000) continue;
   			chi_square_altitude_R[int_altitude][int_R] = _chi_square;
   			h2_altitude_R->SetBinContent(int_altitude, int_R, _chi_square);
   		}

   		

   	}
   	
   	for(int i=0; i<301; i++){
   		for(int j=0; j<171; j++){
   			//if(chi_square_list[i][j]<20)cout<<chi_square_list[i][j]<<endl;
   		}
   	}
	

	res_file->Close();
    
    
    
       		
   		
   		// delay direct_HV reflected_H   20.27979  29.24219  39.17188  ==>> normalized to 18: 21.471340, 30.960331, 41.473446
   		// direct_HV reflected_H   15.93579  24.02539  33.1875
   		// delay direct_HV   13.74463  21.34863  30.09375
   		// direct_HV   9.30396   15.78906  23.57812
    
    
	gStyle->SetPalette(55);
	

	Double_t levels[] = {0, 21.471340, 30.960331, 41.473446, 100, 1.0e3, 1.0e4, 1.0e5, 1.0e12};
	h2_theta_R->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
	
	TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
	h2_theta_R->Draw("colorz");
	//h2_theta_R->GetZaxis()->SetRangeUser(4.5, 7.1);
	
	
	h2_theta_phi->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
	TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
	h2_theta_phi->Draw("colorz");
	//h2_theta_phi->GetZaxis()->SetRangeUser(4.5, 7.1);
	//double contours[3] = {21.471340, 30.960331, 41.473446};
	
   	//h2_theta_phi->SetContour(3, contours);
   	//h2_theta_phi->Draw("cont3 same");
   	h2_theta_phi->SetLineColor(kRed);
   	
   	
   	h2_altitude_R->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
   	TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
   	h2_altitude_R->Draw("colorz");
	
}// end of main function
              

     

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
     	     						

