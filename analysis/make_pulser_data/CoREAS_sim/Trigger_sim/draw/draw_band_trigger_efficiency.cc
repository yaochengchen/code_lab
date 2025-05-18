

void draw_band_trigger_efficiency(){

	TFile* f0 = new TFile("../efficiency_vs_strength.root");
  
  	TGraph * gr_eff_vs_strength[1000][8][8] = {NULL};
	
	for(int step=500; step<901; step++){
		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){
				char *txtname = Form("gr_Step%d_Board%d_Band%d", step, board, band); 
				gr_eff_vs_strength[step-300][board-1][band-1] = (TGraph*)f0->Get(txtname); 
			}
		}
	}


double Sgn[18]={26,38,54,77,100,127,137,154,173,196,220,243,349,506,715,997,1418,1965};

		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){ 
				
				
				//gr_eff_vs_strength[200][board-1][band-1]->Draw("apl");
				
				int n = gr_eff_vs_strength[200][board-1][band-1]->GetN();
				double *x = gr_eff_vs_strength[200][board-1][band-1]->GetX();
				double *y = gr_eff_vs_strength[200][board-1][band-1]->GetY();
				
				for(int i=0; i<n; i++){
					x[i] /= 2.;// from V_pp to V_p
				}
				
				for(int i=11; i<n; i++){
					if(y[i]<0.995){y[i]=0.999;}
				}
				
				for(int i=0; i<3; i++){
					if(y[i]>0.1){y[i]=0.0001;}
				}
				gr_eff_vs_strength[200][board-1][band-1] = new TGraph(n, x, y);
				
			}
		}
		

		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){ 
				
				if(board==1&&band==1){
					gr_eff_vs_strength[200][board-1][band-1]->Draw("al");
					gr_eff_vs_strength[200][board-1][band-1]->GetXaxis()->SetTitle("V_{peak} (mV)");
					gr_eff_vs_strength[200][board-1][band-1]->GetYaxis()->SetTitle("Efficiency");
					gr_eff_vs_strength[200][board-1][band-1]->GetXaxis()->SetRangeUser(0,500);
					gr_eff_vs_strength[200][board-1][band-1]->GetYaxis()->SetRangeUser(0,1);
				}
				else{
					gr_eff_vs_strength[200][board-1][band-1]->Draw("l same");
				}
				gr_eff_vs_strength[200][board-1][band-1]->SetLineColor(1);
			}
		}




		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){ 
				
				
				//gr_eff_vs_strength[200][board-1][band-1]->Draw("apl");
				
				int n = gr_eff_vs_strength[400][board-1][band-1]->GetN();
				double *x = gr_eff_vs_strength[400][board-1][band-1]->GetX();
				double *y = gr_eff_vs_strength[400][board-1][band-1]->GetY();
				
				for(int i=0; i<n; i++){
					x[i] /= 2.;// from V_pp to V_p
				}
				
				for(int i=13; i<n; i++){
					if(y[i]<0.995){y[i]=0.999;}
				}
				
				for(int i=0; i<3; i++){
					if(y[i]>0.001){y[i]=0.0001;}
				}
				gr_eff_vs_strength[400][board-1][band-1] = new TGraph(n, x, y);
				
			}
		}
		

		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){ 
				
				
				gr_eff_vs_strength[400][board-1][band-1]->Draw("l same");
				
				gr_eff_vs_strength[400][board-1][band-1]->SetLineColor(2);
			}
		}


		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){ 
				
				
				//gr_eff_vs_strength[200][board-1][band-1]->Draw("apl");
				
				int n = gr_eff_vs_strength[600][board-1][band-1]->GetN();
				double *x = gr_eff_vs_strength[600][board-1][band-1]->GetX();
				double *y = gr_eff_vs_strength[600][board-1][band-1]->GetY();
				
				for(int i=0; i<n; i++){
					x[i] /= 2.;// from V_pp to V_p
				}
				
				for(int i=13; i<n; i++){
					if(y[i]<0.995){y[i]=0.999;}
				}
				
				for(int i=0; i<3; i++){
					if(y[i]>0.001){y[i]=0.0001;}
				}
				gr_eff_vs_strength[600][board-1][band-1] = new TGraph(n, x, y);
				
			}
		}
		

		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){ 
				
				
				gr_eff_vs_strength[600][board-1][band-1]->Draw("l same");
				
				gr_eff_vs_strength[600][board-1][band-1]->SetLineColor(3);
			}
		}
		
		TLegend *leg_R = new TLegend(0.65,0.6,0.9,0.9);
		leg_R->AddEntry(gr_eff_vs_strength[200][0][0],"Threshold = 500","l");
		leg_R->AddEntry(gr_eff_vs_strength[400][0][0],"Threshold = 700","l");
		leg_R->AddEntry(gr_eff_vs_strength[600][0][0],"Threshold = 900","l");

		leg_R->Draw();
		
		
}
