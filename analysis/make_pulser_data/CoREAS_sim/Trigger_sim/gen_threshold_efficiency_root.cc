
#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h> 
#include <iomanip>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TMultiGraph.h"

using namespace std;


double Sgn[18]={26,38,54,77,100,127,137,154,173,196,220,243,349,506,715,997,1418,1965};


void get_band_efficiency_vs_threshold(TGraph *gr[8][8][18]){
	for(int board=1; board<9; board++){
		for(int band=1; band<9; band++){
			for(int strength=0; strength<18; strength++){
  				char *txtname = Form("/home/cyc/Verilog/TAROGE3_FPGA/soft/V_vs_Sgn/V_pp_%dmV/Board%d_Band%d_efficiency_scan.dat", int(Sgn[strength]), board, band);//
  
       			FILE *fpin;
    			fpin = fopen(txtname,"r");
    			float a;
    			float b;
    			double voltage[20]={0};
    			double efficiency_count[20]={0};
    			double efficiency[20]={0};
    			if(!fpin) {cout<<txtname<<" open !ok"<<endl; return 0;}
 
        		int p = fscanf(fpin,"%f %f", &a, &b); 
    			int N =0;  
    			while(1){
        			voltage[N]=a; efficiency_count[N]=b;
        			p = fscanf(fpin,"%f %f", &a, &b);
        			if(p!=2) break;
        			N+=1;// because of break, N should +1;
            	}
    			fclose(fpin);
    
    			for(int i=N;i>0;i--){
       				if(efficiency_count[i]>efficiency_count[i-1]) efficiency_count[i-1] = efficiency_count[i];
                }
                        
    			for(int i=0;i<N+1;i++){
       				efficiency[i] = efficiency_count[i]/16543.;
       				if(efficiency[i]>1) efficiency[i] = 1;
                }
       
    			gr[board-1][band-1][strength] = new TGraph(N+1,voltage,efficiency);
    		}
    	}
    }
}

void get_band_threshold_vs_step(TGraph *gr[8][8]){
	for(int board=1; board<9; board++){
		for(int band=1; band<9; band++){
  				char *txtname = Form("/home/cyc/Verilog/TAROGE3_FPGA/soft/V_vs_Sgn/Board%d_Band%d_efficiency_scan.dat", board, band);//
  
       			FILE *fpin;
    			fpin = fopen(txtname,"r");
    			float a; //a is not important
    			float b;
    			double step[10010]={0};
    			double voltage[10010]={0};
    			if(!fpin) {cout<<txtname<<" open !ok"<<endl; return 0;}
 
        		int p = fscanf(fpin,"%f %f", &a, &b); 
    			int N =0;  
    			while(1){
        			step[N]=N; voltage[N]=b;
        			p = fscanf(fpin,"%f %f", &a, &b);
        			if(p!=2) break;
        			N+=1;// because of break, N should +1;
            	}
    			fclose(fpin);
    
        
    
    			gr[board-1][band-1] = new TGraph(N+1, step, voltage);
    		
    	}
    }
}

void gen_threshold_efficiency_root(){

	TFile *efficiency_vs_strength = new TFile("efficiency_vs_strength.root","RECREATE");
	
	TGraph * gr_eff_vs_strength[1000][8][8] = {NULL};
	
	TGraph * gr_eff_vs_thred[8][8][18] = {NULL};
	get_band_efficiency_vs_threshold(gr_eff_vs_thred);
	
	TGraph * gr_thred_vs_step[8][8] = {NULL};
	get_band_threshold_vs_step(gr_thred_vs_step);

	double efficiency_list[18]={0};
	double strength_list[18]={0};
	
	for(int step=300; step<1300; step++){
		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){
			
				double threshold = gr_thred_vs_step[board-1][band-1]->Eval(step);
				//cout<<endl<<endl<<board<<"  "<<band<<"  "<<threshold<<endl;
				
				for(int strength=0; strength<18; strength++){
    				double efficiency = gr_eff_vs_thred[board-1][band-1][strength]->Eval(threshold);
    				if(efficiency>1.){efficiency=1.;}
    				if(efficiency<0.){efficiency=0.;}
    				strength_list[strength] = Sgn[strength];
    				efficiency_list[strength] = efficiency;
    				//cout<<efficiency<<"   ";
    			}
    			gr_eff_vs_strength[step-300][board-1][band-1] = new TGraph(18, strength_list, efficiency_list);
    			//delete gr_eff_vs_strength[step-400][board-1][band-1];
    			//TGraph *cyc = new TGraph(18, strength_list, efficiency_list);
    			//delete cyc;
    			char *txtname = Form("gr_Step%d_Board%d_Band%d", step, board, band);
    			gr_eff_vs_strength[step-300][board-1][band-1]->SetName(txtname);
				gr_eff_vs_strength[step-300][board-1][band-1]->SetTitle("Band efficiency vs V_pp signal strength for given step threshold");
    			gr_eff_vs_strength[step-300][board-1][band-1]->Write();
    		}
    	}
    }

    efficiency_vs_strength->Close();
    delete efficiency_vs_strength;

}


