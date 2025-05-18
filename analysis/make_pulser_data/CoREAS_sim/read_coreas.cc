#include <iostream>
#include <vector>
#include <fstream>
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TApplication.h"

#include <unistd.h>
#include <fcntl.h>


#include "constant_value.h"

using namespace std;


int readLongProfile(const char* fname_long, const char* particleName, vector<double> &vX, vector<double> &vN );
bool plotLongProfile(const char* fname_long);

bool saveCoreasToROOT(const char* stnName, const char* dir_cr, int run, const char* dir_out, bool bOverwrite);
bool loopSaveROOT(const char* stnName, const char* dir_cr, int run1, int run2, const char* dir_out, bool bOverwrite);

bool decodeListLine(string str, TVector3& vec_pos, string& str_pos, double& radius, double& phi, double& min_gamma, double &max_gamma);
double getTokenInt(string str, int nth);
double getTokenDouble(string str, int nth);

Bool_t FileExist(const string& fname);



int main(int argc, char *argv[])
{ 
	TApplication *app = new TApplication("App", &argc, argv);
	//parse string  

	// The recognized options are removed from the argument array. The original list of argument options can be retrieved via the Argc() and Argv() method
	printf("argc: %d\n", app->Argc() );
	for(int c=0; c<app->Argc(); c++) printf("%d %s\n", c, app->Argv()[c]);
	if(app->Argc() < 7)
	{
		//help();
		printf("argc = %d < 7   loopSaveROOT(const char* stnName, const char* dir_cr, int run1, int run2, const char* dir_out, bool bOverwrite)\n", argc);
		
		return 1;
	}

	printf("run planefit3d(%s, %s, %s, %s, %s, %s)\n", app->Argv()[1], app->Argv()[2], app->Argv()[3], app->Argv()[4], app->Argv()[5], app->Argv()[6]);
	
	//bool loopSaveROOT(const char* stnName, const char* dir_cr, int run1, int run2, const char* dir_out, bool bOverwrite)
	loopSaveROOT( 
		app->Argv()[1], app->Argv()[2], 
		strtol( app->Argv()[3], NULL, 10), 
		strtol( app->Argv()[4], NULL, 10),
		app->Argv()[5], (bool) strtol( app->Argv()[6], NULL, 10 )  ); //, bSim); //
	//If you manage the “infinite loop” yourself, don’t call “TApplication::Run” at all. Inside of your “infinite loop” you may call “gSystem->ProcessEvents();” in regular intervals.
	//app->Run();

	delete app;

	return 0;
}


//bool readCorsikaInput(const char* fname){

//	return true;
//}

//loop over all [run number] sub-folders under dir_cr
//loopSaveROOT("t4", "/media/cyc/For_Linux/CoREAS_Sim/proton", 31009, 31009, "/media/cyc/For_Linux/CoREAS_Sim/proton/cyc", 1)
bool loopSaveROOT(const char* stnName, const char* dir_cr, int run1, int run2, const char* dir_out, bool bOverwrite)
{

	for(int r=run1; r<=run2; r++)
	{
		//if(r%1000>19){continue;}
		//cout<<"run "<< r<<endl;
		saveCoreasToROOT(stnName, dir_cr, r, dir_out, bOverwrite);
	}
	return true;
}

/*
	convert CoREAS output to ROOT file

	bOverwrite: =false: skip if file already exists
*/
bool saveCoreasToROOT(const char* stnName, const char* dir_cr_, int run, const char* dir_out, bool bOverwrite)
{


	const char* dir_cr = Form("%s/%d/", dir_cr_, run); //old ver: %d


	const string strCoreasDir = Form("SIM%06d_coreas/", 1);
	
	
	const string wf_files = dir_cr + strCoreasDir;
	// 0 means exist, -1 means not
	if(access(wf_files.c_str(), F_OK)!=0)
	{
		cout<<wf_files<<"  do not exist! skip it"<<endl;
		return 0;
	}
		
		//"SIM000001_coreas/"; //where E-field are stored
	const string strLongFile = Form("DAT%06d.long", 1);
	//"DAT000001.long"; 	//output longitudinal profile
	const string strListFile = Form("SIM%06d.list", 1);
	//"SIM000001.list"; //CoREAS input: position where E-field is sampled
	const string strReasFile = Form("SIM%06d.reas", 1);
	//"SIM000001.reas"; //CoREAS input & output

	TTree *tr_cr=NULL, *tr_wf=NULL, *tr_ = NULL, *tr_long=NULL, *tr_parno=NULL, *tr_edepo=NULL, *tr_ghfit=NULL;

	bool bWaveform = false;
	int positionID=0, particleID=0, NXN=0, NXE=0;
	double radius, phi, min_gamma=0., max_gamma=0., depth, depth_prev;
	double energy_eV, zenith_deg, cherenkov_angle_deg, azimuth_corsika_deg, Xmax, alpha_deg, distance_Xmax_m,BfieldMag_uT, BfieldIncline_deg, samplingPeriod_s; 
	double peak_amplitude_north, peak_amplitude_west, peak_amplitude_vert;
	vector<double> vTime, vE_north, vE_west, vE_vert;
	string str, str_pos, strWfFile;
	TVector3 vec_pos, vec_core;

	const char* fname_reas = Form("%s%s", dir_cr, strReasFile.c_str() );
	const char* fname_list = Form("%s%s", dir_cr, strListFile.c_str() );
	const char* fname_long = Form("%s%s", dir_cr, strLongFile.c_str() );


	const char* fname_out = Form("%s/Coreas-%s-r%06d.root", dir_out, "t4" , run);

	if( !bOverwrite && FileExist(fname_out ) )
	{
		printf("output file already exists: %s  skip\n", fname_out );
		return false;
	}

	if(!FileExist(fname_reas) || !FileExist(fname_list) )
	{
		printf("no input: %s or %s\n", fname_reas, fname_list);
		return false;
	}
	//read text 
	ifstream 
		inFile_list( fname_list ), 
		inFile_wf, 
		inFile_reas(fname_reas), 
		inFile_long(fname_long);

	//if(!inFile_list.good() || !inFile_reas.good() )
	//{
	//	printf("no input: %s or %s\n", fname_reas, fname_list);
	//	return false;
	//}


	TFile *f_out = new TFile(fname_out, "recreate");

	tr_cr = new TTree("CRTree","from .reas file");
	tr_cr->Branch("run", &run, "run/I");
	tr_cr->Branch("particleID", &particleID, "particleID/I" ); //14=proton
    tr_cr->Branch("energy_eV", &energy_eV, "energy_eV/D"); //eV
    tr_cr->Branch("zenith_deg", &zenith_deg, "zenith_deg/D");  //deg

    tr_cr->Branch("azimuth_corsika_deg", &azimuth_corsika_deg, "azimuth_corsika_deg/D"); //deg //corsika geomanetic north, 0: shower propagates to north, 90: to west


    tr_cr->Branch("Xmax", &Xmax, "Xmax/D" ); //g/cm2 slant depth
    tr_cr->Branch("alpha_deg", &alpha_deg, "alpha_deg/D"); //GeomagneticAngle: deg  pxB
    tr_cr->Branch("distance_Xmax_m", &distance_Xmax_m, "distance_Xmax_m/D" );   //geometrical distance of shower maximum from core in cm  -> m
    tr_cr->Branch("samplingPeriod_s", &samplingPeriod_s, "samplingPeriod_s/D");
    tr_cr->Branch("vecCore", &vec_core ); //cm to m
    tr_cr->Branch("BfieldMag_uT", &BfieldMag_uT, "BfieldMag_uT/D"); //0.6227455483
    tr_cr->Branch("BfieldIncline_deg", &BfieldIncline_deg, "BfieldIncline_deg/D"); //in degrees, >0: in northern hemisphere, <0: in southern hemisphere

   //decodeReasFile( fname_reas, tr_cr);

    
    //CR parameters should be consistent between .reas & .inp
    while( getline( inFile_reas, str ) ) 
    {
		//cm to m
		if(str.find("CoreCoordinateNorth") != std::string::npos) vec_core.SetX( 0.01 * getTokenDouble(str, 2));
		else if(str.find("CoreCoordinateWest") != std::string::npos) vec_core.SetY( 0.01 * getTokenDouble(str, 2)); //geomanetic basis (0,0,h)
		else if(str.find("CoreCoordinateVertical") != std::string::npos) vec_core.SetZ( 0.01 * getTokenDouble(str, 2)); 
		else if(str.find("TimeResolution") != std::string::npos) samplingPeriod_s = getTokenDouble(str, 2);
		else if(str.find("ShowerZenithAngle") != std::string::npos) zenith_deg = getTokenDouble(str, 2);
		else if(str.find("ShowerAzimuthAngle") != std::string::npos) azimuth_corsika_deg = getTokenDouble(str, 2); //momentum to arriving direction
		else if(str.find("PrimaryParticleEnergy") != std::string::npos) energy_eV = getTokenDouble(str, 2);
		else if(str.find("PrimaryParticleType") != std::string::npos) particleID = getTokenInt(str, 2);
		else if(str.find("DepthOfShowerMaximum") != std::string::npos) Xmax = getTokenDouble(str, 2); //g/cm2
		else if(str.find("DistanceOfShowerMaximum") != std::string::npos) distance_Xmax_m = 0.01* getTokenDouble(str, 2); //cm to m
		else if(str.find("MagneticFieldStrength") != std::string::npos) BfieldMag_uT = getTokenDouble(str, 2) * 100; //gauss to micro Tesla
		else if(str.find("MagneticFieldInclinationAngle") != std::string::npos) BfieldIncline_deg = getTokenDouble(str, 2);
		else if(str.find("GeomagneticAngle") != std::string::npos) alpha_deg = getTokenDouble(str, 2);
		else{
			//printf("skip %s\n", str.c_str());
		}
	}

	inFile_reas.close();




	tr_cr->Fill();
	tr_cr->Write(); 

    //read .long
	tr_long = new TTree("tr_long","longitudinal profile"); //temporary tree which will not be saved
	tr_long->ReadFile( fname_long , "depth/D:gamma:positron:electron:mu_bar:mu:hadron:charged:nuclei:cherenkov", ' ');
	tr_long->SetBranchAddress("depth", &depth);

	cout<<"Number of entries: "<< tr_long->GetEntries() <<endl;



	//first half: longitudinal profile
	//second half: deposited energy

	//separate 
	tr_long->GetEntry(0);
	depth_prev = depth;

	for(int n=1; n<tr_long->GetEntries(); n++)
	{
		tr_long->GetEntry(n);
		if( depth < depth_prev ){ //X should increase monotonically 
			NXN = n; //particle number
			NXE = tr_long->GetEntries() - NXN; //energy posit
			printf("break at entry %d: %.1f < %.1f:\t NXN = %d \t NXE = %d\n", n, depth, depth_prev, NXN, NXE);

			break;
		}else	depth_prev = depth;
	}

	//each entry contains one depth
	//##### alternative  one entry contains vector of all depths
	tr_parno = tr_long->CopyTree( Form("Entry$<%d", NXN) );
	tr_edepo = tr_long->CopyTree( Form("Entry$>=%d", NXN)  );

	tr_parno->SetName("LongProfTree");
	tr_edepo->SetName("EDepositTree");

	f_out->cd();
	tr_parno->Write();
	tr_edepo->Write();


	double P1, P2, P3, P4, P5, P6, chi2dof, avgDevPercent;
	tr_ghfit = new TTree("FitTree", "Gaisser-Hillas curve on all charged particles");
	tr_ghfit->Branch("P1", &P1, "P1/D");
	tr_ghfit->Branch("P2", &P2, "P2/D");
	tr_ghfit->Branch("P3", &P3, "P3/D");
	tr_ghfit->Branch("P4", &P4, "P4/D");
	tr_ghfit->Branch("P5", &P5, "P5/D");
	tr_ghfit->Branch("P6", &P6, "P6/D");
	tr_ghfit->Branch("chi2dof", &chi2dof, "chi2dof/D");
	tr_ghfit->Branch("avgDevPercent", &avgDevPercent, "avgDevPercent/D"); 

    while( getline( inFile_long, str ) ) {

		if(str.find("PARAMETERS") != std::string::npos){
			P1 =  getTokenDouble(str, 2);
			P2 =  getTokenDouble(str, 3);
			P3 =  getTokenDouble(str, 4);
			P4 =  getTokenDouble(str, 5);
			P5 =  getTokenDouble(str, 6);
			P6 =  getTokenDouble(str, 7);
		} 		else if(str.find("CHI**2/DOF") != std::string::npos) chi2dof = getTokenDouble(str, 2);
		else if(str.find("AV. DEVIATION IN %") != std::string::npos ) avgDevPercent = getTokenDouble(str, 5);
		else{
			//printf("skip %s\n", str.c_str());
		}
	}

	inFile_reas.close();

	f_out->cd();
	tr_ghfit->Fill();
	tr_ghfit->Show(0);
	tr_ghfit->Write();



	//waveforms
	tr_wf = new TTree("EfieldTree","");
	tr_wf->Branch("bWaveform", &bWaveform, "bWaveform/O" );	//if waveform data is available
	tr_wf->Branch("positionID", &positionID, "positionID/I");
	tr_wf->Branch("radius", &radius, "radius/D");
	tr_wf->Branch("phi", &phi, "phi/D");
	tr_wf->Branch("index", &str_pos[0], "index/C");
	tr_wf->Branch("vecPos", &vec_pos ); //TVector3 in m

	tr_wf->Branch("minGamma", &min_gamma, "minGamma/D");
	tr_wf->Branch("maxGamma", &max_gamma, "maxGamma/D");

	//t (s)  North, west , vertical, cgs unit to SI
	tr_wf->Branch("vTime", &vTime );
	tr_wf->Branch("vE_north", &vE_north );
	tr_wf->Branch("vE_west", &vE_west );
	tr_wf->Branch("vE_vert", &vE_vert );


	//fixed or variable length?



	tr_wf->Branch("peak_amplitude_north", &peak_amplitude_north, "peak_amplitude_north/D"); //E-field
	tr_wf->Branch("peak_amplitude_west", &peak_amplitude_west, "peak_amplitude_west/D"); //E-field
	tr_wf->Branch("peak_amplitude_vert", &peak_amplitude_vert, "peak_amplitude_vert/D"); //E-field


	tr_ = new TTree("tr_","temp");

	//read .reas file for config

	//read .list for locations
		//read list of files

	//CORSIKA: The zenith angle θ of a particle trajectory is measured between the particle momentum vector and the negativez-axis, and the azimuthal angleφbetween the positivex-axis and thex-y-component of the particle momentum vector (i.e. with respect tonorth) proceeding counterclockwise


	printf("zenith: %.3f deg  azimuth: %.3f deg (cor)\n", zenith_deg, azimuth_corsika_deg);
 


	while( getline( inFile_list, str ) ) {
		
		//vectors in Corsika frame! (gN, gW, z)

		//bool decodeListLine(string str, TVector3& vec_pos, string& str_pos, double& radius, double& phi, double& min_gamma, double &max_gamma) 
		decodeListLine( str, vec_pos, str_pos, radius, phi, min_gamma, max_gamma); //converted from cm to m
		
		if((run>=46001)&&(run<=46019)){if(radius>=1720){continue;}}
		if((run>=47001)&&(run<=47019)){if(radius>=2040){continue;}}
		if((run>=66039)&&(run<=66057)){if(radius>=2240){continue;}}
		if((run>=67001)&&(run<=67019)){if(radius>=2040){continue;}}
		if((run>=57001)&&(run<=57019)){if(radius>=2040){continue;}}
		if(run==57020){continue;}
		if((run>=76039)&&(run<=76057)){if(radius>=2040){continue;}}
		if((run>=77001)&&(run<=77019)){if(radius>=2040){continue;}}
		
		

		peak_amplitude_north=0.;
		peak_amplitude_west=0.;
		peak_amplitude_vert=0.;

		//get waveform accordingly
		bWaveform = false;

		strWfFile = dir_cr + strCoreasDir + "raw_" + str_pos +".dat";
		inFile_wf.open( strWfFile.c_str() );
		if(!inFile_wf.good() ){
			printf("(saveCoreasToROOT) Error: waveform file not found: %s\n", str_pos.c_str() );

			bWaveform = false;

			vTime.assign(  vTime.size(), 0. );
			vE_north.assign( vE_north.size(), 0. );
			vE_west.assign( vE_west.size(), 0. );
			vE_vert.assign( vE_vert.size(), 0. );

		}else{

			//CoREAS: /SIM[]_coreas/*.dat:  columns in the ﬁle denote the absolute time stamp and the north-, west-, and vertical component of the electric ﬁeld in cgs unit
			tr_->ReadFile( strWfFile.c_str(), "times/D:E_north:E_west:E_vert"); 
			tr_->Draw("times:E_north:E_west:E_vert","","goff"); //

			if(tr_->GetEntries()>0 ) bWaveform = true;


			//vector::assign(): Any elements held in the container before the call are destroyed and replaced by newly constructed elements (no assignments of elements take place). This causes an automatic reallocation of the allocated storage space if -and only if- the new vector size surpasses the current vector capacity.

			vTime.assign( tr_->GetV1(), tr_->GetV1() + tr_->GetSelectedRows() ); //begin, end
			vE_north.assign( tr_->GetV2(), tr_->GetV2() + tr_->GetSelectedRows() );
			vE_west.assign( tr_->GetV3(), tr_->GetV3() + tr_->GetSelectedRows() );
			vE_vert.assign( tr_->GetV4(), tr_->GetV4() + tr_->GetSelectedRows() );

			tr_->Reset();

			//cgs to SI units and find peak E-field
			for(std::vector<double>::iterator it = vE_north.begin(); it != vE_north.end(); ++it)
			{
				(*it) *= E_conversion_factor;
				if( *it > peak_amplitude_north ) peak_amplitude_north = *it;
			}
			for(std::vector<double>::iterator it = vE_west.begin(); it != vE_west.end(); ++it)
			{
				(*it) *= E_conversion_factor;
				if( *it > peak_amplitude_west ) peak_amplitude_west = *it;
			}
			for(std::vector<double>::iterator it = vE_vert.begin(); it != vE_vert.end(); ++it)
			{
				(*it) *= E_conversion_factor;
				if( *it > peak_amplitude_vert ) peak_amplitude_vert = *it;
			}


			inFile_wf.close();


		}

		tr_wf->Fill();
		positionID++;
		
	}

	inFile_list.close();

	cout<<"number of entries: "<< tr_wf->GetEntries() <<endl; 

	tr_wf->Write();

	delete tr_;
	delete tr_long;

	f_out->Close();

	printf("saveCoreasToROOT for run [%d] done\n", run);


	return true;
}


double getTokenInt(string str, int nth){

  	char *token = strtok( &str[0], " "); 		  
    // Keep printing tokens while one of the 
    // delimiters present in str[]. 
    int n=0, num=0;
   

    while (token != NULL){ 

    	if(n==nth)  num = atoi(token);
		// a null pointer may be specified, in which case the function continues scanning where a previous successful call to the function ended.
        token = strtok(NULL, " "); 

        n++;
    }

	return num;
}

double getTokenDouble(string str, int nth){


  	char *token = strtok( &str[0], " "); 
		  
    // Keep printing tokens while one of the 
    // delimiters present in str[]. 
    int n=0;
    double  num=0.;

    while (token != NULL){ 

    	if(n==nth)  num = atof(token);
		// a null pointer may be specified, in which case the function continues scanning where a previous successful call to the function ended.
        token = strtok(NULL, " "); 

        n++;
    }

	return num;
}

/*
	CoREAS manual:
	 The columns signify the position to north, the position to west and the height asl, all in cm, followed by a unique
	name for the observer. 

	Within the CoREAS code, the same coordinate conventions as in CORSIKA are being used. For spatial coordinates, this means that a right-handed coordinate system of x, y and z is used where x denotes the geomagnetic north direction,
	y denotes the west direction and z denotes the vertical direction.

	
	example
	//AntennaPosition = 1370.02280194 -2111.95284188 3000.0 pos_25_0 gamma 2.0 1.0e35
*/
bool decodeListLine(string str, TVector3& vec_pos, string& str_pos, double& radius, double& phi, double& min_gamma, double &max_gamma){



  	char *token = strtok( &str[0], " "); 
		  
    // Keep printing tokens while one of the 
    // delimiters present in str[]. 
    int n=0;
    double  num;

    while (token != NULL){ 
    
    	//if(n>=6)printf("%s\n", token);

    	if(n>=2 && n<=4) vec_pos[n-2] = 0.01* atof( token ) ; //##### cm to m
    	else if(n==5){
    		//break into indices
    		int m=0;
    		str_pos = token;

    		size_t pos1 = str_pos.find("_");
    		size_t pos2 = str_pos.find("_", pos1+1);

    		//pos, len
    	//	printf("%s  %s\n", (str_pos.substr(pos1+1, pos2-pos1-1) ).c_str(), (str_pos.substr(pos2+1 ) ).c_str()  );
    		radius = atof( (str_pos.substr(pos1+1, pos2-pos1-1) ).c_str() );
    		phi = atof( (str_pos.substr(pos2+1 ) ).c_str() ); //to the end of string

    		//not worked for nested loop as strtok() use same memory location
    		//	toktok = strtok(NULL, " "); 


    	}else if(n==7) min_gamma =  atof( token ) ;
    	else if(n==8) max_gamma =  atof( token ) ;
    	else if(n>8){

    		printf("(decodeListLine) Warning: number of tokens more than expected %d\n", n);
    		return false;
    	}
       
        // a null pointer may be specified, in which case the function continues scanning where a previous successful call to the function ended.
        token = strtok(NULL, " "); 

       

        n++;
    } 



	return true;
}






/*
read .long file
longitudinal profile & energy deposit


*/
bool plotLongProfile(const char* fname_long){

	const int NParType = 5;
	const char* parName[NParType] = { "#gamma", "e^{-}", "e^{+}", "#mu^{+}+#mu^{-}", "charged"};
	const Color_t color[NParType] = { kBlack, kRed, kBlue, kGreen+2, kOrange+7};

	int NXN=0, NXE=0; 

	TGraph* gr_N[NParType] = {NULL}, *gr_E[NParType]={NULL};
	TMultiGraph *mgr = NULL;
	TCanvas* c = NULL;
	TLegend *lg = NULL;




	mgr = new TMultiGraph("mgr","proton shower, E= 2.54#times10^{17} eV, #theta=70#circ,;slant depth X (g/cm^{2});number N");
	lg = new TLegend(0.6, 0.6, 0.8, 0.85);

	TTree* tr = new TTree("tr","longitudinal profile");
	tr->ReadFile( fname_long, "depth/D:gamma:positron:electron:mu_bar:mu:hadron:charged:nuclei:cherenkov", ' ');

	cout<<"Number of entries: "<< tr->GetEntries() <<endl;

	//first half: longitudinal profile
	//second half: deposited energy
	tr->Draw("depth:gamma:positron:electron:mu+mu_bar:charged", "", "para goff");



	//separate  
	for(int n=1; n<tr->GetSelectedRows(); n++){
		if( tr->GetVal(0)[n] < tr->GetVal(0)[n-1] ){ //X should increase monotonically 
			NXN = n; //particle number
			NXE = tr->GetSelectedRows() - NXN; //energy posit
			printf("break at %d: %.1f < %.1f:\t NXN = %d \t NXE = %d\n", n, tr->GetVal(0)[n], tr->GetVal(0)[n-1], NXN, NXE);

			break;
		}
	}


	//longitudinal profile only
	for(int i=0;i<NParType;i++){

		//gr_N[i] = new TGraph( tr->GetSelectedRows(), tr->GetVal(0), tr->GetVal( i+1 ) );
		gr_N[i]  = new TGraph( NXN, tr->GetVal(0), tr->GetVal( i+1 ) );
		gr_N[i]->SetLineColor( color[i] );	gr_N[i]->SetLineWidth(2);
		mgr->Add( gr_N[i], "l" );
		lg->AddEntry( gr_N[i], parName[i], "l");
	}

	//check   
	for(int n=0; n<gr_N[0]->GetN(); n++){
		printf("%d\t%.2f\t%.3f\n", n, gr_N[0]->GetX()[n], gr_N[0]->GetY()[n]);
	}


	c= new TCanvas("c","par num",800, 600);

	mgr->Draw("al");
	mgr->SetMinimum(1.);
	mgr->GetYaxis()->SetTitleOffset(1.2);
	gPad->SetTicks(1,1);	gPad->SetLogy(1);

	lg->Draw("same");
	return true;
}




/*
for reading multiple files

read .long file
longitudinal profile & energy deposit


*/
int readLongProfile(const char* fname_long, const char* particleName, vector<double> &vX, vector<double> &vN ){


	int NXN=0; 

	TTree tr("tr","");
	NXN = tr.ReadFile( fname_long, "depth/D:gamma:positron:electron:mu:mu_bar:hadron:charged:nuclei:cherenkov", ' ');

	if(NXN <=0){
		printf("Error: no entry is read\n");
		return 0;
	}

	cout<<"Number of entries: "<< tr.GetEntries() <<endl;

	//first half: longitudinal profile
	//second half: deposited energy
	NXN = tr.Draw( Form("depth:%s", particleName), "", "goff");

	if(NXN <=0){
		printf("Error: no entries was selected\n");
		return 0;
	}

	//separate  
	for(int n=1; n<tr.GetSelectedRows(); n++){
		if( tr.GetVal(0)[n] < tr.GetVal(0)[n-1] ){ //X should increase monotonically 
			NXN = n;
			printf("break at %d: %.1f < %.1f:\t NXN = %d\n", n, tr.GetVal(0)[n], tr.GetVal(0)[n-1], NXN);

			break;
		}
	}

	//longitudinal profile only
	//assign from array
	vX.assign( tr.GetVal(0), tr.GetVal(0) + NXN );
	vN.assign( tr.GetVal(1), tr.GetVal(1) + NXN );

	//check   
	//for(int n=0; n<NXN; n++){
	//	printf("%d\t%.2f\t%.3f\n", n, vX[n], vN[n] );
	//}


	return NXN;
}

//UsefulFunction.cpp
Bool_t FileExist(const string& fname){

	if(fname.size()==0) return false;
	//if(!fname) return false;

	//Bool_t TSystem::AccessPathName( const char * path,EAccessMode mode =kFileExists )	
	//Returns FALSE if one can access a file using the specified access mode.
	//The file name must not contain any special shell characters line ~ or $, in those cases first call ExpandPathName(). Attention, bizarre convention of return value!!
	if( gSystem->AccessPathName( fname.c_str(), kFileExists) ) return false;
	else  return true;
}
