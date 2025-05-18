#include "TTree.h"
#include <iostream>
#include "TGraph.h" 
#include "TVirtualFFT.h" 
#include "TFile.h"
#include <fstream>
#include <math.h>
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include <iostream>
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TChain.h"
#include <thread>
#include <unistd.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TSpline.h"

using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8

#define WINDOW_SPACING 150
#define HALF_WINDOW_SPACING 75
#define HALF_WINDOW_WIDTH 76






class Demo{
// 宣告 public 成員
public:
    int a;
    int b;
    int DoSomething();
};	
