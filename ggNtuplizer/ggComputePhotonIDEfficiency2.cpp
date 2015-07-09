#include "TVector3.h"
#include "TH2D.h"
#include "TCut.h"
#include <math.h>
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TColor.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"
#include "TString.h"
#include "TColor.h"
#include "TLegend.h"
#include "TObject.h"
#include <string.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "TStyle.h"
#include <vector>

bool doLogPlot = false;
bool addHashErrors = false;

const int nWP = 3;
enum WpType { WP_LOOSE = 0,
	      WP_MEDIUM, 
	      WP_TIGHT};
const TString wpName[nWP] = 
  {"Loose", "Medium", "Tight"};

enum IsoType {CH_ISO, NH_ISO, PH_ISO};

float computeIso(IsoType isoType, float isoVal, float eta, float rho); 

// Use here one of the WpType values
const WpType wp = WP_MEDIUM;

const TString treename = "ggNtuplizer/EventTree";
// const TString fname1 = "~/workspace/ntuples/ggNtuple_GJets_HT-100to200_singleFile.root";
//const TString fname1 = "/tmp/rslu/job_phys14_gjet_pt40_20bx25.root";
//const TString fname1 = "/media/nik/NIKBERRYHD/job_spring14_DYJets_20bx25.root";
//const TString fname1 = "/nfs/data/eepgnnb/EGamma/job_spring14_DYJets_20bx25.root"; small sample
const TString fname1 = "/nfs/data/eepgnnb/EGamma/job_phys14_zjets_20bx25.root";
//const TString fname1 = "job_spring15_DYJetsToLL_m50_25ns.root";

bool verbose = false;
bool smallEventCount = false;

const int nEtaBins = 2;

const float ptMin = 20;
const float ptMax = 200;
const float barrelEtaLimit = 1.479;
const float endcapUpperEtaLimit = 2.5;
const float endcapLowerEtaLimit = 1.479;

//bool passesZmass(vector< vector<float> > photonCollection);

bool passWorkingPoint(WpType wp, bool isBarrel, float pt,
		      float hOverE, float full5x5_sigmaIetaIeta, 
		      float chIso, float nhIso, float phIso);

bool isMatched(float pEta, float pPhi,
	       std::vector<int> *mcPID,
	       std::vector<int> *mcMomPID,
	       std::vector<int> *mcStatus,
	       std::vector<float> *mcEta,
	       std::vector<float> *mcPhi);

void set_plot_style();
TText* doPreliminary(double x_pos,double y_pos);
void doPlot(TH1D* plot, Color_t fillColour);
TLegend* doLegend(TH1D* plot);
TH1D* doHashErrors(TH1D* plot);

void ggComputePhotonIDEfficiency2()
{

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  set_plot_style();
 // gROOT->SetBatch();
  //
  // Find the tree
  //
  TFile *file1 = new TFile(fname1);
  if( !file1 )
    assert(0);
  TTree *tree = (TTree*)file1->Get(treename);
  if( !tree )
    assert(0);
  std::cout << "Found the tree" << std::endl; fflush(stdout);
  
  //Set histos
  //Event variables
  TH1D* nPhoSigAll = new TH1D("nPhoSigAll", "", 65, 0, 13);
        nPhoSigAll->GetXaxis()->SetTitle("N_{#gamma}");
        nPhoSigAll->GetYaxis()->SetTitle("Number of events");

  TH1D* nPhoSigEB = new TH1D("nPhoSigEB", "", 65, 0, 13);
        nPhoSigEB->GetXaxis()->SetTitle("N_{#gamma}");
        nPhoSigEB->GetYaxis()->SetTitle("Number of events");

  TH1D* nPhoSigEE = new TH1D("nPhoSigEE", "", 65, 0, 13);
        nPhoSigEE->GetXaxis()->SetTitle("N_{#gamma}");
        nPhoSigEE->GetYaxis()->SetTitle("Number of events");  

  TH1D* nPhoSigAll_noconv = new TH1D("nPhoSigAll_noconv", "", 65, 0, 13);
        nPhoSigAll_noconv->GetXaxis()->SetTitle("N_{#gamma}");
        nPhoSigAll_noconv->GetYaxis()->SetTitle("Number of events");

  TH1D* nPhoSigEB_noconv = new TH1D("nPhoSigEB_noconv", "", 65, 0, 13);
        nPhoSigEB_noconv->GetXaxis()->SetTitle("N_{#gamma}");
        nPhoSigEB_noconv->GetYaxis()->SetTitle("Number of events");

  TH1D* nPhoSigEE_noconv = new TH1D("nPhoSigEE_noconv", "", 65, 0, 13);
        nPhoSigEE_noconv->GetXaxis()->SetTitle("N_{#gamma}");
        nPhoSigEE_noconv->GetYaxis()->SetTitle("Number of events");               

  TH1D* nVtxSigAll = new TH1D("nVtxSigAll", "", 210, 0, 42);
        nVtxSigAll->GetXaxis()->SetTitle("N_{nVtx}");
        nVtxSigAll->GetYaxis()->SetTitle("Number of events");

  TH1D* nVtxSigEB = new TH1D("nVtxSigEB", "", 210, 0, 42);
        nVtxSigEB->GetXaxis()->SetTitle("N_{nVtx}(barrel)");
        nVtxSigEB->GetYaxis()->SetTitle("Number of events");

  TH1D* nVtxSigEE = new TH1D("nVtxSigEE", "", 210, 0, 42);
        nVtxSigEE->GetXaxis()->SetTitle("N_{nVtx}(endcap)");
        nVtxSigEE->GetYaxis()->SetTitle("Number of events");

  TH1D* nVtxSigAll_noconv = new TH1D("nVtxSigAll_noconv", "", 210, 0, 42);
        nVtxSigAll_noconv->GetXaxis()->SetTitle("N_{nVtx}");
        nVtxSigAll_noconv->GetYaxis()->SetTitle("Number of events");

  TH1D* nVtxSigEB_noconv = new TH1D("nVtxSigEB_noconv", "", 210, 0, 42);
        nVtxSigEB_noconv->GetXaxis()->SetTitle("N_{nVtx}(barrel)");
        nVtxSigEB_noconv->GetYaxis()->SetTitle("Number of events");

  TH1D* nVtxSigEE_noconv = new TH1D("nVtxSigEE_noconv", "", 210, 0, 42);
        nVtxSigEE_noconv->GetXaxis()->SetTitle("N_{nVtx}(endcap)");
        nVtxSigEE_noconv->GetYaxis()->SetTitle("Number of events");        

  TH1D* PUSigAll = new TH1D("PUSigAll", "", 90, 0, 90);
        PUSigAll->GetXaxis()->SetTitle("PU");
        PUSigAll->GetYaxis()->SetTitle("Number of events");

  TH1D* PUSigEB = new TH1D("PUSigEB", "", 90, 0, 90);
        PUSigEB->GetXaxis()->SetTitle("PU");
        PUSigEB->GetYaxis()->SetTitle("Number of events");

  TH1D* PUSigEE = new TH1D("PUSigEE", "", 90, 0, 90);
        PUSigEE->GetXaxis()->SetTitle("PU");
        PUSigEE->GetYaxis()->SetTitle("Number of events");

  TH1D* PUSigEB_noconv = new TH1D("PUSigEB_noconv", "", 90, 0, 90);
        PUSigEB_noconv->GetXaxis()->SetTitle("PU");
        PUSigEB_noconv->GetYaxis()->SetTitle("Number of events");

  TH1D* PUSigEE_noconv = new TH1D("PUSigEE_noconv", "", 90, 0, 90);
        PUSigEE_noconv->GetXaxis()->SetTitle("PU");
        PUSigEE_noconv->GetYaxis()->SetTitle("Number of events");

  TH1D* nTrksSigAll = new TH1D("nTrksSigAll", "", 240, 0, 240);
        nTrksSigAll->GetXaxis()->SetTitle("N_{Tracks}");
        nTrksSigAll->GetYaxis()->SetTitle("Number of events");  

  TH1D* nTrksSigEB = new TH1D("nTrksSigEB", "", 240, 0, 240);
        nTrksSigEB->GetXaxis()->SetTitle("N_{Tracks}");
        nTrksSigEB->GetYaxis()->SetTitle("Number of events"); 
        
  TH1D* nTrksSigEE = new TH1D("nTrksSigEE", "", 240, 0, 240);
        nTrksSigEE->GetXaxis()->SetTitle("N_{Tracks}");
        nTrksSigEE->GetYaxis()->SetTitle("Number of events"); 

  TH1D* nTrksSigAll_noconv = new TH1D("nTrksSigAll_noconv", "", 240, 0, 240);
        nTrksSigAll_noconv->GetXaxis()->SetTitle("N_{Tracks}");
        nTrksSigAll_noconv->GetYaxis()->SetTitle("Number of events");          

  TH1D* nTrksSigEB_noconv = new TH1D("nTrksSigEB_noconv", "", 240, 0, 240);
        nTrksSigEB_noconv->GetXaxis()->SetTitle("N_{Tracks}");
        nTrksSigEB_noconv->GetYaxis()->SetTitle("Number of events"); 

  TH1D* nTrksSigEE_noconv = new TH1D("nTrksSigEE_noconv", "", 240, 0, 240);
        nTrksSigEE_noconv->GetXaxis()->SetTitle("N_{Tracks}");
        nTrksSigEE_noconv->GetYaxis()->SetTitle("Number of events");           

  //Full detector
  TH1D *phoETSigAll = new TH1D("phoETSigAll","",185, 15, 200);	
        phoETSigAll->GetXaxis()->SetTitle("E_{T} [GeV]");
        phoETSigAll->GetYaxis()->SetTitle("Number of events");    

  TH1D *phoETSigAll_noconv = new TH1D("phoETSigAll_noconv","",185, 15, 200);  
        phoETSigAll_noconv->GetXaxis()->SetTitle("E_{T} [GeV]");
        phoETSigAll_noconv->GetYaxis()->SetTitle("Number of events"); 

  TH1D *phoEtaSigAll = new TH1D("phoETaSigAll","",200, -5, 5);
        phoEtaSigAll->GetXaxis()->SetTitle("#eta_{SC} [GeV]");
        phoEtaSigAll->GetYaxis()->SetTitle("Number of events");

  TH2D *phoETvsEtaSigAll = new TH2D("phoETvsEtaSigAll","", 200, -5, 5, 185, 15, 200); 
        phoETvsEtaSigAll->GetXaxis()->SetTitle("#eta_{SC}");
        phoETvsEtaSigAll->GetYaxis()->SetTitle("E_{T} [GeV]"); 

  TH1D* phoPhiSigAll = new TH1D("phoPhiSigAll", "", 200, -4, 4);
        phoPhiSigAll->GetXaxis()->SetTitle("#phi");
        phoPhiSigAll->GetYaxis()->SetTitle("Number of events");

  
  //barrel histos
  TH1D* phoETSigEB = new TH1D("phoETSigEB", "", 185, 15, 200);
  	    phoETSigEB->GetXaxis()->SetTitle("E_{T} [GeV]");
	      phoETSigEB->GetYaxis()->SetTitle("Number of events");

  TH1D* phoETSigEB_noconv = new TH1D("phoETSigEB_noconv", "", 185, 15, 200); 
    	  phoETSigEB_noconv->GetXaxis()->SetTitle("E_{T} [GeV]");
	      phoETSigEB_noconv->GetYaxis()->SetTitle("Number of events");

  TH1D* phoEtaSigEB = new TH1D("phoEtaSigEB", "", 200, -3, 3);
        phoEtaSigEB->GetXaxis()->SetTitle("#eta_{SC}");
        phoEtaSigEB->GetYaxis()->SetTitle("Number of events");

  TH1D* phoEtaSigEB_noconv = new TH1D("phoEtaSigEB_noconv", "", 200, -3, 3);
        phoEtaSigEB_noconv->GetXaxis()->SetTitle("#eta_{SC}");
        phoEtaSigEB_noconv->GetYaxis()->SetTitle("Number of events");         
  
  TH1D* phoPhiSigEB = new TH1D("phoPhiSigEB", "", 200, -4, 4);
        phoPhiSigEB->GetXaxis()->SetTitle("#phi");
        phoPhiSigEB->GetYaxis()->SetTitle("Number of events");

  TH1D* phoPhiSigEB_noconv = new TH1D("phoPhiSigEB_noconv", "", 200, -4, 4);
        phoPhiSigEB_noconv->GetXaxis()->SetTitle("#phi");
        phoPhiSigEB_noconv->GetYaxis()->SetTitle("Number of events");

  //endcap histos
  TH1D* phoETSigEE = new TH1D("phoETSigEE", "", 185, 15, 200);
    	  phoETSigEE->GetXaxis()->SetTitle("E_{T} [GeV]");
	      phoETSigEE->GetYaxis()->SetTitle("Number of events");

  TH1D* phoETSigEE_noconv = new TH1D("phoETSigEE_noconv", "", 185, 15, 200);  
    	  phoETSigEE_noconv->GetXaxis()->SetTitle("E_{T} [GeV]");
	      phoETSigEE_noconv->GetYaxis()->SetTitle("Number of events");

  TH1D* phoEtaSigEE = new TH1D("phoEtaSigEE", "", 200, -3, 3);
        phoEtaSigEE->GetXaxis()->SetTitle("#eta_{SC}");
        phoEtaSigEE->GetYaxis()->SetTitle("Number of events");

  TH1D* phoEtaSigEE_noconv = new TH1D("phoEtaSigEE_noconv", "", 200, -3, 3);
        phoEtaSigEE_noconv->GetXaxis()->SetTitle("#eta_{SC}");
        phoEtaSigEE_noconv->GetYaxis()->SetTitle("Number of events");      

  TH1D* phoPhiSigEE = new TH1D("phoPhiSigEE", "", 200, -4, 4);
        phoPhiSigEE->GetXaxis()->SetTitle("#phi");
        phoPhiSigEE->GetYaxis()->SetTitle("Number of events");

  TH1D* phoPhiSigEE_noconv = new TH1D("phoPhiSigEE_noconv", "", 200, -4, 4);
        phoPhiSigEE_noconv->GetXaxis()->SetTitle("#phi");
        phoPhiSigEE_noconv->GetYaxis()->SetTitle("Number of events"); 

  TH2D* FRvsET = new TH2D("FRvsET", "", 185, 15, 200, 185, 15, 200);
        FRvsET->GetXaxis()->SetTitle("E_{T} [GeV]");
        FRvsET->GetYaxis()->SetTitle("Fake Rate");     
             
  
  // Event-level variables:   
  int nPho;
  float rho;
  std::vector<double> *nPU = 0;
  int nTrksPV;  
  int nVtx;
  //std::vector<int> *nVtx = 0;

  // Per-photon variables
  // Kinematics
  std::vector <float> *phoEt = 0;     
  std::vector <float> *phoSCEta = 0;    
  std::vector <float> *phoPhi = 0;  
    
  // ID related variables
  std::vector <float> *phoSigmaIEtaIEta_2012 = 0; 
  std::vector <float> *phoHoverE = 0;    
  std::vector <float> *phoPFPhoIso = 0;    
  std::vector <float> *phoPFChIso = 0;    
  std::vector <float> *phoPFNeuIso = 0;    
  std::vector <int> *phohasPixelSeed = 0; 
  std::vector <int> *phoEleVeto = 0;
  
  // Electron variables
  std::vector <int> *eleConvVeto = 0;  
  
  // MC variables  
  std::vector <int> *mcPID = 0;     
  std::vector <int> *mcMomPID = 0;     
  std::vector <int> *mcStatus = 0;     
  std::vector <float> *mcEta = 0;     
  std::vector <float> *mcPhi = 0;     

  // Declare branches
  TBranch *b_nPho = 0;
  TBranch *b_rho = 0;
  TBranch *b_nPU = 0;
  TBranch *b_nTrksPV = 0;
  TBranch *b_nVtx = 0;
  TBranch *b_phoEt = 0;
  TBranch *b_phoSCEta = 0;
  TBranch *b_phoPhi = 0;
  TBranch *b_phoSigmaIEtaIEta_2012 = 0;
  TBranch *b_phoHoverE = 0;
  TBranch *b_phoPFPhoIso = 0;
  TBranch *b_phoPFChIso = 0;
  TBranch *b_phoPFNeuIso = 0;
  TBranch *b_phohasPixelSeed = 0;
  TBranch *b_phoEleVeto = 0;
  //
  TBranch *b_mcPID;
  TBranch *b_mcMomPID;
  TBranch *b_mcStatus;
  TBranch *b_mcEta;
  TBranch *b_mcPhi;

  TBranch *b_eleConvVeto = 0;

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nPho", &nPho, &b_nPho);
  tree->SetBranchAddress("rho", &rho, &b_rho);
  tree->SetBranchAddress("nPU", &nPU, &b_nPU);
  tree->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
  tree->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
  tree->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
  tree->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
  tree->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
  tree->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012, &b_phoSigmaIEtaIEta_2012);
  tree->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
  tree->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
  tree->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
  tree->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
  tree->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
  tree->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
  //
  tree->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
  tree->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
  tree->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
  tree->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
  tree->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);

  tree->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   //
  // The first loop is to fill the weight histograms Weight histograms
  //
  TH2D *hSignal = new TH2D("hSignal","",200,-5,5,185,15,200);
  hSignal->GetXaxis()->SetTitle("#eta_{SC}");
  hSignal->GetYaxis()->SetTitle("E_{T} [GeV]");
  TH2D *hBackground = new TH2D("hBackground","",200,-5,5,185,15,200);
  hBackground->GetXaxis()->SetTitle("#eta_{SC}");
  hBackground->GetYaxis()->SetTitle("E_{T} [GeV]"); 
 
  UInt_t maxEvents = tree->GetEntries();
  if( smallEventCount )
    maxEvents = 1000;
  if(verbose)
    std::cout << "Start loop over events for WEIGHTS, total events = " << tree->GetEntries() << std::endl; 
           
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){
    
    if( ievent%100000 == 0){
    	std::cout << "." ; fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);
    
    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    b_nPU->GetEntry(tentry);
    b_nTrksPV->GetEntry(tentry);
    b_nVtx->GetEntry(tentry);

    // Get data for all photons in this event, only vars of interest
    b_phoEt->GetEntry(tentry);
    b_phoSCEta->GetEntry(tentry);
    b_phoPhi->GetEntry(tentry);
    b_phohasPixelSeed->GetEntry(tentry);
    b_phoEleVeto->GetEntry(tentry);

    b_mcPID->GetEntry(tentry);
    b_mcMomPID->GetEntry(tentry);
    b_mcStatus->GetEntry(tentry);
    b_mcEta->GetEntry(tentry);
    b_mcPhi->GetEntry(tentry);

    b_eleConvVeto->GetEntry(tentry);
    

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){
   
      
      // Preselection
      if( !(phoEt->at(ipho) > ptMin && phoEt->at(ipho) < ptMax ) ) continue;

      // Match to MC truth
      bool isTrue = isMatched( phoSCEta->at(ipho), phoPhi->at(ipho),
			       mcPID, mcMomPID, mcStatus,
			       mcEta, mcPhi);

      if( isTrue ) {
	hSignal->Fill( phoSCEta->at(ipho), phoEt->at(ipho) );

      }else{
	hBackground->Fill( phoSCEta->at(ipho), phoEt->at(ipho) );
      }
    }// end loop over photons
  } // end loop over events


TCanvas* can1 = new TCanvas("can1");
	hSignal->Draw("colz");
can1->SaveAs("Plots/Signal.pdf");

TCanvas* can2 = new TCanvas("can2");
	hBackground->Draw("colz");
can2->SaveAs("Plots/Background.pdf"); 

  //
  // Prepare for computing efficiencies
  //
  double sumSignalDenomEB = 0;
  double sumSignalNumEB   = 0;
  double sumSignalDenomEE = 0;
  double sumSignalNumEE   = 0;
  double sumBackDenomEB = 0;
  double sumBackNumEB   = 0;
  double sumBackDenomEE = 0;
  double sumBackNumEE   = 0;

  double sumSignalDenomEBErr2 = 0;
  double sumSignalNumEBErr2   = 0;
  double sumSignalDenomEEErr2 = 0;
  double sumSignalNumEEErr2   = 0;
  double sumBackDenomEBErr2 = 0;
  double sumBackNumEBErr2   = 0;
  double sumBackDenomEEErr2 = 0;
  double sumBackNumEEErr2   = 0;

  int nSigBarrel = 0;
  float nSigBarrelWeighted = 0;
  int nBgBarrel  = 0;
  float nBgBarrelWeighted  = 0;

  // 
  // The second loop over events is for efficiency computation
  //
  
  if(verbose)
    std::cout << "Start loop over events for EFFICIENCY caluclation, total events = " << tree->GetEntries() << std::endl; 
           
 vector< vector<float> > photonCollection;	   
 vector<float> photon;	   
	   
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      std::cout << "." ; fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);

    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    b_nPU->GetEntry(tentry);
    b_nTrksPV->GetEntry(tentry);
    b_nVtx->GetEntry(tentry);

    if(verbose)
      std::cout << "Event " << ievent << ", number of photons " << nPho << std::endl;

    // Get data for all photons in this event, only vars of interest
    b_phoEt->GetEntry(tentry);
    b_phoSCEta->GetEntry(tentry);
    b_phoPhi->GetEntry(tentry);
    b_phoSigmaIEtaIEta_2012->GetEntry(tentry);
    b_phoHoverE->GetEntry(tentry);
    b_phoPFPhoIso->GetEntry(tentry);
    b_phoPFChIso->GetEntry(tentry);
    b_phoPFNeuIso->GetEntry(tentry);
    b_phohasPixelSeed->GetEntry(tentry);
    b_phoEleVeto->GetEntry(tentry);
    b_eleConvVeto->GetEntry(tentry);

    b_mcPID->GetEntry(tentry);
    b_mcMomPID->GetEntry(tentry);
    b_mcStatus->GetEntry(tentry);
    b_mcEta->GetEntry(tentry);
    b_mcPhi->GetEntry(tentry);

//std::cout << "nPU at ievent = " << nPU->at(ievent) << std::endl;
//std::cout << "nVtx at ievent = " << nVtx->at(ievent) << std::endl;
  // PUSigAll->Fill(nPU->at(ievent));

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){
      
      // Preselection
      if( !(phoEt->at(ipho) > ptMin && phoEt->at(ipho) < ptMax ) ) continue;
      if( !( phoEleVeto->at(ipho) == 1) ) continue;
//      if( !( phohasPixelSeed->at(ipho) == 0 ) ) continue;
     // if( nPho < 2 ) continue; 
         
//      bool withinZMassRange = passesZmass(photonCollection);


      bool isBarrel = (fabs(phoSCEta->at(ipho)) < barrelEtaLimit);
      bool isEndcap = (endcapLowerEtaLimit < fabs(phoSCEta->at(ipho)) && fabs(phoSCEta->at(ipho)) < endcapUpperEtaLimit);

      // Correct for pile-up
      float chIsoWithEA = computeIso(CH_ISO, 
				     phoPFChIso->at(ipho), 
				     phoSCEta->at(ipho), rho);
      float nhIsoWithEA = computeIso(NH_ISO, 
				     phoPFNeuIso->at(ipho), 
				     phoSCEta->at(ipho), rho);
      float phIsoWithEA = computeIso(PH_ISO, 
				     phoPFPhoIso->at(ipho), 
				     phoSCEta->at(ipho), rho);

      // Compute ID decision
      bool pass = passWorkingPoint( wp, isBarrel, phoEt->at(ipho),
				    phoHoverE->at(ipho),
				    phoSigmaIEtaIEta_2012->at(ipho),
				    chIsoWithEA, nhIsoWithEA, phIsoWithEA);

      // Match to MC truth
      bool isTrue = isMatched( phoSCEta->at(ipho), phoPhi->at(ipho),
			       mcPID, mcMomPID, mcStatus,
			       mcEta, mcPhi);
      
      // We reweight signal and background (separately) to have 
      // a flat pt and eta distribution. This step is a matter of definition.
      double binContent = 0;
      TH2D *hEvents = hSignal;
      if( !isTrue )
	hEvents = hBackground;
      binContent = hEvents->GetBinContent
	( hEvents->FindBin( phoSCEta->at(ipho), phoEt->at(ipho) ) );
      double weight = 1;
      if( binContent == 0 ){
      	std::cout << "Error! Zero! pt= "<< phoEt->at(ipho) << ", eta= " << phoSCEta->at(ipho) << std::endl;
      }else{
      	weight = 1./binContent;
      }
      

  if( isTrue ) {
	 phoETSigAll->Fill(phoEt->at(ipho), weight);
	 phoEtaSigAll->Fill(phoSCEta->at(ipho), weight);
   phoPhiSigAll->Fill(phoPhi->at(ipho), weight);
	 phoETvsEtaSigAll->Fill(phoSCEta->at(ipho), phoEt->at(ipho), weight);

    std::vector<int> nVtxAllVec; 
    nVtxAllVec.push_back(nVtx);
    for (unsigned int j = 0; j < nVtxAllVec.size(); ++j)
    {
    // std::cout << "nVtx at event = " << nVtxVec[j] << std::endl;
      nVtxSigAll->Fill(nVtxAllVec[j], weight);

    }
    std::vector<int> nPhoAllVec; 
    nPhoAllVec.push_back(nPho);
    for (unsigned int j = 0; j < nPhoAllVec.size(); ++j)
    {
    // std::cout << "nPho at event = " << nPhoVec[j] << std::endl;
      nPhoSigAll->Fill(nPhoAllVec[j], weight);

    }

    std::vector<int> nTrksAllVec; 
    nTrksAllVec.push_back(nTrksPV);
    for (unsigned int j = 0; j < nTrksAllVec.size(); ++j)
    {
    // std::cout << "nTrks at event = " << nTrksVec[j] << std::endl;
      nTrksSigAll->Fill(nTrksAllVec[j], weight);

    }

	 if( isBarrel ) {
	  nSigBarrel++;
	  nSigBarrelWeighted += weight;
	  sumSignalDenomEB += weight;
	  sumSignalDenomEBErr2 += weight*weight;

	  if( pass ) {
	    sumSignalNumEB += weight;
	    sumSignalNumEBErr2 += weight*weight;
	    
	    phoETSigEB->Fill(phoEt->at(ipho), weight);
      phoEtaSigEB->Fill(phoSCEta->at(ipho), weight);
      phoPhiSigEB->Fill(phoPhi->at(ipho), weight);
      PUSigEB->Fill(nPU->at(ipho), weight);

      std::vector<double> nVtxSigEBVec;
      nVtxSigEBVec.push_back(nVtx);
      for (unsigned int j = 0; j < nVtxSigEBVec.size(); ++j)
      {
      // std::cout << "nVtx at event = " << nVtxVec[j] << std::endl;
      nVtxSigEB->Fill(nVtxSigEBVec[j], weight);

      }

      std::vector<double> nPhoSigEBVec;
      nPhoSigEBVec.push_back(nPho);
      for (unsigned int j = 0; j < nPhoSigEBVec.size(); ++j)
      {
      // std::cout << "nPho at event = " << nPhoVec[j] << std::endl;
      nPhoSigEB->Fill(nPhoSigEBVec[j], weight);

      }

      std::vector<double> nTrksSigEBVec;
      nTrksSigEBVec.push_back(nTrksPV);
      for (unsigned int j = 0; j < nTrksSigEBVec.size(); ++j)
      {
      // std::cout << "nTrks at event = " << nTrksVec[j] << std::endl;
      nTrksSigEB->Fill(nTrksSigEBVec[j], weight);

      }

	  }
	 } else if( isEndcap ){

	  sumSignalDenomEE += weight;
	  sumSignalDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	    sumSignalNumEE += weight;
	    sumSignalNumEEErr2 += weight*weight;
	    
	    phoETSigEE->Fill(phoEt->at(ipho), weight);
      phoEtaSigEE->Fill(phoSCEta->at(ipho), weight);
      phoPhiSigEE->Fill(phoPhi->at(ipho), weight);
      PUSigEE->Fill(nPU->at(ipho), weight);

      std::vector<double> nVtxSigEEVec;
      nVtxSigEEVec.push_back(nVtx);
      for (unsigned int j = 0; j < nVtxSigEEVec.size(); ++j)
      {
      // std::cout << "nVtx at event = " << nVtxVec[j] << std::endl;
        nVtxSigEE->Fill(nVtxSigEEVec[j], weight);

      }

      std::vector<double> nPhoSigEEVec;
      nPhoSigEEVec.push_back(nPho);
      for (unsigned int j = 0; j < nPhoSigEEVec.size(); ++j)
      {
      // std::cout << "nPho at event = " << nPhoVec[j] << std::endl;
        nPhoSigEE->Fill(nPhoSigEEVec[j], weight);

      }

      std::vector<double> nTrksSigEEVec;
      nTrksSigEEVec.push_back(nTrksPV);
      for (unsigned int j = 0; j < nTrksSigEEVec.size(); ++j)
      {
      // std::cout << "nTrks at event = " << nTrksVec[j] << std::endl;
        nTrksSigEE->Fill(nTrksSigEEVec[j], weight);

      }

	  }

	 }// end barrel / endcap
    } // end if signal
	

      if( !isTrue ) {
      
 	
	if( isBarrel ) {
	  nBgBarrel++;
	  nBgBarrelWeighted += weight;
	  
	  sumBackDenomEB += weight;
	  sumBackDenomEBErr2 += weight*weight;
	  if( pass ) {
	    sumBackNumEB += weight;
	    sumBackNumEBErr2 += weight*weight;
	  }
	} else if( isEndcap ){
	  sumBackDenomEE += weight;
	  sumBackDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	    sumBackNumEE += weight;
	    sumBackNumEEErr2 += weight*weight;
	  }

	}// end barrel / endcap
      } // end if background
    } // end loop over photons

}// end loop over events  


///=========Third loop for full selection sans Electron Conversion Veto===========

for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      std::cout << "." ; fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);

    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    b_nPU->GetEntry(tentry);

    if(verbose)
      std::cout << "Event " << ievent << ", number of photons " << nPho << std::endl;

    // Get data for all photons in this event, only vars of interest
    b_phoEt->GetEntry(tentry);
    b_phoSCEta->GetEntry(tentry);
    b_phoPhi->GetEntry(tentry);
    b_phoSigmaIEtaIEta_2012->GetEntry(tentry);
    b_phoHoverE->GetEntry(tentry);
    b_phoPFPhoIso->GetEntry(tentry);
    b_phoPFChIso->GetEntry(tentry);
    b_phoPFNeuIso->GetEntry(tentry);
    b_phohasPixelSeed->GetEntry(tentry);
    b_phoEleVeto->GetEntry(tentry);
    b_eleConvVeto->GetEntry(tentry);

    b_mcPID->GetEntry(tentry);
    b_mcMomPID->GetEntry(tentry);
    b_mcStatus->GetEntry(tentry);
    b_mcEta->GetEntry(tentry);
    b_mcPhi->GetEntry(tentry);

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){
      
      // Preselection
      if( !(phoEt->at(ipho) > ptMin && phoEt->at(ipho) < ptMax ) ) continue;
      if( !( phoEleVeto->at(ipho) == 0) ) continue;
////      if( !( phohasPixelSeed->at(ipho) == 0 ) ) continue;
//      if( nPho < 2 ) continue; 

//      bool withinZMassRange = passesZmass();

      bool isBarrel = (fabs(phoSCEta->at(ipho)) < barrelEtaLimit);
      bool isEndcap = (endcapLowerEtaLimit < fabs(phoSCEta->at(ipho)) && fabs(phoSCEta->at(ipho)) < endcapUpperEtaLimit);

      // Correct for pile-up
      float chIsoWithEA = computeIso(CH_ISO, 
				     phoPFChIso->at(ipho), 
				     phoSCEta->at(ipho), rho);
      float nhIsoWithEA = computeIso(NH_ISO, 
				     phoPFNeuIso->at(ipho), 
				     phoSCEta->at(ipho), rho);
      float phIsoWithEA = computeIso(PH_ISO, 
				     phoPFPhoIso->at(ipho), 
				     phoSCEta->at(ipho), rho);

      // Compute ID decision
      bool pass = passWorkingPoint( wp, isBarrel, phoEt->at(ipho),
				    phoHoverE->at(ipho),
				    phoSigmaIEtaIEta_2012->at(ipho),
				    chIsoWithEA, nhIsoWithEA, phIsoWithEA);

      // Match to MC truth
      bool isTrue = isMatched( phoSCEta->at(ipho), phoPhi->at(ipho),
			       mcPID, mcMomPID, mcStatus,
			       mcEta, mcPhi);
      
      // We reweight signal and background (separately) to have 
      // a flat pt and eta distribution. This step is a matter of definition.
      double binContent = 0;
      TH2D *hEvents = hSignal;
      if( !isTrue )
	hEvents = hBackground;
      binContent = hEvents->GetBinContent( hEvents->FindBin( phoSCEta->at(ipho), phoEt->at(ipho) ) );
      double weight = 1;
      if( binContent == 0 ){
      	std::cout << "Error! Zero! pt= "<< phoEt->at(ipho) << ", eta= " << phoSCEta->at(ipho) << std::endl;
      }else{
      	weight = 1./binContent;
      }
      

  if( isTrue ) {
	//phoEtCheckHist_conv->Fill(phoEt->at(ipho), weight);
	//phoSCEtaCheckHist_conv->Fill(phoSCEta->at(ipho), weight);
	//phoTotCheck_conv->Fill(phoSCEta->at(ipho), phoEt->at(ipho), weight);
    phoETSigAll_noconv->Fill(phoEt->at(ipho), weight);

    std::vector<int> nVtxAllVec; 
    nVtxAllVec.push_back(nVtx);
    for (unsigned int j = 0; j < nVtxAllVec.size(); ++j)
    {
    // std::cout << "nVtx at event = " << nVtxVec[j] << std::endl;
      nVtxSigAll_noconv->Fill(nVtxAllVec[j], weight);

    }
    std::vector<int> nPhoAllVec; 
    nPhoAllVec.push_back(nPho);
    for (unsigned int j = 0; j < nPhoAllVec.size(); ++j)
    {
    // std::cout << "nPho at event = " << nPhoVec[j] << std::endl;
      nPhoSigAll_noconv->Fill(nPhoAllVec[j], weight);

    }

    std::vector<int> nTrksAllVec; 
    nTrksAllVec.push_back(nTrksPV);
    for (unsigned int j = 0; j < nTrksAllVec.size(); ++j)
    {
    // std::cout << "nTrks at event = " << nTrksVec[j] << std::endl;
      nTrksSigAll_noconv->Fill(nTrksAllVec[j], weight);

    }  
	  	

	if( isBarrel ) {
	  //nSigBarrel++;
	  //nSigBarrelWeighted += weight;
	  //sumSignalDenomEB += weight;
	  //sumSignalDenomEBErr2 += weight*weight;


	  if( pass ) {
	    //sumSignalNumEB += weight;
	    //sumSignalNumEBErr2 += weight*weight;
	    
	    phoETSigEB_noconv->Fill(phoEt->at(ipho), weight);
      phoEtaSigEB_noconv->Fill(phoSCEta->at(ipho), weight);
      phoPhiSigEB_noconv->Fill(phoPhi->at(ipho), weight);
      PUSigEB_noconv->Fill(nPU->at(ipho), weight);
      //nTrksSigEB_noconv->Fill(nTrksPV->at(ipho), weight);

	  }
	} else if( isEndcap ){
	  //sumSignalDenomEE += weight;
	  //sumSignalDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	   // sumSignalNumEE += weight;
	   // sumSignalNumEEErr2 += weight*weight;
	    
	    phoETSigEE_noconv->Fill(phoEt->at(ipho), weight);
      phoEtaSigEE_noconv->Fill(phoSCEta->at(ipho), weight);
      phoPhiSigEE_noconv->Fill(phoPhi->at(ipho), weight);
      PUSigEE_noconv->Fill(nPU->at(ipho), weight);
      //nTrksSigEE_noconv->Fill(nTrksPV->at(ipho), weight);

	  }

	}// end barrel / endcap
      } // end if signal
	

      if( !isTrue ) {
      
 	
	if( isBarrel ) {
	  //nBgBarrel++;
	  //nBgBarrelWeighted += weight;
	  
	  //sumBackDenomEB += weight;
	  //sumBackDenomEBErr2 += weight*weight;
	  if( pass ) {
	  //  sumBackNumEB += weight;
	  //  sumBackNumEBErr2 += weight*weight;
	  }
	} else if( isEndcap ){
	//  sumBackDenomEE += weight;
	//  sumBackDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	//    sumBackNumEE += weight;
	//    sumBackNumEEErr2 += weight*weight;
	  }

	}// end barrel / endcap
      } // end if background
    } // end loop over photons
    
 //   N_Pho_noconv->Fill(nPho->at(ipho), weight);
 //   PU_noconv->Fill(nPU->at(ipho), weight);
 //   N_Trks_noconv->Fill(nTrksPV->at(ipho), weight);

}// end loop over events  

///=========Fourth loop for creating 2D plot for fake rate with a variable===========

for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      std::cout << "." ; fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);

    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    b_nPU->GetEntry(tentry);
    b_nTrksPV->GetEntry(tentry);

    if(verbose)
      std::cout << "Event " << ievent << ", number of photons " << nPho << std::endl;

    // Get data for all photons in this event, only vars of interest
    b_phoEt->GetEntry(tentry);
    b_phoSCEta->GetEntry(tentry);
    b_phoPhi->GetEntry(tentry);
    b_phoSigmaIEtaIEta_2012->GetEntry(tentry);
    b_phoHoverE->GetEntry(tentry);
    b_phoPFPhoIso->GetEntry(tentry);
    b_phoPFChIso->GetEntry(tentry);
    b_phoPFNeuIso->GetEntry(tentry);
    b_phohasPixelSeed->GetEntry(tentry);
    b_phoEleVeto->GetEntry(tentry);
    b_eleConvVeto->GetEntry(tentry);

    b_mcPID->GetEntry(tentry);
    b_mcMomPID->GetEntry(tentry);
    b_mcStatus->GetEntry(tentry);
    b_mcEta->GetEntry(tentry);
    b_mcPhi->GetEntry(tentry);

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){
      
      // Preselection
      if( !(phoEt->at(ipho) > ptMin && phoEt->at(ipho) < ptMax ) ) continue;
//      if( !( phoEleVeto->at(ipho) == 0) ) continue;
//      if( !( phohasPixelSeed->at(ipho) == 0 ) ) continue;
//      if( nPho < 2 ) continue; 

//      bool withinZMassRange = passesZmass();

      bool isBarrel = (fabs(phoSCEta->at(ipho)) < barrelEtaLimit);
      bool isEndcap = (endcapLowerEtaLimit < fabs(phoSCEta->at(ipho)) && fabs(phoSCEta->at(ipho)) < endcapUpperEtaLimit);

      // Correct for pile-up
      float chIsoWithEA = computeIso(CH_ISO, 
             phoPFChIso->at(ipho), 
             phoSCEta->at(ipho), rho);
      float nhIsoWithEA = computeIso(NH_ISO, 
             phoPFNeuIso->at(ipho), 
             phoSCEta->at(ipho), rho);
      float phIsoWithEA = computeIso(PH_ISO, 
             phoPFPhoIso->at(ipho), 
             phoSCEta->at(ipho), rho);

      // Compute ID decision
      bool pass = passWorkingPoint( wp, isBarrel, phoEt->at(ipho),
            phoHoverE->at(ipho),
            phoSigmaIEtaIEta_2012->at(ipho),
            chIsoWithEA, nhIsoWithEA, phIsoWithEA);

      // Match to MC truth
      bool isTrue = isMatched( phoSCEta->at(ipho), phoPhi->at(ipho),
             mcPID, mcMomPID, mcStatus,
             mcEta, mcPhi);
      
      // We reweight signal and background (separately) to have 
      // a flat pt and eta distribution. This step is a matter of definition.
      double binContent = 0;
      TH2D *hEvents = hSignal;
      if( !isTrue )
  hEvents = hBackground;
      binContent = hEvents->GetBinContent( hEvents->FindBin( phoSCEta->at(ipho), phoEt->at(ipho) ) );
      double weight = 1;
      if( binContent == 0 ){
        std::cout << "Error! Zero! pt= "<< phoEt->at(ipho) << ", eta= " << phoSCEta->at(ipho) << std::endl;
      }else{
        weight = 1./binContent;
      }
      
float ET, ET_noconv;

  if( isTrue && pass && isBarrel && (phoEleVeto->at(ipho) == 0)){

        ET = phoEt->at(ipho);
  }

  if( isTrue && pass && isBarrel && (phoEleVeto->at(ipho) == 1)){
    
        ET_noconv = phoEt->at(ipho);
  }

  float denom = ET + ET_noconv;
  float fakerate = ET/denom;
  //std::cout << "fakerate " << fakerate << std::endl;
  //std::cout << "ET " << ET << " ET_noconv " << ET_noconv << " ET+ET_noconv " << ET+ET_noconv <<
               //" ET/(ET+ET_noconv) " << ET/(ET+ET_noconv) << " ET/denom " << ET/denom << std::endl;
  FRvsET->Fill(ET, fakerate, weight);

  if( isTrue ) {
      

  if( isBarrel ) {
    //nSigBarrel++;
    //nSigBarrelWeighted += weight;
    //sumSignalDenomEB += weight;
    //sumSignalDenomEBErr2 += weight*weight;

    if( pass ) {
      //sumSignalNumEB += weight;
      //sumSignalNumEBErr2 += weight*weight;
    

    }
  } else if( isEndcap ){
    //sumSignalDenomEE += weight;
    //sumSignalDenomEEErr2 += weight*weight;    
    if( pass ) {
     // sumSignalNumEE += weight;
     // sumSignalNumEEErr2 += weight*weight;
      

    }

  }// end barrel / endcap
      } // end if signal
  

      if( !isTrue ) {
      
  
  if( isBarrel ) {
    //nBgBarrel++;
    //nBgBarrelWeighted += weight;
    
    //sumBackDenomEB += weight;
    //sumBackDenomEBErr2 += weight*weight;
    if( pass ) {
    //  sumBackNumEB += weight;
    //  sumBackNumEBErr2 += weight*weight;
    }
  } else if( isEndcap ){
  //  sumBackDenomEE += weight;
  //  sumBackDenomEEErr2 += weight*weight;    
    if( pass ) {
  //    sumBackNumEE += weight;
  //    sumBackNumEEErr2 += weight*weight;
    }

  }// end barrel / endcap
      } // end if background
    } // end loop over photons
    
 //   N_Pho_noconv->Fill(nPho->at(ipho), weight);
 //   PU_noconv->Fill(nPU->at(ipho), weight);
 //   N_Trks_noconv->Fill(nTrksPV->at(ipho), weight);

}// end loop over events  

TH1D* FRSigEB = (TH1D*)phoETSigEB->Clone("FRSigEB");
      FRSigEB->GetXaxis()->SetTitle("Fake Rate");  
      FRSigEB->GetYaxis()->SetTitle("Number of events"); 
      FRSigEB->Divide(phoETSigEB_noconv);

TH1D* FRSigEE = (TH1D*)phoETSigEE->Clone("FRSigEE");
      FRSigEE->GetXaxis()->SetTitle("Fake Rate");  
      FRSigEE->GetYaxis()->SetTitle("Number of events"); 
      FRSigEE->Divide(phoETSigEE_noconv);  

TH1D* FRSigAll = (TH1D*)phoETSigAll->Clone("FRSigAll");
      FRSigAll->GetXaxis()->SetTitle("Fake Rate");  
      FRSigAll->GetYaxis()->SetTitle("Number of events"); 
      FRSigAll->Divide(phoETSigAll_noconv);        

  std::cout << "" << std::endl;
  std::cout << "***********************EB Signal Events*************************" << std::endl;
  std::cout << "Number of signal events = " << phoETSigEB->Integral() << std::endl;
  std::cout << "Number of events passing selection but not conversion Veto = " << phoETSigEB_noconv->Integral() << std::endl;
  std::cout << "Number of events in ratio of the two = " << FRSigEB->Integral() << std::endl;
  std::cout << "Ratio of events with Conversion Veto to events without, R = " << phoETSigEB->Integral()/phoETSigEB_noconv->Integral() << " (barrel)" << std::endl;

  std::cout << "***********************EE Signal Events*************************" << std::endl;
  std::cout << "Number of signal events = " << phoETSigEE->Integral() << std::endl;
  std::cout << "Number of events passing selection but not conversion Veto = " << phoETSigEE_noconv->Integral() << std::endl;
  std::cout << "Number of events in ratio of the two = " << FRSigEE->Integral() << std::endl;
  std::cout << "Ratio of events with Conversion Veto to events without, R = " << phoETSigEE->Integral()/phoETSigEE_noconv->Integral() << " (endcaps)" << std::endl;

  std::vector<TH1D*> barrelHistos;

  barrelHistos.push_back(phoETSigEB);
  barrelHistos.push_back(phoETSigEB_noconv);
  barrelHistos.push_back(phoEtaSigEB);
  barrelHistos.push_back(phoEtaSigEB_noconv);
  barrelHistos.push_back(phoPhiSigEB);
  barrelHistos.push_back(phoPhiSigEB_noconv);
  barrelHistos.push_back(PUSigEB);
  barrelHistos.push_back(PUSigEB_noconv); 
  barrelHistos.push_back(nTrksSigEB);
  barrelHistos.push_back(nTrksSigEB_noconv);
  barrelHistos.push_back(FRSigEB);
  barrelHistos.push_back(nVtxSigEB);
  barrelHistos.push_back(nVtxSigEB_noconv);
  barrelHistos.push_back(nPhoSigEB);
  barrelHistos.push_back(nPhoSigEB_noconv);
  std::vector<TH1D*> endcapHistos;

  endcapHistos.push_back(phoETSigEE);
  endcapHistos.push_back(phoETSigEE_noconv);
  endcapHistos.push_back(phoEtaSigEE);
  endcapHistos.push_back(phoEtaSigEE_noconv);
  endcapHistos.push_back(phoPhiSigEE);
  endcapHistos.push_back(phoPhiSigEE_noconv);
  endcapHistos.push_back(PUSigEE);
  endcapHistos.push_back(PUSigEE_noconv);
  endcapHistos.push_back(nTrksSigEE);
  endcapHistos.push_back(nTrksSigEE_noconv);
  endcapHistos.push_back(FRSigEE);
  endcapHistos.push_back(nVtxSigEE);
  endcapHistos.push_back(nVtxSigEE_noconv);  
  endcapHistos.push_back(nPhoSigEE);
  endcapHistos.push_back(nPhoSigEE_noconv);
  std::vector<TH1D*> fullHistos;

  fullHistos.push_back(phoETSigAll);
//  endcapHistos.push_back(phoETSigAll_noconv);
  fullHistos.push_back(phoEtaSigAll);
//  endcapHistos.push_back(phoEtaSigAll_noconv);
  fullHistos.push_back(phoPhiSigAll);
//  endcapHistos.push_back(phoPhiSigAll_noconv);
  fullHistos.push_back(PUSigAll);
  fullHistos.push_back(nTrksSigAll);
  fullHistos.push_back(nTrksSigAll_noconv);  
  fullHistos.push_back(FRSigAll);
  fullHistos.push_back(nVtxSigAll);
  fullHistos.push_back(nVtxSigAll_noconv);  
  fullHistos.push_back(nPhoSigAll);
  fullHistos.push_back(nPhoSigAll_noconv);  

// Plotting section
for (unsigned int i = 0; i < barrelHistos.size(); ++i)
{
  doPlot(barrelHistos[i], kAzure+8);
}

for (unsigned int i = 0; i < endcapHistos.size(); ++i)
{
  doPlot(endcapHistos[i], kTeal+8);
}

for (unsigned int i = 0; i < fullHistos.size(); ++i)
{
  doPlot(fullHistos[i], kRed+1);
}

  TCanvas *c2 = new TCanvas("c2","c2",10,10,900, 600);

          phoETvsEtaSigAll->Draw("colz");

  c2->SaveAs("Plots/phoETvsEtaSigAll.pdf");

  TCanvas* c3 = new TCanvas("c3");
          FRvsET->Draw();
  c3->SaveAs("Plots/FRvsET.pdf");        


  std::cout << "barrel signal a=" << sumSignalNumEB << "   b=" << sumSignalDenomEB << std::endl;
  std::cout << "barrel background a=" << sumBackNumEB << "    b=" << sumBackDenomEB << std::endl;

  std::cout << "Efficiencies for the working point " << wpName[wp].Data() << std::endl;

  // Compute signal efficiencies
  double effSignalEB = sumSignalNumEB / sumSignalDenomEB;
  double effSignalEBErr = sqrt( sumSignalDenomEBErr2 
				* effSignalEB*(1-effSignalEB)
				/(sumSignalDenomEB*sumSignalDenomEB) );
  std::cout << "Signal barrel efficiency: " << effSignalEB*100 << " +- " << effSignalEBErr*100 << " %%" << std::endl; //  std::cout << "Signal barrel efficiency: %5.1f +- %5.1f %%" << std::endl; 
	 

  double effSignalEE = sumSignalNumEE / sumSignalDenomEE;
  double effSignalEEErr = sqrt( sumSignalDenomEEErr2 
				* effSignalEE*(1-effSignalEE)
				/(sumSignalDenomEE*sumSignalDenomEE) );
  std::cout << "Signal endcap efficiency: " << effSignalEE*100 << " +- " << effSignalEEErr*100 << " %%" << std::endl; //same as above 
	 

  // Compute background efficiencies
  double effBackEB = sumBackNumEB / sumBackDenomEB;
  double effBackEBErr = sqrt( sumBackDenomEBErr2 
				* effBackEB*(1-effBackEB)
				/(sumBackDenomEB*sumBackDenomEB) );
  std::cout << "Background barrel efficiency: " << effBackEB*100 << " +- " << effBackEBErr*100 << " %%" << std::endl;  
	 

  double effBackEE = sumBackNumEE / sumBackDenomEE;
  double effBackEEErr = sqrt( sumBackDenomEEErr2 
				* effBackEE*(1-effBackEE)
				/(sumBackDenomEE*sumBackDenomEE) );
  std::cout << "Background endcap efficiency: " << effBackEE*100 << " +- " << effBackEEErr*100 << " %%" << std::endl; 
	

  std::cout << "" << std::endl;
  std::cout << " signal photons: " << nSigBarrel << "    weighted:   " << nSigBarrelWeighted << std::endl;
  std::cout << " backgr photons: " << nBgBarrel << "    weighted:   " << nBgBarrelWeighted << std::endl;
  
  std::cout << "Integral of fake rate " << FRvsET->Integral() << std::endl;

}


// Effective area corrections
const int nEtaBinsEA = 7;
const float etaBinLimits[nEtaBinsEA+1] = {
  0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};

//Old Effective areas
/* const float areaPhotons[nEtaBinsEA] = {
  0.0894, 0.0750, 0.0423, 0.0561, 0.0882, 0.1144, 0.1684
};
const float areaNeutralHadrons[nEtaBinsEA] = {
  0.049, 0.0108, 0.0019, 0.0037, 0.0062, 0.0130, 0.1699
};
const float areaChargedHadrons[nEtaBinsEA] = {
  0.0089, 0.0062, 0.0086, 0.0041, 0.0113, 0.0085, 0.0039
}; */

//New effective areas
const float areaPhotons[nEtaBinsEA] = {
  0.0982, 0.0857, 0.0484, 0.0668, 0.0868, 0.0982, 0.1337
};
const float areaNeutralHadrons[nEtaBinsEA] = {
  0.0126, 0.0237, 0.0, 0.0, 0.0, 0.0, 0.0769
};
const float areaChargedHadrons[nEtaBinsEA] = {
  0.0080, 0.0079, 0.0080, 0.0048, 0.0029, 0.0036, 0.0016
};

float computeIso(IsoType isoType, float isoVal, float eta, float rho){

  const float *areas = 0;
  if( isoType == CH_ISO )
    areas = areaChargedHadrons;
  else if (isoType == NH_ISO)
    areas = areaNeutralHadrons;
  else if (isoType == PH_ISO)
    areas = areaPhotons;
  else
    assert(0);
  // Compute isolation with effective area correction for PU
  // Find eta bin first. If eta>2.5, the last eta bin is used.
  int etaBin = 0; 
  while ( etaBin < nEtaBinsEA-1 
	  && fabs( eta ) > etaBinLimits[nEtaBinsEA+1] )
    { ++etaBin; };

  float isoValWithEA =  std::max( (float)0.0, isoVal - rho * areas[etaBin] );
  return isoValWithEA;				  
}


// // ID cuts and code - almost latest version
// const float hOverECut[2][nWP] = 
//   { { 0.032, 0.020, 0.012 },
//     { 0.023, 0.011, 0.011 } };
// 
// const float sieieCut[2][nWP] = 
//   { {0.0100, 0.0099, 0.0098},
//     {0.0270, 0.0269, 0.0264} };
// 
// const float chIsoCut[2][nWP] = 
//   { {2.94, 2.62, 1.91},
//     {3.07, 1.40, 1.26} };
// 
// const float nhIso_A[2][nWP] = 
//   { {3.16, 2.69, 2.55},
//     {17.16, 4.92, 2.71} };
// 
// const float nhIso_B[2][nWP] = 
//   { {0.0023, 0.0023, 0.0023},
//     {0.0116, 0.0116, 0.0116} };
// 
// const float phIso_A[2][nWP] = 
//   { {4.43, 1.35, 1.29},
//     {2.11, 2.11, 1.91} };
// 
// const float phIso_B[2][nWP] = 
//   { {0.0004, 0.0004, 0.0004},
//     {0.0037, 0.0037, 0.0037} };

// ID cuts and code - newest
const float hOverECut[2][nWP] =
  { { 0.028, 0.012, 0.010 },
    { 0.093, 0.023, 0.015 } };

const float sieieCut[2][nWP] =
  { {0.0107, 0.0100, 0.0100},
    {0.0272, 0.0267, 0.0265} };

const float chIsoCut[2][nWP] =
  { {2.67, 1.79, 1.66},
    {1.79, 1.09, 1.04} };

const float nhIso_A[2][nWP] =
  { {7.23, 0.16, 0.14},
    {8.89, 4.31, 3.89} };

//----------add new ones of these (exponentials) -----------
const float nhIso_B[2][nWP] =
  { {0.0172, 0.0172, 0.0172},
    {0.01725, 0.0172, 0.0172} };

const float phIso_A[2][nWP] =
  { {2.11, 1.90, 1.40},
    {3.09, 1.90, 1.40} };

const float phIso_B[2][nWP] =
  { {0.0014, 0.0014, 0.0014},
    {0.0091, 0.0091, 0.0091} };

bool passWorkingPoint(WpType iwp, bool isBarrel, float pt,
		      float hOverE, float full5x5_sigmaIetaIeta, 
		      float chIso, float nhIso, float phIso){

  int ieta = 0;
  if( !isBarrel )
    ieta = 1;

/*   bool result = 1
    && hOverE < hOverECut[ieta][iwp]
    && full5x5_sigmaIetaIeta > 0 // in case miniAOD sets this to zero due to pre-selection of storage
    && full5x5_sigmaIetaIeta < sieieCut[ieta][iwp]
    && chIso < chIsoCut[ieta][iwp]
    && nhIso < nhIso_A[ieta][iwp] + pt * nhIso_B[ieta][iwp]
    && phIso < phIso_A[ieta][iwp] + pt * phIso_B[ieta][iwp] ; */
    
    bool result = 1
    && hOverE < hOverECut[ieta][iwp]
    && full5x5_sigmaIetaIeta > 0 // in case miniAOD sets this to zero due to pre-selection of storage
    && full5x5_sigmaIetaIeta < sieieCut[ieta][iwp]
    && chIso < chIsoCut[ieta][iwp]
    && nhIso < nhIso_A[ieta][iwp] + pt * nhIso_B[ieta][iwp]
    && phIso < phIso_A[ieta][iwp] + pt * phIso_B[ieta][iwp] ;  
  
  return result;
}


//
// MC truth matching
//

// ggNtuple arrays needed:
//   mcStatus, mcPID, mcMomPID, mcEta, mcPhi 
bool isMatched(float pEta, float pPhi,
	       std::vector<int> *mcPID,
	       std::vector<int> *mcMomPID,
	       std::vector<int> *mcStatus,
	       std::vector<float> *mcEta,
	       std::vector<float> *mcPhi){

  bool isMatched = false;
  // double genPt = -1;   
  if(verbose) std::cout << "Check match for photon eta= " << pEta << ", phi= " << pPhi << std::endl;
  for(unsigned int imc = 0; imc < (*mcPID).size(); imc++){//uint? imc...

    if(verbose) std::cout << "   Check next particle: pid= " << (*mcPID)[imc] << "  status= " << (*mcStatus)[imc] << "   mom= " << (*mcMomPID)[imc] << "  (eta,phi)=(" << (*mcEta)[imc] << 
    			", " << (*mcPhi)[imc] << ")" << std::endl;
		         
    double msts = (*mcStatus)[imc];
    if((msts != 1)||((*mcPID)[imc] != 11)) continue; //22
    if(verbose) std::cout << "      passed pid and status" << std::endl;

/*    if((fabs((*mcMomPID)[imc]) !=21) 
       &&(fabs((*mcMomPID)[imc]) !=1)
       &&(fabs((*mcMomPID)[imc]) !=2)
       &&(fabs((*mcMomPID)[imc]) !=3)
       &&(fabs((*mcMomPID)[imc]) !=4)
       &&(fabs((*mcMomPID)[imc]) !=5)
       &&(fabs((*mcMomPID)[imc]) !=6))continue;    
    if(verbose) std::cout << "       passed mother" << std::endl;*/

    if(fabs((*mcMomPID)[imc]) != 23) continue;
    if(verbose) std::cout << "       passed mother" << std::endl;
    
    double meta = (*mcEta)[imc];
    double mphi = (*mcPhi)[imc];
    
    TVector3 mcphoton;
    TVector3 mcElectron;
    TVector3 recoPHOTOn;
  //  mcphoton.SetPtEtaPhi(1.0,meta,mphi);
    mcElectron.SetPtEtaPhi(1.0,meta,mphi);
    recoPHOTOn.SetPtEtaPhi(1.0,pEta,pPhi);
    
//    double DR = mcphoton.DrEtaPhi(recoPHOTOn);
    double DR = mcElectron.DrEtaPhi(recoPHOTOn);
    if(verbose) std::cout << "        dR= " << DR << std::endl;
    if(DR < 0.1 ){
      isMatched = true;
      if(verbose) std::cout << "         PASSED ALL" << std::endl;
      // genPt = (*mcPt)[imc];
      break;
    }  
  }

  return isMatched;
}

/*bool passesZmass(vector< vector<float> > photonCollection){

	bool passesZmass = false;

	TLorentzVector electron1;
	TLorentzVector electron2;
	
	electron1.SetPtEtaPhiE(photonCollection[0][0],photonCollection[0][1],photonCollection[0][2],photonCollection[0][0]);
	electron2.SetPtEtaPhiE(photonCollection[1][0],photonCollection[1][1],photonCollection[1][2],photonCollection[1][0]);
	
	float mass_ll = (electron1 + electron2).M();
	std::cout << "Inv Mass = " << mass_ll << std::endl;
		
	if(mass_ll > 60 && mass_ll < 120)
		passesZmass = true;
		
		return passesZmass; 
		//TH1F * Histo_mll = new TH1F("Histo_mll","mll",60,60,120);

}*/

void set_plot_style() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}            

TText* doPreliminary(double x_pos,double y_pos){

          ostringstream stream;
          stream  << "CMS Preliminary";                                                                             

          TLatex* text = new TLatex(x_pos, y_pos, stream.str().c_str());
          text->SetNDC(true);
          text->SetTextFont(62);
          text->SetTextSize(0.03);  // for thesis

          return text;
}   

void doPlot(TH1D* plot, Color_t fillColour){

TString name = string(plot->GetName());

  TCanvas* can = new TCanvas("Plot","Plot",635, 600); 
      
      plot->SetFillColor(fillColour); 
      plot->SetLineColor(kBlack);  
      plot->Rebin(5);
      plot->SetMaximum(plot->GetBinContent(plot->GetMaximumBin())*1.3);
      plot->Draw();

      if(addHashErrors){
            TH1D* hashErrors = doHashErrors(plot);
            hashErrors->Draw("same");
      }

      TText* prelim = doPreliminary(0.14, 0.84);
      prelim->Draw("same");

      TLegend* leg = doLegend(plot);
      leg->Draw("same");

  if(doLogPlot){    
      can->SetLogy();
      can->SaveAs("Plots/Log/"+name+".pdf");
  }else{   
      can->SaveAs("Plots/"+name+".pdf");
  }    

  delete can;
  delete prelim;
  delete leg;

}      

TLegend* doLegend(TH1D* plot){

      TLegend* leg = new TLegend(0.65,0.82,0.85,0.88);
      leg->SetTextSize(0.04);
      leg->SetBorderSize(0);
      leg->SetFillColor(10);
      leg->AddEntry(plot , "DY 20bx25", "f");
      return leg;
}           

TH1D* doHashErrors(TH1D* plot){
  TH1D * hashErrors = plot;

  hashErrors->SetFillColor(kBlack);
  hashErrors->SetFillStyle(3354);
  hashErrors->SetMarkerSize(0.);
  hashErrors->SetStats(0);

  return hashErrors;
}
