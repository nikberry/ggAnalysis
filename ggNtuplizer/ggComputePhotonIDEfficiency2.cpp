#include "TVector3.h"
#include "TH2D.h"
#include "TCut.h"
#include <math.h>
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TH1D.h"
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
void doPlot(TH1D* plot);
TLegend* doLegend(TH1D* plot);

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
  
  TH1D *phoEtCheckHist_conv = new TH1D("phoEtCheckHistConv","",185, 15, 200);	
        phoEtCheckHist_conv->GetXaxis()->SetTitle("E_{T} [GeV]");
        phoEtCheckHist_conv->GetYaxis()->SetTitle("Number of events");
        phoEtCheckHist_conv->SetFillColor(kAzure+8);
        phoEtCheckHist_conv->SetLineColor(kBlack);      
  TH1D *phoSCEtaCheckHist_conv = new TH1D("phoSCEtaCheckHistConv","",200, -5, 5);
        phoSCEtaCheckHist_conv->GetXaxis()->SetTitle("#eta_{SC} [GeV]");
        phoSCEtaCheckHist_conv->GetYaxis()->SetTitle("Number of events");
        phoSCEtaCheckHist_conv->SetFillColor(kAzure+8);
        phoSCEtaCheckHist_conv->SetLineColor(kBlack);
  TH2D *phoTotCheck_conv = new TH2D("phoTotCheckConv","", 200, -5, 5, 185, 15, 200); 
        phoTotCheck_conv->GetXaxis()->SetTitle("#eta_{SC}");
        phoTotCheck_conv->GetYaxis()->SetTitle("E_{T} [GeV]"); 
  
  //barrel histos
  TH1D* phoETSigEB = new TH1D("phoETSigEB", "", 185, 15, 200);
  	    phoETSigEB->GetXaxis()->SetTitle("E_{T} [GeV]");
	      phoETSigEB->GetYaxis()->SetTitle("Number of events");

  TH1D* phoETSigEB_noconv = new TH1D("phoETSigEB_noconv", "", 185, 15, 200); 
    	  phoETSigEB_noconv->GetXaxis()->SetTitle("E_{T} [GeV]");
	      phoETSigEB_noconv->GetYaxis()->SetTitle("Number of events");
        phoETSigEB_noconv->SetFillColor(kAzure+8);
        phoETSigEB_noconv->SetLineColor(kBlack);
  
  //endcap histos
  TH1D* phoETSigEE = new TH1D("phoETSigEE", "", 185, 15, 200);
    	  phoETSigEE->GetXaxis()->SetTitle("E_{T} [GeV]");
	      phoETSigEE->GetYaxis()->SetTitle("Number of events");
        phoETSigEE->SetFillColor(kTeal+8);
        phoETSigEE->SetLineColor(kBlack);
  TH1D* phoETSigEE_noconv = new TH1D("phoETSigEE_noconv", "", 185, 15, 200);  
    	  phoETSigEE_noconv->GetXaxis()->SetTitle("E_{T} [GeV]");
	      phoETSigEE_noconv->GetYaxis()->SetTitle("Number of events");
        phoETSigEE_noconv->SetFillColor(kTeal+8);
        phoETSigEE_noconv->SetLineColor(kBlack);
  
  // Event-level variables:
  int nPho;
  float rho;
  
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
  TBranch *b_phoEt = 0;
  TBranch *b_phoSCEta = 0;
  TBranch *b_phoPhi = 0;
  TBranch *b_phoSigmaIEtaIEta_2012 = 0;
  TBranch *b_phoHoverE = 0;
  TBranch *b_phoPFPhoIso = 0;
  TBranch *b_phoPFChIso = 0;
  TBranch *b_phoPFNeuIso = 0;
  TBranch *b_phohasPixelSeed = 0;
  TBranch * b_phoEleVeto = 0;
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
      if( !( phoEleVeto->at(ipho) == 1) ) continue;
//      if( !( phohasPixelSeed->at(ipho) == 0 ) ) continue;
     // if( nPho < 2 ) continue; 
      
  //    photon.push_back(phoEt->at(ipho));
  //    photon.push_back(phoSCEta->at(ipho));
  //    photon.push_back(phoPhi->at(ipho));
  //    photonCollection.push_back(photon);
   
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
	 phoEtCheckHist_conv->Fill(phoEt->at(ipho), weight);
	 phoSCEtaCheckHist_conv->Fill(phoSCEta->at(ipho), weight);
	 phoTotCheck_conv->Fill(phoSCEta->at(ipho), phoEt->at(ipho), weight);

	  	

	 if( isBarrel ) {
	  nSigBarrel++;
	  nSigBarrelWeighted += weight;
	  sumSignalDenomEB += weight;
	  sumSignalDenomEBErr2 += weight*weight;

	  if( pass ) {
	    sumSignalNumEB += weight;
	    sumSignalNumEBErr2 += weight*weight;
	    
	    phoETSigEB->Fill(phoEt->at(ipho), weight);

	  }
	 } else if( isEndcap ){
	  sumSignalDenomEE += weight;
	  sumSignalDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	    sumSignalNumEE += weight;
	    sumSignalNumEEErr2 += weight*weight;
	    
	    phoETSigEE->Fill(phoEt->at(ipho), weight);

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

	  	

	if( isBarrel ) {
	  //nSigBarrel++;
	  //nSigBarrelWeighted += weight;
	  //sumSignalDenomEB += weight;
	  //sumSignalDenomEBErr2 += weight*weight;

	  if( pass ) {
	    //sumSignalNumEB += weight;
	    //sumSignalNumEBErr2 += weight*weight;
	    
	    phoETSigEB_noconv->Fill(phoEt->at(ipho), weight);

	  }
	} else if( isEndcap ){
	  //sumSignalDenomEE += weight;
	  //sumSignalDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	   // sumSignalNumEE += weight;
	   // sumSignalNumEEErr2 += weight*weight;
	    
	    phoETSigEE_noconv->Fill(phoEt->at(ipho), weight);

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
    

}// end loop over events  

TH1D* FRSigEB = (TH1D*)phoETSigEB->Clone("FR");
      FRSigEB->SetFillColor(kAzure+8);
  FRSigEB->Divide(phoETSigEB_noconv);

TH1D* FRSigEE = (TH1D*)phoETSigEE->Clone("FR");
      FRSigEB->SetFillColor(kTeal+8);
  FRSigEE->Divide(phoETSigEE_noconv);  

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

// Plotting section

doPlot(phoETSigEB);


TCanvas* signalEBnoconv = new TCanvas("signalEBnoconv");
  phoETSigEB_noconv->Draw();
signalEBnoconv->SaveAs("Plots/phoETSigEB_noconv.pdf");

TCanvas* ratioSigEB = new TCanvas("ratioSigEB");
  FRSigEB->Draw();
ratioSigEB->SaveAs("Plots/ratioSigEB.pdf");

TCanvas* signalEE = new TCanvas("signalEE");
  phoETSigEE->Draw();
signalEE->SaveAs("Plots/phoETSigEE.pdf");

TCanvas* signalEEnoconv = new TCanvas("signalEEnoconv");
  phoETSigEE_noconv->Draw();
signalEEnoconv->SaveAs("Plots/phoETSigEE_noconv.pdf");

TCanvas* ratioSigEE = new TCanvas("ratioSigEE");
  FRSigEE->Draw();
ratioSigEE->SaveAs("Plots/ratioSigEE.pdf");


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
  
  TCanvas *c2 = new TCanvas("c2","c2",10,10,900, 600);
  c2->Divide(3,1);
  c2->cd(1);
  phoEtCheckHist_conv->Draw();
  c2->cd(2);
  phoSCEtaCheckHist_conv->Draw();
  c2->cd(3);
  phoTotCheck_conv->Draw("colz");

  c2->SaveAs("Plots/phoCheckHists.pdf");


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
      if(verbose) std::cout << "         PASSE ALL" << std::endl;
      // genPt = (*mcPt)[imc];
      break;
    }  
  }

  return isMatched;
}

/*bool passesZmass(vector< vector<float> > photonCollection){

	bool passesZmass = false;

	TLorentzVector photon1;
	TLorentzVector photon2;
	
	photon1.SetPtEtaPhiE(photonCollection[0][0],photonCollection[0][1],photonCollection[0][2],photonCollection[0][0]);
	photon2.SetPtEtaPhiE(photonCollection[1][0],photonCollection[1][1],photonCollection[1][2],photonCollection[1][0]);
	
	float mass_ll = (photon1 + photon2).M();
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

void doPlot(TH1D* plot){

TString name = string(plot->GetName());

  TCanvas* can = new TCanvas("Plot","Plot",635, 600); 
      plot->SetMaximum(plot->GetBinContent(plot->GetMaximumBin())*1.3);
      plot->SetFillColor(kAzure+8);
      //plot->SetFillColor(kTeal+8); 
      plot->SetLineColor(kBlack);  
      plot->Draw();

      TText* prelim = doPreliminary(0.14, 0.82);
      prelim->Draw("same");

      TLegend* leg = doLegend(plot);
      leg->Draw("same");

  if(doLogPlot){    
      can->SetLogy();
      can->SaveAs("Plots/Log/"+name+".pdf");
  }else{   
      can->SaveAs("Plots/"+name+".pdf");
  }    

}      

TLegend* doLegend(TH1D* plot){

      TLegend* leg = new TLegend(0.70,0.48,0.90,0.88);
      leg->SetTextSize(0.04);
      leg->SetBorderSize(0);
      leg->SetFillColor(10);
      leg->AddEntry(plot , "2012 data", "lpe");
      return leg;
}           

