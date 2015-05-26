#include "TVector3.h"
#include "TH2D.h"
#include "TCut.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TH1D.h"
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

const int nWP = 3;
enum WpType { WP_LOOSE = 0,
	      WP_MEDIUM, 
	      WP_TIGHT};
const TString wpName[nWP] = 
  {"Loose", "Medium", "Tight"};

enum IsoType {CH_ISO, NH_ISO, PH_ISO};

float computeIso(IsoType isoType, float isoVal, float eta, float rho); 

// Use here one of the WpType values
const WpType wp = WP_TIGHT;

const TString treename = "ggNtuplizer/EventTree";
// const TString fname1 = "~/workspace/ntuples/ggNtuple_GJets_HT-100to200_singleFile.root";
//const TString fname1 = "/tmp/rslu/job_phys14_gjet_pt40_20bx25.root";
//const TString fname1 = "/media/nik/NIKBERRYHD/job_spring14_DYJets_20bx25.root";
const TString fname1 = "/nfs/data/eepgnnb/EGamma/job_spring14_DYJets_20bx25.root";

bool verbose = false;
bool smallEventCount = false;

const int nEtaBins = 2;

const float ptMin = 30;
const float ptMax = 200;
const float barrelEtaLimit = 1.479;
const float endcapUpperEtaLimit = 2.5;
const float endcapLowerEtaLimit = 1.479;

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

void ggComputePhotonIDEfficiency()
{

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 
//  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  set_plot_style();
  gROOT->SetBatch();
  //
  // Find the tree
  //
  TFile *file1 = new TFile(fname1);
  if( !file1 )
    assert(0);
  TTree *tree = (TTree*)file1->Get(treename);
  if( !tree )
    assert(0);
  printf("Found the tree\n"); fflush(stdout);


//Histos for conv safe veto
  TH1D *phoEtCheckHist_conv = new TH1D("phoEtCheckHistConv","",185, 15, 200);
  TH1D *phoSCEtaCheckHist_conv = new TH1D("phoSCEtaCheckHistConv","",200, -5, 5);
  TH2D *phoTotCheck_conv = new TH2D("phoTotCheckConv","", 200, -5, 5, 185, 15, 200);

  TH2D* FRvsEta =  new TH2D("FRvsEta", "", 100, -5, 5, 200, 0, 5); 
//  TH2D* FRvsPt = new TH2D("FRvsPt", "", 200, 0, 200, 200, 0, 200);
//  TH2D* FRvsPU = new TH2D("FRvsPU", "", 200, 0, 200, 200, 0, 200); 

  TH1D* phoET_Sig_conv = new TH1D("phoETSigConv", "", 200, 0, 200);
  	phoET_Sig_conv->GetXaxis()->SetTitle("E_T");
  	phoET_Sig_conv->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_Sig_conv = new TH1D("phoEtaSigConv", "", 100, -5, 5);
  	phoEta_Sig_conv->GetXaxis()->SetTitle("#eta");
  	phoEta_Sig_conv->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_Sig_conv = new TH1D("phoPhiSigConv", "", 100, -4, 4);
        phoPhi_Sig_conv->GetXaxis()->SetTitle("#phi");
        phoPhi_Sig_conv->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_Bkgd_conv = new TH1D("phoEtaBkgdConv", "", 100, -5, 5);
//  TH1D* h_nPho = new TH1D("nPho", "", 100, 0, 10);

//Histos for pixel veto
  TH1D *phoEtCheckHist_pix = new TH1D("phoEtCheckHistPix","",185, 15, 200);
  TH1D *phoSCEtaCheckHist_pix = new TH1D("phoSCEtaCheckHistPix","",200, -5, 5);
  TH2D *phoTotCheck_pix = new TH2D("phoTotCheckPix","", 200, -5, 5, 185, 15, 200);
  TH1D* phoET_Sig_pix = new TH1D("phoETSigPix", "", 200, 0, 200);
        phoET_Sig_pix->GetXaxis()->SetTitle("E_{T}");
        phoET_Sig_pix->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_Sig_pix = new TH1D("phoEtaSigPix", "", 100, -5, 5);
        phoEta_Sig_pix->GetXaxis()->SetTitle("#eta");
        phoEta_Sig_pix->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_Sig_pix = new TH1D("phoPhiSigPix", "", 100, -4, 4);
        phoPhi_Sig_pix->GetXaxis()->SetTitle("#phi");
        phoPhi_Sig_pix->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_Bkgd_pix = new TH1D("phoEtaBkgdPix", "", 100, -5, 5);

//Histos for inverted cut
  TH1D* phoET_SigInv_conv = new TH1D("phoETSigInvConv", "", 200, 0, 200);
        phoET_SigInv_conv->GetXaxis()->SetTitle("E_{T}");
        phoET_SigInv_conv->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_SigInv_conv = new TH1D("phoEtaSigInvConv", "", 100, -5, 5);
        phoEta_SigInv_conv->GetXaxis()->SetTitle("#eta");
        phoEta_SigInv_conv->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_SigInv_conv = new TH1D("phoPhiSigInvConv", "", 100, -4, 4);
        phoPhi_SigInv_conv->GetXaxis()->SetTitle("#phi");
        phoPhi_SigInv_conv->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_BkgdInv_conv = new TH1D("phoEtaBkgdInvConv", "", 100, -5, 5); 


  TH1D* phoET_SigInv_pix = new TH1D("phoETSigInvPix", "", 200, 0, 200);
        phoET_SigInv_pix->GetXaxis()->SetTitle("E_{T}");
        phoET_SigInv_pix->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_SigInv_pix = new TH1D("phoEtaSigInvPix", "", 100, -5, 5);
        phoEta_SigInv_pix->GetXaxis()->SetTitle("#eta");
        phoEta_SigInv_pix->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_SigInv_pix = new TH1D("phoPhiSigInvPix", "", 100, -4, 4);
        phoPhi_SigInv_pix->GetXaxis()->SetTitle("#phi");
        phoPhi_SigInv_pix->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_BkgdInv_pix = new TH1D("phoEtaBkgdInvPix", "", 100, -5, 5);

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
  hSignal->GetYaxis()->SetTitle("E_{T}");
  TH2D *hBackground = new TH2D("hBackground","",200,-5,5,185,15,200);
  hBackground->GetXaxis()->SetTitle("#eta_{SC}");
  hBackground->GetYaxis()->SetTitle("E_{T}"); 
 
  UInt_t maxEvents = tree->GetEntries();
  if( smallEventCount )
    maxEvents = 1000;
  if(verbose)
    printf("Start loop over events for WEIGHTS, total events = %lld\n", 
           tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){
    
    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
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
//      if( !( phoEleVeto->at(ipho) == 1) ) continue;
 //     if( !( phohasPixelSeed->at(ipho) == 0 ) ) continue;
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

//TCanvas* canFRvsEta = new TCanvas("FRvsEta");
//	FRvsEta->Draw();
//canFRvsEta->SaveAs("Plots/FRvsEta.pdf");
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


  TH1D* phoET_Sig_conv_barrel = new TH1D("phoETSigConvBarrel", "", 200, 0, 200);
        phoET_Sig_conv_barrel->GetXaxis()->SetTitle("E_{T}");
        phoET_Sig_conv_barrel->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_Sig_conv_barrel = new TH1D("phoEtaSigConvBarrel", "", 100, -5, 5);
        phoEta_Sig_conv_barrel->GetXaxis()->SetTitle("#eta");
        phoEta_Sig_conv_barrel->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_Sig_conv_barrel = new TH1D("phoPhiSigConvBarrel", "", 100, -4, 4);
        phoPhi_Sig_conv_barrel->GetXaxis()->SetTitle("#phi");
        phoPhi_Sig_conv_barrel->GetYaxis()->SetTitle("Events");
//  TH1D* phoEta_Bkgd_conv = new TH1D("phoEtaBkgdConv", "", 100, -5, 5);

  TH1D* phoET_SigFull_conv_barrel = new TH1D("phoETSigFullConvEB", "", 200, 0, 200);
        phoET_SigFull_conv_barrel->GetXaxis()->SetTitle("E_{T}");
        phoET_SigFull_conv_barrel->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_SigFull_conv_barrel = new TH1D("phoEtaSigFullConvEB", "", 100, -5, 5);
        phoEta_SigFull_conv_barrel->GetXaxis()->SetTitle("#eta");
        phoEta_SigFull_conv_barrel->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_SigFull_conv_barrel = new TH1D("phoPhiSigFullConvEB", "", 100, -4, 4);
        phoPhi_SigFull_conv_barrel->GetXaxis()->SetTitle("#phi");
        phoPhi_SigFull_conv_barrel->GetYaxis()->SetTitle("Events");
//  TH1D* phoEta_BkgdFull_conv_barrel = new TH1D("phoEtaBkgdConv", "", 100, -5, 5);

  TH1D* phoET_SigFull_conv_endcap = new TH1D("phoETSigFullConvEE", "", 200, 0, 200);
        phoET_SigFull_conv_endcap->GetXaxis()->SetTitle("E_{T}");
        phoET_SigFull_conv_endcap->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_SigFull_conv_endcap = new TH1D("phoEtaSigFullConvEE", "", 100, -5, 5);
        phoEta_SigFull_conv_endcap->GetXaxis()->SetTitle("#eta");
        phoEta_SigFull_conv_endcap->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_SigFull_conv_endcap = new TH1D("phoPhiSigFullConvEE", "", 100, -4, 4);
        phoPhi_SigFull_conv_endcap->GetXaxis()->SetTitle("#phi");
        phoPhi_SigFull_conv_endcap->GetYaxis()->SetTitle("Events");


  // 
  // The second loop over events is for efficiency computation
  //
  if(verbose)
    printf("Start loop over events for EFFICIENCY caluclation, total events = %lld\n", 
           tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);

    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of photons %u\n", ievent, nPho);

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
 //     if( !( phohasPixelSeed->at(ipho) == 0 ) ) continue;

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
      	printf("Error! Zero! pt=%f, eta=%f\n", phoEt->at(ipho), phoSCEta->at(ipho));
      }else{
      	weight = 1./binContent;
      }
      

      if( isTrue ) {
	phoEtCheckHist_conv->Fill(phoEt->at(ipho), weight);
	phoSCEtaCheckHist_conv->Fill(phoSCEta->at(ipho), weight);
	phoTotCheck_conv->Fill(phoSCEta->at(ipho), phoEt->at(ipho), weight);
	phoET_Sig_conv->Fill(phoEt->at(ipho), weight);
	phoEta_Sig_conv->Fill(phoSCEta->at(ipho), weight);
	phoPhi_Sig_conv->Fill(phoPhi->at(ipho), weight);
//	FRvsEta->Fill(phoSCEta->at(ipho), phoSCEta->at(ipho)/(phoSCEta->at(ipho)+phoSCEta->at(ipho)), weight);
//	FRvsET->Fill();
//	FRvsPU->Fill();
	  	

	if( isBarrel ) {
	  nSigBarrel++;
	  nSigBarrelWeighted += weight;
	  sumSignalDenomEB += weight;
	  sumSignalDenomEBErr2 += weight*weight;

	  phoET_Sig_conv_barrel->Fill(phoEt->at(ipho), weight);
	  if( pass ) {
	    sumSignalNumEB += weight;
	    sumSignalNumEBErr2 += weight*weight;

	    phoET_SigFull_conv_barrel->Fill(phoEt->at(ipho), weight);
	    phoEta_SigFull_conv_barrel->Fill(phoSCEta->at(ipho), weight);
	    phoPhi_SigFull_conv_barrel->Fill(phoPhi->at(ipho), weight);
	  }
	} else if( isEndcap ){
	  sumSignalDenomEE += weight;
	  sumSignalDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	    sumSignalNumEE += weight;
	    sumSignalNumEEErr2 += weight*weight;

            phoET_SigFull_conv_endcap->Fill(phoEt->at(ipho), weight);
            phoEta_SigFull_conv_endcap->Fill(phoSCEta->at(ipho), weight);
            phoPhi_SigFull_conv_endcap->Fill(phoPhi->at(ipho), weight);

	  }

	}// end barrel / endcap
      } // end if signal
	

      if( !isTrue ) {

	phoEta_Bkgd_conv->Fill(phoSCEta->at(ipho), weight);
 	
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

  TCanvas* canoe = new TCanvas("canoe");

  phoET_SigFull_conv_barrel->SetLineColor(kGreen);
  phoET_SigFull_conv_barrel->SetFillColor(kGreen);
  phoET_SigFull_conv_barrel->Draw("e");

  canoe->SaveAs("Plots/phoETSigFullConvBarrel.pdf"); 

  TCanvas* canoe2 = new TCanvas("canoe2");

  phoEta_SigFull_conv_barrel->SetLineColor(kGreen);
  phoEta_SigFull_conv_barrel->SetFillColor(kGreen);
  phoEta_SigFull_conv_barrel->Draw("e");

  canoe2->SaveAs("Plots/phoEtaSigFullConvBarrel.pdf");

  TCanvas* canoe3 = new TCanvas("canoe3");

  phoPhi_SigFull_conv_barrel->SetLineColor(kGreen);
  phoPhi_SigFull_conv_barrel->SetFillColor(kGreen);
  phoPhi_SigFull_conv_barrel->Draw("e");

  canoe3->SaveAs("Plots/phoPhiSigFullConvBarrel.pdf");

  TCanvas* canoe4 = new TCanvas("canoe4");

  phoET_SigFull_conv_endcap->SetLineColor(kGreen);
  phoET_SigFull_conv_endcap->SetFillColor(kGreen);
  phoET_SigFull_conv_endcap->Draw("e");

  canoe4->SaveAs("Plots/phoETSigFullConvEndcap.pdf");

  TCanvas* canoe5 = new TCanvas("canoe5");

  phoEta_SigFull_conv_endcap->SetLineColor(kGreen);
  phoEta_SigFull_conv_endcap->SetFillColor(kGreen);
  phoEta_SigFull_conv_endcap->Draw("e");

  canoe5->SaveAs("Plots/phoEtaSigFullConvEndcap.pdf");

  TCanvas* canoe6 = new TCanvas("canoe6");

  phoPhi_SigFull_conv_endcap->SetLineColor(kGreen);
  phoPhi_SigFull_conv_endcap->SetFillColor(kGreen);
  phoPhi_SigFull_conv_endcap->Draw("e");

  canoe6->SaveAs("Plots/phoPhiSigFullConvEndcap.pdf");

  TCanvas* stages =  new TCanvas("stages");
  stages->Divide(3,1);
  stages->cd(1);
  phoET_Sig_conv->Draw("e");
  stages->cd(2);
  phoET_Sig_conv_barrel->Draw("e");
  stages->cd(3);
  phoET_SigFull_conv_barrel->Draw("e");
  stages->SaveAs("Plots/selectionstages.pdf");

 
  // 
  // The third loop over events is for pixel veto plots  
  // //
  if(verbose)
    printf("Start loop over events for pixel veto plots, total events = %lld\n", 
           tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);

    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of photons %u\n", ievent, nPho);

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
 //     if( !( phoEleVeto->at(ipho) == 1) ) continue;
      if( !( phohasPixelSeed->at(ipho) == 0 ) ) continue;

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
      	printf("Error! Zero! pt=%f, eta=%f\n", phoEt->at(ipho), phoSCEta->at(ipho));
      }else{
      	weight = 1./binContent;
      }
      

      if( isTrue ) {
	phoEtCheckHist_pix->Fill(phoEt->at(ipho), weight);
	phoSCEtaCheckHist_pix->Fill(phoSCEta->at(ipho), weight);
	phoTotCheck_pix->Fill(phoSCEta->at(ipho), phoEt->at(ipho), weight);
	phoET_Sig_pix->Fill(phoEt->at(ipho), weight);
	phoEta_Sig_pix->Fill(phoSCEta->at(ipho), weight);
	phoPhi_Sig_pix->Fill(phoPhi->at(ipho), weight);
	  	

	if( isBarrel ) {
	  nSigBarrel++;
	  nSigBarrelWeighted += weight;
	  sumSignalDenomEB += weight;
	  sumSignalDenomEBErr2 += weight*weight;
	  if( pass ) {
	    sumSignalNumEB += weight;
	    sumSignalNumEBErr2 += weight*weight;
	  }
	} else if( isEndcap ){
	  sumSignalDenomEE += weight;
	  sumSignalDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	    sumSignalNumEE += weight;
	    sumSignalNumEEErr2 += weight*weight;
	  }

	}// end barrel / endcap
      } // end if signal
	

      if( !isTrue ) {

	phoEta_Bkgd_pix->Fill(phoSCEta->at(ipho), weight);
 	
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

  TH1D* phoET_SigInv_conv_barrel = new TH1D("phoETSigInvConv", "", 200, 0, 200);
        phoET_SigInv_conv_barrel->GetXaxis()->SetTitle("E_{T}");
        phoET_SigInv_conv_barrel->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_SigInv_conv_barrel = new TH1D("phoEtaSigInvConv", "", 100, -5, 5);
        phoEta_SigInv_conv_barrel->GetXaxis()->SetTitle("#eta");
        phoEta_SigInv_conv_barrel->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_SigInv_conv_barrel = new TH1D("phoPhiSigInvConv", "", 100, -4, 4);
        phoPhi_SigInv_conv_barrel->GetXaxis()->SetTitle("#phi");
        phoPhi_SigInv_conv_barrel->GetYaxis()->SetTitle("Events");

  TH1D* phoET_SigInvFull_conv_barrel = new TH1D("phoETSigInvConv", "", 200, 0, 200);
        phoET_SigInvFull_conv_barrel->GetXaxis()->SetTitle("E_{T}");
        phoET_SigInvFull_conv_barrel->GetYaxis()->SetTitle("Events");
  TH1D* phoEta_SigInvFull_conv_barrel = new TH1D("phoEtaSigInvConv", "", 100, -5, 5);
        phoEta_SigInvFull_conv_barrel->GetXaxis()->SetTitle("#eta");
        phoEta_SigInvFull_conv_barrel->GetYaxis()->SetTitle("Events");
  TH1D* phoPhi_SigInvFull_conv_barrel = new TH1D("phoPhiSigInvConv", "", 100, -4, 4);
        phoPhi_SigInvFull_conv_barrel->GetXaxis()->SetTitle("#phi");
        phoPhi_SigInvFull_conv_barrel->GetYaxis()->SetTitle("Events");

 // 
  // The fourth loop over events is for inverted conv veto (control sample) and fakerate computation and TH2D filling
  //
  if(verbose)
    printf("Start loop over events for EFFICIENCY caluclation, total events = %lld\n", 
           tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);

    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of photons %u\n", ievent, nPho);

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
 //     if( !( phohasPixelSeed->at(ipho) == 0 ) ) continue;

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
      	printf("Error! Zero! pt=%f, eta=%f\n", phoEt->at(ipho), phoSCEta->at(ipho));
      }else{
      	weight = 1./binContent;
      }
      

//  fakeRate->Fill(phoEt->at(ipho), phoEt->at(ipho)/(phoEt->at(ipho)+pho), weight);

      if( isTrue ) {
	phoET_SigInv_conv->Fill(phoEt->at(ipho), weight);
	phoEta_SigInv_conv->Fill(phoSCEta->at(ipho), weight);
	phoPhi_SigInv_conv->Fill(phoPhi->at(ipho), weight);
	  	

	if( isBarrel ) {
	  nSigBarrel++;
	  nSigBarrelWeighted += weight;
	  sumSignalDenomEB += weight;
	  sumSignalDenomEBErr2 += weight*weight;

	  phoET_SigInv_conv_barrel->Fill(phoEt->at(ipho), weight);

	  if( pass ) {
	    sumSignalNumEB += weight;
	    sumSignalNumEBErr2 += weight*weight;

	    phoET_SigInvFull_conv_barrel->Fill(phoEt->at(ipho), weight);

	  }
	} else if(isEndcap){
	  sumSignalDenomEE += weight;
	  sumSignalDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	    sumSignalNumEE += weight;
	    sumSignalNumEEErr2 += weight*weight;
	  }

	}// end barrel / endcap
      } // end if signal
	

      if( !isTrue ) {

	//phoEta_BkgdInv_conv->Fill(phoSCEta->at(ipho), weight);
 	
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

  TCanvas* control = new TCanvas("control");
  control->Divide(3,1);
  control->cd(1);
  phoET_SigInv_conv->SetLineColor(kRed);
  phoET_SigInv_conv->SetFillColor(kRed);
  phoET_SigInv_conv->Draw("e");
  control->cd(2);
  phoET_SigInv_conv_barrel->SetLineColor(kRed);
  phoET_SigInv_conv_barrel->SetFillColor(kRed);
  phoET_SigInv_conv_barrel->Draw("e");
  control->cd(3);
  phoET_SigInvFull_conv_barrel->SetLineColor(kRed);
  phoET_SigInvFull_conv_barrel->SetFillColor(kRed);
  phoET_SigInvFull_conv_barrel->Draw("e");

  control->SaveAs("Plots/phoETSigInvConv.pdf");

  TCanvas* sigcontrolcomp = new TCanvas("sigcontrolcomp");
  phoET_Sig_conv->SetLineColor(kCyan+1);
  phoET_Sig_conv->Draw();
  phoET_SigInv_conv->SetLineColor(kRed);
  phoET_SigInv_conv->Draw("same");
  sigcontrolcomp->SaveAs("Plots/SigControlComparison.pdf");


  TH1D* ratio = (TH1D*)phoET_Sig_conv->Clone();
  ratio->Divide(phoET_SigInv_conv);
  printf("Number of events in ratio: %f", ratio->Integral()); 
  TCanvas* ratiocan = new TCanvas("ratiocan");
  ratio->Draw();
  ratiocan->SaveAs("Plots/ratio.pdf");

  TH1D* numFR = (TH1D*)phoET_Sig_conv->Clone();
  TH1D* denomFR = (TH1D*)phoET_Sig_conv->Clone();
  denomFR->Add(phoET_SigInv_conv);
  TH1D* f = (TH1D*)numFR->Clone();
  f->Divide(denomFR);
  f->SetLineColor(kMagenta+1);
  TCanvas* FR = new TCanvas("FR");
  f->Draw("e");
  FR->SaveAs("Plots/FR.pdf");

  printf("barrel signal a=%f   b=%f\n", sumSignalNumEB, sumSignalDenomEB);
  printf("barrel back a=%f    b=%f\n", sumBackNumEB, sumBackDenomEB);

  printf("\nEfficiencies for the working point %s\n", wpName[wp].Data());

  // Compute signal efficiencies
  double effSignalEB = sumSignalNumEB / sumSignalDenomEB;
  double effSignalEBErr = sqrt( sumSignalDenomEBErr2 
				* effSignalEB*(1-effSignalEB)
				/(sumSignalDenomEB*sumSignalDenomEB) );
  printf("Signal barrel efficiency: %5.1f +- %5.1f %%\n", 
	 effSignalEB*100, effSignalEBErr*100 );

  double effSignalEE = sumSignalNumEE / sumSignalDenomEE;
  double effSignalEEErr = sqrt( sumSignalDenomEEErr2 
				* effSignalEE*(1-effSignalEE)
				/(sumSignalDenomEE*sumSignalDenomEE) );
  printf("Signal endcap efficiency: %5.1f +- %5.1f %%\n", 
	 effSignalEE*100, effSignalEEErr*100 );

  // Compute background efficiencies
  double effBackEB = sumBackNumEB / sumBackDenomEB;
  double effBackEBErr = sqrt( sumBackDenomEBErr2 
				* effBackEB*(1-effBackEB)
				/(sumBackDenomEB*sumBackDenomEB) );
  printf("Background barrel efficiency: %5.1f +- %5.1f %%\n", 
	 effBackEB*100, effBackEBErr*100 );

  double effBackEE = sumBackNumEE / sumBackDenomEE;
  double effBackEEErr = sqrt( sumBackDenomEEErr2 
				* effBackEE*(1-effBackEE)
				/(sumBackDenomEE*sumBackDenomEE) );
  printf("Background endcap efficiency: %5.1f +- %5.1f %%\n", 
	 effBackEE*100, effBackEEErr*100 );

  printf("\n");
  printf(" signal photons: %d    weighted:   %.1f\n", nSigBarrel, nSigBarrelWeighted);
  printf(" backgr photons: %d    weighted:   %.1f\n", nBgBarrel, nBgBarrelWeighted);
  
  TCanvas *c2 = new TCanvas("c2","c2",10,10,900, 600);
  c2->Divide(3,1);
  c2->cd(1);
  phoEtCheckHist_conv->Draw();
  c2->cd(2);
  phoSCEtaCheckHist_conv->Draw();
  c2->cd(3);
  phoTotCheck_conv->Draw("colz");

  c2->SaveAs("Plots/phoCheckHists.pdf");
  c2->SaveAs("Plots/phoCheckHists.png");

//Try various ways even if it works.
TH1D* Num = (TH1D*)phoEta_Sig_conv->Clone();
TH1D* Denom = (TH1D*)phoEta_Bkgd_conv->Clone();
Denom->Add(Num);
Denom->SetLineColor(kRed);
Num->Divide(Denom);
Num->SetFillColor(kCyan+1);
Num->SetLineColor(kBlack);
printf("Number of events: %f \n", Num->Integral());

  TCanvas *c3 = new TCanvas("c3","c3",10,10,900, 600);
 	c3->Divide(3,1);
	c3->cd(1);
		phoEta_Sig_conv->Draw(); 
	c3->cd(2);
		phoEta_Bkgd_conv->Draw();
	c3->cd(3);
		phoEta_Bkgd_conv->Draw();
		Denom->Draw("same");
  c3->SaveAs("Plots/tests.pdf");

  TCanvas* c4 = new TCanvas("c4");
  	Num->Draw();
  c4->SaveAs("Plots/fakeRateEta.pdf");

  TCanvas* c5 = new TCanvas("c5");
  	FRvsEta->Draw("colz");	
  c5->SaveAs("Plots/FRvsEta.pdf");

  phoEta_Sig_pix->SetLineColor(kRed);
  phoET_Sig_pix->SetLineColor(kRed);
  phoPhi_Sig_pix->SetLineColor(kRed);

  TCanvas* c6 = new TCanvas("c6","",10,10,1400, 600);
  	c6->Divide(3,1);
	c6->cd(1);
		phoEta_Sig_pix->Draw();
		phoEta_Sig_conv->Draw("same");
	c6->cd(2);
		phoET_Sig_pix->Draw();
		phoET_Sig_conv->Draw("same");
	c6->cd(3);
		phoPhi_Sig_pix->Draw();
		phoPhi_Sig_conv->Draw("same");
  c6->SaveAs("Plots/PVvsCV.pdf");

  TCanvas* c7 = new TCanvas("c7","",10,10,1400, 600);
	c7->Divide(3,1);  
	c7->cd(1);
  		phoEta_SigInv_conv->Draw("e");
  	c7->cd(2);
		phoET_SigInv_conv->Draw("e");
	c7->cd(3);
		phoPhi_SigInv_conv->Draw("e");
  c7->SaveAs("Plots/phoInvVeto.pdf");


/*  TCanvas* c8 = new TCanvas("c8","",10,10,1400, 600);
  TH1D* efficiency_conv_Num = (TH1D*)phoET_Sig_conv->Clone();
  TH1D* efficiency_conv_Denom = (TH1D*)phoET_SigInv_conv->Clone();
  TH1D* efficiency_conv = (TH1D*)efficiency_conv_Num->Divide(efficiency_conv_Denom);


  TH1D* efficiency_pix_Num = (TH1D*)phoET_Sig_pix->Clone();
  TH1D* efficiency_pix_Denom = (TH1D*)phoET_SigInv_pix->Clone();
  TH1D* efficiency_pix = (TH1D*)efficiency_pix_Num->Divide(efficiency_pix_Denom);

        c8->Divide(3,1);
                c8->cd(1);
                        phoET_Sig_conv->Draw();
                        phoET_Sig_pix->Draw("same");
                c8->cd(2);
                        phoET_SigInv_conv->Draw();
                        phoET_SigInv_pix->Draw("same");
                c8->cd(3);
                        efficiency_conv->Draw();
                        efficiency_pix->Draw("same");
  c8->SaveAs("Plots/Efficiency.pdf");

*/
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
  { { 0.553, 0.058, 0.019 },
    { 0.062, 0.020, 0.016 } };

const float sieieCut[2][nWP] =
  { {0.0099, 0.0099, 0.0099},
    {0.0284, 0.0268, 0.0263} };

const float chIsoCut[2][nWP] =
  { {2.49, 1.91, 1.61},
    {1.04, 0.82, 0.69} };

const float nhIso_A[2][nWP] =
  { {15.43, 4.66, 3.98},
    {19.71, 14.65, 4.52} };

const float nhIso_B[2][nWP] =
  { {0.007, 0.007, 0.007},
    {0.0129, 0.0129, 0.0129} };

const float phIso_A[2][nWP] =
  { {9.42, 4.29, 3.01},
    {11.88, 4.06, 3.61} };

const float phIso_B[2][nWP] =
  { {0.0033, 0.0033, 0.0033},
    {0.0108, 0.0108, 0.0108} };

bool passWorkingPoint(WpType iwp, bool isBarrel, float pt,
		      float hOverE, float full5x5_sigmaIetaIeta, 
		      float chIso, float nhIso, float phIso){

  int ieta = 0;
  if( !isBarrel )
    ieta = 1;

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
  if(verbose) printf("Check match for photon eta= %f   phi= %f\n", pEta, pPhi);
  for(unsigned int imc = 0; imc < (*mcPID).size(); imc++){//uint? imc...

    if(verbose) printf("   Check next particle: pid= %d  status= %d   mom= %d  (eta,phi)=(%f,%f)\n",
		       (*mcPID)[imc], (*mcStatus)[imc], (*mcMomPID)[imc], (*mcEta)[imc], (*mcPhi)[imc]);
    double msts = (*mcStatus)[imc];
    if((msts != 1)||((*mcPID)[imc] != 22))continue;
    if(verbose)printf("      passed pid and status\n");

    if((fabs((*mcMomPID)[imc]) !=21) 
       &&(fabs((*mcMomPID)[imc]) !=1)
       &&(fabs((*mcMomPID)[imc]) !=2)
       &&(fabs((*mcMomPID)[imc]) !=3)
       &&(fabs((*mcMomPID)[imc]) !=4)
       &&(fabs((*mcMomPID)[imc]) !=5)
       &&(fabs((*mcMomPID)[imc]) !=6))continue;    
    if(verbose) printf("       passed mother\n");
    
    double meta = (*mcEta)[imc];
    double mphi = (*mcPhi)[imc];
    
    TVector3 mcphoton;
    TVector3 recoPHOTOn;
    mcphoton.SetPtEtaPhi(1.0,meta,mphi);
    recoPHOTOn.SetPtEtaPhi(1.0,pEta,pPhi);
    
    double DR = mcphoton.DrEtaPhi(recoPHOTOn);
    if(verbose) printf("        dR=%f\n", DR);
    if(DR < 0.1 ){
      isMatched = true;
      if(verbose)printf("         PASSE ALL\n");
      // genPt = (*mcPt)[imc];
      break;
    }  
  }

  return isMatched;
}

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
