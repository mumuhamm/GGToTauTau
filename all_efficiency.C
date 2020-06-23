// Author Muhammad Alibordi
// GG To TauTau analysis

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TH1.h>
#include <string>
#include <TMath.h>
#include <vector>

#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TBuffer.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
using namespace RooFit;
using namespace std;

void all_efficiency(){
    
    TFile *fIn1 = new TFile("trktrig_eff_UPCSMuon_MC.root.root");
    TFile *fIn2 = new TFile("trktrig_eff_UPCSMuon_Open.root");
    //HLT_HIUPC_SingleMu0_NotMBHF2AND_v*
    TH1F* h_pt_denom = (TH1F*) fIn1->Get("eff/h_pt_probe_all");
    TH1F* h_pt_num = (TH1F*) fIn1->Get("eff/h_pt_probe_pass");
    //HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v*
    TH1F* h_pt_denom_open = (TH1F*) fIn2->Get("eff/h_pt_probe_all");
    TH1F* h_pt_num_open = (TH1F*) fIn2->Get("eff/h_pt_probe_pass");
    
    TCanvas* c_pt = new TCanvas("c_pt","c_pt");
    c_pt->cd();
    TH2F* h_pt_axis = new TH2F("h_pt_axis",";p_{T} (GeV);#epsilon_{trigger}(p_{T})",100,0,20,20,0,1.2);
    h_pt_axis->GetYaxis()->SetTitleOffset(0.98);
    h_pt_axis->Draw();
    
    TEfficiency* h_pt_eff = new TEfficiency(*h_pt_num, *h_pt_denom);
    h_pt_eff->SetLineColor(kRed);
    h_pt_eff->SetMarkerColor(kRed);
    h_pt_eff->Draw("pe same");
    
    TEfficiency* h_pt_eff_open = new TEfficiency(*h_pt_num_open, *h_pt_denom_open);
    h_pt_eff_open->SetLineColor(kBlue);
    h_pt_eff_open->SetMarkerColor(kBlue);
    h_pt_eff_open->Draw("pe same");
    
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(h_pt_eff,"HLT_HIUPC_SingleMu0_NotMBHF2AND","lep");
    leg->AddEntry(h_pt_eff_open,"HLT_HIUPC_SingleMuOpen_NotMBHF2AND","lep");
    leg->Draw();
    c_pt->Update();

}
