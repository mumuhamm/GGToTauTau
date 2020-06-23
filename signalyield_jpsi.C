/// \author 07/2018 - Muhammad Alibordi 
// Test of RooKeyPDF has ability to discriminate the background and

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit ;

TH1* makeTH1() ;
TTree* makeTTree() ;
void signalyield_jpsi()
{

  Int_t nbins = 100;
  //TChain* chain_data = new TChain("eff/eff_muon");
 // chain_data->Add("/Users/md/Documents/GGToTauTau/muonstudy/trktrig_eff_UPCSMuon_Open_Jpsi.root");
    TFile *fIn1 = new TFile("abordi_trktrig_eff.root");
    TH1F* Mass = (TH1F*) fIn1->Get("eff/h_trigmll_pass");
  //Int_t nevt = (int)chain_data->GetEntries();
  //std::cout<<"Number of total events"<<nevt<<"\n";
    
  RooRealVar *trkdimuon = new RooRealVar("trkdimuon", "M_{J/#psi} GeV/c^{2}",2.9,3.34);
  RooDataHist *data = new RooDataHist("data", "data", *trkdimuon, Import(*Mass));
    
  
  
    
    RooRealVar *mu= new RooRealVar("mu", "mu", 3.016, 3.0, 3.2);
    RooRealVar *lambda = new RooRealVar("lambda", "lambda", 1.044, 0.0, 1);
    RooRealVar *gamma = new RooRealVar("gamma", "gamma", 2.0, 0.0, 0.1);
    RooRealVar *delta = new RooRealVar("delta", "delta", 1, 0, 10);
    
    RooRealVar *mean = new RooRealVar("mean","mean",3.016, 2.9, 3.2) ;
    RooRealVar *sigma = new RooRealVar("sigma1","sigma1",0.0,0.,3) ;
    
    RooJohnson *john = new RooJohnson("john", "john", *trkdimuon, *mu, *lambda, *gamma, *delta);
    
    RooGaussian *gauss = new RooGaussian("gauss1","gauss1",*trkdimuon,*mean,*sigma) ;
    RooRealVar *nSig = new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC",500,0.,10e5);
    RooExtendPdf *sigpdf_inter = new RooExtendPdf("sigpdf_inter", "extended signal p.d.f",*john,*nSig);
    RooFitResult* fitRes = sigpdf_inter->fitTo(*data,Save(), Extended(1));//data_SigReg
    fitRes->Print("v");
  
    RooPlot* jpsimass = trkdimuon->frame(Title("M_{J/#psi} (GeV/c^{2}) "),Bins(100));
    data->plotOn(jpsimass,DataError(RooAbsData::SumW2));
    sigpdf_inter->plotOn(jpsimass) ;
    sigpdf_inter->paramOn(jpsimass);
    RooPlot* pullframe = trkdimuon->frame(RooFit::Title("Mass pull"));
    RooHist* hpull1 = jpsimass->pullHist();
    
    pullframe->addPlotable(hpull1,"P0") ;
    pullframe->SetMinimum(-3) ;
    pullframe->SetMaximum(+3) ;
    pullframe->SetYTitle("pull");
    pullframe->SetMarkerStyle(20);
    pullframe->SetNdivisions(10);
    TCanvas * pull_Can = new TCanvas("pull_Can", "mass Pull", 800, 200);
    pullframe->Draw();
    Double_t chisquare_mass = jpsimass->chiSquare();
    cout<<"Chi square of mass fit is :"<< chisquare_mass<< endl;
    sigpdf_inter->plotOn(jpsimass, Components(*gauss), LineColor(3), LineWidth(1), LineStyle(4));

    TCanvas *c = new TCanvas("c", "c",0,0,600,600);
    
    jpsimass->Draw();
    

}
