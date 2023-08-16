#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooAbsData.h"
#include "RooGaussian.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooUnblindUniform.h"
#include "RooExtendPdf.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TKey.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooDoubleCBFast.h"

#include "../interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "../../tdrStyle/tdrstyle.C"
#include "../../tdrStyle/CMS_lumi.C"

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

bool BLIND = false;
bool VETO = false;
bool runFtestCheckWithToys=false;

// for blind=true
// these values are configured in main
float mN = 2.75;
float sigma = 0.025;
int nsigma = 10;
float mN_low  = 2.625;
float mN_high = 2.875;
float mN_low_unblind = 2.7;
float mN_high_unblind = 2.8;
float jpsi_veto_low = 3.097-0.07;
float jpsi_veto_high = 3.097+0.07;
float d0_veto_low = 1.76-0.05;
float d0_veto_high = 1.76+0.05;

// not configured in main
int nBinsForFit  = 60; // kept baseline values for Hgg 
int nBinsForPlot = 60;  // ""

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
/* Construct the pdf model using the PdfModelBuilder */
  
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Polynomial") return pdfsModel.getPolynomial(Form("%s_poly%d",ext,order),order); 
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries){
/* Basic fitting routine, fit is not extended */

  int ntries=0;
  RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  
//  mass->setRange("unblindReg_1",mN_low,mN_low);
  data->Print("v");
  params_test->Print("v");
  int stat=1;
  double minnll=10e8;
//  data->var("hnl_mass")->setRange("left",mN_low,jpsi_veto_low);	
//  data->var("hnl_mass")->setRange("right",jpsi_veto_high,mN_high);	
  while (stat!=0){
    if (ntries>=MaxTries) break;
    RooFitResult *fitTest;
    if (VETO){
	std::cout << "vetoing resonance!___________________" << std::endl;
//	 RooFitResult *fitTest_left = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Range("left")); 
//	 RooFitResult *fitTest_right = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Range("right")); 
	 fitTest = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Range("reso,blind_left,blind_right"), RooFit::SumW2Error(kTRUE));
	}
    else  fitTest = pdf->fitTo(*data,RooFit::Save(1),RooFit::Minimizer("Minuit2","minimize"),RooFit::Range("blind_left,blind_right"));
    stat = fitTest->status();
    minnll = fitTest->minNll();
    if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
    ntries++; 
  }
  *stat_t = stat;
  *NLL = minnll;
}

double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooDataSet *data, std::string name){
/* Get probability of the f-test. Currently toys are not used and only simple TMath::Prob is used. */
 
  double prob_asym = TMath::Prob(chi2,ndof);
  if (!runFtestCheckWithToys) return prob_asym;

  int ndata = data->sumEntries();
   
  RooFitResult *fitNullData;// fit the pdfs to the data and keep this fit Result (for randomizing)
  RooFitResult *fitTestData;
  
  if (VETO){
  fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Range("reso,blind_left,blind_right")); //FIXME
  fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Range("reso,blind_left,blind_right")); 
  }else{

  fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Range("blind_left,blind_right")); //FIXME
  fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
    ,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1),RooFit::Range("blind_left,blind_right")); 


  }
  // Ok we want to check the distribution in toys then 
  // Step 1, cache the parameters of each pdf so as not to upset anything 
  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);
 
  int ntoys = 5000;
  TCanvas *can = new TCanvas();
  can->SetLogy();
  TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
  TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
  TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

  TGraph *gChi2 = new TGraph();
  gChi2->SetLineColor(kGreen+2);
  double w = toyhist.GetBinWidth(1);

  int ipoint=0;

  for (int b=0;b<0/*toyhist.GetNbinsX()*/;b++){
    double x = toyhist.GetBinCenter(b+1);
    if (x>0){
      gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
      ipoint++;
    }
  }
  int npass =0; int nsuccesst =0;
  mass->setBins(nBinsForFit);
  for (int itoy = 0 ; itoy < ntoys ; itoy++){

    params_null->assignValueOnly(preParams_null);
    params_test->assignValueOnly(preParams_test);
    RooDataSet *binnedtoy = pdfNull->generate(RooArgSet(*mass),ndata,0,1);

    int stat_n=1;
    int stat_t=1;
    int ntries = 0;
    double nllNull,nllTest;
    // Iterate on the fit 
    int MaxTries = 2;
    while (stat_n!=0){
      if (ntries>=MaxTries) break;
      RooFitResult *fitNull;
      if (VETO) fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Range("reso,blind_left,blind_right"));
      else fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Range("blind_left,blind_right"));
      //,RooFit::Optimize(0));

      nllNull = fitNull->minNll();
      stat_n = fitNull->status();
      if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
      ntries++; 
    }
  
    ntries = 0;
    while (stat_t!=0){
      if (ntries>=MaxTries) break;
      RooFitResult *fitTest;
      if (VETO) fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Range("reso,blind_left,blind_right"));
      else fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Range("blind_left,blind_right"));
      nllTest = fitTest->minNll();
      stat_t = fitTest->status();
      if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars()); 
      ntries++; 
    }
       
    toyhistStatN.Fill(stat_n);
    toyhistStatT.Fill(stat_t);

    if (stat_t !=0 || stat_n !=0) continue;
    nsuccesst++;
    double chi2_t = 2*(nllNull-nllTest);
    if (chi2_t >= chi2) npass++;
        toyhist.Fill(chi2_t);

  } // end loop over toys

  double prob=0;
  if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
  toyhist.Scale(1./(w*toyhist.Integral()));
  toyhist.Draw();
  TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
  lData.SetLineWidth(2);
  lData.Draw();
  gChi2->Draw("L");
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
  std::cout << "debug get probability f test: " << name.c_str() << std::endl; 
  can->SaveAs(name.c_str());

  TCanvas *stas =new TCanvas();
  toyhistStatN.SetLineColor(2);
  toyhistStatT.SetLineColor(1); 
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.87); leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(&toyhistStatN,"Null Hyp","L");
  leg->AddEntry(&toyhistStatT,"Test Hyp","L");
  toyhistStatN.Draw();
  toyhistStatT.Draw("same");
  leg->Draw();
  stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));
  //reassign params
  params_null->assignValueOnly(preParams_null);
  params_test->assignValueOnly(preParams_test);

  delete can; delete stas;
  delete gChi2;
  delete leg;
  delete lat;

  // Still return the asymptotic prob (usually its close to the toys one)
  return prob_asym;

}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name){
/* Get goodness of fit, based on chi-square, using binned dataset and fitted pdf 
   use toys or chi-square distributions depending on avg number of events in bin */

  double prob;
  int ntoys = 500;
  // Routine to calculate the goodness of fit. 
  name+="_gofTest.pdf";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  //norm.removeRange();

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
 // pdf->fitTo(*data,RooFit::Range("left,right"));
  if (VETO){
  data->plotOn(plot_chi2,Binning(nBinsForFit),Name("data"),RooFit::CutRange("reso,blind_left,blind_right"));

  pdf->plotOn(plot_chi2,Name("pdf"),/*RooFit::Range("reso,blind_left,blind_right"),*/RooFit::NormRange("reso,blind_left,blind_right"));
  }else{

  data->plotOn(plot_chi2,Binning(nBinsForFit),Name("data"),RooFit::CutRange("blind_left,blind_right"));

  pdf->plotOn(plot_chi2,Name("pdf"),/*RooFit::Range("full,blind_left,blind_right"),*/RooFit::NormRange("blind_left,blind_right"));

  }
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
 
  if ((double)data->sumEntries()/nBinsForFit < 5 ){

    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    // store pre-fit params 
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();
 
    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
      std::cout << "[INFO] " <<Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << std::flush;
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataSet *binnedtoy = pdf->generate(RooArgSet(*mass),nToyEvents,0,1);
      //RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);
//      pdf->fitTo(*binnedtoy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(0)); 

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t);
      if (VETO){ 
             pdf->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1)/*,RooFit::Range("reso,blind_left,blind_right")*/);
      pdf->plotOn(plot_t/*RooFit::Range("reso,blind_left,blind_right"),RooFit::NormRange("reso,blind_left,blind_right")*/);//,RooFit::NormRange("fitdata_1,fitdata_2"));
      }else{
	  pdf->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1)
                                              ,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1)/*,RooFit::Range("blind_left,blind_right")*/);
      pdf->plotOn(plot_t/*RooFit::Range("blind_left,blind_right"),RooFit::NormRange("blind_left,blind_right")*/);//,RooFit::NormRange("fitdata_1,fitdata_2"));
	}
      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForFit-np));
      delete plot_t;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;

    TCanvas *can = new TCanvas();
    double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    TArrow lData(chi2*(nBinsForFit-np),toyhist.GetMaximum(),chi2*(nBinsForFit-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    // back to best fit   
    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2*(nBinsForFit-np),nBinsForFit-np);
  }
  std::cout << "[INFO] GOF Chi2 in Observed =  " << chi2*(nBinsForFit-np) << std::endl;
  std::cout << "[INFO] GOF p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name,vector<string> flashggCats_, int status, double *prob){
/* Plot single pdf vs data, with pulls */
    
  // Chi2 taken from full range fit

  RooPlot *plot_chi2 = mass->frame();
  if (VETO){
	pdf->plotOn(plot_chi2/*,RooFit::Range("reso,blind_left,blind_right")*/,RooFit::NormRange("reso,blind_left,blind_right"));
        data->plotOn(plot_chi2,Binning(nBinsForFit),RooFit::Range("reso,blind_left,blind_right"));
  }
  else{
	 pdf->plotOn(plot_chi2/*,RooFit::Range("full,blind_left,blind_right")*/,RooFit::NormRange("blind_left,blind_right"));
         data->plotOn(plot_chi2,Binning(nBinsForFit),RooFit::Range("blind_left,blind_right"));
  }
  double integral = pdf->getNorm(*mass);
  std::cout << "___________________________________________"<<  integral << std::endl;
  int np = pdf->getParameters(*data)->getSize()+1; //Because this pdf has no extend
  double chi2 = plot_chi2->chiSquare(np);
 
  *prob = getGoodnessOfFit(mass,pdf,data,name);
  RooPlot *plot = mass->frame();
  if (VETO) {
   data->plotOn(plot,Binning(nBinsForPlot),CutRange("reso,blind_left,blind_right"));
//   data->plotOn(plot,/*Binning(nBinsForPlot),*/RooFit::CutRange("blind_right"));
 //  data->plotOn(plot,/*Binning(nBinsForPlot),*/RooFit::CutRange("reso"));
/*   data->plotOn(plot,Binning(nBinsForPlot),RooFit::CutRange("left"));*/
  //  data->plotOn(plot,Binning(nBinsForPlot),Invisible());
/*   data->plotOn(plot,Binning(nBinsForPlot),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(nBinsForPlot),Invisible());*/
  }
  else data->plotOn(plot,RooFit::CutRange("blind_left,blind_right"));

  TCanvas *canv = new TCanvas();
  canv->Divide(1,2);
  canv->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.10);
  gPad->SetPad(0.01,0.2,0.99,0.99);
  plot->GetXaxis()->SetTitleSize(0.04);
  plot->GetXaxis()->SetTitle("m_{\ell#pi} (GeV)");
  plot->GetYaxis()->SetTitleSize(0.04);;
  //plot->GetYaxis()->SetTitle("Entries");
  plot->GetYaxis()->SetTitleOffset(1.35);

  if (VETO)  pdf->plotOn(plot,RooFit::NormRange("reso,blind_left,blind_right")/*,RooFit::Range("reso,blind_left,blind_right")*/);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  else  pdf->plotOn(plot,RooFit::NormRange("blind_left,blind_right")/*,RooFit::Range("full")*/);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  pdf->paramOn(plot,RooFit::Layout(0.34,0.85,0.89),RooFit::Format("NEA",AutoPrecision(1)));
  plot->getAttText()->SetTextSize(0.025);
  plot->SetMaximum(plot->GetMaximum()*1.4);
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->Draw();
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.034);
  lat->DrawLatex(0.15,0.94,Form("#chi^{2}/ndof = %.3f, Prob = %.2f, Fit Status = %d ",chi2,*prob,status));

  TLatex *lab = new TLatex();
  lab->SetNDC();
  lab->SetTextFont(42);
  lab->SetTextSize(0.034);
  std::string delimiter = "/";
  size_t pos = 0;
  std::string s = name;
  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    //std::cout << token << std::endl;
    s.erase(0, pos + delimiter.length());
  }
  //delimiter = ".";
  //token = s.substr(0, s.find(delimiter));
  //std::cout << s << std::endl;
  lab->DrawLatex(0.55,0.25,s.c_str());

  // pulls
  canv->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetPad(0.01,0.01,0.99,0.2);
  gPad->SetGridy();
  RooHist* hpull = plot->pullHist();
  RooPlot *plot2 = mass->frame();
  plot2->GetYaxis()->SetNdivisions(504);
  plot2->GetYaxis()->SetLabelSize(0.17);
  plot2->GetYaxis()->SetTitleSize(0.17);
  plot2->GetYaxis()->SetTitleOffset(0.24);
  //plot2->GetYaxis()->SetRangeUser(-3.0,3.0);
  plot2->SetMinimum(-3.);
  plot2->SetMaximum(3.);
  plot2->GetYaxis()->SetTitle("Pulls") ;
  plot2->GetXaxis()->SetTitle("");
  plot2->GetXaxis()->SetLabelOffset(5);
  plot2->addPlotable(hpull,"P"); 
  plot2->Draw();

  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  //plot_chi2->Draw();
  //canv->SaveAs((name+"debug").c_str());

  delete canv;
  delete lat;
}
void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
/* Plot MultiPdf vs data */
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TLegend *leg = new TLegend(0.6,0.65,0.95,0.90);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",mN_low,mN_low_unblind);
  mass->setRange("unblindReg_2",mN_high_unblind,mN_high);
  if (VETO) {
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("reso,blind_left,blind_right"));
//    data->plotOn(plot,Binning(nBinsForPlot),CutRange("blind_left"));
//    data->plotOn(plot,Binning(nBinsForPlot),CutRange("blind_right"));
    data->plotOn(plot,Binning(nBinsForPlot),Invisible());
  }
  else{

  //  data->plotOn(plot,Binning(nBinsForPlot),CutRange("full"));
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("blind_left,blind_right"));
  //  data->plotOn(plot,Binning(nBinsForPlot),CutRange("blind_right"));
    data->plotOn(plot,Binning(nBinsForPlot),Invisible());
  }  
  TCanvas *canv = new TCanvas();
  //TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
  //pad1->SetBottomMargin(0.18);
  //pad1->Draw();
  //pad1->cd();

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - %s",flashggCats_[cat].c_str()),"LEP");
  int style=1;
  RooAbsPdf *pdf;
  RooCurve *nomBkgCurve;
  int bestcol= -1;
  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    if (VETO){
  //      pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::Range("left"));  
    //    pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::Range("right"));  
        pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::Range("reso,blind_left,blind_right"));  
	pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style),RooFit::Range("full"),RooFit::NormRange("reso,blind_left,blind_right"));
    }else{
        pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::Range("blind_left,blind_right"));  
	pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style),RooFit::Range("full"),RooFit::NormRange("blind_left,blind_right"));
	}
  //  pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style),RooFit::Range("full"),RooFit::NormRange("full"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) {
      ext=" (Best Fit Pdf) ";
      pdf= pdfs->getCurrentPdf();
      nomBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
      bestcol = col;
    }
    string pdfName = pdfs->getCurrentPdf()->GetName(); 
    std::string token;
    std::string delimiter = "_";
    size_t pos = 0;
    while ((pos = pdfName.find(delimiter)) != std::string::npos) {
      token = pdfName.substr(0, pos);
      pdfName.erase(0, pos + delimiter.length());
    }
    leg->AddEntry(pdfLeg,Form("%s%s",pdfName.c_str(),ext.c_str()),"L");
  }
  plot->SetTitle(Form("Category %s",flashggCats_[cat].c_str()));
  plot->SetMaximum(plot->GetMaximum()*1.4);
  plot->GetXaxis()->SetTitle("m_{l\pi} (GeV)");
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
/* Plot several Pdfs vs data, without ratio plot, (used for the "truth") */
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",mN_low,mN_low_unblind);
  mass->setRange("unblindReg_2",mN_high_unblind,mN_high);
  if (BLIND) {
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(nBinsForPlot),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(nBinsForPlot),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForPlot));

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  if(flashggCats_.size() >0){
  leg->AddEntry(datLeg,Form("Data - %s",flashggCats_[cat].c_str()),"LEP");
  } else {
  leg->AddEntry(datLeg,Form("Data - %d",cat),"LEP");
  }
  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    if (VETO)it->second->plotOn(plot,LineColor(col),LineStyle(style),RooFit::Range("full"),RooFit::NormRange("reso,blind_left,blind_right"));
    else it->second->plotOn(plot,LineColor(col),LineStyle(style),RooFit::Range("full"),RooFit::NormRange("blind_left,blind_right"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
    i++;
  }
  plot->SetMaximum(plot->GetMaximum()*1.4);
  plot->SetTitle(Form(" %s",flashggCats_[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( canv, 0, 0);
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
}

void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}

int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooCategory *cat, bool silent=false){
/* Get index of the best fit pdf (minimum NLL, including correction) among functions in the multipdf.
   All fits are performed again. */

  double global_minNll = 1E10;
  int best_index = 0;
  int number_of_indeces = cat->numTypes();
    
  RooArgSet snap,clean;
  RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
  params->remove(*cat);  // pdf_index is RooCategory, removed from parameters
  params->snapshot(snap);
  params->snapshot(clean);
  if (!silent) {
    //params->Print("V");
  }
 
  // Uncomment to try to make converge a failed fit
  //bkg->setDirtyInhibit(1);
  //RooAbsReal *nllm = bkg->createNLL(*data);
  //RooMinimizer minim(*nllm);
  //minim.setStrategy(1);
  
  for (int id=0;id<number_of_indeces;id++){    
    params->assignValueOnly(clean);
    cat->setIndex(id);

    //RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

    if (!silent) {
      /*
      std::cout << "BEFORE  MAKING FIT" << std::endl;
      params->Print("V");
      std::cout << "-----------------------" << std::endl;    
      */
    }
    
    //minim.minimize("Minuit2","minimize");
    double minNll=0; //(nllm->getVal())+bkg->getCorrection();
    int fitStatus=1;    
    runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/7);
    // Add the penalty

    minNll=minNll+bkg->getCorrection();

    if (!silent) {
      /*
      std::cout << "After Minimization ------------------  " <<std::endl;
      std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
      bkg->Print("v");
      bkg->getCurrentPdf()->getParameters(*data)->Print("V");
      std::cout << " ------------------------------------  " << std::endl;
  
      params->Print("V");
      */
      std::cout << "[INFO] AFTER FITTING" << std::endl;
      std::cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
      std::cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<std::endl;
      std::cout << "[INFO] NLL + c = " <<  minNll << std::endl;
      std::cout << "-----------------------" << std::endl;
    }
      
    if (minNll < global_minNll){
      global_minNll = minNll;
      snap.assignValueOnly(*params);
      best_index=id;
    }
  } // end loop over pdf_index 
  cat->setIndex(best_index);
  params->assignValueOnly(snap);

  std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
  std::cout << "[INFO] Best fit parameters " << std::endl;
  params->Print("V");
  
  return best_index;
}

int main(int argc, char* argv[]){
 
  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  string year_ = "2016";
  //int year_ = 2017;

  string fileName;
  int ncats;
  int singleCategory;
  int catOffset;
  string datfile;
  string outDir;
  string outfilename;
  bool is2011=false;
  bool verbose=false;
  bool saveMultiPdf=false;
  int isFlashgg_ =1;
  string flashggCatsStr_;
  vector<string> flashggCats_;
  bool isData_ =0;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("singleCat", po::value<int>(&singleCategory)->default_value(-1),                           "Run A single Category")
    ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&outfilename),                                           "Save a MultiPdf model with the appropriate pdfs")
    ("runFtestCheckWithToys",                                                                   "When running the F-test, use toys to calculate pvals (and make plots) ")
    ("is2011",                                                                                  "Run 2011 config")
    ("is2012",                                                                                  "Run 2012 config")
    ("blind",                                                                                   "blind plots")
    ("isFlashgg",  po::value<int>(&isFlashgg_)->default_value(1),                               "Use Flashgg output ")
    ("isData",  po::value<bool>(&isData_)->default_value(0),                                    "Use Data not MC ")
    ("flashggCats,f", po::value<string>(&flashggCatsStr_)->default_value("UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag"),                  "Flashgg category names to consider")
    ("year", po::value<string>(&year_)->default_value("2016"),                                  "Dataset year")
    ("catOffset", po::value<int>(&catOffset)->default_value(0),                                 "Category numbering scheme offset")
    ("mN", po::value<float>(&mN)->default_value(2.75),                                          "Mass of the peak, for center of window")
    ("sigma", po::value<float>(&sigma)->default_value(0.025),                                   "Sigma of the peak, for size of window")
    ("nsigma", po::value<int>(&nsigma)->default_value(10),                                       "Sigma multiplier, for size of window")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
//  if (vm.count("blind")) BLIND=true;
  saveMultiPdf = vm.count("saveMultiPdf");

  if (vm.count("verbose")) verbose=true;
  if (vm.count("runFtestCheckWithToys")) runFtestCheckWithToys=true;

  std::cout << "DEBUG mN=" << mN << std::endl; 

  double veto_low[3] = {1.814,3.047,3.636};
  double veto_high[3] = {1.914,3.147,3.736};
  double vetolow;
  double vetohigh;
  mN_low  = mN - nsigma * sigma;
  mN_high = mN + nsigma * sigma;
  mN_low_unblind =  mN - 2 * sigma;
  mN_high_unblind = mN + 2 * sigma;
/*  if (mN ==3){
	
	jpsi_veto_low = 3.097-0.05
	jpsi_veto_high = 3.097+0.05

  }
  if (mN ==2){
	
	d0_veto_low = 1.76-0.05
	d0_veto_high = 1.76+0.05

  }*/

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }
  split(flashggCats_,flashggCatsStr_,boost::is_any_of(","));
  
  int startingCategory=0;
  if (singleCategory >-1){
    ncats=singleCategory+1;  
    startingCategory=singleCategory;
  }
  if (isFlashgg_==1){
    ncats= flashggCats_.size();
  }

  if(verbose) std::cout << "[INFO] SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
  outputfile = new TFile(outfilename.c_str(),"RECREATE");
  outputws = new RooWorkspace(); outputws->SetName("multipdf");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile;
  inFile = TFile::Open(fileName.c_str());
  std::cout<< "check if file exists" << fileName << std::endl;
  RooWorkspace *inWS;
  if(isFlashgg_){
    if (isData_){
      inWS = (RooWorkspace*)inFile->Get("wS");
    } else {
      inWS = (RooWorkspace*)inFile->Get("wS");
    }
  } else {
    inWS = (RooWorkspace*)inFile->Get("wS");
  }
  if (verbose){
	 std::cout << "[INFO]  inWS open " << inWS << std::endl;
	 inWS->Print();
	}
  if (saveMultiPdf){
    transferMacros(inFile,outputfile);

    RooRealVar *intL; 
    RooRealVar *sqrts;

    if (isFlashgg_){
      //intL  = (RooRealVar*)inWS->var("IntLumi");
      intL  = intLumi_;
      sqrts = (RooRealVar*)inWS->var("SqrtS");
      if (!sqrts){ sqrts = new RooRealVar("SqrtS","SqrtS",13); }
    std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;


    } else {
      //intL  = (RooRealVar*)inWS->var("IntLumi");
      intL  = intLumi_;
      sqrts = (RooRealVar*)inWS->var("Sqrts");
    }
    outputws->import(*intL);
    outputws->import(*sqrts);
    std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;
  }

  // Set up which families of functions you want to test
  vector<string> functionClasses;
//  functionClasses.push_back("Bernstein");
  functionClasses.push_back("Exponential");
  functionClasses.push_back("PowerLaw");
  functionClasses.push_back("Laurent");
//  functionClasses.push_back("Chebychev");
//  functionClasses.push_back("Polynomial");
  map<string,string> namingMap;
//  namingMap.insert(pair<string,string>("Bernstein","pol"));
  namingMap.insert(pair<string,string>("Exponential","exp"));
  namingMap.insert(pair<string,string>("PowerLaw","pow"));
  namingMap.insert(pair<string,string>("Laurent","lau"));
//  namingMap.insert(pair<string,string>("Chebychev","che"));
 // namingMap.insert(pair<string,string>("Polynomial","pol"));

  FILE *resFile ;
  if  (singleCategory >-1) resFile = fopen(Form("%s/fTestResults_%s.txt",outDir.c_str(),flashggCats_[singleCategory].c_str()),"w");
  else resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
  vector<map<string,int> > choices_vec;
  vector<map<string,std::vector<int> > > choices_envelope_vec;
  vector<map<string,RooAbsPdf*> > pdfs_vec;

  PdfModelBuilder pdfsModel;
/*  inWS->var("hnl_mass")->setMin(mN_low);
  inWS->var("hnl_mass")->setMax(mN_high);*/
 /* inWS->var("hnl_mass")->setRange("blind_left",mN_low,mN_low_unblind);
  inWS->var("hnl_mass")->setRange("blind_right",mN_high_unblind,mN_high);*/
  std:: cout <<  mN_low << " " <<  jpsi_veto_low  << " " <<  mN_high << " " << jpsi_veto_high << std::endl;

  //loop to treat SM resonances
  for (int sm_idx= 0; sm_idx<3;sm_idx++){

   if ((mN_high <veto_low[sm_idx] || mN_low>veto_high[sm_idx])){ //define only blinded and unblinded ranges
   inWS->var("hnl_mass")->setRange("full",mN_low,mN_high);
   inWS->var("hnl_mass")->setRange("blind_left",mN_low,mN_low_unblind);
   inWS->var("hnl_mass")->setRange("blind_right",mN_high_unblind,mN_high);
  std:: cout <<"no restriction " <<  mN_low << " " <<  veto_low[sm_idx]  << " " <<  mN_high << " " << veto_high[sm_idx] << std::endl;
 
  }else if( mN_low<veto_low[sm_idx] && mN_high>veto_low[sm_idx] && mN_high<veto_high[sm_idx]){

  std:: cout <<"right edge " <<  mN_low << " " <<  veto_low[sm_idx]  << " " <<  mN_high << " " << veto_high[sm_idx] << std::endl;
   inWS->var("hnl_mass")->setRange("full",mN_low,veto_low[sm_idx]);	
   inWS->var("hnl_mass")->setRange("blind_left",mN_low,mN_low_unblind);
   inWS->var("hnl_mass")->setRange("blind_right",mN_high_unblind,veto_low[sm_idx]);
  break;
  }else if( mN_low<veto_low[sm_idx] && mN_high>veto_high[sm_idx]){
   
  std:: cout <<"in the middle " <<  mN_low << " " <<  veto_low[sm_idx]  << " " <<  mN_high << " " << veto_high[sm_idx] << std::endl;
   if (veto_high[sm_idx]>mN_high_unblind){
   inWS->var("hnl_mass")->setRange("reso",veto_high[sm_idx],mN_high);
   inWS->var("hnl_mass")->setRange("blind_left",mN_low,mN_low_unblind);
   inWS->var("hnl_mass")->setRange("blind_right",mN_high_unblind,veto_low[sm_idx]);
   }else if(veto_high[sm_idx]<mN_low_unblind){

   inWS->var("hnl_mass")->setRange("reso",mN_low,veto_low[sm_idx]);	
   inWS->var("hnl_mass")->setRange("blind_left",veto_high[sm_idx],mN_low_unblind);
   inWS->var("hnl_mass")->setRange("blind_right",mN_high_unblind,mN_high);

   }
   VETO = true;
   vetolow = veto_low[sm_idx];
   vetohigh =veto_high[sm_idx];
   break;
  }else if(mN_low>veto_low[sm_idx] && mN_low<veto_high[sm_idx] && mN_high>veto_high[sm_idx]){

  std:: cout <<"left edge" <<  mN_low << " " <<  veto_low[sm_idx]  << " " <<  mN_high << " " << veto_high[sm_idx] << std::endl;
   inWS->var("hnl_mass")->setRange("full",veto_high[sm_idx],mN_high);	
   inWS->var("hnl_mass")->setMin(veto_high[sm_idx]);	
   inWS->var("hnl_mass")->setRange("blind_left",veto_high[sm_idx],mN_low_unblind);
   inWS->var("hnl_mass")->setRange("blind_right",mN_high_unblind,mN_high);
  break;	
  }
  }
  inWS->var("hnl_mass")->Print();
/*   if ((mN_high <d0_veto_low || mN_low>d0_veto_high)){ //define only blinded and unblinded ranges
   inWS->var("hnl_mass")->setRange("full",mN_low,mN_high);
  
//  inWS->var("hnl_mass")->setRange("blinded_left",mN_low,mN_low_unblind);
//  inWS->var("hnl_mass")->setRange("blinded_right",mN_high_unblind,mN_high);
//  inWS->var("hnl_mass")->setRange("blinded",RooFit::Range("blinded_right","blinded_left"));

  }else if( mN_low<d0_veto_low && mN_high>d0_veto_low && mN_high<d0_veto_high) {

   inWS->var("hnl_mass")->setRange("full",mN_low,d0_veto_low);	

  }else if ( mN_low<d0_veto_low && mN_high>d0_veto_high){

   inWS->var("hnl_mass")->setRange("left",mN_low,d0_veto_low);	
   inWS->var("hnl_mass")->setRange("right",d0_veto_high,mN_high);
//   VETO = true;	


  }else if(mN_low>d0_veto_low && mN_low<d0_veto_high && mN_high>d0_veto_high){

   inWS->var("hnl_mass")->setRange("full",d0_veto_low,mN_high);	

  }
  if ((mN_high <jpsi_veto_low || mN_low>jpsi_veto_high)){ //define only blinded and unblinded ranges
  inWS->var("hnl_mass")->setRange("full",mN_low,mN_high);
	
 }else if( mN_low<jpsi_veto_low && mN_high>jpsi_veto_low && mN_high<jpsi_veto_high) {

   inWS->var("hnl_mass")->setRange("full",mN_low,jpsi_veto_low);	

  }else if ( mN_low<jpsi_veto_low && mN_high>jpsi_veto_high){

   inWS->var("hnl_mass")->setRange("left",mN_low,jpsi_veto_low);	
   inWS->var("hnl_mass")->setRange("right",jpsi_veto_high,mN_high);	
  // mass->setRange("left",mN_low,jpsi_veto_low);	
  // mass->setRange("right",jpsi_veto_high,mN_high);	
   VETO = true;	


  }else if(mN_low>jpsi_veto_low && mN_low<jpsi_veto_high && mN_high>jpsi_veto_high){

   inWS->var("hnl_mass")->setRange("full",jpsi_veto_low,mN_high);	

 }


*/
	
//  inWS->var("hnl_mass")->setRange("veto_d0",mN_low,mN_high);
  RooRealVar *mass = (RooRealVar*)inWS->var("hnl_mass"); 
  std:: cout << "[INFO] Got mass from ws " << mass << std::endl;
  pdfsModel.setObsVar(mass);
  double upperEnvThreshold = 0.1; // upper threshold on prob_ftest to include function in envelope (looser than truth function)
  double minGofThreshold = 0.01;  // minimal goodness of fit to include function in envelope

  fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
  fprintf(resFile,"\\hline\n");

  std::string ext = is2011 ? "7TeV" : "8TeV";
        if( isFlashgg_ ){
          if( year_ == "all" ){ ext = "13TeV"; }
          else{ ext = Form("%s_13TeV",year_.c_str()); }
        }

  std::cout << "[INFO] Number of categories to process: " << ncats << std::endl;
  for (int cat=startingCategory; cat<ncats; cat++){

    map<string,int> choices;
    map<string,std::vector<int> > choices_envelope;
    map<string,RooAbsPdf*> pdfs;
    map<string,RooAbsPdf*> allPdfs;
    string catname;
    if (isFlashgg_){
      catname = Form("%s",flashggCats_[cat].c_str());
    } else {
      catname = Form("cat%d",cat);
    }
    
    // Option 1: Use as input an unbinned RooDataSet and bin it
    /*
    RooDataSet *dataFull;
    RooDataSet *dataFull0;
    if (isData_) {
    dataFull = (RooDataSet*)inWS->data(Form("Data_13TeV_%s",catname.c_str()));
    if (verbose) std::cout << "[INFO] opened data for  "  << Form("Data_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
    }
    else 
    {dataFull = (RooDataSet*)inWS->data(Form("data_mass_%s",catname.c_str()));
    if (verbose) std::cout << "[INFO] opened data for  "  << Form("data_mass_%s",catname.c_str()) <<" - " << dataFull <<std::endl;
    }

    mass->setBins(nBinsForFit);
    RooDataSet *data;
    string thisdataBinned_name;

    if ( isFlashgg_){
      thisdataBinned_name =Form("CAT_roohist_data_mass_%s",flashggCats_[cat].c_str());
    } else {
      thisdataBinned_name= Form("CAT_roohist_data_mass_cat%d",cat);
    }
    RooDataHist thisdataBinned(thisdataBinned_name.c_str(),"data",*mass,*dataFull);
    data = (RooDataSet*)&thisdataBinned; 
    // both "data" and "thisdataBinned" are binned (number of bins is given by "mass" bins) 
    // "dataFull" is unbinned
    */

    // Option 2 (equivalent): Use as input a binned RooDataHist as is
    string data_name = Form("CAT_roohist_data_mass_%s",flashggCats_[cat].c_str()); 
    RooDataSet  *data_in       = (RooDataSet*)inWS->data("data_obs");
  //   RooDataSet* data_small;
   // if (BLIND) data_small =(RooDataSet*)  data_in->reduce(RooArgSet{*mass},("hnl_mass> "+std::to_string(mN_low)+" && hnl_mass<"+ std::to_string(mN_high)+" && hnl_mass>"+ std::to_string(mN_high_unblind)+" && hnl_mass<"+ std::to_string(mN_low_unblind)).c_str());
    std::cout << "entries dataset " << data_in->sumEntries() << std::endl;
    RooDataSet* data_small;

    data_small = (RooDataSet*)data_in->reduce(RooArgSet(*mass));

	
   // std::unique_ptr<RooAbsData> data_small; 
  // if (VETO)  data_small = (RooDataSet*)  data_small->reduce(("(hnl_mass> "+std::to_string(mN_low)+" && hnl_mass<"+ std::to_string(jpsi_veto_low) +") || (hnl_mass> "+std::to_string(jpsi_veto_high)+" && hnl_mass<"+ std::to_string(mN_high) +")").c_str());
   //else  data_small = (RooDataSet*)  data_small->reduce(RooFit::CutRange("full"));
   // if (VETO)  data_small = (RooDataSet*)  data_in->reduce(RooArgSet{*mass},("(hnl_mass> "+std::to_string(mN_low)+" && hnl_mass<"+ std::to_string(jpsi_veto_low) +") || (hnl_mass> "+std::to_string(jpsi_veto_high)+" && hnl_mass<"+ std::to_string(mN_high) +")" ).c_str());
  //  if (VETO)  data_small = (RooDataSet*)  data_in->reduce(RooArgSet{*mass},("hnl_mass> "+std::to_string(mN_low)+" && hnl_mass<"+ std::to_string(mN_high)).c_str());
   // std::cout << " (hnl_mass> "+std::to_string(jpsi_veto_high)+" && hnl_mass<"+ std::to_string(mN_high) <<std::endl;
    std::cout << "entries dataset " << data_small->sumEntries() << std::endl;
    TH1D* hist =(TH1D*) data_small->createHistogram("hnl_mass",nBinsForFit);
   /* if(VETO){
    for (int i=1; i<nBinsForFit;i++){
		if (hist->GetBinLowEdge(i)>jpsi_veto_low && (hist->GetBinLowEdge(i)+hist->GetBinWidth(i))<jpsi_veto_high){
			 hist->SetBinContent(i,0.);
			 hist->SetBinError(i,0.);

	}
    }
    }*/
    hist->SaveAs("testhist.root");
    std::cout << "entries histo " << hist->GetEntries() << std::endl;
    RooDataHist *thisdataBinned= new RooDataHist("data","data",*mass,hist);
 
  //  RooDataSet  *data       = (RooDataSet*)thisdataBinned;
    RooDataSet  *data       = (RooDataSet*)data_small;
    std::cout << "entries dataset " << data->sumEntries() << std::endl;
    
  


    
    RooArgList storedPdfs("store");

    fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
    fprintf(resFile,"\\hline\n");

    double MinimimNLLSoFar=1e10;
    int simplebestFitPdfIndex = 0;

    // Standard F-Test to find the truth functions
    for (vector<string>::iterator funcType=functionClasses.begin(); 
        funcType!=functionClasses.end(); funcType++){

      std::cout << "======================================= " << std::endl;
      std::cout << "====> FAMILY " << funcType->c_str() << std::endl;
      std::cout << "======================================= " << std::endl;

      double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.; 
      int order=1; int prev_order=0; int cache_order=0;

      RooAbsPdf *prev_pdf=NULL;
      RooAbsPdf *cache_pdf=NULL;
      std::vector<int> pdforders;

      std::cout << "===> F-TEST for Truth determination" << std::endl;

      int counter =0;
      while (prob<0.05 && order < 7){ 
        cout << "==> " << *funcType << " " << order << endl;
        RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("ftest_pdf_%d_%s",(cat+catOffset),ext.c_str()));
        if (!bkgPdf){
          // assume this order is not allowed
          order++;
        }
        else {

          //bkgPdf->Print();
          int fitStatus = 0;
          runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/7);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
          if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
       
          chi2 = 2.*(prevNll-thisNll);
          if (chi2<0. && order>1) chi2=0.;
          if (prev_pdf!=NULL){
            prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data
                ,Form("%s/Ftest_from_%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,(cat+catOffset)));
            std::cout << "[INFO] F-test Prob == " << prob << std::endl;
          } else {
            prob = 0;
          }
          double gofProb=0;
          // otherwise we get it later ...
          if (!saveMultiPdf) plot(mass,bkgPdf,data,Form("%s/%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb);
          cout << "[INFO] function type, order, prevNLL, thisNLL, chi2, prob " << endl;
          cout << "[INFO] " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
          prevNll=thisNll;
          cache_order=prev_order;
          cache_pdf=prev_pdf;
          prev_order=order;
          prev_pdf=bkgPdf;
          order++;
        }
        counter++;
      } // end condition for performing f-test

      // next line is commented, as we want to save only the final result (that takes into account both GOF and F-test results
      //fprintf(resFile,"%15s & %d & %5.3f & %5.3f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      choices.insert(pair<string,int>(*funcType,cache_order));
      pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));

      int truthOrder = cache_order;

      // Now run loop to determine functions inside envelope
      std::cout << "===> F-TEST and GOF for ENVELOPE determination" << std::endl;
      if (saveMultiPdf){
        chi2=0.;
        thisNll=0.;
        prevNll=0.;
        prob=0.;
        order=1;
        prev_order=0;
        cache_order=0;
        std::cout << "[INFO] Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;

        while (prob<upperEnvThreshold){
          cout << "==> " << *funcType << " " << order << endl;
          RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,Form("env_pdf_%d_%s",(cat+catOffset),ext.c_str()));
          if (!bkgPdf ){
            // assume this order is not allowed
            if (order >6) { std::cout << " [WARNING] could not add ] " << std::endl; break ;}
            order++;
          }
          else {
            // Fit and chi-square calculation is repeated
            //RooFitResult *fitRes;
            int fitStatus=0;
            runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/7);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
            //thisNll = fitRes->minNll();
            if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
            double myNll = 2.*thisNll;
            chi2 = 2.*(prevNll-thisNll);
            if (chi2<0. && order>1) chi2=0.;
            prob = TMath::Prob(chi2,order-prev_order); 

            cout << "[INFO] function type, order, prevNLL, thisNLL, chi2, prob " << endl;
            cout << "[INFO] " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
            prevNll=thisNll;
            cache_order=prev_order;
            cache_pdf=prev_pdf;

            // Calculate goodness of fit (will use toys for lowstats)
            double gofProb =0; 
            plot(mass,bkgPdf,data,Form("%s/%s%d_%s",outDir.c_str(),funcType->c_str(),order,catname.c_str()),flashggCats_,fitStatus,&gofProb);

            if ((prob < upperEnvThreshold) ) { // Looser requirements for the envelope

              if (gofProb > minGofThreshold || order == truthOrder ) {  // Good looking fit or one of our regular truth functions
              //if (gofProb > minGofThreshold) { // minimal requirement on the goodness of fit

                std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
                  << " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() <<  std::endl;
                allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
                storedPdfs.add(*bkgPdf);
                pdforders.push_back(order);

                // Keep track but we shall redo this later
                if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
                  simplebestFitPdfIndex = storedPdfs.getSize()-1;
                  MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
                }
            //  }
            }
            prev_order=order;
            prev_pdf=bkgPdf;
            order++;
          }
	}
        } // end while

        fprintf(resFile,"%15s & %d & %5.3f & %5.3f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
        choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
      }
    } // end loop over families

    fprintf(resFile,"\\hline\n");
    choices_vec.push_back(choices);
    choices_envelope_vec.push_back(choices_envelope);
    pdfs_vec.push_back(pdfs);

    plot(mass,pdfs,data,Form("%s/truths_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat);

    if (saveMultiPdf){
      // Put selectedModels into a MultiPdf
      string catindexname;
      string catname;
      if (isFlashgg_){
        catindexname = Form("pdfindex_%s_%s",flashggCats_[cat].c_str(),ext.c_str());
        catname = Form("%s",flashggCats_[cat].c_str());
      } else {
        catindexname = Form("pdfindex_%d_%s",(cat+catOffset),ext.c_str());
        catname = Form("cat%d",(cat+catOffset));
      }
      RooCategory catIndex(catindexname.c_str(),"c");
      RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg",data->sumEntries(),0,3*data->sumEntries());
      RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_%s_%s_bkgshape",catname.c_str(),ext.c_str()),"all pdfs",catIndex,storedPdfs);
       //nBackground.removeRange(); // bug in roofit will break combine until dev branch brought in
      //double check the best pdf!
      int bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,!verbose);
     // RooArgSet nset(*mass);
     
      RooRealVar nGaus(Form("gaus_veto_norm",catname.c_str(),ext.c_str()),"ngaus",data->sumEntries(),0,3*data->sumEntries());
      RooRealVar* mean = new RooRealVar(Form("sm_mean_",catname.c_str()),Form("sm_mean_",catname.c_str()),(vetolow+vetohigh)/2,(vetolow+vetohigh)/2);
      RooRealVar* sigma = new RooRealVar(Form("sm_sigma_",catname.c_str()),Form("sm_sigma_",catname.c_str()),(vetolow+vetohigh)/2*0.00697693+0.00185739,(vetolow+vetohigh)/2*0.00697693+0.00185739);
      RooRealVar* n1 = new RooRealVar(Form("sm_n1_",catname.c_str()),Form("sm_n1_",catname.c_str()),2.5,2.5);
      RooRealVar* alpha1 = new RooRealVar(Form("sm_alpha1_",catname.c_str()),Form("sm_alpha1_",catname.c_str()),1.2,1.2);
      RooRealVar* n2 = new RooRealVar(Form("sm_n2_",catname.c_str()),Form("sm_n2_",catname.c_str()),4,4);
      RooRealVar* alpha2 = new RooRealVar(Form("sm_alpha2_",catname.c_str()),Form("sm_alpha2_",catname.c_str()),1.5,1.5);
      //RooAbsPdf* gaus_veto = new RooGaussian("gaus_veto","gaus_veto",*mass,*mean,*sigma);
      RooDoubleCBFast gaus_veto("gaus_veto","gaus_veto",*mass,*mean,*sigma,*alpha1,*n1,*alpha2,*n2);	
    
      if (VETO){
      inWS->var("hnl_mass")->setRange("left",mN_low,mN_low_unblind);
      inWS->var("hnl_mass")->setRange("right",mN_high_unblind,mN_high);

      RooAbsPdf* bestFit =(RooAbsPdf*) storedPdfs.at(bestFitPdfIndex);
      RooAddPdf* bestFit_ext =new RooAddPdf("extendedBestFit","extendedBestFit",
					    RooArgList(*bestFit,gaus_veto),RooArgList(nBackground,nGaus));
      bestFit_ext->fitTo(*data,Range("left,right"));
	

      }

      catIndex.setIndex(bestFitPdfIndex);
      std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
      std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
      storedPdfs.Print();
      std::cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;

      std::cout << "[INFO] Simple check of index "<< simplebestFitPdfIndex <<std::endl;
      std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;

      mass->setBins(nBinsForFit);
      //RooDataHist dataBinned(Form("roohist_data_mass_%s",catname.c_str()),"data",*mass,*dataFull);

      // Save it (also a binned version of the dataset
      outputws->import(*pdf);
      if (VETO){
	 outputws->import(gaus_veto);
	 outputws->import(nGaus);
	}
      outputws->import(nBackground);
      outputws->import(catIndex);
      //outputws->import(dataBinned);
      outputws->import(*data);
      plot(mass,pdf,&catIndex,data,Form("%s/multipdf_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,bestFitPdfIndex);

    } // end if saveMultiPdf

  } // end loop over categories 

  if (saveMultiPdf){
    outputfile->cd();
    outputws->Write();
    outputfile->Close();  
  }

  // Write recommended options to screen and to file
/*  FILE *dfile = fopen(datfile.c_str(),"w");
  cout << "[RESULT] Recommended options based on truth" << endl;

  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",(cat+catOffset)); 
    for (map<string,int>::iterator it=choices_vec[cat-startingCategory].begin(); it!=choices_vec[cat-startingCategory].end(); it++){
      cout << "\t" << it->first << " - " << it->second << endl;
      fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
    }
    fprintf(dfile,"\n");
  }


  cout << "[RESULT] Recommended options for envelope" << endl;
  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",(cat+catOffset)); 
    for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat-startingCategory].begin(); it!=choices_envelope_vec[cat-startingCategory].end(); it++){
      std::vector<int> ords = it->second;
      for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
        cout << "\t" << it->first << " - " << *ordit << endl;
        fprintf(dfile,"envel=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
      }
    }
    fprintf(dfile,"\n");
  }*/

  inFile->Close();
}
