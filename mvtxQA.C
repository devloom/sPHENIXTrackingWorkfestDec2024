#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>

Double_t langaufun(Double_t *x, Double_t *par)
{
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) 
  {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}



void GetRunList(std::vector<int>& vect)
{
  std::string line;
  ifstream irunlist;
  irunlist.open("lists/runlist.txt");
  if(!irunlist)
  {
    std::cout << "failed to open runlist " << std::endl;
    exit(1);
  }
  while (std::getline(irunlist, line)) {
    vect.push_back(std::stoi(line));
  }
  irunlist.close(); 
}

// Note that many runs from 53332-53485 are not in lustre (just in hpss) so the QAhtml is missing from here
//, TH1F& h1mvtx_chip_bc, TH1F& h1mvtx_strobe_bc
int QA(int runnum, TFile &outputFile, std::vector<float>& deadstavechipeff, std::vector<float>& hitpeak, 
      float &reduced_chi2_chip, float &reduced_chi2_strobe)
{
  char dstpath[128];

  // Check if the QAhtml hit histograms exist for the ana441 build
  sprintf(dstpath,"/sphenix/data/data02/sphnxpro/QAhtml/aggregated/HIST_DST_TRKR_HIT_run2pp_ana441_2024p007-%08i-9000.root",runnum);
  TFile *qafile = new TFile(dstpath);
  if (!qafile || qafile->IsZombie()) 
  {
    std::cout << "QA file for ana441 does not exist for run: " << runnum << std::endl;
    std::cout << "Trying older build." << std::endl;
    
    // Try the QAhtml hit histograms for the older (new_p007) build
    sprintf(dstpath,"/sphenix/data/data02/sphnxpro/QAhtml/aggregated/HIST_DST_TRKR_HIT_run2pp_new_new_2024p007-%08i-9000.root",runnum);
    qafile = new TFile(dstpath);
    if (!qafile || qafile->IsZombie()) 
    {
      std::cout << "QA file for older build (new_p007) does not exist for run: " << runnum << std::endl;
      delete qafile;
      return 1;
    }
  }


  char hmvtx[128];
  TH2F *h_MvtxRawHitQA_nhits_stave_chip[3];
  int totalchips[3] = {9,9,9};
  int totalstaves[3] = {12,16,20};
  int totalstaveschips[3] = {9*12,9*16,9*20};

  TH1F *h_MvtxRawHitQA_nhits_layer[3];


  //chi2 test to see how uniform the chip_bc and strobe_bc are
  TH1F *h_MvtxRawHitQA_chip_bc = (TH1F*)qafile->Get("h_MvtxRawHitQA_chip_bc");
  TH1F *h_MvtxRawHitQA_strobe_bc = (TH1F*)qafile->Get("h_MvtxRawHitQA_strobe_bc");
  int chip_bc_entries = h_MvtxRawHitQA_chip_bc->GetEntries();
  int strobe_bc_entries = h_MvtxRawHitQA_strobe_bc->GetEntries();
  int chip_bins = std::count_if(
      1 + h_MvtxRawHitQA_chip_bc->GetArray(),
      h_MvtxRawHitQA_chip_bc->GetArray() + h_MvtxRawHitQA_chip_bc->GetNbinsX(),
      [](double content) { return content > 0; }
  );
  int strobe_bins = std::count_if(
      1 + h_MvtxRawHitQA_strobe_bc->GetArray(),
      h_MvtxRawHitQA_strobe_bc->GetArray() + h_MvtxRawHitQA_strobe_bc->GetNbinsX(),
      [](double content) { return content > 0; }
  );

  float chip_bc_expected = chip_bc_entries/(float)chip_bins;
  float strobe_bc_expected = strobe_bc_entries/(float)strobe_bins;

  float chi2_chip_bc = 0.;
  float chi2_strobe_bc = 0.;
  for (int ichip = 1; ichip < h_MvtxRawHitQA_chip_bc->GetNbinsX(); ichip++)
  {
    if (h_MvtxRawHitQA_chip_bc->GetBinContent(ichip) == 0){continue;}
    float chip_bc_observed = h_MvtxRawHitQA_chip_bc->GetBinContent(ichip);
    
    chi2_chip_bc += TMath::Power(chip_bc_observed - chip_bc_expected,2)/chip_bc_expected;
  }
  for (int istrobe = 1; istrobe < h_MvtxRawHitQA_strobe_bc->GetNbinsX(); istrobe++)
  {
    if (h_MvtxRawHitQA_strobe_bc->GetBinContent(istrobe) == 0){continue;}
    float strobe_bc_observed = h_MvtxRawHitQA_strobe_bc->GetBinContent(istrobe);
    chi2_strobe_bc += TMath::Power(strobe_bc_observed - strobe_bc_expected,2)/strobe_bc_expected;
  }

  reduced_chi2_chip = chi2_chip_bc/chip_bins;
  reduced_chi2_strobe = chi2_strobe_bc/strobe_bins;

  for (int i = 0; i < 3; i++)
  {
    sprintf(hmvtx,"h_MvtxRawHitQA_nhits_stave_chip_layer%d",i);
    h_MvtxRawHitQA_nhits_stave_chip[i] = (TH2F*)qafile->Get(hmvtx); 

    sprintf(hmvtx,"h_MvtxRawHitQA_nhits_layer%d",i);
    h_MvtxRawHitQA_nhits_layer[i] = (TH1F*)qafile->Get(hmvtx); 
    sprintf(hmvtx,"h_MvtxRawHitQA_nhits_layer%d_%d",i,runnum);
    h_MvtxRawHitQA_nhits_layer[i]->SetName(hmvtx);

    // get peak position of hit distributions in each layer
    TF1 *ffitData = new TF1("ffitData",langaufun,100,8000,4);
	  ffitData->SetParameter(0,200);
	  ffitData->SetParameter(1,1200);
	  ffitData->SetParameter(2,1e6);
    ffitData->SetParameter(3,300);
	  ffitData->SetParLimits(0,10,500);
	  ffitData->SetParLimits(1,0,2500);
	  ffitData->SetParLimits(3,200,500);

    h_MvtxRawHitQA_nhits_layer[i]->Fit(ffitData,"QR");

    outputFile.cd();
    h_MvtxRawHitQA_nhits_layer[i]->Write();
    hitpeak.push_back(ffitData->GetParameter(1));

    // get number of dead staves or chips in each layer
    int emptystaveschips = 0;
    Float_t *bins = h_MvtxRawHitQA_nhits_stave_chip[i]->GetArray();
    for (int x = 0; x < totalchips[i]; x++)
    {
      for (int y = 0; y < totalstaves[i]; y++)
      {
        if (h_MvtxRawHitQA_nhits_stave_chip[i]->GetBinContent(x+1,y+1) == 0)
        {
          emptystaveschips++;
        }
      }
    }

    deadstavechipeff.push_back(((float)emptystaveschips)/totalstaveschips[i]);
  }

  //gIntt.Draw("ape");
  return 0;
}



void mvtxQA()
{
  // Get runlist
  std::vector<int> runlistVect;
  GetRunList(runlistVect);

  int nlayers = 3;

  TGraphErrors *gEffStavesChips[3];
  TGraphErrors *gHitPeaks[3];
  TGraphErrors *gReducedChi2Chip = new TGraphErrors();
  gReducedChi2Chip->SetName("gReducedChi2Chip");
  TGraphErrors *gReducedChi2Strobe = new TGraphErrors();
  gReducedChi2Strobe->SetName("gReducedChi2Strobe");

  for (int i = 0; i < nlayers; i++)
  {
    gEffStavesChips[i] = new TGraphErrors();
    gEffStavesChips[i]->SetName(Form("gEffStavesChips_%d",i));
    gHitPeaks[i] = new TGraphErrors();
    gHitPeaks[i]->SetName(Form("gHitPeaks_%d",i));
  }

  TFile *fHitDists = new TFile("mvtx/mvtxHitDistributions.root","RECREATE");

  std::ofstream outFile("mvtx/goodmvtxruns.txt");

  int totalruns = 0;
  int deadstavecut = 0;
  int hitpeakcut = 0;
  int bc_chi2cut = 0;

  for (const auto& run : runlistVect)
  {
    //if (run != 53879){continue;}
    if (run == 51771){continue;}
    std::vector<float> deadstavechipeff;
    std::vector<float> hitpeak;
    float reduced_chi2_chip;
    float reduced_chi2_strobe;
    int status = QA(run,*fHitDists,deadstavechipeff,hitpeak,reduced_chi2_chip,reduced_chi2_strobe);
    if (status == 1){continue;}
    gEffStavesChips[0]->SetPoint(gEffStavesChips[0]->GetN(),run,1-deadstavechipeff[0]);
    gEffStavesChips[1]->SetPoint(gEffStavesChips[1]->GetN(),run,1-deadstavechipeff[1]);
    gEffStavesChips[2]->SetPoint(gEffStavesChips[2]->GetN(),run,1-deadstavechipeff[2]);
    gHitPeaks[0]->SetPoint(gHitPeaks[0]->GetN(),run,hitpeak[0]);
    gHitPeaks[1]->SetPoint(gHitPeaks[1]->GetN(),run,hitpeak[1]);
    gHitPeaks[2]->SetPoint(gHitPeaks[2]->GetN(),run,hitpeak[2]);
    gReducedChi2Chip->SetPoint(gReducedChi2Chip->GetN(),run,reduced_chi2_chip);
    gReducedChi2Strobe->SetPoint(gReducedChi2Strobe->GetN(),run,reduced_chi2_strobe);
    std::cout << run << " " << status << std::endl;
    totalruns++;
    if (deadstavechipeff[0] < 0.05 && deadstavechipeff[1] < 0.05 && deadstavechipeff[2] < 0.05)
    {
      deadstavecut++;
      if (hitpeak[0] > 450 && hitpeak[1] > 450 && hitpeak[2] > 450)
      {
        hitpeakcut++;
        if (reduced_chi2_chip > 0 && reduced_chi2_chip < 2.5 && reduced_chi2_strobe > 0 && reduced_chi2_strobe < 2.5)
        {
          bc_chi2cut++;
          outFile << run << "\n";
        }
      }
    }
    
    
    deadstavechipeff.clear();
    hitpeak.clear();

  }

  std::cout << totalruns << " " << deadstavecut << " " << hitpeakcut << " " << bc_chi2cut << std::endl;

  outFile.close();

  fHitDists->Write();
  fHitDists->Close();
  //gEffStavesChips[0]->Draw("ape");

  TFile *fMvtxQA = new TFile("mvtx/mvtxqa.root","RECREATE");
  for (auto graph : gEffStavesChips)
  {
    graph->Write();
    delete graph;
  }
  for (auto graph : gHitPeaks)
  {
    graph->Write();
    delete graph;
  }
  gReducedChi2Chip->Write();
  gReducedChi2Strobe->Write();

  fMvtxQA->Write();
  fMvtxQA->Close();

}
