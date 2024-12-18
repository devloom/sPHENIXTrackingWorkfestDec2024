#include <iostream>
#include <fstream>
#include <vector>
#include <string>


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

int QA(int runnum, TGraphErrors& gInttRMSX, TGraphErrors& gInttRMSY)
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
  
  TH2I *h_InttRawHitQA_intt[8][14];
  char hintt[128];
  for (int i = 0; i < 8; i++)
  {
    for (int j = 0; j < 14; j++)
    {
      sprintf(hintt,"h_InttRawHitQA_intt%d_%d",i,j);
      h_InttRawHitQA_intt[i][j] = (TH2I*)qafile->Get(hintt);
      std::cout << i << " " << j << " " << h_InttRawHitQA_intt[i][j]->GetRMS(1) << " " << h_InttRawHitQA_intt[i][j]->GetRMS(2) << std::endl;
      double effectiveEntries = h_InttRawHitQA_intt[i][j]->GetEffectiveEntries();
      double RMSX = h_InttRawHitQA_intt[i][j]->GetRMS(1);
      double errorRMSX = RMSX / sqrt(2 * effectiveEntries);
      double RMSY = h_InttRawHitQA_intt[i][j]->GetRMS(2);
      double errorRMSY = RMSY / sqrt(2 * effectiveEntries);
      gInttRMSX.SetPoint(14*i+j,14*i+j,RMSX);
      gInttRMSX.SetPointError(14*i+j,0,errorRMSX);
      gInttRMSY.SetPoint(14*i+j,14*i+j,RMSY);
      gInttRMSY.SetPointError(14*i+j,0,errorRMSY);
    }
  }

  //gIntt.Draw("ape");
  return 0;
}



void inttQA()
{
  // Get runlist
  std::vector<int> runlistVect;
  GetRunList(runlistVect);

  // output files
  std::ofstream siliconrunqa("runqa_status_highstats.txt");
  if (!siliconrunqa.is_open()) 
  {
    std::cerr << "Error: Could not open the file for writing!" << std::endl;
    return 1;
  }

  std::vector<TGraphErrors*> gInttRMSXVect;
  std::vector<TGraphErrors*> gInttRMSYVect;
  for (const auto& run : runlistVect)
  {
    if (run == 51771){continue;}
    TGraphErrors *gInttRMSX = new TGraphErrors();
    TGraphErrors *gInttRMSY = new TGraphErrors();
    int status = QA(run,*gInttRMSX,*gInttRMSY);
    gInttRMSX->SetMarkerStyle(20);
    gInttRMSY->SetMarkerStyle(24);
    gInttRMSX->SetMarkerSize(0.8);
    gInttRMSY->SetMarkerSize(0.8);
    gInttRMSX->SetMarkerColor(kBlue);
    gInttRMSY->SetMarkerColor(kRed);
    gInttRMSX->SetLineColor(kBlue);
    gInttRMSY->SetLineColor(kRed);
    gInttRMSX->GetYaxis()->SetTitle("RMS");
    gInttRMSY->GetYaxis()->SetTitle("RMS");
    gInttRMSX->GetYaxis()->SetRangeUser(0,40);
    gInttRMSY->GetYaxis()->SetRangeUser(0,40);
    gInttRMSX->SetName(Form("gInttRMSX_%d",run));
    gInttRMSY->SetName(Form("gInttRMSY_%d",run));
    gInttRMSXVect.push_back(gInttRMSX);
    gInttRMSYVect.push_back(gInttRMSY);
    std::cout << run << " " << status << std::endl;
    //siliconrunqa << run << " " << status << std::endl;  
  }

  siliconrunqa.close();

  TFile *fInttQA = new TFile("inttqa.root","RECREATE");
  for (const auto& graphX : gInttRMSXVect)
  {
    graphX->Write();
    delete graphX;
    
  }
  for (const auto& graphY : gInttRMSYVect)
  {
    graphY->Write(); 
    delete graphY;
  }
  fInttQA->Write();
  fInttQA->Close();

  gInttRMSXVect.clear();
  gInttRMSYVect.clear();


  /*
  TCanvas* canvas = new TCanvas("canvas", "Custom X-Axis", 800, 600);
  gInttRMSYVect[0]->Draw("ape");
  gInttRMSXVect[0]->Draw("pesame");

  TAxis* xAxis = gInttRMSYVect[0]->GetXaxis();
  xAxis->SetLimits(0, 112);
  xAxis->SetRangeUser(0, 112);

  xAxis->SetNdivisions(8, false);
  xAxis->SetLabelSize(0);
  
  for (int i = 0; i <= 7; i++) 
  {
    double tickStart = i * 14;
    double tickEnd = (i + 1) * 14;
    double midpoint = (tickStart + tickEnd) / 2.0;

    TText* label = new TText(midpoint, -1.5, Form("intt%d", i));
    label->SetTextSize(0.03);
    label->SetTextAlign(22);
    label->Draw();
  }
  
  canvas->Update();

  TLatex latexrun;
  latexrun.SetTextSize(0.04);
  latexrun.SetTextAlign(13);
  latexrun.DrawLatexNDC(0.2, 0.7, Form("Run Number: %d", runlistVect[0]));
  latexrun.SetTextColor(kRed);
  latexrun.DrawLatexNDC(0.25, 0.65, "Channel");
  latexrun.SetTextColor(kBlue);
  latexrun.DrawLatexNDC(0.25, 0.6, "Chip");
  */
  

}
