void plotintt()
{
  std::vector<int> vect;
  std::string line;
  ifstream irunlist;
  irunlist.open("../runlist_highstats.txt");
  if(!irunlist)
  {
    std::cout << "failed to open runlist " << std::endl;
    exit(1);
  }
  while (std::getline(irunlist, line)) {
    vect.push_back(std::stoi(line));
  }
  irunlist.close(); 

  TFile *infile = new TFile("../inttqa.root");

  /*
  // ***************** Draw full size example *********************** //
  TCanvas *cexamp = new TCanvas("cexamp","cexamp",800,600);
  cexamp->cd();
  TGraphErrors *gRMSexamp = (TGraphErrors*)infile->Get("gInttRMSX_51944");
  if (!gRMSexamp){exit(1);}
  TGraphErrors *gRMSYexamp = (TGraphErrors*)infile->Get("gInttRMSY_51944");
  if (!gRMSYexamp){exit(1);}
  gRMSexamp->GetYaxis()->SetTitleSize(0.08);
  gRMSYexamp->GetYaxis()->SetTitleSize(0.08);
  gRMSexamp->GetYaxis()->SetTitleOffset(0.5);
  gRMSYexamp->GetYaxis()->SetTitleOffset(0.5);
  gRMSexamp->GetYaxis()->SetRangeUser(-5,50);
  gRMSYexamp->GetYaxis()->SetRangeUser(-5,50);

  gRMSexamp->Draw("ape");
  gRMSYexamp->Draw("pesame");
  
  TAxis* xAxisexamp = gRMSexamp->GetXaxis();
  xAxisexamp->SetLimits(0, 112);
  xAxisexamp->SetRangeUser(0, 112);
  xAxisexamp->SetNdivisions(8, false);
  xAxisexamp->SetLabelSize(0); // Suppress default labels
  
  for (int i = 0; i <= 7; i++) 
  {
    double tickStart = i * 14;
    double tickEnd = (i + 1) * 14;
    double midpoint = (tickStart + tickEnd) / 2.0;
    
    // Create and customize the label
    TText* labelexamp = new TText(midpoint+1.0, -2.0, Form("intt%d", i)); // Adjust y-coordinate (-2.0)
    labelexamp->SetTextSize(0.05); // Increase text size
    labelexamp->SetTextAlign(33);  // Align top-right for better positioning
    labelexamp->SetTextAngle(45);  // Angle the text at 45 degrees
    labelexamp->Draw();
  }
  
  cexamp->Update();
  
  TLatex latexexamp;
  latexexamp.SetTextSize(0.05);
  latexexamp.SetTextAlign(13);
  latexexamp.DrawLatexNDC(0.2, 0.6, "Run Number: 51944");
  latexexamp.SetTextColor(kRed);
  latexexamp.DrawLatexNDC(0.25, 0.5, "Channel");
  latexexamp.SetTextColor(kBlue);
  latexexamp.DrawLatexNDC(0.25, 0.4, "Chip");   
  // *************************************** //
  */
    

  
  TCanvas *c[24];
  //TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
  //c1->Divide(6,6);

  std::cout << vect.size() << std::endl;

  for (int i = 0; i < vect.size(); i++)
  {
    std::cout << vect[i] << std::endl;
    if (i%36==0)
    {
      c[i/36] = new TCanvas(Form("c%d",i/36),Form("c%d",i/36),1800,1200);
      c[i/36]->Divide(6,6);
    }
    
    TGraphErrors *gRMSX = (TGraphErrors*)infile->Get(Form("gInttRMSX_%d",vect[i]));
    if (!gRMSX){continue;}
    TGraphErrors *gRMSY = (TGraphErrors*)infile->Get(Form("gInttRMSY_%d",vect[i]));
    if (!gRMSY){continue;}
    c[i/36]->cd(i%36+1);
    gRMSX->GetYaxis()->SetTitleSize(0.08);
    gRMSY->GetYaxis()->SetTitleSize(0.08);
    gRMSX->GetYaxis()->SetTitleOffset(0.5);
    gRMSY->GetYaxis()->SetTitleOffset(0.5);
    gRMSX->GetYaxis()->SetRangeUser(-5,50);
    gRMSY->GetYaxis()->SetRangeUser(-5,50);

    
    gRMSX->Draw("ape");
    gRMSY->Draw("pesame");
    
    TAxis* xAxis = gRMSX->GetXaxis();
    xAxis->SetLimits(0, 112);
    xAxis->SetRangeUser(0, 112);
    xAxis->SetNdivisions(8, false);
    xAxis->SetLabelSize(0); // Suppress default labels
    
    for (int i = 0; i <= 7; i++) 
    {
      double tickStart = i * 14;
      double tickEnd = (i + 1) * 14;
      double midpoint = (tickStart + tickEnd) / 2.0;
      
      // Create and customize the label
      TText* label = new TText(midpoint+1.0, 1.0, Form("intt%d", i)); // Adjust y-coordinate (-2.0)
      label->SetTextSize(0.1); // Increase text size
      label->SetTextAlign(33);  // Align top-right for better positioning
      label->SetTextAngle(45);  // Angle the text at 45 degrees
      label->Draw();
    }
  
    c[i/36]->Update();
    
    TLatex latexrun;
    latexrun.SetTextSize(0.1);
    latexrun.SetTextAlign(13);
    latexrun.DrawLatexNDC(0.2, 0.6, Form("Run Number: %d", vect[i]));
    latexrun.SetTextColor(kRed);
    latexrun.DrawLatexNDC(0.25, 0.5, "Channel");
    latexrun.SetTextColor(kBlue);
    latexrun.DrawLatexNDC(0.25, 0.4, "Chip"); 

    //c[i / 36]->SaveAs(Form("images/canvas_%d.png", i / 36));
  
  }
  
}
