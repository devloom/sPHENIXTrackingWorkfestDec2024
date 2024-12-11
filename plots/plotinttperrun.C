void plotinttperrun()
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

  std::cout << vect.size() << std::endl;

  TGraphErrors *gFeeEff = new TGraphErrors();
  TGraphErrors *gFeeEff_chip = new TGraphErrors();
  TGraphErrors *gFeeEff_channel = new TGraphErrors();

  for (int i = 0; i < vect.size(); i++)
  {
    if (vect[i] != 51944){continue;}
    //std::cout << vect[i] << std::endl;    
    TGraphErrors *gRMS_Chip = (TGraphErrors*)infile->Get(Form("gInttRMSX_%d",vect[i]));
    if (!gRMS_Chip){continue;}
    TGraphErrors *gRMS_Channel = (TGraphErrors*)infile->Get(Form("gInttRMSY_%d",vect[i]));
    if (!gRMS_Channel){continue;}

    // if RMS = 0, ignore point when calculating bands (do robust LTS regression after removing RMS=0 points)

    TGraphErrors *gRMS_Chip_copy = new TGraphErrors();
    TGraphErrors *gRMS_Channel_copy = new TGraphErrors();

    for (int j = 0; j < gRMS_Chip->GetN(); j++) 
    {
      double x_chip, y_chip;
      double x_channel, y_channel;
      gRMS_Chip->GetPoint(j, x_chip, y_chip);
      gRMS_Channel->GetPoint(j, x_channel, y_channel);

      if (y_chip != 0) 
      {
        double ex_chip = gRMS_Chip->GetErrorX(j);
        double ey_chip = gRMS_Chip->GetErrorY(j);

        gRMS_Chip_copy->SetPoint(gRMS_Chip_copy->GetN(), x_chip, y_chip);
        gRMS_Chip_copy->SetPointError(gRMS_Chip_copy->GetN() - 1, ex_chip, ey_chip);
      }
      if (y_channel != 0) 
      {
        double ex_channel = gRMS_Channel->GetErrorX(j);
        double ey_channel = gRMS_Channel->GetErrorY(j);

        gRMS_Channel_copy->SetPoint(gRMS_Channel_copy->GetN(), x_channel, y_channel);
        gRMS_Channel_copy->SetPointError(gRMS_Channel_copy->GetN() - 1, ex_channel, ey_channel);
      }


    }

    // if the copy histograms are empty that means that all FEEs are dead (RMS are all zero)
    if (gRMS_Chip_copy->GetN()==0){continue;}
    if (gRMS_Channel_copy->GetN()==0){continue;}
   
    TF1 *ffit_Chip = new TF1("ffit_Chip", "pol0", 0, 112);
    TF1 *ffit_Channel = new TF1("ffit_Channel", "pol0", 0, 112);
    ffit_Chip->SetLineColor(kBlue);
    ffit_Channel->SetLineColor(kRed);
    gRMS_Chip_copy->Fit(ffit_Chip, "RSCQ rob=0.90");
    gRMS_Channel_copy->Fit(ffit_Channel, "RSCQ rob=0.90");
    float ffit_chip_par = ffit_Chip->GetParameter(0);
    float ffit_channel_par = ffit_Channel->GetParameter(0);

    //get the 75 percent of points with the closest residuals to fit
    std::vector<std::tuple<double, double, double>> pointsTrimmedVec_Chip;
    std::vector<std::tuple<double, double, double>> pointsTrimmedVec_Channel;
    for (int k = 0; k < gRMS_Chip_copy->GetN(); k++)
    {
      double x_chip_resid,y_chip_resid;
      gRMS_Chip_copy->GetPoint(k,x_chip_resid,y_chip_resid);
      double ey_chip_resid = gRMS_Chip_copy->GetErrorY(k);
      pointsTrimmedVec_Chip.push_back(std::make_tuple(x_chip_resid,y_chip_resid,ey_chip_resid));
    }
    for (int k = 0; k < gRMS_Channel_copy->GetN(); k++)
    {
      double x_channel_resid,y_channel_resid;
      gRMS_Channel_copy->GetPoint(k,x_channel_resid,y_channel_resid);
      double ey_channel_resid = gRMS_Channel_copy->GetErrorY(k);
      pointsTrimmedVec_Channel.push_back(std::make_tuple(x_channel_resid,y_channel_resid,ey_channel_resid));
    }
    
    std::sort(pointsTrimmedVec_Chip.begin(), pointsTrimmedVec_Chip.end(), [ffit_chip_par](const auto& a, const auto& b) {
        return fabs(std::get<1>(a)-ffit_chip_par) < fabs(std::get<1>(b)-ffit_chip_par);
    });

    std::sort(pointsTrimmedVec_Channel.begin(), pointsTrimmedVec_Channel.end(), [ffit_channel_par](const auto& a, const auto& b) {
        return fabs(std::get<1>(a)-ffit_channel_par) < fabs(std::get<1>(b)-ffit_channel_par);
    });

    pointsTrimmedVec_Chip.erase(pointsTrimmedVec_Chip.begin() + (pointsTrimmedVec_Chip.size() - pointsTrimmedVec_Chip.size() / 10), pointsTrimmedVec_Chip.end());
    pointsTrimmedVec_Channel.erase(pointsTrimmedVec_Channel.begin() + (pointsTrimmedVec_Channel.size() - pointsTrimmedVec_Channel.size() / 10), pointsTrimmedVec_Channel.end());

    float sumSquared_chip = 0;
    float sumSquared_channel = 0;
    TGraphErrors *gRMS_Chip_trimmed = new TGraphErrors();
    TGraphErrors *gRMS_Channel_trimmed = new TGraphErrors();
    for (const auto& tuple_chip : pointsTrimmedVec_Chip)
    {
      gRMS_Chip_trimmed->SetPoint(gRMS_Chip_trimmed->GetN(),std::get<0>(tuple_chip),std::get<1>(tuple_chip));
      gRMS_Chip_trimmed->SetPointError(gRMS_Chip_trimmed->GetN() - 1, 0,std::get<2>(tuple_chip));
      sumSquared_chip += (std::get<1>(tuple_chip)-ffit_chip_par)*(std::get<1>(tuple_chip)-ffit_chip_par);
    }
    for (const auto& tuple_channel : pointsTrimmedVec_Channel)
    {
      gRMS_Channel_trimmed->SetPoint(gRMS_Channel_trimmed->GetN(),std::get<0>(tuple_channel),std::get<1>(tuple_channel));
      gRMS_Channel_trimmed->SetPointError(gRMS_Channel_trimmed->GetN() - 1, 0,std::get<2>(tuple_channel));
      sumSquared_channel += (std::get<1>(tuple_channel)-ffit_channel_par)*(std::get<1>(tuple_channel)-ffit_channel_par);
    }

    float stddev_chip = std::sqrt(sumSquared_chip / gRMS_Chip_trimmed->GetN());
    float stddev_channel = std::sqrt(sumSquared_channel / gRMS_Channel_trimmed->GetN());

    int badFees_chip = 0;
    int badFees_channel = 0;
    for (int p = 0; p < gRMS_Chip->GetN(); p++)
    {
      double x_chip_eff, y_chip_eff;
      gRMS_Chip->GetPoint(p,x_chip_eff,y_chip_eff);
      if (y_chip_eff > ffit_chip_par+3*stddev_chip || y_chip_eff < ffit_chip_par-3*stddev_chip)
      {
	badFees_chip++;
      }
    }

    for (int q = 0; q < gRMS_Channel->GetN(); q++)
    {
      double x_channel_eff, y_channel_eff;
      gRMS_Channel->GetPoint(q,x_channel_eff,y_channel_eff);
      if (y_channel_eff > ffit_channel_par+3*stddev_channel || y_channel_eff < ffit_channel_par-3*stddev_channel)
      {
	badFees_channel++;
      }
    }

    

    //std::cout << "Chip: " << ffit_chip_par << " +/- " << stddev_chip << std::endl;
    //std::cout << "Channel: " << ffit_channel_par << " +/- " << stddev_channel << std::endl;
    float fees_eff_chip = 1-((float)badFees_chip)/gRMS_Chip->GetN();
    float fees_eff_channel = 1-((float)badFees_channel)/gRMS_Channel->GetN();
    //std::cout << "FEE efficiency (chip): " << fees_eff_chip << std::endl;
    //std::cout << "FEE efficiency (channel): " << fees_eff_channel << std::endl;

    
    TLine *upperband_chip = new TLine(0,ffit_chip_par+3*stddev_chip,112,ffit_chip_par+3*stddev_chip);
    upperband_chip->SetLineWidth(2);
    upperband_chip->SetLineColor(kBlack);
    TLine *lowerband_chip = new TLine(0,ffit_chip_par-3*stddev_chip,112,ffit_chip_par-3*stddev_chip);
    lowerband_chip->SetLineWidth(2);
    lowerband_chip->SetLineColor(kBlack);
    TLine *upperband_channel = new TLine(0,ffit_channel_par+3*stddev_channel,112,ffit_channel_par+3*stddev_channel);
    upperband_channel->SetLineWidth(2);
    upperband_channel->SetLineColor(kBlack);
    TLine *lowerband_channel = new TLine(0,ffit_channel_par-3*stddev_channel,112,ffit_channel_par-3*stddev_channel);
    lowerband_channel->SetLineWidth(2);
    lowerband_channel->SetLineColor(kBlack);



    // ********* Plot fitted FEE RMS distribution *********** //
    TCanvas* canvas = new TCanvas("canvas", "Custom X-Axis", 800, 600);
    gRMS_Chip->GetYaxis()->SetRangeUser(-5,50);
    gRMS_Chip->Draw("ape");
    ffit_Chip->Draw("same");
    upperband_chip->Draw("same");
    lowerband_chip->Draw("same");
    
    gRMS_Channel->Draw("pesame");
    ffit_Channel->Draw("same");
    upperband_channel->Draw("same");
    lowerband_channel->Draw("same");
    
    TAxis* xAxis = gRMS_Chip->GetXaxis();
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
      TText* label = new TText(midpoint+1.0, -2.0, Form("intt%d", i)); // Adjust y-coordinate (-2.0)
      label->SetTextSize(0.04); // Increase text size
      label->SetTextAlign(33);  // Align top-right for better positioning
      label->SetTextAngle(45);  // Angle the text at 45 degrees
      label->Draw();
    }
    
    canvas->Update();
    
    TLatex latexrun;
    latexrun.SetTextSize(0.04);
    latexrun.SetTextAlign(13);
    latexrun.DrawLatexNDC(0.2, 0.6, Form("Run Number: %d", vect[i]));
    latexrun.SetTextColor(kRed);
    latexrun.DrawLatexNDC(0.25, 0.55, "Channel");
    latexrun.SetTextColor(kBlue);
    latexrun.DrawLatexNDC(0.25, 0.5, "Chip"); 

    // ********************************** //



    gFeeEff_chip->SetPoint(gFeeEff_chip->GetN(),vect[i],fees_eff_chip);
    gFeeEff_channel->SetPoint(gFeeEff_channel->GetN(),vect[i],fees_eff_channel);
    gFeeEff->SetPoint(gFeeEff->GetN(),vect[i],(fees_eff_chip+fees_eff_channel)/2);

    if ((fees_eff_chip+fees_eff_channel)/2 < 0.815)
    {
      std::cout << "Bad run: " << vect[i] << std::endl;
    }
    
  }
  
  /*
  gFeeEff->GetYaxis()->SetRangeUser(0,1);
  gFeeEff->GetYaxis()->SetTitle("Good INTT FEEs");
  gFeeEff->SetMarkerStyle(20);
  gFeeEff->SetMarkerColor(kBlack);

  gFeeEff_chip->GetYaxis()->SetRangeUser(0,1);
  gFeeEff_chip->GetYaxis()->SetTitle("Good INTT FEEs");
  gFeeEff_chip->SetMarkerStyle(20);
  gFeeEff_chip->SetMarkerColor(kBlue);

  gFeeEff_channel->GetYaxis()->SetRangeUser(0,1);
  gFeeEff_channel->SetMarkerStyle(20);
  gFeeEff_channel->SetMarkerColor(kRed);

  //gFeeEff_chip->Draw("ape");
  //gFeeEff_channel->Draw("pesame");
  gFeeEff->Draw("ape");
  */



}
