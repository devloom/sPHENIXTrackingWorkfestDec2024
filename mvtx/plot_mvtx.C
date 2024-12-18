void plot_mvtx()
{
    TFile *f = new TFile("mvtxqa.root");
    TGraphErrors *g0 = (TGraphErrors*)f->Get("gReducedChi2Chip");
    TGraphErrors *g1 = (TGraphErrors*)f->Get("gReducedChi2Strobe");

    g0->SetMarkerStyle(20);
    g0->SetMarkerColor(kBlack);
    g0->GetXaxis()->SetTitle("Run number");
    g0->GetYaxis()->SetTitle("Bunch crossing #chi^{2}/NDF");
    g0->GetYaxis()->SetRangeUser(-100,3500);
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kRed);

    g0->Draw("ape");
    g1->Draw("pesame");

    TLatex latexrun;
    latexrun.SetTextSize(0.04);
    latexrun.SetTextAlign(13);
    latexrun.DrawLatexNDC(0.3, 0.5, "Chip BC");
    latexrun.SetTextColor(kRed);
    latexrun.DrawLatexNDC(0.3, 0.45, "Strobe BC");

    /*
    TGraphErrors *g0 = (TGraphErrors*)f->Get("gHitPeaks_0");
    TGraphErrors *g1 = (TGraphErrors*)f->Get("gHitPeaks_1");
    TGraphErrors *g2 = (TGraphErrors*)f->Get("gHitPeaks_2");

    g0->SetMarkerStyle(20);
    g0->SetMarkerColor(kBlack);
    g0->GetXaxis()->SetTitle("Run number");
    g0->GetYaxis()->SetTitle("Total raw hit MPV");
    g0->GetYaxis()->SetRangeUser(-100,2250);
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kRed);
    g2->SetMarkerStyle(20); 
    g2->SetMarkerColor(kBlue);

    g0->Draw("ape");
    g1->Draw("pesame");
    g2->Draw("pesame");

    TLatex latexrun;
    latexrun.SetTextSize(0.04);
    latexrun.SetTextAlign(13);
    latexrun.DrawLatexNDC(0.7, 0.85, "Layer 0");
    latexrun.SetTextColor(kRed);
    latexrun.DrawLatexNDC(0.7, 0.8, "Layer 1");
    latexrun.SetTextColor(kBlue);
    latexrun.DrawLatexNDC(0.7, 0.75, "Layer 2");  */


}