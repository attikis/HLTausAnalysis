void HistoExtracter()
{
  TCanvas *c1 = new TCanvas("c1","",800,800);
  TFile f1("JetProbability_Histograms_TTBar_11May001.root");

  TF1 *myexp = new TF1("myexp","[0]*exp(-x/[1])",0.1,0.2);
  myexp->SetParameters(2,2);

  //  TH1F *h1 = (TH1F*)f1.Get("hBjetTkabsd0");
  TH1F *h1 = (TH1F*)f1.Get("hGenMatchedPixP");
  h1->SetLineColor(kRed);
  //h1->SetMarkerStyle(20);
  //h1->SetMarkerColor(kRed);
  h1->Draw();
  h1->Scale(1/h1->Integral());
//  h1->Fit(myexp,"R")
  c1->Update();

  //  TH1F *h2 = (TH1F*)f1.Get("hPjetTkabsd0");
  TH1F *h2 = (TH1F*)f1.Get("hGenMatchedPixB");
  //h2->SetMarkerStyle(20);
  h2->Scale(1/h2->Integral());
  h2->Draw("same");

  c1->Update();
  TFile test("AbsD0.root","RECREATE");
  // h1->Write();
  //h2->Write();
  c1->Write();
  test.Close();
}
