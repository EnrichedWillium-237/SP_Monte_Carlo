// style.h
// little plotting classes to make Will's life easier 

void SetTPaveTxt( TPaveText * txtemplate, int txtsize ) {
  txtemplate->SetFillColor(0);
  txtemplate->SetBorderSize(0);
  txtemplate->SetTextFont(43);
  txtemplate->SetTextSize(txtsize);
}

void SetLegend( TLegend * legtemplate, int legsize ) {
  legtemplate->SetFillColor(0);
  legtemplate->SetBorderSize(0);
  legtemplate->SetTextFont(43);
  legtemplate->SetTextSize(legsize);
}

void GraphToHist( TGraphErrors * gin, TH1D * hout) {
  int num = gin->GetN();
  Double_t x[100], y[100], yerr[100];
  for (int i = 0; i<num; i++) {
    gin->GetPoint(i, x[i], y[i]);
    yerr[i] = gin->GetErrorY(i);
    hout->SetBinContent(i+1, y[i]);
    hout->SetBinError(i+1, yerr[i]);
  }
}
