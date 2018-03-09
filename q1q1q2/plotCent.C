# include "TCanvas.h"
# include "TGraphErrors.h"
# include "TH1.h"
# include "TH2.h"
# include "TLegend.h"
# include "TPaveText.h"
# include "TStyle.h"
# include <fstream>
# include <iostream>

void SetTPaveTxt( TPaveText * txtemplate, int txtsize ) {
    txtemplate->SetFillColor(0);
    txtemplate->SetBorderSize(0);
    txtemplate->SetTextFont(43);
    txtemplate->SetTextAlign(12);
    txtemplate->SetTextSize(txtsize);
}

void SetLegend( TLegend * legtemplate, int legsize ) {
    legtemplate->SetFillColor(0);
    legtemplate->SetBorderSize(0);
    legtemplate->SetTextFont(43);
    legtemplate->SetTextSize(legsize);
}

TH1D * h0;
TGraphErrors * gdata;
TGraphErrors * gtoyOdd;
TGraphErrors * gtoyOddZero;
TGraphErrors * gtoyEven;

void plotCent() {

    int ncbins = 11;
    double centval[] = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 45, 55, 65};

    double q1q2data_val[] = {0.006, 0.014, 0.021, 0.027, 0.032, 0.035, 0.037, 0.038, 0.036, 0.029, 0.02};
    double q1q2data_err[] = {0.000175587, 0.000291156, 0.00038145, 0.000526683, 0.000594694, 0.000659344, 0.000741557, 0.000713392, 0.000628951, 0.000501184, 0.000360879};
    double q1q2toyOdd_val[] = {0.00298, 0.01273, 0.02519, 0.03314, 0.04285, 0.04817, 0.05375, 0.06006, 0.06541, 0.07803, 0.10551};
    double q1q2toyOdd_err[] = {0.00217, 0.00178, 0.00203, 0.00182, 0.00265, 0.00231, 0.00209, 0.0013, 0.00167, 0.00219, 0.00211};
    double q1q2toyOddZero_val[] = {0.00833, 0.01314, 0.0179, 0.02547, 0.03601, 0.0371, 0.04235, 0.04349, 0.0443, 0.04349, 0.03428};
    double q1q2toyOddZero_err[] = {0.00256, 0.00228, 0.00227, 0.0022, 0.00211, 0.00191, 0.00189, 0.00215, 0.0025, 0.00215, 0.00245};
    double q1q2toyEven_val[] = {0.01277, 0.03364, 0.04995, 0.06268, 0.0812, 0.09138, 0.10785, 0.12489, 0.16144, 0.22623, 0.32419};
    double q1q2toyEven_err[] = {0.0023, 0.00291, 0.00142, 0.0016, 0.00238, 0.00327, 0.00186, 0.00154, 0.00184, 0.00269, 0.00173};

    gdata = new TGraphErrors(ncbins, centval, q1q2data_val, 0, q1q2data_err);
    gdata->SetMarkerStyle(20);
    gdata->SetMarkerSize(1.3);
    gdata->SetMarkerColor(kBlack);
    gdata->SetLineColor(kBlack);

    gtoyOdd = new TGraphErrors(ncbins, centval, q1q2toyOdd_val, 0, q1q2toyOdd_err);
    gtoyOdd->SetMarkerStyle(20);
    gtoyOdd->SetMarkerSize(1.3);
    gtoyOdd->SetMarkerColor(kRed);
    gtoyOdd->SetLineColor(kRed);

    gtoyOddZero = new TGraphErrors(ncbins, centval, q1q2toyOddZero_val, 0, q1q2toyOddZero_err);
    gtoyOddZero->SetMarkerStyle(20);
    gtoyOddZero->SetMarkerSize(1.3);
    gtoyOddZero->SetMarkerColor(kGreen+2);
    gtoyOddZero->SetLineColor(kGreen+2);

    gtoyEven = new TGraphErrors(ncbins, centval, q1q2toyEven_val, 0, q1q2toyEven_err);
    gtoyEven->SetMarkerStyle(20);
    gtoyEven->SetMarkerSize(1.3);
    gtoyEven->SetMarkerColor(kBlue);
    gtoyEven->SetLineColor(kBlue);


    TCanvas * c0 = new TCanvas("c0", "c0", 650, 600);
    h0 = new TH1D("h0", "", 100, 0, 75);
    h0->GetYaxis()->SetRangeUser(0, 0.1);
    h0->SetXTitle("Centrality (%)");
    h0->SetYTitle("<Q_{1}^{2}Q_{2A}^{*}/|Q_{2A}|>");
    h0->GetXaxis()->CenterTitle();
    h0->GetYaxis()->CenterTitle();
    h0->GetXaxis()->SetNdivisions(509);
    h0->GetYaxis()->SetNdivisions(508);
    h0->GetXaxis()->SetDecimals();
    h0->GetYaxis()->SetDecimals();
    h0->GetXaxis()->SetLabelSize(0.04);
    h0->GetYaxis()->SetLabelSize(0.04);
    h0->Draw();
    gdata->Draw("same p");
    gtoyOdd->Draw("same p");
    //gtoyEven->Draw("same p");
    gtoyOddZero->Draw("same p");
    TLegend * leg0 = new TLegend(0.20, 0.76, 0.40, 0.92);
    SetLegend(leg0, 22);
    leg0->AddEntry(gdata,"PbPb data","p");
    leg0->AddEntry(gtoyOdd,"Toy MC v_{1}^{odd} = max","p");
    leg0->AddEntry(gtoyOddZero,"Toy MC v_{1}^{odd} = 0","p");
    //leg0->AddEntry(gtoyEven,"Toy MC v_{1}^{even} + v_{1}^{non-flow}","p");
    leg0->Draw();
    c0->Print("q1q2_cent.pdf","pdf");
    c0->Print("q1q2_cent.png","png");

}
