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
TH2D * h2D;
TGraphErrors * gv2_000;
TGraphErrors * gv2_005;
TGraphErrors * gv2_010;
TGraphErrors * gv2_030;
TGraphErrors * gv2_060;
TGraphErrors * gv2_090;
TGraphErrors * gv2_120;

void q1q1q2() {

    double v1in[]           = {0.000,    0.005,    0.010,    0.015,    0.020,    0.033,   0.063};
    double v2in[]           = {0.000,    0.005,    0.030,    0.060,    0.090,    0.120};

    double q1q1q2_000[]     = {0.00153,  0.00273,  0.00158, -0.00223, -0.00031, -0.00219,  0.00418};
    double q1q1q2_000_err[] = {0.00282,  0.00157,  0.00194,  0.00197,  0.00147,  0.00184,  0.00248};
    double q1q1q2_005[]     = {0.00293, -0.00106,  0.00102, -0.00154,  0.00095,  0.00636,  0.01762};
    double q1q1q2_005_err[] = {0.00184,  0.00217,  0.00241,  0.00111,  0.00149,  0.00148,  0.00202};
    double q1q1q2_010[]     = {-0.00312, -0.00116, 0.00279,  0.00819,  0.00488,  0.01659,  0.04774};
    double q1q1q2_010_err[] = {0.00310,  0.00213,  0.00285,  0.00217,  0.00339,  0.00226,  0.00250};
    double q1q1q2_030[]     = {0.01070,  0.00768,  0.01102,  0.01253,  0.02106,  0.04470,  0.12605};
    double q1q1q2_030_err[] = {0.00205,  0.00124,  0.00178,  0.00205,  0.00111,  0.00184,  0.00136};
    double q1q1q2_060[]     = {0.02290,  0.02610,  0.02439,  0.03684,  0.04771,  0.08710,  0.21926};
    double q1q1q2_060_err[] = {0.00166,  0.00178,  0.00246,  0.00144,  0.00200,  0.00179,  0.00847};
    double q1q1q2_090[]     = {0.04352,  0.04240,  0.04536,  0.05561,  0.07031,  0.11309,  0.27176};
    double q1q1q2_090_err[] = {0.00180,  0.00226,  0.00236,  0.00192,  0.00185,  0.00107,  0.00926};
    double q1q1q2_120[]     = {0.05908,  0.05554,  0.06469,  0.07703,  0.08079,  0.13416,  0.30154};
    double q1q1q2_120_err[] = {0.00324,  0.00234,  0.00289,  0.00305,  0.00183,  0.00152,  0.00190};

    // h2D = new TH2D("h2D", "h2D", 5, 0, 0.020, 5, 0, 0.12);
    // h2D->SetOption("colz");
    // h2D->SetXTitle("Input v_{1}");
    // h2D->SetYTitle("Input v_{2}");
    // h2D->GetXaxis()->SetNdivisions(509);
    // h2D->GetYaxis()->SetNdivisions(506);
    // h2D->GetXaxis()->SetDecimals();
    // h2D->GetYaxis()->SetDecimals();
    // h2D->GetXaxis()->SetLabelSize(0.04);
    // h2D->GetYaxis()->SetLabelSize(0.04);
    // h2D->GetXaxis()->CenterTitle();
    // h2D->GetYaxis()->CenterTitle();
    // for (int i = 1; i<=5; i++) {
    //     for (int j = 1; j<=5; j++) {
    //         if (j == 1) {h2D->SetBinContent(i, j, q1q1q2_000[i-1]); h2D->SetBinError(i, j, q1q1q2_000_err[i-1]);}
    //         if (j == 2) {h2D->SetBinContent(i, j, q1q1q2_030[i-1]); h2D->SetBinError(i, j, q1q1q2_030_err[i-1]);}
    //         if (j == 3) {h2D->SetBinContent(i, j, q1q1q2_060[i-1]); h2D->SetBinError(i, j, q1q1q2_060_err[i-1]);}
    //         if (j == 4) {h2D->SetBinContent(i, j, q1q1q2_090[i-1]); h2D->SetBinError(i, j, q1q1q2_090_err[i-1]);}
    //         if (j == 5) {h2D->SetBinContent(i, j, q1q1q2_120[i-1]); h2D->SetBinError(i, j, q1q1q2_120_err[i-1]);}
    //     }
    // }


    gv2_000 = new TGraphErrors(7, v1in, q1q1q2_000, 0, q1q1q2_000_err);
    gv2_000->SetMarkerStyle(20);
    gv2_000->SetMarkerSize(1.4);
    gv2_000->SetMarkerColor(kBlack);
    gv2_000->SetLineColor(kBlack);

    gv2_005 = new TGraphErrors(7, v1in, q1q1q2_005, 0, q1q1q2_005_err);
    gv2_005->SetMarkerStyle(20);
    gv2_005->SetMarkerSize(1.4);
    gv2_005->SetMarkerColor(kOrange+7);
    gv2_005->SetLineColor(kOrange+7);

    gv2_010 = new TGraphErrors(7, v1in, q1q1q2_010, 0, q1q1q2_010_err);
    gv2_010->SetMarkerStyle(20);
    gv2_010->SetMarkerSize(1.4);
    gv2_010->SetMarkerColor(kCyan+2);
    gv2_010->SetLineColor(kCyan+2);

    gv2_030 = new TGraphErrors(7, v1in, q1q1q2_030, 0, q1q1q2_030_err);
    gv2_030->SetMarkerStyle(20);
    gv2_030->SetMarkerSize(1.4);
    gv2_030->SetMarkerColor(kRed);
    gv2_030->SetLineColor(kRed);

    gv2_060 = new TGraphErrors(7, v1in, q1q1q2_060, 0, q1q1q2_060_err);
    gv2_060->SetMarkerStyle(20);
    gv2_060->SetMarkerSize(1.4);
    gv2_060->SetMarkerColor(kBlue);
    gv2_060->SetLineColor(kBlue);

    gv2_090 = new TGraphErrors(7, v1in, q1q1q2_090, 0, q1q1q2_090_err);
    gv2_090->SetMarkerStyle(20);
    gv2_090->SetMarkerSize(1.4);
    gv2_090->SetMarkerColor(kGreen+2);
    gv2_090->SetLineColor(kGreen+2);

    gv2_120 = new TGraphErrors(7, v1in, q1q1q2_120, 0, q1q1q2_120_err);
    gv2_120->SetMarkerStyle(20);
    gv2_120->SetMarkerSize(1.4);
    gv2_120->SetMarkerColor(kMagenta);
    gv2_120->SetLineColor(kMagenta);



    TCanvas * c0 = new TCanvas("c0", "c0", 650, 600);
    h0 = new TH1D("h0", "", 100, -0.005, 0.07);
    h0->GetYaxis()->SetRangeUser(-0.01, 0.35);
    h0->SetXTitle("Input v_{1}");
    // h0->SetYTitle("<Q_{1}^{2}Q_{2A}^{*}>");
    h0->SetYTitle("<cos[2(#Psi_{1} - #Psi_{2})]>");
    h0->GetXaxis()->CenterTitle();
    h0->GetYaxis()->CenterTitle();
    h0->GetXaxis()->SetNdivisions(509);
    h0->GetYaxis()->SetNdivisions(508);
    h0->GetXaxis()->SetDecimals();
    h0->GetYaxis()->SetDecimals();
    h0->GetXaxis()->SetLabelSize(0.04);
    h0->GetYaxis()->SetLabelSize(0.04);
    h0->Draw();
    gv2_000->Draw("same p");
    gv2_005->Draw("same p");
    gv2_010->Draw("same p");
    gv2_030->Draw("same p");
    gv2_060->Draw("same p");
    gv2_090->Draw("same p");
    gv2_120->Draw("same p");
    TLegend * leg0 = new TLegend(0.20, 0.48, 0.40, 0.92);
    SetLegend(leg0, 24);
    leg0->AddEntry(gv2_000,"v_{2} = 0.000","p");
    leg0->AddEntry(gv2_005,"v_{2} = 0.005","p");
    leg0->AddEntry(gv2_010,"v_{2} = 0.010","p");
    leg0->AddEntry(gv2_030,"v_{2} = 0.030","p");
    leg0->AddEntry(gv2_060,"v_{2} = 0.060","p");
    leg0->AddEntry(gv2_090,"v_{2} = 0.090","p");
    leg0->AddEntry(gv2_120,"v_{2} = 0.120","p");
    leg0->Draw();
    c0->Print("figq1q12.pdf","pdf");


    // TCanvas * c1 = new TCanvas("c1", "c1", 650, 600);
    // TPad * pad1 = (TPad *) c1->cd();
    // pad1->SetRightMargin(0.15);
    // h2D->Draw();
    // TPaveText * tx1_000 = new TPaveText(0.2, 0.2, 0.33, 0.25, "NDC");
    // SetTPaveTxt(tx1_000, 20);
    // tx1_000->AddText("v_{2} = 0.00");
    // tx1_000->Draw();
    // TPaveText * tx1_030 = new TPaveText(0.2, 0.37, 0.33, 0.42, "NDC");
    // SetTPaveTxt(tx1_030, 20);
    // tx1_030->AddText("v_{2} = 0.03");
    // tx1_030->Draw();
    // TPaveText * tx1_060 = new TPaveText(0.2, 0.52, 0.33, 0.57, "NDC");
    // SetTPaveTxt(tx1_060, 20);
    // tx1_060->AddText("v_{2} = 0.06");
    // tx1_060->Draw();
    // TPaveText * tx1_090 = new TPaveText(0.2, 0.68, 0.33, 0.74, "NDC");
    // SetTPaveTxt(tx1_090, 20);
    // tx1_090->AddText("v_{2} = 0.09");
    // tx1_090->Draw();
    // TPaveText * tx1_120 = new TPaveText(0.2, 0.84, 0.33, 0.89, "NDC");
    // SetTPaveTxt(tx1_120, 20);
    // tx1_120->AddText("v_{2} = 0.12");
    // tx1_120->Draw();



}
