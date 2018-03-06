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
TGraphErrors * gv2_030;
TGraphErrors * gv2_060;
TGraphErrors * gv2_090;
TGraphErrors * gv2_120;

void q1q1q2() {

    double v1in[] = {0.000, 0.005, 0.010, 0.015, 0.020};
    double v2in[] = {0.000, 0.030, 0.060, 0.090, 0.120};

    double q1q1q2_000[] = {0.00065, -0.00284, 0.00061, -0.00156, -0.00070};
    double q1q1q2_000_err[] = {0.00142, 0.00167, 0.00128, 0.00134, 0.00174};
    double q1q1q2_030[] = {0.00936, 0.00579, 0.00927, 0.01406, 0.02224};
    double q1q1q2_030_err[] = {0.00189, 0.00135, 0.00173, 0.00138, 0.00156};
    double q1q1q2_060[] = {0.02189, 0.02316, 0.03038, 0.03947, 0.04683};
    double q1q1q2_060_err[] = {0.00166, 0.00168, 0.00129, 0.00154, 0.00137};
    double q1q1q2_090[] = {0.03912, 0.04178, 0.04522, 0.05425, 0.07089};
    double q1q1q2_090_err[] = {0.00107, 0.00154, 0.00156, 0.00168, 0.00128};
    double q1q1q2_120[] = {0.05466, 0.06092, 0.06581, 0.07078, 0.08746};
    double q1q1q2_120_err[] = {0.00128, 0.00164, 0.00128, 0.00165, 0.00164};

    h2D = new TH2D("h2D", "h2D", 5, 0, 0.020, 5, 0, 0.12);
    h2D->SetOption("colz");
    h2D->SetXTitle("Input v_{1}");
    h2D->SetYTitle("Input v_{2}");
    h2D->GetXaxis()->SetNdivisions(509);
    h2D->GetYaxis()->SetNdivisions(506);
    h2D->GetXaxis()->SetDecimals();
    h2D->GetYaxis()->SetDecimals();
    h2D->GetXaxis()->SetLabelSize(0.04);
    h2D->GetYaxis()->SetLabelSize(0.04);
    h2D->GetXaxis()->CenterTitle();
    h2D->GetYaxis()->CenterTitle();
    for (int i = 1; i<=5; i++) {
        for (int j = 1; j<=5; j++) {
            if (j == 1) {h2D->SetBinContent(i, j, q1q1q2_000[i-1]); h2D->SetBinError(i, j, q1q1q2_000_err[i-1]);}
            if (j == 2) {h2D->SetBinContent(i, j, q1q1q2_030[i-1]); h2D->SetBinError(i, j, q1q1q2_030_err[i-1]);}
            if (j == 3) {h2D->SetBinContent(i, j, q1q1q2_060[i-1]); h2D->SetBinError(i, j, q1q1q2_060_err[i-1]);}
            if (j == 4) {h2D->SetBinContent(i, j, q1q1q2_090[i-1]); h2D->SetBinError(i, j, q1q1q2_090_err[i-1]);}
            if (j == 5) {h2D->SetBinContent(i, j, q1q1q2_120[i-1]); h2D->SetBinError(i, j, q1q1q2_120_err[i-1]);}
        }
    }


    gv2_000 = new TGraphErrors(5, v1in, q1q1q2_000, 0, q1q1q2_000_err);
    gv2_000->SetMarkerStyle(20);
    gv2_000->SetMarkerSize(1.3);
    gv2_000->SetMarkerColor(kBlack);
    gv2_000->SetLineColor(kBlack);

    gv2_030 = new TGraphErrors(5, v1in, q1q1q2_030, 0, q1q1q2_030_err);
    gv2_030->SetMarkerStyle(20);
    gv2_030->SetMarkerSize(1.3);
    gv2_030->SetMarkerColor(kRed);
    gv2_030->SetLineColor(kRed);

    gv2_060 = new TGraphErrors(5, v1in, q1q1q2_060, 0, q1q1q2_060_err);
    gv2_060->SetMarkerStyle(20);
    gv2_060->SetMarkerSize(1.3);
    gv2_060->SetMarkerColor(kBlue);
    gv2_060->SetLineColor(kBlue);

    gv2_090 = new TGraphErrors(5, v1in, q1q1q2_090, 0, q1q1q2_090_err);
    gv2_090->SetMarkerStyle(20);
    gv2_090->SetMarkerSize(1.3);
    gv2_090->SetMarkerColor(kGreen+2);
    gv2_090->SetLineColor(kGreen+2);

    gv2_120 = new TGraphErrors(5, v1in, q1q1q2_120, 0, q1q1q2_120_err);
    gv2_120->SetMarkerStyle(20);
    gv2_120->SetMarkerSize(1.3);
    gv2_120->SetMarkerColor(kMagenta);
    gv2_120->SetLineColor(kMagenta);



    TCanvas * c0 = new TCanvas("c0", "c0", 650, 600);
    h0 = new TH1D("h0", "", 100, -0.001, 0.022);
    h0->GetYaxis()->SetRangeUser(-0.01, 0.12);
    h0->SetXTitle("Input v_{1}");
    h0->SetYTitle("<Q_{1}^{2}Q_{2A}^{*}>");
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
    gv2_030->Draw("same p");
    gv2_060->Draw("same p");
    gv2_090->Draw("same p");
    gv2_120->Draw("same p");
    TLegend * leg0 = new TLegend(0.20, 0.67, 0.40, 0.92);
    SetLegend(leg0, 22);
    leg0->AddEntry(gv2_000,"Input v_{2} = 0.00","p");
    leg0->AddEntry(gv2_030,"Input v_{2} = 0.03","p");
    leg0->AddEntry(gv2_060,"Input v_{2} = 0.06","p");
    leg0->AddEntry(gv2_090,"Input v_{2} = 0.09","p");
    leg0->AddEntry(gv2_120,"Input v_{2} = 0.12","p");
    leg0->Draw();


    TCanvas * c1 = new TCanvas("c1", "c1", 650, 600);
    TPad * pad1 = (TPad *) c1->cd();
    pad1->SetRightMargin(0.15);
    h2D->Draw();
    TPaveText * tx1_000 = new TPaveText(0.2, 0.2, 0.33, 0.25, "NDC");
    SetTPaveTxt(tx1_000, 20);
    tx1_000->AddText("v_{2} = 0.00");
    tx1_000->Draw();
    TPaveText * tx1_030 = new TPaveText(0.2, 0.37, 0.33, 0.42, "NDC");
    SetTPaveTxt(tx1_030, 20);
    tx1_030->AddText("v_{2} = 0.03");
    tx1_030->Draw();
    TPaveText * tx1_060 = new TPaveText(0.2, 0.52, 0.33, 0.57, "NDC");
    SetTPaveTxt(tx1_060, 20);
    tx1_060->AddText("v_{2} = 0.06");
    tx1_060->Draw();
    TPaveText * tx1_090 = new TPaveText(0.2, 0.68, 0.33, 0.74, "NDC");
    SetTPaveTxt(tx1_090, 20);
    tx1_090->AddText("v_{2} = 0.09");
    tx1_090->Draw();
    TPaveText * tx1_120 = new TPaveText(0.2, 0.84, 0.33, 0.89, "NDC");
    SetTPaveTxt(tx1_120, 20);
    tx1_120->AddText("v_{2} = 0.12");
    tx1_120->Draw();

}
