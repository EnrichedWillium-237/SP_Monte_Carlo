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

    double v1in = {0.000, 0.005, 0.010, 0.015, 0.020};
    double v2in = {0.000, 0.030, 0.060, 0.090, 0.120};

    double q1q1q2_000[] = {-0.00132, 0.00308, -0.00029, -0.00055, -0.00027};
    double q1q1q2_000_err[] = {4.74055, 1.19974, 11.89683, 4.93847, 4.61481};
    double q1q1q2_030[] = {0.00844, 0.00759, 0.01145, 0.01517, 0.02197};
    double q1q1q2_030_err[] = {0.30493, 0.30238, 0.17674, 0.14561, 0.09440};
    double q1q1q2_060[] = {0.02249, , 0.02582, 0.03333, 0.04988};
    double q1q1q2_060_err[] = {0.09789, , 0.06921, 0.07048, 0.04119};
    double q1q1q2_090[] = {0.03642, 0.04188, 0.05067, 0.05756, 0.07073};
    double q1q1q2_090_err[] = {0.05602, 0.06465, 0.05271, 0.03401, 0.03487};
    double q1q1q2_120[] = {0.05501, 0.05517, , , };
    double q1q1q2_120_err[] = {0.03935, 0.03328, , , };


    h2D = new TH2D


    int numEven = 4;
    double v1inEven[numEven] = {0.005, 0.010, 0.015, 0.020};
    double q1q1q2Even[numEven] = {0.02646, 0.03366, 0.03746, 0.03099};
    double q1q1q2EvenErr[numEven] = {0.00185, 0.00198, 0.00176, 0.00395};
    gEven = new TGraphErrors(numEven, v1inEven, q1q1q2Even, 0, q1q1q2EvenErr);
    gEven->SetMarkerColor(kBlue);
    gEven->SetLineColor(kBlue);

    double q1q1q2Even1[numEven] = {0.01599, 0.01870, 0.01590, 0.02407};
    double q1q1q2EvenErr1[numEven] = {0.00158, 0.00198, 0.00196, 0.00135};
    gEven1 = new TGraphErrors(numEven, v1inEven, q1q1q2Even1, 0, q1q1q2EvenErr1);
    gEven1->SetMarkerColor(kRed);
    gEven1->SetLineColor(kRed);

    double q1q1q2Even2[numEven] = {0.04316, 0.04655, 0.05874, 0.06804};
    double q1q1q2EvenErr2[numEven] = {0.00149, 0.00245, 0.00254, 0.00226};
    gEven2 = new TGraphErrors(numEven, v1inEven, q1q1q2Even2, 0, q1q1q2EvenErr2);
    gEven2->SetMarkerColor(kGreen+2);
    gEven2->SetLineColor(kGreen+2);


    // double numOdd = 2;
    // double q1q1q2Odd0[numOdd] = {, };
    // double q1q1q2Odd0Err[numOdd] = {, };
    // gOdd0 = new TGraphErrors();


    TCanvas * c0 = new TCanvas("c0", "c0", 650, 600);
    c0->cd();
    h0 = new TH1D("h0", "", 100, 0, 0.022);
    h0->GetYaxis()->SetRangeUser(0.0, 0.05);
    h0->SetXTitle("Input v_{1}^{even}");
    h0->SetYTitle("Q_{1}^{2}Q_{2}^{*}");
    h0->GetXaxis()->CenterTitle();
    h0->GetYaxis()->CenterTitle();
    h0->GetXaxis()->SetNdivisions(509);
    h0->GetYaxis()->SetNdivisions(508);
    h0->Draw();
    gEven->Draw("same p");
    gEven1->Draw("same p");
    gEven2->Draw("same p");
    TLegend * leg0 = new TLegend(0.20, 0.74, 0.40, 0.93);
    SetLegend(leg0, 22);
    leg0->AddEntry(gEven1,"v_{2}^{in} = 0.04","p");
    leg0->AddEntry(gEven,"v_{2}^{in} = 0.07","p");
    leg0->AddEntry(gEven2,"v_{2}^{in} = 0.09","p");
    leg0->Draw();

}
