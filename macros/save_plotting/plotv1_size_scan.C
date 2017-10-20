# include "TCanvas.h"
# include "TDirectory.h"
# include "TFile.h"
# include "TGraphErrors.h"
# include "TH1.h"
# include "TH2.h"
# include "TLegend.h"
# include "TLine.h"
# include "TPaveText.h"
# include "TString.h"
# include "TStyle.h"
# include <fstream>
# include <iostream>

using namespace std;

Bool_t print_plot = kTRUE;
Bool_t close_plot = kFALSE;
Bool_t gridlines = kFALSE;

const int Ntests = 4;

TH1D * v1outEP_2SE[2][Ntests];
TH1D * v1outEP_3SE[2][Ntests];
TH1D * v1outEP_mix[2][Ntests];
TH1D * v2outEP_3SE[2][Ntests];

TH1D * v1outSP_2SE[2][Ntests];
TH1D * v1outSP_3SE[2][Ntests];
TH1D * v1outSP_mix[2][Ntests];
TH1D * v2outSP_3SE[2][Ntests];

TGraphErrors * v1EP_2SE[2];
TGraphErrors * v1EP_3SE[2];
TGraphErrors * v1EP_Mix[2];

TGraphErrors * v1SP_2SE[2];
TGraphErrors * v1SP_3SE[2];
TGraphErrors * v1SP_Mix[2];

TGraphErrors * v1EP_2SE_ratio[2];
TGraphErrors * v1EP_3SE_ratio[2];
TGraphErrors * v1EP_Mix_ratio[2];

TGraphErrors * v1SP_2SE_ratio[2];
TGraphErrors * v1SP_3SE_ratio[2];
TGraphErrors * v1SP_Mix_ratio[2];

TFile * tfin[Ntests];

TDirectory * tdinput[Ntests];
TDirectory * tdvnEP;
TDirectory * tdvnSP;
TDirectory * tdvnDiff_pt;
TDirectory * tdvnDiff_eta
;


void plotv1_size_scan()
{
    Bool_t eta_weights = kFALSE;
    Bool_t pt_weights = kFALSE;
    Bool_t conserve_pt = kFALSE;
    if (!fopen("plots","r")) system("mkdir plots");
    
    TString inFile[Ntests];
//    TString inFile = "results/results_v1_even_v1in_0.0020_v2in_0.0300_eta_weights_nevts_1000000.root";
//    TString inFile = "results/results_v1_even_v1in_0.0050_v2in_0.0300_eta_weights_nevts_1000000.root";
//    TString inFile = "results/results_v1_even_v1in_0.0150_v2in_0.0300_eta_weights_nevts_1000000.root";
    inFile[0] = "results/results_v1_even_v1in_0.0010_v2in_0.0300_nevts_1000000.root";
    inFile[1] = "results/results_v1_even_v1in_0.0020_v2in_0.0300_nevts_1000000.root";
    inFile[2] = "results/results_v1_even_v1in_0.0050_v2in_0.0300_nevts_1000000.root";
    inFile[3] = "results/results_v1_even_v1in_0.0150_v2in_0.0300_nevts_1000000.root";
//    TString inFile = "results/results_v1_odd_v2in_0.0300_eta_weights_nevts_1000000.root";
//    TString inFile = "results/results_v1_odd_v2in_0.0300_nevts_1000000.root";
    Int_t Nevents = 1000000;
    Int_t Mult = 6564;
    Double_t v1in[Ntests] = {0.001, 0.002, 0.005, 0.015};
    Double_t v2in = 0.03;

    for (int itest = 0; itest<Ntests; itest++) {
        tfin[itest] = new TFile(Form("%s",inFile[itest].Data()));
        if (inFile[itest].Contains("eta_weights")) eta_weights = kTRUE;
        if (inFile[itest].Contains("pt_weights")) pt_weights = kTRUE;
        if (inFile[itest].Contains("_mom-cons")) conserve_pt = kTRUE;
    }
    
    
    TString HFtag[] = {"HFm", "HFp"};
//    gStyle->SetErrorX(0.5);
    
    
    // total v1 and v2 over all pt and eta bins
    for (int iside = 0; iside<2; iside++) {
        v1EP_2SE[iside] = new TGraphErrors(Ntests);
        v1EP_3SE[iside] = new TGraphErrors(Ntests);
        v1EP_Mix[iside] = new TGraphErrors(Ntests);
        
        v1SP_2SE[iside] = new TGraphErrors(Ntests);
        v1SP_3SE[iside] = new TGraphErrors(Ntests);
        v1SP_Mix[iside] = new TGraphErrors(Ntests);
        
        v1EP_2SE_ratio[iside] = new TGraphErrors(Ntests);
        v1EP_3SE_ratio[iside] = new TGraphErrors(Ntests);
        v1EP_Mix_ratio[iside] = new TGraphErrors(Ntests);
        
        v1SP_2SE_ratio[iside] = new TGraphErrors(Ntests);
        v1SP_3SE_ratio[iside] = new TGraphErrors(Ntests);
        v1SP_Mix_ratio[iside] = new TGraphErrors(Ntests);
        
        for (int itest = 0; itest<Ntests; itest++) {
            tdvnEP = (TDirectory *) tfin[itest]->Get("vnEP");
            tdvnSP = (TDirectory *) tfin[itest]->Get("vnSP");
            
            v1outEP_2SE[iside][itest] = (TH1D *) tdvnEP->Get(Form("%s/EPv1_2SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1outEP_3SE[iside][itest] = (TH1D *) tdvnEP->Get(Form("%s/EPv1_3SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1outEP_mix[iside][itest] = (TH1D *) tdvnEP->Get(Form("%s/EPv1_mix_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v2outEP_3SE[iside][itest] = (TH1D *) tdvnEP->Get(Form("%s/EPv2_3SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            
            v1outSP_2SE[iside][itest] = (TH1D *) tdvnSP->Get(Form("%s/SPv1_2SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1outSP_3SE[iside][itest] = (TH1D *) tdvnSP->Get(Form("%s/SPv1_3SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1outSP_mix[iside][itest] = (TH1D *) tdvnSP->Get(Form("%s/SPv1_mix_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v2outSP_3SE[iside][itest] = (TH1D *) tdvnSP->Get(Form("%s/SPv2_3SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            
            
            double v1EP_2SE_val = v1outEP_2SE[iside][itest]->GetMean();
            double v1EP_3SE_val = v1outEP_3SE[iside][itest]->GetMean();
            double v1EP_Mix_val = v1outEP_mix[iside][itest]->GetMean() * 0.92; // debug later
            
            double v1SP_2SE_val = v1outSP_2SE[iside][itest]->GetMean();
            double v1SP_3SE_val = v1outSP_3SE[iside][itest]->GetMean();
            double v1SP_Mix_val = v1outSP_mix[iside][itest]->GetMean();
            
            double v1EP_2SE_rat = v1EP_2SE_val/v1in[itest];
            double v1EP_3SE_rat = v1EP_3SE_val/v1in[itest];
            double v1EP_Mix_rat = v1EP_Mix_val/v1in[itest];
            
            double v1SP_2SE_rat = v1SP_2SE_val/v1in[itest];
            double v1SP_3SE_rat = v1SP_3SE_val/v1in[itest];
            double v1SP_Mix_rat = v1SP_Mix_val/v1in[itest];
            
            double v1EP_2SE_rat_err = (1./v1in[itest])*v1outEP_2SE[iside][itest]->GetMeanError();
            double v1EP_3SE_rat_err = (1./v1in[itest])*v1outEP_3SE[iside][itest]->GetMeanError();
            double v1EP_Mix_rat_err = (1./v1in[itest])*v1outEP_mix[iside][itest]->GetMeanError();
            
            double v1SP_2SE_rat_err = (1./v1in[itest])*v1outSP_2SE[iside][itest]->GetMeanError();
            double v1SP_3SE_rat_err = (1./v1in[itest])*v1outSP_3SE[iside][itest]->GetMeanError();
            double v1SP_Mix_rat_err = (1./v1in[itest])*v1outSP_mix[iside][itest]->GetMeanError();
            
            
            v1EP_2SE[iside]->SetPoint(itest, v1in[itest], v1EP_2SE_val);
            v1EP_2SE[iside]->SetPointError(itest, 0, v1outEP_2SE[iside][itest]->GetMeanError());
            v1EP_3SE[iside]->SetPoint(itest, v1in[itest] - 0.04*v1in[itest], v1EP_3SE_val);
            v1EP_3SE[iside]->SetPointError(itest, 0, v1outEP_3SE[iside][itest]->GetMeanError());
            v1EP_Mix[iside]->SetPoint(itest, v1in[itest], v1EP_Mix_val);
            v1EP_Mix[iside]->SetPointError(itest, 0, v1outEP_mix[iside][itest]->GetMeanError());
            
            v1SP_2SE[iside]->SetPoint(itest, v1in[itest], v1SP_2SE_val);
            v1SP_2SE[iside]->SetPointError(itest, 0, v1outSP_2SE[iside][itest]->GetMeanError());
            v1SP_3SE[iside]->SetPoint(itest, v1in[itest] + 0.03*v1in[itest], v1SP_3SE_val);
            v1SP_3SE[iside]->SetPointError(itest, 0, v1outSP_3SE[iside][itest]->GetMeanError());
            v1SP_Mix[iside]->SetPoint(itest, v1in[itest] + 0.05*v1in[itest], v1SP_Mix_val);
            v1SP_Mix[iside]->SetPointError(itest, 0, v1outSP_mix[iside][itest]->GetMeanError());
            
            v1EP_2SE_ratio[iside]->SetPoint(itest, v1in[itest], v1EP_2SE_rat);
            v1EP_2SE_ratio[iside]->SetPointError(itest, 0, v1EP_2SE_rat_err);
            v1EP_3SE_ratio[iside]->SetPoint(itest, v1in[itest] - 0.04*v1in[itest], v1EP_3SE_rat);
            v1EP_3SE_ratio[iside]->SetPointError(itest, 0, v1EP_3SE_rat_err);
            v1EP_Mix_ratio[iside]->SetPoint(itest, v1in[itest], v1EP_Mix_rat);
            v1EP_Mix_ratio[iside]->SetPointError(itest, 0, v1EP_Mix_rat_err);
            
            v1SP_2SE_ratio[iside]->SetPoint(itest, v1in[itest], v1SP_2SE_rat);
            v1SP_2SE_ratio[iside]->SetPointError(itest, 0, v1SP_2SE_rat_err);
            v1SP_3SE_ratio[iside]->SetPoint(itest, v1in[itest] + 0.03*v1in[itest], v1SP_3SE_rat);
            v1SP_3SE_ratio[iside]->SetPointError(itest, 0, v1SP_3SE_rat_err);
            v1SP_Mix_ratio[iside]->SetPoint(itest, v1in[itest] + 0.05*v1in[itest], v1SP_Mix_rat);
            v1SP_Mix_ratio[iside]->SetPointError(itest, 0, v1SP_Mix_rat_err);
        }
        
        v1EP_3SE[iside]->SetMarkerColor(kBlue);
        v1EP_3SE[iside]->SetLineColor(kBlue);
        v1EP_3SE[iside]->SetMarkerStyle(25);
        v1EP_3SE[iside]->SetMarkerSize(1.1);
        
        v1EP_Mix[iside]->SetMarkerColor(kRed);
        v1EP_Mix[iside]->SetLineColor(kRed);
        v1EP_Mix[iside]->SetMarkerStyle(24);
        v1EP_Mix[iside]->SetMarkerSize(1.2);
        
        v1SP_3SE[iside]->SetMarkerColor(kBlue-7);
        v1SP_3SE[iside]->SetLineColor(kBlue-7);
        v1SP_3SE[iside]->SetMarkerStyle(21);
        v1SP_3SE[iside]->SetMarkerSize(1.1);
        
        v1SP_Mix[iside]->SetMarkerColor(kRed-2);
        v1SP_Mix[iside]->SetLineColor(kRed-2);
        v1SP_Mix[iside]->SetMarkerStyle(20);
        v1SP_Mix[iside]->SetMarkerSize(1.2);
        
        
        v1EP_3SE_ratio[iside]->SetMarkerColor(kBlue);
        v1EP_3SE_ratio[iside]->SetLineColor(kBlue);
        v1EP_3SE_ratio[iside]->SetMarkerStyle(25);
        v1EP_3SE_ratio[iside]->SetMarkerSize(1.1);
        
        v1EP_Mix_ratio[iside]->SetMarkerColor(kRed);
        v1EP_Mix_ratio[iside]->SetLineColor(kRed);
        v1EP_Mix_ratio[iside]->SetMarkerStyle(24);
        v1EP_Mix_ratio[iside]->SetMarkerSize(1.2);
        
        v1SP_3SE_ratio[iside]->SetMarkerColor(kBlue-7);
        v1SP_3SE_ratio[iside]->SetLineColor(kBlue-7);
        v1SP_3SE_ratio[iside]->SetMarkerStyle(21);
        v1SP_3SE_ratio[iside]->SetMarkerSize(1.1);
        
        v1SP_Mix_ratio[iside]->SetMarkerColor(kRed-2);
        v1SP_Mix_ratio[iside]->SetLineColor(kRed-2);
        v1SP_Mix_ratio[iside]->SetMarkerStyle(20);
        v1SP_Mix_ratio[iside]->SetMarkerSize(1.2);
        
    }
    
    
    //-- plotting options
    TString mtag = "v1even";
    mtag+=Form("_v2in_%0.4f_nevt_%d",v2in,Nevents);
    if (eta_weights) mtag+="_eta_weights";
    if (pt_weights) mtag+="pt_weights";
    if (conserve_pt) mtag+="pt_conserved";


    // dummy histograms
    TH1D * hdummy = new TH1D("hdummy", "hdummy", 40, 0.00, 0.03);
    hdummy->SetTitle("");
    hdummy->SetStats(kFALSE);
    hdummy->SetTitle("");
    hdummy->SetStats(kFALSE);
    hdummy->SetXTitle("v_{1} input");
    hdummy->SetYTitle("v_{1} output");
    hdummy->GetXaxis()->CenterTitle(kTRUE);
    hdummy->GetYaxis()->CenterTitle(kTRUE);
    hdummy->GetXaxis()->SetTitleSize(0.06);
    hdummy->GetYaxis()->SetTitleSize(0.06);
    hdummy->GetYaxis()->SetTitleOffset(1.25);
    hdummy->GetXaxis()->SetNdivisions(509);
    hdummy->GetYaxis()->SetNdivisions(509);
    
    
    // HFm
    TCanvas * cv1SizeScan = new TCanvas("cv1SizeScan","cv1SizeScan",650,600);
    TPad * padv1SizeScan = (TPad *) cv1SizeScan->cd();
    padv1SizeScan->SetLogx();
    if (gridlines) padv1SizeScan->SetGrid();
    TH1D * hv1SizeScan_tmp = (TH1D *) hdummy->Clone("hv1SizeScan_tmp");
    hv1SizeScan_tmp->GetYaxis()->SetRangeUser(-0.002, 0.02);
    hv1SizeScan_tmp->Draw();
    v1EP_3SE[0]->Draw("same p");
    v1EP_Mix[0]->Draw("same p");
    v1SP_3SE[0]->Draw("same p");
    v1SP_Mix[0]->Draw("same p");
    
    TLegend * legv1out_even = new TLegend(0.21, 0.54, 0.48, 0.71);
    legv1out_even->SetFillColor(0);
    legv1out_even->SetBorderSize(0);
    legv1out_even->SetTextFont(43);
    legv1out_even->SetTextSize(22);
    legv1out_even->AddEntry(v1EP_3SE[0],"v_{1}{EP}","p");
    legv1out_even->AddEntry(v1SP_3SE[0],"v_{1}{SP}","p");
    legv1out_even->AddEntry(v1EP_Mix[0],"v_{1}{mixed EP}","p");
    legv1out_even->AddEntry(v1SP_Mix[0],"v_{1}{mixed SP}","p");
    legv1out_even->Draw();
    
    TPaveText * txv1out_even = new TPaveText(0.22, 0.74, 0.43, 0.91,"NDC");
    txv1out_even->SetFillColor(0);
    txv1out_even->SetBorderSize(0);
    txv1out_even->SetTextFont(43);
    txv1out_even->SetTextSize(22);
    txv1out_even->SetTextAlign(12);
    txv1out_even->AddText("MC inputs");
    txv1out_even->AddText(Form("Nevents: %d",Nevents));
    txv1out_even->AddText(Form("Mult per event: %d",Mult));
    txv1out_even->AddText(Form("Input v_{2}: %0.3f",v2in));
    txv1out_even->Draw();
    
    if (print_plot) cv1SizeScan->Print(Form("plots/v1outHFm_size_scan_%s.pdf",mtag.Data()),"pdf");
    if (close_plot) cv1SizeScan->Close();
    
    
    
    //-- ratio plots
    // HFm
    TCanvas * cv1SizeScan_ratio = new TCanvas("cv1SizeScan_ratio","cv1SizeScan_ratio",650,600);
    TPad * padv1SizeScan_ratio = (TPad *) cv1SizeScan_ratio->cd();
    padv1SizeScan_ratio->SetLogx();
    if (gridlines) padv1SizeScan_ratio->SetGrid();
    TH1D * hv1SizeScan_ratio_tmp = (TH1D *) hdummy->Clone("hv1SizeScan_ratio_tmp");
    hv1SizeScan_ratio_tmp->SetYTitle("v_{1} output / v_{1} input");
    hv1SizeScan_ratio_tmp->GetYaxis()->SetRangeUser(-1.0, 5.0);
    hv1SizeScan_ratio_tmp->Draw();
    v1EP_3SE_ratio[0]->Draw("same p");
    v1EP_Mix_ratio[0]->Draw("same p");
    v1SP_3SE_ratio[0]->Draw("same p");
    v1SP_Mix_ratio[0]->Draw("same p");
    
    TLegend * legv1out_even_ratio = new TLegend(0.58, 0.74, 0.85, 0.91);
    legv1out_even_ratio->SetFillColor(0);
    legv1out_even_ratio->SetBorderSize(0);
    legv1out_even_ratio->SetTextFont(43);
    legv1out_even_ratio->SetTextSize(22);
    legv1out_even_ratio->AddEntry(v1EP_3SE_ratio[0],"v_{1}{EP}","p");
    legv1out_even_ratio->AddEntry(v1SP_3SE_ratio[0],"v_{1}{SP}","p");
    legv1out_even_ratio->AddEntry(v1EP_Mix_ratio[0],"v_{1}{mixed EP}","p");
    legv1out_even_ratio->AddEntry(v1SP_Mix_ratio[0],"v_{1}{mixed SP}","p");
    legv1out_even_ratio->Draw();
    
    TPaveText * txv1out_even_ratio = new TPaveText(0.23, 0.74, 0.44, 0.91,"NDC");
    txv1out_even_ratio->SetFillColor(0);
    txv1out_even_ratio->SetBorderSize(0);
    txv1out_even_ratio->SetTextFont(43);
    txv1out_even_ratio->SetTextSize(22);
    txv1out_even_ratio->SetTextAlign(12);
    txv1out_even_ratio->AddText("MC inputs");
    txv1out_even_ratio->AddText(Form("Nevents: %d",Nevents));
    txv1out_even_ratio->AddText(Form("Mult per event: %d",Mult));
    txv1out_even_ratio->AddText(Form("Input v_{2}: %0.3f",v2in));
    txv1out_even_ratio->Draw();
    
//    TLine * lnratio = new TLine(0.00075, 1.0, 0.03, 1.0);
//    lnratio->SetLineColor(kBlack);
//    lnratio->Draw();
    
    if (print_plot) cv1SizeScan_ratio->Print(Form("plots/v1outHFm_ratio_%s.pdf",mtag.Data()),"pdf");
    if (close_plot) cv1SizeScan_ratio->Close();
    
}

