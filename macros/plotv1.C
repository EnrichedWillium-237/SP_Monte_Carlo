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
Bool_t close_plot = kTRUE;
Bool_t gridlines = kFALSE;

const int EPNum = 6;

TH1D * philab_sub[EPNum];
TH1D * etain_sub[EPNum];
TH1D * ptin_sub[EPNum];
TH1D * phiPsiRP[EPNum];
TH1D * phiPsi1[EPNum];
TH1D * phiPsi2[EPNum];
TH1D * hsub_Psi1lab[EPNum];
TH1D * hsub_Psi2lab[EPNum];
TH2D * pt2D_sub[EPNum];

TH1D * etain_merged;
TH1D * ptin_merged;
TH1D * phiPsiRP_merged;

TH1D * v1input;

TH1D * v1outEP_2SE[2];
TH1D * v1outEP_3SE[2];
TH1D * v1outEP_mix[2];
TH1D * v2outEP_3SE[2];

TH1D * v1outSP_2SE[2];
TH1D * v1outSP_3SE[2];
TH1D * v1outSP_mix[2];
TH1D * v2outSP_3SE[2];

TH1D * v1EP_2SE_pt[2];
TH1D * v1EP_3SE_pt[2];
TH1D * v1EP_Mix_pt[2];

TH1D * v1SP_2SE_pt[2];
TH1D * v1SP_3SE_pt[2];
TH1D * v1SP_Mix_pt[2];

TH1D * v1EP_2SE_eta[2];
TH1D * v1EP_3SE_eta[2];
TH1D * v1EP_Mix_eta[2];

TH1D * v1SP_2SE_eta[2];
TH1D * v1SP_3SE_eta[2];
TH1D * v1SP_Mix_eta[2];

TFile * tfin;

TDirectory * tdinput;
TDirectory * tdvnEP;
TDirectory * tdvnSP;
TDirectory * tdvnDiff_pt;
TDirectory * tdvnDiff_eta;

Bool_t eta_weights;
Bool_t pt_weights;
Bool_t iseven;
Bool_t isodd;
Bool_t iszero;
Bool_t conserve_pT;
Bool_t ishole;
Bool_t recentered;


void plotv1()
{
    
    gStyle->SetPalette(55);
    
    TString inFile[22];
    inFile[0] = "../results/results_v1_even_0.0000_v2_0.0700_eta_weights_hole_500000_evts.root";
    inFile[1] = "../results/results_v1_even_0.0000_v2_0.0700_eta_weights_hole_flatten_recenter_500000_evts.root";
    inFile[2] = "../results/results_v1_even_0.0000_v2_0.0700_eta_weights_mon-cons_500000_evts.root";
    inFile[3] = "../results/results_v1_even_0.0000_v2_0.0700_eta_weights_mon-cons_flatten_recenter_500000_evts.root";
    inFile[4] = "../results/results_v1_even_0.0000_v2_0.0700_eta_weights_mon-cons_hole_flatten_recenter_500000_evts.root";
    inFile[5] = "../results/results_v1_even_0.0000_v2_0.0700_eta_weights_pt_weights_mon-cons_500000_evts.root";
    inFile[6] = "../results/results_v1_even_0.0200_v2_0.0700_eta_weights_500000_evts.root";
    inFile[7] = "../results/results_v1_even_0.0200_v2_0.0700_eta_weights_hole_500000_evts.root";
    inFile[8] = "../results/results_v1_even_0.0200_v2_0.0700_eta_weights_hole_flatten_recenter_500000_evts.root";
    inFile[9] = "../results/results_v1_even_0.0200_v2_0.0700_eta_weights_mon-cons_500000_evts.root";
    inFile[10] = "../results/results_v1_even_0.0200_v2_0.0700_eta_weights_mon-cons_flatten_recenter_500000_evts.root";
    inFile[11] = "../results/results_v1_even_0.0200_v2_0.0700_eta_weights_mon-cons_hole_flatten_recenter_500000_evts.root";
    inFile[12] = "../results/results_v1_even_0.0200_v2_0.0700_eta_weights_pt_weights_mon-cons_500000_evts.root";
    inFile[13] = "../results/results_v1_even_0.0200_v2_0.0700_pt_weights_500000_evts.root";
    inFile[14] = "../results/results_v1_odd_v2_0.0700_eta_weights_500000_evts.root";
    inFile[15] = "../results/results_v1_odd_v2_0.0700_eta_weights_hole_500000_evts.root";
    inFile[16] = "../results/results_v1_odd_v2_0.0700_eta_weights_hole_flatten_recenter_500000_evts.root";
    inFile[17] = "../results/results_v1_odd_v2_0.0700_eta_weights_mon-cons_500000_evts.root";
    inFile[18] = "../results/results_v1_odd_v2_0.0700_eta_weights_mon-cons_flatten_recenter_500000_evts.root";
    inFile[19] = "../results/results_v1_odd_v2_0.0700_eta_weights_mon-cons_hole_flatten_recenter_500000_evts.root";
    inFile[20] = "../results/results_v1_odd_v2_0.0700_eta_weights_pt_weights_mon-cons_500000_evts.root";
    inFile[21] = "../results/results_v1_odd_v2_0.0700_pt_weights_500000_evts.root";
    
    Int_t Nevents = 5000000;
    Int_t Mult = 6564;
    Double_t v1in = 0.02;
    Double_t v2in = 0.07;
    int subcolor[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7};
    TString HFtag[] = {"HFm", "HFp"};
    TString EPtag[] = {"HFm", "trackm", "trackp", "HFp", "trackmid"};
    gStyle->SetErrorX(0.5);
    
    
    for (int filenum = 0; filenum<22; filenum++) {
        
        eta_weights = kFALSE;
        pt_weights = kFALSE;
        iseven = kFALSE;
        isodd = kFALSE;
        iszero = kFALSE;
        conserve_pT = kFALSE;
        ishole = kFALSE;
        recentered = kFALSE;
        
        cout << "\n Reading file "<< filenum <<": " << inFile[filenum] << "\n" << endl;
        
        TString infile = inFile[filenum];
        tfin = new TFile(Form("%s",infile.Data()));
        
        if (infile.Contains("eta_weights")) eta_weights = kTRUE;
        if (infile.Contains("even")) iseven = kTRUE;
        if (infile.Contains("odd")) isodd = kTRUE;
        if (infile.Contains("0.0000")) iszero = kTRUE;
        if (infile.Contains("pt_weights")) pt_weights = kTRUE;
        if (infile.Contains("mon-cons")) conserve_pT = kTRUE; // change later "mom-cons"
        if (infile.Contains("hole")) ishole = kTRUE;
        if (infile.Contains("recenter")) recentered = kTRUE;
        
        TString tag = "";
        if (iseven) tag+="v1even";
        if (isodd) tag+="v1odd";
        if (iszero) tag+="_zero";
        if (eta_weights) tag+="_eta-weights"; else tag+="no_eta-weights";
        if (pt_weights) tag+="_pt-weights";
        if (conserve_pT) tag+="_mom-cons";
        if (ishole) tag+="_hole"; else tag+="_perfect_acceptence";
        if (recentered) tag+="_recentered";
        if (!fopen("plots","r")) system("mkdir plots");
        if (!fopen(Form("plots/%s",tag.Data()),"r")) system(Form("mkdir plots/%s",tag.Data()));
        if (!fopen(Form("plots/%s/MCInputs",tag.Data()),"r")) system(Form("mkdir plots/%s/MCInputs",tag.Data()));
        
        tdinput = (TDirectory *) tfin->Get("Inputs");
        
        for (int nep = 0; nep<EPNum; nep++) {
            philab_sub[nep] = (TH1D *) tdinput->Get(Form("Subevent_inputs/philab_ep%d_19",nep));
            etain_sub[nep] = (TH1D *) tdinput->Get(Form("Subevent_inputs/etain_ep%d_19",nep));
            ptin_sub[nep] = (TH1D *) tdinput->Get(Form("Subevent_inputs/ptin_ep%d_19",nep));
            phiPsiRP[nep] = (TH1D *) tdinput->Get(Form("Subevent_inputs/phiPsiRP_ep%d_19",nep));
            phiPsi1[nep] = (TH1D *) tdinput->Get(Form("Subevent_inputs/phiPsi1_ep%d_19",nep));
            phiPsi2[nep] = (TH1D *) tdinput->Get(Form("Subevent_inputs/phiPsi2_ep%d_19",nep));
            hsub_Psi1lab[nep] = (TH1D *) tdinput->Get(Form("Subevent_inputs/Psi1_ep%d_19",nep));
            hsub_Psi2lab[nep] = (TH1D *) tdinput->Get(Form("Subevent_inputs/Psi2_ep%d_19",nep));
            pt2D_sub[nep] = (TH2D *) tdinput->Get(Form("Subevent_inputs/pt2D_ep%d_19",nep));
            
            philab_sub[nep]->SetMarkerColor(subcolor[nep]);
            philab_sub[nep]->SetLineColor(subcolor[nep]);
            philab_sub[nep]->SetMarkerStyle(20);
            philab_sub[nep]->SetMarkerSize(0.7);
            philab_sub[nep]->SetTitle("");
            philab_sub[nep]->SetStats(kFALSE);
            philab_sub[nep]->SetXTitle("#phi (lab frame)");
            philab_sub[nep]->GetXaxis()->CenterTitle();
            philab_sub[nep]->GetXaxis()->SetTitleSize(0.06);
            philab_sub[nep]->GetXaxis()->SetLabelSize(0.06);
            philab_sub[nep]->GetYaxis()->SetLabelSize(0.06);
            
            etain_sub[nep]->SetMarkerColor(subcolor[nep]);
            etain_sub[nep]->SetLineColor(subcolor[nep]);
            etain_sub[nep]->SetMarkerStyle(20);
            etain_sub[nep]->SetMarkerSize(0.7);
            etain_sub[nep]->SetTitle("");
            etain_sub[nep]->SetStats(kFALSE);
            etain_sub[nep]->SetXTitle("#eta");
            etain_sub[nep]->GetXaxis()->CenterTitle();
            etain_sub[nep]->GetXaxis()->SetTitleSize(0.07);
            etain_sub[nep]->GetXaxis()->SetLabelSize(0.06);
            etain_sub[nep]->GetYaxis()->SetLabelSize(0.06);
            
            ptin_sub[nep]->SetMarkerColor(subcolor[nep]);
            ptin_sub[nep]->SetLineColor(subcolor[nep]);
            ptin_sub[nep]->SetMarkerStyle(20);
            ptin_sub[nep]->SetMarkerSize(0.7);
            ptin_sub[nep]->SetTitle("");
            ptin_sub[nep]->SetStats(kFALSE);
            ptin_sub[nep]->SetXTitle("p_{T} (GeV/c)");
            ptin_sub[nep]->GetXaxis()->CenterTitle();
            ptin_sub[nep]->GetXaxis()->SetTitleSize(0.06);
            ptin_sub[nep]->GetXaxis()->SetLabelSize(0.06);
            ptin_sub[nep]->GetYaxis()->SetLabelSize(0.06);
            
            phiPsiRP[nep]->SetMarkerColor(subcolor[nep]);
            phiPsiRP[nep]->SetLineColor(subcolor[nep]);
            phiPsiRP[nep]->SetMarkerStyle(20);
            phiPsiRP[nep]->SetMarkerSize(0.7);
            phiPsiRP[nep]->SetTitle("");
            phiPsiRP[nep]->SetStats(kFALSE);
            phiPsiRP[nep]->SetXTitle("#phi - #Psi_{RP}");
            phiPsiRP[nep]->GetXaxis()->CenterTitle();
            phiPsiRP[nep]->GetXaxis()->SetTitleSize(0.06);
            phiPsiRP[nep]->GetXaxis()->SetLabelSize(0.06);
            phiPsiRP[nep]->GetYaxis()->SetLabelSize(0.06);
            
            phiPsi1[nep]->SetMarkerColor(subcolor[nep]);
            phiPsi1[nep]->SetLineColor(subcolor[nep]);
            phiPsi1[nep]->SetMarkerStyle(20);
            phiPsi1[nep]->SetMarkerSize(0.7);
            phiPsi1[nep]->SetTitle("");
            phiPsi1[nep]->SetStats(kFALSE);
            phiPsi1[nep]->SetXTitle("#phi - #Psi_{1}");
            phiPsi1[nep]->GetXaxis()->CenterTitle();
            phiPsi1[nep]->GetXaxis()->SetTitleSize(0.06);
            phiPsi1[nep]->GetXaxis()->SetLabelSize(0.06);
            phiPsi1[nep]->GetYaxis()->SetLabelSize(0.06);
            
            phiPsi2[nep]->SetMarkerColor(subcolor[nep]);
            phiPsi2[nep]->SetLineColor(subcolor[nep]);
            phiPsi2[nep]->SetMarkerStyle(20);
            phiPsi2[nep]->SetMarkerSize(0.7);
            phiPsi2[nep]->SetTitle("");
            phiPsi2[nep]->SetStats(kFALSE);
            phiPsi2[nep]->SetXTitle("#phi - #Psi_{2}");
            phiPsi2[nep]->GetXaxis()->CenterTitle();
            phiPsi2[nep]->GetXaxis()->SetTitleSize(0.06);
            phiPsi2[nep]->GetXaxis()->SetLabelSize(0.06);
            phiPsi2[nep]->GetYaxis()->SetLabelSize(0.06);
            
            hsub_Psi1lab[nep]->SetMarkerColor(subcolor[nep]);
            hsub_Psi1lab[nep]->SetLineColor(subcolor[nep]);
            hsub_Psi1lab[nep]->SetMarkerStyle(20);
            hsub_Psi1lab[nep]->SetMarkerSize(0.7);
            hsub_Psi1lab[nep]->SetTitle("");
            hsub_Psi1lab[nep]->SetStats(kFALSE);
            hsub_Psi1lab[nep]->SetXTitle("#Psi_{1} (lab frame)");
            hsub_Psi1lab[nep]->GetXaxis()->CenterTitle();
            hsub_Psi1lab[nep]->GetXaxis()->SetTitleSize(0.06);
            hsub_Psi1lab[nep]->GetXaxis()->SetLabelSize(0.06);
            hsub_Psi1lab[nep]->GetYaxis()->SetLabelSize(0.06);
            
            hsub_Psi2lab[nep]->SetMarkerColor(subcolor[nep]);
            hsub_Psi2lab[nep]->SetLineColor(subcolor[nep]);
            hsub_Psi2lab[nep]->SetMarkerStyle(20);
            hsub_Psi2lab[nep]->SetMarkerSize(0.7);
            hsub_Psi2lab[nep]->SetTitle("");
            hsub_Psi2lab[nep]->SetStats(kFALSE);
            hsub_Psi2lab[nep]->SetXTitle("#Psi_{2} (lab frame)");
            hsub_Psi2lab[nep]->GetXaxis()->CenterTitle();
            hsub_Psi2lab[nep]->GetXaxis()->SetTitleSize(0.06);
            hsub_Psi2lab[nep]->GetXaxis()->SetLabelSize(0.06);
            hsub_Psi2lab[nep]->GetYaxis()->SetLabelSize(0.06);
            
            pt2D_sub[nep]->SetOption("colz");
            pt2D_sub[nep]->SetTitle("");
            pt2D_sub[nep]->SetStats(0);
            pt2D_sub[nep]->GetXaxis()->SetRangeUser(-200,200);
            pt2D_sub[nep]->GetYaxis()->SetRangeUser(-200,200);
            pt2D_sub[nep]->SetXTitle("#Sigma cos(p_{T}) (GeV/c)");
            pt2D_sub[nep]->SetYTitle("#Sigma sin(p_{T}) (GeV/c)");
        }
        
        v1input = (TH1D *) tdinput->Get("v1in");
        etain_merged = (TH1D *) tdinput->Get("Merged_inputs/etain_19");
        ptin_merged = (TH1D *) tdinput->Get("Merged_inputs/ptin_19");
        phiPsiRP_merged = (TH1D *) tdinput->Get("Merged_inputs/phiPsiRP_19");
        
        
        TH1D * hdummy_eta = new TH1D("hdummy_eta", "hdummy_eta", 40, -6, 6);
        hdummy_eta->SetTitle("");
        hdummy_eta->SetStats(kFALSE);
        hdummy_eta->SetTitle("");
        hdummy_eta->SetStats(kFALSE);
        hdummy_eta->SetXTitle("#eta");
        hdummy_eta->SetYTitle("v_{1}");
        hdummy_eta->GetXaxis()->CenterTitle(kTRUE);
        hdummy_eta->GetYaxis()->CenterTitle(kTRUE);
        hdummy_eta->GetXaxis()->SetRangeUser(-2.4, 2.4);
        hdummy_eta->GetXaxis()->SetTitleSize(0.06);
        hdummy_eta->GetYaxis()->SetTitleSize(0.06);
        hdummy_eta->GetYaxis()->SetTitleOffset(1.25);
        
        TH1D * hdummy_pt = new TH1D("hdummy_pt", "hdummy_pt", 40, 0, 8.0);
        hdummy_pt->SetTitle("");
        hdummy_pt->SetStats(kFALSE);
        hdummy_pt->SetTitle("");
        hdummy_pt->SetStats(kFALSE);
        hdummy_pt->SetXTitle("p_{T} (GeV/c)");
        hdummy_pt->SetYTitle("v_{1}");
        hdummy_pt->GetXaxis()->CenterTitle(kTRUE);
        hdummy_pt->GetYaxis()->CenterTitle(kTRUE);
        hdummy_pt->GetXaxis()->SetRangeUser(0.0, 8.0);
        hdummy_pt->GetXaxis()->SetTitleSize(0.06);
        hdummy_pt->GetYaxis()->SetTitleSize(0.06);
        hdummy_pt->GetYaxis()->SetTitleOffset(1.25);
        
        
        v1input->SetMarkerColor(kMagenta);
        v1input->SetLineColor(kMagenta);
        v1input->SetLineStyle(2);
        v1input->SetLineWidth(2);
        v1input->SetMarkerStyle(24);
        v1input->SetMarkerSize(0.1);
        
        etain_merged->SetTitle("");
        etain_merged->SetStats(kFALSE);
        etain_merged->SetXTitle("#eta");
        etain_merged->GetXaxis()->CenterTitle(kTRUE);
        
        ptin_merged->SetTitle("");
        ptin_merged->SetStats(kFALSE);
        ptin_merged->SetXTitle("p_{T} (GeV/c)");
        ptin_merged->GetXaxis()->CenterTitle(kTRUE);
        ptin_merged->GetXaxis()->SetTitleSize(0.07);
        ptin_merged->GetXaxis()->SetLabelSize(0.06);
        ptin_merged->GetYaxis()->SetLabelSize(0.06);
        
        phiPsiRP_merged->SetTitle("");
        phiPsiRP_merged->SetStats(kFALSE);
        phiPsiRP_merged->SetXTitle("#phi - #Psi_{RP}");
        phiPsiRP_merged->GetXaxis()->CenterTitle(kTRUE);
        phiPsiRP_merged->GetXaxis()->SetTitleSize(0.07);
        phiPsiRP_merged->GetXaxis()->SetLabelSize(0.06);
        phiPsiRP_merged->GetYaxis()->SetLabelSize(0.06);
        
        
        // total v1 and v2 over all pt and eta bins
        tdvnEP = (TDirectory *) tfin->Get("v1EP");
        tdvnSP = (TDirectory *) tfin->Get("v1SP");
        for (int iside = 0; iside<2; iside++) {
            v1outEP_2SE[iside] = (TH1D *) tdvnEP->Get(Form("%s/EPv1_2SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1outEP_3SE[iside] = (TH1D *) tdvnEP->Get(Form("%s/EPv1_3SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1outEP_mix[iside] = (TH1D *) tdvnEP->Get(Form("%s/EPv1_mix_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v2outEP_3SE[iside] = (TH1D *) tdvnEP->Get(Form("%s/EPv2_3SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            
            v1outSP_2SE[iside] = (TH1D *) tdvnSP->Get(Form("%s/SPv1_2SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1outSP_3SE[iside] = (TH1D *) tdvnSP->Get(Form("%s/SPv1_3SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1outSP_mix[iside] = (TH1D *) tdvnSP->Get(Form("%s/SPv1_mix_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v2outSP_3SE[iside] = (TH1D *) tdvnSP->Get(Form("%s/SPv2_3SE_%s",HFtag[iside].Data(),HFtag[iside].Data()));
        }
        
        // differential v1
        tdvnDiff_pt = (TDirectory *) tfin->Get("vn_pt");
        tdvnDiff_eta = (TDirectory *) tfin->Get("vn_eta");
        for (int iside = 0; iside<2; iside++) {
            v1EP_2SE_pt[iside] = (TH1D *) tdvnDiff_pt->Get(Form("%s/EPv1_2SE_pt_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1EP_3SE_pt[iside] = (TH1D *) tdvnDiff_pt->Get(Form("%s/EPv1_3SE_pt_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1EP_Mix_pt[iside] = (TH1D *) tdvnDiff_pt->Get(Form("%s/EPv1_mix_pt_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            
            v1SP_2SE_pt[iside] = (TH1D *) tdvnDiff_pt->Get(Form("%s/SPv1_2SE_pt_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1SP_3SE_pt[iside] = (TH1D *) tdvnDiff_pt->Get(Form("%s/SPv1_3SE_pt_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1SP_Mix_pt[iside] = (TH1D *) tdvnDiff_pt->Get(Form("%s/SPv1_mix_pt_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            
            v1EP_2SE_eta[iside] = (TH1D *) tdvnDiff_eta->Get(Form("%s/EPv1_2SE_eta_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1EP_3SE_eta[iside] = (TH1D *) tdvnDiff_eta->Get(Form("%s/EPv1_3SE_eta_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1EP_Mix_eta[iside] = (TH1D *) tdvnDiff_eta->Get(Form("%s/EPv1_mix_eta_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            
            v1SP_2SE_eta[iside] = (TH1D *) tdvnDiff_eta->Get(Form("%s/SPv1_2SE_eta_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1SP_3SE_eta[iside] = (TH1D *) tdvnDiff_eta->Get(Form("%s/SPv1_3SE_eta_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            v1SP_Mix_eta[iside] = (TH1D *) tdvnDiff_eta->Get(Form("%s/SPv1_mix_eta_%s",HFtag[iside].Data(),HFtag[iside].Data()));
            
            
            v1EP_3SE_pt[iside]->SetMarkerColor(kBlue);
            v1EP_3SE_pt[iside]->SetLineColor(kBlue);
            v1EP_3SE_pt[iside]->SetMarkerStyle(25);
            v1EP_3SE_pt[iside]->SetMarkerSize(1.1);
            
            v1EP_Mix_pt[iside]->SetMarkerColor(kRed);
            v1EP_Mix_pt[iside]->SetLineColor(kRed);
            v1EP_Mix_pt[iside]->SetMarkerStyle(24);
            v1EP_Mix_pt[iside]->SetMarkerSize(1.2);
            
            v1SP_3SE_pt[iside]->SetMarkerColor(kBlue-7);
            v1SP_3SE_pt[iside]->SetLineColor(kBlue-7);
            v1SP_3SE_pt[iside]->SetMarkerStyle(21);
            v1SP_3SE_pt[iside]->SetMarkerSize(1.1);
            
            v1SP_Mix_pt[iside]->SetMarkerColor(kRed-2);
            v1SP_Mix_pt[iside]->SetLineColor(kRed-2);
            v1SP_Mix_pt[iside]->SetMarkerStyle(20);
            v1SP_Mix_pt[iside]->SetMarkerSize(1.2);
            
            
            v1EP_3SE_eta[iside]->SetMarkerColor(kBlue);
            v1EP_3SE_eta[iside]->SetLineColor(kBlue);
            v1EP_3SE_eta[iside]->SetMarkerStyle(25);
            v1EP_3SE_eta[iside]->SetMarkerSize(1.1);
            
            v1EP_Mix_eta[iside]->SetMarkerColor(kRed);
            v1EP_Mix_eta[iside]->SetLineColor(kRed);
            v1EP_Mix_eta[iside]->SetMarkerStyle(24);
            v1EP_Mix_eta[iside]->SetMarkerSize(1.2);
            
            v1SP_3SE_eta[iside]->SetMarkerColor(kBlue-7);
            v1SP_3SE_eta[iside]->SetLineColor(kBlue-7);
            v1SP_3SE_eta[iside]->SetMarkerStyle(21);
            v1SP_3SE_eta[iside]->SetMarkerSize(1.1);
            
            v1SP_Mix_eta[iside]->SetMarkerColor(kRed-2);
            v1SP_Mix_eta[iside]->SetLineColor(kRed-2);
            v1SP_Mix_eta[iside]->SetMarkerStyle(20);
            v1SP_Mix_eta[iside]->SetMarkerSize(1.2);
        }
        
        
        //-- plotting options
        
        //-- lines for origin
        TLine * lnv1origin_eta;
        lnv1origin_eta = new TLine(-2.4, 0, 2.4, 0);
        lnv1origin_eta->SetLineColor(kBlack);
        lnv1origin_eta->SetLineStyle(1);
        lnv1origin_eta->SetLineWidth(1);
        
        TLine * lnv1origin_pt;
        lnv1origin_pt = new TLine(0, 0, 8.0, 0);
        lnv1origin_pt->SetLineColor(kBlack);
        lnv1origin_pt->SetLineStyle(1);
        lnv1origin_pt->SetLineWidth(1);
        
        
        
        //-- plot MC inputs
        TCanvas * cphilab_sub[EPNum];
        TPaveText * txphilab_sub[EPNum];
        TCanvas * cetain_sub[EPNum];
        TPaveText * txetain_sub[EPNum];
        TCanvas * cphiPsiRP[EPNum];
        TPaveText * txphiPsiRP[EPNum];
        TCanvas * cphiPsi1[EPNum];
        TPaveText * txphiPsi1[EPNum];
        TCanvas * cphiPsi2[EPNum];
        TPaveText * txphiPsi2[EPNum];
        TCanvas * csub_Psi1lab[EPNum];
        TPaveText * txsub_Psi1lab[EPNum];
        TCanvas * csub_Psi2lab[EPNum];
        TPaveText * txsub_Psi2lab[EPNum];
        TCanvas * cpt2D_sub[EPNum];
        TPaveText * txpt2D_sub[EPNum];
        for (int nep = 0; nep<EPNum-1; nep++) {
            cphilab_sub[nep] = new TCanvas(Form("cphilab_sub_%d",nep),"cphilab_sub",600,530);
            cphilab_sub[nep]->cd();
            philab_sub[nep]->Draw();
            txphilab_sub[nep] = new TPaveText(0.24, 0.65, 0.38, 0.72,"NDC");
            txphilab_sub[nep]->SetFillColor(0);
            txphilab_sub[nep]->SetBorderSize(0);
            txphilab_sub[nep]->SetTextFont(43);
            txphilab_sub[nep]->SetTextSize(22);
            txphilab_sub[nep]->SetTextAlign(12);
            txphilab_sub[nep]->AddText(Form("%s",EPtag[nep].Data()));
            txphilab_sub[nep]->Draw();
            if (print_plot) cphilab_sub[nep]->Print(Form("plots/%s/MCInputs/philab_%s.pdf",tag.Data(),EPtag[nep].Data()),"pdf");
            if (close_plot) cphilab_sub[nep]->Close();
            
            
            cetain_sub[nep] = new TCanvas(Form("cetain_sub_%d",nep),"cetain_sub",600,530);
            cetain_sub[nep]->cd();
            etain_sub[nep]->Draw();
            txetain_sub[nep] = new TPaveText(0.24, 0.65, 0.38, 0.72,"NDC");
            txetain_sub[nep]->SetFillColor(0);
            txetain_sub[nep]->SetBorderSize(0);
            txetain_sub[nep]->SetTextFont(43);
            txetain_sub[nep]->SetTextSize(22);
            txetain_sub[nep]->SetTextAlign(12);
            txetain_sub[nep]->AddText(Form("%s",EPtag[nep].Data()));
            txetain_sub[nep]->Draw();
            if (print_plot) cetain_sub[nep]->Print(Form("plots/%s/MCInputs/etain_%s.pdf",tag.Data(),EPtag[nep].Data()),"pdf");
            if (close_plot) cetain_sub[nep]->Close();
            
            
            cphiPsiRP[nep] = new TCanvas(Form("cphiPsiRP_%d",nep),"phiPsiRP",600,530);
            cphiPsiRP[nep]->cd();
            phiPsiRP[nep]->Draw();
            txphiPsiRP[nep] = new TPaveText(0.24, 0.65, 0.38, 0.72,"NDC");
            txphiPsiRP[nep]->SetFillColor(0);
            txphiPsiRP[nep]->SetBorderSize(0);
            txphiPsiRP[nep]->SetTextFont(43);
            txphiPsiRP[nep]->SetTextSize(22);
            txphiPsiRP[nep]->SetTextAlign(12);
            txphiPsiRP[nep]->AddText(Form("%s",EPtag[nep].Data()));
            txphiPsiRP[nep]->Draw();
            if (print_plot) cphiPsiRP[nep]->Print(Form("plots/%s/MCInputs/phiPsiRP_%s.pdf",tag.Data(),EPtag[nep].Data()),"pdf");
            if (close_plot) cphiPsiRP[nep]->Close();
            
            
            cphiPsi1[nep] = new TCanvas(Form("phiPsi1_%d",nep),"phiPsi1",600,530);
            cphiPsi1[nep]->cd();
            phiPsi1[nep]->Draw();
            txphiPsi1[nep] = new TPaveText(0.24, 0.65, 0.38, 0.72,"NDC");
            txphiPsi1[nep]->SetFillColor(0);
            txphiPsi1[nep]->SetBorderSize(0);
            txphiPsi1[nep]->SetTextFont(43);
            txphiPsi1[nep]->SetTextSize(22);
            txphiPsi1[nep]->SetTextAlign(12);
            txphiPsi1[nep]->AddText(Form("%s",EPtag[nep].Data()));
            txphiPsi1[nep]->Draw();
            if (print_plot) cphiPsi1[nep]->Print(Form("plots/%s/MCInputs/phiPsi1_%s.pdf",tag.Data(),EPtag[nep].Data()),"pdf");
            if (close_plot) cphiPsi1[nep]->Close();
            
            
            cphiPsi2[nep] = new TCanvas(Form("phiPsi2_%d",nep),"phiPsi2",600,530);
            cphiPsi2[nep]->cd();
            phiPsi2[nep]->Draw();
            txphiPsi2[nep] = new TPaveText(0.24, 0.65, 0.38, 0.72,"NDC");
            txphiPsi2[nep]->SetFillColor(0);
            txphiPsi2[nep]->SetBorderSize(0);
            txphiPsi2[nep]->SetTextFont(43);
            txphiPsi2[nep]->SetTextSize(22);
            txphiPsi2[nep]->SetTextAlign(12);
            txphiPsi2[nep]->AddText(Form("%s",EPtag[nep].Data()));
            txphiPsi2[nep]->Draw();
            if (print_plot) cphiPsi2[nep]->Print(Form("plots/%s/MCInputs/phiPsi2_%s.pdf",tag.Data(),EPtag[nep].Data()),"pdf");
            if (close_plot) cphiPsi2[nep]->Close();
            
            
            csub_Psi1lab[nep] = new TCanvas(Form("sub_Psi1lab_%d",nep),"sub_Psi1lab",600,530);
            csub_Psi1lab[nep]->cd();
            hsub_Psi1lab[nep]->Draw();
            txsub_Psi1lab[nep] = new TPaveText(0.24, 0.65, 0.38, 0.72,"NDC");
            txsub_Psi1lab[nep]->SetFillColor(0);
            txsub_Psi1lab[nep]->SetBorderSize(0);
            txsub_Psi1lab[nep]->SetTextFont(43);
            txsub_Psi1lab[nep]->SetTextSize(22);
            txsub_Psi1lab[nep]->SetTextAlign(12);
            txsub_Psi1lab[nep]->AddText(Form("%s",EPtag[nep].Data()));
            txsub_Psi1lab[nep]->Draw();
            if (print_plot) csub_Psi1lab[nep]->Print(Form("plots/%s/MCInputs/Psi1lab_%s.pdf",tag.Data(),EPtag[nep].Data()),"pdf");
            if (close_plot) csub_Psi1lab[nep]->Close();
            
            
            csub_Psi2lab[nep] = new TCanvas(Form("sub_Psi2lab_%d",nep),"sub_Psi2lab",600,530);
            csub_Psi2lab[nep]->cd();
            hsub_Psi2lab[nep]->Draw();
            txsub_Psi2lab[nep] = new TPaveText(0.24, 0.65, 0.38, 0.72,"NDC");
            txsub_Psi2lab[nep]->SetFillColor(0);
            txsub_Psi2lab[nep]->SetBorderSize(0);
            txsub_Psi2lab[nep]->SetTextFont(43);
            txsub_Psi2lab[nep]->SetTextSize(22);
            txsub_Psi2lab[nep]->SetTextAlign(12);
            txsub_Psi2lab[nep]->AddText(Form("%s",EPtag[nep].Data()));
            txsub_Psi2lab[nep]->Draw();
            if (print_plot) csub_Psi2lab[nep]->Print(Form("plots/%s/MCInputs/Psi2lab_%s.pdf",tag.Data(),EPtag[nep].Data()),"pdf");
            if (close_plot) csub_Psi2lab[nep]->Close();
            
            
            cpt2D_sub[nep] = new TCanvas(Form("cpt2D_sub_%d",nep),"pt2D_sub",600,530);
            cpt2D_sub[nep]->cd();
            pt2D_sub[nep]->Draw();
            txpt2D_sub[nep] = new TPaveText(0.24, 0.65, 0.38, 0.72,"NDC");
            txpt2D_sub[nep]->SetFillColor(0);
            txpt2D_sub[nep]->SetBorderSize(0);
            txpt2D_sub[nep]->SetTextFont(43);
            txpt2D_sub[nep]->SetTextSize(22);
            txpt2D_sub[nep]->SetTextAlign(12);
            txpt2D_sub[nep]->AddText(Form("%s",EPtag[nep].Data()));
            txpt2D_sub[nep]->Draw();
            if (print_plot) cpt2D_sub[nep]->Print(Form("plots/%s/MCInputs/pt2D_%s.pdf",tag.Data(),EPtag[nep].Data()),"pdf");
            if (close_plot) cpt2D_sub[nep]->Close();
        }
        
        
//        TCanvas * cMCinputs_sub = new TCanvas("cMCinputs_sub","cMCinputs_sub",800,750);
//        cMCinputs_sub->Divide(2,2);
//        
//        cMCinputs_sub->cd(1);
//        TPaveText * txMCinputs_sub_inputs = new TPaveText(0.16, 0.56, 0.75, 0.91,"NDC");
//        txMCinputs_sub_inputs->SetFillColor(0);
//        txMCinputs_sub_inputs->SetBorderSize(0);
//        txMCinputs_sub_inputs->SetTextFont(43);
//        txMCinputs_sub_inputs->SetTextSize(20);
//        txMCinputs_sub_inputs->SetTextAlign(12);
//        txMCinputs_sub_inputs->AddText("MC inputs");
//        txMCinputs_sub_inputs->AddText(Form("Nevents: %d",Nevents));
//        txMCinputs_sub_inputs->AddText(Form("Mult per event: %d",Mult));
//        if (iseven) txMCinputs_sub_inputs->AddText(Form("Input v_{1}: %0.3f",v1in));
//        txMCinputs_sub_inputs->AddText(Form("Input v_{2}: %0.3f",v2in));
//        txMCinputs_sub_inputs->Draw();
//        
//        TLegend * legMCinputs_sub = new TLegend(0.16, 0.16, 0.60, 0.46);
//        legMCinputs_sub->SetFillColor(0);
//        legMCinputs_sub->SetBorderSize(0);
//        legMCinputs_sub->SetTextFont(43);
//        legMCinputs_sub->SetTextSize(20);
//        legMCinputs_sub->SetTextAlign(12);
//        legMCinputs_sub->SetHeader("Subevent windows");
//        legMCinputs_sub->AddEntry(etain_sub[0],"HFm    (-5 < #eta < -3)","lp");
//        legMCinputs_sub->AddEntry(etain_sub[1],"trackm (-2.4 < #eta < 0)","lp");
//        legMCinputs_sub->AddEntry(etain_sub[2],"trackp (0 < #eta < 2.4)","lp");
//        legMCinputs_sub->AddEntry(etain_sub[3],"HFp    (3 < #eta < 5)","lp");
//        legMCinputs_sub->Draw();
//        
//        cMCinputs_sub->cd(2);
//        etain_sub[1]->Draw();
//        for (int nep = 0; nep<4; nep++) {
//            etain_sub[nep]->Draw("same");
//        }
//        
//        cMCinputs_sub->cd(3);
//        phiPsiRP_merged->GetXaxis()->SetRangeUser(-4,4);
//        phiPsiRP_merged->Draw("same");
//        
//        TPad * padMCinputs_sub = (TPad *) cMCinputs_sub->cd(4);
//        padMCinputs_sub->SetLogy();
//        ptin_merged->GetXaxis()->SetRangeUser(0,5);
//        ptin_merged->Draw("same");
//        
//        
//        //    TPaveText * txMCinputs_sub_outputs = new TPaveText(0.13, 0.19, 0.83, 0.88,"NDC");
//        //    txMCinputs_sub_outputs->SetFillColor(0);
//        //    txMCinputs_sub_outputs->SetBorderSize(0);
//        //    txMCinputs_sub_outputs->SetTextFont(43);
//        //    txMCinputs_sub_outputs->SetTextSize(20);
//        //    txMCinputs_sub_outputs->SetTextAlign(12);
//        //    txMCinputs_sub_outputs->AddText("MC outputs");
//        //    txMCinputs_sub_outputs->AddText("  #it{Event plane}");
//        //    txMCinputs_sub_outputs->AddText(Form("v_{1}{EP}           = %0.4f #pm %0.4f",v1outEP_3SE[0]->GetMean(),v1outEP_2SE[0]->GetMeanError()));
//        //    txMCinputs_sub_outputs->AddText(Form("v_{1}{mixed EP} = %0.4f #pm %0.4f",v1outEP_mix[0]->GetMean(),v1outEP_mix[0]->GetMeanError()));
//        //    txMCinputs_sub_outputs->AddText(Form("v_{2}{EP}           = %0.4f #pm %0.4f",v2outEP_3SE[0]->GetMean(),v2outEP_3SE[0]->GetMeanError()));
//        //    txMCinputs_sub_outputs->AddText("  #it{Scalar-product}");
//        //    txMCinputs_sub_outputs->AddText(Form("v_{1}{SP}           = %0.4f #pm %0.4f",v1outSP_3SE[0]->GetMean(),v1outSP_2SE[0]->GetMeanError()));
//        //    txMCinputs_sub_outputs->AddText(Form("v_{1}{mixed SP} = %0.4f #pm %0.4f",v1outSP_mix[0]->GetMean(),v1outSP_mix[0]->GetMeanError()));
//        //    txMCinputs_sub_outputs->AddText(Form("v_{2}{SP}           = %0.4f #pm %0.4f",v2outSP_3SE[0]->GetMean(),v2outSP_3SE[0]->GetMeanError()));
//        //    txMCinputs_sub_outputs->Draw();
//        
//        if (print_plot) cMCinputs_sub->Print(Form("plots/%s/MCinputs_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//        if (close_plot) cMCinputs_sub->Close();
//        
//        
//        
//        //-- v1(eta)
//        if (isodd) {
//            // HFm
//            TCanvas * cv1outHFm_odd_eta = new TCanvas("cv1outHFm_odd_eta","cv1outHFm_odd_eta",600,550);
//            TPad * padv1outHFm_odd_eta = (TPad *) cv1outHFm_odd_eta->cd();
//            if (gridlines) padv1outHFm_odd_eta->SetGrid();
//            TH1D * hv1outHFm_odd_eta_tmp = (TH1D *) hdummy_eta->Clone("v1outHFm_odd_eta_tmp");
//            hv1outHFm_odd_eta_tmp->GetYaxis()->SetRangeUser(-0.03,0.03);
//            hv1outHFm_odd_eta_tmp->Draw();
//            v1input->Draw("hist same L");
//            v1EP_3SE_eta[0]->Draw("same");
//            v1SP_3SE_eta[0]->Draw("same");
//            v1EP_Mix_v1_eta[0]->Scale(0.92); // debug later
//            v1EP_Mix_v1_eta[0]->Draw("same");
//            v1SP_Mix_v1_eta[0]->Draw("same");
//            
//            TLegend * legv1outHFm_odd_eta = new TLegend(0.22, 0.20, 0.48, 0.45);
//            legv1outHFm_odd_eta->SetFillColor(0);
//            legv1outHFm_odd_eta->SetBorderSize(0);
//            legv1outHFm_odd_eta->SetTextFont(43);
//            legv1outHFm_odd_eta->SetTextSize(22);
//            legv1outHFm_odd_eta->AddEntry(v1EP_3SE_eta[0],"v_{1}{EP}","p");
//            legv1outHFm_odd_eta->AddEntry(v1SP_3SE_eta[0],"v_{1}{SP}","p");
//            legv1outHFm_odd_eta->AddEntry(v1EP_Mix_v1_eta[0],"v_{1}{mixed EP}","p");
//            legv1outHFm_odd_eta->AddEntry(v1SP_Mix_v1_eta[0],"v_{1}{mixed SP}","p");
//            legv1outHFm_odd_eta->AddEntry(v1input,"Input v_{1}","l");
//            legv1outHFm_odd_eta->Draw();
//            
//            TPaveText * txv1outHFm_odd_eta = new TPaveText(0.62, 0.68, 0.82, 0.91,"NDC");
//            txv1outHFm_odd_eta->SetFillColor(0);
//            txv1outHFm_odd_eta->SetBorderSize(0);
//            txv1outHFm_odd_eta->SetTextFont(43);
//            txv1outHFm_odd_eta->SetTextSize(22);
//            txv1outHFm_odd_eta->SetTextAlign(12);
//            txv1outHFm_odd_eta->AddText("MC inputs");
//            txv1outHFm_odd_eta->AddText(Form("Nevents: %d",Nevents));
//            txv1outHFm_odd_eta->AddText(Form("Mult per event: %d",Mult));
//            if (!isodd) txv1outHFm_odd_eta->AddText(Form("Input v_{1}: %0.3f",v1in));
//            txv1outHFm_odd_eta->AddText(Form("Input v_{2}: %0.3f",v2in));
//            txv1outHFm_odd_eta->Draw();
//            
//            lnv1origin_eta->Draw();
//            
//            if (print_plot) cv1outHFm_odd_eta->Print(Form("plots/%s/v1outHFm_eta_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//            if (close_plot) cv1outHFm_odd_eta->Close();
//            
//            
//            // HFp
//            TCanvas * cv1outHFp_odd_eta = new TCanvas("cv1outHFp_odd_eta","cv1outHFp_odd_eta",600,550);
//            TPad * padv1outHFp_odd_eta = (TPad *) cv1outHFp_odd_eta->cd();
//            if (gridlines) padv1outHFp_odd_eta->SetGrid();
//            TH1D * hv1outHFp_odd_eta_tmp = (TH1D *) hdummy_eta->Clone("v1outHFp_odd_eta_tmp");
//            hv1outHFp_odd_eta_tmp->GetYaxis()->SetRangeUser(-0.03,0.03);
//            hv1outHFp_odd_eta_tmp->Draw();
//            v1input->Draw("hist same L");
//            v1EP_3SE_eta[1]->Draw("same");
//            v1SP_3SE_eta[1]->Draw("same");
//            v1EP_Mix_v1_eta[1]->Scale(0.92); // debug later
//            v1EP_Mix_v1_eta[1]->Draw("same");
//            v1SP_Mix_v1_eta[1]->Draw("same");
//            
//            TLegend * legv1outHFp_odd_eta = new TLegend(0.22, 0.20, 0.48, 0.45);
//            legv1outHFp_odd_eta->SetFillColor(0);
//            legv1outHFp_odd_eta->SetBorderSize(0);
//            legv1outHFp_odd_eta->SetTextFont(43);
//            legv1outHFp_odd_eta->SetTextSize(22);
//            legv1outHFp_odd_eta->AddEntry(v1EP_3SE_eta[1],"v_{1}{EP}","p");
//            legv1outHFp_odd_eta->AddEntry(v1SP_3SE_eta[1],"v_{1}{SP}","p");
//            legv1outHFp_odd_eta->AddEntry(v1EP_Mix_v1_eta[1],"v_{1}{mixed EP}","p");
//            legv1outHFp_odd_eta->AddEntry(v1SP_Mix_v1_eta[1],"v_{1}{mixed SP}","p");
//            legv1outHFp_odd_eta->AddEntry(v1input,"Input v_{1}","l");
//            legv1outHFp_odd_eta->Draw();
//            
//            TPaveText * txv1outHFp_odd_eta = new TPaveText(0.62, 0.68, 0.82, 0.91,"NDC");
//            txv1outHFp_odd_eta->SetFillColor(0);
//            txv1outHFp_odd_eta->SetBorderSize(0);
//            txv1outHFp_odd_eta->SetTextFont(43);
//            txv1outHFp_odd_eta->SetTextSize(22);
//            txv1outHFp_odd_eta->SetTextAlign(12);
//            txv1outHFp_odd_eta->AddText("MC inputs");
//            txv1outHFp_odd_eta->AddText(Form("Nevents: %d",Nevents));
//            txv1outHFp_odd_eta->AddText(Form("Mult per event: %d",Mult));
//            if (!isodd) txv1outHFp_odd_eta->AddText(Form("Input v_{1}: %0.3f",v1in));
//            txv1outHFp_odd_eta->AddText(Form("Input v_{2}: %0.3f",v2in));
//            txv1outHFp_odd_eta->Draw();
//            
//            lnv1origin_eta->Draw();
//            
//            if (print_plot) cv1outHFp_odd_eta->Print(Form("plots/%s/v1outHFp_eta_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//            if (close_plot) cv1outHFp_odd_eta->Close();
//        }
//        
//        if (!isodd) {
//            // v1(eta) HFm
//            TCanvas * cv1outHFm_even_eta = new TCanvas("cv1outHFm_even_eta","cv1outHFm_even_eta",600,550);
//            TPad * padv1outHFm_even_eta = (TPad *) cv1outHFm_even_eta->cd();
//            if (gridlines) padv1outHFm_even_eta->SetGrid();
//            TH1D * hv1outHFm_even_eta_tmp = (TH1D *) hdummy_eta->Clone("v1outHFm_even_eta_tmp");
//            if (iseven) hv1outHFm_even_eta_tmp->GetYaxis()->SetRangeUser(0.0125,0.025);
//            hv1outHFm_even_eta_tmp->Draw();
//            v1input->Draw("hist same L");
//            v1EP_3SE_eta[0]->Draw("same");
//            v1SP_3SE_eta[0]->Draw("same");
//            v1EP_Mix_v1_eta[0]->Scale(0.91); // debug later
//            v1EP_Mix_v1_eta[0]->Draw("same");
//            v1SP_Mix_v1_eta[0]->Draw("same");
//            v1SP_Mix_v1_eta[0]->Scale(0.98); // debug later
//            
//            TLegend * legv1outHFm_even_eta = new TLegend(0.22, 0.66, 0.48, 0.91);
//            legv1outHFm_even_eta->SetFillColor(0);
//            legv1outHFm_even_eta->SetBorderSize(0);
//            legv1outHFm_even_eta->SetTextFont(43);
//            legv1outHFm_even_eta->SetTextSize(22);
//            legv1outHFm_even_eta->AddEntry(v1EP_3SE_eta[0],"v_{1}{EP}","p");
//            legv1outHFm_even_eta->AddEntry(v1SP_3SE_eta[0],"v_{1}{SP}","p");
//            legv1outHFm_even_eta->AddEntry(v1EP_Mix_v1_eta[0],"v_{1}{mixed EP}","p");
//            legv1outHFm_even_eta->AddEntry(v1SP_Mix_v1_eta[0],"v_{1}{mixed SP}","p");
//            legv1outHFm_even_eta->AddEntry(v1input,"Input v_{1}","l");
//            legv1outHFm_even_eta->Draw();
//            
//            TPaveText * txv1outHFm_even_eta = new TPaveText(0.58, 0.67, 0.78, 0.90,"NDC");
//            txv1outHFm_even_eta->SetFillColor(0);
//            txv1outHFm_even_eta->SetBorderSize(0);
//            txv1outHFm_even_eta->SetTextFont(43);
//            txv1outHFm_even_eta->SetTextSize(22);
//            txv1outHFm_even_eta->SetTextAlign(12);
//            txv1outHFm_even_eta->AddText("MC inputs");
//            txv1outHFm_even_eta->AddText(Form("Nevents: %d",Nevents));
//            txv1outHFm_even_eta->AddText(Form("Mult per event: %d",Mult));
//            if (iseven) txv1outHFm_even_eta->AddText(Form("Input v_{1}: %0.3f",v1in));
//            txv1outHFm_even_eta->AddText(Form("Input v_{2}: %0.3f",v2in));
//            txv1outHFm_even_eta->Draw();
//            
//            lnv1origin_eta->Draw();
//            
//            if (print_plot) cv1outHFm_even_eta->Print(Form("plots/%s/v1outHFm_eta_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//            if (close_plot) cv1outHFm_even_eta->Close();
//            
//            
//            // v1(eta) HFp
//            TCanvas * cv1outHFp_even_eta = new TCanvas("cv1outHFp_even_eta","cv1outHFp_even_eta",600,550);
//            TPad * padv1outHFp_even_eta = (TPad *) cv1outHFp_even_eta->cd();
//            if (gridlines) padv1outHFp_even_eta->SetGrid();
//            TH1D * hv1outHFp_even_eta_tmp = (TH1D *) hdummy_eta->Clone("v1outHFp_even_eta_tmp");
//            if (iseven) hv1outHFp_even_eta_tmp->GetYaxis()->SetRangeUser(0.0125,0.025);
//            hv1outHFp_even_eta_tmp->Draw();
//            v1input->Draw("hist same L");
//            v1EP_3SE_eta[1]->Draw("same");
//            v1SP_3SE_eta[1]->Draw("same");
//            v1EP_Mix_v1_eta[1]->Scale(0.92); // debug later
//            v1EP_Mix_v1_eta[1]->Draw("same");
//            v1SP_Mix_v1_eta[1]->Draw("same");
//            
//            TLegend * legv1outHFp_even_eta = new TLegend(0.22, 0.66, 0.48, 0.91);
//            legv1outHFp_even_eta->SetFillColor(0);
//            legv1outHFp_even_eta->SetBorderSize(0);
//            legv1outHFp_even_eta->SetTextFont(43);
//            legv1outHFp_even_eta->SetTextSize(22);
//            legv1outHFp_even_eta->AddEntry(v1EP_3SE_eta[1],"v_{1}{EP}","p");
//            legv1outHFp_even_eta->AddEntry(v1SP_3SE_eta[1],"v_{1}{SP}","p");
//            legv1outHFp_even_eta->AddEntry(v1EP_Mix_v1_eta[1],"v_{1}{mixed EP}","p");
//            legv1outHFp_even_eta->AddEntry(v1SP_Mix_v1_eta[1],"v_{1}{mixed SP}","p");
//            legv1outHFp_even_eta->AddEntry(v1input,"Input v_{1}","l");
//            legv1outHFp_even_eta->Draw();
//            
//            TPaveText * txv1outHFp_even_eta = new TPaveText(0.58, 0.67, 0.78, 0.90,"NDC");
//            txv1outHFp_even_eta->SetFillColor(0);
//            txv1outHFp_even_eta->SetBorderSize(0);
//            txv1outHFp_even_eta->SetTextFont(43);
//            txv1outHFp_even_eta->SetTextSize(22);
//            txv1outHFp_even_eta->SetTextAlign(12);
//            txv1outHFp_even_eta->AddText("MC inputs");
//            txv1outHFp_even_eta->AddText(Form("Nevents: %d",Nevents));
//            txv1outHFp_even_eta->AddText(Form("Mult per event: %d",Mult));
//            if (iseven) txv1outHFp_even_eta->AddText(Form("Input v_{1}: %0.3f",v1in));
//            txv1outHFp_even_eta->AddText(Form("Input v_{2}: %0.3f",v2in));
//            txv1outHFp_even_eta->Draw();
//            
//            lnv1origin_eta->Draw();
//            
//            if (print_plot) cv1outHFp_even_eta->Print(Form("plots/%s/v1outHFp_eta_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//            if (close_plot) cv1outHFp_even_eta->Close();
//            
//            
//            
//            // v1(pt) HFm
//            TCanvas * cv1outHFm_even_pt = new TCanvas("cv1outHFm_even_pt","cv1outHFm_even_pt",600,550);
//            TPad * padv1outHFm_even_pt = (TPad *) cv1outHFm_even_pt->cd();
//            if (gridlines) padv1outHFm_even_pt->SetGrid();
//            TH1D * hv1outHFm_even_pt_tmp = (TH1D *) hdummy_pt->Clone("v1outHFm_even_pt_tmp");
//            if (iseven) hv1outHFm_even_pt_tmp->GetYaxis()->SetRangeUser(0.0110,0.0235);
//            hv1outHFm_even_pt_tmp->Draw();
//            TLine * lnv1in = new TLine(0.0, v1in, 8.0, v1in);
//            lnv1in->SetLineColor(kMagenta);
//            lnv1in->SetLineStyle(2);
//            lnv1in->SetLineWidth(2);
//            lnv1in->Draw();
//            v1EP_3SE_pt[0]->Draw("same");
//            v1SP_3SE_pt[0]->Draw("same");
//            v1EP_Mix_v1_pt[0]->Scale(0.90); // debug later
//            v1EP_Mix_v1_pt[0]->Draw("same");
//            v1SP_Mix_v1_pt[0]->Scale(0.98);
//            v1SP_Mix_v1_pt[0]->Draw("same");
//            
//            TLegend * legv1outHFm_even_pt = new TLegend(0.22, 0.66, 0.48, 0.91);
//            legv1outHFm_even_pt->SetFillColor(0);
//            legv1outHFm_even_pt->SetBorderSize(0);
//            legv1outHFm_even_pt->SetTextFont(43);
//            legv1outHFm_even_pt->SetTextSize(22);
//            legv1outHFm_even_pt->AddEntry(v1EP_3SE_pt[0],"v_{1}{EP}","p");
//            legv1outHFm_even_pt->AddEntry(v1SP_3SE_pt[0],"v_{1}{SP}","p");
//            legv1outHFm_even_pt->AddEntry(v1EP_Mix_v1_pt[0],"v_{1}{mixed EP}","p");
//            legv1outHFm_even_pt->AddEntry(v1SP_Mix_v1_pt[0],"v_{1}{mixed SP}","p");
//            legv1outHFm_even_pt->AddEntry(lnv1in,"Input v_{1}","l");
//            legv1outHFm_even_pt->Draw();
//            
//            TPaveText * txv1outHFm_even_pt = new TPaveText(0.58, 0.67, 0.78, 0.90,"NDC");
//            txv1outHFm_even_pt->SetFillColor(0);
//            txv1outHFm_even_pt->SetBorderSize(0);
//            txv1outHFm_even_pt->SetTextFont(43);
//            txv1outHFm_even_pt->SetTextSize(22);
//            txv1outHFm_even_pt->SetTextAlign(12);
//            txv1outHFm_even_pt->AddText("MC inputs");
//            txv1outHFm_even_pt->AddText(Form("Nevents: %d",Nevents));
//            txv1outHFm_even_pt->AddText(Form("Mult per event: %d",Mult));
//            if (iseven) txv1outHFm_even_pt->AddText(Form("Input v_{1}: %0.3f",v1in));
//            txv1outHFm_even_pt->AddText(Form("Input v_{2}: %0.3f",v2in));
//            txv1outHFm_even_pt->Draw();
//            
//            lnv1origin_pt->Draw();
//            
//            if (print_plot) cv1outHFm_even_pt->Print(Form("plots/%s/v1outHFm_pt_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//            if (close_plot) cv1outHFm_even_pt->Close();
//            
//            
//            
//            // v1(pt) HFp
//            TCanvas * cv1outHFp_even_pt = new TCanvas("cv1outHFp_even_pt","cv1outHFp_even_pt",600,550);
//            TPad * padv1outHFp_even_pt = (TPad *) cv1outHFp_even_pt->cd();
//            if (gridlines) padv1outHFp_even_pt->SetGrid();
//            TH1D * hv1outHFp_even_pt_tmp = (TH1D *) hdummy_pt->Clone("v1outHFp_even_pt_tmp");
//            if (iseven) hv1outHFp_even_pt_tmp->GetYaxis()->SetRangeUser(0.0110,0.0235);
//            hv1outHFp_even_pt_tmp->Draw();
//            lnv1in->Draw();
//            v1EP_3SE_pt[1]->Draw("same");
//            v1SP_3SE_pt[1]->Draw("same");
//            v1EP_Mix_v1_pt[1]->Scale(0.90); // debug later
//            v1EP_Mix_v1_pt[1]->Draw("same");
//            v1SP_Mix_v1_pt[1]->Scale(0.98);
//            v1SP_Mix_v1_pt[1]->Draw("same");
//            
//            TLegend * legv1outHFp_even_pt = new TLegend(0.22, 0.66, 0.48, 0.91);
//            legv1outHFp_even_pt->SetFillColor(0);
//            legv1outHFp_even_pt->SetBorderSize(0);
//            legv1outHFp_even_pt->SetTextFont(43);
//            legv1outHFp_even_pt->SetTextSize(22);
//            legv1outHFp_even_pt->AddEntry(v1EP_3SE_pt[1],"v_{1}{EP}","p");
//            legv1outHFp_even_pt->AddEntry(v1SP_3SE_pt[1],"v_{1}{SP}","p");
//            legv1outHFp_even_pt->AddEntry(v1EP_Mix_v1_pt[1],"v_{1}{mixed EP}","p");
//            legv1outHFp_even_pt->AddEntry(v1SP_Mix_v1_pt[1],"v_{1}{mixed SP}","p");
//            legv1outHFp_even_pt->AddEntry(lnv1in,"Input v_{1}","l");
//            legv1outHFp_even_pt->Draw();
//            
//            TPaveText * txv1outHFp_even_pt = new TPaveText(0.58, 0.67, 0.78, 0.90,"NDC");
//            txv1outHFp_even_pt->SetFillColor(0);
//            txv1outHFp_even_pt->SetBorderSize(0);
//            txv1outHFp_even_pt->SetTextFont(43);
//            txv1outHFp_even_pt->SetTextSize(22);
//            txv1outHFp_even_pt->SetTextAlign(12);
//            txv1outHFp_even_pt->AddText("MC inputs");
//            txv1outHFp_even_pt->AddText(Form("Nevents: %d",Nevents));
//            txv1outHFp_even_pt->AddText(Form("Mult per event: %d",Mult));
//            if (iseven) txv1outHFp_even_pt->AddText(Form("Input v_{1}: %0.3f",v1in));
//            txv1outHFp_even_pt->AddText(Form("Input v_{2}: %0.3f",v2in));
//            txv1outHFp_even_pt->Draw();
//            
//            lnv1origin_pt->Draw();
//            
//            if (print_plot) cv1outHFp_even_pt->Print(Form("plots/%s/v1outHFp_pt_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//            if (close_plot) cv1outHFp_even_pt->Close();
//        }
//        
//        
//        TCanvas * cv1input = new TCanvas("cv1input","cv1input",600,550);
//        TPad * padv1input = (TPad *) cv1input->cd();
//        if (gridlines) padv1input->SetGrid();
//        TH1D * hv1input_tmp = (TH1D *) hdummy_eta->Clone("v1input_tmp");
//        hv1input_tmp->GetXaxis()->SetRangeUser(-6.,6.);
//        hv1input_tmp->GetYaxis()->SetRangeUser(-0.05,0.05);
//        hv1input_tmp->SetXTitle("#eta");
//        hv1input_tmp->SetYTitle("MC input v_{1}");
//        hv1input_tmp->Draw();
//        
//        TBox * bv1input = new TBox(-2.4, -0.05, 2.4, 0.05);
//        bv1input->SetFillColor(kBlue-7);
//        bv1input->SetFillStyle(3001);
//        bv1input->Draw("same");
//        
//        TH1D * hv1in_tmp = (TH1D *) v1input->Clone("hv1in_tmp");
//        hv1in_tmp->SetMarkerColor(kBlack);
//        hv1in_tmp->SetLineColor(kBlack);
//        hv1in_tmp->SetLineStyle(1);
//        hv1in_tmp->SetMarkerStyle(20);
//        hv1in_tmp->SetMarkerSize(1.2);
//        hv1in_tmp->Draw("same");
//        
//        if (print_plot) cv1input->Print(Form("plots/%s/v1outHFm_eta_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//        if (close_plot) cv1input->Close();
//        
//        
//        
//        
//        TCanvas * cetasubs = new TCanvas("cetasubs","cetasubs",600,550);
//        TH1D * hetasubs_tmp = (TH1D *) etain_sub[1]->Clone("etasubs_tmp");
//        hetasubs_tmp->SetXTitle("#eta");
//        hetasubs_tmp->SetYTitle("");
//        hetasubs_tmp->Draw();
//        hetasubs_tmp->SetTitle("");
//        hetasubs_tmp->SetStats(kFALSE);
//        hetasubs_tmp->GetXaxis()->CenterTitle(kTRUE);
//        hetasubs_tmp->GetYaxis()->CenterTitle(kTRUE);
//        hetasubs_tmp->GetXaxis()->SetRangeUser(-6, 6);
//        hetasubs_tmp->GetXaxis()->SetTitleSize(0.05);
//        hetasubs_tmp->GetYaxis()->SetTitleSize(0.05);
//        hetasubs_tmp->GetYaxis()->SetTitleOffset(1.25);
//        hetasubs_tmp->Draw();
//        
//        TBox * betasubHFm = new TBox(-5.0, 0, -3.0, 1.02e7);
//        betasubHFm->SetFillColor(kRed-9);
//        betasubHFm->SetFillStyle(3001);
//        betasubHFm->Draw("same");
//        TBox * betasubtrackm = new TBox(-2.4, 0, 0, 1.02e7);
//        betasubtrackm->SetFillColor(kBlue-9);
//        betasubtrackm->SetFillStyle(3001);
//        betasubtrackm->Draw("same");
//        TBox * betasubtrackp = new TBox(0, 0, 2.4, 1.02e7);
//        betasubtrackp->SetFillColor(kGreen-9);
//        betasubtrackp->SetFillStyle(3001);
//        betasubtrackp->Draw("same");
//        TBox * betasubHFp = new TBox(3.0, 0, 5.0, 1.02e7);
//        betasubHFp->SetFillColor(kMagenta-9);
//        betasubHFp->SetFillStyle(3001);
//        betasubHFp->Draw("same");
//        
//        etain_sub[1]->Draw("same");
//        for (int nep = 0; nep<4; nep++) {
//            etain_sub[nep]->Draw("same");
//        }
//        
//        if (print_plot) cetasubs->Print(Form("plots/%s/etasubs_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//        if (close_plot) cetasubs->Close();
//        
//        
//        
//        TCanvas * cphilab = new TCanvas("cphilab","cphilab",600,550);
//        cphilab->cd();
//        philab_sub[1]->Draw();
//        
//        if (print_plot) cphilab->Print(Form("plots/%s/philab_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//        if (close_plot) cphilab->Close();
//        
//        
//        
//        TCanvas * cpt2D = new TCanvas("cpt2D","cpt2D",800,700);
//        cpt2D->Divide(2,2);
//        TLine * ln2Dx = new TLine(-200, 0, 200, 0);
//        ln2Dx->SetLineWidth(2);
//        ln2Dx->SetLineStyle(2);
//        TLine * ln2Dy = new TLine(0, -200, 0, 200);
//        ln2Dy->SetLineWidth(2);
//        ln2Dy->SetLineStyle(2);
//        TPaveText * tx2D0 = new TPaveText(0.27, 0.81, 0.37, 0.92,"NDC");
//        tx2D0->SetFillColor(0);
//        tx2D0->SetBorderSize(0);
//        tx2D0->SetTextFont(43);
//        tx2D0->SetTextSize(24);
//        tx2D0->SetTextAlign(12);
//        TPaveText * tx2D1 = (TPaveText *) tx2D0->Clone();
//        TPaveText * tx2D2 = (TPaveText *) tx2D0->Clone();
//        TPaveText * tx2D3 = (TPaveText *) tx2D0->Clone();
//        
//        cpt2D->cd(1);
//        pt2D_sub[0]->Draw();
//        ln2Dx->Draw();
//        ln2Dy->Draw();
//        tx2D0->SetTextColor(kRed);
//        tx2D0->AddText("HF-");
//        tx2D0->Draw();
//        cpt2D->cd(2);
//        pt2D_sub[1]->Draw();
//        ln2Dx->Draw();
//        ln2Dy->Draw();
//        tx2D1->SetTextColor(kBlue);
//        tx2D1->AddText("tracker-");
//        tx2D1->Draw();
//        cpt2D->cd(3);
//        pt2D_sub[2]->Draw();
//        ln2Dx->Draw();
//        ln2Dy->Draw();
//        tx2D2->SetTextColor(kGreen+2);
//        tx2D2->AddText("tracker+");
//        tx2D2->Draw();
//        cpt2D->cd(4);
//        pt2D_sub[3]->Draw();
//        ln2Dx->Draw();
//        ln2Dy->Draw();
//        tx2D3->SetTextColor(kMagenta);
//        tx2D3->AddText("HF+");
//        tx2D3->Draw();
//        
//        if (print_plot) cpt2D->Print(Form("plots/%s/pt2D_%s.pdf",tag.Data(),mtag.Data()),"pdf");
//        if (close_plot) cpt2D->Close();
        
    }
}

