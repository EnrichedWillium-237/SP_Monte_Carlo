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

# include "style.h"

using namespace std;

Bool_t close_plots = kFALSE;
Bool_t gridlines = kFALSE;

static const int nptbins = 21;
static const double ptbins[] = {
     0.00,  0.20,  0.40,  0.60,  0.80,  1.00,  1.20,  1.40,
     1.60,  1.80,  2.00,  2.25,  2.50,  2.75,  3.00,  3.50,
     4.00,  4.50,  5.00,  6.00,  7.00,  8.00};

static const int netabins = 12;
static const double etabins[] = {
    -2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,
     0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const Int_t numEP = 6;
//--                             HF-  track- track+  HF+ trackmid
static const double ecutmin[] = {-5.0, -2.4,  2.0,  3.0, -0.5};
static const double ecutmax[] = {-3.0, -2.0,  2.4,  5.0,  0.5};
static const double pcutmin[] = { 0.3,  0.3,  0.3,  0.3,  0.3};
static const double pcutmax[] = {30.0,  3.0,  3.0, 30.0, 3.9};

static const int nv1Ebins = 28;
static const double v1Ebins[] = {
    -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0,
    -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,
     2.4,  2.8,  3.2,  3.6,  4.0,  4.4,  4.8,  5.2,  5.6};

const std::string SideName[] = { "min", "pos" };
const std::string EPName[] = { "HFm", "trackm", "trackp", "HFp", "trkmid" };

//-- initial MC throw

TH1D * inputParms;
TH1D * init_v1in;

//-- event plane and scalar-product vn

TH1D * rescor1_2SE;
TH1D * rescor1_3SE[2];
TH1D * rescor1_mix[2];
TH1D * rescor1_mix_3[2];
TH1D * rescor2HF_2SE;
TH1D * rescor2HF_3SE[2];
TH1D * rescor2Trk_3SE[2];
TH1D * rescor3HF_2SE;
TH1D * rescor3HF_3SE[2];
TH1D * rescor3Trk_3SE[2];
TH1D * EPv1_2SE[2];
TH1D * EPv1_3SE[2];
TH1D * EPv1_mix[2];
TH1D * EPv1_mix_3[2];
TH1D * EPv2_2SE[2];
TH1D * EPv2_3SE[2];
TH1D * EPv3_2SE[2];
TH1D * EPv3_3SE[2];

TH1D * SPdenom1_2SE;
TH1D * SPdenom1_3SE[2];
TH1D * SPdenom1_mix[2];
TH1D * SPdenom1_mix_3[2];
TH1D * SPdenom2HF_2SE;
TH1D * SPdenom2_3SE[2];
TH1D * SPdenom3HF_2SE;
TH1D * SPdenom3HF_3SE[2];
TH1D * SPv1_2SE[2];
TH1D * SPv1_3SE[2];
TH1D * SPv1_mix[2];
TH1D * SPv1_mix_3[2];
TH1D * SPv2_2SE[2];
TH1D * SPv2_3SE[2];
TH1D * SPv3_2SE[2];
TH1D * SPv3_3SE[2];

//-- differential vn(pT)

TH1D * DiffEPv1_2SE_pt[2];
TH1D * DiffEPv1_3SE_pt[2];
TH1D * DiffEPv1_mix_pt[2];
TH1D * DiffEPv1_mix_3_pt[2];
TH1D * DiffEPv2_2SE_pt[2];
TH1D * DiffEPv2_3SE_pt[2];
TH1D * DiffEPv3_2SE_pt[2];
TH1D * DiffEPv3_3SE_pt[2];

TH1D * DiffSPv1_2SE_pt[2];
TH1D * DiffSPv1_3SE_pt[2];
TH1D * DiffSPv1_mix_pt[2];
TH1D * DiffSPv1_mix_3_pt[2];
TH1D * DiffSPv2_2SE_pt[2];
TH1D * DiffSPv2_3SE_pt[2];
TH1D * DiffSPv3_2SE_pt[2];
TH1D * DiffSPv3_3SE_pt[2];

TH1D * Diffqcnt_pt[2];
TH1D * Diffncnt_pt[2];

//-- differential v1(eta)

TH1D * DiffEPv1_2SE_eta[2];
TH1D * DiffEPv1_3SE_eta[2];
TH1D * DiffEPv1_mix_eta[2];
TH1D * DiffEPv1_mix_3_eta[2];

TH1D * DiffSPv1_2SE_eta[2];
TH1D * DiffSPv1_3SE_eta[2];
TH1D * DiffSPv1_mix_eta[2];
TH1D * DiffSPv1_mix_3_eta[2];

TH1D * Diffqcnt_eta[2];
TH1D * Diffncnt_eta[2];

void v1_eta()
{
    gStyle->SetPalette(55);
    gStyle->SetErrorX(0.5);

    int subcolor[] = {kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+7};

    Int_t NumEvnts = 0;
    Double_t v1in = 0.;
    Double_t v2in = 0.;
    Double_t v3in = 0.;
    Bool_t isodd = kFALSE;
    Bool_t eta_weights = kFALSE;
    Bool_t pt_weights = kFALSE;
    Bool_t conserve_pT = kFALSE;
    Bool_t recenter = kFALSE;
    Bool_t flatten = kFALSE;
    Double_t etaTrkMin = 0.;
    Double_t etaTrkMax = 0.;


    //results_v1_even_0.0000_v2_0.0700_eta_weights_mom-cons_2.0_to_2.4_1000000_evts.root
    //results_v1_even_0.0000_v2_0.0700_eta_weights_pt_weights_mom-cons_2.0_to_2.4_1000000_evts.root
    //results_v1_even_0.0200_v2_0.0700_eta_weights_2.0_to_2.4_1000000_evts.root
    results_v1_odd_v2_0.0700_eta_weights_2.0_to_2.4_1000000_evts.root


    //-- retrieve data from MC sample

    init_v1in = (TH1D *) tfin->Get("Inputs/v1in");
    inputParms = (TH1D *) tfin->Get("Inputs/Input_Parameters");
    NumEvnts = inputParms->GetBinContent(1);
    if ( inputParms->GetBinContent(2)>0 ) isodd = kTRUE;
    if (!isodd) v1in = init_v1in->GetBinContent(1);
    v2in = inputParms->GetBinContent(3);
    v3in = inputParms->GetBinContent(4);
    if ( inputParms->GetBinContent(5)>0 ) eta_weights = kTRUE;
    if ( inputParms->GetBinContent(6)>0 ) pt_weights = kTRUE;
    if ( inputParms->GetBinContent(7)>0 ) conserve_pT = kTRUE;
    if ( inputParms->GetBinContent(8)>0 ) recenter = kTRUE;
    if ( inputParms->GetBinContent(9)>0 ) flatten = kTRUE;
    etaTrkMin = inputParms->GetBinContent(10);
    etaTrkMax = inputParms->GetBinContent(11);

    rescor1_2SE = (TH1D *) tfin->Get("v1EP/HFp/rescor1_2SE");
    rescor2HF_2SE = (TH1D *) tfin->Get("v1EP/HFp/rescor2HF_2SE");
    rescor3HF_2SE = (TH1D *) tfin->Get("v1EP/HFp/rescor3HF_2SE");
    for (int iside = 0; iside<2; iside++) {
        rescor1_3SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/rescor1_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        rescor1_mix[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/rescor1_mix_%s",SideName[iside].data(),SideName[iside].data()));
        rescor1_mix_3[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/rescor1_mix_3_%s",SideName[iside].data(),SideName[iside].data()));
        rescor2HF_3SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/rescor2HF_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        rescor3HF_3SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/rescor3HF_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        EPv1_2SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/EPv1_2SE_%s",SideName[iside].data(),SideName[iside].data()));
        EPv1_3SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/EPv1_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        EPv1_mix[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/EPv1_mix_%s",SideName[iside].data(),SideName[iside].data()));
        EPv1_mix_3[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/EPv1_mix_3_%s",SideName[iside].data(),SideName[iside].data()));
        EPv2_2SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/EPv2_2SE_%s",SideName[iside].data(),SideName[iside].data()));
        EPv2_3SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/EPv2_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        EPv3_2SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/EPv3_2SE_%s",SideName[iside].data(),SideName[iside].data()));
        EPv3_3SE[iside] = (TH1D *) tfin->Get(Form("v1EP/%s/EPv3_3SE_%s",SideName[iside].data(),SideName[iside].data()));
    }

    SPdenom1_2SE = (TH1D *) tfin->Get("v1SP/HFp/SPdenom1_2SE");
    SPdenom2HF_2SE = (TH1D *) tfin->Get("v1SP/HFp/SPdenom2HF_2SE");
    SPdenom3HF_2SE = (TH1D *) tfin->Get("v1SP/HFp/SPdenom3HF_2SE");
    for (int iside = 0; iside<2; iside++) {
        SPdenom1_3SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPdenom1_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        SPdenom1_mix[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPdenom1_mix_%s",SideName[iside].data(),SideName[iside].data()));
        SPdenom1_mix_3[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPdenom1_mix_3_%s",SideName[iside].data(),SideName[iside].data()));
        SPdenom2_3SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPdenom2_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        SPdenom3HF_3SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPdenom3HF_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        SPv1_2SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPv1_2SE_%s",SideName[iside].data(),SideName[iside].data()));
        SPv1_3SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPv1_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        SPv1_mix[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPv1_mix_%s",SideName[iside].data(),SideName[iside].data()));
        SPv1_mix_3[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPv1_mix_3_%s",SideName[iside].data(),SideName[iside].data()));
        SPv2_2SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPv2_2SE_%s",SideName[iside].data(),SideName[iside].data()));
        SPv2_3SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPv2_3SE_%s",SideName[iside].data(),SideName[iside].data()));
        SPv3_2SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPv3_2SE_%s",SideName[iside].data(),SideName[iside].data()));
        SPv3_3SE[iside] = (TH1D *) tfin->Get(Form("v1SP/%s/SPv3_3SE_%s",SideName[iside].data(),SideName[iside].data()));
    }

    for (int iside = 0; iside<2; iside++) {
        DiffEPv1_2SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/EPv1_2SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv1_3SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/EPv1_3SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv1_mix_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/EPv1_mix_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv1_mix_3_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/EPv1_mix_3_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv2_2SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/EPv2_2SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv2_3SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/EPv2_3SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv3_2SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/EPv3_2SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv3_3SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/EPv3_3SE_pt_%s",SideName[iside].data(),SideName[iside].data()));

        DiffSPv1_2SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/SPv1_2SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv1_3SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/SPv1_3SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv1_mix_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/SPv1_mix_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv1_mix_3_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/SPv1_mix_3_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv2_2SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/SPv2_2SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv2_3SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/SPv2_3SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv3_2SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/SPv3_2SE_pt_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv3_3SE_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/SPv3_3SE_pt_%s",SideName[iside].data(),SideName[iside].data()));

        Diffqcnt_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/Qcnt_pt_%s",SideName[iside].data(),SideName[iside].data()));
        Diffncnt_pt[iside] = (TH1D *) tfin->Get(Form("vn_pt/%s/Ncnt_pt_%s",SideName[iside].data(),SideName[iside].data()));

        DiffEPv1_2SE_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/EPv1_2SE_eta_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv1_3SE_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/EPv1_3SE_eta_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv1_mix_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/EPv1_mix_eta_%s",SideName[iside].data(),SideName[iside].data()));
        DiffEPv1_mix_3_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/EPv1_mix_3_eta_%s",SideName[iside].data(),SideName[iside].data()));

        DiffSPv1_2SE_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/SPv1_2SE_eta_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv1_3SE_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/SPv1_3SE_eta_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv1_mix_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/SPv1_mix_eta_%s",SideName[iside].data(),SideName[iside].data()));
        DiffSPv1_mix_3_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/SPv1_mix_3_eta_%s",SideName[iside].data(),SideName[iside].data()));

        Diffqcnt_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/Qcnt_eta__%s",SideName[iside].data(),SideName[iside].data()));
        Diffncnt_eta[iside] = (TH1D *) tfin->Get(Form("vn_eta/%s/Ncnt_eta_%s",SideName[iside].data(),SideName[iside].data()));
    }

    TString mtag = "v1";
    if (isodd) mtag+="_odd";
    else mtag+=Form("_even_%0.4f",v1in);
    mtag+=Form("_v2_%0.4f_v3_%0.4f",v2in,v3in);
    if (eta_weights) mtag+="_eta_weights";
    if (pt_weights) mtag+="_pt_weights";
    if (conserve_pT) mtag+="_momcons";

    init_v1in->SetMarkerColor(kMagenta);
    init_v1in->SetLineColor(kMagenta);
    init_v1in->SetMarkerStyle(25);
    init_v1in->SetMarkerSize(1.2);
    init_v1in->SetLineWidth(2);
    init_v1in->SetLineStyle(2);


    //-- plot v1(eta)

    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen(Form("plots/%s",mtag.Data()),"r")) system(Form("mkdir plots/%s",mtag.Data()));


    //-- Event plane v1(eta) using HFm
    TCanvas * cv1EP_eta_HFm = new TCanvas("cv1EP_eta_HFm","cv1EP_eta_HFm",650,600);
    TPad * padv1EP_eta_HFm = (TPad *) cv1EP_eta_HFm->cd();
    if (gridlines) padv1EP_eta_HFm->SetGrid();
    DiffEPv1_2SE_eta[0]->SetTitle("");
    DiffEPv1_2SE_eta[0]->SetStats(0);
    DiffEPv1_2SE_eta[0]->SetXTitle("#eta");
    DiffEPv1_2SE_eta[0]->SetYTitle("v_{1}");
    DiffEPv1_2SE_eta[0]->SetMarkerColor(kBlue);
    DiffEPv1_2SE_eta[0]->SetLineColor(kBlue);
    DiffEPv1_2SE_eta[0]->SetMarkerStyle(20);
    DiffEPv1_2SE_eta[0]->SetMarkerSize(1.2);
    if (isodd) DiffEPv1_2SE_eta[0]->GetYaxis()->SetRangeUser(-0.015,0.015);
    else DiffEPv1_2SE_eta[0]->GetYaxis()->SetRangeUser(-0.025,0.025);
    DiffEPv1_2SE_eta[0]->Draw("same");

    init_v1in->Draw("hist same c");

    DiffEPv1_3SE_eta[0]->SetMarkerColor(kRed);
    DiffEPv1_3SE_eta[0]->SetLineColor(kRed);
    DiffEPv1_3SE_eta[0]->SetMarkerStyle(21);
    DiffEPv1_3SE_eta[0]->SetMarkerSize(1.1);
    DiffEPv1_3SE_eta[0]->Draw("same");

    DiffEPv1_mix_eta[0]->SetMarkerColor(kGreen+2);
    DiffEPv1_mix_eta[0]->SetLineColor(kGreen+2);
    DiffEPv1_mix_eta[0]->SetMarkerStyle(24);
    DiffEPv1_mix_eta[0]->SetMarkerSize(1.2);
    DiffEPv1_mix_eta[0]->Draw("same");

    TLegend * legv1EP_eta_HFm = new TLegend(0.2, 0.2, 0.4, 0.4);
    SetLegend(legv1EP_eta_HFm, 18);
    legv1EP_eta_HFm->SetHeader("HF- Event plane");
    legv1EP_eta_HFm->AddEntry(init_v1in,"Input v_{1}","l");
    legv1EP_eta_HFm->AddEntry(DiffEPv1_2SE_eta[0],"2 subevent v_{1}","p");
    legv1EP_eta_HFm->AddEntry(DiffEPv1_3SE_eta[0],"3 subevent v_{1}","p");
    legv1EP_eta_HFm->AddEntry(DiffEPv1_mix_eta[0],"mixed","p");
    legv1EP_eta_HFm->Draw();

    TPaveText * txv1EP_eta_HFm = new TPaveText(0.45, 0.17, 0.65, 0.41,"NDC");
    SetTPaveTxt(txv1EP_eta_HFm, 18);
    txv1EP_eta_HFm->AddText(Form("N_{events}: %d",NumEvnts));
    if (!isodd) txv1EP_eta_HFm->AddText(Form("Input v_{1}: %0.3f",v1in));
    txv1EP_eta_HFm->AddText(Form("Input v_{2}: %0.3f",v2in));
    txv1EP_eta_HFm->AddText(Form("Input v_{3}: %0.3f",v3in));
    if (eta_weights) txv1EP_eta_HFm->AddText("#eta-weights");
    if (pt_weights) txv1EP_eta_HFm->AddText("p_{T}-weights");
    if (conserve_pT) txv1EP_eta_HFm->AddText("p_{T} conserved");
    txv1EP_eta_HFm->Draw();

    cv1EP_eta_HFm->Print(Form("plots/%s/v1EP_eta_HFm.png",mtag.Data()),"png");
    if (close_plots) cv1EP_eta_HFm->Close();


    //-- Event plane v1(eta) using HFp
    TCanvas * cv1EP_eta_HFp = new TCanvas("cv1EP_eta_HFp","cv1EP_eta_HFp",650,600);
    TPad * padv1EP_eta_HFp = (TPad *) cv1EP_eta_HFp->cd();
    if (gridlines) padv1EP_eta_HFp->SetGrid();
    DiffEPv1_2SE_eta[1]->SetTitle("");
    DiffEPv1_2SE_eta[1]->SetStats(0);
    DiffEPv1_2SE_eta[1]->SetXTitle("#eta");
    DiffEPv1_2SE_eta[1]->SetYTitle("v_{1}");
    DiffEPv1_2SE_eta[1]->SetMarkerColor(kBlue);
    DiffEPv1_2SE_eta[1]->SetLineColor(kBlue);
    DiffEPv1_2SE_eta[1]->SetMarkerStyle(20);
    DiffEPv1_2SE_eta[1]->SetMarkerSize(1.2);
    if (isodd) DiffEPv1_2SE_eta[1]->GetYaxis()->SetRangeUser(-0.015,0.015);
    else DiffEPv1_2SE_eta[1]->GetYaxis()->SetRangeUser(v1in-0.01,v1in+0.01);
    DiffEPv1_2SE_eta[1]->Draw("same");

    init_v1in->Draw("hist same c");

    DiffEPv1_3SE_eta[1]->SetMarkerColor(kRed);
    DiffEPv1_3SE_eta[1]->SetLineColor(kRed);
    DiffEPv1_3SE_eta[1]->SetMarkerStyle(21);
    DiffEPv1_3SE_eta[1]->SetMarkerSize(1.1);
    DiffEPv1_3SE_eta[1]->Draw("same");

    DiffEPv1_mix_eta[1]->SetMarkerColor(kGreen+2);
    DiffEPv1_mix_eta[1]->SetLineColor(kGreen+2);
    DiffEPv1_mix_eta[1]->SetMarkerStyle(24);
    DiffEPv1_mix_eta[1]->SetMarkerSize(1.2);
    DiffEPv1_mix_eta[1]->Draw("same");

    TLegend * legv1EP_eta_HFp = new TLegend(0.2, 0.2, 0.4, 0.4);
    SetLegend(legv1EP_eta_HFp, 18);
    legv1EP_eta_HFp->SetHeader("HF+ Event plane");
    legv1EP_eta_HFp->AddEntry(init_v1in,"Input v_{1}","l");
    legv1EP_eta_HFp->AddEntry(DiffEPv1_2SE_eta[1],"2 subevent v_{1}","p");
    legv1EP_eta_HFp->AddEntry(DiffEPv1_3SE_eta[1],"3 subevent v_{1}","p");
    legv1EP_eta_HFp->AddEntry(DiffEPv1_mix_eta[1],"mixed","p");
    legv1EP_eta_HFp->Draw();

    TPaveText * txv1EP_eta_HFp = new TPaveText(0.45, 0.17, 0.65, 0.41,"NDC");
    SetTPaveTxt(txv1EP_eta_HFp, 18);
    txv1EP_eta_HFp->AddText(Form("N_{events}: %d",NumEvnts));
    if (!isodd) txv1EP_eta_HFp->AddText(Form("Input v_{1}: %0.3f",v1in));
    txv1EP_eta_HFp->AddText(Form("Input v_{2}: %0.3f",v2in));
    txv1EP_eta_HFp->AddText(Form("Input v_{3}: %0.3f",v3in));
    if (eta_weights) txv1EP_eta_HFp->AddText("#eta-weights");
    if (pt_weights) txv1EP_eta_HFp->AddText("p_{T}-weights");
    if (conserve_pT) txv1EP_eta_HFp->AddText("p_{T} conserved");
    txv1EP_eta_HFp->Draw();

    cv1EP_eta_HFp->Print(Form("plots/%s/v1EP_eta_HFp.png",mtag.Data()),"png");
    if (close_plots) cv1EP_eta_HFp->Close();



    //-- Scalar-product v1(eta) using HFm
    TCanvas * cv1SP_eta_HFm = new TCanvas("cv1SP_eta_HFm","cv1SP_eta_HFm",650,600);
    TPad * padv1SP_eta_HFm = (TPad *) cv1SP_eta_HFm->cd();
    if (gridlines) padv1SP_eta_HFm->SetGrid();
    DiffSPv1_2SE_eta[0]->SetTitle("");
    DiffSPv1_2SE_eta[0]->SetStats(0);
    DiffSPv1_2SE_eta[0]->SetXTitle("#eta");
    DiffSPv1_2SE_eta[0]->SetYTitle("v_{1}");
    DiffSPv1_2SE_eta[0]->SetMarkerColor(kBlue);
    DiffSPv1_2SE_eta[0]->SetLineColor(kBlue);
    DiffSPv1_2SE_eta[0]->SetMarkerStyle(20);
    DiffSPv1_2SE_eta[0]->SetMarkerSize(1.2);
    if (isodd) DiffSPv1_2SE_eta[0]->GetYaxis()->SetRangeUser(-0.015,0.015);
    else DiffSPv1_2SE_eta[0]->GetYaxis()->SetRangeUser(-0.025,0.025);
    DiffSPv1_2SE_eta[0]->Draw("same");


    init_v1in->Draw("hist same c");

    DiffSPv1_3SE_eta[0]->SetMarkerColor(kRed);
    DiffSPv1_3SE_eta[0]->SetLineColor(kRed);
    DiffSPv1_3SE_eta[0]->SetMarkerStyle(21);
    DiffSPv1_3SE_eta[0]->SetMarkerSize(1.1);
    DiffSPv1_3SE_eta[0]->Draw("same");

    DiffSPv1_mix_eta[0]->SetMarkerColor(kGreen+2);
    DiffSPv1_mix_eta[0]->SetLineColor(kGreen+2);
    DiffSPv1_mix_eta[0]->SetMarkerStyle(24);
    DiffSPv1_mix_eta[0]->SetMarkerSize(1.2);
    DiffSPv1_mix_eta[0]->Draw("same");

    TLegend * legv1SP_eta_HFm = new TLegend(0.2, 0.2, 0.4, 0.4);
    SetLegend(legv1SP_eta_HFm, 18);
    legv1SP_eta_HFm->SetHeader("HF- Scalar-product");
    legv1SP_eta_HFm->AddEntry(init_v1in,"Input v_{1}","l");
    legv1SP_eta_HFm->AddEntry(DiffSPv1_2SE_eta[0],"2 subevent v_{1}","p");
    legv1SP_eta_HFm->AddEntry(DiffSPv1_3SE_eta[0],"3 subevent v_{1}","p");
    legv1SP_eta_HFm->AddEntry(DiffSPv1_mix_eta[0],"mixed","p");
    legv1SP_eta_HFm->Draw();

    TPaveText * txv1SP_eta_HFm = new TPaveText(0.45, 0.17, 0.65, 0.41,"NDC");
    SetTPaveTxt(txv1SP_eta_HFm, 18);
    txv1SP_eta_HFm->AddText(Form("N_{events}: %d",NumEvnts));
    if (!isodd) txv1SP_eta_HFm->AddText(Form("Input v_{1}: %0.3f",v1in));
    txv1SP_eta_HFm->AddText(Form("Input v_{2}: %0.3f",v2in));
    txv1SP_eta_HFm->AddText(Form("Input v_{3}: %0.3f",v3in));
    if (eta_weights) txv1SP_eta_HFm->AddText("#eta-weights");
    if (pt_weights) txv1SP_eta_HFm->AddText("p_{T}-weights");
    if (conserve_pT) txv1SP_eta_HFm->AddText("p_{T} conserved");
    txv1SP_eta_HFm->Draw();

    cv1SP_eta_HFm->Print(Form("plots/%s/v1SP_eta_HFm.png",mtag.Data()),"png");
    if (close_plots) cv1SP_eta_HFm->Close();


    //-- Scalar-product v1(eta) using HFp
    TCanvas * cv1SP_eta_HFp = new TCanvas("cv1SP_eta_HFp","cv1SP_eta_HFp",650,600);
    TPad * padv1SP_eta_HFp = (TPad *) cv1SP_eta_HFp->cd();
    if (gridlines) padv1SP_eta_HFp->SetGrid();
    DiffSPv1_2SE_eta[1]->SetTitle("");
    DiffSPv1_2SE_eta[1]->SetStats(0);
    DiffSPv1_2SE_eta[1]->SetXTitle("#eta");
    DiffSPv1_2SE_eta[1]->SetYTitle("v_{1}");
    DiffSPv1_2SE_eta[1]->SetMarkerColor(kBlue);
    DiffSPv1_2SE_eta[1]->SetLineColor(kBlue);
    DiffSPv1_2SE_eta[1]->SetMarkerStyle(20);
    DiffSPv1_2SE_eta[1]->SetMarkerSize(1.2);
    if (isodd) DiffSPv1_2SE_eta[1]->GetYaxis()->SetRangeUser(-0.015,0.015);
    else DiffSPv1_2SE_eta[1]->GetYaxis()->SetRangeUser(v1in-0.01,v1in+0.01);
    DiffSPv1_2SE_eta[1]->Draw("same");

    init_v1in->Draw("hist same c");

    DiffSPv1_3SE_eta[1]->SetMarkerColor(kRed);
    DiffSPv1_3SE_eta[1]->SetLineColor(kRed);
    DiffSPv1_3SE_eta[1]->SetMarkerStyle(21);
    DiffSPv1_3SE_eta[1]->SetMarkerSize(1.1);
    DiffSPv1_3SE_eta[1]->Draw("same");

    DiffSPv1_mix_eta[1]->SetMarkerColor(kGreen+2);
    DiffSPv1_mix_eta[1]->SetLineColor(kGreen+2);
    DiffSPv1_mix_eta[1]->SetMarkerStyle(24);
    DiffSPv1_mix_eta[1]->SetMarkerSize(1.2);
    DiffSPv1_mix_eta[1]->Draw("same");

    TLegend * legv1SP_eta_HFp = new TLegend(0.2, 0.2, 0.4, 0.4);
    SetLegend(legv1SP_eta_HFp, 18);
    legv1SP_eta_HFp->SetHeader("HF+ Scalar-product");
    legv1SP_eta_HFp->AddEntry(init_v1in,"Input v_{1}","l");
    legv1SP_eta_HFp->AddEntry(DiffSPv1_2SE_eta[1],"2 subevent v_{1}","p");
    legv1SP_eta_HFp->AddEntry(DiffSPv1_3SE_eta[1],"3 subevent v_{1}","p");
    legv1SP_eta_HFp->AddEntry(DiffSPv1_mix_eta[1],"mixed","p");
    legv1SP_eta_HFp->Draw();

    TPaveText * txv1SP_eta_HFp = new TPaveText(0.45, 0.17, 0.65, 0.41,"NDC");
    SetTPaveTxt(txv1SP_eta_HFp, 18);
    txv1SP_eta_HFp->AddText(Form("N_{events}: %d",NumEvnts));
    if (!isodd) txv1SP_eta_HFp->AddText(Form("Input v_{1}: %0.3f",v1in));
    txv1SP_eta_HFp->AddText(Form("Input v_{2}: %0.3f",v2in));
    txv1SP_eta_HFp->AddText(Form("Input v_{3}: %0.3f",v3in));
    if (eta_weights) txv1SP_eta_HFp->AddText("#eta-weights");
    if (pt_weights) txv1SP_eta_HFp->AddText("p_{T}-weights");
    if (conserve_pT) txv1SP_eta_HFp->AddText("p_{T} conserved");
    txv1SP_eta_HFp->Draw();

    cv1SP_eta_HFp->Print(Form("plots/%s/v1SP_eta_HFp.png",mtag.Data()),"png");
    if (close_plots) cv1SP_eta_HFp->Close();



    //-- Event plane mixed v1(eta) for both HFp and HFm
    TCanvas * cv1EPmix_eta_HFpm = new TCanvas("cv1EPmix_eta_HFpm","cv1EPmix_eta_HFpm",1200,600);
    cv1EPmix_eta_HFpm->Divide(2,1);
    TPad * padv1EPmix_eta_HFpm = (TPad *) cv1EPmix_eta_HFpm->cd(1);
    if (gridlines) padv1EPmix_eta_HFpm->SetGrid();
    DiffEPv1_mix_eta[0]->SetTitle("");
    DiffEPv1_mix_eta[0]->SetStats(0);
    DiffEPv1_mix_eta[0]->SetXTitle("#eta");
    DiffEPv1_mix_eta[0]->SetYTitle("v_{1}");
    DiffEPv1_mix_eta[0]->SetMarkerColor(kGreen+2);
    DiffEPv1_mix_eta[0]->SetLineColor(kGreen+2);
    DiffEPv1_mix_eta[0]->SetMarkerStyle(24);
    DiffEPv1_mix_eta[0]->SetMarkerSize(1.2);
    if (isodd) DiffEPv1_mix_eta[0]->GetYaxis()->SetRangeUser(-0.015,0.015);
    else DiffEPv1_mix_eta[0]->GetYaxis()->SetRangeUser(v1in-0.01,v1in+0.01);
    DiffEPv1_mix_eta[0]->Draw("same");

    init_v1in->Draw("hist same c");

    DiffEPv1_mix_3_eta[0]->SetMarkerColor(kBlue);
    DiffEPv1_mix_3_eta[0]->SetLineColor(kBlue);
    DiffEPv1_mix_3_eta[0]->SetMarkerStyle(24);
    DiffEPv1_mix_3_eta[0]->SetMarkerSize(1.2);
    DiffEPv1_mix_3_eta[0]->Draw("same");

    TLegend * legv1EPmix_eta_HFpm = new TLegend(0.21, 0.74, 0.41, 0.90);
    SetLegend(legv1EPmix_eta_HFpm, 18);
    legv1EPmix_eta_HFpm->AddEntry(init_v1in,"Input v_{1}","l");
    legv1EPmix_eta_HFpm->AddEntry(DiffEPv1_mix_eta[0],"v_{1}{#Psi_{1},#Psi_{2}}","p");
    legv1EPmix_eta_HFpm->AddEntry(DiffEPv1_mix_3_eta[0],"v_{1}{#Psi_{1},#Psi_{2},#Psi_{3}}","p");
    legv1EPmix_eta_HFpm->Draw();

    TPaveText * txv1EPmix_eta_HFpm_side = new TPaveText(0.81, 0.84, 0.89, 0.91,"NDC");
    SetTPaveTxt(txv1EPmix_eta_HFpm_side, 20);
    txv1EPmix_eta_HFpm_side->AddText("HF-");
    txv1EPmix_eta_HFpm_side->Draw();


    TPad * padv1EPmix_eta_HFpm_2 = (TPad *) cv1EPmix_eta_HFpm->cd(2);
    if (gridlines) padv1EPmix_eta_HFpm_2->SetGrid();
    DiffEPv1_mix_eta[1]->SetTitle("");
    DiffEPv1_mix_eta[1]->SetStats(0);
    DiffEPv1_mix_eta[1]->SetXTitle("#eta");
    DiffEPv1_mix_eta[1]->SetYTitle("v_{1}");
    DiffEPv1_mix_eta[1]->SetMarkerColor(kGreen+2);
    DiffEPv1_mix_eta[1]->SetLineColor(kGreen+2);
    DiffEPv1_mix_eta[1]->SetMarkerStyle(24);
    DiffEPv1_mix_eta[1]->SetMarkerSize(1.2);
    if (isodd) DiffEPv1_mix_eta[1]->GetYaxis()->SetRangeUser(-0.015,0.015);
    else DiffEPv1_mix_eta[1]->GetYaxis()->SetRangeUser(v1in-0.01,v1in+0.01);
    DiffEPv1_mix_eta[1]->Draw("same");

    init_v1in->Draw("hist same c");

    DiffEPv1_mix_3_eta[1]->SetMarkerColor(kBlue);
    DiffEPv1_mix_3_eta[1]->SetLineColor(kBlue);
    DiffEPv1_mix_3_eta[1]->SetMarkerStyle(24);
    DiffEPv1_mix_3_eta[1]->SetMarkerSize(1.2);
    DiffEPv1_mix_3_eta[1]->Draw("same");

    TPaveText * txv1EPmix_eta_HFpm = new TPaveText(0.22, 0.66, 0.42, 0.90,"NDC");
    SetTPaveTxt(txv1EPmix_eta_HFpm, 18);
    txv1EPmix_eta_HFpm->AddText(Form("N_{events}: %d",NumEvnts));
    if (!isodd) txv1EPmix_eta_HFpm->AddText(Form("Input v_{1}: %0.3f",v1in));
    txv1EPmix_eta_HFpm->AddText(Form("Input v_{2}: %0.3f",v2in));
    txv1EPmix_eta_HFpm->AddText(Form("Input v_{3}: %0.3f",v3in));
    if (eta_weights) txv1EPmix_eta_HFpm->AddText("#eta-weights");
    if (pt_weights) txv1EPmix_eta_HFpm->AddText("p_{T}-weights");
    if (conserve_pT) txv1EPmix_eta_HFpm->AddText("p_{T} conserved");
    txv1EPmix_eta_HFpm->Draw();

    TPaveText * txv1EPmix_eta_HFpm_side_2 = new TPaveText(0.81, 0.84, 0.89, 0.91,"NDC");
    SetTPaveTxt(txv1EPmix_eta_HFpm_side_2, 20);
    txv1EPmix_eta_HFpm_side_2->AddText("HF+");
    txv1EPmix_eta_HFpm_side_2->Draw();

    cv1EPmix_eta_HFpm->Print(Form("plots/%s/v1EPmix_eta_HFpm.png",mtag.Data()),"png");
    if (close_plots) cv1EPmix_eta_HFpm->Close();



    //-- Scalar-product mixed v1(eta) for both HFp and HFm
    TCanvas * cv1SPmix_eta_HFpm = new TCanvas("cv1SPmix_eta_HFpm","cv1SPmix_eta_HFpm",1200,600);
    cv1SPmix_eta_HFpm->Divide(2,1);
    TPad * padv1SPmix_eta_HFpm = (TPad *) cv1SPmix_eta_HFpm->cd(1);
    if (gridlines) padv1SPmix_eta_HFpm->SetGrid();
    DiffSPv1_mix_eta[0]->SetTitle("");
    DiffSPv1_mix_eta[0]->SetStats(0);
    DiffSPv1_mix_eta[0]->SetXTitle("#eta");
    DiffSPv1_mix_eta[0]->SetYTitle("v_{1}");
    DiffSPv1_mix_eta[0]->SetMarkerColor(kGreen+2);
    DiffSPv1_mix_eta[0]->SetLineColor(kGreen+2);
    DiffSPv1_mix_eta[0]->SetMarkerStyle(24);
    DiffSPv1_mix_eta[0]->SetMarkerSize(1.2);
    if (isodd) DiffSPv1_mix_eta[0]->GetYaxis()->SetRangeUser(-0.015,0.015);
    else DiffSPv1_mix_eta[0]->GetYaxis()->SetRangeUser(v1in-0.01,v1in+0.01);
    DiffSPv1_mix_eta[0]->Draw("same");

    init_v1in->Draw("hist same c");

    DiffSPv1_mix_3_eta[0]->SetMarkerColor(kBlue);
    DiffSPv1_mix_3_eta[0]->SetLineColor(kBlue);
    DiffSPv1_mix_3_eta[0]->SetMarkerStyle(24);
    DiffSPv1_mix_3_eta[0]->SetMarkerSize(1.2);
    DiffSPv1_mix_3_eta[0]->Draw("same");

    TLegend * legv1SPmix_eta_HFpm = new TLegend(0.21, 0.74, 0.41, 0.90);
    SetLegend(legv1SPmix_eta_HFpm, 18);
    legv1SPmix_eta_HFpm->AddEntry(init_v1in,"Input v_{1}","l");
    legv1SPmix_eta_HFpm->AddEntry(DiffSPv1_mix_eta[0],"v_{1}{Q_{1},Q_{2}}","p");
    legv1SPmix_eta_HFpm->AddEntry(DiffSPv1_mix_3_eta[0],"v_{1}{Q_{1},Q_{2},Q_{3}}","p");
    legv1SPmix_eta_HFpm->Draw();

    TPaveText * txv1SPmix_eta_HFpm_side = new TPaveText(0.81, 0.84, 0.89, 0.91,"NDC");
    SetTPaveTxt(txv1SPmix_eta_HFpm_side, 20);
    txv1SPmix_eta_HFpm_side->AddText("HF-");
    txv1SPmix_eta_HFpm_side->Draw();


    TPad * padv1SPmix_eta_HFpm_2 = (TPad *) cv1SPmix_eta_HFpm->cd(2);
    if (gridlines) padv1SPmix_eta_HFpm_2->SetGrid();
    DiffSPv1_mix_eta[1]->SetTitle("");
    DiffSPv1_mix_eta[1]->SetStats(0);
    DiffSPv1_mix_eta[1]->SetXTitle("#eta");
    DiffSPv1_mix_eta[1]->SetYTitle("v_{1}");
    DiffSPv1_mix_eta[1]->SetMarkerColor(kGreen+2);
    DiffSPv1_mix_eta[1]->SetLineColor(kGreen+2);
    DiffSPv1_mix_eta[1]->SetMarkerStyle(24);
    DiffSPv1_mix_eta[1]->SetMarkerSize(1.2);
    if (isodd) DiffSPv1_mix_eta[1]->GetYaxis()->SetRangeUser(-0.015,0.015);
    else DiffSPv1_mix_eta[1]->GetYaxis()->SetRangeUser(v1in-0.01,v1in+0.01);
    DiffSPv1_mix_eta[1]->Draw("same");

    init_v1in->Draw("hist same c");

    DiffSPv1_mix_3_eta[1]->SetMarkerColor(kBlue);
    DiffSPv1_mix_3_eta[1]->SetLineColor(kBlue);
    DiffSPv1_mix_3_eta[1]->SetMarkerStyle(24);
    DiffSPv1_mix_3_eta[1]->SetMarkerSize(1.2);
    DiffSPv1_mix_3_eta[1]->Draw("same");

    TPaveText * txv1SPmix_eta_HFpm = new TPaveText(0.22, 0.66, 0.42, 0.90,"NDC");
    SetTPaveTxt(txv1SPmix_eta_HFpm, 18);
    txv1SPmix_eta_HFpm->AddText(Form("N_{events}: %d",NumEvnts));
    if (!isodd) txv1SPmix_eta_HFpm->AddText(Form("Input v_{1}: %0.3f",v1in));
    txv1SPmix_eta_HFpm->AddText(Form("Input v_{2}: %0.3f",v2in));
    txv1SPmix_eta_HFpm->AddText(Form("Input v_{3}: %0.3f",v3in));
    if (eta_weights) txv1SPmix_eta_HFpm->AddText("#eta-weights");
    if (pt_weights) txv1SPmix_eta_HFpm->AddText("p_{T}-weights");
    if (conserve_pT) txv1SPmix_eta_HFpm->AddText("p_{T} conserved");
    txv1SPmix_eta_HFpm->Draw();

    TPaveText * txv1SPmix_eta_HFpm_side_2 = new TPaveText(0.81, 0.84, 0.89, 0.91,"NDC");
    SetTPaveTxt(txv1SPmix_eta_HFpm_side_2, 20);
    txv1SPmix_eta_HFpm_side_2->AddText("HF+");
    txv1SPmix_eta_HFpm_side_2->Draw();

    cv1SPmix_eta_HFpm->Print(Form("plots/%s/v1SPmix_eta_HFpm.png",mtag.Data()),"png");
    if (close_plots) cv1SPmix_eta_HFpm->Close();

}
