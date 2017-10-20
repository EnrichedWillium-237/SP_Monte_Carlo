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

Bool_t print_plots = kTRUE;
Bool_t close_plots = kTRUE;
Bool_t gridlines = kFALSE;

static const int EPNum = 6;

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
TH2D * hpt2D_merged;

TH1D * v1input;

void v1plot()
{

    gStyle->SetPalette(55);

    Int_t Nevents = 1e6;
    Int_t Mult = 6564;
    Double_t v1in = 0.02;
    Double_t v2in = 0.07;
    int subcolor[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7};
    TString HFtag[] = {"HFm", "HFp"};
    TString EPtag[] = {"HFm", "trackm", "trackp", "HFp", "trackmid"};
    gStyle->SetErrorX(0.5);

    TString inFile[20];
    inFile[1]  = "../results/results_fixPsip2_v1_even_0.0000_v2_0.0700_1000000_evts22.root";
    inFile[2]  = "../results/results_fixPsip2_v1_even_0.0000_v2_0.0700_mom-cons_1000000_evts22.root";
    inFile[3]  = "../results/results_fixPsip2_v1_even_0.0000_v2_0.0700_pt_weights_1000000_evts22.root";
    inFile[4]  = "../results/results_fixPsip2_v1_even_0.0000_v2_0.0700_pt_weights_mom-cons_1000000_evts22.root";
    inFile[5]  = "../results/results_fixPsip2_v1_even_0.0200_v2_0.0700_1000000_evts22.root";
    inFile[6]  = "../results/results_fixPsip2_v1_even_0.0200_v2_0.0700_eta_weights_1000000_evts22.root";
    inFile[7]  = "../results/results_fixPsip2_v1_even_0.0200_v2_0.0700_eta_weights_mom-cons_1000000_evts22.root";
    inFile[8]  = "../results/results_fixPsip2_v1_even_0.0200_v2_0.0700_mom-cons_1000000_evts22.root";
    inFile[9]  = "../results/results_fixPsip2_v1_even_0.0200_v2_0.0700_pt_weights_1000000_evts22.root";
    inFile[10] = "../results/results_fixPsip2_v1_even_0.0200_v2_0.0700_pt_weights_mom-cons_1000000_evts22.root";
    inFile[11] = "../results/results_fixPsip2_v1_odd_v2_0.0700_1000000_evts22.root";
    inFile[12] = "../results/results_fixPsip2_v1_odd_v2_0.0700_eta_weights_1000000_evts22.root";

    if (!fopen("plots","r")) system("mkdir plots");

    TFile * tfin01 = new TFile(Form("%s",inFile[1].Data()));
    TDirectory * tdv1_01 = (TDirectory *) tfin01->Get("vn_pt");
    TH1D * v1SP_null_pt = (TH1D *) tdv1_01->Get("HFm/SPv1_3SE_pt_HFm");
    TH1D * v1mix_pt_null = (TH1D *) tdv1_01->Get("HFm/SPv1_mix_pt_HFm");

    TFile * tfin02 = new TFile(Form("%s",inFile[2].Data()));
    TDirectory * tdv1_02 = (TDirectory *) tfin02->Get("vn_pt");
    TH1D * v1SP_null_pt_momcons = (TH1D *) tdv1_02->Get("HFm/SPv1num_pt_HFm");
    TH1D * v1mix_null_pt_momcons = (TH1D *) tdv1_02->Get("HFm/SPv1numMix_pt_HFm");
    TH1D * v1SP_null_pt_momcons_denom = (TH1D *) tfin02->Get("v1SP/HFm/SPdenom1_3SE_HFm");
    v1SP_null_pt_momcons->Scale(1/v1SP_null_pt_momcons_denom->GetMean());
    TH1D * v1mix_null_pt_momcons_denom = (TH1D *) tfin02->Get("v1SP/HFm/SPdenom1_mix_HFm");
    v1mix_null_pt_momcons->Scale(1/v1mix_null_pt_momcons_denom->GetMean());

    TFile * tfin08 = new TFile(Form("%s",inFile[8].Data()));
    TDirectory * tdinput08 = (TDirectory *) tfin08->Get("Inputs");
    TH2D * hpt2D_merged_momcons = (TH2D *) tdinput08->Get("Merged_inputs/pt2D_19");

    TFile * tfin11 = new TFile(Form("%s",inFile[11].Data()));
    TH1D * v1SPodd_HFp_noeweight = (TH1D *) tfin11->Get("vn_eta/HFp/SPv1_2SE_eta_HFp");
    TH1D * v1SPodd_HFm_noeweight = (TH1D *) tfin11->Get("vn_eta/HFm/SPv1_2SE_eta_HFm");
    TH1D * v1Mixodd_HFp_noeweight = (TH1D *) tfin11->Get("vn_eta/HFp/SPv1_mix_eta_HFp");
    TH1D * v1Mixodd_HFm_noeweight = (TH1D *) tfin11->Get("vn_eta/HFm/SPv1_mix_eta_HFm");


    TFile * tfin12 = new TFile(Form("%s",inFile[12].Data()));
    TDirectory * tdinput12 = (TDirectory *) tfin12->Get("Inputs");
    v1input = (TH1D *) tdinput12->Get("v1in");
    etain_merged = (TH1D *) tdinput12->Get("Merged_inputs/etain_19");
    ptin_merged = (TH1D *) tdinput12->Get("Merged_inputs/ptin_19");
    phiPsiRP_merged = (TH1D *) tdinput12->Get("Merged_inputs/phiPsiRP_19");
    hpt2D_merged = (TH2D *) tdinput12->Get("Merged_inputs/pt2D_19");
    TH1D * v1SPodd_HFp_eweight = (TH1D *) tfin12->Get("vn_eta/HFp/EPv1_2SE_eta_HFp");
    TH1D * v1SPodd_HFm_eweight = (TH1D *) tfin12->Get("vn_eta/HFm/EPv1_2SE_eta_HFm");
    TH1D * v1Mixodd_HFp_eweight = (TH1D *) tfin12->Get("vn_eta/HFp/SPv1_mix_eta_HFp");
    TH1D * v1Mixodd_HFm_eweight = (TH1D *) tfin12->Get("vn_eta/HFm/SPv1_mix_eta_HFm");

    for (int nep = 0; nep<EPNum; nep++) {
        if (nep == EPNum) continue;
        philab_sub[nep] = (TH1D *) tdinput12->Get(Form("Subevent_inputs/philab_ep%d_19",nep));
        etain_sub[nep] = (TH1D *) tdinput12->Get(Form("Subevent_inputs/etain_ep%d_19",nep));
        ptin_sub[nep] = (TH1D *) tdinput12->Get(Form("Subevent_inputs/ptin_ep%d_19",nep));
        phiPsiRP[nep] = (TH1D *) tdinput12->Get(Form("Subevent_inputs/phiPsiRP_ep%d_19",nep));
        phiPsi1[nep] = (TH1D *) tdinput12->Get(Form("Subevent_inputs/phiPsi1_ep%d_19",nep));
        phiPsi2[nep] = (TH1D *) tdinput12->Get(Form("Subevent_inputs/phiPsi2_ep%d_19",nep));
        hsub_Psi1lab[nep] = (TH1D *) tdinput12->Get(Form("Subevent_inputs/Psi1_ep%d_19",nep));
        hsub_Psi2lab[nep] = (TH1D *) tdinput12->Get(Form("Subevent_inputs/Psi2_ep%d_19",nep));
        pt2D_sub[nep] = (TH2D *) tdinput12->Get(Form("Subevent_inputs/pt2D_ep%d_19",nep));

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
    }



    v1input->SetMarkerColor(kBlack);
    v1input->SetLineColor(kBlack);
    v1input->SetMarkerStyle(20);
    v1input->SetMarkerSize(0.7);
    v1input->SetTitle("");
    v1input->SetStats(kFALSE);
    v1input->SetXTitle("#eta");
    v1input->SetYTitle("Input v_{1}");
    v1input->GetXaxis()->CenterTitle();
    v1input->GetXaxis()->SetTitleSize(0.06);
    v1input->GetYaxis()->SetTitleSize(0.06);
    v1input->GetXaxis()->SetTitleOffset(1.10);
    v1input->GetYaxis()->SetTitleOffset(1.30);
    v1input->GetXaxis()->SetLabelSize(0.06);
    v1input->GetYaxis()->SetLabelSize(0.05);

    etain_merged->SetMarkerColor(kBlack);
    etain_merged->SetLineColor(kBlack);
    etain_merged->SetMarkerStyle(20);
    etain_merged->SetMarkerSize(0.7);
    etain_merged->SetTitle("");
    etain_merged->SetStats(kFALSE);
    etain_merged->SetXTitle("#eta");
    etain_merged->GetXaxis()->CenterTitle();
    etain_merged->GetXaxis()->SetTitleSize(0.06);
    etain_merged->GetXaxis()->SetTitleOffset(1.10);
    etain_merged->GetXaxis()->SetLabelSize(0.06);
    etain_merged->GetYaxis()->SetLabelSize(0.05);

    ptin_merged->SetMarkerColor(kBlack);
    ptin_merged->SetLineColor(kBlack);
    ptin_merged->SetMarkerStyle(20);
    ptin_merged->SetMarkerSize(0.7);
    ptin_merged->SetTitle("");
    ptin_merged->SetStats(kFALSE);
    ptin_merged->SetXTitle("p_{T} (GeV/c)");
    ptin_merged->GetXaxis()->CenterTitle();
    ptin_merged->GetXaxis()->SetTitleSize(0.06);
    ptin_merged->GetXaxis()->SetTitleOffset(1.10);
    ptin_merged->GetXaxis()->SetLabelSize(0.06);
    ptin_merged->GetYaxis()->SetLabelSize(0.05);

    phiPsiRP_merged->SetMarkerColor(kBlack);
    phiPsiRP_merged->SetLineColor(kBlack);
    phiPsiRP_merged->SetMarkerStyle(20);
    phiPsiRP_merged->SetMarkerSize(0.7);
    phiPsiRP_merged->SetTitle("");
    phiPsiRP_merged->SetStats(kFALSE);
    phiPsiRP_merged->SetXTitle("#phi - #Psi_{RP}");
    phiPsiRP_merged->GetXaxis()->CenterTitle();
    phiPsiRP_merged->GetXaxis()->SetTitleSize(0.06);
    phiPsiRP_merged->GetXaxis()->SetTitleOffset(1.10);
    phiPsiRP_merged->GetXaxis()->SetLabelSize(0.06);
    phiPsiRP_merged->GetYaxis()->SetLabelSize(0.05);

    hpt2D_merged->SetTitle("");
    hpt2D_merged->SetStats(0);
    hpt2D_merged->GetXaxis()->SetRangeUser(-250,250);
    hpt2D_merged->GetYaxis()->SetRangeUser(-250,250);
    hpt2D_merged->SetXTitle("#Sigma cos(p_{T})  (GeV/c)");
    hpt2D_merged->SetYTitle("#Sigma sin(p_{T})  (GeV/c)");
    hpt2D_merged->GetXaxis()->CenterTitle();
    hpt2D_merged->GetYaxis()->CenterTitle();
    hpt2D_merged->GetXaxis()->SetTitleSize(0.06);
    hpt2D_merged->GetYaxis()->SetTitleSize(0.06);
    hpt2D_merged->GetXaxis()->SetTitleOffset(1.10);
    hpt2D_merged->GetYaxis()->SetTitleOffset(1.25);
    hpt2D_merged->GetXaxis()->SetLabelSize(0.04);
    hpt2D_merged->GetYaxis()->SetLabelSize(0.04);

    hpt2D_merged_momcons->SetTitle("");
    hpt2D_merged_momcons->SetStats(0);
    hpt2D_merged_momcons->GetXaxis()->SetRangeUser(-250,250);
    hpt2D_merged_momcons->GetYaxis()->SetRangeUser(-250,250);
    hpt2D_merged_momcons->SetXTitle("#Sigma cos(p_{T})  (GeV/c)");
    hpt2D_merged_momcons->SetYTitle("#Sigma sin(p_{T})  (GeV/c)");
    hpt2D_merged_momcons->GetXaxis()->CenterTitle();
    hpt2D_merged_momcons->GetYaxis()->CenterTitle();
    hpt2D_merged_momcons->GetXaxis()->SetTitleSize(0.06);
    hpt2D_merged_momcons->GetYaxis()->SetTitleSize(0.06);
    hpt2D_merged_momcons->GetXaxis()->SetTitleOffset(1.10);
    hpt2D_merged_momcons->GetYaxis()->SetTitleOffset(1.25);
    hpt2D_merged_momcons->GetXaxis()->SetLabelSize(0.04);
    hpt2D_merged_momcons->GetYaxis()->SetLabelSize(0.04);


    v1SP_null_pt->SetMarkerColor(kBlue);
    v1SP_null_pt->SetLineColor(kBlue);
    v1SP_null_pt->SetMarkerStyle(21);
    v1SP_null_pt->SetMarkerSize(1.1);

    v1mix_pt_null->SetMarkerColor(kBlue);
    v1mix_pt_null->SetLineColor(kBlue);
    v1mix_pt_null->SetMarkerStyle(25);
    v1mix_pt_null->SetMarkerSize(1.1);

    v1SP_null_pt_momcons->SetMarkerColor(kBlue);
    v1SP_null_pt_momcons->SetLineColor(kBlue);
    v1SP_null_pt_momcons->SetMarkerStyle(21);
    v1SP_null_pt_momcons->SetMarkerSize(1.1);

    v1mix_null_pt_momcons->SetMarkerColor(kBlue);
    v1mix_null_pt_momcons->SetLineColor(kBlue);
    v1mix_null_pt_momcons->SetMarkerStyle(25);
    v1mix_null_pt_momcons->SetMarkerSize(1.1);

    v1SPodd_HFp_noeweight->SetMarkerColor(kBlue);
    v1SPodd_HFp_noeweight->SetLineColor(kBlue);
    v1SPodd_HFp_noeweight->SetMarkerStyle(21);
    v1SPodd_HFp_noeweight->SetMarkerSize(1.1);

    v1SPodd_HFm_noeweight->SetMarkerColor(kBlue);
    v1SPodd_HFm_noeweight->SetLineColor(kBlue);
    v1SPodd_HFm_noeweight->SetMarkerStyle(21);
    v1SPodd_HFm_noeweight->SetMarkerSize(1.1);

    v1Mixodd_HFp_noeweight->SetMarkerColor(kBlue);
    v1Mixodd_HFp_noeweight->SetLineColor(kBlue);
    v1Mixodd_HFp_noeweight->SetMarkerStyle(25);
    v1Mixodd_HFp_noeweight->SetMarkerSize(1.1);

    v1Mixodd_HFm_noeweight->SetMarkerColor(kBlue);
    v1Mixodd_HFm_noeweight->SetLineColor(kBlue);
    v1Mixodd_HFm_noeweight->SetMarkerStyle(25);
    v1Mixodd_HFm_noeweight->SetMarkerSize(1.1);

    v1SPodd_HFp_eweight->SetMarkerColor(kBlue);
    v1SPodd_HFp_eweight->SetLineColor(kBlue);
    v1SPodd_HFp_eweight->SetMarkerStyle(21);
    v1SPodd_HFp_eweight->SetMarkerSize(1.1);

    v1SPodd_HFm_eweight->SetMarkerColor(kBlue);
    v1SPodd_HFm_eweight->SetLineColor(kBlue);
    v1SPodd_HFm_eweight->SetMarkerStyle(21);
    v1SPodd_HFm_eweight->SetMarkerSize(1.1);

    v1Mixodd_HFp_eweight->SetMarkerColor(kBlue);
    v1Mixodd_HFp_eweight->SetLineColor(kBlue);
    v1Mixodd_HFp_eweight->SetMarkerStyle(25);
    v1Mixodd_HFp_eweight->SetMarkerSize(1.1);

    v1Mixodd_HFm_eweight->SetMarkerColor(kBlue);
    v1Mixodd_HFm_eweight->SetLineColor(kBlue);
    v1Mixodd_HFm_eweight->SetMarkerStyle(25);
    v1Mixodd_HFm_eweight->SetMarkerSize(1.1);



    //-- Monte Carlo inputs for merged events
    TCanvas * cMCinputs_merged = new TCanvas("cMCinputs_merged","cMCinputs_merged",850,750);
    cMCinputs_merged->Divide(2,2);

    cMCinputs_merged->cd(1);
    v1input->Draw();
    TPaveText * txMCinputs_merged_1 = new TPaveText(0.2, 0.75, 0.59, 0.89, "NDC");
    txMCinputs_merged_1->SetFillColor(0);
    txMCinputs_merged_1->SetBorderSize(0);
    txMCinputs_merged_1->SetTextFont(43);
    txMCinputs_merged_1->SetTextSize(20);
    txMCinputs_merged_1->SetTextAlign(12);
    txMCinputs_merged_1->AddText("Rapidty-odd v_{1}");
    txMCinputs_merged_1->AddText("Constant v_{2} = 0.07");
    txMCinputs_merged_1->Draw();

    cMCinputs_merged->cd(2);
    etain_merged->Draw();

    TPad * padMCinputs_merged = (TPad *) cMCinputs_merged->cd(3);
    padMCinputs_merged->SetLogy();
    ptin_merged->Draw();

    cMCinputs_merged->cd(4);
    phiPsiRP_merged->GetYaxis()->SetRangeUser(0, 3800e3);
    phiPsiRP_merged->Draw();

    if (print_plots) cMCinputs_merged->Print("plots/MCinputs_merged.png","png");
    if (close_plots) cMCinputs_merged->Close();



    //-- 2D plots before and after pT-conservation

    TCanvas * cpt2D_merged = new TCanvas("cpt2D_merged","cpt2D_merged",500,500);
    hpt2D_merged->Draw("col");
    if (print_plots) cpt2D_merged->Print("plots/pt2D_merged.png","png");
    if (close_plots) cpt2D_merged->Close();

    TCanvas * cpt2D_merged_momcons = new TCanvas("cpt2D_merged_momcons","cpt2D_merged_momcons",500,500);
    hpt2D_merged_momcons->Draw("col");
    if (print_plots) cpt2D_merged_momcons->Print("plots/pt2D_merged_momcons.png","png");
    if (close_plots) cpt2D_merged_momcons->Close();



    //-- plot v1 = 0 with and without mom-cons

    TCanvas * cv1null_pt = new TCanvas("cv1null_pt","cv1null_pt",500,500);
    cv1null_pt->cd();
    TH1D * hv1null_pt = new TH1D("hv1null_pt", "", 50, 0, 8);
    hv1null_pt->SetTitle("");
    hv1null_pt->SetStats(kFALSE);
    hv1null_pt->SetXTitle("p_{T} (GeV/c)");
    hv1null_pt->SetYTitle("v_{1}");
    hv1null_pt->GetXaxis()->CenterTitle(kTRUE);
    hv1null_pt->GetYaxis()->CenterTitle(kTRUE);
    hv1null_pt->GetYaxis()->SetRangeUser(-0.08, 0.08);
    hv1null_pt->GetXaxis()->SetTitleSize(0.05);
    hv1null_pt->GetYaxis()->SetTitleSize(0.05);
    hv1null_pt->GetYaxis()->SetTitleOffset(1.25);
    hv1null_pt->Draw();
    v1SP_null_pt->Draw("same");
    v1mix_pt_null->Draw("same");

    TPaveText * txv1null_pt = new TPaveText(0.64, 0.74, 0.91, 0.91, "NDC");
    txv1null_pt->SetFillColor(0);
    txv1null_pt->SetBorderSize(0);
    txv1null_pt->SetTextFont(43);
    txv1null_pt->SetTextSize(20);
    txv1null_pt->SetTextAlign(12);
    txv1null_pt->AddText("10^{6} MC events");
    txv1null_pt->AddText("Input v_{1} = 0");
    txv1null_pt->AddText("Input v_{1} = 0.07");
    txv1null_pt->Draw();

    TPaveText * txv1null_pt_1 = new TPaveText(0.21, 0.85, 0.48, 0.92, "NDC");
    txv1null_pt_1->SetFillColor(0);
    txv1null_pt_1->SetBorderSize(0);
    txv1null_pt_1->SetTextFont(43);
    txv1null_pt_1->SetTextSize(20);
    txv1null_pt_1->SetTextAlign(12);
    txv1null_pt_1->AddText("p_{T} not conserved");
    txv1null_pt_1->Draw();

    if (print_plots) cv1null_pt->Print("plots/v1null_pt_pt.png","png");
    if (close_plots) cv1null_pt->Close();



    TCanvas * cv1null_pt_momcons = new TCanvas("cv1null_pt_momcons","cv1null_pt_momcons",500,500);
    cv1null_pt_momcons->cd();
    TH1D * hv1null_pt_momcons = (TH1D *) hv1null_pt->Clone("hv1null_pt_momcons");
    // hv1null_pt_momcons->GetYaxis()->SetRangeUser(-0.5, 0.2);
    hv1null_pt_momcons->GetYaxis()->SetRangeUser(-0.08, 0.08);
    hv1null_pt_momcons->Draw();
    v1SP_null_pt_momcons->Draw("same");
    v1mix_null_pt_momcons->Draw("same");

    TLegend * legv1null_pt_momcons = new TLegend(0.59, 0.80, 0.80, 0.91);
    legv1null_pt_momcons->SetFillColor(0);
    legv1null_pt_momcons->SetBorderSize(0);
    legv1null_pt_momcons->SetTextFont(43);
    legv1null_pt_momcons->SetTextSize(20);
    legv1null_pt_momcons->AddEntry(v1SP_null_pt_momcons,"Scalar-product","p");
    legv1null_pt_momcons->AddEntry(v1mix_null_pt_momcons,"Mixed harmonic","p");
    legv1null_pt_momcons->Draw();

    TPaveText * txv1null_pt_momcons_1 = new TPaveText(0.21, 0.85, 0.48, 0.92, "NDC");
    txv1null_pt_momcons_1->SetFillColor(0);
    txv1null_pt_momcons_1->SetBorderSize(0);
    txv1null_pt_momcons_1->SetTextFont(43);
    txv1null_pt_momcons_1->SetTextSize(20);
    txv1null_pt_momcons_1->SetTextAlign(12);
    txv1null_pt_momcons_1->AddText("p_{T} conserved");
    txv1null_pt_momcons_1->Draw();

    if (print_plots) cv1null_pt_momcons->Print("plots/v1null_pt_momcons.png","png");
    if (close_plots) cv1null_pt_momcons->Close();



    //-- v1odd with and without eta-weighting

    TCanvas * cv1odd_HFp_noeweight = new TCanvas("cv1odd_HFp_noeweight","cv1odd_HFp_noeweight",500,500);
    cv1odd_HFp_noeweight->cd();
    TH1D * hv1odd_HFp_noeweight = new TH1D("hv1odd_HFp_noeweight", "", 50, -2.4, 2.4);
    hv1odd_HFp_noeweight->SetTitle("");
    hv1odd_HFp_noeweight->SetStats(kFALSE);
    hv1odd_HFp_noeweight->SetXTitle("#eta");
    hv1odd_HFp_noeweight->SetYTitle("v_{1}");
    hv1odd_HFp_noeweight->GetXaxis()->CenterTitle(kTRUE);
    hv1odd_HFp_noeweight->GetYaxis()->CenterTitle(kTRUE);
    hv1odd_HFp_noeweight->GetYaxis()->SetRangeUser(-0.01, 0.01);
    hv1odd_HFp_noeweight->GetXaxis()->SetTitleSize(0.05);
    hv1odd_HFp_noeweight->GetYaxis()->SetTitleSize(0.05);
    hv1odd_HFp_noeweight->GetYaxis()->SetTitleOffset(1.25);
    hv1odd_HFp_noeweight->Draw();
    v1input->Draw("same");
    v1SPodd_HFp_noeweight->Draw("same");
    v1Mixodd_HFp_noeweight->Draw("same");

    TPaveText * txv1odd_HFp_noeweight = new TPaveText(0.22, 0.84, 0.33, 0.91, "NDC");
    txv1odd_HFp_noeweight->SetFillColor(0);
    txv1odd_HFp_noeweight->SetBorderSize(0);
    txv1odd_HFp_noeweight->SetTextFont(43);
    txv1odd_HFp_noeweight->SetTextSize(26);
    txv1odd_HFp_noeweight->SetTextColor(kRed);
    txv1odd_HFp_noeweight->AddText("v_{1}^{HF+}");
    txv1odd_HFp_noeweight->Draw();

    TLegend * legv1odd_HFp_noeweight = new TLegend(0.2, 0.18, 0.5, 0.35);
    legv1odd_HFp_noeweight->SetFillColor(0);
    legv1odd_HFp_noeweight->SetBorderSize(0);
    legv1odd_HFp_noeweight->SetTextFont(43);
    legv1odd_HFp_noeweight->SetTextSize(20);
    legv1odd_HFp_noeweight->AddEntry(v1input,"Input v_{1}","p");
    legv1odd_HFp_noeweight->AddEntry(v1SPodd_HFp_noeweight,"Scalar-product","p");
    legv1odd_HFp_noeweight->AddEntry(v1Mixodd_HFp_noeweight,"Mixed harmonic","p");
    legv1odd_HFp_noeweight->Draw();

    if (print_plots) cv1odd_HFp_noeweight->Print("plots/v1odd_HFp_noeweight.png","png");
    if (close_plots) cv1odd_HFp_noeweight->Close();



    TCanvas * cv1odd_HFm_noeweight = new TCanvas("cv1odd_HFm_noeweight","cv1odd_HFm_noeweight",500,500);
    cv1odd_HFm_noeweight->cd();
    TH1D * hv1odd_HFm_noeweight = (TH1D *) hv1odd_HFp_noeweight->Clone("hv1odd_HFm_noeweight");
    hv1odd_HFm_noeweight->Draw();
    v1input->Draw("same");
    v1SPodd_HFm_noeweight->Draw("same");
    v1Mixodd_HFm_noeweight->Draw("same");

    TPaveText * txv1odd_HFm_noeweight = new TPaveText(0.22, 0.84, 0.33, 0.91, "NDC");
    txv1odd_HFm_noeweight->SetFillColor(0);
    txv1odd_HFm_noeweight->SetBorderSize(0);
    txv1odd_HFm_noeweight->SetTextFont(43);
    txv1odd_HFm_noeweight->SetTextSize(26);
    txv1odd_HFm_noeweight->SetTextColor(kCyan+2);
    txv1odd_HFm_noeweight->AddText("v_{1}^{HF-}");
    txv1odd_HFm_noeweight->Draw();

    if (print_plots) cv1odd_HFm_noeweight->Print("plots/v1odd_HFm_noeweight.png","png");
    if (close_plots) cv1odd_HFm_noeweight->Close();



    TCanvas * cv1odd_HFp_eweight = new TCanvas("cv1odd_HFp_eweight","cv1odd_HFp_eweight",500,500);
    cv1odd_HFp_eweight->cd();
    TH1D * hv1odd_HFp_eweight = (TH1D *) hv1odd_HFp_noeweight->Clone("hv1odd_HFp_eweight");
    hv1odd_HFp_eweight->Draw();
    v1input->Draw("same");
    v1SPodd_HFp_eweight->Draw("same");
    v1Mixodd_HFp_eweight->Draw("same");

    TPaveText * txv1odd_HFp_eweight = new TPaveText(0.22, 0.84, 0.33, 0.91, "NDC");
    txv1odd_HFp_eweight->SetFillColor(0);
    txv1odd_HFp_eweight->SetBorderSize(0);
    txv1odd_HFp_eweight->SetTextFont(43);
    txv1odd_HFp_eweight->SetTextSize(26);
    txv1odd_HFp_eweight->SetTextColor(kRed);
    txv1odd_HFp_eweight->AddText("v_{1}^{HF+}");
    txv1odd_HFp_eweight->Draw();

    if (print_plots) cv1odd_HFp_eweight->Print("plots/v1odd_HFp_eweight.png","png");
    if (close_plots) cv1odd_HFp_eweight->Close();



    TCanvas * cv1odd_HFm_eweight = new TCanvas("cv1odd_HFm_eweight","cv1odd_HFm_eweight",500,500);
    cv1odd_HFm_eweight->cd();
    TH1D * hv1odd_HFm_eweight = (TH1D *) hv1odd_HFp_noeweight->Clone("hv1odd_HFm_eweight");
    hv1odd_HFm_eweight->Draw();
    v1input->Draw("same");
    v1SPodd_HFm_eweight->Draw("same");
    v1Mixodd_HFm_eweight->Draw("same");

    TPaveText * txv1odd_HFm_eweight = new TPaveText(0.22, 0.84, 0.33, 0.91, "NDC");
    txv1odd_HFm_eweight->SetFillColor(0);
    txv1odd_HFm_eweight->SetBorderSize(0);
    txv1odd_HFm_eweight->SetTextFont(43);
    txv1odd_HFm_eweight->SetTextSize(26);
    txv1odd_HFm_eweight->SetTextColor(kCyan+2);
    txv1odd_HFm_eweight->AddText("v_{1}^{HF-}");
    txv1odd_HFm_eweight->Draw();

    if (print_plots) cv1odd_HFm_eweight->Print("plots/v1odd_HFm_eweight.png","png");
    if (close_plots) cv1odd_HFm_eweight->Close();



    //-- compare v1non-flow with data

    # include "dataCompare.h"
    TGraphErrors * gdata_momcons_pt_compare = new TGraphErrors(ndataMomCons, dataMomCons_xval, dataMomCons_yval_c40t50, 0, dataMomCons_yvalErr_c40t50);
    gdata_momcons_pt_compare->SetMarkerColor(kBlue);
    gdata_momcons_pt_compare->SetLineColor(kBlue);
    gdata_momcons_pt_compare->SetMarkerStyle(21);
    gdata_momcons_pt_compare->SetMarkerSize(1.1);

    TCanvas * cdata_momcons_pt_compare = new TCanvas("cdata_momcons_pt_compare","cdata_momcons_pt_compare",500,500);
    cdata_momcons_pt_compare->cd();
    TH1D * hdata_momcons_pt_compare = (TH1D *) hv1null_pt->Clone("hdata_momcons_pt_compare");
    // hdata_momcons_pt_compare->GetYaxis()->SetRangeUser(-0.5, 0.2);
    hdata_momcons_pt_compare->GetXaxis()->SetRangeUser(0.0, 4.0);
    hdata_momcons_pt_compare->GetYaxis()->SetRangeUser(-0.3, 0.15);
    hdata_momcons_pt_compare->SetYTitle("v_{1}^{non-flow}");
    hdata_momcons_pt_compare->GetYaxis()->SetTitleOffset(1.6);
    hdata_momcons_pt_compare->GetYaxis()->CenterTitle(kFALSE);
    hdata_momcons_pt_compare->Draw();
    v1SP_null_pt_momcons->SetMarkerColor(kOrange+7);
    v1SP_null_pt_momcons->SetLineColor(kOrange+7);
    v1SP_null_pt_momcons->SetMarkerStyle(20);
    v1SP_null_pt_momcons->SetMarkerSize(1.2);
    v1SP_null_pt_momcons->Draw("same");
    gdata_momcons_pt_compare->Draw("same p");

    TLegend * legdata_momcons_pt_compare = new TLegend(0.40, 0.75, 0.78, 0.92);
    legdata_momcons_pt_compare->SetFillColor(0);
    legdata_momcons_pt_compare->SetBorderSize(0);
    legdata_momcons_pt_compare->SetTextFont(43);
    legdata_momcons_pt_compare->SetTextSize(20);
    legdata_momcons_pt_compare->SetHeader("40-50% Centrality");
    legdata_momcons_pt_compare->AddEntry(gdata_momcons_pt_compare,"CMS PbPb 5.02 TeV","p");
    legdata_momcons_pt_compare->AddEntry(v1SP_null_pt_momcons,"Toy MC","p");
    legdata_momcons_pt_compare->Draw();

    if (print_plots) cdata_momcons_pt_compare->Print("plots/data_momcons_pt_compare.png","png");
    //if (close_plots) cdata_momcons_pt_compare->Close();


}
