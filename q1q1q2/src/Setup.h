# include "TCanvas.h"
# include "TDatime.h"
# include "TFile.h"
# include "TH1.h"
# include "TH2.h"
# include "TLegend.h"
# include "TPaveText.h"
# include "TStopwatch.h"
# include "TString.h"
# include "TStyle.h"
# include <cmath>
# include <complex>
# include <fstream>
# include <iostream>

# include "MCEvent.h"

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
static const double ecutmin[] = {-5.0, -0.4,  0.0,  4.0, -0.5};
static const double ecutmax[] = {-4.0,  0.0,  0.4,  5.0,  0.5};
static const double pcutmin[] = { 0.3,  0.3,  0.3,  0.3,  0.3};
static const double pcutmax[] = {30.0,  3.0,  3.0, 30.0,  3.0};

static const int nv1Ebins = 28;
static const double v1Ebins[] = {
    -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0,
    -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,
     2.4,  2.8,  3.2,  3.6,  4.0,  4.4,  4.8,  5.2,  5.6};

static const double v1oddin[] = {
     0.00378,  0.00350,  0.00322,  0.00294,  0.00266,  0.00238,
     0.00210,  0.00182,  0.00154,  0.00126,  0.00098,  0.00070,
     0.00042,  0.00014, -0.00014, -0.00042, -0.00070, -0.00098,
    -0.00126, -0.00154, -0.00182, -0.00210, -0.00238, -0.00266,
    -0.00294, -0.00322, -0.00350, -0.00378};

static const int v1Emult[] = {
    129,  149,  170,  191,  213,  234,
    253,  270,  281,  286,  285,  280,
    273,  268,  268,  273,  280,  285,
    286,  281,  270,  253,  234,  213,
    191,  170,  149,  129};

static const Double_t phiMinHole[] = {0.0,  2.2,  0.0,  0.0,  0.0};
static const Double_t phiMaxHole[] = {0.0,  2.6,  0.0,  0.0,  0.0};
static const Double_t etaMinHole[] = {0.0, -2.0,  0.0,  0.0,  0.0};
static const Double_t etaMaxHole[] = {0.0, -0.4,  0.0,  0.0,  0.0};

static const Int_t multMax = 8000;


//-- initial MC throw

TH1D * inputParms;
TH1D * hinit_v1in;
TH1D * hinit_philab[nv1Ebins];
TH1D * hinit_eta[nv1Ebins];
TH1D * hinit_pt[nv1Ebins];
TH1D * hinit_phiPsiRP[nv1Ebins];

//-- merge MC inputs into a single event

TH1D * hphilab;
TH1D * hetain;
TH1D * hptin;
TH1D * hphiPsiRP;
TH1D * hptAve2;
TH1D * hweights;
TH1D * hwpt;
TH2D * hpt2D;

//-- divide into subevents

TH1D * hsub_philab[numEP];
TH1D * hsub_etain[numEP];
TH1D * hsub_ptin[numEP];
TH1D * hsub_phiPsiRP[numEP];
TH1D * hsub_phiPsi1[numEP];
TH1D * hsub_phiPsi2[numEP];
TH1D * hsub_Psi1[numEP];
TH1D * hsub_Psi2[numEP];
TH1D * hsub_Psi1lab[numEP];
TH1D * hsub_Psi2lab[numEP];

TH1D * hq112;

TFile * tfout;
Int_t iseed = 0;
Int_t counter = 0;


Double_t bounds(int ord, double ang) {
    while (ang >  TMath::Pi()/ord) ang-=TMath::TwoPi()/ord;
    while (ang < -TMath::Pi()/ord) ang+=TMath::TwoPi()/ord;
    return ang;
}

void Setup()
{

    hq112 = new TH1D("q112", "", 200, -0.1, 0.2);

}


void WriteToFile( TString mtag ) {
    if (!fopen("results","r")) system("mkdir results");
    tfout = new TFile(Form("results/results_%s.root",mtag.Data()),"recreate");

    TDirectory * tdinput = (TDirectory *) tfout->mkdir("Inputs");
    tdinput->cd();
    inputParms->Write();
    hinit_v1in->Write();

    TDirectory * tdinit = (TDirectory *) tdinput->mkdir("Initial_throw");
    tdinit->cd();
    for (int vbin = 0; vbin<nv1Ebins; vbin++) hinit_philab[vbin]->Write();
    for (int vbin = 0; vbin<nv1Ebins; vbin++) hinit_eta[vbin]->Write();
    for (int vbin = 0; vbin<nv1Ebins; vbin++) hinit_pt[vbin]->Write();
    for (int vbin = 0; vbin<nv1Ebins; vbin++) hinit_phiPsiRP[vbin]->Write();

    TDirectory * tdmerged = (TDirectory *) tdinput->mkdir("Merged_inputs");
    tdmerged->cd();
    hphilab->Write();
    hetain->Write();
    hptin->Write();
    hphiPsiRP->Write();
    hptAve2->Write();
    hweights->Write();
    hwpt->Write();
    hpt2D->Write();

    TDirectory * tdsubevnt = (TDirectory *) tdinput->mkdir("Subevent_inputs");
    tdsubevnt->cd();
    for (int nep = 0; nep<numEP; nep++) hsub_philab[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_etain[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_ptin[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_phiPsiRP[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_phiPsi1[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_phiPsi2[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_phiPsi3[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi1[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi2[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi3[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi1lab[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi2lab[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi3lab[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_ptAve[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_ptAve2[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_weights[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_wpt[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_pt2D[nep]->Write();

    TDirectory * tdqvec = (TDirectory *) tfout->mkdir("Q-Vectors");
    TDirectory * tdqvec_side[2];
    TDirectory * tdqvecNorm = (TDirectory *) tfout->mkdir("Normalized_Q-Vectors");
    TDirectory * tdqvecNorm_side[2];
    TDirectory * tdEPCorr = (TDirectory *) tfout->mkdir("EP_correlations");
    TDirectory * tdEPCos = (TDirectory *) tfout->mkdir("EP_cosines");
    TDirectory * tdvnep = (TDirectory *) tfout->mkdir("v1EP");
    TDirectory * tdvnep_side[2];
    TDirectory * tdvnsp = (TDirectory *) tfout->mkdir("v1SP");
    TDirectory * tdvnsp_side[2];
    TDirectory * tdvnDiff_pt = (TDirectory *) tfout->mkdir("vn_pt");
    TDirectory * tdvnDiff_pt_side[2];
    TDirectory * tdvnDiff_eta = (TDirectory *) tfout->mkdir("vn_eta");
    TDirectory * tdvnDiff_eta_side[2];

    for (int iside = 0; iside<2; iside++) {
        tdqvec_side[iside] = (TDirectory *) tdqvec->mkdir(Form("%s",SideName[iside].data()));
        tdqvec_side[iside]->cd();
        hQ1nA_final[iside]->Write();
        hQ1AB_final[iside]->Write();
        hQ1AC_final[iside]->Write();
        hQ1BC_final[iside]->Write();
        hQ2AB_final[iside]->Write();
        hQ2AC_final[iside]->Write();
        hQ2BC_final[iside]->Write();
        hQ3AB_final[iside]->Write();
        hQ3AC_final[iside]->Write();
        hQ3BC_final[iside]->Write();
        hQmixnAC_final[iside]->Write();
        hQmixABC_final[iside]->Write();
        hQmixnAC_3_final[iside]->Write();
        hQmixABC_3_final[iside]->Write();
        hQ1Q1Q2[iside]->Write();

        tdqvecNorm_side[iside] = (TDirectory *) tdqvecNorm->mkdir(Form("%s",SideName[iside].data()));
        tdqvecNorm_side[iside]->cd();
        hQ1nAnorm_final[iside]->Write();
        hQ1ABnorm_final[iside]->Write();
        hQ1ACnorm_final[iside]->Write();
        hQ1BCnorm_final[iside]->Write();
        hQ2ABnorm_final[iside]->Write();
        hQ2ACnorm_final[iside]->Write();
        hQ2BCnorm_final[iside]->Write();
        hQ3ABnorm_final[iside]->Write();
        hQ3ACnorm_final[iside]->Write();
        hQ3BCnorm_final[iside]->Write();
        hQmixnACnorm_final[iside]->Write();
        hQmixABCnorm_final[iside]->Write();
        hQmixnAC_3norm_final[iside]->Write();
        hQmixABC_3norm_final[iside]->Write();
        hQ1Q1Q2norm[iside]->Write();
        hq112->Write();

        tdvnep_side[iside] = (TDirectory *) tdvnep->mkdir(Form("%s",SideName[iside].data()));
        tdvnep_side[iside]->cd();
        rescor1_2SE->Write();
        rescor1_3SE[iside]->Write();
        rescor1_mix[iside]->Write();
        rescor1_mix3[iside]->Write();
        rescor2HF_2SE->Write();
        rescor2HF_3SE[iside]->Write();
        rescor2Trk_3SE[iside]->Write();
        rescor3HF_2SE->Write();
        rescor3HF_3SE[iside]->Write();
        rescor3Trk_3SE[iside]->Write();
        EPv1obs[iside]->Write();
        EPv1obsMix[iside]->Write();
        EPv1obsMix_3[iside]->Write();
        EPv2obs[iside]->Write();
        EPv3obs[iside]->Write();
        EPv1_2SE[iside]->Write();
        EPv1_3SE[iside]->Write();
        EPv1_mix[iside]->Write();
        EPv1_mix_3[iside]->Write();
        EPv2_2SE[iside]->Write();
        EPv2_3SE[iside]->Write();
        EPv3_2SE[iside]->Write();
        EPv3_3SE[iside]->Write();

        tdvnsp_side[iside] = (TDirectory *) tdvnsp->mkdir(Form("%s",SideName[iside].data()));
        tdvnsp_side[iside]->cd();
        SPdenom1_2SE->Write();
        SPdenom1_3SE[iside]->Write();
        SPdenom1_mix[iside]->Write();
        SPdenom1_mix_3[iside]->Write();
        SPdenom2HF_2SE->Write();
        SPdenom2HF_3SE[iside]->Write();
        SPdenom2Trk_3SE[iside]->Write();
        SPdenom3HF_2SE->Write();
        SPdenom3HF_3SE[iside]->Write();
        SPdenom3Trk_3SE[iside]->Write();
        SPv1num[iside]->Write();
        SPv1numMix[iside]->Write();
        SPv1numMix_3[iside]->Write();
        SPv2num[iside]->Write();
        SPv3num[iside]->Write();
        SPv1_2SE[iside]->Write();
        SPv1_3SE[iside]->Write();
        SPv1_mix[iside]->Write();
        SPv1_mix_3[iside]->Write();
        SPv2_2SE[iside]->Write();
        SPv2_3SE[iside]->Write();
        SPv3_2SE[iside]->Write();
        SPv3_3SE[iside]->Write();

        tdvnDiff_pt_side[iside] = (TDirectory *) tdvnDiff_pt->mkdir(Form("%s",SideName[iside].data()));
        tdvnDiff_pt_side[iside]->cd();
        DiffEPv1_2SE_pt[iside]->Write();
        DiffEPv1_3SE_pt[iside]->Write();
        DiffEPv1_mix_pt[iside]->Write();
        DiffEPv1_mix_3_pt[iside]->Write();
        DiffEPv2_2SE_pt[iside]->Write();
        DiffEPv2_3SE_pt[iside]->Write();
        DiffEPv3_2SE_pt[iside]->Write();
        DiffEPv3_3SE_pt[iside]->Write();
        DiffSPv1_2SE_pt[iside]->Write();
        DiffSPv1_3SE_pt[iside]->Write();
        DiffSPv1_mix_pt[iside]->Write();
        DiffSPv1_mix_3_pt[iside]->Write();
        DiffSPv2_2SE_pt[iside]->Write();
        DiffSPv2_3SE_pt[iside]->Write();
        DiffSPv3_2SE_pt[iside]->Write();
        DiffSPv3_3SE_pt[iside]->Write();
        DiffEPv1obs_pt[iside]->Write();
        DiffEPv1obsMix_pt[iside]->Write();
        DiffEPv1obsMix_3_pt[iside]->Write();
        DiffEPv2obs_pt[iside]->Write();
        DiffEPv3obs_pt[iside]->Write();
        DiffSPv1num_pt[iside]->Write();
        DiffSPv1numMix_pt[iside]->Write();
        DiffSPv1numMix_3_pt[iside]->Write();
        DiffSPv2num_pt[iside]->Write();
        DiffSPv3num_pt[iside]->Write();
        Diffqcnt_pt[iside]->Write();
        Diffncnt_pt[iside]->Write();

        tdvnDiff_eta_side[iside] = (TDirectory *) tdvnDiff_eta->mkdir(Form("%s",SideName[iside].data()));
        tdvnDiff_eta_side[iside]->cd();
        DiffEPv1_2SE_eta[iside]->Write();
        DiffEPv1_3SE_eta[iside]->Write();
        DiffEPv1_mix_eta[iside]->Write();
        DiffEPv1_mix_3_eta[iside]->Write();
        DiffSPv1_2SE_eta[iside]->Write();
        DiffSPv1_3SE_eta[iside]->Write();
        DiffSPv1_mix_eta[iside]->Write();
        DiffSPv1_mix_3_eta[iside]->Write();
        DiffEPv1obs_eta[iside]->Write();
        DiffEPv1obsMix_eta[iside]->Write();
        DiffEPv1obsMix_3_eta[iside]->Write();
        DiffSPv1num_eta[iside]->Write();
        DiffSPv1numMix_eta[iside]->Write();
        DiffSPv1numMix_3_eta[iside]->Write();
        Diffqcnt_eta[iside]->Write();
        Diffncnt_eta[iside]->Write();
    }

    tdEPCorr->cd();
    corr_HFm1_HFp1->Write();
    corr_trackmid1_HFm1->Write();
    corr_trackmid1_HFp1->Write();
    corr_trackm1_HFm1->Write();
    corr_trackm1_trackp1->Write();
    corr_trackm1_HFp1->Write();
    corr_trackp1_HFm1->Write();
    corr_trackp1_trackm1->Write();
    corr_trackp1_HFp1->Write();
    corr_HFm1_trackm2->Write();
    corr_HFm1_trackp2->Write();
    corr_HFp1_trackm2->Write();
    corr_HFp1_trackp2->Write();
    corr_HFm1_HFp2->Write();
    corr_HFp1_HFm2->Write();
    corr_trackmid1_HFm2->Write();
    corr_trackmid1_HFp2->Write();
    corr_trackm1_HFm2->Write();
    corr_trackm1_trackp2->Write();
    corr_trackm1_HFp2->Write();
    corr_trackp1_HFm2->Write();
    corr_trackp1_trackm2->Write();
    corr_trackp1_HFp2->Write();
    corr_HFm2_HFp2->Write();
    corr_trackmid2_HFm2->Write();
    corr_trackmid2_HFp2->Write();
    corr_trackm2_HFm2->Write();
    corr_trackm2_trackp2->Write();
    corr_trackm2_HFp2->Write();
    corr_trackp2_HFm2->Write();
    corr_trackp2_HFp2->Write();

    tdEPCos->cd();
    cos_HFm1_HFp1->Write();
    cos_trackmid1_HFm1->Write();
    cos_trackmid1_HFp1->Write();
    cos_trackm1_HFm1->Write();
    cos_trackm1_trackp1->Write();
    cos_trackm1_HFp1->Write();
    cos_trackp1_HFm1->Write();
    cos_trackp1_trackm1->Write();
    cos_trackp1_HFp1->Write();
    cos_HFm1_trackm2->Write();
    cos_HFm1_trackp2->Write();
    cos_HFp1_trackm2->Write();
    cos_HFp1_trackp2->Write();
    cos_HFm1_HFp2->Write();
    cos_HFp1_HFm2->Write();
    cos_trackmid1_HFm2->Write();
    cos_trackmid1_HFp2->Write();
    cos_trackm1_HFm2->Write();
    cos_trackm1_trackp2->Write();
    cos_trackm1_HFp2->Write();
    cos_trackp1_HFm2->Write();
    cos_trackp1_trackm2->Write();
    cos_trackp1_HFp2->Write();
    cos_HFm2_HFp2->Write();
    cos_trackmid2_HFm2->Write();
    cos_trackmid2_HFp2->Write();
    cos_trackm2_HFm2->Write();
    cos_trackm2_trackp2->Write();
    cos_trackm2_HFp2->Write();
    cos_trackp2_HFm2->Write();
    cos_trackp2_HFp2->Write();
}


void ComputeVN( Int_t nevents, Int_t evtmult, Bool_t isodd, Double_t setv1, Double_t setv2, Double_t setv3, Bool_t eta_weights, Bool_t pt_weights, Bool_t conserve_pT, Bool_t addholes, Bool_t flatten, Bool_t recenter, Int_t iseed, TString mtag, Int_t ntries ) {

    double v1in = setv1;
    double v2in = setv2;
    double v3in = setv3;

    for (int iside = 0; iside<2; iside++) {
        for (int pbin = 0; pbin<nptbins; pbin++) {
            DiffEPv1obs_pt[iside]->SetBinContent(pbin+1, EPv1obs_pt[iside][pbin]->GetMean());
            DiffEPv1obs_pt[iside]->SetBinError(pbin+1, EPv1obs_pt[iside][pbin]->GetMeanError());
            DiffEPv1obsMix_pt[iside]->SetBinContent(pbin+1, EPv1obsMix_pt[iside][pbin]->GetMean());
            DiffEPv1obsMix_pt[iside]->SetBinError(pbin+1, EPv1obsMix_pt[iside][pbin]->GetMeanError());
            DiffEPv1obsMix_3_pt[iside]->SetBinContent(pbin+1, EPv1obsMix_3_pt[iside][pbin]->GetMean());
            DiffEPv1obsMix_3_pt[iside]->SetBinError(pbin+1, EPv1obsMix_3_pt[iside][pbin]->GetMeanError());
            DiffEPv2obs_pt[iside]->SetBinContent(pbin+1, EPv2obs_pt[iside][pbin]->GetMean());
            DiffEPv2obs_pt[iside]->SetBinError(pbin+1, EPv2obs_pt[iside][pbin]->GetMeanError());
            DiffEPv3obs_pt[iside]->SetBinContent(pbin+1, EPv3obs_pt[iside][pbin]->GetMean());
            DiffEPv3obs_pt[iside]->SetBinError(pbin+1, EPv3obs_pt[iside][pbin]->GetMeanError());
            DiffSPv1num_pt[iside]->SetBinContent(pbin+1, SPv1num_pt[iside][pbin]->GetMean());
            DiffSPv1num_pt[iside]->SetBinError(pbin+1, SPv1num_pt[iside][pbin]->GetMeanError());
            DiffSPv1numMix_pt[iside]->SetBinContent(pbin+1, SPv1numMix_pt[iside][pbin]->GetMean());
            DiffSPv1numMix_pt[iside]->SetBinError(pbin+1, SPv1numMix_pt[iside][pbin]->GetMeanError());
            DiffSPv1numMix_3_pt[iside]->SetBinContent(pbin+1, SPv1numMix_3_pt[iside][pbin]->GetMean());
            DiffSPv1numMix_3_pt[iside]->SetBinError(pbin+1, SPv1numMix_3_pt[iside][pbin]->GetMeanError());
            DiffSPv2num_pt[iside]->SetBinContent(pbin+1, SPv2num_pt[iside][pbin]->GetMean());
            DiffSPv2num_pt[iside]->SetBinError(pbin+1, SPv2num_pt[iside][pbin]->GetMeanError());
            DiffSPv3num_pt[iside]->SetBinContent(pbin+1, SPv3num_pt[iside][pbin]->GetMean());
            DiffSPv3num_pt[iside]->SetBinError(pbin+1, SPv3num_pt[iside][pbin]->GetMeanError());
            Diffqcnt_pt[iside]->SetBinContent(pbin+1, hqcnt_pt[iside][pbin]->GetMean());
            Diffqcnt_pt[iside]->SetBinError(pbin+1, hqcnt_pt[iside][pbin]->GetMeanError());
            Diffncnt_pt[iside]->SetBinContent(pbin+1, hncnt_pt[iside][pbin]->GetMean());
            Diffncnt_pt[iside]->SetBinError(pbin+1, hncnt_pt[iside][pbin]->GetMeanError());

            DiffEPv1_2SE_pt[iside]->SetBinContent(pbin+1, EPv1_2SE_pt[iside][pbin]->GetMean());
            DiffEPv1_2SE_pt[iside]->SetBinError(pbin+1, EPv1_2SE_pt[iside][pbin]->GetMeanError());
            DiffEPv1_3SE_pt[iside]->SetBinContent(pbin+1, EPv1_3SE_pt[iside][pbin]->GetMean());
            DiffEPv1_3SE_pt[iside]->SetBinError(pbin+1, EPv1_3SE_pt[iside][pbin]->GetMeanError());
            DiffEPv1_mix_pt[iside]->SetBinContent(pbin+1, EPv1_mix_pt[iside][pbin]->GetMean());
            DiffEPv1_mix_pt[iside]->SetBinError(pbin+1, EPv1_mix_pt[iside][pbin]->GetMeanError());
            DiffEPv1_mix_3_pt[iside]->SetBinContent(pbin+1, EPv1_mix_3_pt[iside][pbin]->GetMean());
            DiffEPv1_mix_3_pt[iside]->SetBinError(pbin+1, EPv1_mix_3_pt[iside][pbin]->GetMeanError());
            DiffEPv2_2SE_pt[iside]->SetBinContent(pbin+1, EPv2_2SE_pt[iside][pbin]->GetMean());
            DiffEPv2_2SE_pt[iside]->SetBinError(pbin+1, EPv2_2SE_pt[iside][pbin]->GetMeanError());
            DiffEPv2_3SE_pt[iside]->SetBinContent(pbin+1, EPv2_3SE_pt[iside][pbin]->GetMean());
            DiffEPv2_3SE_pt[iside]->SetBinError(pbin+1, EPv2_3SE_pt[iside][pbin]->GetMeanError());
            DiffEPv3_2SE_pt[iside]->SetBinContent(pbin+1, EPv3_2SE_pt[iside][pbin]->GetMean());
            DiffEPv3_2SE_pt[iside]->SetBinError(pbin+1, EPv3_2SE_pt[iside][pbin]->GetMeanError());
            DiffEPv3_3SE_pt[iside]->SetBinContent(pbin+1, EPv3_3SE_pt[iside][pbin]->GetMean());
            DiffEPv3_3SE_pt[iside]->SetBinError(pbin+1, EPv3_3SE_pt[iside][pbin]->GetMeanError());

            DiffSPv1_2SE_pt[iside]->SetBinContent(pbin+1, SPv1_2SE_pt[iside][pbin]->GetMean());
            DiffSPv1_2SE_pt[iside]->SetBinError(pbin+1, SPv1_2SE_pt[iside][pbin]->GetMeanError());
            DiffSPv1_3SE_pt[iside]->SetBinContent(pbin+1, SPv1_3SE_pt[iside][pbin]->GetMean());
            DiffSPv1_3SE_pt[iside]->SetBinError(pbin+1, SPv1_3SE_pt[iside][pbin]->GetMeanError());
            DiffSPv1_mix_pt[iside]->SetBinContent(pbin+1, SPv1_mix_pt[iside][pbin]->GetMean());
            DiffSPv1_mix_pt[iside]->SetBinError(pbin+1, SPv1_mix_pt[iside][pbin]->GetMeanError());
            DiffSPv1_mix_3_pt[iside]->SetBinContent(pbin+1, SPv1_mix_3_pt[iside][pbin]->GetMean());
            DiffSPv1_mix_3_pt[iside]->SetBinError(pbin+1, SPv1_mix_3_pt[iside][pbin]->GetMeanError());
            DiffSPv2_2SE_pt[iside]->SetBinContent(pbin+1, SPv2_2SE_pt[iside][pbin]->GetMean());
            DiffSPv2_2SE_pt[iside]->SetBinError(pbin+1, SPv2_2SE_pt[iside][pbin]->GetMeanError());
            DiffSPv2_3SE_pt[iside]->SetBinContent(pbin+1, SPv2_3SE_pt[iside][pbin]->GetMean());
            DiffSPv2_3SE_pt[iside]->SetBinError(pbin+1, SPv2_3SE_pt[iside][pbin]->GetMeanError());
            DiffSPv3_2SE_pt[iside]->SetBinContent(pbin+1, SPv3_2SE_pt[iside][pbin]->GetMean());
            DiffSPv3_2SE_pt[iside]->SetBinError(pbin+1, SPv3_2SE_pt[iside][pbin]->GetMeanError());
            DiffSPv3_3SE_pt[iside]->SetBinContent(pbin+1, SPv3_3SE_pt[iside][pbin]->GetMean());
            DiffSPv3_3SE_pt[iside]->SetBinError(pbin+1, SPv3_3SE_pt[iside][pbin]->GetMeanError());
        }


        for (int ebin = 0; ebin<netabins; ebin++) {
            hQ1nA_final[iside]->SetBinContent(ebin+1, hQ1nA[iside][ebin]->GetMean());
            hQ1nA_final[iside]->SetBinError(ebin+1, hQ1nA[iside][ebin]->GetMeanError());
            hQmixnAC_final[iside]->SetBinContent(ebin+1, hQmixnAC[iside][ebin]->GetMean());
            hQmixnAC_final[iside]->SetBinError(ebin+1, hQmixnAC[iside][ebin]->GetMeanError());
            hQmixnAC_3_final[iside]->SetBinContent(ebin+1, hQmixnAC_3[iside][ebin]->GetMean());
            hQmixnAC_3_final[iside]->SetBinError(ebin+1, hQmixnAC_3[iside][ebin]->GetMeanError());
            hQ1nAnorm_final[iside]->SetBinContent(ebin+1, hQ1nAnorm[iside][ebin]->GetMean());
            hQ1nAnorm_final[iside]->SetBinError(ebin+1, hQ1nAnorm[iside][ebin]->GetMeanError());
            hQmixnACnorm_final[iside]->SetBinContent(ebin+1, hQmixnACnorm[iside][ebin]->GetMean());
            hQmixnACnorm_final[iside]->SetBinError(ebin+1, hQmixnACnorm[iside][ebin]->GetMeanError());
            hQmixnAC_3norm_final[iside]->SetBinContent(ebin+1, hQmixnAC_3norm[iside][ebin]->GetMean());
            hQmixnAC_3norm_final[iside]->SetBinError(ebin+1, hQmixnAC_3norm[iside][ebin]->GetMeanError());

            DiffEPv1obs_eta[iside]->SetBinContent(ebin+1, EPv1obs_eta[iside][ebin]->GetMean());
            DiffEPv1obs_eta[iside]->SetBinError(ebin+1, EPv1obs_eta[iside][ebin]->GetMeanError());
            DiffEPv1obsMix_eta[iside]->SetBinContent(ebin+1, EPv1obsMix_eta[iside][ebin]->GetMean());
            DiffEPv1obsMix_eta[iside]->SetBinError(ebin+1, EPv1obsMix_eta[iside][ebin]->GetMeanError());
            DiffEPv1obsMix_3_eta[iside]->SetBinContent(ebin+1, EPv1obsMix_3_eta[iside][ebin]->GetMean());
            DiffEPv1obsMix_3_eta[iside]->SetBinError(ebin+1, EPv1obsMix_3_eta[iside][ebin]->GetMeanError());
            DiffSPv1num_eta[iside]->SetBinContent(ebin+1, SPv1num_eta[iside][ebin]->GetMean());
            DiffSPv1num_eta[iside]->SetBinError(ebin+1, SPv1num_eta[iside][ebin]->GetMeanError());
            DiffSPv1numMix_eta[iside]->SetBinContent(ebin+1, SPv1numMix_eta[iside][ebin]->GetMean());
            DiffSPv1numMix_eta[iside]->SetBinError(ebin+1, SPv1numMix_eta[iside][ebin]->GetMeanError());
            DiffSPv1numMix_3_eta[iside]->SetBinContent(ebin+1, SPv1numMix_3_eta[iside][ebin]->GetMean());
            DiffSPv1numMix_3_eta[iside]->SetBinError(ebin+1, SPv1numMix_3_eta[iside][ebin]->GetMeanError());
            Diffqcnt_eta[iside]->SetBinContent(ebin+1, hqcnt_eta[iside][ebin]->GetMean());
            Diffqcnt_eta[iside]->SetBinError(ebin+1, hqcnt_eta[iside][ebin]->GetMeanError());
            Diffncnt_eta[iside]->SetBinContent(ebin+1, hncnt_eta[iside][ebin]->GetMean());
            Diffncnt_eta[iside]->SetBinError(ebin+1, hncnt_eta[iside][ebin]->GetMeanError());

            DiffEPv1_2SE_eta[iside]->SetBinContent(ebin+1, EPv1_2SE_eta[iside][ebin]->GetMean());
            DiffEPv1_2SE_eta[iside]->SetBinError(ebin+1, EPv1_2SE_eta[iside][ebin]->GetMeanError());
            DiffEPv1_3SE_eta[iside]->SetBinContent(ebin+1, EPv1_3SE_eta[iside][ebin]->GetMean());
            DiffEPv1_3SE_eta[iside]->SetBinError(ebin+1, EPv1_3SE_eta[iside][ebin]->GetMeanError());
            DiffEPv1_mix_eta[iside]->SetBinContent(ebin+1, EPv1_mix_eta[iside][ebin]->GetMean());
            DiffEPv1_mix_eta[iside]->SetBinError(ebin+1, EPv1_mix_eta[iside][ebin]->GetMeanError());
            DiffEPv1_mix_3_eta[iside]->SetBinContent(ebin+1, EPv1_mix_3_eta[iside][ebin]->GetMean());
            DiffEPv1_mix_3_eta[iside]->SetBinError(ebin+1, EPv1_mix_3_eta[iside][ebin]->GetMeanError());

            DiffSPv1_2SE_eta[iside]->SetBinContent(ebin+1, SPv1_2SE_eta[iside][ebin]->GetMean());
            DiffSPv1_2SE_eta[iside]->SetBinError(ebin+1, SPv1_2SE_eta[iside][ebin]->GetMeanError());
            DiffSPv1_3SE_eta[iside]->SetBinContent(ebin+1, SPv1_3SE_eta[iside][ebin]->GetMean());
            DiffSPv1_3SE_eta[iside]->SetBinError(ebin+1, SPv1_3SE_eta[iside][ebin]->GetMeanError());
            DiffSPv1_mix_eta[iside]->SetBinContent(ebin+1, SPv1_mix_eta[iside][ebin]->GetMean());
            DiffSPv1_mix_eta[iside]->SetBinError(ebin+1, SPv1_mix_eta[iside][ebin]->GetMeanError());
            DiffSPv1_mix_3_eta[iside]->SetBinContent(ebin+1, SPv1_mix_3_eta[iside][ebin]->GetMean());
            DiffSPv1_mix_3_eta[iside]->SetBinError(ebin+1, SPv1_mix_3_eta[iside][ebin]->GetMeanError());
        }
    }

    ofstream fout;
    if (!fopen("logs","r")) system(Form("mkdir logs"));
    TString tag = Form("data_final_tally_%s",mtag.Data());
    TString foutname = "logs/"+tag+".dat";
    fout.open(foutname.Data());
    fout << "====================" << endl;
    fout << "nevents:    " << nevents*ntries << endl;
    fout << "event multiplicity: " << evtmult << endl;
    if (isodd) fout << "Rapidity-odd v1 input " << endl;
    else fout << "Rapidity-even v1:  " << setv1 << endl;
    fout << "Input v2:     " << setv2 << endl;
    fout << "Input v3:     " << setv3 << endl;
    if (eta_weights) fout << "Using eta-dependent weights..." <<endl;
    if (pt_weights) fout << "Using pt-dependent weights..." <<endl;
    if (conserve_pT) fout << "Conserving transverse momentum event by event... " << endl;
    if (addholes) fout << "Holes in detector acceptance... " << endl;
    if (recenter) fout << "Recentered... " << endl;
    if (flatten) fout << "Flattened... " << endl;
    fout << Form("tracker subevent %0.1f < |eta| < %0.1f",ecutmin[2],ecutmax[2]) << endl;
    fout << "iseed:    " << iseed << endl;
    fout << "   -------   " << endl;


    fout<<"\n  ---- Final Tally (negative eta) ---- \n"<<endl;
    fout<<"rescor1_2SE:     "<<Form("%.5f",rescor1_2SE->GetMean())<<" +/- "<<Form("%.5f",rescor1_2SE->GetMeanError())<<endl;
    fout<<"rescor1_3SE:     "<<Form("%.5f",rescor1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_3SE[0]->GetMeanError())<<endl;
    fout<<"rescor2HF_2SE:   "<<Form("%.5f",rescor2HF_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor2HF_2SE->GetMeanError())<<endl;
    fout<<"rescor2HF_3SE:   "<<Form("%.5f",rescor2HF_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor2HF_3SE[0]->GetMeanError())<<endl;
    fout<<"rescor2Trk_3SE:  "<<Form("%.5f",rescor2Trk_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor2Trk_3SE[0]->GetMeanError())<<endl;
    fout<<"rescor3HF_2SE:   "<<Form("%.5f",rescor3HF_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor3HF_2SE->GetMeanError())<<endl;
    fout<<"rescor3HF_3SE:   "<<Form("%.5f",rescor3HF_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor3HF_3SE[0]->GetMeanError())<<endl;
    fout<<"rescor3Trk_3SE:  "<<Form("%.5f",rescor3Trk_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor3Trk_3SE[0]->GetMeanError())<<endl;
    fout<<"rescor1_mix:     "<<Form("%.5f",rescor1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix[0]->GetMeanError())<<endl;
    fout<<"rescor1_mix3:    "<<Form("%.5f",rescor1_mix3[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix3[0]->GetMeanError())<<endl;
    fout<<"EPv1obs:         "<<Form("%.5f",EPv1obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obs[0]->GetMeanError())<<endl;
    fout<<"EPv2obs:         "<<Form("%.5f",EPv2obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2obs[0]->GetMeanError())<<endl;
    fout<<"EPv3obs:         "<<Form("%.5f",EPv3obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv3obs[0]->GetMeanError())<<endl;
    fout<<"EPv1obs_mix:     "<<Form("%.5f",EPv1obsMix[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix[0]->GetMeanError())<<endl;
    fout<<"EPv1obsMix_3:     "<<Form("%.5f",EPv1obsMix_3[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix_3[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"SPdenom1_2SE:    "<<Form("%.5f",SPdenom1_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_2SE->GetMeanError())<<endl;
    fout<<"SPdenom1_3SE:    "<<Form("%.5f",SPdenom1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_3SE[0]->GetMeanError())<<endl;
    fout<<"SPdenom2HF_2SE:  "<<Form("%.5f",SPdenom2HF_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom2HF_2SE->GetMeanError())<<endl;
    fout<<"SPdenom2HF_3SE:  "<<Form("%.5f",SPdenom2HF_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2HF_3SE[0]->GetMeanError())<<endl;
    fout<<"SPdenom2Trk_3SE: "<<Form("%.5f",SPdenom2Trk_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2Trk_3SE[0]->GetMeanError())<<endl;
    fout<<"SPdenom3HF_2SE:  "<<Form("%.5f",SPdenom3HF_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom3HF_2SE->GetMeanError())<<endl;
    fout<<"SPdenom3HF_3SE:  "<<Form("%.5f",SPdenom3HF_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom3HF_3SE[0]->GetMeanError())<<endl;
    fout<<"SPdenom3Trk_3SE: "<<Form("%.5f",SPdenom3Trk_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom3Trk_3SE[0]->GetMeanError())<<endl;
    fout<<"SPdenom1_mix:    "<<Form("%.5f",SPdenom1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix[0]->GetMeanError())<<endl;
    fout<<"SPdenom1_mix_3:  "<<Form("%.5f",SPdenom1_mix_3[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix_3[0]->GetMeanError())<<endl;
    fout<<"SPv1num:         "<<Form("%.5f",SPv1num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1num[0]->GetMeanError())<<endl;
    fout<<"SPv2num:         "<<Form("%.5f",SPv2num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2num[0]->GetMeanError())<<endl;
    fout<<"SPv3num:         "<<Form("%.5f",SPv3num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv3num[0]->GetMeanError())<<endl;
    fout<<"SPv1numMix:      "<<Form("%.5f",SPv1numMix[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix[0]->GetMeanError())<<endl;
    fout<<"SPv1numMix_3:    "<<Form("%.5f",SPv1numMix_3[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix_3[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"Q1Q1Q2:          "<<Form("%.5f",hQ1Q1Q2[0]->GetMean())<<" +/- "<<Form("%.5f",hQ1Q1Q2[0]->GetMeanError())<<endl;
    fout<<"Q1Q1Q2norm:      "<<Form("%.5f",hQ1Q1Q2norm[0]->GetMean())<<" +/- "<<Form("%.5f",hQ1Q1Q2norm[0]->GetMeanError())<<endl;
    fout<<"q112:            "<<Form("%.5f",hq112->GetMean())<<" +/- "<<Form("%.5f",hq112->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"EPv1_2SE:        "<<Form("%.5f",EPv1_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_2SE[0]->GetMeanError())<<endl;
    fout<<"EPv1_3SE:        "<<Form("%.5f",EPv1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_3SE[0]->GetMeanError())<<endl;
    fout<<"EPv1_mix:        "<<Form("%.5f",EPv1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix[0]->GetMeanError())<<endl;
    fout<<"EPv1_mix_3:      "<<Form("%.5f",EPv1_mix_3[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix_3[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"SPv1_2SE:        "<<Form("%.5f",SPv1_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_2SE[0]->GetMeanError())<<endl;
    fout<<"SPv1_3SE:        "<<Form("%.5f",SPv1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_3SE[0]->GetMeanError())<<endl;
    fout<<"SPv1_mix:        "<<Form("%.5f",SPv1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix[0]->GetMeanError())<<endl;
    fout<<"SPv1_mix_3:      "<<Form("%.5f",SPv1_mix_3[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix_3[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"EPv2_2SE:        "<<Form("%.5f",EPv2_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2_2SE[0]->GetMeanError())<<endl;
    fout<<"EPv2_3SE:        "<<Form("%.5f",EPv2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2_3SE[0]->GetMeanError())<<endl;
    fout<<"SPv2_2SE:        "<<Form("%.5f",SPv2_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2_2SE[0]->GetMeanError())<<endl;
    fout<<"SPv2_3SE:        "<<Form("%.5f",SPv2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2_3SE[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"EPv3_2SE:        "<<Form("%.5f",EPv3_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv3_2SE[0]->GetMeanError())<<endl;
    fout<<"EPv3_3SE:        "<<Form("%.5f",EPv3_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv3_3SE[0]->GetMeanError())<<endl;
    fout<<"SPv3_2SE:        "<<Form("%.5f",SPv3_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv3_2SE[0]->GetMeanError())<<endl;
    fout<<"SPv3_3SE:        "<<Form("%.5f",SPv3_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv3_3SE[0]->GetMeanError())<<endl;


    fout<<"\n  ---- Final Tally (positive eta) ---- \n"<<endl;
    fout<<"rescor1_2SE:     "<<Form("%.5f",rescor1_2SE->GetMean())<<" +/- "<<Form("%.5f",rescor1_2SE->GetMeanError())<<endl;
    fout<<"rescor1_3SE:     "<<Form("%.5f",rescor1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_3SE[1]->GetMeanError())<<endl;
    fout<<"rescor2HF_2SE:   "<<Form("%.5f",rescor2HF_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor2HF_2SE->GetMeanError())<<endl;
    fout<<"rescor2HF_3SE:   "<<Form("%.5f",rescor2HF_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor2HF_3SE[1]->GetMeanError())<<endl;
    fout<<"rescor2Trk_3SE:  "<<Form("%.5f",rescor2Trk_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor2Trk_3SE[1]->GetMeanError())<<endl;
    fout<<"rescor3HF_2SE:   "<<Form("%.5f",rescor3HF_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor3HF_2SE->GetMeanError())<<endl;
    fout<<"rescor3HF_3SE:   "<<Form("%.5f",rescor3HF_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor3HF_3SE[1]->GetMeanError())<<endl;
    fout<<"rescor3Trk_3SE:  "<<Form("%.5f",rescor3Trk_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor3Trk_3SE[1]->GetMeanError())<<endl;
    fout<<"rescor1_mix:     "<<Form("%.5f",rescor1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix[1]->GetMeanError())<<endl;
    fout<<"rescor1_mix3:    "<<Form("%.5f",rescor1_mix3[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix3[1]->GetMeanError())<<endl;
    fout<<"EPv1obs:         "<<Form("%.5f",EPv1obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obs[1]->GetMeanError())<<endl;
    fout<<"EPv2obs:         "<<Form("%.5f",EPv2obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2obs[1]->GetMeanError())<<endl;
    fout<<"EPv3obs:         "<<Form("%.5f",EPv3obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv3obs[1]->GetMeanError())<<endl;
    fout<<"EPv1obsMix:      "<<Form("%.5f",EPv1obsMix[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix[1]->GetMeanError())<<endl;
    fout<<"EPv1obsMix_3:    "<<Form("%.5f",EPv1obsMix_3[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix_3[1]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"SPdenom1_2SE:    "<<Form("%.5f",SPdenom1_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_2SE->GetMeanError())<<endl;
    fout<<"SPdenom1_3SE:    "<<Form("%.5f",SPdenom1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_3SE[1]->GetMeanError())<<endl;
    fout<<"SPdenom2HF_2SE:  "<<Form("%.5f",SPdenom2HF_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom2HF_2SE->GetMeanError())<<endl;
    fout<<"SPdenom2HF_3SE:  "<<Form("%.5f",SPdenom2HF_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2HF_3SE[1]->GetMeanError())<<endl;
    fout<<"SPdenom2Trk_3SE: "<<Form("%.5f",SPdenom2Trk_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2Trk_3SE[1]->GetMeanError())<<endl;
    fout<<"SPdenom3HF_2SE:  "<<Form("%.5f",SPdenom3HF_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom3HF_2SE->GetMeanError())<<endl;
    fout<<"SPdenom3HF_3SE:  "<<Form("%.5f",SPdenom3HF_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom3HF_3SE[1]->GetMeanError())<<endl;
    fout<<"SPdenom3Trk_3SE: "<<Form("%.5f",SPdenom3Trk_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom3Trk_3SE[1]->GetMeanError())<<endl;
    fout<<"SPdenom1_mix:    "<<Form("%.5f",SPdenom1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix[1]->GetMeanError())<<endl;
    fout<<"SPdenom1_mix_3:  "<<Form("%.5f",SPdenom1_mix_3[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix_3[1]->GetMeanError())<<endl;
    fout<<"SPv1num:         "<<Form("%.5f",SPv1num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1num[1]->GetMeanError())<<endl;
    fout<<"SPv2num:         "<<Form("%.5f",SPv2num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2num[1]->GetMeanError())<<endl;
    fout<<"SPv3num:         "<<Form("%.5f",SPv3num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv3num[1]->GetMeanError())<<endl;
    fout<<"SPv1numMix:      "<<Form("%.5f",SPv1numMix[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix[1]->GetMeanError())<<endl;
    fout<<"SPv1numMix_3:    "<<Form("%.5f",SPv1numMix_3[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix_3[1]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"Q1Q1Q2:          "<<Form("%.5f",hQ1Q1Q2[1]->GetMean())<<" +/- "<<Form("%.5f",hQ1Q1Q2[1]->GetMeanError())<<endl;
    fout<<"Q1Q1Q2norm:      "<<Form("%.5f",hQ1Q1Q2norm[1]->GetMean())<<" +/- "<<Form("%.5f",hQ1Q1Q2norm[1]->GetMeanError())<<endl;
    fout<<"q112:            "<<Form("%.5f",hq112->GetMean())<<" +/- "<<Form("%.5f",hq112->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"EPv1_2SE:        "<<Form("%.5f",EPv1_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_2SE[1]->GetMeanError())<<endl;
    fout<<"EPv1_3SE:        "<<Form("%.5f",EPv1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_3SE[1]->GetMeanError())<<endl;
    fout<<"EPv1_mix:        "<<Form("%.5f",EPv1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix[1]->GetMeanError())<<endl;
    fout<<"EPv1_mix_3:      "<<Form("%.5f",EPv1_mix_3[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix_3[1]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"SPv1_2SE:        "<<Form("%.5f",SPv1_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_2SE[1]->GetMeanError())<<endl;
    fout<<"SPv1_3SE:        "<<Form("%.5f",SPv1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_3SE[1]->GetMeanError())<<endl;
    fout<<"SPv1_mix:        "<<Form("%.5f",SPv1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix[1]->GetMeanError())<<endl;
    fout<<"SPv1_mix_3:      "<<Form("%.5f",SPv1_mix_3[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix_3[1]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"EPv2_2SE:        "<<Form("%.5f",EPv2_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2_2SE[1]->GetMeanError())<<endl;
    fout<<"EPv2_3SE:        "<<Form("%.5f",EPv2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2_3SE[1]->GetMeanError())<<endl;
    fout<<"SPv2_2SE:        "<<Form("%.5f",SPv2_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2_2SE[1]->GetMeanError())<<endl;
    fout<<"SPv2_3SE:        "<<Form("%.5f",SPv2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2_3SE[1]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"EPv3_2SE:        "<<Form("%.5f",EPv3_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv3_2SE[1]->GetMeanError())<<endl;
    fout<<"EPv3_3SE:        "<<Form("%.5f",EPv3_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv3_3SE[1]->GetMeanError())<<endl;
    fout<<"SPv3_2SE:        "<<Form("%.5f",SPv3_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv3_2SE[1]->GetMeanError())<<endl;
    fout<<"SPv3_3SE:        "<<Form("%.5f",SPv3_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv3_3SE[1]->GetMeanError())<<endl;
    fout<<"\n --- Comparisons to input vn --- \n"<<endl;
    if (!isodd) {
        fout<<Form("EPv1_2SE/v1in  (HFm): %0.4f",EPv1_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_2SE[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("EPv1_3SE/v1in  (HFm): %0.4f",EPv1_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_3SE[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("EPv1_mix/v1in  (HFm): %0.4f",EPv1_mix[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_mix[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("EPv1_mix_3/v1in  (HFm): %0.4f",EPv1_mix_3[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix_3[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_mix_3[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix_3[1]->GetMeanError()/v1in)<<endl;
        cout<<"   -----        "<<endl;
        fout<<Form("SPv1_2SE/v1in  (HFm): %0.4f",SPv1_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_2SE[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("SPv1_3SE/v1in  (HFm): %0.4f",SPv1_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_3SE[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("SPv1_mix/v1in  (HFm): %0.4f",SPv1_mix[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_mix[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("SPv1_mix_3/v1in  (HFm): %0.4f",SPv1_mix_3[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix_3[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_mix_3[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix_3[1]->GetMeanError()/v1in)<<endl;
        fout<<"   -----        "<<endl;
    }
    fout<<Form("EPv2_2SE/v2in  (HFm): %0.4f",EPv2_2SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_2SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",EPv2_2SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_2SE[1]->GetMeanError()/v2in)<<endl;
    fout<<Form("EPv2_3SE/v2in  (HFm): %0.4f",EPv2_3SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_3SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",EPv2_3SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_3SE[1]->GetMeanError()/v2in)<<endl;
    fout<<Form("SPv2_2SE/v2in  (HFm): %0.4f",SPv2_2SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_2SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",SPv2_2SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_2SE[1]->GetMeanError()/v2in)<<endl;
    fout<<Form("SPv2_3SE/v2in  (HFm): %0.4f",SPv2_3SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_3SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",SPv2_3SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_3SE[1]->GetMeanError()/v2in)<<endl;
    fout<<"   -----        "<<endl;
    fout<<Form("EPv3_2SE/v3in  (HFm): %0.4f",EPv3_2SE[0]->GetMean()/v3in)<<Form(" +/- %0.4f",EPv3_2SE[0]->GetMeanError()/v3in)<<Form("\t (HFp): %0.4f",EPv3_2SE[1]->GetMean()/v3in)<<Form(" +/- %0.4f",EPv3_2SE[1]->GetMeanError()/v3in)<<endl;
    fout<<Form("EPv3_3SE/v3in  (HFm): %0.4f",EPv3_3SE[0]->GetMean()/v3in)<<Form(" +/- %0.4f",EPv3_3SE[0]->GetMeanError()/v3in)<<Form("\t (HFp): %0.4f",EPv3_3SE[1]->GetMean()/v3in)<<Form(" +/- %0.4f",EPv3_3SE[1]->GetMeanError()/v3in)<<endl;
    fout<<Form("SPv3_2SE/v3in  (HFm): %0.4f",SPv3_2SE[0]->GetMean()/v3in)<<Form(" +/- %0.4f",SPv3_2SE[0]->GetMeanError()/v3in)<<Form("\t (HFp): %0.4f",SPv3_2SE[1]->GetMean()/v3in)<<Form(" +/- %0.4f",SPv3_2SE[1]->GetMeanError()/v3in)<<endl;
    fout<<Form("SPv3_3SE/v3in  (HFm): %0.4f",SPv3_3SE[0]->GetMean()/v3in)<<Form(" +/- %0.4f",SPv3_3SE[0]->GetMeanError()/v3in)<<Form("\t (HFp): %0.4f",SPv3_3SE[1]->GetMean()/v3in)<<Form(" +/- %0.4f",SPv3_3SE[1]->GetMeanError()/v3in)<<endl;
    fout<<"   -----        "<<endl;
    double q1q1q2ave = 0.5*(hQ1Q1Q2norm[0]->GetMean() + hQ1Q1Q2norm[1]->GetMean());
    double q1q1q2aveErr = 0.5*sqrt( pow(hQ1Q1Q2norm[0]->GetMeanError(),2) + pow(hQ1Q1Q2norm[1]->GetMeanError(),2) );
    if (!isodd) cout<<"v1: "<<setv1<<"\tv2: "<<setv2<<"\tQ1Q1Q2norm:      "<<Form("%.5f",q1q1q2ave)<<" +/- "<<Form("%.5f",q1q1q2aveErr)<<endl;
    fout<<"\n ...done \n"<<endl;
    fout.close();



    cout << "====================" << endl;
    cout << "nevents:    " << nevents*ntries << endl;
    cout << "event multiplicity: " << evtmult << endl;
    if (isodd) cout << "Rapidity-odd v1 input " << endl;
    else cout << "Rapidity-even v1:  " << setv1 << endl;
    cout << "Input v2:     " << setv2 << endl;
    cout << "Input v3:     " << setv3 << endl;
    if (eta_weights) cout << "Using eta-dependent weights..." <<endl;
    if (pt_weights) cout << "Using pt-dependent weights..." <<endl;
    if (conserve_pT) cout << "Conserving transverse momentum event by event... " << endl;
    if (addholes) cout << "Holes in detector acceptance... " << endl;
    if (recenter) cout << "Recentered... " << endl;
    if (flatten) cout << "Flattened... " << endl;
    cout << Form("tracker subevent %0.1f < |eta| < %0.1f",ecutmin[2],ecutmax[2]) << endl;
    cout << "iseed:    " << iseed << endl;
    cout << "   -------   " << endl;


    cout<<"\n  ---- Final Tally (negative eta) ---- \n"<<endl;
    cout<<"rescor1_2SE:     "<<Form("%.5f",rescor1_2SE->GetMean())<<" +/- "<<Form("%.5f",rescor1_2SE->GetMeanError())<<endl;
    cout<<"rescor1_3SE:     "<<Form("%.5f",rescor1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_3SE[0]->GetMeanError())<<endl;
    cout<<"rescor2HF_2SE:   "<<Form("%.5f",rescor2HF_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor2HF_2SE->GetMeanError())<<endl;
    cout<<"rescor2HF_3SE:   "<<Form("%.5f",rescor2HF_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor2HF_3SE[0]->GetMeanError())<<endl;
    cout<<"rescor2Trk_3SE:  "<<Form("%.5f",rescor2Trk_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor2Trk_3SE[0]->GetMeanError())<<endl;
    cout<<"rescor3HF_2SE:   "<<Form("%.5f",rescor3HF_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor3HF_2SE->GetMeanError())<<endl;
    cout<<"rescor3HF_3SE:   "<<Form("%.5f",rescor3HF_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor3HF_3SE[0]->GetMeanError())<<endl;
    cout<<"rescor3Trk_3SE:  "<<Form("%.5f",rescor3Trk_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor3Trk_3SE[0]->GetMeanError())<<endl;
    cout<<"rescor1_mix:     "<<Form("%.5f",rescor1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix[0]->GetMeanError())<<endl;
    cout<<"rescor1_mix3:    "<<Form("%.5f",rescor1_mix3[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix3[0]->GetMeanError())<<endl;
    cout<<"EPv1obs:         "<<Form("%.5f",EPv1obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obs[0]->GetMeanError())<<endl;
    cout<<"EPv2obs:         "<<Form("%.5f",EPv2obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2obs[0]->GetMeanError())<<endl;
    cout<<"EPv3obs:         "<<Form("%.5f",EPv3obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv3obs[0]->GetMeanError())<<endl;
    cout<<"EPv1obs_mix:     "<<Form("%.5f",EPv1obsMix[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix[0]->GetMeanError())<<endl;
    cout<<"EPv1obsMix_3:     "<<Form("%.5f",EPv1obsMix_3[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix_3[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"SPdenom1_2SE:    "<<Form("%.5f",SPdenom1_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_2SE->GetMeanError())<<endl;
    cout<<"SPdenom1_3SE:    "<<Form("%.5f",SPdenom1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_3SE[0]->GetMeanError())<<endl;
    cout<<"SPdenom2HF_2SE:  "<<Form("%.5f",SPdenom2HF_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom2HF_2SE->GetMeanError())<<endl;
    cout<<"SPdenom2HF_3SE:  "<<Form("%.5f",SPdenom2HF_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2HF_3SE[0]->GetMeanError())<<endl;
    cout<<"SPdenom2Trk_3SE: "<<Form("%.5f",SPdenom2Trk_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2Trk_3SE[0]->GetMeanError())<<endl;
    cout<<"SPdenom3HF_2SE:  "<<Form("%.5f",SPdenom3HF_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom3HF_2SE->GetMeanError())<<endl;
    cout<<"SPdenom3HF_3SE:  "<<Form("%.5f",SPdenom3HF_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom3HF_3SE[0]->GetMeanError())<<endl;
    cout<<"SPdenom3Trk_3SE: "<<Form("%.5f",SPdenom3Trk_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom3Trk_3SE[0]->GetMeanError())<<endl;
    cout<<"SPdenom1_mix:    "<<Form("%.5f",SPdenom1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix[0]->GetMeanError())<<endl;
    cout<<"SPdenom1_mix_3:  "<<Form("%.5f",SPdenom1_mix_3[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix_3[0]->GetMeanError())<<endl;
    cout<<"SPv1num:         "<<Form("%.5f",SPv1num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1num[0]->GetMeanError())<<endl;
    cout<<"SPv2num:         "<<Form("%.5f",SPv2num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2num[0]->GetMeanError())<<endl;
    cout<<"SPv3num:         "<<Form("%.5f",SPv3num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv3num[0]->GetMeanError())<<endl;
    cout<<"SPv1numMix:      "<<Form("%.5f",SPv1numMix[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix[0]->GetMeanError())<<endl;
    cout<<"SPv1numMix_3:    "<<Form("%.5f",SPv1numMix_3[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix_3[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"Q1Q1Q2:          "<<Form("%.5f",hQ1Q1Q2[0]->GetMean())<<" +/- "<<Form("%.5f",hQ1Q1Q2[0]->GetMeanError())<<endl;
    cout<<"Q1Q1Q2norm:      "<<Form("%.5f",hQ1Q1Q2norm[0]->GetMean())<<" +/- "<<Form("%.5f",hQ1Q1Q2norm[0]->GetMeanError())<<endl;
    cout<<"q112:            "<<Form("%.5f",hq112->GetMean())<<" +/- "<<Form("%.5f",hq112->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"EPv1_2SE:        "<<Form("%.5f",EPv1_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_2SE[0]->GetMeanError())<<endl;
    cout<<"EPv1_3SE:        "<<Form("%.5f",EPv1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_3SE[0]->GetMeanError())<<endl;
    cout<<"EPv1_mix:        "<<Form("%.5f",EPv1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix[0]->GetMeanError())<<endl;
    cout<<"EPv1_mix_3:      "<<Form("%.5f",EPv1_mix_3[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix_3[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"SPv1_2SE:        "<<Form("%.5f",SPv1_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_2SE[0]->GetMeanError())<<endl;
    cout<<"SPv1_3SE:        "<<Form("%.5f",SPv1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_3SE[0]->GetMeanError())<<endl;
    cout<<"SPv1_mix:        "<<Form("%.5f",SPv1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix[0]->GetMeanError())<<endl;
    cout<<"SPv1_mix_3:      "<<Form("%.5f",SPv1_mix_3[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix_3[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"EPv2_2SE:        "<<Form("%.5f",EPv2_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2_2SE[0]->GetMeanError())<<endl;
    cout<<"EPv2_3SE:        "<<Form("%.5f",EPv2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2_3SE[0]->GetMeanError())<<endl;
    cout<<"SPv2_2SE:        "<<Form("%.5f",SPv2_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2_2SE[0]->GetMeanError())<<endl;
    cout<<"SPv2_3SE:        "<<Form("%.5f",SPv2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2_3SE[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"EPv3_2SE:        "<<Form("%.5f",EPv3_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv3_2SE[0]->GetMeanError())<<endl;
    cout<<"EPv3_3SE:        "<<Form("%.5f",EPv3_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv3_3SE[0]->GetMeanError())<<endl;
    cout<<"SPv3_2SE:        "<<Form("%.5f",SPv3_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv3_2SE[0]->GetMeanError())<<endl;
    cout<<"SPv3_3SE:        "<<Form("%.5f",SPv3_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv3_3SE[0]->GetMeanError())<<endl;


    cout<<"\n  ---- Final Tally (positive eta) ---- \n"<<endl;
    cout<<"rescor1_2SE:     "<<Form("%.5f",rescor1_2SE->GetMean())<<" +/- "<<Form("%.5f",rescor1_2SE->GetMeanError())<<endl;
    cout<<"rescor1_3SE:     "<<Form("%.5f",rescor1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_3SE[1]->GetMeanError())<<endl;
    cout<<"rescor2HF_2SE:   "<<Form("%.5f",rescor2HF_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor2HF_2SE->GetMeanError())<<endl;
    cout<<"rescor2HF_3SE:   "<<Form("%.5f",rescor2HF_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor2HF_3SE[1]->GetMeanError())<<endl;
    cout<<"rescor2Trk_3SE:  "<<Form("%.5f",rescor2Trk_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor2Trk_3SE[1]->GetMeanError())<<endl;
    cout<<"rescor3HF_2SE:   "<<Form("%.5f",rescor3HF_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor3HF_2SE->GetMeanError())<<endl;
    cout<<"rescor3HF_3SE:   "<<Form("%.5f",rescor3HF_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor3HF_3SE[1]->GetMeanError())<<endl;
    cout<<"rescor3Trk_3SE:  "<<Form("%.5f",rescor3Trk_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor3Trk_3SE[1]->GetMeanError())<<endl;
    cout<<"rescor1_mix:     "<<Form("%.5f",rescor1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix[1]->GetMeanError())<<endl;
    cout<<"rescor1_mix3:    "<<Form("%.5f",rescor1_mix3[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix3[1]->GetMeanError())<<endl;
    cout<<"EPv1obs:         "<<Form("%.5f",EPv1obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obs[1]->GetMeanError())<<endl;
    cout<<"EPv2obs:         "<<Form("%.5f",EPv2obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2obs[1]->GetMeanError())<<endl;
    cout<<"EPv3obs:         "<<Form("%.5f",EPv3obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv3obs[1]->GetMeanError())<<endl;
    cout<<"EPv1obsMix:      "<<Form("%.5f",EPv1obsMix[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix[1]->GetMeanError())<<endl;
    cout<<"EPv1obsMix_3:    "<<Form("%.5f",EPv1obsMix_3[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix_3[1]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"SPdenom1_2SE:    "<<Form("%.5f",SPdenom1_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_2SE->GetMeanError())<<endl;
    cout<<"SPdenom1_3SE:    "<<Form("%.5f",SPdenom1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_3SE[1]->GetMeanError())<<endl;
    cout<<"SPdenom2HF_2SE:  "<<Form("%.5f",SPdenom2HF_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom2HF_2SE->GetMeanError())<<endl;
    cout<<"SPdenom2HF_3SE:  "<<Form("%.5f",SPdenom2HF_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2HF_3SE[1]->GetMeanError())<<endl;
    cout<<"SPdenom2Trk_3SE: "<<Form("%.5f",SPdenom2Trk_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2Trk_3SE[1]->GetMeanError())<<endl;
    cout<<"SPdenom3HF_2SE:  "<<Form("%.5f",SPdenom3HF_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom3HF_2SE->GetMeanError())<<endl;
    cout<<"SPdenom3HF_3SE:  "<<Form("%.5f",SPdenom3HF_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom3HF_3SE[1]->GetMeanError())<<endl;
    cout<<"SPdenom3Trk_3SE: "<<Form("%.5f",SPdenom3Trk_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom3Trk_3SE[1]->GetMeanError())<<endl;
    cout<<"SPdenom1_mix:    "<<Form("%.5f",SPdenom1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix[1]->GetMeanError())<<endl;
    cout<<"SPdenom1_mix_3:  "<<Form("%.5f",SPdenom1_mix_3[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix_3[1]->GetMeanError())<<endl;
    cout<<"SPv1num:         "<<Form("%.5f",SPv1num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1num[1]->GetMeanError())<<endl;
    cout<<"SPv2num:         "<<Form("%.5f",SPv2num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2num[1]->GetMeanError())<<endl;
    cout<<"SPv3num:         "<<Form("%.5f",SPv3num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv3num[1]->GetMeanError())<<endl;
    cout<<"SPv1numMix:      "<<Form("%.5f",SPv1numMix[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix[1]->GetMeanError())<<endl;
    cout<<"SPv1numMix_3:    "<<Form("%.5f",SPv1numMix_3[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix_3[1]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"Q1Q1Q2:          "<<Form("%.5f",hQ1Q1Q2[1]->GetMean())<<" +/- "<<Form("%.5f",hQ1Q1Q2[1]->GetMeanError())<<endl;
    cout<<"Q1Q1Q2norm:      "<<Form("%.5f",hQ1Q1Q2norm[1]->GetMean())<<" +/- "<<Form("%.5f",hQ1Q1Q2norm[1]->GetMeanError())<<endl;
    cout<<"q112:            "<<Form("%.5f",hq112->GetMean())<<" +/- "<<Form("%.5f",hq112->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"EPv1_2SE:        "<<Form("%.5f",EPv1_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_2SE[1]->GetMeanError())<<endl;
    cout<<"EPv1_3SE:        "<<Form("%.5f",EPv1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_3SE[1]->GetMeanError())<<endl;
    cout<<"EPv1_mix:        "<<Form("%.5f",EPv1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix[1]->GetMeanError())<<endl;
    cout<<"EPv1_mix_3:      "<<Form("%.5f",EPv1_mix_3[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix_3[1]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"SPv1_2SE:        "<<Form("%.5f",SPv1_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_2SE[1]->GetMeanError())<<endl;
    cout<<"SPv1_3SE:        "<<Form("%.5f",SPv1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_3SE[1]->GetMeanError())<<endl;
    cout<<"SPv1_mix:        "<<Form("%.5f",SPv1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix[1]->GetMeanError())<<endl;
    cout<<"SPv1_mix_3:      "<<Form("%.5f",SPv1_mix_3[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix_3[1]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"EPv2_2SE:        "<<Form("%.5f",EPv2_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2_2SE[1]->GetMeanError())<<endl;
    cout<<"EPv2_3SE:        "<<Form("%.5f",EPv2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2_3SE[1]->GetMeanError())<<endl;
    cout<<"SPv2_2SE:        "<<Form("%.5f",SPv2_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2_2SE[1]->GetMeanError())<<endl;
    cout<<"SPv2_3SE:        "<<Form("%.5f",SPv2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2_3SE[1]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"EPv3_2SE:        "<<Form("%.5f",EPv3_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv3_2SE[1]->GetMeanError())<<endl;
    cout<<"EPv3_3SE:        "<<Form("%.5f",EPv3_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv3_3SE[1]->GetMeanError())<<endl;
    cout<<"SPv3_2SE:        "<<Form("%.5f",SPv3_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv3_2SE[1]->GetMeanError())<<endl;
    cout<<"SPv3_3SE:        "<<Form("%.5f",SPv3_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv3_3SE[1]->GetMeanError())<<endl;
    cout<<"\n --- Comparisons to input vn --- \n"<<endl;
    if (!isodd) {
        cout<<Form("EPv1_2SE/v1in  (HFm): %0.4f",EPv1_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_2SE[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("EPv1_3SE/v1in  (HFm): %0.4f",EPv1_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_3SE[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("EPv1_mix/v1in  (HFm): %0.4f",EPv1_mix[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_mix[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("EPv1_mix_3/v1in  (HFm): %0.4f",EPv1_mix_3[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix_3[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_mix_3[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix_3[1]->GetMeanError()/v1in)<<endl;
        cout<<"   -----        "<<endl;
        cout<<Form("SPv1_2SE/v1in  (HFm): %0.4f",SPv1_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_2SE[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("SPv1_3SE/v1in  (HFm): %0.4f",SPv1_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_3SE[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("SPv1_mix/v1in  (HFm): %0.4f",SPv1_mix[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_mix[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("SPv1_mix_3/v1in  (HFm): %0.4f",SPv1_mix_3[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix_3[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_mix_3[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix_3[1]->GetMeanError()/v1in)<<endl;
        cout<<"   -----        "<<endl;
    }
    cout<<Form("EPv2_2SE/v2in  (HFm): %0.4f",EPv2_2SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_2SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",EPv2_2SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_2SE[1]->GetMeanError()/v2in)<<endl;
    cout<<Form("EPv2_3SE/v2in  (HFm): %0.4f",EPv2_3SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_3SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",EPv2_3SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_3SE[1]->GetMeanError()/v2in)<<endl;
    cout<<Form("SPv2_2SE/v2in  (HFm): %0.4f",SPv2_2SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_2SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",SPv2_2SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_2SE[1]->GetMeanError()/v2in)<<endl;
    cout<<Form("SPv2_3SE/v2in  (HFm): %0.4f",SPv2_3SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_3SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",SPv2_3SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_3SE[1]->GetMeanError()/v2in)<<endl;
    cout<<"   -----        "<<endl;
    cout<<Form("EPv3_2SE/v3in  (HFm): %0.4f",EPv3_2SE[0]->GetMean()/v3in)<<Form(" +/- %0.4f",EPv3_2SE[0]->GetMeanError()/v3in)<<Form("\t (HFp): %0.4f",EPv3_2SE[1]->GetMean()/v3in)<<Form(" +/- %0.4f",EPv3_2SE[1]->GetMeanError()/v3in)<<endl;
    cout<<Form("EPv3_3SE/v3in  (HFm): %0.4f",EPv3_3SE[0]->GetMean()/v3in)<<Form(" +/- %0.4f",EPv3_3SE[0]->GetMeanError()/v3in)<<Form("\t (HFp): %0.4f",EPv3_3SE[1]->GetMean()/v3in)<<Form(" +/- %0.4f",EPv3_3SE[1]->GetMeanError()/v3in)<<endl;
    cout<<Form("SPv3_2SE/v3in  (HFm): %0.4f",SPv3_2SE[0]->GetMean()/v3in)<<Form(" +/- %0.4f",SPv3_2SE[0]->GetMeanError()/v3in)<<Form("\t (HFp): %0.4f",SPv3_2SE[1]->GetMean()/v3in)<<Form(" +/- %0.4f",SPv3_2SE[1]->GetMeanError()/v3in)<<endl;
    cout<<Form("SPv3_3SE/v3in  (HFm): %0.4f",SPv3_3SE[0]->GetMean()/v3in)<<Form(" +/- %0.4f",SPv3_3SE[0]->GetMeanError()/v3in)<<Form("\t (HFp): %0.4f",SPv3_3SE[1]->GetMean()/v3in)<<Form(" +/- %0.4f",SPv3_3SE[1]->GetMeanError()/v3in)<<endl;
    cout<<"   -----         "<<endl;
    cout<<"v1: "<<setv1<<"\tv2: "<<setv2<<"\tq112: "<<Form("%.5f",hq112->GetMean())<<" +/- "<<Form("%.5f",hq112->GetMeanError())<<endl;
    cout<<"\n ...done \n"<<endl;

}
