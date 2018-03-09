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
# include "HiEvtPlaneFlatten.h"

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

//-- The following values assume a centrality of 20-30%
// static const double v1oddin[] = {
    // -0.0373, -0.0324, -0.0281, -0.0232,	-0.0189, -0.0140,
    // -0.0097, -0.0049, -0.0006,	0.0043,  0.0054,  0.0049,
    //  0.0031,  0.0012, -0.0012, -0.0031,	-0.0049, -0.0054,
    // -0.0043,  0.0006,  0.0049,  0.0097,	 0.0140,  0.0189,
    //  0.0232,  0.0281,  0.0324,  0.0373};
//-- Values assuming a linear v1odd with no turn-around behavior (no decorrelation correction)
// static const double v1oddin[] = {
//      0.027,  0.025,  0.023,  0.021,  0.019,  0.017,
//      0.015,  0.013,  0.011,  0.009,  0.007,  0.005,
//      0.003,  0.001, -0.001, -0.003, -0.005, -0.007,
//     -0.009, -0.011, -0.013, -0.015, -0.017, -0.019,
//     -0.021, -0.023, -0.025, -0.027};
//-- Values assuming a linear v1odd with no turn-around behavior (decorrelation correction applied)
static const double v1oddin[] = {
     0.00378,  0.00350,  0.00322,  0.00294,  0.00266,  0.00238,
     0.00210,  0.00182,  0.00154,  0.00126,  0.00098,  0.00070,
     0.00042,  0.00014, -0.00014, -0.00042, -0.00070, -0.00098,
    -0.00126, -0.00154, -0.00182, -0.00210, -0.00238, -0.00266,
    -0.00294, -0.00322, -0.00350, -0.00378};

static const double v1Ebins[] = {
        -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0,
        -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,
         2.4,  2.8,  3.2,  3.6,  4.0,  4.4,  4.8,  5.2,  5.6};

// 0 to 5%
//static const double v1Emult[] = {
//         369, 355, 386, 439, 498, 554, 600, 635, 657, 670, 674, 674, 671, 669, 667, 668, 668, 668, 665, 655, 636, 606, 562, 507, 442, 375, 318, 288};

// 5 to 10%
// static const double v1Emult[] = {
//          312, 297, 320, 363, 411, 456, 493, 520, 537, 547, 551, 551, 550, 549, 548, 548, 546, 541, 532, 515, 488, 452, 405, 353, 300, 260, 247, 287};
//
// // 10 to 20%
// static const double v1Emult[] = {
//          235, 225, 242, 274, 310, 344, 371, 391, 403, 409, 411, 409, 407, 405, 405, 405, 406, 408, 407, 402, 391, 374, 348, 315, 276, 236, 204, 191};
//
// 20 - 30%
static const int v1Emult[] = {
         129, 149, 170, 191, 213, 234, 253, 270, 281, 286, 285, 280, 273, 268, 268, 273, 280,  285, 286, 281, 270, 253, 234, 213, 191, 170, 149, 129};
//
// // 30 to 40%
// static const double v1Emult[] = {
//          128, 106, 107, 120, 138, 154, 168, 176, 181, 181, 180, 177, 174, 173, 172, 174, 176, 179, 181, 180, 176, 168, 155, 140, 123, 111, 111, 134};
//
// // 40 to 50%
// static const double v1Emult[] = {
//          86, 71, 70, 77, 87, 96, 104, 109, 111, 111, 110, 108, 106, 105, 105, 106, 108, 110, 111, 111, 109, 104, 96, 87, 77, 70, 70, 83};
//
// // 50 to 60%
// static const double v1Emult[] = {
//          52, 42, 41, 45, 50, 56, 60, 62, 64, 63, 62, 61, 60, 59, 59, 60, 61, 62, 63, 64, 63, 60, 56, 50, 45, 41, 41, 49};
//
// // 60 to 70%
// static const double v1Emult[] = {
//          34, 27, 25, 25, 27, 30, 31, 32, 33, 33, 32, 31, 30, 30, 30, 30, 31, 32, 33, 33, 32, 31, 29, 27, 24, 22, 22, 26};


static const Double_t phiMinHole[] = {0.0,  2.2,  0.0,  0.0,  0.0};
static const Double_t phiMaxHole[] = {0.0,  2.6,  0.0,  0.0,  0.0};
static const Double_t etaMinHole[] = {0.0, -2.0,  0.0,  0.0,  0.0};
static const Double_t etaMaxHole[] = {0.0, -0.4,  0.0,  0.0,  0.0};

const std::string SideName[] = { "min", "pos" };
const std::string EPName[] = { "HFm", "trackm", "trackp", "HFp", "trkmid" };
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
TH1D * hsub_phiPsi3[numEP];
TH1D * hsub_Psi1[numEP];
TH1D * hsub_Psi2[numEP];
TH1D * hsub_Psi3[numEP];
TH1D * hsub_Psi1lab[numEP];
TH1D * hsub_Psi2lab[numEP];
TH1D * hsub_Psi3lab[numEP];
TH1D * hsub_ptAve[numEP];
TH1D * hsub_ptAve2[numEP];
TH1D * hsub_weights[numEP];
TH1D * hsub_wpt[numEP];
TH2D * hsub_pt2D[numEP];

//-- Q-vector correlations

TH1D * hQ1nA[2][netabins];
TH1D * hQ1AB[2][netabins];
TH1D * hQ1AC[2][netabins];
TH1D * hQ1BC[2][netabins];

TH1D * hQ2nA[2][netabins];
TH1D * hQ2AB[2][netabins];
TH1D * hQ2AC[2][netabins];
TH1D * hQ2BC[2][netabins];

TH1D * hQ3nA[2][netabins];
TH1D * hQ3AB[2][netabins];
TH1D * hQ3AC[2][netabins];
TH1D * hQ3BC[2][netabins];

TH1D * hQmixnAC[2][netabins];
TH1D * hQmixABC[2][netabins];

TH1D * hQmixnAC_3[2][netabins];
TH1D * hQmixABC_3[2][netabins];

TH1D * hQ1nA_final[2];
TH1D * hQ1AB_final[2];
TH1D * hQ1AC_final[2];
TH1D * hQ1BC_final[2];

TH1D * hQ2AB_final[2];
TH1D * hQ2AC_final[2];
TH1D * hQ2BC_final[2];

TH1D * hQ3AB_final[2];
TH1D * hQ3AC_final[2];
TH1D * hQ3BC_final[2];

TH1D * hQmixnAC_final[2];
TH1D * hQmixABC_final[2];

TH1D * hQmixnAC_3_final[2];
TH1D * hQmixABC_3_final[2];

TH1D * hQ1Q1Q2[2];

//-- normalized Q-vector correlations

TH1D * hQ1nAnorm[2][netabins];
TH1D * hQ1ABnorm[2][netabins];
TH1D * hQ1ACnorm[2][netabins];
TH1D * hQ1BCnorm[2][netabins];

TH1D * hQ2nAnorm[2][netabins];
TH1D * hQ2ABnorm[2][netabins];
TH1D * hQ2ACnorm[2][netabins];
TH1D * hQ2BCnorm[2][netabins];

TH1D * hQ3nAnorm[2][netabins];
TH1D * hQ3ABnorm[2][netabins];
TH1D * hQ3ACnorm[2][netabins];
TH1D * hQ3BCnorm[2][netabins];

TH1D * hQmixnACnorm[2][netabins];
TH1D * hQmixABCnorm[2][netabins];

TH1D * hQmixnAC_3norm[2][netabins];
TH1D * hQmixABC_3norm[2][netabins];

TH1D * hQ1nAnorm_final[2];
TH1D * hQ1ABnorm_final[2];
TH1D * hQ1ACnorm_final[2];
TH1D * hQ1BCnorm_final[2];

TH1D * hQ2ABnorm_final[2];
TH1D * hQ2ACnorm_final[2];
TH1D * hQ2BCnorm_final[2];

TH1D * hQ3ABnorm_final[2];
TH1D * hQ3ACnorm_final[2];
TH1D * hQ3BCnorm_final[2];

TH1D * hQmixnACnorm_final[2];
TH1D * hQmixABCnorm_final[2];

TH1D * hQmixnAC_3norm_final[2];
TH1D * hQmixABC_3norm_final[2];

TH1D * hQ1Q1Q2norm[2];
TH1D * hq112;

//-- event plane and scalar-product vn

TH1D * rescor1_2SE;
TH1D * rescor1_3SE[2];
TH1D * rescor1_mix[2];
TH1D * rescor1_mix3[2];
TH1D * rescor2HF_2SE;
TH1D * rescor2HF_3SE[2];
TH1D * rescor2Trk_3SE[2];
TH1D * rescor3HF_2SE;
TH1D * rescor3HF_3SE[2];
TH1D * rescor3Trk_3SE[2];
TH1D * EPv1obs[2];
TH1D * EPv1obsMix[2];
TH1D * EPv1obsMix_3[2];
TH1D * EPv2obs[2];
TH1D * EPv3obs[2];
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
TH1D * SPdenom2HF_3SE[2];
TH1D * SPdenom2Trk_3SE[2];
TH1D * SPdenom3HF_2SE;
TH1D * SPdenom3HF_3SE[2];
TH1D * SPdenom3Trk_3SE[2];
TH1D * SPv1num[2];
TH1D * SPv1numMix[2];
TH1D * SPv1numMix_3[2];
TH1D * SPv2num[2];
TH1D * SPv3num[2];
TH1D * SPv1_2SE[2];
TH1D * SPv1_3SE[2];
TH1D * SPv1_mix[2];
TH1D * SPv1_mix_3[2];
TH1D * SPv2_2SE[2];
TH1D * SPv2_3SE[2];
TH1D * SPv3_2SE[2];
TH1D * SPv3_3SE[2];

//-- differential vn(pT)

TH1D * EPv1obs_pt[2][nptbins];
TH1D * EPv1obsMix_pt[2][nptbins];
TH1D * EPv1obsMix_3_pt[2][nptbins];
TH1D * EPv2obs_pt[2][nptbins];
TH1D * EPv3obs_pt[2][nptbins];
TH1D * EPv1_2SE_pt[2][nptbins];
TH1D * EPv1_3SE_pt[2][nptbins];
TH1D * EPv1_mix_pt[2][nptbins];
TH1D * EPv1_mix_3_pt[2][nptbins];
TH1D * EPv2_2SE_pt[2][nptbins];
TH1D * EPv2_3SE_pt[2][nptbins];
TH1D * EPv3_2SE_pt[2][nptbins];
TH1D * EPv3_3SE_pt[2][nptbins];

TH1D * SPv1num_pt[2][nptbins];
TH1D * SPv1numMix_pt[2][nptbins];
TH1D * SPv1numMix_3_pt[2][nptbins];
TH1D * SPv2num_pt[2][nptbins];
TH1D * SPv3num_pt[2][nptbins];
TH1D * SPv1_2SE_pt[2][nptbins];
TH1D * SPv1_3SE_pt[2][nptbins];
TH1D * SPv1_mix_pt[2][nptbins];
TH1D * SPv1_mix_3_pt[2][nptbins];
TH1D * SPv2_2SE_pt[2][nptbins];
TH1D * SPv2_3SE_pt[2][nptbins];
TH1D * SPv3_2SE_pt[2][nptbins];
TH1D * SPv3_3SE_pt[2][nptbins];

TH1D * hqcnt_pt[2][nptbins];
TH1D * hncnt_pt[2][nptbins];

TH1D * DiffEPv1obs_pt[2];
TH1D * DiffEPv1obsMix_pt[2];
TH1D * DiffEPv1obsMix_3_pt[2];
TH1D * DiffEPv2obs_pt[2];
TH1D * DiffEPv3obs_pt[2];
TH1D * DiffEPv1_2SE_pt[2];
TH1D * DiffEPv1_3SE_pt[2];
TH1D * DiffEPv1_mix_pt[2];
TH1D * DiffEPv1_mix_3_pt[2];
TH1D * DiffEPv2_2SE_pt[2];
TH1D * DiffEPv2_3SE_pt[2];
TH1D * DiffEPv3_2SE_pt[2];
TH1D * DiffEPv3_3SE_pt[2];

TH1D * DiffSPv1num_pt[2];
TH1D * DiffSPv1numMix_pt[2];
TH1D * DiffSPv1numMix_3_pt[2];
TH1D * DiffSPv2num_pt[2];
TH1D * DiffSPv3num_pt[2];
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

TH1D * EPv1obs_eta[2][netabins];
TH1D * EPv1obsMix_eta[2][netabins];
TH1D * EPv1obsMix_3_eta[2][netabins];
TH1D * EPv1_2SE_eta[2][netabins];
TH1D * EPv1_3SE_eta[2][netabins];
TH1D * EPv1_mix_eta[2][netabins];
TH1D * EPv1_mix_3_eta[2][netabins];

TH1D * SPv1num_eta[2][netabins];
TH1D * SPv1numMix_eta[2][netabins];
TH1D * SPv1numMix_3_eta[2][netabins];
TH1D * SPv1_2SE_eta[2][netabins];
TH1D * SPv1_3SE_eta[2][netabins];
TH1D * SPv1_mix_eta[2][netabins];
TH1D * SPv1_mix_3_eta[2][netabins];

TH1D * hqcnt_eta[2][netabins];
TH1D * hncnt_eta[2][netabins];

TH1D * DiffEPv1obs_eta[2];
TH1D * DiffEPv1obsMix_eta[2];
TH1D * DiffEPv1obsMix_3_eta[2];
TH1D * DiffEPv1_2SE_eta[2];
TH1D * DiffEPv1_3SE_eta[2];
TH1D * DiffEPv1_mix_eta[2];
TH1D * DiffEPv1_mix_3_eta[2];

TH1D * DiffSPv1num_eta[2];
TH1D * DiffSPv1numMix_eta[2];
TH1D * DiffSPv1numMix_3_eta[2];
TH1D * DiffSPv1_2SE_eta[2];
TH1D * DiffSPv1_3SE_eta[2];
TH1D * DiffSPv1_mix_eta[2];
TH1D * DiffSPv1_mix_3_eta[2];

TH1D * Diffqcnt_eta[2];
TH1D * Diffncnt_eta[2];

//-- event plane correlations

TH2D * corr_HFm1_HFp1;
TH2D * corr_trackmid1_HFm1;
TH2D * corr_trackmid1_HFp1;
TH2D * corr_trackm1_HFm1;
TH2D * corr_trackm1_trackp1;
TH2D * corr_trackm1_HFp1;
TH2D * corr_trackp1_HFm1;
TH2D * corr_trackp1_trackm1;
TH2D * corr_trackp1_HFp1;

TH2D * corr_HFm1_trackm2;
TH2D * corr_HFm1_trackp2;
TH2D * corr_HFp1_trackm2;
TH2D * corr_HFp1_trackp2;
TH2D * corr_HFm1_HFp2;
TH2D * corr_HFp1_HFm2;
TH2D * corr_trackmid1_HFm2;
TH2D * corr_trackmid1_HFp2;
TH2D * corr_trackm1_HFm2;
TH2D * corr_trackm1_trackp2;
TH2D * corr_trackm1_HFp2;
TH2D * corr_trackp1_HFm2;
TH2D * corr_trackp1_trackm2;
TH2D * corr_trackp1_HFp2;

TH2D * corr_HFm2_HFp2;
TH2D * corr_trackmid2_HFm2;
TH2D * corr_trackmid2_HFp2;
TH2D * corr_trackm2_HFm2;
TH2D * corr_trackm2_trackp2;
TH2D * corr_trackm2_HFp2;
TH2D * corr_trackp2_HFm2;
TH2D * corr_trackp2_HFp2;

//-- cosines of correlations

TH1D * cos_HFm1_HFp1;
TH1D * cos_trackmid1_HFm1;
TH1D * cos_trackmid1_HFp1;
TH1D * cos_trackm1_HFm1;
TH1D * cos_trackm1_trackp1;
TH1D * cos_trackm1_HFp1;
TH1D * cos_trackp1_HFm1;
TH1D * cos_trackp1_trackm1;
TH1D * cos_trackp1_HFp1;

TH1D * cos_HFm1_trackm2;
TH1D * cos_HFm1_trackp2;
TH1D * cos_HFp1_trackm2;
TH1D * cos_HFp1_trackp2;
TH1D * cos_HFm1_HFp2;
TH1D * cos_HFp1_HFm2;
TH1D * cos_trackmid1_HFm2;
TH1D * cos_trackmid1_HFp2;
TH1D * cos_trackm1_HFm2;
TH1D * cos_trackm1_trackp2;
TH1D * cos_trackm1_HFp2;
TH1D * cos_trackp1_HFm2;
TH1D * cos_trackp1_trackm2;
TH1D * cos_trackp1_HFp2;

TH1D * cos_HFm2_HFp2;
TH1D * cos_trackmid2_HFm2;
TH1D * cos_trackmid2_HFp2;
TH1D * cos_trackm2_HFm2;
TH1D * cos_trackm2_trackp2;
TH1D * cos_trackm2_HFp2;
TH1D * cos_trackp2_HFm2;
TH1D * cos_trackp2_HFp2;

HiEvtPlaneFlatten * flat1[6];
HiEvtPlaneFlatten * flat2[6];
HiEvtPlaneFlatten * flat3[6];

TF1 * fnTheory;

TFile * tfout;
Int_t iseed = 0;
Int_t counter = 0;


Double_t bounds(int ord, double ang) {
    while (ang >  TMath::Pi()/ord) ang-=TMath::TwoPi()/ord;
    while (ang < -TMath::Pi()/ord) ang+=TMath::TwoPi()/ord;
    return ang;
}

typedef complex<double> comp;
comp zero (0,0);

void Setup()
{
    inputParms = new TH1D("Input_Parameters", "", 11, 0, 11);
    const char * inputLabel[] = {"n_{evts}", "odd v_{1}^{in}", "v_{2}^{in}", "v_{3}^{in}", "#eta-weights", "p_{T}-weights", "mom-cons", "recenter", "flatten", "#eta_{trk}^{min}", "#eta_{trk}^{max}"};
    for (int i = 1; i<=11; i++) inputParms->SetBinContent(i, 0);
    for (int i = 1; i<=11; i++) inputParms->GetXaxis()->SetBinLabel(i, inputLabel[i-1]);

    hinit_v1in = new TH1D("v1in", "", nv1Ebins, v1Ebins);
    rescor1_2SE = new TH1D("rescor1_2SE", "", 200, 0, 1);
    rescor2HF_2SE = new TH1D("rescor2HF_2SE", "", 200, 0, 1);
    rescor3HF_2SE = new TH1D("rescor3HF_2SE", "", 200, 0, 1);
    SPdenom1_2SE = new TH1D("SPdenom1_2SE", "", 200, 0, 0.1);
    SPdenom2HF_2SE = new TH1D("SPdenom2HF_2SE", "", 200, 0, 0.1);
    SPdenom3HF_2SE = new TH1D("SPdenom3HF_2SE", "", 200, 0, 0.1);

    hq112 = new TH1D("q112", "", 200, -0.1, 0.2);
    for (int iside = 0; iside<2; iside++) {
        for (int ebin = 0; ebin<netabins; ebin++) {

            hQ1nA[iside][ebin] = new TH1D(Form("Q1nA_%s_%d",SideName[iside].data(),ebin), "", 150, -1e-4, 1e-4);
            hQ2nA[iside][ebin] = new TH1D(Form("Q2nA_%s_%d",SideName[iside].data(),ebin), "", 150, -1e-3, 1e-3);
            hQ3nA[iside][ebin] = new TH1D(Form("Q3nA_%s_%d",SideName[iside].data(),ebin), "", 150, -1e-3, 1e-3);
            hQmixnAC[iside][ebin] = new TH1D(Form("QmixnAC_%s_%d",SideName[iside].data(),ebin), "", 150, -5e-5, 5e-5);
            hQmixnAC_3[iside][ebin] = new TH1D(Form("QmixnAC_3_%s_%d",SideName[iside].data(),ebin), "", 150, -5e-5, 5e-5);

            hQ1nAnorm[iside][ebin] = new TH1D(Form("Q1nAnorm_%s_%d",SideName[iside].data(),ebin), "", 150, -1e-4, 1e-4);
            hQ2nAnorm[iside][ebin] = new TH1D(Form("Q2nAnorm_%s_%d",SideName[iside].data(),ebin), "", 150, -1e-3, 1e-3);
            hQ3nAnorm[iside][ebin] = new TH1D(Form("Q3nAnorm_%s_%d",SideName[iside].data(),ebin), "", 150, -1e-3, 1e-3);
            hQmixnACnorm[iside][ebin] = new TH1D(Form("QmixnACnorm_%s_%d",SideName[iside].data(),ebin), "", 150, -5e-5, 5e-5);
            hQmixnAC_3norm[iside][ebin] = new TH1D(Form("QmixnAC_3norm_%s_%d",SideName[iside].data(),ebin), "", 150, -5e-5, 5e-5);
        }

        hQ1nA_final[iside] = new TH1D(Form("Q1nA_%s",SideName[iside].data()), "", netabins, etabins);
        hQ1AB_final[iside] = new TH1D(Form("Q1AB_%s",SideName[iside].data()), "", 200, -1e-4, 1e-4);
        hQ1AC_final[iside] = new TH1D(Form("Q1AC_%s",SideName[iside].data()), "", 200, -1e-4, 1e-4);
        hQ1BC_final[iside] = new TH1D(Form("Q1BC_%s",SideName[iside].data()), "", 200, -1e-4, 1e-4);
        hQ2AB_final[iside] = new TH1D(Form("Q2AB_%s",SideName[iside].data()), "", 200, -5e-3, 5e-3);
        hQ2AC_final[iside] = new TH1D(Form("Q2AC_%s",SideName[iside].data()), "", 200, -5e-3, 5e-3);
        hQ2BC_final[iside] = new TH1D(Form("Q2BC_%s",SideName[iside].data()), "", 200, -5e-3, 5e-3);
        hQ3AB_final[iside] = new TH1D(Form("Q3AB_%s",SideName[iside].data()), "", 200, -5e-3, 5e-3);
        hQ3AC_final[iside] = new TH1D(Form("Q3AC_%s",SideName[iside].data()), "", 200, -5e-3, 5e-3);
        hQ3BC_final[iside] = new TH1D(Form("Q3BC_%s",SideName[iside].data()), "", 200, -5e-3, 5e-3);
        hQmixnAC_final[iside] = new TH1D(Form("QmixnAC_%s",SideName[iside].data()), "", netabins, etabins);
        hQmixABC_final[iside] = new TH1D(Form("QmixABC_%s",SideName[iside].data()), "", 200, -5e-5, 5e-5);
        hQmixnAC_3_final[iside] = new TH1D(Form("QmixnAC_3_final_%s",SideName[iside].data()), "", netabins, etabins);
        hQmixABC_3_final[iside] = new TH1D(Form("QmixABC_3_final_%s",SideName[iside].data()), "", 200, -5e-5, 5e-5);
        hQ1Q1Q2[iside] = new TH1D(Form("Q1Q1Q2_%s",SideName[iside].data()), "", 200, -0.001, 0.001);

        hQ1nAnorm_final[iside] = new TH1D(Form("Q1nAnorm_%s",SideName[iside].data()), "", netabins, etabins);
        hQ1ABnorm_final[iside] = new TH1D(Form("Q1ABnorm_%s",SideName[iside].data()), "", 200, -0.1, 0.1);
        hQ1ACnorm_final[iside] = new TH1D(Form("Q1ACnorm_%s",SideName[iside].data()), "", 200, -0.1, 0.1);
        hQ1BCnorm_final[iside] = new TH1D(Form("Q1BCnorm_%s",SideName[iside].data()), "", 200, -0.1, 0.1);
        hQ2ABnorm_final[iside] = new TH1D(Form("Q2ABnorm_%s",SideName[iside].data()), "", 200, -1, 1);
        hQ2ACnorm_final[iside] = new TH1D(Form("Q2ACnorm_%s",SideName[iside].data()), "", 200, -1, 1);
        hQ2BCnorm_final[iside] = new TH1D(Form("Q2BCnorm_%s",SideName[iside].data()), "", 200, -1, 1);
        hQ3ABnorm_final[iside] = new TH1D(Form("Q3ABnorm_%s",SideName[iside].data()), "", 200, -1, 1);
        hQ3ACnorm_final[iside] = new TH1D(Form("Q3ACnorm_%s",SideName[iside].data()), "", 200, -1, 1);
        hQ3BCnorm_final[iside] = new TH1D(Form("Q3BCnorm_%s",SideName[iside].data()), "", 200, -1, 1);
        hQmixnACnorm_final[iside] = new TH1D(Form("QmixnACnorm_%s",SideName[iside].data()), "", netabins, etabins);
        hQmixABCnorm_final[iside] = new TH1D(Form("QmixABCnorm_%s",SideName[iside].data()), "", 200, -5e-5, 5e-5);
        hQmixnAC_3norm_final[iside] = new TH1D(Form("QmixnAC_3norm_final_%s",SideName[iside].data()), "", netabins, etabins);
        hQmixABC_3norm_final[iside] = new TH1D(Form("QmixABC_3norm_final_%s",SideName[iside].data()), "", 200, -5e-5, 5e-5);
        hQ1Q1Q2norm[iside] = new TH1D(Form("Q1Q1Q2norm_%s",SideName[iside].data()), "", 200, -0.12, 0.12);

        rescor1_3SE[iside] = new TH1D(Form("rescor1_3SE_%s",SideName[iside].data()), "", 200, 0, 1);
        rescor1_mix[iside] = new TH1D(Form("rescor1_mix_%s",SideName[iside].data()), "", 200, 0, 1);
        rescor1_mix3[iside] = new TH1D(Form("rescor1_mix3_%s",SideName[iside].data()), "", 200, 0, 1);
        rescor2HF_3SE[iside] = new TH1D(Form("rescor2HF_3SE_%s",SideName[iside].data()), "", 200, 0, 1);
        rescor2Trk_3SE[iside] = new TH1D(Form("rescor2Trk_3SE_%s",SideName[iside].data()), "", 200, 0, 1);
        rescor3HF_3SE[iside] = new TH1D(Form("rescor3HF_3SE_%s",SideName[iside].data()), "", 200, 0, 1);
        rescor3Trk_3SE[iside] = new TH1D(Form("rescor3Trk_3SE_%s",SideName[iside].data()), "", 200, 0, 1);
        EPv1obs[iside] = new TH1D(Form("EPv1obs_%s",SideName[iside].data()), "", 200, -0.1, 0.1);
        EPv1obsMix[iside] = new TH1D(Form("EPv1obsMix_%s",SideName[iside].data()), "", 200, -0.1, 0.1);
        EPv1obsMix_3[iside] = new TH1D(Form("EPv1obsMix_3_%s",SideName[iside].data()), "", 200, -0.1, 0.1);
        EPv2obs[iside] = new TH1D(Form("EPv2obs_%s",SideName[iside].data()), "", 200, -0.1, 0.1);
        EPv3obs[iside] = new TH1D(Form("EPv3obs_%s",SideName[iside].data()), "", 200, -0.1, 0.1);
        EPv1_2SE[iside] = new TH1D(Form("EPv1_2SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        EPv1_3SE[iside] = new TH1D(Form("EPv1_3SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        EPv1_mix[iside] = new TH1D(Form("EPv1_mix_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        EPv1_mix_3[iside] = new TH1D(Form("EPv1_mix_3_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        EPv2_2SE[iside] = new TH1D(Form("EPv2_2SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        EPv2_3SE[iside] = new TH1D(Form("EPv2_3SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        EPv3_2SE[iside] = new TH1D(Form("EPv3_2SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        EPv3_3SE[iside] = new TH1D(Form("EPv3_3SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);

        SPdenom1_3SE[iside] = new TH1D(Form("SPdenom1_3SE_%s",SideName[iside].data()), "", 200, 0, 0.1);
        SPdenom1_mix[iside] = new TH1D(Form("SPdenom1_mix_%s",SideName[iside].data()), "", 200, 0, 0.1);
        SPdenom1_mix_3[iside] = new TH1D(Form("SPdenom1_mix_3_%s",SideName[iside].data()), "", 200, 0, 0.1);
        SPdenom2HF_3SE[iside] = new TH1D(Form("SPdenom2HF_3SE_%s",SideName[iside].data()), "", 200, 0, 0.1);
        SPdenom2Trk_3SE[iside] = new TH1D(Form("SPdenom2Trk_3SE_%s",SideName[iside].data()), "", 200, 0, 0.1);
        SPdenom3HF_3SE[iside] = new TH1D(Form("SPdenom3HF_3SE_%s",SideName[iside].data()), "", 200, 0, 0.1);
        SPdenom3Trk_3SE[iside] = new TH1D(Form("SPdenom3Trk_3SE_%s",SideName[iside].data()), "", 200, 0, 0.1);
        SPv1num[iside] = new TH1D(Form("SPv1num_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv1numMix[iside] = new TH1D(Form("SPv1numMix_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv1numMix_3[iside] = new TH1D(Form("SPv1numMix_3_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv2num[iside] = new TH1D(Form("SPv2num_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv3num[iside] = new TH1D(Form("SPv3num_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv1_2SE[iside] = new TH1D(Form("SPv1_2SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv1_3SE[iside] = new TH1D(Form("SPv1_3SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv1_mix[iside] = new TH1D(Form("SPv1_mix_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv1_mix_3[iside] = new TH1D(Form("SPv1_mix_3_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv2_2SE[iside] = new TH1D(Form("SPv2_2SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv2_3SE[iside] = new TH1D(Form("SPv2_3SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv3_2SE[iside] = new TH1D(Form("SPv3_2SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);
        SPv3_3SE[iside] = new TH1D(Form("SPv3_3SE_%s",SideName[iside].data()), "", 200, -0.15, 0.15);

        for (int pbin = 0; pbin<nptbins; pbin++) {
            EPv1obs_pt[iside][pbin] = new TH1D(Form("EPv1obs_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            EPv1obsMix_pt[iside][pbin] = new TH1D(Form("EPv1obsMix_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            EPv1obsMix_3_pt[iside][pbin] = new TH1D(Form("EPv1obsMix_3_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            EPv2obs_pt[iside][pbin] = new TH1D(Form("EPv2obs_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            EPv3obs_pt[iside][pbin] = new TH1D(Form("EPv3obs_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            EPv1_2SE_pt[iside][pbin] = new TH1D(Form("EPv1_2SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            EPv1_3SE_pt[iside][pbin] = new TH1D(Form("EPv1_3SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            EPv1_mix_pt[iside][pbin] = new TH1D(Form("EPv1_mix_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            EPv1_mix_3_pt[iside][pbin] = new TH1D(Form("EPv1_mix_3_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            EPv2_2SE_pt[iside][pbin] = new TH1D(Form("EPv2_2SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            EPv2_3SE_pt[iside][pbin] = new TH1D(Form("EPv2_3SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            EPv3_2SE_pt[iside][pbin] = new TH1D(Form("EPv3_2SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            EPv3_3SE_pt[iside][pbin] = new TH1D(Form("EPv3_3SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);

            SPv1num_pt[iside][pbin] = new TH1D(Form("SPv1num_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            SPv1numMix_pt[iside][pbin] = new TH1D(Form("SPv1numMix_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            SPv1numMix_3_pt[iside][pbin] = new TH1D(Form("SPv1numMix_3_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            SPv2num_pt[iside][pbin] = new TH1D(Form("SPv2num_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            SPv3num_pt[iside][pbin] = new TH1D(Form("SPv3num_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.1, 0.1);
            SPv1_2SE_pt[iside][pbin] = new TH1D(Form("SPv1_2SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            SPv1_3SE_pt[iside][pbin] = new TH1D(Form("SPv1_3SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            SPv1_mix_pt[iside][pbin] = new TH1D(Form("SPv1_mix_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            SPv1_mix_3_pt[iside][pbin] = new TH1D(Form("SPv1_mix_3_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            SPv2_2SE_pt[iside][pbin] = new TH1D(Form("SPv2_2SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            SPv2_3SE_pt[iside][pbin] = new TH1D(Form("SPv2_3SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            SPv3_2SE_pt[iside][pbin] = new TH1D(Form("SPv3_2SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);
            SPv3_3SE_pt[iside][pbin] = new TH1D(Form("SPv3_3SE_pt_%s_%d",SideName[iside].data(),pbin), "", 200, -0.15, 0.15);

            hqcnt_pt[iside][pbin] = new TH1D(Form("qcnt_pt_%s_%d",SideName[iside].data(),pbin), "", 1000, 0, 1e6);
            hncnt_pt[iside][pbin] = new TH1D(Form("ncnt_pt_%s_%d",SideName[iside].data(),pbin), "", 1000, 0, 1e9);
        }

        DiffEPv1obs_pt[iside] = new TH1D(Form("EPv1obs_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv1obsMix_pt[iside] = new TH1D(Form("EPv1obsMix_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv1obsMix_3_pt[iside] = new TH1D(Form("EPv1obsMix_3_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv2obs_pt[iside] = new TH1D(Form("EPv2obs_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv3obs_pt[iside] = new TH1D(Form("EPv3obs_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv1_2SE_pt[iside] = new TH1D(Form("EPv1_2SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv1_3SE_pt[iside] = new TH1D(Form("EPv1_3SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv1_mix_pt[iside] = new TH1D(Form("EPv1_mix_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv1_mix_3_pt[iside] = new TH1D(Form("EPv1_mix_3_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv2_2SE_pt[iside] = new TH1D(Form("EPv2_2SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv2_3SE_pt[iside] = new TH1D(Form("EPv2_3SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv3_2SE_pt[iside] = new TH1D(Form("EPv3_2SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffEPv3_3SE_pt[iside] = new TH1D(Form("EPv3_3SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);

        DiffSPv1num_pt[iside] = new TH1D(Form("SPv1num_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv1numMix_pt[iside] = new TH1D(Form("SPv1numMix_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv1numMix_3_pt[iside] = new TH1D(Form("SPv1numMix_3_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv2num_pt[iside] = new TH1D(Form("SPv2num_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv3num_pt[iside] = new TH1D(Form("SPv3num_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv1_2SE_pt[iside] = new TH1D(Form("SPv1_2SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv1_3SE_pt[iside] = new TH1D(Form("SPv1_3SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv1_mix_pt[iside] = new TH1D(Form("SPv1_mix_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv1_mix_3_pt[iside] = new TH1D(Form("SPv1_mix_3_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv2_2SE_pt[iside] = new TH1D(Form("SPv2_2SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv2_3SE_pt[iside] = new TH1D(Form("SPv2_3SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv3_2SE_pt[iside] = new TH1D(Form("SPv3_2SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        DiffSPv3_3SE_pt[iside] = new TH1D(Form("SPv3_3SE_pt_%s",SideName[iside].data()), "", nptbins, ptbins);

        Diffqcnt_pt[iside] = new TH1D(Form("Qcnt_pt_%s",SideName[iside].data()), "", nptbins, ptbins);
        Diffncnt_pt[iside] = new TH1D(Form("Ncnt_pt_%s",SideName[iside].data()), "", nptbins, ptbins);

        for (int ebin = 0; ebin<netabins; ebin++) {
            EPv1obs_eta[iside][ebin] = new TH1D(Form("EPv1obs_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.1, 0.1);
            EPv1obsMix_eta[iside][ebin] = new TH1D(Form("EPv1obsMix_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.1, 0.1);
            EPv1obsMix_3_eta[iside][ebin] = new TH1D(Form("EPv1obsMix_3_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.1, 0.1);
            EPv1_2SE_eta[iside][ebin] = new TH1D(Form("EPv1_2SE_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.15, 0.15);
            EPv1_3SE_eta[iside][ebin] = new TH1D(Form("EPv1_3SE_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.15, 0.15);
            EPv1_mix_eta[iside][ebin] = new TH1D(Form("EPv1_mix_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.15, 0.15);
            EPv1_mix_3_eta[iside][ebin] = new TH1D(Form("EPv1_mix_3_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.15, 0.15);

            SPv1num_eta[iside][ebin] = new TH1D(Form("SPv1num_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.1, 0.1);
            SPv1numMix_eta[iside][ebin] = new TH1D(Form("SPv1numMix_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.1, 0.1);
            SPv1numMix_3_eta[iside][ebin] = new TH1D(Form("SPv1num_Mix_3_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.1, 0.1);
            SPv1_2SE_eta[iside][ebin] = new TH1D(Form("SPv1_2SE_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.15, 0.15);
            SPv1_3SE_eta[iside][ebin] = new TH1D(Form("SPv1_3SE_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.15, 0.15);
            SPv1_mix_eta[iside][ebin] = new TH1D(Form("SPv1_mix_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.15, 0.15);
            SPv1_mix_3_eta[iside][ebin] = new TH1D(Form("SPv1_mix_3_eta_%s_%d",SideName[iside].data(),ebin), "", 200, -0.15, 0.15);

            hqcnt_eta[iside][ebin] = new TH1D(Form("qcnt_eta_%s_%d",SideName[iside].data(),ebin), "", 1000, 0, 1e6);
            hncnt_eta[iside][ebin] = new TH1D(Form("ncnt_eta_%s_%d",SideName[iside].data(),ebin), "", 1000, 0, 1e9);
        }

        DiffEPv1obs_eta[iside] = new TH1D(Form("EPv1obs_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffEPv1obsMix_eta[iside] = new TH1D(Form("EPv1obsMix_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffEPv1obsMix_3_eta[iside] = new TH1D(Form("EPv1obsMix_3_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffEPv1_2SE_eta[iside] = new TH1D(Form("EPv1_2SE_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffEPv1_3SE_eta[iside] = new TH1D(Form("EPv1_3SE_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffEPv1_mix_eta[iside] = new TH1D(Form("EPv1_mix_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffEPv1_mix_3_eta[iside] = new TH1D(Form("EPv1_mix_3_eta_%s",SideName[iside].data()), "", netabins, etabins);

        DiffSPv1num_eta[iside] = new TH1D(Form("SPv1num_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffSPv1numMix_eta[iside] = new TH1D(Form("SPv1numMix_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffSPv1numMix_3_eta[iside] = new TH1D(Form("SPv1numMix_3_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffSPv1_2SE_eta[iside] = new TH1D(Form("SPv1_2SE_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffSPv1_3SE_eta[iside] = new TH1D(Form("SPv1_3SE_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffSPv1_mix_eta[iside] = new TH1D(Form("SPv1_mix_eta_%s",SideName[iside].data()), "", netabins, etabins);
        DiffSPv1_mix_3_eta[iside] = new TH1D(Form("SPv1_mix_3_eta_%s",SideName[iside].data()), "", netabins, etabins);

        Diffqcnt_eta[iside] = new TH1D(Form("Qcnt_eta_%s",SideName[iside].data()), "Qcnt_eta", netabins, etabins);
        Diffncnt_eta[iside] = new TH1D(Form("Ncnt_eta_%s",SideName[iside].data()), "Ncnt_eta", netabins, etabins);
    }

    corr_HFm1_HFp1 = new TH2D("HFp1_vs_HFm1", "HFp1_vs_HFm1", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
    corr_HFm1_HFp1->Sumw2();
    corr_HFm1_HFp1->SetOption("colz");
    corr_trackmid1_HFm1 = (TH2D *) corr_HFm1_HFp1->Clone("HFm1_vs_trackmid1");
    corr_trackmid1_HFp1 = (TH2D *) corr_HFm1_HFp1->Clone("HFp1_vs_trackmid1");
    corr_trackm1_HFm1 = (TH2D *) corr_HFm1_HFp1->Clone("HFm1_vs_trackm1");
    corr_trackm1_trackp1 = (TH2D *) corr_HFm1_HFp1->Clone("trackp1_vs_trackm1");
    corr_trackm1_HFp1 = (TH2D *) corr_HFm1_HFp1->Clone("HFp1_vs_trackm1");
    corr_trackp1_HFm1 = (TH2D *) corr_HFm1_HFp1->Clone("HFm1_vs_trackp1");
    corr_trackp1_trackm1 = (TH2D *) corr_HFm1_HFp1->Clone("trackm1_vs_trackp1");
    corr_trackp1_HFp1 = (TH2D *) corr_HFm1_HFp1->Clone("HFp1_vs_trackp1");

    corr_HFm1_trackm2 = new TH2D("trackm2_vs_HFm1", "trackm2_vs_HFm1", 100, -TMath::Pi(), TMath::Pi(), 100, -0.5*TMath::Pi(), 0.5*TMath::Pi());
    corr_HFm1_trackm2->Sumw2();
    corr_HFm1_trackm2->SetOption("colz");
    corr_HFm1_trackp2 = (TH2D *) corr_HFm1_trackm2->Clone("trackp2_vs_HFm1");
    corr_HFp1_trackm2 = (TH2D *) corr_HFm1_trackm2->Clone("trackm2_vs_HFp1");
    corr_HFp1_trackp2 = (TH2D *) corr_HFm1_trackm2->Clone("trackp2_vs_HFp1");
    corr_HFm1_HFp2 = (TH2D *) corr_HFm1_trackm2->Clone("HFp2_vs_HFm1");
    corr_HFp1_HFm2 = (TH2D *) corr_HFm1_trackm2->Clone("HFm2_vs_HFp1");
    corr_trackmid1_HFm2 = (TH2D *) corr_HFm1_trackm2->Clone("HFm2_vs_trackmid1");
    corr_trackmid1_HFp2 = (TH2D *) corr_HFm1_trackm2->Clone("HFp2_vs_trackmid1");
    corr_trackm1_HFm2 = (TH2D *) corr_HFm1_trackm2->Clone("HFm2_vs_trackm1");
    corr_trackm1_trackp2 = (TH2D *) corr_HFm1_trackm2->Clone("trackp2_vs_trackm1");
    corr_trackm1_HFp2 = (TH2D *) corr_HFm1_trackm2->Clone("HFp2_vs_trackm1");
    corr_trackp1_HFm2 = (TH2D *) corr_HFm1_trackm2->Clone("HFm2_vs_trackp1");
    corr_trackp1_trackm2 = (TH2D *) corr_HFm1_trackm2->Clone("trackm2_vs_trackp1");
    corr_trackp1_HFp2 = (TH2D *) corr_HFm1_trackm2->Clone("HFp2_vs_trackp1");

    corr_HFm2_HFp2 = new TH2D("HFp2_vs_HFm2", "HFp2_vs_HFm2", 100, -0.5*TMath::Pi(), 0.5*TMath::Pi(), 100, -0.5*TMath::Pi(), 0.5*TMath::Pi());
    corr_HFm2_HFp2->Sumw2();
    corr_HFm2_HFp2->SetOption("colz");
    corr_trackmid2_HFm2 = (TH2D *) corr_HFm2_HFp2->Clone("HFm2_vs_trackmid2");
    corr_trackmid2_HFp2 = (TH2D *) corr_HFm2_HFp2->Clone("HFp2_vs_trackmid2");
    corr_trackm2_HFm2 = (TH2D *) corr_HFm2_HFp2->Clone("HFm2_vs_trackm2");
    corr_trackm2_trackp2 = (TH2D *) corr_HFm2_HFp2->Clone("trackp2_vs_trackm2");
    corr_trackm2_HFp2 = (TH2D *) corr_HFm2_HFp2->Clone("HFp2_vs_trackm2");
    corr_trackp2_HFm2 = (TH2D *) corr_HFm2_HFp2->Clone("HFm2_vs_trackp2");
    corr_trackp2_HFp2 = (TH2D *) corr_HFm2_HFp2->Clone("HFp2_vs_trackp2");

    //-- cosines of correlations

    cos_HFm1_HFp1 = new TH1D("cos_HFm1_HFp1", "cos_HFm1_HFp1", 150, -1.2, 1.2);
    cos_trackmid1_HFm1 = new TH1D("cos_trackmid1_HFm1", "cos_trackmid1_HFm1", 150, -1.2, 1.2);
    cos_trackmid1_HFp1 = new TH1D("cos_trackmid1_HFp1", "cos_trackmid1_HFp1", 150, -1.2, 1.2);
    cos_trackm1_HFm1 = new TH1D("cos_trackm1_HFm1", "cos_trackm1_HFm1", 150, -1.2, 1.2);
    cos_trackm1_trackp1 = new TH1D("cos_trackm1_trackp1", "cos_trackm1_trackp1", 150, -1.2, 1.2);
    cos_trackm1_HFp1 = new TH1D("cos_trackm1_HFp1", "cos_trackm1_HFp1", 150, -1.2, 1.2);
    cos_trackp1_HFm1 = new TH1D("cos_trackp1_HFm1", "cos_trackp1_HFm1", 150, -1.2, 1.2);
    cos_trackp1_trackm1 = new TH1D("cos_trackp1_trackm1", "cos_trackp1_trackm1", 150, -1.2, 1.2);
    cos_trackp1_HFp1 = new TH1D("cos_trackp1_HFp1", "cos_trackp1_HFp1", 150, -1.2, 1.2);

    cos_HFm1_trackm2 = new TH1D("cos_HFm1_trackm2", "cos_HFm1_trackm2", 150, -1.2, 1.2);
    cos_HFm1_trackp2 = new TH1D("cos_HFm1_trackp2", "cos_HFm1_trackp2", 150, -1.2, 1.2);
    cos_HFp1_trackm2 = new TH1D("cos_HFp1_trackm2", "cos_HFp1_trackm2", 150, -1.2, 1.2);
    cos_HFp1_trackp2 = new TH1D("cos_HFp1_trackp2", "cos_HFp1_trackp2", 150, -1.2, 1.2);
    cos_HFm1_HFp2 = new TH1D("cos_HFm1_HFp2", "cos_HFm1_HFp2", 150, -1.2, 1.2);
    cos_HFp1_HFm2 = new TH1D("cos_HFp1_HFm2", "cos_HFp1_HFm2", 150, -1.2, 1.2);
    cos_trackmid1_HFm2 = new TH1D("cos_trackmid1_HFm2", "cos_trackmid1_HFm2", 150, -1.2, 1.2);
    cos_trackmid1_HFp2 = new TH1D("cos_trackmid1_HFp2", "cos_trackmid1_HFp2", 150, -1.2, 1.2);
    cos_trackm1_HFm2 = new TH1D("cos_trackm1_HFm2", "cos_trackm1_HFm2", 150, -1.2, 1.2);
    cos_trackm1_trackp2 = new TH1D("cos_trackm1_trackp2", "cos_trackm1_trackp2", 150, -1.2, 1.2);
    cos_trackm1_HFp2 = new TH1D("cos_trackm1_HFp2", "cos_trackm1_HFp2", 150, -1.2, 1.2);
    cos_trackp1_HFm2 = new TH1D("cos_trackp1_HFm2", "cos_trackp1_HFm2", 150, -1.2, 1.2);
    cos_trackp1_trackm2 = new TH1D("cos_trackp1_trackm2", "cos_trackp1_trackm2", 150, -1.2, 1.2);
    cos_trackp1_HFp2 = new TH1D("cos_trackp1_HFp2", "cos_trackp1_HFp2", 150, -1.2, 1.2);

    cos_HFm2_HFp2 = new TH1D("cos_HFm2_HFp2", "cos_HFm2_HFp2", 150, -1.2, 1.2);
    cos_trackmid2_HFm2 = new TH1D("cos_trackmid2_HFm2", "cos_trackmid2_HFm2", 150, -1.2, 1.2);
    cos_trackmid2_HFp2 = new TH1D("cos_trackmid2_HFp2", "cos_trackmid2_HFp2", 150, -1.2, 1.2);
    cos_trackm2_HFm2 = new TH1D("cos_trackm2_HFm2", "cos_trackm2_HFm2", 150, -1.2, 1.2);
    cos_trackm2_trackp2 = new TH1D("cos_trackm2_trackp2", "cos_trackm2_trackp2", 150, -1.2, 1.2);
    cos_trackm2_HFp2 = new TH1D("cos_trackm2_HFp2", "cos_trackm2_HFp2", 150, -1.2, 1.2);
    cos_trackp2_HFm2 = new TH1D("cos_trackp2_HFm2", "cos_trackp2_HFm2", 150, -1.2, 1.2);
    cos_trackp2_HFp2 = new TH1D("cos_trackp2_HFp2", "cos_trackp2_HFp2", 150, -1.2, 1.2);
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
