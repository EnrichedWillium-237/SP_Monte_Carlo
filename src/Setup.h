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
static const double ecutmin[] = {-5.0, -2.4,  2.0,  3.0, -0.5};
static const double ecutmax[] = {-3.0, -2.0,  2.4,  5.0,  0.5};
static const double pcutmin[] = { 0.3,  0.3,  0.3,  0.3,  0.3};
static const double pcutmax[] = {30.0,  3.0,  3.0, 30.0, 99.9};

static const int nv1Ebins = 28;
static const double v1Ebins[] = {
    -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0,
    -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,
     2.4,  2.8,  3.2,  3.6,  4.0,  4.4,  4.8,  5.2,  5.6};

//-- The following values assume a centrality of 20-30%
static const double v1oddin[] = {
    -0.0373, -0.0324, -0.0281, -0.0232,	-0.0189, -0.0140,
    -0.0097, -0.0049, -0.0006,	0.0043,  0.0054,  0.0049,
     0.0031,  0.0012, -0.0012, -0.0031,	-0.0049, -0.0054,
    -0.0043,  0.0006,  0.0049,  0.0097,	 0.0140,  0.0189,
     0.0232,  0.0281,  0.0324,  0.0373};
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

const std::string HFName[] = { "HFm", "HFp" };
const std::string EPName[] = { "HFm", "trackm", "trackp", "HFp", "trkmid" };
static const Int_t multMax = 8000;


//-- initial MC throw

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
TH1D * hQmixnAA[2][netabins];
TH1D * hQmixnAB[2][netabins];
TH1D * hQmixnAC[2][netabins];
TH1D * hQmixABC[2][netabins];
TH1D * hQmixAABB[2][netabins];
TH1D * hQmixAACC[2][netabins];
TH1D * hQmixBBCC[2][netabins];

TH1D * hQ1nA_final[2];
TH1D * hQ1AB_final[2];
TH1D * hQ1AC_final[2];
TH1D * hQ1BC_final[2];
TH1D * hQ2nA_final[2];
TH1D * hQ2AB_final[2];
TH1D * hQ2AC_final[2];
TH1D * hQ2BC_final[2];
TH1D * hQmixnAA_final[2];
TH1D * hQmixnAB_final[2];
TH1D * hQmixnAC_final[2];
TH1D * hQmixABC_final[2];
TH1D * hQmixAABB_final[2];
TH1D * hQmixAACC_final[2];
TH1D * hQmixBBCC_final[2];

//-- normalized Q-vector correlations

TH1D * hQ1nAnorm[2][netabins];
TH1D * hQ1ABnorm[2][netabins];
TH1D * hQ1ACnorm[2][netabins];
TH1D * hQ1BCnorm[2][netabins];
TH1D * hQ2nAnorm[2][netabins];
TH1D * hQ2ABnorm[2][netabins];
TH1D * hQ2ACnorm[2][netabins];
TH1D * hQ2BCnorm[2][netabins];
TH1D * hQmixnAAnorm[2][netabins];
TH1D * hQmixnABnorm[2][netabins];
TH1D * hQmixnACnorm[2][netabins];
TH1D * hQmixABCnorm[2][netabins];
TH1D * hQmixAABBnorm[2][netabins];
TH1D * hQmixAACCnorm[2][netabins];
TH1D * hQmixBBCCnorm[2][netabins];

TH1D * hQ1nAnorm_final[2];
TH1D * hQ1ABnorm_final[2];
TH1D * hQ1ACnorm_final[2];
TH1D * hQ1BCnorm_final[2];
TH1D * hQ2nAnorm_final[2];
TH1D * hQ2ABnorm_final[2];
TH1D * hQ2ACnorm_final[2];
TH1D * hQ2BCnorm_final[2];
TH1D * hQmixnAAnorm_final[2];
TH1D * hQmixnABnorm_final[2];
TH1D * hQmixnACnorm_final[2];
TH1D * hQmixABCnorm_final[2];
TH1D * hQmixAABBnorm_final[2];
TH1D * hQmixAACCnorm_final[2];
TH1D * hQmixBBCCnorm_final[2];

//-- event plane and scalar-product vn

TH1D * rescor1_2SE;     // rescor1 2 subevent
TH1D * rescor1_3SE[2];  // rescor1 3 subevent
TH1D * rescor1_mix[2];  // rescor1 mixed harmonic
TH1D * rescor2_2SE;     // rescor2 2 subevent
TH1D * rescor2_3SE[2];  // rescor2 3 subevent
TH1D * EPv1obs[2];
TH1D * EPv1obsMix[2];
TH1D * EPv2obs[2];
TH1D * EPv1_2SE[2];
TH1D * EPv1_3SE[2];
TH1D * EPv1_mix[2];
TH1D * EPv2_2SE[2];
TH1D * EPv2_3SE[2];

TH1D * SPdenom1_2SE;
TH1D * SPdenom1_3SE[2];
TH1D * SPdenom1_mix[2];
TH1D * SPdenom2_2SE;
TH1D * SPdenom2_3SE[2];
TH1D * SPv1num[2];
TH1D * SPv1numMix[2];
TH1D * SPv2num[2];
TH1D * SPv1_2SE[2];
TH1D * SPv1_3SE[2];
TH1D * SPv1_mix[2];
TH1D * SPv2_2SE[2];
TH1D * SPv2_3SE[2];

//-- differential vn(pT)

TH1D * EPv1obs_pt[2][nptbins];
TH1D * EPv1obsMix_pt[2][nptbins];
TH1D * EPv2obs_pt[2][nptbins];
TH1D * EPv1_2SE_pt[2][nptbins];
TH1D * EPv1_3SE_pt[2][nptbins];
TH1D * EPv1_mix_pt[2][nptbins];
TH1D * EPv2_2SE_pt[2][nptbins];
TH1D * EPv2_3SE_pt[2][nptbins];

TH1D * SPv1num_pt[2][nptbins];
TH1D * SPv1numMix_pt[2][nptbins];
TH1D * SPv2num_pt[2][nptbins];
TH1D * SPv1_2SE_pt[2][nptbins];
TH1D * SPv1_3SE_pt[2][nptbins];
TH1D * SPv1_mix_pt[2][nptbins];
TH1D * SPv2_2SE_pt[2][nptbins];
TH1D * SPv2_3SE_pt[2][nptbins];

TH1D * hqcnt_pt[2][nptbins];
TH1D * hncnt_pt[2][nptbins];

TH1D * DiffEPv1obs_pt[2];
TH1D * DiffEPv1obsMix_pt[2];
TH1D * DiffEPv2obs_pt[2];
TH1D * DiffEPv1_2SE_pt[2];
TH1D * DiffEPv1_3SE_pt[2];
TH1D * DiffEPv1_mix_pt[2];
TH1D * DiffEPv2_2SE_pt[2];
TH1D * DiffEPv2_3SE_pt[2];

TH1D * DiffSPv1num_pt[2];
TH1D * DiffSPv1numMix_pt[2];
TH1D * DiffSPv2num_pt[2];
TH1D * DiffSPv1_2SE_pt[2];
TH1D * DiffSPv1_3SE_pt[2];
TH1D * DiffSPv1_mix_pt[2];
TH1D * DiffSPv2_2SE_pt[2];
TH1D * DiffSPv2_3SE_pt[2];

TH1D * Diffqcnt_pt[2];
TH1D * Diffncnt_pt[2];

//-- differential vn(eta)

TH1D * EPv1obs_eta[2][netabins];
TH1D * EPv1obsMix_eta[2][netabins];
TH1D * EPv2obs_eta[2][netabins];
TH1D * EPv1_2SE_eta[2][netabins];
TH1D * EPv1_3SE_eta[2][netabins];
TH1D * EPv1_mix_eta[2][netabins];
TH1D * EPv2_2SE_eta[2][netabins];
TH1D * EPv2_3SE_eta[2][netabins];

TH1D * SPv1num_eta[2][netabins];
TH1D * SPv1numMix_eta[2][netabins];
TH1D * SPv2num_eta[2][netabins];
TH1D * SPv1_2SE_eta[2][netabins];
TH1D * SPv1_3SE_eta[2][netabins];
TH1D * SPv1_mix_eta[2][netabins];
TH1D * SPv2_2SE_eta[2][netabins];
TH1D * SPv2_3SE_eta[2][netabins];

TH1D * hqcnt_eta[2][netabins];
TH1D * hncnt_eta[2][netabins];

TH1D * DiffEPv1obs_eta[2];
TH1D * DiffEPv1obsMix_eta[2];
TH1D * DiffEPv2obs_eta[2];
TH1D * DiffEPv1_2SE_eta[2];
TH1D * DiffEPv1_3SE_eta[2];
TH1D * DiffEPv1_mix_eta[2];
TH1D * DiffEPv2_2SE_eta[2];
TH1D * DiffEPv2_3SE_eta[2];

TH1D * DiffSPv1num_eta[2];
TH1D * DiffSPv1numMix_eta[2];
TH1D * DiffSPv2num_eta[2];
TH1D * DiffSPv1_2SE_eta[2];
TH1D * DiffSPv1_3SE_eta[2];
TH1D * DiffSPv1_mix_eta[2];
TH1D * DiffSPv2_2SE_eta[2];
TH1D * DiffSPv2_3SE_eta[2];

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
    hinit_v1in = new TH1D("v1in", "v1in", nv1Ebins, v1Ebins);
    rescor1_2SE = new TH1D("rescor1_2SE", "rescor1_2SE", 200, 0, 1);
    rescor2_2SE = new TH1D("rescor1_3SE", "rescor1_3SE", 200, 0, 1);
    SPdenom1_2SE = new TH1D("SPdenom1_2SE", "SPdenom1_2SE", 200, 0, 1);
    SPdenom2_2SE = new TH1D("SPdenom2_2SE", "SPdenom2_2SE", 200, 0, 1);

    for (int iside = 0; iside<2; iside++) {
        for (int ebin = 0; ebin<netabins; ebin++) {

            hQ1nA[iside][ebin] = new TH1D(Form("Q1nA_%s_%d",HFName[iside].data(),ebin), "Q1nA", 150, -1e-3, 1e-3);
            hQ2nA[iside][ebin] = new TH1D(Form("Q2nA_%s_%d",HFName[iside].data(),ebin), "Q2nA", 150, -5e-3, 5e-3);
            hQmixnAA[iside][ebin] = new TH1D(Form("QmixnAA_%s_%d",HFName[iside].data(),ebin), "QmixnAA", 150, -1e-4, 1e-4);
            hQmixnAB[iside][ebin] = new TH1D(Form("QmixnAB_%s_%d",HFName[iside].data(),ebin), "QmixnAB", 150, -1e-4, 1e-4);
            hQmixnAC[iside][ebin] = new TH1D(Form("QmixnAC_%s_%d",HFName[iside].data(),ebin), "QmixnAC", 150, -1e-4, 1e-4);

            hQ1nAnorm[iside][ebin] = new TH1D(Form("Q1nAnorm_%s_%d",HFName[iside].data(),ebin), "Q1nAnorm", 150, -1e-3, 1e-3);
            hQ2nAnorm[iside][ebin] = new TH1D(Form("Q2nAnorm_%s_%d",HFName[iside].data(),ebin), "Q2nAnorm", 150, -5e-3, 5e-3);
            hQmixnAAnorm[iside][ebin] = new TH1D(Form("QmixnAAnorm_%s_%d",HFName[iside].data(),ebin), "QmixnAAnorm", 150, -1e-4, 1e-4);
            hQmixnABnorm[iside][ebin] = new TH1D(Form("QmixnABnorm_%s_%d",HFName[iside].data(),ebin), "QmixnABnorm", 150, -1e-4, 1e-4);
            hQmixnACnorm[iside][ebin] = new TH1D(Form("QmixnACnorm_%s_%d",HFName[iside].data(),ebin), "QmixnACnorm", 150, -1e-4, 1e-4);
        }

        hQ1nA_final[iside] = new TH1D(Form("Q1nA_%s",HFName[iside].data()), "Q1nA", netabins, etabins);
        hQ1AB_final[iside] = new TH1D(Form("Q1AB_%s",HFName[iside].data()), "Q1AB", 200, -1e-3, 1e-3);
        hQ1AC_final[iside] = new TH1D(Form("Q1AC_%s",HFName[iside].data()), "Q1AC", 200, -1e-3, 1e-3);
        hQ1BC_final[iside] = new TH1D(Form("Q1BC_%s",HFName[iside].data()), "Q1BC", 200, -1e-3, 1e-3);
        hQ2nA_final[iside] = new TH1D(Form("Q2nA_%s",HFName[iside].data()), "Q2nA", netabins, etabins);
        hQ2AB_final[iside] = new TH1D(Form("Q2AB_%s",HFName[iside].data()), "Q2AB", 200, -1e-2, 1e-2);
        hQ2AC_final[iside] = new TH1D(Form("Q2AC_%s",HFName[iside].data()), "Q2AC", 200, -1e-2, 1e-2);
        hQ2BC_final[iside] = new TH1D(Form("Q2BC_%s",HFName[iside].data()), "Q2BC", 200, -1e-2, 1e-2);
        hQmixnAA_final[iside] = new TH1D(Form("QmixnAA_%s",HFName[iside].data()), "QmixnAA", netabins, etabins);
        hQmixnAB_final[iside] = new TH1D(Form("QmixnAB_%s",HFName[iside].data()), "QmixnAB", netabins, etabins);
        hQmixnAC_final[iside] = new TH1D(Form("QmixnAC_%s",HFName[iside].data()), "QmixnAC", netabins, etabins);
        hQmixABC_final[iside] = new TH1D(Form("QmixABC_%s",HFName[iside].data()), "QmixABC", 200, -1e-4, 1e-4);
        hQmixAABB_final[iside] = new TH1D(Form("QmixAABB_%s",HFName[iside].data()), "QmixAABB", 200, -1e-5, 1e-5);
        hQmixAACC_final[iside] = new TH1D(Form("QmixAACC_%s",HFName[iside].data()), "QmixAACC", 200, -1e-5, 1e-5);
        hQmixBBCC_final[iside] = new TH1D(Form("QmixBBCC_%s",HFName[iside].data()), "QmixBBCC", 200, -1e-5, 1e-5);

        hQ1nAnorm_final[iside] = new TH1D(Form("Q1nAnorm_%s",HFName[iside].data()), "Q1nAnorm", netabins, etabins);
        hQ1ABnorm_final[iside] = new TH1D(Form("Q1ABnorm_%s",HFName[iside].data()), "Q1ABnorm", 200, -0.5, 0.5);
        hQ1ACnorm_final[iside] = new TH1D(Form("Q1ACnorm_%s",HFName[iside].data()), "Q1ACnorm", 200, -0.5, 0.5);
        hQ1BCnorm_final[iside] = new TH1D(Form("Q1BCnorm_%s",HFName[iside].data()), "Q1BCnorm", 200, -0.5, 0.5);
        hQ2nAnorm_final[iside] = new TH1D(Form("Q2nAnorm_%s",HFName[iside].data()), "Q2nAnorm", netabins, etabins);
        hQ2ABnorm_final[iside] = new TH1D(Form("Q2ABnorm_%s",HFName[iside].data()), "Q2ABnorm", 200, -1, 1);
        hQ2ACnorm_final[iside] = new TH1D(Form("Q2ACnorm_%s",HFName[iside].data()), "Q2ACnorm", 200, -1, 1);
        hQ2BCnorm_final[iside] = new TH1D(Form("Q2BCnorm_%s",HFName[iside].data()), "Q2BCnorm", 200, -1, 1);
        hQmixnAAnorm_final[iside] = new TH1D(Form("QmixnAAnorm_%s",HFName[iside].data()), "QmixnAAnorm", netabins, etabins);
        hQmixnABnorm_final[iside] = new TH1D(Form("QmixnABnorm_%s",HFName[iside].data()), "QmixnABnorm", netabins, etabins);
        hQmixnACnorm_final[iside] = new TH1D(Form("QmixnACnorm_%s",HFName[iside].data()), "QmixnACnorm", netabins, etabins);
        hQmixABCnorm_final[iside] = new TH1D(Form("QmixABCnorm_%s",HFName[iside].data()), "QmixABCnorm", netabins, etabins);
        hQmixAABBnorm_final[iside] = new TH1D(Form("QmixAABBnorm_%s",HFName[iside].data()), "QmixAABBnorm", 200, -1, 1);
        hQmixAACCnorm_final[iside] = new TH1D(Form("QmixAACCnorm_%s",HFName[iside].data()), "QmixAACCnorm", 200, -1, 1);
        hQmixBBCCnorm_final[iside] = new TH1D(Form("QmixBBCCnorm_%s",HFName[iside].data()), "QmixBBCCnorm", 200, -1, 1);

        rescor1_3SE[iside] = new TH1D(Form("rescor1_3SE_%s",HFName[iside].data()), "rescor1_3SE", 200, 0, 1);
        rescor1_mix[iside] = new TH1D(Form("rescor1_mix_%s",HFName[iside].data()), "rescor1_mix", 200, 0, 1);
        rescor2_3SE[iside] = new TH1D(Form("rescor2_3SE_%s",HFName[iside].data()), "rescor2_3SE", 200, 0, 1);
        EPv1obs[iside] = new TH1D(Form("EPv1obs_%s",HFName[iside].data()), "EPv1obs", 200, -0.1, 0.1);
        EPv1obsMix[iside] = new TH1D(Form("EPv1obsMix_%s",HFName[iside].data()), "EPv1obsMix", 200, -0.1, 0.1);
        EPv2obs[iside] = new TH1D(Form("EPv2obs_%s",HFName[iside].data()), "EPv2obs", 200, -0.1, 0.1);
        EPv1_2SE[iside] = new TH1D(Form("EPv1_2SE_%s",HFName[iside].data()), "EPv1_2SE", 200, -0.15, 0.15);
        EPv1_3SE[iside] = new TH1D(Form("EPv1_3SE_%s",HFName[iside].data()), "EPv1_3SE", 200, -0.15, 0.15);
        EPv1_mix[iside] = new TH1D(Form("EPv1_mix_%s",HFName[iside].data()), "EPv1_mix", 200, -0.15, 0.15);
        EPv2_2SE[iside] = new TH1D(Form("EPv2_2SE_%s",HFName[iside].data()), "EPv2_2SE", 200, -0.15, 0.15);
        EPv2_3SE[iside] = new TH1D(Form("EPv2_3SE_%s",HFName[iside].data()), "EPv2_3SE", 200, -0.15, 0.15);

        SPdenom1_3SE[iside] = new TH1D(Form("SPdenom1_3SE_%s",HFName[iside].data()), "SPdenom1_3SE", 200, -0.1, 0.1);
        SPdenom1_mix[iside] = new TH1D(Form("SPdenom1_mix_%s",HFName[iside].data()), "SPdenom1_mix", 200, -0.1, 0.1);
        SPdenom2_3SE[iside] = new TH1D(Form("SPdenom2_3SE_%s",HFName[iside].data()), "SPdenom2_3SE", 200, -0.1, 0.1);
        SPv1num[iside] = new TH1D(Form("SPv1num_%s",HFName[iside].data()), "SPv1num", 200, -0.1, 0.1);
        SPv1numMix[iside] = new TH1D(Form("SPv1numMix_%s",HFName[iside].data()), "SPv1numMix", 200, -0.1, 0.1);
        SPv2num[iside] = new TH1D(Form("SPv2num_%s",HFName[iside].data()), "SPv2num", 200, -0.1, 0.1);
        SPv1_2SE[iside] = new TH1D(Form("SPv1_2SE_%s",HFName[iside].data()), "SPv1_2SE", 200, -0.1, 0.1);
        SPv1_3SE[iside] = new TH1D(Form("SPv1_3SE_%s",HFName[iside].data()), "SPv1_3SE", 200, -0.1, 0.1);
        SPv1_mix[iside] = new TH1D(Form("SPv1_mix_%s",HFName[iside].data()), "SPv1_mix", 200, -0.1, 0.1);
        SPv2_2SE[iside] = new TH1D(Form("SPv2_2SE_%s",HFName[iside].data()), "SPv2_2SE", 200, -0.1, 0.1);
        SPv2_3SE[iside] = new TH1D(Form("SPv2_3SE_%s",HFName[iside].data()), "SPv2_3SE", 200, -0.1, 0.1);

        for (int pbin = 0; pbin<nptbins; pbin++) {
            EPv1obs_pt[iside][pbin] = new TH1D(Form("EPv1obs_pt_%s_%d",HFName[iside].data(),pbin), "EPv1obs_pt", 200, -0.1, 0.1);
            EPv1obsMix_pt[iside][pbin] = new TH1D(Form("EPv1obsMix_pt_%s_%d",HFName[iside].data(),pbin), "EPv1obsMix_pt", 200, -0.1, 0.1);
            EPv2obs_pt[iside][pbin] = new TH1D(Form("EPv2obs_pt_%s_%d",HFName[iside].data(),pbin), "EPv2obs_pt", 200, -0.1, 0.1);
            EPv1_2SE_pt[iside][pbin] = new TH1D(Form("EPv1_2SE_pt_%s_%d",HFName[iside].data(),pbin), "EPv1_2SE_pt", 200, -0.15, 0.15);
            EPv1_3SE_pt[iside][pbin] = new TH1D(Form("EPv1_3SE_pt_%s_%d",HFName[iside].data(),pbin), "EPv1_3SE_pt", 200, -0.15, 0.15);
            EPv1_mix_pt[iside][pbin] = new TH1D(Form("EPv1_mix_pt_%s_%d",HFName[iside].data(),pbin), "EPv1_mix_pt", 200, -0.15, 0.15);
            EPv2_2SE_pt[iside][pbin] = new TH1D(Form("EPv2_2SE_pt_%s_%d",HFName[iside].data(),pbin), "EPv2_2SE_pt", 200, -0.15, 0.15);
            EPv2_3SE_pt[iside][pbin] = new TH1D(Form("EPv2_3SE_pt_%s_%d",HFName[iside].data(),pbin), "EPv2_3SE_pt", 200, -0.15, 0.15);

            SPv1num_pt[iside][pbin] = new TH1D(Form("SPv1num_pt_%s_%d",HFName[iside].data(),pbin), "SPv1num_pt", 200, -0.1, 0.1);
            SPv1numMix_pt[iside][pbin] = new TH1D(Form("SPv1numMix_pt_%s_%d",HFName[iside].data(),pbin), "SPv1numMix_pt", 200, -0.1, 0.1);
            SPv2num_pt[iside][pbin] = new TH1D(Form("SPv2num_pt_%s_%d",HFName[iside].data(),pbin), "SPv2num_pt", 200, -0.1, 0.1);
            SPv1_2SE_pt[iside][pbin] = new TH1D(Form("SPv1_2SE_pt_%s_%d",HFName[iside].data(),pbin), "SPv1_2SE_pt", 200, -0.15, 0.15);
            SPv1_3SE_pt[iside][pbin] = new TH1D(Form("SPv1_3SE_pt_%s_%d",HFName[iside].data(),pbin), "SPv1_3SE_pt", 200, -0.15, 0.15);
            SPv1_mix_pt[iside][pbin] = new TH1D(Form("SPv1_mix_pt_%s_%d",HFName[iside].data(),pbin), "SPv1_mix_pt", 200, -0.15, 0.15);
            SPv2_2SE_pt[iside][pbin] = new TH1D(Form("SPv2_2SE_pt_%s_%d",HFName[iside].data(),pbin), "SPv2_2SE_pt", 200, -0.15, 0.15);
            SPv2_3SE_pt[iside][pbin] = new TH1D(Form("SPv2_3SE_pt_%s_%d",HFName[iside].data(),pbin), "SPv2_3SE_pt", 200, -0.15, 0.15);

            hqcnt_pt[iside][pbin] = new TH1D(Form("qcnt_pt_%s_%d",HFName[iside].data(),pbin), "qcnt_pt", 1000, 0, 1e6);
            hncnt_pt[iside][pbin] = new TH1D(Form("ncnt_pt_%s_%d",HFName[iside].data(),pbin), "ncnt_pt", 1000, 0, 1e9);
        }

        DiffEPv1obs_pt[iside] = new TH1D(Form("EPv1obs_pt_%s",HFName[iside].data()), "EPv1obs_pt", nptbins, ptbins);
        DiffEPv1obsMix_pt[iside] = new TH1D(Form("EPv1obsMix_pt_%s",HFName[iside].data()), "EPv1obsMix_pt", nptbins, ptbins);
        DiffEPv2obs_pt[iside] = new TH1D(Form("EPv2obs_pt_%s",HFName[iside].data()), "EPv2obs_pt", nptbins, ptbins);
        DiffEPv1_2SE_pt[iside] = new TH1D(Form("EPv1_2SE_pt_%s",HFName[iside].data()), "EPv1_2SE_pt", nptbins, ptbins);
        DiffEPv1_3SE_pt[iside] = new TH1D(Form("EPv1_3SE_pt_%s",HFName[iside].data()), "EPv1_3SE_pt", nptbins, ptbins);
        DiffEPv1_mix_pt[iside] = new TH1D(Form("EPv1_mix_pt_%s",HFName[iside].data()), "EPv1_mix_pt", nptbins, ptbins);
        DiffEPv2_2SE_pt[iside] = new TH1D(Form("EPv2_2SE_pt_%s",HFName[iside].data()), "EPv2_2SE_pt", nptbins, ptbins);
        DiffEPv2_3SE_pt[iside] = new TH1D(Form("EPv2_3SE_pt_%s",HFName[iside].data()), "EPv2_3SE_pt", nptbins, ptbins);

        DiffSPv1num_pt[iside] = new TH1D(Form("SPv1num_pt_%s",HFName[iside].data()), "SPv1num_pt", nptbins, ptbins);
        DiffSPv1numMix_pt[iside] = new TH1D(Form("SPv1numMix_pt_%s",HFName[iside].data()), "SPv1numMix_pt", nptbins, ptbins);
        DiffSPv2num_pt[iside] = new TH1D(Form("SPv2num_pt_%s",HFName[iside].data()), "SPv2num_pt", nptbins, ptbins);
        DiffSPv1_2SE_pt[iside] = new TH1D(Form("SPv1_2SE_pt_%s",HFName[iside].data()), "SPv1_2SE_pt", nptbins, ptbins);
        DiffSPv1_3SE_pt[iside] = new TH1D(Form("SPv1_3SE_pt_%s",HFName[iside].data()), "SPv1_3SE_pt", nptbins, ptbins);
        DiffSPv1_mix_pt[iside] = new TH1D(Form("SPv1_mix_pt_%s",HFName[iside].data()), "SPv1_mix_pt", nptbins, ptbins);
        DiffSPv2_2SE_pt[iside] = new TH1D(Form("SPv2_2SE_pt_%s",HFName[iside].data()), "SPv2_2SE_pt", nptbins, ptbins);
        DiffSPv2_3SE_pt[iside] = new TH1D(Form("SPv2_3SE_pt_%s",HFName[iside].data()), "SPv2_3SE_pt", nptbins, ptbins);

        Diffqcnt_pt[iside] = new TH1D(Form("Qcnt_pt_%s",HFName[iside].data()), "Qcnt_pt", nptbins, ptbins);
        Diffncnt_pt[iside] = new TH1D(Form("Ncnt_pt_%s",HFName[iside].data()), "Ncnt_pt", nptbins, ptbins);

        for (int ebin = 0; ebin<netabins; ebin++) {
            EPv1obs_eta[iside][ebin] = new TH1D(Form("EPv1obs_eta_%s_%d",HFName[iside].data(),ebin), "EPv1obs_eta", 200, -0.1, 0.1);
            EPv1obsMix_eta[iside][ebin] = new TH1D(Form("EPv1obsMix_eta_%s_%d",HFName[iside].data(),ebin), "EPv1obsMix_eta", 200, -0.1, 0.1);
            EPv2obs_eta[iside][ebin] = new TH1D(Form("EPv2obs_eta_%s_%d",HFName[iside].data(),ebin), "EPv2obs_eta", 200, -0.1, 0.1);
            EPv1_2SE_eta[iside][ebin] = new TH1D(Form("EPv1_2SE_eta_%s_%d",HFName[iside].data(),ebin), "EPv1_2SE_eta", 200, -0.15, 0.15);
            EPv1_3SE_eta[iside][ebin] = new TH1D(Form("EPv1_3SE_eta_%s_%d",HFName[iside].data(),ebin), "EPv1_3SE_eta", 200, -0.15, 0.15);
            EPv1_mix_eta[iside][ebin] = new TH1D(Form("EPv1_mix_eta_%s_%d",HFName[iside].data(),ebin), "EPv1_mix_eta", 200, -0.15, 0.15);
            EPv2_2SE_eta[iside][ebin] = new TH1D(Form("EPv2_2SE_eta_%s_%d",HFName[iside].data(),ebin), "EPv2_2SE_eta", 200, -0.15, 0.15);
            EPv2_3SE_eta[iside][ebin] = new TH1D(Form("EPv2_3SE_eta_%s_%d",HFName[iside].data(),ebin), "EPv2_3SE_eta", 200, -0.15, 0.15);

            SPv1num_eta[iside][ebin] = new TH1D(Form("SPv1num_eta_%s_%d",HFName[iside].data(),ebin), "SPv1num_eta", 200, -0.1, 0.1);
            SPv1numMix_eta[iside][ebin] = new TH1D(Form("SPv1numMix_eta_%s_%d",HFName[iside].data(),ebin), "SPv1numMix_eta", 200, -0.1, 0.1);
            SPv2num_eta[iside][ebin] = new TH1D(Form("SPv2num_eta_%s_%d",HFName[iside].data(),ebin), "SPv2num_eta", 200, -0.1, 0.1);
            SPv1_2SE_eta[iside][ebin] = new TH1D(Form("SPv1_2SE_eta_%s_%d",HFName[iside].data(),ebin), "SPv1_2SE_eta", 200, -0.15, 0.15);
            SPv1_3SE_eta[iside][ebin] = new TH1D(Form("SPv1_3SE_eta_%s_%d",HFName[iside].data(),ebin), "SPv1_3SE_eta", 200, -0.15, 0.15);
            SPv1_mix_eta[iside][ebin] = new TH1D(Form("SPv1_mix_eta_%s_%d",HFName[iside].data(),ebin), "SPv1_mix_eta", 200, -0.15, 0.15);
            SPv2_2SE_eta[iside][ebin] = new TH1D(Form("SPv2_2SE_eta_%s_%d",HFName[iside].data(),ebin), "SPv2_2SE_eta", 200, -0.15, 0.15);
            SPv2_3SE_eta[iside][ebin] = new TH1D(Form("SPv2_3SE_eta_%s_%d",HFName[iside].data(),ebin), "SPv2_3SE_eta", 200, -0.15, 0.15);

            hqcnt_eta[iside][ebin] = new TH1D(Form("qcnt_eta_%s_%d",HFName[iside].data(),ebin), "qcnt_eta", 1000, 0, 1e6);
            hncnt_eta[iside][ebin] = new TH1D(Form("ncnt_eta_%s_%d",HFName[iside].data(),ebin), "ncnt_eta", 1000, 0, 1e9);
        }

        DiffEPv1obs_eta[iside] = new TH1D(Form("EPv1obs_eta_%s",HFName[iside].data()), "EPv1obs_eta", netabins, etabins);
        DiffEPv1obsMix_eta[iside] = new TH1D(Form("EPv1obsMix_eta_%s",HFName[iside].data()), "EPv1obsMix_eta", netabins, etabins);
        DiffEPv2obs_eta[iside] = new TH1D(Form("EPv2obs_eta_%s",HFName[iside].data()), "EPv2obs_eta", netabins, etabins);
        DiffEPv1_2SE_eta[iside] = new TH1D(Form("EPv1_2SE_eta_%s",HFName[iside].data()), "EPv1_2SE_eta", netabins, etabins);
        DiffEPv1_3SE_eta[iside] = new TH1D(Form("EPv1_3SE_eta_%s",HFName[iside].data()), "EPv1_3SE_eta", netabins, etabins);
        DiffEPv1_mix_eta[iside] = new TH1D(Form("EPv1_mix_eta_%s",HFName[iside].data()), "EPv1_mix_eta", netabins, etabins);
        DiffEPv2_2SE_eta[iside] = new TH1D(Form("EPv2_2SE_eta_%s",HFName[iside].data()), "EPv2_2SE_eta", netabins, etabins);
        DiffEPv2_3SE_eta[iside] = new TH1D(Form("EPv2_3SE_eta_%s",HFName[iside].data()), "EPv2_3SE_eta", netabins, etabins);

        DiffSPv1num_eta[iside] = new TH1D(Form("SPv1num_eta_%s",HFName[iside].data()), "SPv1num_eta", netabins, etabins);
        DiffSPv1numMix_eta[iside] = new TH1D(Form("SPv1numMix_eta_%s",HFName[iside].data()), "SPv1numMix_eta", netabins, etabins);
        DiffSPv2num_eta[iside] = new TH1D(Form("SPv2num_eta_%s",HFName[iside].data()), "SPv2num_eta", netabins, etabins);
        DiffSPv1_2SE_eta[iside] = new TH1D(Form("SPv1_2SE_eta_%s",HFName[iside].data()), "SPv1_2SE_eta", netabins, etabins);
        DiffSPv1_3SE_eta[iside] = new TH1D(Form("SPv1_3SE_eta_%s",HFName[iside].data()), "SPv1_3SE_eta", netabins, etabins);
        DiffSPv1_mix_eta[iside] = new TH1D(Form("SPv1_mix_eta_%s",HFName[iside].data()), "SPv1_mix_eta", netabins, etabins);
        DiffSPv2_2SE_eta[iside] = new TH1D(Form("SPv2_2SE_eta_%s",HFName[iside].data()), "SPv2_2SE_eta", netabins, etabins);
        DiffSPv2_3SE_eta[iside] = new TH1D(Form("SPv2_3SE_eta_%s",HFName[iside].data()), "SPv2_3SE_eta", netabins, etabins);

        Diffqcnt_eta[iside] = new TH1D(Form("Qcnt_eta_%s",HFName[iside].data()), "Qcnt_eta", netabins, etabins);
        Diffncnt_eta[iside] = new TH1D(Form("Ncnt_eta_%s",HFName[iside].data()), "Ncnt_eta", netabins, etabins);
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
    tfout = new TFile(Form("results/results_fixPsip2_%s22.root",mtag.Data()),"recreate");

    TDirectory * tdinput = (TDirectory *) tfout->mkdir("Inputs");
    tdinput->cd();
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
    for (int nep = 0; nep<numEP; nep++) hsub_Psi1[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi2[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi1lab[nep]->Write();
    for (int nep = 0; nep<numEP; nep++) hsub_Psi2lab[nep]->Write();
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
        tdqvec_side[iside] = (TDirectory *) tdqvec->mkdir(Form("%s",HFName[iside].data()));
        tdqvec_side[iside]->cd();
        hQ1nA_final[iside]->Write();
        hQ1AB_final[iside]->Write();
        hQ1AC_final[iside]->Write();
        hQ1BC_final[iside]->Write();
        hQ2nA_final[iside]->Write();
        hQ2AB_final[iside]->Write();
        hQ2AC_final[iside]->Write();
        hQ2BC_final[iside]->Write();
        hQmixnAA_final[iside]->Write();
        hQmixnAB_final[iside]->Write();
        hQmixnAC_final[iside]->Write();
        hQmixABC_final[iside]->Write();
        hQmixAABB_final[iside]->Write();
        hQmixAACC_final[iside]->Write();
        hQmixBBCC_final[iside]->Write();

        tdqvecNorm_side[iside] = (TDirectory *) tdqvecNorm->mkdir(Form("%s",HFName[iside].data()));
        tdqvecNorm_side[iside]->cd();
        hQ1nAnorm_final[iside]->Write();
        hQ1ABnorm_final[iside]->Write();
        hQ1ACnorm_final[iside]->Write();
        hQ1BCnorm_final[iside]->Write();
        hQ2nAnorm_final[iside]->Write();
        hQ2ABnorm_final[iside]->Write();
        hQ2ACnorm_final[iside]->Write();
        hQ2BCnorm_final[iside]->Write();
        hQmixnAAnorm_final[iside]->Write();
        hQmixnABnorm_final[iside]->Write();
        hQmixnACnorm_final[iside]->Write();
        hQmixABCnorm_final[iside]->Write();
        hQmixAABBnorm_final[iside]->Write();
        hQmixAACCnorm_final[iside]->Write();
        hQmixBBCCnorm_final[iside]->Write();

        tdvnep_side[iside] = (TDirectory *) tdvnep->mkdir(Form("%s",HFName[iside].data()));
        tdvnep_side[iside]->cd();
        rescor1_3SE[iside]->Write();
        rescor1_mix[iside]->Write();
        rescor2_3SE[iside]->Write();
        EPv1obs[iside]->Write();
        EPv1obsMix[iside]->Write();
        EPv2obs[iside]->Write();
        EPv1_2SE[iside]->Write();
        EPv1_3SE[iside]->Write();
        EPv1_mix[iside]->Write();
        EPv2_2SE[iside]->Write();
        EPv2_3SE[iside]->Write();

        tdvnsp_side[iside] = (TDirectory *) tdvnsp->mkdir(Form("%s",HFName[iside].data()));
        tdvnsp_side[iside]->cd();
        SPdenom1_3SE[iside]->Write();
        SPdenom1_mix[iside]->Write();
        SPdenom2_3SE[iside]->Write();
        SPv1num[iside]->Write();
        SPv1numMix[iside]->Write();
        SPv2num[iside]->Write();
        SPv1_2SE[iside]->Write();
        SPv1_3SE[iside]->Write();
        SPv1_mix[iside]->Write();
        SPv2_2SE[iside]->Write();
        SPv2_3SE[iside]->Write();

        tdvnDiff_pt_side[iside] = (TDirectory *) tdvnDiff_pt->mkdir(Form("%s",HFName[iside].data()));
        tdvnDiff_pt_side[iside]->cd();
        DiffEPv1_2SE_pt[iside]->Write();
        DiffEPv1_3SE_pt[iside]->Write();
        DiffEPv1_mix_pt[iside]->Write();
        DiffEPv2_2SE_pt[iside]->Write();
        DiffEPv2_3SE_pt[iside]->Write();
        DiffSPv1_2SE_pt[iside]->Write();
        DiffSPv1_3SE_pt[iside]->Write();
        DiffSPv1_mix_pt[iside]->Write();
        DiffSPv2_2SE_pt[iside]->Write();
        DiffSPv2_3SE_pt[iside]->Write();
        DiffEPv1obs_pt[iside]->Write();
        DiffEPv1obsMix_pt[iside]->Write();
        DiffEPv2obs_pt[iside]->Write();
        DiffSPv1num_pt[iside]->Write();
        DiffSPv1numMix_pt[iside]->Write();
        DiffSPv2num_pt[iside]->Write();
        Diffqcnt_pt[iside]->Write();
        Diffncnt_pt[iside]->Write();

        tdvnDiff_eta_side[iside] = (TDirectory *) tdvnDiff_eta->mkdir(Form("%s",HFName[iside].data()));
        tdvnDiff_eta_side[iside]->cd();
        DiffEPv1_2SE_eta[iside]->Write();
        DiffEPv1_3SE_eta[iside]->Write();
        DiffEPv1_mix_eta[iside]->Write();
        DiffEPv2_2SE_eta[iside]->Write();
        DiffEPv2_3SE_eta[iside]->Write();
        DiffSPv1_2SE_eta[iside]->Write();
        DiffSPv1_3SE_eta[iside]->Write();
        DiffSPv1_mix_eta[iside]->Write();
        DiffSPv2_2SE_eta[iside]->Write();
        DiffSPv2_3SE_eta[iside]->Write();
        DiffEPv1obs_eta[iside]->Write();
        DiffEPv1obsMix_eta[iside]->Write();
        DiffEPv2obs_eta[iside]->Write();
        DiffSPv1num_eta[iside]->Write();
        DiffSPv1numMix_eta[iside]->Write();
        DiffSPv2num_eta[iside]->Write();
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


void ComputeVN( Int_t nevents, Int_t evtmult, Bool_t isodd, Double_t setv1, Double_t setv2, Bool_t eta_weights, Bool_t pt_weights, Bool_t conserve_pT, Bool_t addholes, Bool_t flatten, Bool_t recenter, Int_t iseed, TString mtag, Int_t ntries ) {

    double v1in = setv1;
    double v2in = setv2;

    for (int iside = 0; iside<2; iside++) {
        for (int pbin = 0; pbin<nptbins; pbin++) {
            DiffEPv1obs_pt[iside]->SetBinContent(pbin+1, EPv1obs_pt[iside][pbin]->GetMean());
            DiffEPv1obs_pt[iside]->SetBinError(pbin+1, EPv1obs_pt[iside][pbin]->GetMeanError());
            DiffEPv1obsMix_pt[iside]->SetBinContent(pbin+1, EPv1obsMix_pt[iside][pbin]->GetMean());
            DiffEPv1obsMix_pt[iside]->SetBinError(pbin+1, EPv1obsMix_pt[iside][pbin]->GetMeanError());
            DiffEPv2obs_pt[iside]->SetBinContent(pbin+1, EPv2obs_pt[iside][pbin]->GetMean());
            DiffEPv2obs_pt[iside]->SetBinError(pbin+1, EPv2obs_pt[iside][pbin]->GetMeanError());
            DiffSPv1num_pt[iside]->SetBinContent(pbin+1, SPv1num_pt[iside][pbin]->GetMean());
            DiffSPv1num_pt[iside]->SetBinError(pbin+1, SPv1num_pt[iside][pbin]->GetMeanError());
            DiffSPv1numMix_pt[iside]->SetBinContent(pbin+1, SPv1numMix_pt[iside][pbin]->GetMean());
            DiffSPv1numMix_pt[iside]->SetBinError(pbin+1, SPv1numMix_pt[iside][pbin]->GetMeanError());
            DiffSPv2num_pt[iside]->SetBinContent(pbin+1, SPv2num_pt[iside][pbin]->GetMean());
            DiffSPv2num_pt[iside]->SetBinError(pbin+1, SPv2num_pt[iside][pbin]->GetMeanError());
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
            DiffEPv2_2SE_pt[iside]->SetBinContent(pbin+1, EPv2_2SE_pt[iside][pbin]->GetMean());
            DiffEPv2_2SE_pt[iside]->SetBinError(pbin+1, EPv2_2SE_pt[iside][pbin]->GetMeanError());
            DiffEPv2_3SE_pt[iside]->SetBinContent(pbin+1, EPv2_3SE_pt[iside][pbin]->GetMean());
            DiffEPv2_3SE_pt[iside]->SetBinError(pbin+1, EPv2_3SE_pt[iside][pbin]->GetMeanError());

            DiffSPv1_2SE_pt[iside]->SetBinContent(pbin+1, SPv1_2SE_pt[iside][pbin]->GetMean());
            DiffSPv1_2SE_pt[iside]->SetBinError(pbin+1, SPv1_2SE_pt[iside][pbin]->GetMeanError());
            DiffSPv1_3SE_pt[iside]->SetBinContent(pbin+1, SPv1_3SE_pt[iside][pbin]->GetMean());
            DiffSPv1_3SE_pt[iside]->SetBinError(pbin+1, SPv1_3SE_pt[iside][pbin]->GetMeanError());
            DiffSPv1_mix_pt[iside]->SetBinContent(pbin+1, SPv1_mix_pt[iside][pbin]->GetMean());
            DiffSPv1_mix_pt[iside]->SetBinError(pbin+1, SPv1_mix_pt[iside][pbin]->GetMeanError());
            DiffSPv2_2SE_pt[iside]->SetBinContent(pbin+1, SPv2_2SE_pt[iside][pbin]->GetMean());
            DiffSPv2_2SE_pt[iside]->SetBinError(pbin+1, SPv2_2SE_pt[iside][pbin]->GetMeanError());
            DiffSPv2_3SE_pt[iside]->SetBinContent(pbin+1, SPv2_3SE_pt[iside][pbin]->GetMean());
            DiffSPv2_3SE_pt[iside]->SetBinError(pbin+1, SPv2_3SE_pt[iside][pbin]->GetMeanError());
        }


        for (int ebin = 0; ebin<netabins; ebin++) {
            hQ1nA_final[iside]->SetBinContent(ebin+1, hQ1nA[iside][ebin]->GetMean());
            hQ1nA_final[iside]->SetBinError(ebin+1, hQ1nA[iside][ebin]->GetMeanError());
            hQ2nA_final[iside]->SetBinContent(ebin+1, hQ2nA[iside][ebin]->GetMean());
            hQ2nA_final[iside]->SetBinError(ebin+1, hQ2nA[iside][ebin]->GetMeanError());
            hQmixnAA_final[iside]->SetBinContent(ebin+1, hQmixnAA[iside][ebin]->GetMean());
            hQmixnAA_final[iside]->SetBinError(ebin+1, hQmixnAA[iside][ebin]->GetMeanError());
            hQmixnAB_final[iside]->SetBinContent(ebin+1, hQmixnAB[iside][ebin]->GetMean());
            hQmixnAB_final[iside]->SetBinError(ebin+1, hQmixnAB[iside][ebin]->GetMeanError());
            hQmixnAC_final[iside]->SetBinContent(ebin+1, hQmixnAC[iside][ebin]->GetMean());
            hQmixnAC_final[iside]->SetBinError(ebin+1, hQmixnAC[iside][ebin]->GetMeanError());

            hQ1nAnorm_final[iside]->SetBinContent(ebin+1, hQ1nAnorm[iside][ebin]->GetMean());
            hQ1nAnorm_final[iside]->SetBinError(ebin+1, hQ1nAnorm[iside][ebin]->GetMeanError());
            hQ2nAnorm_final[iside]->SetBinContent(ebin+1, hQ2nAnorm[iside][ebin]->GetMean());
            hQ2nAnorm_final[iside]->SetBinError(ebin+1, hQ2nAnorm[iside][ebin]->GetMeanError());
            hQmixnAAnorm_final[iside]->SetBinContent(ebin+1, hQmixnAAnorm[iside][ebin]->GetMean());
            hQmixnAAnorm_final[iside]->SetBinError(ebin+1, hQmixnAAnorm[iside][ebin]->GetMeanError());
            hQmixnABnorm_final[iside]->SetBinContent(ebin+1, hQmixnABnorm[iside][ebin]->GetMean());
            hQmixnABnorm_final[iside]->SetBinError(ebin+1, hQmixnABnorm[iside][ebin]->GetMeanError());
            hQmixnACnorm_final[iside]->SetBinContent(ebin+1, hQmixnACnorm[iside][ebin]->GetMean());
            hQmixnACnorm_final[iside]->SetBinError(ebin+1, hQmixnACnorm[iside][ebin]->GetMeanError());


            DiffEPv1obs_eta[iside]->SetBinContent(ebin+1, EPv1obs_eta[iside][ebin]->GetMean());
            DiffEPv1obs_eta[iside]->SetBinError(ebin+1, EPv1obs_eta[iside][ebin]->GetMeanError());
            DiffEPv1obsMix_eta[iside]->SetBinContent(ebin+1, EPv1obsMix_eta[iside][ebin]->GetMean());
            DiffEPv1obsMix_eta[iside]->SetBinError(ebin+1, EPv1obsMix_eta[iside][ebin]->GetMeanError());
            DiffEPv2obs_eta[iside]->SetBinContent(ebin+1, EPv2obs_eta[iside][ebin]->GetMean());
            DiffEPv2obs_eta[iside]->SetBinError(ebin+1, EPv2obs_eta[iside][ebin]->GetMeanError());
            DiffSPv1num_eta[iside]->SetBinContent(ebin+1, SPv1num_eta[iside][ebin]->GetMean());
            DiffSPv1num_eta[iside]->SetBinError(ebin+1, SPv1num_eta[iside][ebin]->GetMeanError());
            DiffSPv1numMix_eta[iside]->SetBinContent(ebin+1, SPv1numMix_eta[iside][ebin]->GetMean());
            DiffSPv1numMix_eta[iside]->SetBinError(ebin+1, SPv1numMix_eta[iside][ebin]->GetMeanError());
            DiffSPv2num_eta[iside]->SetBinContent(ebin+1, SPv2num_eta[iside][ebin]->GetMean());
            DiffSPv2num_eta[iside]->SetBinError(ebin+1, SPv2num_eta[iside][ebin]->GetMeanError());
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
            DiffEPv2_2SE_eta[iside]->SetBinContent(ebin+1, EPv2_2SE_eta[iside][ebin]->GetMean());
            DiffEPv2_2SE_eta[iside]->SetBinError(ebin+1, EPv2_2SE_eta[iside][ebin]->GetMeanError());
            DiffEPv2_3SE_eta[iside]->SetBinContent(ebin+1, EPv2_3SE_eta[iside][ebin]->GetMean());
            DiffEPv2_3SE_eta[iside]->SetBinError(ebin+1, EPv2_3SE_eta[iside][ebin]->GetMeanError());

            DiffSPv1_2SE_eta[iside]->SetBinContent(ebin+1, SPv1_2SE_eta[iside][ebin]->GetMean());
            DiffSPv1_2SE_eta[iside]->SetBinError(ebin+1, SPv1_2SE_eta[iside][ebin]->GetMeanError());
            DiffSPv1_3SE_eta[iside]->SetBinContent(ebin+1, SPv1_3SE_eta[iside][ebin]->GetMean());
            DiffSPv1_3SE_eta[iside]->SetBinError(ebin+1, SPv1_3SE_eta[iside][ebin]->GetMeanError());
            DiffSPv1_mix_eta[iside]->SetBinContent(ebin+1, SPv1_mix_eta[iside][ebin]->GetMean());
            DiffSPv1_mix_eta[iside]->SetBinError(ebin+1, SPv1_mix_eta[iside][ebin]->GetMeanError());
            DiffSPv2_2SE_eta[iside]->SetBinContent(ebin+1, SPv2_2SE_eta[iside][ebin]->GetMean());
            DiffSPv2_2SE_eta[iside]->SetBinError(ebin+1, SPv2_2SE_eta[iside][ebin]->GetMeanError());
            DiffSPv2_3SE_eta[iside]->SetBinContent(ebin+1, SPv2_3SE_eta[iside][ebin]->GetMean());
            DiffSPv2_3SE_eta[iside]->SetBinError(ebin+1, SPv2_3SE_eta[iside][ebin]->GetMeanError());
        }
    }

    ofstream fout;
    if (!fopen("data","r")) system(Form("mkdir data"));
    TString tag = Form("data_final_tally_%s",mtag.Data());
    TString foutname = "data/"+tag+".dat";
    fout.open(foutname.Data());
    fout << "====================" << endl;
    fout << "nevents:    " << nevents*ntries << endl;
    fout << "event multiplicity: " << evtmult << endl;
    if (isodd) fout << "Rapidity-odd v1 input " << endl;
    else fout << "Rapidity-even v1:  " << setv1 << endl;
    fout << "setv2:     " << setv2 << endl;
    if (eta_weights) fout << "Using eta-dependent weights..." <<endl;
    if (pt_weights) fout << "Using pt-dependent weights..." <<endl;
    if (conserve_pT) fout << "Conserving transverse momentum event by event... " << endl;
    if (addholes) fout << "Holes in detector acceptance... " << endl;
    if (recenter) fout << "Recentered... " << endl;
    if (flatten) fout << "Flattened... " << endl;
    fout << "iseed:    " << iseed << endl;
    fout << "   -------   " << endl;


    fout<<"\n  ---- Final Tally (HFm) ---- \n"<<endl;
    fout<<"rescor1_2SE:     "<<Form("%.5f",rescor1_2SE->GetMean())<<" +/- "<<Form("%.5f",rescor1_2SE->GetMeanError())<<endl;
    fout<<"rescor1_3SE:     "<<Form("%.5f",rescor1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_3SE[0]->GetMeanError())<<endl;
    fout<<"rescor2_2SE:     "<<Form("%.5f",rescor2_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor2_2SE->GetMeanError())<<endl;
    fout<<"rescor2_2SE:     "<<Form("%.5f",rescor2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor2_3SE[0]->GetMeanError())<<endl;
    fout<<"rescor_Mix:      "<<Form("%.5f",rescor1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix[0]->GetMeanError())<<endl;
    fout<<"EPv1obs:         "<<Form("%.5f",EPv1obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obs[0]->GetMeanError())<<endl;
    fout<<"EPv2obs:         "<<Form("%.5f",EPv2obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2obs[0]->GetMeanError())<<endl;
    fout<<"EPv1obs_mix:     "<<Form("%.5f",EPv1obsMix[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"SPdenom1_2SE:    "<<Form("%.5f",SPdenom1_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_2SE->GetMeanError())<<endl;
    fout<<"SPdenom1_3SE:    "<<Form("%.5f",SPdenom1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_3SE[0]->GetMeanError())<<endl;
    fout<<"SPdenom2_2SE:    "<<Form("%.5f",SPdenom2_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom2_2SE->GetMeanError())<<endl;
    fout<<"SPdenom2_3SE:    "<<Form("%.5f",SPdenom2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2_3SE[0]->GetMeanError())<<endl;
    fout<<"SPdenom1_mix:    "<<Form("%.5f",SPdenom1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix[0]->GetMeanError())<<endl;
    fout<<"SPv1num:         "<<Form("%.5f",SPv1num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1num[0]->GetMeanError())<<endl;
    fout<<"SPv2num:         "<<Form("%.5f",SPv2num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2num[0]->GetMeanError())<<endl;
    fout<<"SPv1numMix:      "<<Form("%.5f",SPv1numMix[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"EPv1_2SE:        "<<Form("%.5f",EPv1_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_2SE[0]->GetMeanError())<<endl;
    fout<<"EPv1_3SE:        "<<Form("%.5f",EPv1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_3SE[0]->GetMeanError())<<endl;
    fout<<"EPv1_mix:        "<<Form("%.5f",EPv1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"SPv1_2SE:        "<<Form("%.5f",SPv1_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_2SE[0]->GetMeanError())<<endl;
    fout<<"SPv1_3SE:        "<<Form("%.5f",SPv1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_3SE[0]->GetMeanError())<<endl;
    fout<<"SPv1_mix:        "<<Form("%.5f",SPv1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix[0]->GetMeanError())<<endl;
    fout<<"   -----         "<<endl;
    fout<<"EPv2_2SE:        "<<Form("%.5f",EPv2_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2_2SE[0]->GetMeanError())<<endl;
    fout<<"EPv2_3SE:        "<<Form("%.5f",EPv2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2_3SE[0]->GetMeanError())<<endl;
    fout<<"SPv2_2SE:        "<<Form("%.5f",SPv2_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2_2SE[0]->GetMeanError())<<endl;
    fout<<"SPv2_3SE:        "<<Form("%.5f",SPv2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2_3SE[0]->GetMeanError())<<endl;


    fout<<"\n  ---- Final Tally (HFp) ---- \n"<<endl;
    fout<<"rescor1_2SE:    "<<Form("%.5f",rescor1_2SE->GetMean())<<" +/- "<<Form("%.5f",rescor1_2SE->GetMeanError())<<endl;
    fout<<"rescor1_3SE:    "<<Form("%.5f",rescor1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_3SE[1]->GetMeanError())<<endl;
    fout<<"rescor2_2SE:    "<<Form("%.5f",rescor2_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor2_2SE->GetMeanError())<<endl;
    fout<<"rescor2_2SE:    "<<Form("%.5f",rescor2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor2_3SE[1]->GetMeanError())<<endl;
    fout<<"rescor_mix:     "<<Form("%.5f",rescor1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix[1]->GetMeanError())<<endl;
    fout<<"EPv1obs:        "<<Form("%.5f",EPv1obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obs[1]->GetMeanError())<<endl;
    fout<<"EPv2obs:        "<<Form("%.5f",EPv2obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2obs[1]->GetMeanError())<<endl;
    fout<<"EPv1obsMix:     "<<Form("%.5f",EPv1obsMix[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix[1]->GetMeanError())<<endl;
    fout<<"   -----        "<<endl;
    fout<<"SPdenom1_2SE:   "<<Form("%.5f",SPdenom1_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_2SE->GetMeanError())<<endl;
    fout<<"SPdenom1_3SE:   "<<Form("%.5f",SPdenom1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_3SE[1]->GetMeanError())<<endl;
    fout<<"SPdenom2_2SE:   "<<Form("%.5f",SPdenom2_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom2_2SE->GetMeanError())<<endl;
    fout<<"SPdenom2_3SE:   "<<Form("%.5f",SPdenom2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2_3SE[1]->GetMeanError())<<endl;
    fout<<"SPdenom_mix:    "<<Form("%.5f",SPdenom1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix[1]->GetMeanError())<<endl;
    fout<<"SPv1num:        "<<Form("%.5f",SPv1num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1num[1]->GetMeanError())<<endl;
    fout<<"SPv2num:        "<<Form("%.5f",SPv2num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2num[1]->GetMeanError())<<endl;
    fout<<"SPv1numMix:     "<<Form("%.5f",SPv1numMix[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix[1]->GetMeanError())<<endl;
    fout<<"   -----        "<<endl;
    fout<<"EPv1_2SE:       "<<Form("%.5f",EPv1_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_2SE[1]->GetMeanError())<<endl;
    fout<<"EPv1_3SE:       "<<Form("%.5f",EPv1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_3SE[1]->GetMeanError())<<endl;
    fout<<"EPv1_mix:       "<<Form("%.5f",EPv1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix[1]->GetMeanError())<<endl;
    fout<<"   -----        "<<endl;
    fout<<"SPv1_2SE:       "<<Form("%.5f",SPv1_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_2SE[1]->GetMeanError())<<endl;
    fout<<"SPv1_3SE:       "<<Form("%.5f",SPv1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_3SE[1]->GetMeanError())<<endl;
    fout<<"SPv1_mix:       "<<Form("%.5f",SPv1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix[1]->GetMeanError())<<endl;
    fout<<"   -----        "<<endl;
    fout<<"EPv2_2SE:       "<<Form("%.5f",EPv2_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2_2SE[1]->GetMeanError())<<endl;
    fout<<"EPv2_3SE:       "<<Form("%.5f",EPv2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2_3SE[1]->GetMeanError())<<endl;
    fout<<"SPv2_2SE:       "<<Form("%.5f",SPv2_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2_2SE[1]->GetMeanError())<<endl;
    fout<<"SPv2_3SE:       "<<Form("%.5f",SPv2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2_3SE[1]->GetMeanError())<<endl;
    fout<<"\n --- Comparisons to input vn --- \n"<<endl;
    if (!isodd) {
        fout<<Form("EPv1_2SE/v1in  (HFm): %0.4f",EPv1_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_2SE[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("EPv1_3SE/v1in  (HFm): %0.4f",EPv1_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_3SE[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("EPv1_mix/v1in  (HFm): %0.4f",EPv1_mix[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_mix[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("SPv1_2SE/v1in  (HFm): %0.4f",SPv1_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_2SE[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("SPv1_3SE/v1in  (HFm): %0.4f",SPv1_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_3SE[1]->GetMeanError()/v1in)<<endl;
        fout<<Form("SPv1_mix/v1in  (HFm): %0.4f",SPv1_mix[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_mix[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix[1]->GetMeanError()/v1in)<<endl;
    }
    fout<<Form("EPv2_2SE/v1in  (HFm): %0.4f",EPv2_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv2_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv2_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv2_2SE[1]->GetMeanError()/v1in)<<endl;
    fout<<Form("EPv2_3SE/v1in  (HFm): %0.4f",EPv2_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv2_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv2_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv2_3SE[1]->GetMeanError()/v1in)<<endl;
    fout<<Form("SPv2_2SE/v1in  (HFm): %0.4f",SPv2_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv2_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv2_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv2_2SE[1]->GetMeanError()/v1in)<<endl;
    fout<<Form("SPv2_3SE/v1in  (HFm): %0.4f",SPv2_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv2_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv2_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv2_3SE[1]->GetMeanError()/v1in)<<endl;
    fout<<"\n ...done \n"<<endl;
    fout.close();


    cout<<"\n  ---- Final Tally (HFm) ---- \n"<<endl;
    cout<<"rescor1_2SE:     "<<Form("%.5f",rescor1_2SE->GetMean())<<" +/- "<<Form("%.5f",rescor1_2SE->GetMeanError())<<endl;
    cout<<"rescor1_3SE:     "<<Form("%.5f",rescor1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_3SE[0]->GetMeanError())<<endl;
    cout<<"rescor2_2SE:     "<<Form("%.5f",rescor2_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor2_2SE->GetMeanError())<<endl;
    cout<<"rescor2_2SE:     "<<Form("%.5f",rescor2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",rescor2_3SE[0]->GetMeanError())<<endl;
    cout<<"rescor_Mix:      "<<Form("%.5f",rescor1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix[0]->GetMeanError())<<endl;
    cout<<"EPv1obs:         "<<Form("%.5f",EPv1obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obs[0]->GetMeanError())<<endl;
    cout<<"EPv2obs:         "<<Form("%.5f",EPv2obs[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2obs[0]->GetMeanError())<<endl;
    cout<<"EPv1obs_mix:     "<<Form("%.5f",EPv1obsMix[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"SPdenom1_2SE:    "<<Form("%.5f",SPdenom1_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_2SE->GetMeanError())<<endl;
    cout<<"SPdenom1_3SE:    "<<Form("%.5f",SPdenom1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_3SE[0]->GetMeanError())<<endl;
    cout<<"SPdenom2_2SE:    "<<Form("%.5f",SPdenom2_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom2_2SE->GetMeanError())<<endl;
    cout<<"SPdenom2_3SE:    "<<Form("%.5f",SPdenom2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2_3SE[0]->GetMeanError())<<endl;
    cout<<"SPdenom1_mix:    "<<Form("%.5f",SPdenom1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix[0]->GetMeanError())<<endl;
    cout<<"SPv1num:         "<<Form("%.5f",SPv1num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1num[0]->GetMeanError())<<endl;
    cout<<"SPv2num:         "<<Form("%.5f",SPv2num[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2num[0]->GetMeanError())<<endl;
    cout<<"SPv1numMix:      "<<Form("%.5f",SPv1numMix[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"EPv1_2SE:        "<<Form("%.5f",EPv1_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_2SE[0]->GetMeanError())<<endl;
    cout<<"EPv1_3SE:        "<<Form("%.5f",EPv1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_3SE[0]->GetMeanError())<<endl;
    cout<<"EPv1_mix:        "<<Form("%.5f",EPv1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"SPv1_2SE:        "<<Form("%.5f",SPv1_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_2SE[0]->GetMeanError())<<endl;
    cout<<"SPv1_3SE:        "<<Form("%.5f",SPv1_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_3SE[0]->GetMeanError())<<endl;
    cout<<"SPv1_mix:        "<<Form("%.5f",SPv1_mix[0]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix[0]->GetMeanError())<<endl;
    cout<<"   -----         "<<endl;
    cout<<"EPv2_2SE:        "<<Form("%.5f",EPv2_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2_2SE[0]->GetMeanError())<<endl;
    cout<<"EPv2_3SE:        "<<Form("%.5f",EPv2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",EPv2_3SE[0]->GetMeanError())<<endl;
    cout<<"SPv2_2SE:        "<<Form("%.5f",SPv2_2SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2_2SE[0]->GetMeanError())<<endl;
    cout<<"SPv2_3SE:        "<<Form("%.5f",SPv2_3SE[0]->GetMean())<<" +/- "<<Form("%.5f",SPv2_3SE[0]->GetMeanError())<<endl;


    cout<<"\n  ---- Final Tally (HFp) ---- \n"<<endl;
    cout<<"rescor1_2SE:    "<<Form("%.5f",rescor1_2SE->GetMean())<<" +/- "<<Form("%.5f",rescor1_2SE->GetMeanError())<<endl;
    cout<<"rescor1_3SE:    "<<Form("%.5f",rescor1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_3SE[1]->GetMeanError())<<endl;
    cout<<"rescor2_2SE:    "<<Form("%.5f",rescor2_2SE->GetMean())<<" +/-" <<Form("%.5f",rescor2_2SE->GetMeanError())<<endl;
    cout<<"rescor2_2SE:    "<<Form("%.5f",rescor2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",rescor2_3SE[1]->GetMeanError())<<endl;
    cout<<"rescor_mix:     "<<Form("%.5f",rescor1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",rescor1_mix[1]->GetMeanError())<<endl;
    cout<<"EPv1obs:        "<<Form("%.5f",EPv1obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obs[1]->GetMeanError())<<endl;
    cout<<"EPv2obs:        "<<Form("%.5f",EPv2obs[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2obs[1]->GetMeanError())<<endl;
    cout<<"EPv1obsMix:     "<<Form("%.5f",EPv1obsMix[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1obsMix[1]->GetMeanError())<<endl;
    cout<<"   -----        "<<endl;
    cout<<"SPdenom1_2SE:   "<<Form("%.5f",SPdenom1_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_2SE->GetMeanError())<<endl;
    cout<<"SPdenom1_3SE:   "<<Form("%.5f",SPdenom1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_3SE[1]->GetMeanError())<<endl;
    cout<<"SPdenom2_2SE:   "<<Form("%.5f",SPdenom2_2SE->GetMean())<<" +/- "<<Form("%.5f",SPdenom2_2SE->GetMeanError())<<endl;
    cout<<"SPdenom2_3SE:   "<<Form("%.5f",SPdenom2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom2_3SE[1]->GetMeanError())<<endl;
    cout<<"SPdenom_mix:    "<<Form("%.5f",SPdenom1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",SPdenom1_mix[1]->GetMeanError())<<endl;
    cout<<"SPv1num:        "<<Form("%.5f",SPv1num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1num[1]->GetMeanError())<<endl;
    cout<<"SPv2num:        "<<Form("%.5f",SPv2num[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2num[1]->GetMeanError())<<endl;
    cout<<"SPv1numMix:     "<<Form("%.5f",SPv1numMix[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1numMix[1]->GetMeanError())<<endl;
    cout<<"   -----        "<<endl;
    cout<<"EPv1_2SE:       "<<Form("%.5f",EPv1_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_2SE[1]->GetMeanError())<<endl;
    cout<<"EPv1_3SE:       "<<Form("%.5f",EPv1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_3SE[1]->GetMeanError())<<endl;
    cout<<"EPv1_mix:       "<<Form("%.5f",EPv1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",EPv1_mix[1]->GetMeanError())<<endl;
    cout<<"   -----        "<<endl;
    cout<<"SPv1_2SE:       "<<Form("%.5f",SPv1_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_2SE[1]->GetMeanError())<<endl;
    cout<<"SPv1_3SE:       "<<Form("%.5f",SPv1_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_3SE[1]->GetMeanError())<<endl;
    cout<<"SPv1_mix:       "<<Form("%.5f",SPv1_mix[1]->GetMean())<<" +/- "<<Form("%.5f",SPv1_mix[1]->GetMeanError())<<endl;
    cout<<"   -----        "<<endl;
    cout<<"EPv2_2SE:       "<<Form("%.5f",EPv2_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2_2SE[1]->GetMeanError())<<endl;
    cout<<"EPv2_3SE:       "<<Form("%.5f",EPv2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",EPv2_3SE[1]->GetMeanError())<<endl;
    cout<<"SPv2_2SE:       "<<Form("%.5f",SPv2_2SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2_2SE[1]->GetMeanError())<<endl;
    cout<<"SPv2_3SE:       "<<Form("%.5f",SPv2_3SE[1]->GetMean())<<" +/- "<<Form("%.5f",SPv2_3SE[1]->GetMeanError())<<endl;
    cout<<"\n --- Comparisons to input vn --- \n"<<endl;
    if (!isodd) {
        cout<<Form("EPv1_2SE/v1in  (HFm): %0.4f",EPv1_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_2SE[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("EPv1_3SE/v1in  (HFm): %0.4f",EPv1_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_3SE[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("EPv1_mix/v1in  (HFm): %0.4f",EPv1_mix[0]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",EPv1_mix[1]->GetMean()/v1in)<<Form(" +/- %0.4f",EPv1_mix[1]->GetMeanError()/v1in)<<endl;
        cout<<"   -----        "<<endl;
        cout<<Form("SPv1_2SE/v1in  (HFm): %0.4f",SPv1_2SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_2SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_2SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_2SE[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("SPv1_3SE/v1in  (HFm): %0.4f",SPv1_3SE[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_3SE[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_3SE[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_3SE[1]->GetMeanError()/v1in)<<endl;
        cout<<Form("SPv1_mix/v1in  (HFm): %0.4f",SPv1_mix[0]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix[0]->GetMeanError()/v1in)<<Form("\t (HFp): %0.4f",SPv1_mix[1]->GetMean()/v1in)<<Form(" +/- %0.4f",SPv1_mix[1]->GetMeanError()/v1in)<<endl;
        cout<<"   -----        "<<endl;
    }
    cout<<Form("EPv2_2SE/v2in  (HFm): %0.4f",EPv2_2SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_2SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",EPv2_2SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_2SE[1]->GetMeanError()/v2in)<<endl;
    cout<<Form("EPv2_3SE/v2in  (HFm): %0.4f",EPv2_3SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_3SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",EPv2_3SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",EPv2_3SE[1]->GetMeanError()/v2in)<<endl;
    cout<<Form("SPv2_2SE/v2in  (HFm): %0.4f",SPv2_2SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_2SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",SPv2_2SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_2SE[1]->GetMeanError()/v2in)<<endl;
    cout<<Form("SPv2_3SE/v2in  (HFm): %0.4f",SPv2_3SE[0]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_3SE[0]->GetMeanError()/v2in)<<Form("\t (HFp): %0.4f",SPv2_3SE[1]->GetMean()/v2in)<<Form(" +/- %0.4f",SPv2_3SE[1]->GetMeanError()/v2in)<<endl;
    cout<<"\n ...done \n"<<endl;
}
