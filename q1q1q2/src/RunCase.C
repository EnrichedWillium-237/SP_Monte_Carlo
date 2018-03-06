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

TH1D * hq112;

TFile * tfout;
Int_t iseed = 0;
Int_t counter = 0;


Double_t bounds(int ord, double ang) {
    while (ang >  TMath::Pi()/ord) ang-=TMath::TwoPi()/ord;
    while (ang < -TMath::Pi()/ord) ang+=TMath::TwoPi()/ord;
    return ang;
}

void RunCase( Int_t indx, Int_t nevents, Int_t evtmult, Double_t setv1, Double_t setv2, Int_t iseed, TString mtag, Int_t Ntries )
{

    int cindx = indx;
    TRandom3 * ran = new TRandom3(iseed+cindx);
    MCEvent * event[nv1Ebins];

    for (int vbin = 0; vbin<nv1Ebins; vbin++) {
        hinit_philab[vbin] = new TH1D(Form("init_philab_v%i_%i",vbin,cindx), "", 150, -4, 4);
        hinit_eta[vbin] = new TH1D(Form("init_eta_v%i_%i",vbin,cindx), "", 150, -6, 6);
        hinit_pt[vbin] = new TH1D(Form("init_pt_v%i_%i",vbin,cindx), "", 150, 0, 8);
        hinit_phiPsiRP[vbin] = new TH1D(Form("init_phiPsiRP_v%i_%i",vbin,cindx), "", 150, -4, 4);
        event[vbin] = new MCEvent(setv1, setv2, 0.0, 0.0, 0.0, 0.0);
    }

    Double_t c1[numEP];
    Double_t s1[numEP];
    Double_t c2[numEP];
    Double_t s2[numEP];

    Double_t q112 = 0;

    Double_t q1x[numEP];
    Double_t q1y[numEP];
    Double_t q2x[numEP];
    Double_t q2y[numEP];
    Double_t qcnt[numEP];
    Double_t Psi1[numEP];
    Double_t Psi2[numEP];
    Int_t subcnt[numEP];

    int nevts = 0;

    Double_t phiTracks_init[nv1Ebins][multMax];
    Double_t etaTracks_init[nv1Ebins][multMax];
    Double_t ptTracks_init[nv1Ebins][multMax];

    Double_t phiTracks_merge[multMax];
    Double_t etaTracks_merge[multMax];
    Double_t ptTracks_merge[multMax];

    Double_t phiTracks_sub[numEP][multMax];
    Double_t etaTracks_sub[numEP][multMax];
    Double_t ptTracks_sub[numEP][multMax];

    //-- main event loop

    cout << "Entering primary event loop" << endl;
    for (int vbin = 0; vbin<nv1Ebins; vbin++) {
        double v1tmp;
        v1tmp = setv1;
        int multtmp = v1Emult[vbin];
        event[vbin]->SetSeed(iseed+cindx);
        event[vbin]->Setv1(v1tmp);
        event[vbin]->SetMult(multtmp);
        event[vbin]->SetEventParms();
        hinit_v1in->SetBinContent(vbin+1,v1tmp);
        hinit_v1in->SetBinError(vbin+1,0.00001);
    }

    //-- make initial throws and merge into single event
    for (Int_t ievent = 0; ievent<nevents; ievent++) {
        if (fmod(double(ievent+1), nevents/20.)==0) cout << " event: " << ievent+1
            << "/" << nevents << "\trun: " << cindx+1 << "/" << Ntries <<  endl;

        event[0]->SetPsiRandom();
        Double_t PsiRP = event[0]->GetPsi();
//        Double_t PsiRP = ran->Uniform(-TMath::Pi(), TMath::Pi());
        Int_t jcnt = 0;
        for (int vbin = 0; vbin<nv1Ebins; vbin++) {
            double vmin = v1Ebins[vbin];
            double vmax = v1Ebins[vbin+1];
            event[vbin]->SetPsi(PsiRP);
            event[vbin]->GetThrowPhi(phiTracks_init[vbin]);
            event[vbin]->GetEtaRandom(etaTracks_init[vbin],vmin,vmax);
            event[vbin]->GetPtDist(ptTracks_init[vbin]);
            // event[vbin]->GetPtRandom(ptTracks_init[vbin]);
            for (Int_t j = 0; j<event[vbin]->GetMult(); j++) {
                double ph = bounds(1,phiTracks_init[vbin][j]);
                double eta = etaTracks_init[vbin][j];
                double pt = ptTracks_init[vbin][j];
                phiTracks_merge[jcnt] = ph;
                etaTracks_merge[jcnt] = eta;
                ptTracks_merge[jcnt] = pt;

                jcnt++;
            }
        }

        //-- sort particles into subevents
        for (int nep = 0; nep<numEP; nep++) {
            subcnt[nep] = 0;
        }
        for (Int_t j = 0; j<evtmult; j++) {
            double ph = bounds(1,phiTracks_merge[j]);
            double eta = etaTracks_merge[j];
            double pt = ptTracks_merge[j];
            // HF-, track-, track+, HF+, and trackmid
            bool isinsc = 0;
            for (int nep = 0; nep<numEP-1; nep++) {
                if (eta>=ecutmin[nep] && eta<ecutmax[nep] && pt>=pcutmin[nep] && pt<pcutmax[nep]) {
                    isinsc = 1;
                    phiTracks_sub[nep][subcnt[nep]] = ph;
                    etaTracks_sub[nep][subcnt[nep]] = eta;
                    ptTracks_sub[nep][subcnt[nep]] = pt;
                    ++subcnt[nep];
                }
            }
        }

        //-- calculate event plane angles
        for (int nep = 0; nep<numEP; nep++) {
            q1x[nep] = 0;
            q1y[nep] = 0;
            q2x[nep] = 0;
            q2y[nep] = 0;
            Psi1[nep] = 0;
            Psi2[nep] = 0;
            qcnt[nep] = 0;
            double ptxsub = 0;
            double ptysub = 0;

            for (Int_t j = 0; j<subcnt[nep]; j++) {
                double ph = bounds(1,phiTracks_sub[nep][j]);
                double eta = etaTracks_sub[nep][j];
                double pt = ptTracks_sub[nep][j];

                q1x[nep] += TMath::Cos(ph);
                q1y[nep] += TMath::Sin(ph);
                q2x[nep] += TMath::Cos(2*ph);
                q2y[nep] += TMath::Sin(2*ph);
                qcnt[nep]++;
            }
            Psi1[nep] = TMath::ATan2(q1y[nep],q1x[nep]);
            Psi2[nep] = TMath::ATan2(q2y[nep],q2x[nep])/2.;
        }

        //-- calculate vn
        q112 += TMath::Cos(2*(Psi1[1] - Psi2[3]));
        ++nevts;
    }
    q112/=nevts;
    hq112->Fill(q112);

    ofstream fout;
    if (!fopen("logs","r")) system(Form("mkdir logs"));
    TString tag = Form("logs_%s",mtag.Data());
    TString foutname = "logs/"+tag+".dat";
    fout.open(foutname.Data());
    fout << "====================" << endl;
    fout << "nevents:    " << nevents << endl;
    fout << "Input v1:  " << setv1 << endl;
    fout << "Input v2:  " << setv2 << endl;
    fout << "iseed:    " << iseed << endl;

    fout<<"   -----             "<<endl;
    fout<<"<cos2(Psi1 - Psi2)>: "<<Form("%.6f",q112)<<endl;
    fout<<"   -----             \n"<<endl;
    fout.close();

    cout<<"   -----             "<<endl;
    cout<<"<cos2(Psi1 - Psi2)>: "<<Form("%.6f",q112)<<endl;
    cout<<"   -----             \n"<<endl;

}
