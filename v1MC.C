# include "src/Setup.h"
# include "src/RunCase.C"

using namespace std;

void v1MC( Double_t setv1 = 0.01, Double_t setv2 = 0.06, Double_t setv3 = 0.03, bool odd = false, Int_t NumE = 1e4, bool eweights = false, bool pweights = false, bool mom_cons = false )
{
    int ntries = 10;
    // Double_t v1in =  0.050;
    // Double_t v2in =  0.070;
    // Double_t v3in =  0.030;
    // Bool_t isodd = kFALSE;
    // Int_t NumEvnts = 1e3;
    // Bool_t eta_weights = kFALSE;
    // Bool_t pt_weights = kFALSE;
    // Bool_t conserve_pT = kFALSE;
    Double_t v1in =  setv1;
    Double_t v2in =  setv2;
    Double_t v3in =  setv3;
    Bool_t isodd = kFALSE; if (odd) isodd = kTRUE;
    Int_t NumEvnts = NumE;
    Bool_t eta_weights = kFALSE; if (eweights) eta_weights = kTRUE;
    Bool_t pt_weights = kFALSE; if (pweights) pt_weights = kTRUE;
    Bool_t conserve_pT = kFALSE; if (mom_cons) conserve_pT = kTRUE;
    Bool_t addholes = kFALSE;
    Bool_t flatten = kFALSE;
    Bool_t recenter = kFALSE;
    Int_t nevents = NumEvnts/ntries;

    TString tag = "v1";
    if (isodd) tag+="_odd";
    else tag+=Form("_even_%0.4f",v1in);
    tag+=Form("_v2_%0.4f",v2in);
    if (eta_weights) tag+="_eta_weights";
    if (pt_weights) tag+="_pt_weights";
    if (conserve_pT) tag+="_mom-cons";
    if (addholes) tag+="_hole";
    if (flatten) tag+="_flatten";
    if (recenter) tag+="_recenter";
    tag+=Form("_%0.1f_to_%0.1f",ecutmin[2],ecutmax[2]);
    tag+=Form("_%d_evts",NumEvnts);
    Int_t evtmult = 0;
    for (int vbin = 0; vbin<nv1Ebins; vbin++) {
        int multin = v1Emult[vbin];
        evtmult+=multin;
    }

    TH1::SetDefaultSumw2();
    gStyle->SetPalette(55);
    gStyle->SetErrorX(0.5);
    if (isodd) cout << "Rapidity-odd input v1" << endl;
    else cout << "Rapidity-even input v1: " << v1in << endl;
    cout << "Input v2: " << v2in << endl;
    cout << "Input v3: " << v3in << endl;
    if (eta_weights) cout << "Using eta-dependent weighting... " << endl;
    if (pt_weights) cout << "Using pt-dependent weighting... " << endl;
    if (conserve_pT) cout << "Conserving transverse momentum event by event... " << endl;
    if (addholes) cout << "Holes in detector acceptance... " << endl;
    cout << "Tracker event plane: " << ecutmin[2] << " < |eta| < " << ecutmax[2] << endl;
    cout << "Running " << NumEvnts << " events over " << ntries << " runs (" << nevents << " per run) \n" << endl;

    Setup();

    inputParms->SetBinContent(1, NumEvnts);
    if (isodd) inputParms->SetBinContent(2, 1);
    inputParms->SetBinContent(3, v2in);
    inputParms->SetBinContent(4, v3in);
    if (eta_weights)    inputParms->SetBinContent(5, 1);
    if (pt_weights)     inputParms->SetBinContent(6, 1);
    if (conserve_pT)    inputParms->SetBinContent(7, 1);
    if (recenter)       inputParms->SetBinContent(8, 1);
    if (flatten)        inputParms->SetBinContent(9, 1);
    inputParms->SetBinContent(10, ecutmin[2]);
    inputParms->SetBinContent(11, ecutmax[2]);

    TStopwatch * sw0 = new TStopwatch();
    TStopwatch * sw1 = new TStopwatch();
    sw0->Start();
    for (int cindx = 0; cindx<ntries; cindx++) {
        sw1->Start();
        TDatime ts;
        iseed = ts.GetTime();
        RunCase( cindx, nevents, evtmult, isodd, v1in, v2in, v3in, eta_weights, pt_weights, conserve_pT, addholes, flatten, recenter, iseed, tag.Data(), ntries );

        sw0->Continue();
        sw1->Continue();
        double elapse0 = sw0->RealTime();
        double elapse1 = sw1->RealTime();
        cout << "Elapse time: " << Form("%3.2f",elapse0) << " seconds" << endl;
        cout << "Estimated total run time: " << Form("%3.2f",elapse1*ntries/60./60.) << " hours\n" << endl;
    }
    ComputeVN( nevents, evtmult, isodd, v1in, v2in, v3in, eta_weights, pt_weights, conserve_pT, addholes, flatten, recenter, iseed, tag.Data(), ntries );

    WriteToFile( isodd, v1in, v2in, v3in, eta_weights, pt_weights, conserve_pT, addholes, flatten, recenter );
}
