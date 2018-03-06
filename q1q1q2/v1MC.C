# include "src/RunCase.C"

using namespace std;

void v1MC()
{
    int ntries = 10;
    Double_t v1in = 0.000;
    Double_t v2in = 0.060;
    Bool_t isodd = kFALSE;
    Int_t NumEvnts = 1e4;
    Int_t nevents = NumEvnts/ntries;

    TString tag = "v1";
    if (isodd) tag+="_odd";
    else tag+=Form("_even_%0.4f",v1in);
    tag+=Form("_v2_%0.4f",v2in);
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
    cout << "Input v1: " << v1in << endl;
    cout << "Input v2: " << v2in << endl;
    cout << "Tracker event plane: " << ecutmin[2] << " < |eta| < " << ecutmax[2] << endl;
    cout << "Running " << NumEvnts << " events over " << ntries << " runs (" << nevents << " per run) \n" << endl;

    hinit_v1in = new TH1D("v1in", "", nv1Ebins, v1Ebins);
    hq112 = new TH1D("q112", "", 200, -0.1, 0.2);

    TStopwatch * sw0 = new TStopwatch();
    TStopwatch * sw1 = new TStopwatch();
    sw0->Start();
    for (int cindx = 0; cindx<ntries; cindx++) {
        sw1->Start();
        TDatime ts;
        iseed = ts.GetTime();
        RunCase( cindx, nevents, evtmult, v1in, v2in, iseed, tag.Data(), ntries );

        sw0->Continue();
        sw1->Continue();
        double elapse0 = sw0->RealTime();
        double elapse1 = sw1->RealTime();
        cout << "Elapse time: " << Form("%3.2f",elapse0) << " seconds" << endl;
        cout << "Estimated total run time: " << Form("%3.2f",elapse1*ntries/60./60.) << " hours\n" << endl;
    }

    cout<<"v1: "<<v1in<<"\tv2: "<<v2in<<"\tq112: "<<Form("%.5f",hq112->GetMean())<<" +/- "<<Form("%.5f",hq112->GetMeanError())<<endl;
}
