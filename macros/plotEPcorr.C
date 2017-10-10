
Bool_t print_plot = kTRUE;
Bool_t close_plot = kTRUE;
Bool_t loop_plot_print = kTRUE;
Bool_t xprof_plot = kTRUE;

const double xtlsize1x1 = 0.07;
const double ytlsize1x1 = 0.07;
const double xtloffset1x1 = 1.0;
const double ytloffset1x1 = 1.15;
const double xlbsize1x1 = 0.05;
const double ylbsize1x1 = 0.05;
const double xlboffset1x1 = 0.01;
const double ylboffset1x1 = 0.01;

const double xtlsize_pan = 0.07;
const double ytlsize_pan = 0.07;
const double xtloffset_pan = 1.00;
const double ytloffset_pan = 0.4;
const double xlbsize_pan = 0.05;
const double ylbsize_pan = 0.05;
const double xlboffset_pan = 0.01;
const double ylboffset_pan = 0.01;

double rghtmrg = 0.13;
double rghtmrg_pan = 0.10;

// double zcolmin1v1 = 2400;
// double zcolmax1v1 = 3400;
double zcolmin1v1 = 800;
double zcolmax1v1 = 1600;

// double zcolmin1v2 = 2200;
// double zcolmax1v2 = 3600;
double zcolmin1v2 = 800;
double zcolmax1v2 = 1700;

// double zcolmin2v2 = -1;
// double zcolmax2v2 = 10000;
double zcolmin2v2 = -1;
double zcolmax2v2 = 4000;

double xprofymin = -0.20;
double xprofymax = 0.20;

double xprofymin2ndOrd = -1;
double xprofymax2ndOrd = 1;

Bool_t eta_weights;
Bool_t pt_weights;
Bool_t iseven;
Bool_t isodd;
Bool_t iszero;
Bool_t conserve_pT;
Bool_t ishole;
Bool_t recentered;

# include "EPplot.h"
# include "EPperiodic.h"
# include "EPXprofile.h"

void plotEPcorr()
{

    gStyle->SetPalette(55);
    int NumFile = 2;

    TString inFile[10];

    inFile[0] = "../results/results_v1_odd_v2_0.0700_eta_weights_5000000_evts_02.root";
    inFile[1] = "../results/results_v1_odd_v2_0.0700_eta_weights_2000000_evts_06.root";
    inFile[2] = "../results/results_v1_odd_v2_0.0700_eta_weights_2000000_evts_10.root";
    inFile[3] = "../results/results_v1_odd_v2_0.0700_eta_weights_2000000_evts_14.root";
    inFile[4] = "../results/results_v1_odd_v2_0.0700_eta_weights_2000000_evts_18.root";
    inFile[5] = "../results/results_v1_odd_v2_0.0700_eta_weights_5000000_evts_22.root";
    inFile[8] = "../results/results_v1_odd_v2_0.0700_eta_weights_5000000_evts_0to2p4.root";
    inFile[9] = "../results/results_v1_odd_v2_0.0700_eta_weights_1000000_evts_1p2to2p4.root";

    for (int filenum = NumFile; filenum<=NumFile; filenum++) {

        eta_weights = kFALSE;
        pt_weights = kFALSE;
        iseven = kFALSE;
        isodd = kFALSE;
        iszero = kFALSE;
        conserve_pT = kFALSE;
        ishole = kFALSE;
        recentered = kFALSE;

        cout << "\n Reading file "<< filenum <<": " << inFile[filenum] << "\n" << endl;
        EPplot( inFile[filenum] );

        if (loop_plot_print) EPperiodic();
        if (xprof_plot) EPXprofile();

        TString tag = "";
        if (iseven) tag+="v1even";
        if (isodd) tag+="v1odd";
        if (iszero) tag+="_zero";
        if (eta_weights) tag+="_eta-weights"; else tag+="_no_eta-weights";
        if (pt_weights) tag+="_pt-weights";
        if (conserve_pT) tag+="_mom-cons";
        if (ishole) tag+="_hole"; else tag+="_perfect_acceptence";
        if (recentered) tag+="_recentered";
        if (inFile[filenum].Contains("02")) tag+="_02";
        if (inFile[filenum].Contains("06")) tag+="_06";
        if (inFile[filenum].Contains("10")) tag+="_10";
        if (inFile[filenum].Contains("14")) tag+="_14";
        if (inFile[filenum].Contains("18")) tag+="_18";
        if (inFile[filenum].Contains("22")) tag+="_22";
        if (inFile[filenum].Contains("0to2p4")) tag+="_0to2p4";
        if (inFile[filenum].Contains("1p2to2p4")) tag+="_1p2to2p4";
        if (inFile[filenum].Contains("2p0to2p4")) tag+="_2p0to2p4";

        if (inFile[filenum].Contains("2000000")) {
            zcolmin1v1 = 0;
            zcolmax1v1 = 2000;

            zcolmin1v2 = 0;
            zcolmax1v2 = 2000;

            zcolmin2v2 = -1;
            zcolmax2v2 = 5000;
        }

        if (!fopen("plots","r")) system("mkdir plots");
        if (!fopen(Form("plots/%s",tag.Data()),"r")) system(Form("mkdir plots/%s",tag.Data()));
        if (!fopen(Form("plots/%s/EPCorrelations",tag.Data()),"r")) system(Form("mkdir plots/%s/EPCorrelations",tag.Data()));
        if (!fopen(Form("plots/%s/EPCorrelations/png",tag.Data()),"r")) system(Form("mkdir plots/%s/EPCorrelations/png",tag.Data()));
        if (!fopen(Form("plots/%s/EPCorrelations/pdf",tag.Data()),"r")) system(Form("mkdir plots/%s/EPCorrelations/pdf",tag.Data()));

        TCanvas * ccorr_HFm1_HFp1 = new TCanvas("ccorr_HFm1_HFp1","ccorr_HFm1_HFp1",600,530);
        TPad * padcorr_HFm1_HFp1 = (TPad *) ccorr_HFm1_HFp1->cd();
        padcorr_HFm1_HFp1->SetRightMargin(rghtmrg);
        corr_HFm1_HFp1->Draw();
        if (print_plot) ccorr_HFm1_HFp1->Print(Form("plots/%s/EPCorrelations/png/HFp1_vs_HFm1.png",tag.Data()),"png");
        if (print_plot) ccorr_HFm1_HFp1->Print(Form("plots/%s/EPCorrelations/pdf/HFp1_vs_HFm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_HFm1_HFp1->Close();


        TCanvas * ccorr_trackmid1_HFm1 = new TCanvas("ccorr_trackmid1_HFm1","ccorr_trackmid1_HFm1",600,530);
        TPad * padcorr_trackmid1_HFm1 = (TPad *) ccorr_trackmid1_HFm1->cd();
        padcorr_trackmid1_HFm1->SetRightMargin(rghtmrg);
        corr_trackmid1_HFm1->Draw();
        if (print_plot) ccorr_trackmid1_HFm1->Print(Form("plots/%s/EPCorrelations/png/HFm1_vs_trackmid1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackmid1_HFm1->Print(Form("plots/%s/EPCorrelations/pdf/HFm1_vs_trackmid1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackmid1_HFm1->Close();


        TCanvas * ccorr_trackmid1_HFp1 = new TCanvas("ccorr_trackmid1_HFp1","ccorr_trackmid1_HFp1",600,530);
        TPad * padcorr_trackmid1_HFp1 = (TPad *) ccorr_trackmid1_HFp1->cd();
        padcorr_trackmid1_HFp1->SetRightMargin(rghtmrg);
        corr_trackmid1_HFp1->Draw();
        if (print_plot) ccorr_trackmid1_HFp1->Print(Form("plots/%s/EPCorrelations/png/HFp1_vs_trackmid1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackmid1_HFp1->Print(Form("plots/%s/EPCorrelations/pdf/HFp1_vs_trackmid1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackmid1_HFp1->Close();


        TCanvas * ccorr_trackm1_HFm1 = new TCanvas("ccorr_trackm1_HFm1","ccorr_trackm1_HFm1",600,530);
        TPad * padcorr_trackm1_HFm1 = (TPad *) ccorr_trackm1_HFm1->cd();
        padcorr_trackm1_HFm1->SetRightMargin(rghtmrg);
        corr_trackm1_HFm1->Draw();
        if (print_plot) ccorr_trackm1_HFm1->Print(Form("plots/%s/EPCorrelations/png/HFm1_vs_trackm1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm1_HFm1->Print(Form("plots/%s/EPCorrelations/pdf/HFm1_vs_trackm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm1_HFm1->Close();


        TCanvas * ccorr_trackm1_trackp1 = new TCanvas("ccorr_trackm1_trackp1","ccorr_trackm1_trackp1",600,530);
        TPad * padcorr_trackm1_trackp1 = (TPad *) ccorr_trackm1_trackp1->cd();
        padcorr_trackm1_trackp1->SetRightMargin(rghtmrg);
        corr_trackm1_trackp1->Draw();
        if (print_plot) ccorr_trackm1_trackp1->Print(Form("plots/%s/EPCorrelations/png/trackp1_vs_trackm1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm1_trackp1->Print(Form("plots/%s/EPCorrelations/pdf/trackp1_vs_trackm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm1_trackp1->Close();


        TCanvas * ccorr_trackm1_HFp1 = new TCanvas("ccorr_trackm1_HFp1","ccorr_trackm1_HFp1",600,530);
        TPad * padcorr_trackm1_HFp1 = (TPad *) ccorr_trackm1_HFp1->cd();
        padcorr_trackm1_HFp1->SetRightMargin(rghtmrg);
        corr_trackm1_HFp1->Draw();
        if (print_plot) ccorr_trackm1_HFp1->Print(Form("plots/%s/EPCorrelations/png/HFp1_vs_trackm1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm1_HFp1->Print(Form("plots/%s/EPCorrelations/pdf/HFp1_vs_trackm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm1_HFp1->Close();


        TCanvas * ccorr_trackp1_HFm1 = new TCanvas("ccorr_trackp1_HFm1","ccorr_trackp1_HFm1",600,530);
        TPad * padcorr_trackp1_HFm1 = (TPad *) ccorr_trackp1_HFm1->cd();
        padcorr_trackp1_HFm1->SetRightMargin(rghtmrg);
        corr_trackp1_HFm1->Draw();
        if (print_plot) ccorr_trackp1_HFm1->Print(Form("plots/%s/EPCorrelations/png/HFm1_vs_trackp1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackp1_HFm1->Print(Form("plots/%s/EPCorrelations/pdf/HFm1_vs_trackp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackp1_HFm1->Close();


        TCanvas * ccorr_trackp1_trackm1 = new TCanvas("ccorr_trackp1_trackm1","ccorr_trackp1_trackm1",600,530);
        TPad * padcorr_trackp1_trackm1 = (TPad *) ccorr_trackp1_trackm1->cd();
        padcorr_trackp1_trackm1->SetRightMargin(rghtmrg);
        corr_trackp1_trackm1->Draw();
        if (print_plot) ccorr_trackp1_trackm1->Print(Form("plots/%s/EPCorrelations/png/trackm1_vs_trackp1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackp1_trackm1->Print(Form("plots/%s/EPCorrelations/pdf/trackm1_vs_trackp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackp1_trackm1->Close();


        TCanvas * ccorr_trackp1_HFp1 = new TCanvas("ccorr_trackp1_HFp1","ccorr_trackp1_HFp1",600,530);
        TPad * padcorr_trackp1_HFp1 = (TPad *) ccorr_trackp1_HFp1->cd();
        padcorr_trackp1_HFp1->SetRightMargin(rghtmrg);
        corr_trackp1_HFp1->Draw();
        if (print_plot) ccorr_trackp1_HFp1->Print(Form("plots/%s/EPCorrelations/png/HFp1_vs_trackp1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackp1_HFp1->Print(Form("plots/%s/EPCorrelations/pdf/HFp1_vs_trackp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackp1_HFp1->Close();


        TCanvas * ccorr_HFm1_trackm2 = new TCanvas("ccorr_HFm1_trackm2","ccorr_HFm1_trackm2",600,530);
        TPad * padcorr_HFm1_trackm2 = (TPad *) ccorr_HFm1_trackm2->cd();
        padcorr_HFm1_trackm2->SetRightMargin(rghtmrg);
        corr_HFm1_trackm2->Draw();
        if (print_plot) ccorr_HFm1_trackm2->Print(Form("plots/%s/EPCorrelations/png/trackm2_vs_HFm1.png",tag.Data()),"png");
        if (print_plot) ccorr_HFm1_trackm2->Print(Form("plots/%s/EPCorrelations/pdf/trackm2_vs_HFm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_HFm1_trackm2->Close();


        TCanvas * ccorr_HFm1_trackp2 = new TCanvas("ccorr_HFm1_trackp2","ccorr_HFm1_trackp2",600,530);
        TPad * padcorr_HFm1_trackp2 = (TPad *) ccorr_HFm1_trackp2->cd();
        padcorr_HFm1_trackp2->SetRightMargin(rghtmrg);
        corr_HFm1_trackp2->Draw();
        if (print_plot) ccorr_HFm1_trackp2->Print(Form("plots/%s/EPCorrelations/png/trackp2_vs_HFm1.png",tag.Data()),"png");
        if (print_plot) ccorr_HFm1_trackp2->Print(Form("plots/%s/EPCorrelations/pdf/trackp2_vs_HFm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_HFm1_trackp2->Close();


        TCanvas * ccorr_HFp1_trackm2 = new TCanvas("ccorr_HFp1_trackm2","ccorr_HFp1_trackm2",600,530);
        TPad * padcorr_HFp1_trackm2 = (TPad *) ccorr_HFp1_trackm2->cd();
        padcorr_HFp1_trackm2->SetRightMargin(rghtmrg);
        corr_HFp1_trackm2->Draw();
        if (print_plot) ccorr_HFp1_trackm2->Print(Form("plots/%s/EPCorrelations/png/trackm2_vs_HFp1.png",tag.Data()),"png");
        if (print_plot) ccorr_HFp1_trackm2->Print(Form("plots/%s/EPCorrelations/pdf/trackm2_vs_HFp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_HFp1_trackm2->Close();


        TCanvas * ccorr_HFp1_trackp2 = new TCanvas("ccorr_HFp1_trackp2","ccorr_HFp1_trackp2",600,530);
        TPad * padcorr_HFp1_trackp2 = (TPad *) ccorr_HFp1_trackp2->cd();
        padcorr_HFp1_trackp2->SetRightMargin(rghtmrg);
        corr_HFp1_trackp2->Draw();
        if (print_plot) ccorr_HFp1_trackp2->Print(Form("plots/%s/EPCorrelations/png/trackp2_vs_HFp1.png",tag.Data()),"png");
        if (print_plot) ccorr_HFp1_trackp2->Print(Form("plots/%s/EPCorrelations/pdf/trackp2_vs_HFp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_HFp1_trackp2->Close();


        TCanvas * ccorr_HFm1_HFp2 = new TCanvas("ccorr_HFm1_HFp2","ccorr_HFm1_HFp2",600,530);
        TPad * padcorr_HFm1_HFp2 = (TPad *) ccorr_HFm1_HFp2->cd();
        padcorr_HFm1_HFp2->SetRightMargin(rghtmrg);
        corr_HFm1_HFp2->Draw();
        if (print_plot) ccorr_HFm1_HFp2->Print(Form("plots/%s/EPCorrelations/png/HFp2_vs_HFm1.png",tag.Data()),"png");
        if (print_plot) ccorr_HFm1_HFp2->Print(Form("plots/%s/EPCorrelations/pdf/HFp2_vs_HFm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_HFm1_HFp2->Close();


        TCanvas * ccorr_HFp1_HFm2 = new TCanvas("ccorr_HFp1_HFm2","ccorr_HFp1_HFm2",600,530);
        TPad * padcorr_HFp1_HFm2 = (TPad *) ccorr_HFp1_HFm2->cd();
        padcorr_HFp1_HFm2->SetRightMargin(rghtmrg);
        corr_HFp1_HFm2->Draw();
        if (print_plot) ccorr_HFp1_HFm2->Print(Form("plots/%s/EPCorrelations/png/HFm2_vs_HFp1.png",tag.Data()),"png");
        if (print_plot) ccorr_HFp1_HFm2->Print(Form("plots/%s/EPCorrelations/pdf/HFm2_vs_HFp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_HFp1_HFm2->Close();


        TCanvas * ccorr_trackmid1_HFm2 = new TCanvas("ccorr_trackmid1_HFm2","ccorr_trackmid1_HFm2",600,530);
        TPad * padcorr_trackmid1_HFm2 = (TPad *) ccorr_trackmid1_HFm2->cd();
        padcorr_trackmid1_HFm2->SetRightMargin(rghtmrg);
        corr_trackmid1_HFm2->Draw();
        if (print_plot) ccorr_trackmid1_HFm2->Print(Form("plots/%s/EPCorrelations/png/HFm2_vs_trackmid1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackmid1_HFm2->Print(Form("plots/%s/EPCorrelations/pdf/HFm2_vs_trackmid1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackmid1_HFm2->Close();


        TCanvas * ccorr_trackmid1_HFp2 = new TCanvas("ccorr_trackmid1_HFp2","ccorr_trackmid1_HFp2",600,530);
        TPad * padcorr_trackmid1_HFp2 = (TPad *) ccorr_trackmid1_HFp2->cd();
        padcorr_trackmid1_HFp2->SetRightMargin(rghtmrg);
        corr_trackmid1_HFp2->Draw();
        if (print_plot) ccorr_trackmid1_HFp2->Print(Form("plots/%s/EPCorrelations/png/HFp2_vs_trackmid1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackmid1_HFp2->Print(Form("plots/%s/EPCorrelations/pdf/HFp2_vs_trackmid1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackmid1_HFp2->Close();


        TCanvas * ccorr_trackm1_HFm2 = new TCanvas("ccorr_trackm1_HFm2","ccorr_trackm1_HFm2",600,530);
        TPad * padcorr_trackm1_HFm2 = (TPad *) ccorr_trackm1_HFm2->cd();
        padcorr_trackm1_HFm2->SetRightMargin(rghtmrg);
        corr_trackm1_HFm2->Draw();
        if (print_plot) ccorr_trackm1_HFm2->Print(Form("plots/%s/EPCorrelations/png/HFm2_vs_trackm1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm1_HFm2->Print(Form("plots/%s/EPCorrelations/pdf/HFm2_vs_trackm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm1_HFm2->Close();


        TCanvas * ccorr_trackm1_trackp2 = new TCanvas("ccorr_trackm1_trackp2","ccorr_trackm1_trackp2",600,530);
        TPad * padcorr_trackm1_trackp2 = (TPad *) ccorr_trackm1_trackp2->cd();
        padcorr_trackm1_trackp2->SetRightMargin(rghtmrg);
        corr_trackm1_trackp2->Draw();
        if (print_plot) ccorr_trackm1_trackp2->Print(Form("plots/%s/EPCorrelations/png/trackp2_vs_trackm1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm1_trackp2->Print(Form("plots/%s/EPCorrelations/pdf/trackp2_vs_trackm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm1_trackp2->Close();


        TCanvas * ccorr_trackm1_HFp2 = new TCanvas("ccorr_trackm1_HFp2","ccorr_trackm1_HFp2",600,530);
        TPad * padcorr_trackm1_HFp2 = (TPad *) ccorr_trackm1_HFp2->cd();
        padcorr_trackm1_HFp2->SetRightMargin(rghtmrg);
        corr_trackm1_HFp2->Draw();
        if (print_plot) ccorr_trackm1_HFp2->Print(Form("plots/%s/EPCorrelations/png/HFp2_vs_trackm1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm1_HFp2->Print(Form("plots/%s/EPCorrelations/pdf/HFp2_vs_trackm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm1_HFp2->Close();


        TCanvas * ccorr_trackp1_HFm2 = new TCanvas("ccorr_trackp1_HFm2","ccorr_trackp1_HFm2",600,530);
        TPad * padcorr_trackp1_HFm2 = (TPad *) ccorr_trackp1_HFm2->cd();
        padcorr_trackp1_HFm2->SetRightMargin(rghtmrg);
        corr_trackp1_HFm2->Draw();
        if (print_plot) ccorr_trackp1_HFm2->Print(Form("plots/%s/EPCorrelations/png/HFm2_vs_trackp1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackp1_HFm2->Print(Form("plots/%s/EPCorrelations/pdf/HFm2_vs_trackp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackp1_HFm2->Close();


        TCanvas * ccorr_trackp1_trackm2 = new TCanvas("ccorr_trackp1_trackm2","ccorr_trackp1_trackm2",600,530);
        TPad * padcorr_trackp1_trackm2 = (TPad *) ccorr_trackp1_trackm2->cd();
        padcorr_trackp1_trackm2->SetRightMargin(rghtmrg);
        corr_trackp1_trackm2->Draw();
        if (print_plot) ccorr_trackp1_trackm2->Print(Form("plots/%s/EPCorrelations/png/trackm2_vs_trackp1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackp1_trackm2->Print(Form("plots/%s/EPCorrelations/pdf/trackm2_vs_trackp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackp1_trackm2->Close();


        TCanvas * ccorr_trackp1_HFp2 = new TCanvas("ccorr_trackp1_HFp2","ccorr_trackp1_HFp2",600,530);
        TPad * padcorr_trackp1_HFp2 = (TPad *) ccorr_trackp1_HFp2->cd();
        padcorr_trackp1_HFp2->SetRightMargin(rghtmrg);
        corr_trackp1_HFp2->Draw();
        if (print_plot) ccorr_trackp1_HFp2->Print(Form("plots/%s/EPCorrelations/png/HFp2_vs_trackp1.png",tag.Data()),"png");
        if (print_plot) ccorr_trackp1_HFp2->Print(Form("plots/%s/EPCorrelations/pdf/HFp2_vs_trackp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackp1_HFp2->Close();


        TCanvas * ccorr_HFm2_HFp2 = new TCanvas("ccorr_HFm2_HFp2","ccorr_HFm2_HFp2",600,530);
        TPad * padcorr_HFm2_HFp2 = (TPad *) ccorr_HFm2_HFp2->cd();
        padcorr_HFm2_HFp2->SetRightMargin(rghtmrg);
        corr_HFm2_HFp2->Draw();
        if (print_plot) ccorr_HFm2_HFp2->Print(Form("plots/%s/EPCorrelations/png/HFp2_vs_HFm2.png",tag.Data()),"png");
        if (print_plot) ccorr_HFm2_HFp2->Print(Form("plots/%s/EPCorrelations/pdf/HFp2_vs_HFm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_HFm2_HFp2->Close();


        TCanvas * ccorr_trackmid2_HFm2 = new TCanvas("ccorr_trackmid2_HFm2","ccorr_trackmid2_HFm2",600,530);
        TPad * padcorr_trackmid2_HFm2 = (TPad *) ccorr_trackmid2_HFm2->cd();
        padcorr_trackmid2_HFm2->SetRightMargin(rghtmrg);
        corr_trackmid2_HFm2->Draw();
        if (print_plot) ccorr_trackmid2_HFm2->Print(Form("plots/%s/EPCorrelations/png/HFm2_vs_trackmid2.png",tag.Data()),"png");
        if (print_plot) ccorr_trackmid2_HFm2->Print(Form("plots/%s/EPCorrelations/pdf/HFm2_vs_trackmid2.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackmid2_HFm2->Close();


        TCanvas * ccorr_trackmid2_HFp2 = new TCanvas("ccorr_trackmid2_HFp2","ccorr_trackmid2_HFp2",600,530);
        TPad * padcorr_trackmid2_HFp2 = (TPad *) ccorr_trackmid2_HFp2->cd();
        padcorr_trackmid2_HFp2->SetRightMargin(rghtmrg);
        corr_trackmid2_HFp2->Draw();
        if (print_plot) ccorr_trackmid2_HFp2->Print(Form("plots/%s/EPCorrelations/png/HFp2_vs_trackmid2.png",tag.Data()),"png");
        if (print_plot) ccorr_trackmid2_HFp2->Print(Form("plots/%s/EPCorrelations/pdf/HFp2_vs_trackmid2.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackmid2_HFp2->Close();


        TCanvas * ccorr_trackm2_HFm2 = new TCanvas("ccorr_trackm2_HFm2","ccorr_trackm2_HFm2",600,530);
        TPad * padcorr_trackm2_HFm2 = (TPad *) ccorr_trackm2_HFm2->cd();
        padcorr_trackm2_HFm2->SetRightMargin(rghtmrg);
        corr_trackm2_HFm2->Draw();
        if (print_plot) ccorr_trackm2_HFm2->Print(Form("plots/%s/EPCorrelations/png/HFm2_vs_trackm2.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm2_HFm2->Print(Form("plots/%s/EPCorrelations/pdf/HFm2_vs_trackm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm2_HFm2->Close();


        TCanvas * ccorr_trackm2_trackp2 = new TCanvas("ccorr_trackm2_trackp2","ccorr_trackm2_trackp2",600,530);
        TPad * padcorr_trackm2_trackp2 = (TPad *) ccorr_trackm2_trackp2->cd();
        padcorr_trackm2_trackp2->SetRightMargin(rghtmrg);
        corr_trackm2_trackp2->Draw();
        if (print_plot) ccorr_trackm2_trackp2->Print(Form("plots/%s/EPCorrelations/png/trackp2_vs_trackm2.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm2_trackp2->Print(Form("plots/%s/EPCorrelations/pdf/trackp2_vs_trackm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm2_trackp2->Close();


        TCanvas * ccorr_trackm2_HFp2 = new TCanvas("ccorr_trackm2_HFp2","ccorr_trackm2_HFp2",600,530);
        TPad * padcorr_trackm2_HFp2 = (TPad *) ccorr_trackm2_HFp2->cd();
        padcorr_trackm2_HFp2->SetRightMargin(rghtmrg);
        corr_trackm2_HFp2->Draw();
        if (print_plot) ccorr_trackm2_HFp2->Print(Form("plots/%s/EPCorrelations/png/HFp2_vs_trackm2.png",tag.Data()),"png");
        if (print_plot) ccorr_trackm2_HFp2->Print(Form("plots/%s/EPCorrelations/pdf/HFp2_vs_trackm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackm2_HFp2->Close();


        TCanvas * ccorr_trackp2_HFm2 = new TCanvas("ccorr_trackp2_HFm2","ccorr_trackp2_HFm2",600,530);
        TPad * padcorr_trackp2_HFm2 = (TPad *) ccorr_trackp2_HFm2->cd();
        padcorr_trackp2_HFm2->SetRightMargin(rghtmrg);
        corr_trackp2_HFm2->Draw();
        if (print_plot) ccorr_trackp2_HFm2->Print(Form("plots/%s/EPCorrelations/png/HFm2_vs_trackp2.png",tag.Data()),"png");
        if (print_plot) ccorr_trackp2_HFm2->Print(Form("plots/%s/EPCorrelations/pdf/HFm2_vs_trackp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackp2_HFm2->Close();


        TCanvas * ccorr_trackp2_HFp2 = new TCanvas("ccorr_trackp2_HFp2","ccorr_trackp2_HFp2",600,530);
        TPad * padcorr_trackp2_HFp2 = (TPad *) ccorr_trackp2_HFp2->cd();
        padcorr_trackp2_HFp2->SetRightMargin(rghtmrg);
        corr_trackp2_HFp2->Draw();
        if (print_plot) ccorr_trackp2_HFp2->Print(Form("plots/%s/EPCorrelations/png/HFp2_vs_trackp2.png",tag.Data()),"png");
        if (print_plot) ccorr_trackp2_HFp2->Print(Form("plots/%s/EPCorrelations/pdf/HFp2_vs_trackp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccorr_trackp2_HFp2->Close();


        //----------------

        if (!fopen(Form("plots/%s/EPCosines/",tag.Data()),"r")) system(Form("mkdir plots/%s/EPCosines",tag.Data()));
        if (!fopen(Form("plots/%s/EPCosines/pdf",tag.Data()),"r")) system(Form("mkdir plots/%s/EPCosines/pdf",tag.Data()));

        TPaveText * txCosMean = new TPaveText(0.34, 0.79, 0.70, 0.88, "NDC");
        txCosMean->SetFillColor(0);
        txCosMean->SetBorderSize(0);
        txCosMean->SetTextFont(43);
        txCosMean->SetTextSize(30);


        TCanvas * ccos_HFm1_HFp1 = new TCanvas("ccos_HFm1_HFp1","ccos_HFm1_HFp1",600,530);
        ccos_HFm1_HFp1->cd();
        cos_HFm1_HFp1->Draw();
        TPaveText * txcos_HFm1_HFp1 = (TPaveText *) txCosMean->Clone();
        txcos_HFm1_HFp1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_HFm1_HFp1->GetMean(),cos_HFm1_HFp1->GetMeanError()));
        txcos_HFm1_HFp1->Draw();
        if (print_plot) ccos_HFm1_HFp1->Print(Form("plots/%s/EPCosines/pdf/cos_HFm1_HFp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_HFm1_HFp1->Close();

        TCanvas * ccos_trackmid1_HFm1 = new TCanvas("ccos_trackmid1_HFm1","ccos_trackmid1_HFm1",600,530);
        ccos_trackmid1_HFm1->cd();
        cos_trackmid1_HFm1->Draw();
        TPaveText * txcos_trackmid1_HFm1 = (TPaveText *) txCosMean->Clone();
        txcos_trackmid1_HFm1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackmid1_HFm1->GetMean(),cos_trackmid1_HFm1->GetMeanError()));
        txcos_trackmid1_HFm1->Draw();
        if (print_plot) ccos_trackmid1_HFm1->Print(Form("plots/%s/EPCosines/pdf/cos_trackmid1_HFm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackmid1_HFm1->Close();

        TCanvas * ccos_trackmid1_HFp1 = new TCanvas("ccos_trackmid1_HFp1","ccos_trackmid1_HFp1",600,530);
        ccos_trackmid1_HFp1->cd();
        cos_trackmid1_HFp1->Draw();
        TPaveText * txcos_trackmid1_HFp1 = (TPaveText *) txCosMean->Clone();
        txcos_trackmid1_HFp1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackmid1_HFp1->GetMean(),cos_trackmid1_HFp1->GetMeanError()));
        txcos_trackmid1_HFp1->Draw();
        if (print_plot) ccos_trackmid1_HFp1->Print(Form("plots/%s/EPCosines/pdf/cos_trackmid1_HFp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackmid1_HFp1->Close();

        TCanvas * ccos_trackm1_HFm1 = new TCanvas("ccos_trackm1_HFm1","ccos_trackm1_HFm1",600,530);
        ccos_trackm1_HFm1->cd();
        cos_trackm1_HFm1->Draw();
        TPaveText * txcos_trackm1_HFm1 = (TPaveText *) txCosMean->Clone();
        txcos_trackm1_HFm1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm1_HFm1->GetMean(),cos_trackm1_HFm1->GetMeanError()));
        txcos_trackm1_HFm1->Draw();
        if (print_plot) ccos_trackm1_HFm1->Print(Form("plots/%s/EPCosines/pdf/cos_trackm1_HFm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm1_HFm1->Close();

        TCanvas * ccos_trackm1_trackp1 = new TCanvas("cos_trackm1_trackp1","cos_trackm1_trackp1",600,530);
        ccos_trackm1_trackp1->cd();
        cos_trackm1_trackp1->Draw();
        TPaveText * txcos_trackm1_trackp1 = (TPaveText *) txCosMean->Clone();
        txcos_trackm1_trackp1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm1_trackp1->GetMean(),cos_trackm1_trackp1->GetMeanError()));
        txcos_trackm1_trackp1->Draw();
        if (print_plot) ccos_trackm1_trackp1->Print(Form("plots/%s/EPCosines/pdf/cos_trackm1_trackp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm1_trackp1->Close();

        TCanvas * ccos_trackm1_HFp1 = new TCanvas("ccos_trackm1_HFp1","ccos_trackm1_HFp1",600,530);
        ccos_trackm1_HFp1->cd();
        cos_trackm1_HFp1->Draw();
        TPaveText * txcos_trackm1_HFp1 = (TPaveText *) txCosMean->Clone();
        txcos_trackm1_HFp1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm1_HFp1->GetMean(),cos_trackm1_HFp1->GetMeanError()));
        txcos_trackm1_HFp1->Draw();
        if (print_plot) ccos_trackm1_HFp1->Print(Form("plots/%s/EPCosines/pdf/cos_trackm1_HFp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm1_HFp1->Close();

        TCanvas * ccos_trackp1_HFm1 = new TCanvas("ccos_trackp1_HFm1","ccos_trackp1_HFm1",600,530);
        ccos_trackp1_HFm1->cd();
        cos_trackp1_HFm1->Draw();
        TPaveText * txcos_trackp1_HFm1 = (TPaveText *) txCosMean->Clone();
        txcos_trackp1_HFm1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackp1_HFm1->GetMean(),cos_trackp1_HFm1->GetMeanError()));
        txcos_trackp1_HFm1->Draw();
        if (print_plot) ccos_trackp1_HFm1->Print(Form("plots/%s/EPCosines/pdf/cos_trackp1_HFm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackp1_HFm1->Close();

        TCanvas * ccos_trackp1_trackm1 = new TCanvas("ccos_trackp1_trackm1","ccos_trackp1_trackm1",600,530);
        ccos_trackp1_trackm1->cd();
        cos_trackp1_trackm1->Draw();
        TPaveText * txcos_trackp1_trackm1 = (TPaveText *) txCosMean->Clone();
        txcos_trackp1_trackm1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackp1_trackm1->GetMean(),cos_trackp1_trackm1->GetMeanError()));
        txcos_trackp1_trackm1->Draw();
        if (print_plot) ccos_trackp1_trackm1->Print(Form("plots/%s/EPCosines/pdf/cos_trackp1_trackm1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackp1_trackm1->Close();

        TCanvas * ccos_trackp1_HFp1 = new TCanvas("ccos_trackp1_HFp1","ccos_trackp1_HFp1",600,530);
        ccos_trackp1_HFp1->cd();
        cos_trackp1_HFp1->Draw();
        TPaveText * txcos_trackp1_HFp1 = (TPaveText *) txCosMean->Clone();
        txcos_trackp1_HFp1->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackp1_HFp1->GetMean(),cos_trackp1_HFp1->GetMeanError()));
        txcos_trackp1_HFp1->Draw();
        if (print_plot) ccos_trackp1_HFp1->Print(Form("plots/%s/EPCosines/pdf/cos_trackp1_HFp1.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackp1_HFp1->Close();


        TCanvas * ccos_HFm1_trackm2 = new TCanvas("ccos_HFm1_trackm2","ccos_HFm1_trackm2",600,530);
        ccos_HFm1_trackm2->cd();
        cos_HFm1_trackm2->Draw();
        TPaveText * txcos_HFm1_trackm2 = (TPaveText *) txCosMean->Clone();
        txcos_HFm1_trackm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_HFm1_trackm2->GetMean(),cos_HFm1_trackm2->GetMeanError()));
        txcos_HFm1_trackm2->Draw();
        if (print_plot) ccos_HFm1_trackm2->Print(Form("plots/%s/EPCosines/pdf/cos_HFm1_trackm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_HFm1_trackm2->Close();

        TCanvas * ccos_HFm1_trackp2 = new TCanvas("ccos_HFm1_trackp2","ccos_HFm1_trackp2",600,530);
        ccos_HFm1_trackp2->cd();
        cos_HFm1_trackp2->Draw();
        TPaveText * txcos_HFm1_trackp2 = (TPaveText *) txCosMean->Clone();
        txcos_HFm1_trackp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_HFm1_trackp2->GetMean(),cos_HFm1_trackp2->GetMeanError()));
        txcos_HFm1_trackp2->Draw();
        if (print_plot) ccos_HFm1_trackp2->Print(Form("plots/%s/EPCosines/pdf/cos_HFm1_trackp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_HFm1_trackp2->Close();

        TCanvas * ccos_HFp1_trackm2 = new TCanvas("ccos_HFp1_trackm2","ccos_HFp1_trackm2",600,530);
        ccos_HFp1_trackm2->cd();
        cos_HFp1_trackm2->Draw();
        TPaveText * txcos_HFp1_trackm2 = (TPaveText *) txCosMean->Clone();
        txcos_HFp1_trackm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_HFp1_trackm2->GetMean(),cos_HFp1_trackm2->GetMeanError()));
        txcos_HFp1_trackm2->Draw();
        if (print_plot) ccos_HFp1_trackm2->Print(Form("plots/%s/EPCosines/pdf/cos_HFp1_trackm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_HFp1_trackm2->Close();

        TCanvas * ccos_HFp1_trackp2 = new TCanvas("ccos_HFp1_trackp2","ccos_HFp1_trackp2",600,530);
        ccos_HFp1_trackp2->cd();
        cos_HFp1_trackp2->Draw();
        TPaveText * txcos_HFp1_trackp2 = (TPaveText *) txCosMean->Clone();
        txcos_HFp1_trackp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_HFp1_trackp2->GetMean(),cos_HFp1_trackp2->GetMeanError()));
        txcos_HFp1_trackp2->Draw();
        if (print_plot) ccos_HFp1_trackp2->Print(Form("plots/%s/EPCosines/pdf/cos_HFp1_trackp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_HFp1_trackp2->Close();

        TCanvas * ccos_HFm1_HFp2 = new TCanvas("ccos_HFm1_HFp2","ccos_HFm1_HFp2",600,530);
        ccos_HFm1_HFp2->cd();
        cos_HFm1_HFp2->Draw();
        TPaveText * txcos_HFm1_HFp2 = (TPaveText *) txCosMean->Clone();
        txcos_HFm1_HFp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_HFm1_HFp2->GetMean(),cos_HFm1_HFp2->GetMeanError()));
        txcos_HFm1_HFp2->Draw();
        if (print_plot) ccos_HFm1_HFp2->Print(Form("plots/%s/EPCosines/pdf/cos_HFm1_HFp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_HFm1_HFp2->Close();

        TCanvas * ccos_HFp1_HFm2 = new TCanvas("ccos_HFp1_HFm2","ccos_HFp1_HFm2",600,530);
        ccos_HFp1_HFm2->cd();
        cos_HFp1_HFm2->Draw();
        TPaveText * txcos_HFp1_HFm2 = (TPaveText *) txCosMean->Clone();
        txcos_HFp1_HFm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_HFp1_HFm2->GetMean(),cos_HFp1_HFm2->GetMeanError()));
        txcos_HFp1_HFm2->Draw();
        if (print_plot) ccos_HFp1_HFm2->Print(Form("plots/%s/EPCosines/pdf/cos_HFp1_HFm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_HFp1_HFm2->Close();

        TCanvas * ccos_trackmid1_HFm2 = new TCanvas("ccos_trackmid1_HFm2","ccos_trackmid1_HFm2",600,530);
        ccos_trackmid1_HFm2->cd();
        cos_trackmid1_HFm2->Draw();
        TPaveText * txcos_trackmid1_HFm2 = (TPaveText *) txCosMean->Clone();
        txcos_trackmid1_HFm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackmid1_HFm2->GetMean(),cos_trackmid1_HFm2->GetMeanError()));
        txcos_trackmid1_HFm2->Draw();
        if (print_plot) ccos_trackmid1_HFm2->Print(Form("plots/%s/EPCosines/pdf/cos_trackmid1_HFm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackmid1_HFm2->Close();

        TCanvas * ccos_trackmid1_HFp2 = new TCanvas("ccos_trackmid1_HFp2","ccos_trackmid1_HFp2",600,530);
        ccos_trackmid1_HFp2->cd();
        cos_trackmid1_HFp2->Draw();
        TPaveText * txcos_trackmid1_HFp2 = (TPaveText *) txCosMean->Clone();
        txcos_trackmid1_HFp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackmid1_HFp2->GetMean(),cos_trackmid1_HFp2->GetMeanError()));
        txcos_trackmid1_HFp2->Draw();
        if (print_plot) ccos_trackmid1_HFp2->Print(Form("plots/%s/EPCosines/pdf/cos_trackmid1_HFp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackmid1_HFp2->Close();

        TCanvas * ccos_trackm1_HFm2 = new TCanvas("ccos_trackm1_HFm2","ccos_trackm1_HFm2",600,530);
        ccos_trackm1_HFm2->cd();
        cos_trackm1_HFm2->Draw();
        TPaveText * txcos_trackm1_HFm2 = (TPaveText *) txCosMean->Clone();
        txcos_trackm1_HFm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm1_HFm2->GetMean(),cos_trackm1_HFm2->GetMeanError()));
        txcos_trackm1_HFm2->Draw();
        if (print_plot) ccos_trackm1_HFm2->Print(Form("plots/%s/EPCosines/pdf/cos_trackm1_HFm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm1_HFm2->Close();

        TCanvas * ccos_trackm1_trackp2 = new TCanvas("ccos_trackm1_trackp2","ccos_trackm1_trackp2",600,530);
        ccos_trackm1_trackp2->cd();
        cos_trackm1_trackp2->Draw();
        TPaveText * txcos_trackm1_trackp2 = (TPaveText *) txCosMean->Clone();
        txcos_trackm1_trackp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm1_trackp2->GetMean(),cos_trackm1_trackp2->GetMeanError()));
        txcos_trackm1_trackp2->Draw();
        if (print_plot) ccos_trackm1_trackp2->Print(Form("plots/%s/EPCosines/pdf/cos_trackm1_trackp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm1_trackp2->Close();

        TCanvas * ccos_trackm1_HFp2 = new TCanvas("ccos_trackm1_HFp2","ccos_trackm1_HFp2",600,530);
        ccos_trackm1_HFp2->cd();
        cos_trackm1_HFp2->Draw();
        TPaveText * txcos_trackm1_HFp2 = (TPaveText *) txCosMean->Clone();
        txcos_trackm1_HFp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm1_HFp2->GetMean(),cos_trackm1_HFp2->GetMeanError()));
        txcos_trackm1_HFp2->Draw();
        if (print_plot) ccos_trackm1_HFp2->Print(Form("plots/%s/EPCosines/pdf/cos_trackm1_HFp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm1_HFp2->Close();

        TCanvas * ccos_trackp1_HFm2 = new TCanvas("ccos_trackp1_HFm2","ccos_trackp1_HFm2",600,530);
        ccos_trackp1_HFm2->cd();
        cos_trackp1_HFm2->Draw();
        TPaveText * txcos_trackp1_HFm2 = (TPaveText *) txCosMean->Clone();
        txcos_trackp1_HFm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackp1_HFm2->GetMean(),cos_trackp1_HFm2->GetMeanError()));
        txcos_trackp1_HFm2->Draw();
        if (print_plot) ccos_trackp1_HFm2->Print(Form("plots/%s/EPCosines/pdf/cos_trackp1_HFm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackp1_HFm2->Close();

        TCanvas * ccos_trackp1_trackm2 = new TCanvas("ccos_trackp1_trackm2","ccos_trackp1_trackm2",600,530);
        ccos_trackp1_trackm2->cd();
        cos_trackp1_trackm2->Draw();
        TPaveText * txcos_trackp1_trackm2 = (TPaveText *) txCosMean->Clone();
        txcos_trackp1_trackm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackp1_trackm2->GetMean(),cos_trackp1_trackm2->GetMeanError()));
        txcos_trackp1_trackm2->Draw();
        if (print_plot) ccos_trackp1_trackm2->Print(Form("plots/%s/EPCosines/pdf/cos_trackp1_trackm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackp1_trackm2->Close();

        TCanvas * ccos_trackp1_HFp2 = new TCanvas("ccos_trackp1_HFp2","ccos_trackp1_HFp2",600,530);
        ccos_trackp1_HFp2->cd();
        cos_trackp1_HFp2->Draw();
        TPaveText * txcos_trackp1_HFp2 = (TPaveText *) txCosMean->Clone();
        txcos_trackp1_HFp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackp1_HFp2->GetMean(),cos_trackp1_HFp2->GetMeanError()));
        txcos_trackp1_HFp2->Draw();
        if (print_plot) ccos_trackp1_HFp2->Print(Form("plots/%s/EPCosines/pdf/cos_trackp1_HFp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackp1_HFp2->Close();


        TCanvas * ccos_HFm2_HFp2 = new TCanvas("ccos_HFm2_HFp2","ccos_HFm2_HFp2",600,530);
        ccos_HFm2_HFp2->cd();
        cos_HFm2_HFp2->Draw();
        TPaveText * txcos_HFm2_HFp2 = (TPaveText *) txCosMean->Clone();
        txcos_HFm2_HFp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_HFm2_HFp2->GetMean(),cos_HFm2_HFp2->GetMeanError()));
        txcos_HFm2_HFp2->Draw();
        if (print_plot) ccos_HFm2_HFp2->Print(Form("plots/%s/EPCosines/pdf/cos_HFm2_HFp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_HFm2_HFp2->Close();

        TCanvas * ccos_trackmid2_HFm2 = new TCanvas("ccos_trackmid2_HFm2","ccos_trackmid2_HFm2",600,530);
        ccos_trackmid2_HFm2->cd();
        cos_trackmid2_HFm2->Draw();
        TPaveText * txcos_trackmid2_HFm2 = (TPaveText *) txCosMean->Clone();
        txcos_trackmid2_HFm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackmid2_HFm2->GetMean(),cos_trackmid2_HFm2->GetMeanError()));
        txcos_trackmid2_HFm2->Draw();
        if (print_plot) ccos_trackmid2_HFm2->Print(Form("plots/%s/EPCosines/pdf/cos_trackmid2_HFm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackmid2_HFm2->Close();

        TCanvas * ccos_trackmid2_HFp2 = new TCanvas("ccos_trackmid2_HFp2","ccos_trackmid2_HFp2",600,530);
        ccos_trackmid2_HFp2->cd();
        cos_trackmid2_HFp2->Draw();
        TPaveText * txcos_trackmid2_HFp2 = (TPaveText *) txCosMean->Clone();
        txcos_trackmid2_HFp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackmid2_HFp2->GetMean(),cos_trackmid2_HFp2->GetMeanError()));
        txcos_trackmid2_HFp2->Draw();
        if (print_plot) ccos_trackmid2_HFp2->Print(Form("plots/%s/EPCosines/pdf/cos_trackmid2_HFp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackmid2_HFp2->Close();

        TCanvas * ccos_trackm2_HFm2 = new TCanvas("ccos_trackm2_HFm2","ccos_trackm2_HFm2",600,530);
        ccos_trackm2_HFm2->cd();
        cos_trackm2_HFm2->Draw();
        TPaveText * txcos_trackm2_HFm2 = (TPaveText *) txCosMean->Clone();
        txcos_trackm2_HFm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm2_HFm2->GetMean(),cos_trackm2_HFm2->GetMeanError()));
        txcos_trackm2_HFm2->Draw();
        if (print_plot) ccos_trackm2_HFm2->Print(Form("plots/%s/EPCosines/pdf/cos_trackm2_HFm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm2_HFm2->Close();

        TCanvas * ccos_trackm2_trackp2 = new TCanvas("ccos_trackm2_trackp2","ccos_trackm2_trackp2",600,530);
        ccos_trackm2_trackp2->cd();
        cos_trackm2_trackp2->Draw();
        TPaveText * txcos_trackm2_trackp2 = (TPaveText *) txCosMean->Clone();
        txcos_trackm2_trackp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm2_trackp2->GetMean(),cos_trackm2_trackp2->GetMeanError()));
        txcos_trackm2_trackp2->Draw();
        if (print_plot) ccos_trackm2_trackp2->Print(Form("plots/%s/EPCosines/pdf/cos_trackm2_trackp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm2_trackp2->Close();

        TCanvas * ccos_trackm2_HFp2 = new TCanvas("ccos_trackm2_HFp2","ccos_trackm2_HFp2",600,530);
        ccos_trackm2_HFp2->cd();
        cos_trackm2_HFp2->Draw();
        TPaveText * txcos_trackm2_HFp2 = (TPaveText *) txCosMean->Clone();
        txcos_trackm2_HFp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackm2_HFp2->GetMean(),cos_trackm2_HFp2->GetMeanError()));
        txcos_trackm2_HFp2->Draw();
        if (print_plot) ccos_trackm2_HFp2->Print(Form("plots/%s/EPCosines/pdf/cos_trackm2_HFp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackm2_HFp2->Close();

        TCanvas * ccos_trackp2_HFm2 = new TCanvas("ccos_trackp2_HFm2","ccos_trackp2_HFm2",600,530);
        ccos_trackp2_HFm2->cd();
        cos_trackp2_HFm2->Draw();
        TPaveText * txcos_trackp2_HFm2 = (TPaveText *) txCosMean->Clone();
        txcos_trackp2_HFm2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackp2_HFm2->GetMean(),cos_trackp2_HFm2->GetMeanError()));
        txcos_trackp2_HFm2->Draw();
        if (print_plot) ccos_trackp2_HFm2->Print(Form("plots/%s/EPCosines/pdf/cos_trackp2_HFm2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackp2_HFm2->Close();

        TCanvas * ccos_trackp2_HFp2 = new TCanvas("ccos_trackp2_HFp2","ccos_trackp2_HFp2",600,530);
        ccos_trackp2_HFp2->cd();
        cos_trackp2_HFp2->Draw();
        TPaveText * txcos_trackp2_HFp2 = (TPaveText *) txCosMean->Clone();
        txcos_trackp2_HFp2->AddText(Form("Mean: %0.4f #pm %0.4f",cos_trackp2_HFp2->GetMean(),cos_trackp2_HFp2->GetMeanError()));
        txcos_trackp2_HFp2->Draw();
        if (print_plot) ccos_trackp2_HFp2->Print(Form("plots/%s/EPCosines/pdf/cos_trackp2_HFp2.pdf",tag.Data()),"pdf");
        if (close_plot) ccos_trackp2_HFp2->Close();



        //-- periodicity plots
        if (loop_plot_print) {
            if (!fopen(Form("plots/%s/EPperiodic",tag.Data()),"r")) system(Form("mkdir plots/%s/EPperiodic",tag.Data()));
            if (!fopen(Form("plots/%s/EPperiodic/png",tag.Data()),"r")) system(Form("mkdir plots/%s/EPperiodic/png",tag.Data()));
            if (!fopen(Form("plots/%s/EPperiodic/pdf",tag.Data()),"r")) system(Form("mkdir plots/%s/EPperiodic/pdf",tag.Data()));


            TCanvas * cpanx_corr_HFm1_HFp1 = new TCanvas("cpanx_corr_HFm1_HFp1","cpanx_corr_HFm1_HFp1",1100,400);
            TPad * padpanx_corr_HFm1_HFp1 = (TPad *) cpanx_corr_HFm1_HFp1->cd();
            padpanx_corr_HFm1_HFp1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_HFm1_HFp1->SetLeftMargin(0.1);
            panx_corr_HFm1_HFp1->Draw();
            if (print_plot) cpanx_corr_HFm1_HFp1->Print(Form("plots/%s/EPperiodic/png/panx_corr_HFm1_HFp1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_HFm1_HFp1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_HFm1_HFp1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_HFm1_HFp1->Close();

            TCanvas * cpanx_corr_trackmid1_HFm1 = new TCanvas("cpanx_corr_trackmid1_HFm1","cpanx_corr_trackmid1_HFm1",1100,400);
            TPad * padpanx_corr_trackmid1_HFm1 = (TPad *) cpanx_corr_trackmid1_HFm1->cd();
            padpanx_corr_trackmid1_HFm1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackmid1_HFm1->SetLeftMargin(0.1);
            panx_corr_trackmid1_HFm1->Draw();
            if (print_plot) cpanx_corr_trackmid1_HFm1->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackmid1_HFm1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackmid1_HFm1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackmid1_HFm1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackmid1_HFm1->Close();

            TCanvas * cpanx_corr_trackmid1_HFp1 = new TCanvas("cpanx_corr_trackmid1_HFp1","cpanx_corr_trackmid1_HFp1",1100,400);
            TPad * padpanx_corr_trackmid1_HFp1 = (TPad *) cpanx_corr_trackmid1_HFp1->cd();
            padpanx_corr_trackmid1_HFp1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackmid1_HFp1->SetLeftMargin(0.1);
            panx_corr_trackmid1_HFp1->Draw();
            if (print_plot) cpanx_corr_trackmid1_HFp1->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackmid1_HFp1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackmid1_HFp1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackmid1_HFp1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackmid1_HFp1->Close();

            TCanvas * cpanx_corr_trackm1_HFm1 = new TCanvas("cpanx_corr_trackm1_HFm1","cpanx_corr_trackm1_HFm1",1100,400);
            TPad * padpanx_corr_trackm1_HFm1 = (TPad *) cpanx_corr_trackm1_HFm1->cd();
            padpanx_corr_trackm1_HFm1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm1_HFm1->SetLeftMargin(0.1);
            panx_corr_trackm1_HFm1->Draw();
            if (print_plot) cpanx_corr_trackm1_HFm1->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm1_HFm1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm1_HFm1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm1_HFm1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm1_HFm1->Close();

            TCanvas * cpanx_corr_trackm1_trackp1 = new TCanvas("cpanx_corr_trackm1_trackp1","cpanx_corr_trackm1_trackp1",1100,400);
            TPad * padpanx_corr_trackm1_trackp1 = (TPad *) cpanx_corr_trackm1_trackp1->cd();
            padpanx_corr_trackm1_trackp1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm1_trackp1->SetLeftMargin(0.1);
            panx_corr_trackm1_trackp1->Draw();
            if (print_plot) cpanx_corr_trackm1_trackp1->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm1_trackp1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm1_trackp1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm1_trackp1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm1_trackp1->Close();

            TCanvas * cpanx_corr_trackm1_HFp1 = new TCanvas("cpanx_corr_trackm1_HFp1","cpanx_corr_trackm1_HFp1",1100,400);
            TPad * padpanx_corr_trackm1_HFp1 = (TPad *) cpanx_corr_trackm1_HFp1->cd();
            padpanx_corr_trackm1_HFp1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm1_HFp1->SetLeftMargin(0.1);
            panx_corr_trackm1_HFp1->Draw();
            if (print_plot) cpanx_corr_trackm1_HFp1->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm1_HFp1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm1_HFp1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm1_HFp1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm1_HFp1->Close();

            TCanvas * cpanx_corr_trackp1_HFm1 = new TCanvas("cpanx_corr_trackp1_HFm1","cpanx_corr_trackp1_HFm1",1100,400);
            TPad * padpanx_corr_trackp1_HFm1 = (TPad *) cpanx_corr_trackp1_HFm1->cd();
            padpanx_corr_trackp1_HFm1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackp1_HFm1->SetLeftMargin(0.1);
            panx_corr_trackp1_HFm1->Draw();
            if (print_plot) cpanx_corr_trackp1_HFm1->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackp1_HFm1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackp1_HFm1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackp1_HFm1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackp1_HFm1->Close();

            TCanvas * cpanx_corr_trackp1_trackm1 = new TCanvas("cpanx_corr_trackp1_trackm1","cpanx_corr_trackp1_trackm1",1100,400);
            TPad * padpanx_corr_trackp1_trackm1 = (TPad *) cpanx_corr_trackp1_trackm1->cd();
            padpanx_corr_trackp1_trackm1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackp1_trackm1->SetLeftMargin(0.1);
            panx_corr_trackp1_trackm1->Draw();
            if (print_plot) cpanx_corr_trackp1_trackm1->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackp1_trackm1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackp1_trackm1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackp1_trackm1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackp1_trackm1->Close();

            TCanvas * cpanx_corr_trackp1_HFp1 = new TCanvas("cpanx_corr_trackp1_HFp1","cpanx_corr_trackp1_HFp1",1100,400);
            TPad * padpanx_corr_trackp1_HFp1 = (TPad *) cpanx_corr_trackp1_HFp1->cd();
            padpanx_corr_trackp1_HFp1->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackp1_HFp1->SetLeftMargin(0.1);
            panx_corr_trackp1_HFp1->Draw();
            if (print_plot) cpanx_corr_trackp1_HFp1->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackp1_HFp1.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackp1_HFp1->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackp1_HFp1.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackp1_HFp1->Close();


            TCanvas * cpanx_corr_HFm1_trackm2 = new TCanvas("cpanx_corr_HFm1_trackm2","cpanx_corr_HFm1_trackm2",1100,400);
            TPad * padpanx_corr_HFm1_trackm2 = (TPad *) cpanx_corr_HFm1_trackm2->cd();
            padpanx_corr_HFm1_trackm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_HFm1_trackm2->SetLeftMargin(0.1);
            panx_corr_HFm1_trackm2->Draw();
            if (print_plot) cpanx_corr_HFm1_trackm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_HFm1_trackm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_HFm1_trackm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_HFm1_trackm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_HFm1_trackm2->Close();

            TCanvas * cpanx_corr_HFm1_trackp2 = new TCanvas("cpanx_corr_HFm1_trackp2","cpanx_corr_HFm1_trackp2",1100,400);
            TPad * padpanx_corr_HFm1_trackp2 = (TPad *) cpanx_corr_HFm1_trackp2->cd();
            padpanx_corr_HFm1_trackp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_HFm1_trackp2->SetLeftMargin(0.1);
            panx_corr_HFm1_trackp2->Draw();
            if (print_plot) cpanx_corr_HFm1_trackp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_HFm1_trackp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_HFm1_trackp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_HFm1_trackp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_HFm1_trackp2->Close();

            TCanvas * cpanx_corr_HFp1_trackm2 = new TCanvas("cpanx_corr_HFp1_trackm2","cpanx_corr_HFp1_trackm2",1100,400);
            TPad * padpanx_corr_HFp1_trackm2 = (TPad *) cpanx_corr_HFp1_trackm2->cd();
            padpanx_corr_HFp1_trackm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_HFp1_trackm2->SetLeftMargin(0.1);
            panx_corr_HFp1_trackm2->Draw();
            if (print_plot) cpanx_corr_HFp1_trackm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_HFp1_trackm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_HFp1_trackm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_HFp1_trackm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_HFp1_trackm2->Close();

            TCanvas * cpanx_corr_HFp1_trackp2 = new TCanvas("cpanx_corr_HFp1_trackp2","cpanx_corr_HFp1_trackp2",1100,400);
            TPad * padpanx_corr_HFp1_trackp2 = (TPad *) cpanx_corr_HFp1_trackp2->cd();
            padpanx_corr_HFp1_trackp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_HFp1_trackp2->SetLeftMargin(0.1);
            panx_corr_HFp1_trackp2->Draw();
            if (print_plot) cpanx_corr_HFp1_trackp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_HFp1_trackp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_HFp1_trackp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_HFp1_trackp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_HFp1_trackp2->Close();

            TCanvas * cpanx_corr_HFm1_HFp2 = new TCanvas("cpanx_corr_HFm1_HFp2","cpanx_corr_HFm1_HFp2",1100,400);
            TPad * padpanx_corr_HFm1_HFp2 = (TPad *) cpanx_corr_HFm1_HFp2->cd();
            padpanx_corr_HFm1_HFp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_HFm1_HFp2->SetLeftMargin(0.1);
            panx_corr_HFm1_HFp2->Draw();
            if (print_plot) cpanx_corr_HFm1_HFp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_HFm1_HFp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_HFm1_HFp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_HFm1_HFp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_HFm1_HFp2->Close();

            TCanvas * cpanx_corr_HFp1_HFm2 = new TCanvas("cpanx_corr_HFp1_HFm2","cpanx_corr_HFp1_HFm2",1100,400);
            TPad * padpanx_corr_HFp1_HFm2 = (TPad *) cpanx_corr_HFp1_HFm2->cd();
            padpanx_corr_HFp1_HFm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_HFp1_HFm2->SetLeftMargin(0.1);
            panx_corr_HFp1_HFm2->Draw();
            if (print_plot) cpanx_corr_HFp1_HFm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_HFp1_HFm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_HFp1_HFm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_HFp1_HFm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_HFp1_HFm2->Close();

            TCanvas * cpanx_corr_trackmid1_HFm2 = new TCanvas("cpanx_corr_trackmid1_HFm2","cpanx_corr_trackmid1_HFm2",1100,400);
            TPad * padpanx_corr_trackmid1_HFm2 = (TPad *) cpanx_corr_trackmid1_HFm2->cd();
            padpanx_corr_trackmid1_HFm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackmid1_HFm2->SetLeftMargin(0.1);
            panx_corr_trackmid1_HFm2->Draw();
            if (print_plot) cpanx_corr_trackmid1_HFm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackmid1_HFm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackmid1_HFm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackmid1_HFm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackmid1_HFm2->Close();

            TCanvas * cpanx_corr_trackmid1_HFp2 = new TCanvas("cpanx_corr_trackmid1_HFp2","cpanx_corr_trackmid1_HFp2",1100,400);
            TPad * padpanx_corr_trackmid1_HFp2 = (TPad *) cpanx_corr_trackmid1_HFp2->cd();
            padpanx_corr_trackmid1_HFp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackmid1_HFp2->SetLeftMargin(0.1);
            panx_corr_trackmid1_HFp2->Draw();
            if (print_plot) cpanx_corr_trackmid1_HFp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackmid1_HFp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackmid1_HFp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackmid1_HFp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackmid1_HFp2->Close();

            TCanvas * cpanx_corr_trackm1_HFm2 = new TCanvas("cpanx_corr_trackm1_HFm2","cpanx_corr_trackm1_HFm2",1100,400);
            TPad * padpanx_corr_trackm1_HFm2 = (TPad *) cpanx_corr_trackm1_HFm2->cd();
            padpanx_corr_trackm1_HFm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm1_HFm2->SetLeftMargin(0.1);
            panx_corr_trackm1_HFm2->Draw();
            if (print_plot) cpanx_corr_trackm1_HFm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm1_HFm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm1_HFm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm1_HFm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm1_HFm2->Close();

            TCanvas * cpanx_corr_trackm1_trackp2 = new TCanvas("cpanx_corr_trackm1_trackp2","cpanx_corr_trackm1_trackp2",1100,400);
            TPad * padpanx_corr_trackm1_trackp2 = (TPad *) cpanx_corr_trackm1_trackp2->cd();
            padpanx_corr_trackm1_trackp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm1_trackp2->SetLeftMargin(0.1);
            panx_corr_trackm1_trackp2->Draw();
            if (print_plot) cpanx_corr_trackm1_trackp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm1_trackp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm1_trackp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm1_trackp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm1_trackp2->Close();

            TCanvas * cpanx_corr_trackm1_HFp2 = new TCanvas("cpanx_corr_trackm1_HFp2","cpanx_corr_trackm1_HFp2",1100,400);
            TPad * padpanx_corr_trackm1_HFp2 = (TPad *) cpanx_corr_trackm1_HFp2->cd();
            padpanx_corr_trackm1_HFp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm1_HFp2->SetLeftMargin(0.1);
            panx_corr_trackm1_HFp2->Draw();
            if (print_plot) cpanx_corr_trackm1_HFp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm1_HFp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm1_HFp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm1_HFp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm1_HFp2->Close();

            TCanvas * cpanx_corr_trackp1_HFm2 = new TCanvas("cpanx_corr_trackp1_HFm2","cpanx_corr_trackp1_HFm2",1100,400);
            TPad * padpanx_corr_trackp1_HFm2 = (TPad *) cpanx_corr_trackp1_HFm2->cd();
            padpanx_corr_trackp1_HFm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackp1_HFm2->SetLeftMargin(0.1);
            panx_corr_trackp1_HFm2->Draw();
            if (print_plot) cpanx_corr_trackp1_HFm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackp1_HFm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackp1_HFm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackp1_HFm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackp1_HFm2->Close();

            TCanvas * cpanx_corr_trackp1_trackm2 = new TCanvas("cpanx_corr_trackp1_trackm2","cpanx_corr_trackp1_trackm2",1100,400);
            TPad * padpanx_corr_trackp1_trackm2 = (TPad *) cpanx_corr_trackp1_trackm2->cd();
            padpanx_corr_trackp1_trackm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackp1_trackm2->SetLeftMargin(0.1);
            panx_corr_trackp1_trackm2->Draw();
            if (print_plot) cpanx_corr_trackp1_trackm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackp1_trackm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackp1_trackm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackp1_trackm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackp1_trackm2->Close();

            TCanvas * cpanx_corr_trackp1_HFp2 = new TCanvas("cpanx_corr_trackp1_HFp2","cpanx_corr_trackp1_HFp2",1100,400);
            TPad * padpanx_corr_trackp1_HFp2 = (TPad *) cpanx_corr_trackp1_HFp2->cd();
            padpanx_corr_trackp1_HFp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackp1_HFp2->SetLeftMargin(0.1);
            panx_corr_trackp1_HFp2->Draw();
            if (print_plot) cpanx_corr_trackp1_HFp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackp1_HFp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackp1_HFp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackp1_HFp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackp1_HFp2->Close();


            TCanvas * cpanx_corr_HFm2_HFp2 = new TCanvas("cpanx_corr_HFm2_HFp2","cpanx_corr_HFm2_HFp2",1100,400);
            TPad * padpanx_corr_HFm2_HFp2 = (TPad *) cpanx_corr_HFm2_HFp2->cd();
            padpanx_corr_HFm2_HFp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_HFm2_HFp2->SetLeftMargin(0.1);
            panx_corr_HFm2_HFp2->Draw();
            if (print_plot) cpanx_corr_HFm2_HFp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_HFm2_HFp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_HFm2_HFp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_HFm2_HFp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_HFm2_HFp2->Close();

            TCanvas * cpanx_corr_trackmid2_HFm2 = new TCanvas("cpanx_corr_trackmid2_HFm2","cpanx_corr_trackmid2_HFm2",1100,400);
            TPad * padpanx_corr_trackmid2_HFm2 = (TPad *) cpanx_corr_trackmid2_HFm2->cd();
            padpanx_corr_trackmid2_HFm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackmid2_HFm2->SetLeftMargin(0.1);
            panx_corr_trackmid2_HFm2->Draw();
            if (print_plot) cpanx_corr_trackmid2_HFm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackmid2_HFm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackmid2_HFm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackmid2_HFm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackmid2_HFm2->Close();

            TCanvas * cpanx_corr_trackmid2_HFp2 = new TCanvas("cpanx_corr_trackmid2_HFp2","cpanx_corr_trackmid2_HFp2",1100,400);
            TPad * padpanx_corr_trackmid2_HFp2 = (TPad *) cpanx_corr_trackmid2_HFp2->cd();
            padpanx_corr_trackmid2_HFp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackmid2_HFp2->SetLeftMargin(0.1);
            panx_corr_trackmid2_HFp2->Draw();
            if (print_plot) cpanx_corr_trackmid2_HFp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackmid2_HFp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackmid2_HFp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackmid2_HFp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackmid2_HFp2->Close();

            TCanvas * cpanx_corr_trackm2_HFm2 = new TCanvas("cpanx_corr_trackm2_HFm2","cpanx_corr_trackm2_HFm2",1100,400);
            TPad * padpanx_corr_trackm2_HFm2 = (TPad *) cpanx_corr_trackm2_HFm2->cd();
            padpanx_corr_trackm2_HFm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm2_HFm2->SetLeftMargin(0.1);
            panx_corr_trackm2_HFm2->Draw();
            if (print_plot) cpanx_corr_trackm2_HFm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm2_HFm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm2_HFm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm2_HFm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm2_HFm2->Close();

            TCanvas * cpanx_corr_trackm2_trackp2 = new TCanvas("cpanx_corr_trackm2_trackp2","cpanx_corr_trackm2_trackp2",1100,400);
            TPad * padpanx_corr_trackm2_trackp2 = (TPad *) cpanx_corr_trackm2_trackp2->cd();
            padpanx_corr_trackm2_trackp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm2_trackp2->SetLeftMargin(0.1);
            panx_corr_trackm2_trackp2->Draw();
            if (print_plot) cpanx_corr_trackm2_trackp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm2_trackp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm2_trackp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm2_trackp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm2_trackp2->Close();

            TCanvas * cpanx_corr_trackm2_HFp2 = new TCanvas("cpanx_corr_trackm2_HFp2","cpanx_corr_trackm2_HFp2",1100,400);
            TPad * padpanx_corr_trackm2_HFp2 = (TPad *) cpanx_corr_trackm2_HFp2->cd();
            padpanx_corr_trackm2_HFp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackm2_HFp2->SetLeftMargin(0.1);
            panx_corr_trackm2_HFp2->Draw();
            if (print_plot) cpanx_corr_trackm2_HFp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackm2_HFp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackm2_HFp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackm2_HFp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackm2_HFp2->Close();

            TCanvas * cpanx_corr_trackp2_HFm2 = new TCanvas("cpanx_corr_trackp2_HFm2","cpanx_corr_trackp2_HFm2",1100,400);
            TPad * padpanx_corr_trackp2_HFm2 = (TPad *) cpanx_corr_trackp2_HFm2->cd();
            padpanx_corr_trackp2_HFm2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackp2_HFm2->SetLeftMargin(0.1);
            panx_corr_trackp2_HFm2->Draw();
            if (print_plot) cpanx_corr_trackp2_HFm2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackp2_HFm2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackp2_HFm2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackp2_HFm2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackp2_HFm2->Close();

            TCanvas * cpanx_corr_trackp2_HFp2 = new TCanvas("cpanx_corr_trackp2_HFp2","cpanx_corr_trackp2_HFp2",1100,400);
            TPad * padpanx_corr_trackp2_HFp2 = (TPad *) cpanx_corr_trackp2_HFp2->cd();
            padpanx_corr_trackp2_HFp2->SetRightMargin(rghtmrg_pan);
            padpanx_corr_trackp2_HFp2->SetLeftMargin(0.1);
            panx_corr_trackp2_HFp2->Draw();
            if (print_plot) cpanx_corr_trackp2_HFp2->Print(Form("plots/%s/EPperiodic/png/panx_corr_trackp2_HFp2.png",tag.Data()),"png");
            if (print_plot) cpanx_corr_trackp2_HFp2->Print(Form("plots/%s/EPperiodic/pdf/panx_corr_trackp2_HFp2.pdf",tag.Data()),"pdf");
            if (close_plot) cpanx_corr_trackp2_HFp2->Close();

        }

        //-- draw x-profiles of 2D histos
        if (xprof_plot) EPXprofPrint( tag );
    }

}
