# include "TMath.h"

TH2D * panx_corr_HFm1_HFp1;
TH2D * panx_corr_trackmid1_HFm1;
TH2D * panx_corr_trackmid1_HFp1;
TH2D * panx_corr_trackm1_HFm1;
TH2D * panx_corr_trackm1_trackp1;
TH2D * panx_corr_trackm1_HFp1;
TH2D * panx_corr_trackp1_HFm1;
TH2D * panx_corr_trackp1_trackm1;
TH2D * panx_corr_trackp1_HFp1;

TH2D * panx_corr_HFm1_trackm2;
TH2D * panx_corr_HFm1_trackp2;
TH2D * panx_corr_HFp1_trackm2;
TH2D * panx_corr_HFp1_trackp2;
TH2D * panx_corr_HFm1_HFp2;
TH2D * panx_corr_HFp1_HFm2;
TH2D * panx_corr_trackmid1_HFm2;
TH2D * panx_corr_trackmid1_HFp2;
TH2D * panx_corr_trackm1_HFm2;
TH2D * panx_corr_trackm1_trackp2;
TH2D * panx_corr_trackm1_HFp2;
TH2D * panx_corr_trackp1_HFm2;
TH2D * panx_corr_trackp1_trackm2;
TH2D * panx_corr_trackp1_HFp2;

TH2D * panx_corr_HFm2_HFp2;
TH2D * panx_corr_trackmid2_HFm2;
TH2D * panx_corr_trackmid2_HFp2;
TH2D * panx_corr_trackm2_HFm2;
TH2D * panx_corr_trackm2_trackp2;
TH2D * panx_corr_trackm2_HFp2;
TH2D * panx_corr_trackp2_HFm2;
TH2D * panx_corr_trackp2_HFp2;

void EPperiodic()
{

    int nxbins = corr_HFm1_HFp1->GetNbinsX();
    int nybins = corr_HFm1_HFp1->GetNbinsY();

    double xbinmin1v1 = -TMath::Pi();
    double xbinmax1v1 =  TMath::Pi();
    double ybinmin1v1 = -TMath::Pi();
    double ybinmax1v1 =  TMath::Pi();

    double xbinmin1v2 = -TMath::Pi();
    double xbinmax1v2 =  TMath::Pi();
    double ybinmin1v2 = -0.5*TMath::Pi();
    double ybinmax1v2 =  0.5*TMath::Pi();

    double xbinmin2v2 = -0.5*TMath::Pi();
    double xbinmax2v2 =  0.5*TMath::Pi();
    double ybinmin2v2 = -0.5*TMath::Pi();
    double ybinmax2v2 =  0.5*TMath::Pi();

    panx_corr_HFm1_HFp1 = new TH2D("panx_corr_HFm1_HFp1", "panx_corr_HFm1_HFp1", 300, xbinmin1v1, 5*xbinmax1v1, nybins, ybinmin1v1, ybinmax1v1);
    panx_corr_HFm1_HFp1->SetOption("colz");
    panx_corr_HFm1_HFp1->SetTitle("");
    panx_corr_HFm1_HFp1->SetStats(kFALSE);
    panx_corr_HFm1_HFp1->GetXaxis()->CenterTitle(kTRUE);
    panx_corr_HFm1_HFp1->GetYaxis()->CenterTitle(kTRUE);
    panx_corr_HFm1_HFp1->GetXaxis()->SetTitleSize(xtlsize_pan);
    panx_corr_HFm1_HFp1->GetYaxis()->SetTitleSize(ytlsize_pan);
    panx_corr_HFm1_HFp1->GetXaxis()->SetTitleOffset(xtloffset_pan);
    panx_corr_HFm1_HFp1->GetYaxis()->SetTitleOffset(ytloffset_pan);
    panx_corr_HFm1_HFp1->GetXaxis()->SetLabelSize(xlbsize_pan);
    panx_corr_HFm1_HFp1->GetYaxis()->SetLabelSize(ylbsize_pan);
    panx_corr_HFm1_HFp1->GetXaxis()->SetLabelOffset(xlboffset_pan);
    panx_corr_HFm1_HFp1->GetYaxis()->SetLabelOffset(ylboffset_pan);
    panx_corr_HFm1_HFp1->SetXTitle("#Psi_{1}{HF-}");
    panx_corr_HFm1_HFp1->SetYTitle("#Psi_{1}{HF+}");
    panx_corr_HFm1_HFp1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);

    panx_corr_trackmid1_HFm1 = (TH2D *) panx_corr_HFm1_HFp1->Clone("panx_corr_trackmid1_HFm1");
    panx_corr_trackmid1_HFm1->SetXTitle("#Psi_{1}{trackmid}");
    panx_corr_trackmid1_HFm1->SetYTitle("#Psi_{1}{HF-}");
    panx_corr_trackmid1_HFm1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);

    panx_corr_trackmid1_HFp1 = (TH2D *) panx_corr_HFm1_HFp1->Clone("panx_corr_trackmid1_HFp1");
    panx_corr_trackmid1_HFp1->SetXTitle("#Psi_{1}{trackmid}");
    panx_corr_trackmid1_HFp1->SetYTitle("#Psi_{1}{HF+}");
    panx_corr_trackmid1_HFp1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);

    panx_corr_trackm1_HFm1 = (TH2D *) panx_corr_HFm1_HFp1->Clone("panx_corr_trackm1_HFm1");
    panx_corr_trackm1_HFm1->SetXTitle("#Psi_{1}{track-}");
    panx_corr_trackm1_HFm1->SetYTitle("#Psi_{1}{HF-}");
    panx_corr_trackm1_HFm1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);

    panx_corr_trackm1_trackp1 = (TH2D *) panx_corr_HFm1_HFp1->Clone("panx_corr_trackm1_trackp1");
    panx_corr_trackm1_trackp1->SetXTitle("#Psi_{1}{track-}");
    panx_corr_trackm1_trackp1->SetYTitle("#Psi_{1}{track+}");
    panx_corr_trackm1_trackp1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);

    panx_corr_trackm1_HFp1 = (TH2D *) panx_corr_HFm1_HFp1->Clone("panx_corr_trackm1_HFp1");
    panx_corr_trackm1_HFp1->SetXTitle("#Psi_{1}{track-}");
    panx_corr_trackm1_HFp1->SetYTitle("#Psi_{1}{HF+}");
    panx_corr_trackm1_HFp1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);

    panx_corr_trackp1_HFm1 = (TH2D *) panx_corr_HFm1_HFp1->Clone("panx_corr_trackp1_HFm1");
    panx_corr_trackp1_HFm1->SetXTitle("#Psi_{1}{track+}");
    panx_corr_trackp1_HFm1->SetYTitle("#Psi_{1}{HF-}");
    panx_corr_trackp1_HFm1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);

    panx_corr_trackp1_trackm1 = (TH2D *) panx_corr_HFm1_HFp1->Clone("panx_corr_trackp1_trackm1");
    panx_corr_trackp1_trackm1->SetXTitle("#Psi_{1}{track+}");
    panx_corr_trackp1_trackm1->SetYTitle("#Psi_{1}{track-}");
    panx_corr_trackp1_trackm1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);

    panx_corr_trackp1_HFp1 = (TH2D *) panx_corr_HFm1_HFp1->Clone("panx_corr_trackp1_HFp1");
    panx_corr_trackp1_HFp1->SetXTitle("#Psi_{1}{track+}");
    panx_corr_trackp1_HFp1->SetYTitle("#Psi_{1}{HF+}");
    panx_corr_trackp1_HFp1->GetZaxis()->SetRangeUser(zcolmin1v1,zcolmax1v1);



    panx_corr_HFm1_trackm2 = new TH2D("panx_corr_HFm1_trackm2", "panx_corr_HFm1_trackm2", 300, xbinmin1v2, 5*xbinmax1v2, nybins, ybinmin1v2, ybinmax1v2);
    panx_corr_HFm1_trackm2->SetOption("colz");
    panx_corr_HFm1_trackm2->SetTitle("");
    panx_corr_HFm1_trackm2->SetStats(kFALSE);
    panx_corr_HFm1_trackm2->GetXaxis()->CenterTitle(kTRUE);
    panx_corr_HFm1_trackm2->GetYaxis()->CenterTitle(kTRUE);
    panx_corr_HFm1_trackm2->GetXaxis()->SetTitleSize(xtlsize_pan);
    panx_corr_HFm1_trackm2->GetYaxis()->SetTitleSize(ytlsize_pan);
    panx_corr_HFm1_trackm2->GetXaxis()->SetTitleOffset(xtloffset_pan);
    panx_corr_HFm1_trackm2->GetYaxis()->SetTitleOffset(ytloffset_pan);
    panx_corr_HFm1_trackm2->GetXaxis()->SetLabelSize(xlbsize_pan);
    panx_corr_HFm1_trackm2->GetYaxis()->SetLabelSize(ylbsize_pan);
    panx_corr_HFm1_trackm2->GetXaxis()->SetLabelOffset(xlboffset_pan);
    panx_corr_HFm1_trackm2->GetYaxis()->SetLabelOffset(ylboffset_pan);
    panx_corr_HFm1_trackm2->SetXTitle("#Psi_{1}{HF-}");
    panx_corr_HFm1_trackm2->SetYTitle("#Psi_{2}{track-}");
    panx_corr_HFm1_trackm2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_HFm1_trackp2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_HFm1_trackp2");
    panx_corr_HFm1_trackp2->SetXTitle("#Psi_{1}{HF-}");
    panx_corr_HFm1_trackp2->SetYTitle("#Psi_{2}{track+}");
    panx_corr_HFm1_trackp2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_HFp1_trackm2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_HFp1_trackm2");
    panx_corr_HFp1_trackm2->SetXTitle("#Psi_{1}{HF+}");
    panx_corr_HFp1_trackm2->SetYTitle("#Psi_{2}{track-}");
    panx_corr_HFp1_trackm2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_HFp1_trackp2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_HFp1_trackp2");
    panx_corr_HFp1_trackp2->SetXTitle("#Psi_{1}{HF+}");
    panx_corr_HFp1_trackp2->SetYTitle("#Psi_{2}{track+}");
    panx_corr_HFp1_trackp2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_HFm1_HFp2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_HFm1_HFp2");
    panx_corr_HFm1_HFp2->SetXTitle("#Psi_{1}{HF-}");
    panx_corr_HFm1_HFp2->SetYTitle("#Psi_{2}{HF+}");
    panx_corr_HFm1_HFp2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_HFp1_HFm2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_HFp1_HFm2");
    panx_corr_HFp1_HFm2->SetXTitle("#Psi_{1}{HF+}");
    panx_corr_HFp1_HFm2->SetYTitle("#Psi_{2}{HF-}");
    panx_corr_HFp1_HFm2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_trackmid1_HFm2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_trackmid1_HFm2");
    panx_corr_trackmid1_HFm2->SetXTitle("#Psi_{1}{trackmid}");
    panx_corr_trackmid1_HFm2->SetYTitle("#Psi_{2}{HF-}");
    panx_corr_trackmid1_HFm2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_trackmid1_HFp2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_trackmid1_HFp2");
    panx_corr_trackmid1_HFp2->SetXTitle("#Psi_{1}{trackmid}");
    panx_corr_trackmid1_HFp2->SetYTitle("#Psi_{2}{HF+}");
    panx_corr_trackmid1_HFp2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_trackm1_HFm2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_trackm1_HFm2");
    panx_corr_trackm1_HFm2->SetXTitle("#Psi_{1}{track-}");
    panx_corr_trackm1_HFm2->SetYTitle("#Psi_{2}{HF-}");
    panx_corr_trackm1_HFm2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_trackm1_trackp2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_trackm1_trackp2");
    panx_corr_trackm1_trackp2->SetXTitle("#Psi_{1}{track-}");
    panx_corr_trackm1_trackp2->SetYTitle("#Psi_{2}{track+}");
    panx_corr_trackm1_trackp2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_trackm1_HFp2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_trackm1_HFp2");
    panx_corr_trackm1_HFp2->SetXTitle("#Psi_{1}{track-}");
    panx_corr_trackm1_HFp2->SetYTitle("#Psi_{2}{HF+}");
    panx_corr_trackm1_HFp2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_trackp1_HFm2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_trackp1_HFm2");
    panx_corr_trackp1_HFm2->SetXTitle("#Psi_{1}{track+}");
    panx_corr_trackp1_HFm2->SetYTitle("#Psi_{2}{HF-}");
    panx_corr_trackp1_HFm2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_trackp1_trackm2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_trackp1_trackm2");
    panx_corr_trackp1_trackm2->SetXTitle("#Psi_{1}{track+}");
    panx_corr_trackp1_trackm2->SetYTitle("#Psi_{2}{track-}");
    panx_corr_trackp1_trackm2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);

    panx_corr_trackp1_HFp2 = (TH2D *) panx_corr_HFm1_trackm2->Clone("panx_corr_trackp1_HFp2");
    panx_corr_trackp1_HFp2->SetXTitle("#Psi_{1}{track+}");
    panx_corr_trackp1_HFp2->SetYTitle("#Psi_{2}{tracHF+}");
    panx_corr_trackp1_HFp2->GetZaxis()->SetRangeUser(zcolmin1v2,zcolmax1v2);



    panx_corr_HFm2_HFp2 = new TH2D("panx_corr_HFm2_HFp2", "panx_corr_HFm2_HFp2", 300, xbinmin2v2, 5*xbinmax2v2, nybins, ybinmin2v2, ybinmax2v2);
    panx_corr_HFm2_HFp2->SetOption("colz");
    panx_corr_HFm2_HFp2->SetTitle("");
    panx_corr_HFm2_HFp2->SetStats(kFALSE);
    panx_corr_HFm2_HFp2->GetXaxis()->CenterTitle(kTRUE);
    panx_corr_HFm2_HFp2->GetYaxis()->CenterTitle(kTRUE);
    panx_corr_HFm2_HFp2->GetXaxis()->SetTitleSize(xtlsize_pan);
    panx_corr_HFm2_HFp2->GetYaxis()->SetTitleSize(ytlsize_pan);
    panx_corr_HFm2_HFp2->GetXaxis()->SetTitleOffset(xtloffset_pan);
    panx_corr_HFm2_HFp2->GetYaxis()->SetTitleOffset(ytloffset_pan);
    panx_corr_HFm2_HFp2->GetXaxis()->SetLabelSize(xlbsize_pan);
    panx_corr_HFm2_HFp2->GetYaxis()->SetLabelSize(ylbsize_pan);
    panx_corr_HFm2_HFp2->GetXaxis()->SetLabelOffset(xlboffset_pan);
    panx_corr_HFm2_HFp2->GetYaxis()->SetLabelOffset(ylboffset_pan);
    panx_corr_HFm2_HFp2->SetXTitle("#Psi_{2}{HF-}");
    panx_corr_HFm2_HFp2->SetYTitle("#Psi_{2}{HF+}");
    panx_corr_HFm2_HFp2->GetZaxis()->SetRangeUser(zcolmin2v2,zcolmax2v2);

    panx_corr_trackmid2_HFm2 = (TH2D *) panx_corr_HFm2_HFp2->Clone("panx_corr_trackmid2_HFm2");
    panx_corr_trackmid2_HFm2->SetXTitle("#Psi_{2}{trackmid}");
    panx_corr_trackmid2_HFm2->SetYTitle("#Psi_{2}{HF-}");
    panx_corr_trackmid2_HFm2->GetZaxis()->SetRangeUser(zcolmin2v2,zcolmax2v2);

    panx_corr_trackmid2_HFp2 = (TH2D *) panx_corr_HFm2_HFp2->Clone("panx_corr_trackmid2_HFp2");
    panx_corr_trackmid2_HFp2->SetXTitle("#Psi_{2}{trackmid}");
    panx_corr_trackmid2_HFp2->SetYTitle("#Psi_{2}{HF+}");
    panx_corr_trackmid2_HFp2->GetZaxis()->SetRangeUser(zcolmin2v2,zcolmax2v2);

    panx_corr_trackm2_HFm2 = (TH2D *) panx_corr_HFm2_HFp2->Clone("panx_corr_trackm2_HFm2");
    panx_corr_trackm2_HFm2->SetXTitle("#Psi_{2}{track-}");
    panx_corr_trackm2_HFm2->SetYTitle("#Psi_{2}{HF-}");
    panx_corr_trackm2_HFm2->GetZaxis()->SetRangeUser(zcolmin2v2,zcolmax2v2);

    panx_corr_trackm2_trackp2 = (TH2D *) panx_corr_HFm2_HFp2->Clone("panx_corr_trackm2_trackp2");
    panx_corr_trackm2_trackp2->SetXTitle("#Psi_{2}{track-}");
    panx_corr_trackm2_trackp2->SetYTitle("#Psi_{2}{track+}");
    panx_corr_trackm2_trackp2->GetZaxis()->SetRangeUser(zcolmin2v2,zcolmax2v2);

    panx_corr_trackm2_HFp2 = (TH2D *) panx_corr_HFm2_HFp2->Clone("panx_corr_trackm2_HFp2");
    panx_corr_trackm2_HFp2->SetXTitle("#Psi_{2}{track-}");
    panx_corr_trackm2_HFp2->SetYTitle("#Psi_{2}{HF+}");
    panx_corr_trackm2_HFp2->GetZaxis()->SetRangeUser(zcolmin2v2,zcolmax2v2);

    panx_corr_trackp2_HFm2 = (TH2D *) panx_corr_HFm2_HFp2->Clone("panx_corr_trackp2_HFm2");
    panx_corr_trackp2_HFm2->SetXTitle("#Psi_{2}{track+}");
    panx_corr_trackp2_HFm2->SetYTitle("#Psi_{2}{HF-}");
    panx_corr_trackp2_HFm2->GetZaxis()->SetRangeUser(zcolmin2v2,zcolmax2v2);

    panx_corr_trackp2_HFp2 = (TH2D *) panx_corr_HFm2_HFp2->Clone("panx_corr_trackp2_HFp2");
    panx_corr_trackp2_HFp2->SetXTitle("#Psi_{2}{track+}");
    panx_corr_trackp2_HFp2->SetYTitle("#Psi_{2}{HF+}");
    panx_corr_trackp2_HFp2->GetZaxis()->SetRangeUser(zcolmin2v2,zcolmax2v2);

    for (int np = 0; np<3; np++) {
        for (int xbin = 1; xbin<=nxbins; xbin++) {
            for (int ybin = 1; ybin<=nybins; ybin++) {
                panx_corr_HFm1_HFp1->SetBinContent(xbin+np*100, ybin, corr_HFm1_HFp1->GetBinContent(xbin,ybin));
                panx_corr_HFm1_HFp1->SetBinError(xbin+np*100, ybin, corr_HFm1_HFp1->GetBinError(xbin,ybin));

                panx_corr_trackmid1_HFm1->SetBinContent(xbin+np*100, ybin, corr_trackmid1_HFm1->GetBinContent(xbin,ybin));
                panx_corr_trackmid1_HFm1->SetBinError(xbin+np*100, ybin, corr_trackmid1_HFm1->GetBinError(xbin,ybin));

                panx_corr_trackmid1_HFp1->SetBinContent(xbin+np*100, ybin, corr_trackmid1_HFp1->GetBinContent(xbin,ybin));
                panx_corr_trackmid1_HFp1->SetBinError(xbin+np*100, ybin, corr_trackmid1_HFp1->GetBinError(xbin,ybin));

                panx_corr_trackm1_HFm1->SetBinContent(xbin+np*100, ybin, corr_trackm1_HFm1->GetBinContent(xbin,ybin));
                panx_corr_trackm1_HFm1->SetBinError(xbin+np*100, ybin, corr_trackm1_HFm1->GetBinError(xbin,ybin));

                panx_corr_trackm1_trackp1->SetBinContent(xbin+np*100, ybin, corr_trackm1_trackp1->GetBinContent(xbin,ybin));
                panx_corr_trackm1_trackp1->SetBinError(xbin+np*100, ybin, corr_trackm1_trackp1->GetBinError(xbin,ybin));

                panx_corr_trackm1_HFp1->SetBinContent(xbin+np*100, ybin, corr_trackm1_HFp1->GetBinContent(xbin,ybin));
                panx_corr_trackm1_HFp1->SetBinError(xbin+np*100, ybin, corr_trackm1_HFp1->GetBinError(xbin,ybin));

                panx_corr_trackp1_HFm1->SetBinContent(xbin+np*100, ybin, corr_trackp1_HFm1->GetBinContent(xbin,ybin));
                panx_corr_trackp1_HFm1->SetBinError(xbin+np*100, ybin, corr_trackp1_HFm1->GetBinError(xbin,ybin));

                panx_corr_trackp1_trackm1->SetBinContent(xbin+np*100, ybin, corr_trackp1_trackm1->GetBinContent(xbin,ybin));
                panx_corr_trackp1_trackm1->SetBinError(xbin+np*100, ybin, corr_trackp1_trackm1->GetBinError(xbin,ybin));

                panx_corr_trackp1_HFp1->SetBinContent(xbin+np*100, ybin, corr_trackp1_HFp1->GetBinContent(xbin,ybin));
                panx_corr_trackp1_HFp1->SetBinError(xbin+np*100, ybin, corr_trackp1_HFp1->GetBinError(xbin,ybin));


                panx_corr_HFm1_trackm2->SetBinContent(xbin+np*100, ybin, corr_HFm1_trackm2->GetBinContent(xbin,ybin));
                panx_corr_HFm1_trackm2->SetBinError(xbin+np*100, ybin, corr_HFm1_trackm2->GetBinError(xbin,ybin));

                panx_corr_HFm1_trackp2->SetBinContent(xbin+np*100, ybin, corr_HFm1_trackp2->GetBinContent(xbin,ybin));
                panx_corr_HFm1_trackp2->SetBinError(xbin+np*100, ybin, corr_HFm1_trackp2->GetBinError(xbin,ybin));

                panx_corr_HFp1_trackm2->SetBinContent(xbin+np*100, ybin, corr_HFp1_trackm2->GetBinContent(xbin,ybin));
                panx_corr_HFp1_trackm2->SetBinError(xbin+np*100, ybin, corr_HFp1_trackm2->GetBinError(xbin,ybin));

                panx_corr_HFp1_trackp2->SetBinContent(xbin+np*100, ybin, corr_HFp1_trackp2->GetBinContent(xbin,ybin));
                panx_corr_HFp1_trackp2->SetBinError(xbin+np*100, ybin, corr_HFp1_trackp2->GetBinError(xbin,ybin));

                panx_corr_HFm1_HFp2->SetBinContent(xbin+np*100, ybin, corr_HFm1_HFp2->GetBinContent(xbin,ybin));
                panx_corr_HFm1_HFp2->SetBinError(xbin+np*100, ybin, corr_HFm1_HFp2->GetBinError(xbin,ybin));

                panx_corr_HFp1_HFm2->SetBinContent(xbin+np*100, ybin, corr_HFp1_HFm2->GetBinContent(xbin,ybin));
                panx_corr_HFp1_HFm2->SetBinError(xbin+np*100, ybin, corr_HFm1_HFp1->GetBinError(xbin,ybin));

                panx_corr_trackmid1_HFm2->SetBinContent(xbin+np*100, ybin, corr_trackmid1_HFm2->GetBinContent(xbin,ybin));
                panx_corr_trackmid1_HFm2->SetBinError(xbin+np*100, ybin, corr_trackmid1_HFm2->GetBinError(xbin,ybin));

                panx_corr_trackmid1_HFp2->SetBinContent(xbin+np*100, ybin, corr_HFm1_HFp1->GetBinContent(xbin,ybin));
                panx_corr_trackmid1_HFp2->SetBinError(xbin+np*100, ybin, corr_HFm1_HFp1->GetBinError(xbin,ybin));

                panx_corr_trackm1_HFm2->SetBinContent(xbin+np*100, ybin, corr_trackm1_HFm2->GetBinContent(xbin,ybin));
                panx_corr_trackm1_HFm2->SetBinError(xbin+np*100, ybin, corr_trackm1_HFm2->GetBinError(xbin,ybin));

                panx_corr_trackm1_trackp2->SetBinContent(xbin+np*100, ybin, corr_trackm1_trackp2->GetBinContent(xbin,ybin));
                panx_corr_trackm1_trackp2->SetBinError(xbin+np*100, ybin, corr_trackm1_trackp2->GetBinError(xbin,ybin));

                panx_corr_trackm1_HFp2->SetBinContent(xbin+np*100, ybin, corr_trackm1_HFp2->GetBinContent(xbin,ybin));
                panx_corr_trackm1_HFp2->SetBinError(xbin+np*100, ybin, corr_trackm1_HFp2->GetBinError(xbin,ybin));

                panx_corr_trackp1_HFm2->SetBinContent(xbin+np*100, ybin, corr_trackp1_HFm2->GetBinContent(xbin,ybin));
                panx_corr_trackp1_HFm2->SetBinError(xbin+np*100, ybin, corr_trackp1_HFm2->GetBinError(xbin,ybin));

                panx_corr_trackp1_trackm2->SetBinContent(xbin+np*100, ybin, corr_trackp1_trackm2->GetBinContent(xbin,ybin));
                panx_corr_trackp1_trackm2->SetBinError(xbin+np*100, ybin, corr_trackp1_trackm2->GetBinError(xbin,ybin));

                panx_corr_trackp1_HFp2->SetBinContent(xbin+np*100, ybin, corr_trackp1_HFp2->GetBinContent(xbin,ybin));
                panx_corr_trackp1_HFp2->SetBinError(xbin+np*100, ybin, corr_trackp1_HFp2->GetBinError(xbin,ybin));


                panx_corr_HFm2_HFp2->SetBinContent(xbin+np*100, ybin, corr_HFm2_HFp2->GetBinContent(xbin,ybin));
                panx_corr_HFm2_HFp2->SetBinError(xbin+np*100, ybin, corr_HFm2_HFp2->GetBinError(xbin,ybin));

                panx_corr_trackmid2_HFm2->SetBinContent(xbin+np*100, ybin, corr_trackmid2_HFm2->GetBinContent(xbin,ybin));
                panx_corr_trackmid2_HFm2->SetBinError(xbin+np*100, ybin, corr_trackmid2_HFm2->GetBinError(xbin,ybin));

                panx_corr_trackmid2_HFp2->SetBinContent(xbin+np*100, ybin, corr_trackmid2_HFp2->GetBinContent(xbin,ybin));
                panx_corr_trackmid2_HFp2->SetBinError(xbin+np*100, ybin, corr_trackmid2_HFp2->GetBinError(xbin,ybin));

                panx_corr_trackm2_HFm2->SetBinContent(xbin+np*100, ybin, corr_trackm2_HFm2->GetBinContent(xbin,ybin));
                panx_corr_trackm2_HFm2->SetBinError(xbin+np*100, ybin, corr_trackm2_HFm2->GetBinError(xbin,ybin));

                panx_corr_trackm2_trackp2->SetBinContent(xbin+np*100, ybin, corr_trackm2_trackp2->GetBinContent(xbin,ybin));
                panx_corr_trackm2_trackp2->SetBinError(xbin+np*100, ybin, corr_trackm2_trackp2->GetBinError(xbin,ybin));

                panx_corr_trackm2_HFp2->SetBinContent(xbin+np*100, ybin, corr_trackm2_HFp2->GetBinContent(xbin,ybin));
                panx_corr_trackm2_HFp2->SetBinError(xbin+np*100, ybin, corr_trackm2_HFp2->GetBinError(xbin,ybin));

                panx_corr_trackp2_HFm2->SetBinContent(xbin+np*100, ybin, corr_trackp2_HFm2->GetBinContent(xbin,ybin));
                panx_corr_trackp2_HFm2->SetBinError(xbin+np*100, ybin, corr_trackp2_HFm2->GetBinError(xbin,ybin));

                panx_corr_trackp2_HFp2->SetBinContent(xbin+np*100, ybin, corr_trackp2_HFp2->GetBinContent(xbin,ybin));
                panx_corr_trackp2_HFp2->SetBinError(xbin+np*100, ybin, corr_trackp2_HFp2->GetBinError(xbin,ybin));
            }
        }
    }

}
